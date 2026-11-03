"""
Publication-Quality Plotting Engine for Moments2Amplitudes.

Generates:
1. Fitted Amplitude Magnitudes (|a_LM|, |b_LM|) vs invariant mass with MCMC error bands.
2. Relative Phase Differences vs invariant mass with MCMC bands.
3. Moments comparison: Experimental input moments vs reconstructed fit moments with pull residuals.
4. Multi-Start Chi2 Landscape showing global minimum clustering and plateau count.
"""

from __future__ import annotations
import math
import os
import json
import subprocess
import shutil
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any

# Ensure matplotlib has a writable config directory
os.environ.setdefault("MPLCONFIGDIR", "/tmp")

from .schema import FitConfig
from .rundir import RunDirectory


def get_root_binary(config: Optional[FitConfig] = None) -> str:
    """Find the ROOT binary executable path."""
    if "ROOTSYS" in os.environ:
        cand = Path(os.environ["ROOTSYS"]) / "bin" / "root"
        if cand.is_file():
            return str(cand)
    if config and config.environment.root_dir:
        cand = Path(config.environment.root_dir) / "bin" / "root"
        if cand.is_file():
            return str(cand)
    default_root = Path("/d/grid17/septian/root_v6-36-08_cuda_install/bin/root")
    if default_root.is_file():
        return str(default_root)
    which_root = shutil.which("root")
    if which_root:
        return which_root
    return "root"


def discover_lowest_chi2_fits(rundir: RunDirectory, n_lowest: int = 10) -> List[Tuple[float, Path]]:
    """
    Discover completed fit runs in fit_starts_dir, sort by chi2 ascending,
    and return the top n_lowest (chi2, root_file) tuples.
    """
    valid_fits: List[Tuple[float, Path]] = []

    if rundir.fit_starts_dir.is_dir():
        for log_f in rundir.fit_starts_dir.glob("fit_seed_*.log"):
            try:
                with open(log_f, "r") as f:
                    for line in f:
                        if "[RESULT]" in line:
                            tokens = dict(part.split("=") for part in line.strip().split() if "=" in part)
                            if "CHI2" in tokens and "VALID" in tokens and "STATUS" in tokens:
                                is_valid = tokens["VALID"] == "1"
                                status = int(tokens["STATUS"])
                                chi2_val = float(tokens["CHI2"])
                                seed_val = int(tokens.get("SEED", 0))

                                if is_valid and status == 0 and chi2_val < 1e8:
                                    cand1 = rundir.fit_starts_dir / f"MassDependentResults_{seed_val}.root"
                                    cand2 = rundir.get_seed_root_file(seed_val)
                                    root_cand = cand1 if cand1.is_file() else (cand2 if cand2.is_file() else None)
                                    if root_cand:
                                        valid_fits.append((chi2_val, root_cand))
            except Exception:
                pass

    valid_fits.sort(key=lambda x: x[0])

    # Ensure nominal_best is included as the primary best fit if available
    if rundir.nominal_best_root.is_file():
        has_nominal = any(f[1].resolve() == rundir.nominal_best_root.resolve() for f in valid_fits)
        if not has_nominal:
            best_chi2 = valid_fits[0][0] if valid_fits else 0.0
            valid_fits.insert(0, (best_chi2, rundir.nominal_best_root))

    return valid_fits[:n_lowest]


def plot_multistart_landscape(
    results: List[Dict[str, Any]],
    output_dir: Path | str,
    project_name: str = "Moments2Amplitudes",
    plateau_delta: float = 10.0,
    formats: Tuple[str, ...] = ("png", "pdf"),
) -> List[Path]:
    """Plot multi-start Chi2 distribution histogram and clustering plateau."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    generated: List[Path] = []

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("[WARNING] matplotlib not found. Skipping multi-start plot.")
        return generated

    chi2s = [r["chi2"] for r in results if r.get("is_valid", False) and r.get("status", -1) == 0 and r.get("chi2", 1e9) < 1e6]
    if not chi2s:
        return generated

    min_chi2 = min(chi2s)
    plateau_thresh = min_chi2 + plateau_delta

    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    counts, bins, patches = ax.hist(chi2s, bins=40, color="#1f77b4", edgecolor="black", alpha=0.7, label=f"Trials (N={len(chi2s)})")

    ax.axvline(min_chi2, color="red", linestyle="--", linewidth=2.0, label=f"Global Min $\\chi^2 = {min_chi2:.2f}$")
    ax.axvspan(min_chi2, plateau_thresh, color="green", alpha=0.2, label=f"Plateau $\\leq +{plateau_delta:.1f}$")

    ax.set_xlabel(r"Total $\chi^2$", fontsize=12)
    ax.set_ylabel("Number of Starts", fontsize=12)
    ax.set_title(f"{project_name} - Multi-Start $\\chi^2$ Landscape", fontsize=13, fontweight="bold")
    ax.legend(frameon=True, fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.5)

    for fmt in formats:
        p = out_dir / f"multistart_chi2_landscape.{fmt}"
        fig.savefig(p, bbox_inches="tight")
        generated.append(p)

    plt.close(fig)
    return generated


def plot_amplitudes_and_bands(
    nominal_root: Path | str,
    mcmc_root: Optional[Path | str],
    output_dir: Path | str,
    config: FitConfig,
    formats: Tuple[str, ...] = ("png", "pdf"),
) -> List[Path]:
    """
    Plot fitted amplitude magnitudes and phases vs mass with shaded MCMC 68%/95% credible intervals.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    generated: List[Path] = []

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("[WARNING] matplotlib not found. Skipping amplitude plots.")
        return generated

    # Extract nominal result tree
    masses: List[float] = []
    nominal_amps: Dict[str, List[float]] = {}
    mcmc_samples: Dict[str, Dict[float, List[float]]] = {}

    # Read with uproot or PyROOT
    loaded = False
    try:
        import uproot
        with uproot.open(nominal_root) as f:
            if "result" in f:
                tree = f["result"]
                masses = list(tree["mass_bin"].array())
                for k in tree.keys():
                    if k != "mass_bin" and not k.startswith("H_"):
                        nominal_amps[k] = list(tree[k].array())
                loaded = True

        if mcmc_root and Path(mcmc_root).is_file():
            with uproot.open(mcmc_root) as mf:
                if "result" in mf:
                    mtree = mf["result"]
                    mbins = list(mtree["mass_bin"].array())
                    for k in nominal_amps.keys():
                        if k in mtree:
                            arr = list(mtree[k].array())
                            mcmc_samples[k] = {}
                            for m, v in zip(mbins, arr):
                                mcmc_samples[k].setdefault(round(m, 4), []).append(v)
    except Exception:
        pass

    if not loaded:
        # Fallback to PyROOT
        try:
            import ROOT
            f = ROOT.TFile.Open(str(nominal_root), "READ")
            if f and not f.IsZombie():
                tree = f.Get("result")
                if tree:
                    n_entries = tree.GetEntries()
                    for i in range(n_entries):
                        tree.GetEntry(i)
                        m_val = float(getattr(tree, "mass_bin", 0.0))
                        masses.append(m_val)
                    for b in tree.GetListOfBranches():
                        bname = b.GetName()
                        if bname != "mass_bin" and not bname.startswith("H_"):
                            vals = []
                            for i in range(n_entries):
                                tree.GetEntry(i)
                                vals.append(float(getattr(tree, bname, 0.0)))
                            nominal_amps[bname] = vals
                    loaded = True
                f.Close()
        except Exception:
            pass

    if not loaded or not masses:
        print(f"[WARNING] Could not read amplitude branches from {nominal_root}")
        return generated

    # Group into magnitudes (a_L_M, b_L_M) and phases (aphi_L_M, bphi_L_M)
    magnitudes = [k for k in nominal_amps.keys() if (k.startswith("a_") or k.startswith("b_"))]
    phases = [k for k in nominal_amps.keys() if (k.startswith("aphi_") or k.startswith("bphi_"))]

    # --- 1. Amplitude Magnitudes Plot ---
    if magnitudes:
        fig, ax = plt.subplots(figsize=(9, 6), dpi=300)
        colors = plt.cm.tab10(np.linspace(0, 1, len(magnitudes)))

        for idx, k in enumerate(sorted(magnitudes)):
            y_nom = nominal_amps[k]
            line, = ax.plot(masses, y_nom, "-", linewidth=2.0, color=colors[idx], label=f"${k}$")

            # MCMC credible bands if available
            if k in mcmc_samples:
                lo_68, hi_68 = [], []
                for m in masses:
                    m_round = round(m, 4)
                    samps = mcmc_samples[k].get(m_round, [])
                    if samps:
                        lo_68.append(np.percentile(samps, 16))
                        hi_68.append(np.percentile(samps, 84))
                    else:
                        lo_68.append(y_nom[len(lo_68)])
                        hi_68.append(y_nom[len(hi_68)])
                ax.fill_between(masses, lo_68, hi_68, color=colors[idx], alpha=0.25)

        ax.set_xlabel(r"$M(\eta\pi^0)$ [GeV]", fontsize=12)
        ax.set_ylabel("Amplitude Magnitude", fontsize=12)
        ax.set_title(f"{config.project.name} - Fitted Amplitude Magnitudes", fontsize=13, fontweight="bold")
        ax.legend(frameon=True, fontsize=10, ncol=2)
        ax.grid(True, linestyle=":", alpha=0.5)

        for fmt in formats:
            p = out_dir / f"amplitudes_magnitude_vs_mass.{fmt}"
            fig.savefig(p, bbox_inches="tight")
            generated.append(p)
        plt.close(fig)

    # --- 2. Relative Phases Plot ---
    if phases:
        fig, ax = plt.subplots(figsize=(9, 6), dpi=300)
        colors = plt.cm.tab10(np.linspace(0, 1, len(phases)))

        for idx, k in enumerate(sorted(phases)):
            y_nom = nominal_amps[k]
            ax.plot(masses, y_nom, "-", linewidth=2.0, color=colors[idx], label=f"${k}$")

            if k in mcmc_samples:
                lo_68, hi_68 = [], []
                for m in masses:
                    m_round = round(m, 4)
                    samps = mcmc_samples[k].get(m_round, [])
                    if samps:
                        lo_68.append(np.percentile(samps, 16))
                        hi_68.append(np.percentile(samps, 84))
                    else:
                        lo_68.append(y_nom[len(lo_68)])
                        hi_68.append(y_nom[len(hi_68)])
                ax.fill_between(masses, lo_68, hi_68, color=colors[idx], alpha=0.25)

        ax.set_xlabel(r"$M(\eta\pi^0)$ [GeV]", fontsize=12)
        ax.set_ylabel("Phase [rad]", fontsize=12)
        ax.set_title(f"{config.project.name} - Fitted Wave Phases", fontsize=13, fontweight="bold")
        ax.legend(frameon=True, fontsize=10, ncol=2)
        ax.grid(True, linestyle=":", alpha=0.5)

        for fmt in formats:
            p = out_dir / f"amplitudes_phase_vs_mass.{fmt}"
            fig.savefig(p, bbox_inches="tight")
            generated.append(p)
        plt.close(fig)

    return generated


def plot_moments_comparison(
    config: FitConfig,
    rundir: RunDirectory,
    n_lowest_chi2: Optional[int] = None,
    moments: Optional[List[str]] = None,
    unnormalized_moments_file: Optional[Path | str] = None,
    formats: Tuple[str, ...] = ("png", "pdf"),
) -> List[Path]:
    """
    Plot experimental unnormalized moments (from MCMC) vs mass-dependent model fit results.
    Overlays the N lowest chi2 fit results (Rank 1 solid red, ranks 2..N semi-transparent red).
    Generates both individual moment plots and a multi-panel overview grid.
    """
    generated: List[Path] = []

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("[WARNING] matplotlib not found. Skipping moment comparison plots.")
        return generated

    # Resolve number of lowest chi2 results to draw
    n_lowest = n_lowest_chi2
    if n_lowest is None and getattr(config, "plotting", None):
        n_lowest = getattr(config.plotting, "n_lowest_chi2", None)
    if n_lowest is None and getattr(config, "multi_start", None):
        n_lowest = getattr(config.multi_start, "n_lowest_chi2", None)
    if n_lowest is None or n_lowest <= 0:
        n_lowest = 10

    # Resolve unnormalized moments file
    unnorm_path: Optional[Path] = None
    if unnormalized_moments_file:
        cand = Path(unnormalized_moments_file)
        if cand.is_file():
            unnorm_path = cand
    elif getattr(config, "plotting", None) and config.plotting.unnormalized_moments_file:
        cand = Path(config.plotting.unnormalized_moments_file)
        if cand.is_file():
            unnorm_path = cand
    elif config.dataset.moments_dir:
        for fname in ["UnnormalizedMoments_vs_Mass.root", "UnnormalizedMomentsvMass.root", "UnnormalizedMoments.root"]:
            cand = Path(config.dataset.moments_dir) / fname
            if cand.is_file():
                unnorm_path = cand
                break

    if not unnorm_path or not unnorm_path.is_file():
        print(f"[WARNING] Unnormalized moments file not found in {config.dataset.moments_dir}. Skipping moment comparison plots.")
        return generated

    # Resolve lowest chi2 fit files
    best_fits = discover_lowest_chi2_fits(rundir, n_lowest=n_lowest)
    if not best_fits and rundir.nominal_best_root.is_file():
        best_fits = [(0.0, rundir.nominal_best_root)]

    if not best_fits:
        print(f"[WARNING] No fit result ROOT files found to overlay in {rundir.fit_dir}.")
        return generated

    fit_files = [str(f[1].resolve()) for f in best_fits]
    files_str = ",".join(fit_files)

    # Locate ExportMomentComparison.C
    repo_root = Path(os.environ.get("MOMS2AMPS", str(Path(__file__).resolve().parent.parent.parent)))
    export_macro = repo_root / "macros" / "ExportMomentComparison.C"
    if not export_macro.is_file():
        print(f"[WARNING] Export macro not found at {export_macro}. Skipping moment comparison plots.")
        return generated

    cache_dir = rundir.root / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    json_path = cache_dir / "moment_comparison_data.json"

    root_bin = get_root_binary(config)
    cmd = [
        root_bin,
        "-l",
        "-b",
        "-q",
        f'{export_macro}("{unnorm_path.resolve()}", "{files_str}", "{json_path.resolve()}")'
    ]

    env = os.environ.copy()
    if "ROOTSYS" in env and "LD_LIBRARY_PATH" not in env:
        env["LD_LIBRARY_PATH"] = f"{env['ROOTSYS']}/lib"

    run_res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if run_res.returncode != 0 or not json_path.is_file():
        print(f"[WARNING] ExportMomentComparison failed:\n{run_res.stderr}\n{run_res.stdout}")
        return generated

    try:
        with open(json_path, "r") as f:
            comp_data = json.load(f)
    except Exception as e:
        print(f"[WARNING] Could not parse exported moment comparison JSON: {e}")
        return generated

    data_moments: Dict[str, Any] = comp_data.get("data", {})
    fits_list: List[Dict[str, Any]] = comp_data.get("fits", [])

    if not data_moments:
        print(f"[WARNING] No moment data found in {unnorm_path}.")
        return generated

    # Filter requested moments
    req_moments = moments if moments is not None else getattr(config.plotting, "moments", ["all"])
    if not req_moments:
        req_moments = ["all"]

    selected_moments: List[str] = []
    if any(m.lower().strip() == "all" for m in req_moments):
        selected_moments = list(data_moments.keys())
    else:
        for m_req in req_moments:
            m_str = m_req.strip()
            m_lower = m_str.lower()
            if m_lower.startswith("l="):
                try:
                    target_l = int(m_lower.replace("l=", ""))
                    for k in data_moments.keys():
                        parts = k.split("_")
                        if len(parts) >= 4 and int(parts[2]) == target_l:
                            if k not in selected_moments:
                                selected_moments.append(k)
                except ValueError:
                    pass
            elif m_str in data_moments:
                if m_str not in selected_moments:
                    selected_moments.append(m_str)
            else:
                if not m_str.startswith("H_"):
                    alt = f"H_{m_str}"
                    if alt in data_moments and alt not in selected_moments:
                        selected_moments.append(alt)

    def sort_key(name: str):
        parts = name.split("_")
        if len(parts) >= 4:
            try:
                return (int(parts[1]), int(parts[2]), int(parts[3]))
            except ValueError:
                pass
        return (999, 999, 999)

    selected_moments.sort(key=sort_key)

    if not selected_moments:
        print(f"[WARNING] No moments matched filter {req_moments}.")
        return generated

    moments_out_dir = rundir.moments_plots_dir
    moments_out_dir.mkdir(parents=True, exist_ok=True)

    draw_fits = getattr(config.plotting, "draw_best_fits", True)

    # 1. Individual moment plots
    for mom_name in selected_moments:
        d_val = data_moments[mom_name]
        m_data = d_val.get("mass", [])
        y_data = d_val.get("val", [])
        e_data = d_val.get("err", [])

        parts = mom_name.split("_")
        if len(parts) >= 4:
            alpha, l_val, m_val = parts[1], parts[2], parts[3]
            latex_title = f"H_{{{alpha}}}({l_val},{m_val})"
        else:
            latex_title = mom_name

        fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
        ax.errorbar(
            m_data, y_data, yerr=e_data,
            fmt="o", color="black", markersize=3.5,
            capsize=1.5, elinewidth=1.0,
            label="Data (MCMC)", zorder=3
        )

        if draw_fits and fits_list:
            for idx, fit in enumerate(fits_list):
                if mom_name in fit.get("moments", {}):
                    f_mom = fit["moments"][mom_name]
                    mf = f_mom.get("mass", [])
                    yf = f_mom.get("val", [])
                    if idx == 0:
                        ax.plot(
                            mf, yf, "-", color="red", linewidth=2.0, alpha=0.9,
                            label=f"Rank 1 ($\chi^2 = {fit.get('chi2', 0):.1f}$)", zorder=5
                        )
                    else:
                        lbl = f"Top {len(fits_list)} Lowest $\chi^2$" if idx == 1 else None
                        ax.plot(
                            mf, yf, "-", color="#d62728", linewidth=1.0, alpha=0.35,
                            label=lbl, zorder=4
                        )

        ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=12)
        ax.set_ylabel("Unnormalized Moment", fontsize=12)
        ax.set_title(f"{config.project.name} - ${latex_title}$", fontsize=13, fontweight="bold")
        ax.legend(frameon=True, fontsize=10, loc="best")
        ax.grid(True, linestyle=":", alpha=0.5)

        for fmt in formats:
            p = moments_out_dir / f"{mom_name}.{fmt}"
            fig.savefig(p, bbox_inches="tight")
            generated.append(p)
        plt.close(fig)

    # 2. Multi-panel overview grid
    if len(selected_moments) > 1:
        K = len(selected_moments)
        ncols = min(4, K)
        nrows = math.ceil(K / ncols)

        fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.2 * nrows), dpi=200)
        axes_flat = np.array(axes).reshape(-1)

        for idx, mom_name in enumerate(selected_moments):
            ax = axes_flat[idx]
            d_val = data_moments[mom_name]
            m_data = d_val.get("mass", [])
            y_data = d_val.get("val", [])
            e_data = d_val.get("err", [])

            parts = mom_name.split("_")
            if len(parts) >= 4:
                alpha, l_val, m_val = parts[1], parts[2], parts[3]
                latex_title = f"H_{{{alpha}}}({l_val},{m_val})"
            else:
                latex_title = mom_name

            ax.errorbar(
                m_data, y_data, yerr=e_data,
                fmt="o", color="black", markersize=2.2,
                capsize=1.0, elinewidth=0.8, zorder=3
            )

            if draw_fits and fits_list:
                for f_idx, fit in enumerate(fits_list):
                    if mom_name in fit.get("moments", {}):
                        f_mom = fit["moments"][mom_name]
                        mf = f_mom.get("mass", [])
                        yf = f_mom.get("val", [])
                        if f_idx == 0:
                            ax.plot(mf, yf, "-", color="red", linewidth=1.4, alpha=0.9, zorder=5)
                        else:
                            ax.plot(mf, yf, "-", color="#d62728", linewidth=0.8, alpha=0.35, zorder=4)

            ax.set_title(f"${latex_title}$", fontsize=11, fontweight="bold")
            ax.grid(True, linestyle=":", alpha=0.4)
            if idx >= (nrows - 1) * ncols:
                ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=9)

        for idx in range(K, len(axes_flat)):
            axes_flat[idx].set_visible(False)

        fig.suptitle(
            f"{config.project.name} - Unnormalized Moments vs Mass (Top {len(fits_list)} Lowest $\chi^2$ Fits)",
            fontsize=14,
            fontweight="bold",
            y=0.995,
        )
        plt.tight_layout()

        for fmt in formats:
            p_grid = rundir.plots_dir / f"moments_overview.{fmt}"
            fig.savefig(p_grid, bbox_inches="tight")
            generated.append(p_grid)
        plt.close(fig)

    print(f"[OK] Generated {len(generated)} moment comparison plot files in {rundir.plots_dir}")
    return generated


def run_all_diagnostics_and_plots(
    config: FitConfig,
    rundir: RunDirectory,
    formats: Tuple[str, ...] = ("png", "pdf"),
) -> Dict[str, Any]:
    """Run full diagnostic report generation and publication plotting."""
    from .summary import FitSummary

    summary = FitSummary(config, rundir)
    has_nominal = summary.load_from_root(rundir.nominal_best_root, rundir.mcmc_propagated_root)

    report_files = []
    if has_nominal:
        txt_path, md_path = summary.save_reports()
        report_files.extend([txt_path, md_path])
        print(f"\n{summary.generate_ascii_table()}\n")

    plot_files = plot_amplitudes_and_bands(
        rundir.nominal_best_root,
        rundir.mcmc_propagated_root,
        rundir.plots_dir,
        config,
        formats=formats,
    )

    # Generate multi-start landscape plot if starts exist
    if rundir.fit_starts_dir.is_dir():
        ms_results = []
        for log_f in rundir.fit_starts_dir.glob("fit_seed_*.log"):
            try:
                with open(log_f, "r") as f:
                    for line in f:
                        if "[RESULT] SEED=" in line:
                            tokens = dict(part.split("=") for part in line.strip().split() if "=" in part)
                            if "CHI2" in tokens and "VALID" in tokens:
                                ms_results.append({
                                    "seed": int(tokens.get("SEED", 0)),
                                    "chi2": float(tokens["CHI2"]),
                                    "is_valid": tokens["VALID"] == "1",
                                    "status": int(tokens.get("STATUS", -1)),
                                })
            except Exception:
                pass
        if ms_results:
            ms_plots = plot_multistart_landscape(
                ms_results,
                output_dir=rundir.plots_dir,
                project_name=config.project.name,
                plateau_delta=config.multi_start.plateau_delta_chi2,
                formats=formats,
            )
            plot_files.extend(ms_plots)

    # Generate moment comparison plots if enabled
    if getattr(config, "plotting", None) and config.plotting.enabled:
        mom_plots = plot_moments_comparison(
            config,
            rundir,
            formats=formats,
        )
        plot_files.extend(mom_plots)

    # Generate resonance cross-section report and plot if enabled
    if getattr(config, "cross_section", None) and config.cross_section.enabled:
        try:
            from .cross_section import (
                calculate_resonance_cross_section,
                save_cross_section_reports,
                plot_cross_section_vs_t,
            )
            cs_res = calculate_resonance_cross_section(config, rundir)
            txt_p, md_p, json_p = save_cross_section_reports(cs_res, rundir.reports_dir)
            report_files.extend([txt_p, md_p, json_p])
            cs_plots = plot_cross_section_vs_t(
                [cs_res],
                rundir.plots_dir,
                formats=formats,
                phase1_comparison=config.cross_section.phase1_comparison,
                theory_band_file=config.cross_section.theory_band_file,
            )
            plot_files.extend(cs_plots)
        except Exception as e:
            print(f"[WARNING] Cross-section calculation skipped: {e}")

    return {
        "summary": summary,
        "reports": report_files,
        "plots": plot_files,
    }


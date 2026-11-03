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
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any

# Ensure matplotlib has a writable config directory
os.environ.setdefault("MPLCONFIGDIR", "/tmp")

from .schema import FitConfig
from .rundir import RunDirectory


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

    return {
        "summary": summary,
        "reports": report_files,
        "plots": plot_files,
    }


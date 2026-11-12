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
from concurrent.futures import ProcessPoolExecutor, as_completed

# Ensure matplotlib has a writable config directory
os.environ.setdefault("MPLCONFIGDIR", "/tmp")

from .schema import FitConfig
from .rundir import RunDirectory
from .environment import (
    get_root_binary,
    get_root_subprocess_env,
    resolve_rootsys,
    resolve_brufit,
    is_pixi_or_conda_path,
    get_default_max_workers,
)


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


def format_amplitude_label(name: str) -> str:
    """Format amplitude names into valid LaTeX mathtext (e.g. a_2_-1 -> $a_{2,-1}$)."""
    parts = name.split("_")
    if len(parts) >= 3:
        prefix = parts[0]
        l_str = parts[1]
        m_str = "_".join(parts[2:])
        if prefix in ("a", "b"):
            return f"${prefix}_{{{l_str},{m_str}}}$"
        elif prefix in ("aphi", "bphi"):
            base = "a" if prefix == "aphi" else "b"
            return f"$\\phi({base}_{{{l_str},{m_str}}})$"
    return name.replace("_", "\\_")


def smooth_curve(
    x: List[float] | Any,
    y: List[float] | Any,
    num_points: int = 200,
) -> Tuple[Any, Any]:
    """
    Interpolate (x, y) into a dense, smooth curve for plotting.
    Uses CubicSpline if scipy is available and len(x) >= 4.
    Preserves original points if interpolation fails or has insufficient points.
    """
    import numpy as np

    if len(x) < 4:
        return np.array(x), np.array(y)
    try:
        from scipy.interpolate import CubicSpline
        x_arr = np.array(x, dtype=float)
        y_arr = np.array(y, dtype=float)
        sort_idx = np.argsort(x_arr)
        xs = x_arr[sort_idx]
        ys = y_arr[sort_idx]
        # Remove exact duplicate points along x
        unique_mask = np.concatenate(([True], np.diff(xs) > 1e-7))
        xs = xs[unique_mask]
        ys = ys[unique_mask]
        if len(xs) < 4:
            return x_arr, y_arr
        cs = CubicSpline(xs, ys, bc_type="natural")
        x_dense = np.linspace(xs[0], xs[-1], num_points)
        y_dense = cs(x_dense)
        return x_dense, y_dense
    except Exception:
        return np.array(x), np.array(y)


def is_moment_fitted(mom_name: str, config: FitConfig) -> bool:
    """Check if a moment was included in the fit based on chi2 config."""
    chi2_cfg = getattr(config, "chi2", None)
    if not chi2_cfg:
        return True

    # 1. Check no_use patterns (e.g. 'H_3' excludes all L=3 moments)
    no_use = getattr(chi2_cfg, "no_use", []) or []
    for pat in no_use:
        if pat and (pat in mom_name or mom_name.startswith(pat)):
            return False

    parts = mom_name.split("_")
    if len(parts) >= 4 and parts[0] == "H":
        try:
            alpha, L, M = int(parts[1]), int(parts[2]), int(parts[3])
        except ValueError:
            return True
    elif len(parts) >= 3:
        try:
            alpha, L, M = int(parts[0]), int(parts[1]), int(parts[2])
        except ValueError:
            return True
    else:
        return True

    # Also check no_use for L-prefix match (e.g. pat == 'H_3' and L == 3)
    for pat in no_use:
        if pat and pat.startswith("H_") and len(pat.split("_")) == 2:
            try:
                if int(pat.split("_")[1]) == L:
                    return False
            except ValueError:
                pass

    if getattr(chi2_cfg, "include_all", False):
        return True

    included_l = getattr(chi2_cfg, "included_l", None)
    if included_l is not None and len(included_l) > 0:
        if L not in included_l:
            return False

    excluded_l = getattr(chi2_cfg, "excluded_l", None) or []
    if L in excluded_l:
        return False

    max_l = getattr(chi2_cfg, "max_l", None)
    if max_l is not None and L > max_l:
        return False

    return True


def extract_amplitude_l(name: str) -> Optional[int]:
    """Extract angular momentum L from amplitude branch name (e.g. a_2_-1 -> 2, aphi_0_0 -> 0)."""
    parts = name.split("_")
    if len(parts) >= 3 and parts[0] in ("a", "b", "aphi", "bphi"):
        try:
            return int(parts[1])
        except ValueError:
            return None
    return None


def is_amplitude_mass_dependent(name: str, config: FitConfig) -> bool:
    """
    Return True if amplitude name corresponds to a mass-dependent wave,
    False if it corresponds to a mass-independent (piecewise) wave.
    """
    l_val = extract_amplitude_l(name)
    if l_val is None:
        return True

    waves = getattr(config, "waves", []) or []
    for w in waves:
        if getattr(w, "l", None) == l_val:
            w_type = getattr(w, "type", None)
            w_type_str = str(w_type.value if hasattr(w_type, "value") else w_type).lower()
            if w_type_str in ("piecewise", "mass_independent"):
                return False
            elif w_type_str in ("mass_dependent", "model"):
                return True

    return True


def compute_bin_edges(masses: List[float] | Any, default_step: float = 0.04) -> Any:
    """Compute bin edges for histogram representation given bin centers."""
    import numpy as np
    m_arr = np.array(masses, dtype=float)
    if len(m_arr) == 0:
        return np.array([])
    if len(m_arr) == 1:
        return np.array([m_arr[0] - default_step / 2.0, m_arr[0] + default_step / 2.0])
    inner_edges = (m_arr[:-1] + m_arr[1:]) / 2.0
    first_edge = m_arr[0] - (inner_edges[0] - m_arr[0])
    last_edge = m_arr[-1] + (m_arr[-1] - inner_edges[-1])
    return np.concatenate([[first_edge], inner_edges, [last_edge]])


def extract_credible_bounds(
    masses: List[float],
    mcmc_dict: Dict[float, List[float]],
    nom_vals: List[float],
) -> Tuple[List[float], List[float]]:
    """Extract 16th and 84th percentiles for credible interval band across mass bins."""
    import numpy as np
    lo_68, hi_68 = [], []
    for idx, m in enumerate(masses):
        m_round = round(m, 4)
        samps = mcmc_dict.get(m_round, [])
        if samps:
            lo_68.append(float(np.percentile(samps, 16)))
            hi_68.append(float(np.percentile(samps, 84)))
        else:
            lo_68.append(nom_vals[idx])
            hi_68.append(nom_vals[idx])
    return lo_68, hi_68


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
        # Fallback: export to JSON using ROOT macro ExportAmplitudes.C
        try:
            repo_root = Path(os.environ.get("MOMS2AMPS", str(Path(__file__).resolve().parent.parent.parent)))
            export_macro = repo_root / "macros" / "ExportAmplitudes.C"
            if export_macro.is_file() and Path(nominal_root).is_file():
                cache_dir = out_dir / ".cache"
                cache_dir.mkdir(parents=True, exist_ok=True)
                json_path = cache_dir / "amplitudes_data.json"
                root_bin = get_root_binary(config)
                mcmc_arg = str(Path(mcmc_root).resolve()) if (mcmc_root and Path(mcmc_root).is_file()) else ""
                cmd = [
                    root_bin,
                    "-l",
                    "-b",
                    "-q",
                    f'{export_macro.resolve()}("{Path(nominal_root).resolve()}", "{mcmc_arg}", "{json_path.resolve()}")'
                ]
                env = get_root_subprocess_env(config)
                run_res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
                if run_res.returncode == 0 and json_path.is_file():
                    with open(json_path, "r") as jf:
                        jdata = json.load(jf)
                        raw_masses = jdata.get("masses", [])
                        raw_nominal = jdata.get("nominal", {})
                        if raw_masses and raw_nominal:
                            masses = [float(m) for m in raw_masses]
                            nominal_amps = {k: [float(v) for v in vals] for k, vals in raw_nominal.items()}
                            raw_mcmc = jdata.get("mcmc", {})
                            for k, m_dict in raw_mcmc.items():
                                mcmc_samples[k] = {}
                                for m_str, vals in m_dict.items():
                                    try:
                                        mcmc_samples[k][round(float(m_str), 4)] = [float(v) for v in vals]
                                    except ValueError:
                                        pass
                            loaded = True
        except Exception as e:
            pass

    if not loaded or not masses:
        print(f"[WARNING] Could not read amplitude branches from {nominal_root}")
        return generated

    # Filter active amplitudes: exclude dummy branches that are all zeros and not declared in config.waves
    active_waves_l = [w.l for w in getattr(config, "waves", [])] if getattr(config, "waves", None) else []

    def is_amp_active(k: str) -> bool:
        y_vals = nominal_amps.get(k, [])
        if any(abs(v) > 1e-9 for v in y_vals):
            return True
        l_val = extract_amplitude_l(k)
        if l_val is not None and l_val in active_waves_l:
            return True
        return False

    magnitudes = [k for k in sorted(nominal_amps.keys()) if (k.startswith("a_") or k.startswith("b_")) and is_amp_active(k)]
    phases = [k for k in sorted(nominal_amps.keys()) if (k.startswith("aphi_") or k.startswith("bphi_")) and is_amp_active(k)]

    fit_step = 0.04
    if hasattr(config, "dataset") and getattr(config.dataset, "mass_bins", None):
        fit_step = getattr(config.dataset.mass_bins, "step", 0.04) or 0.04
    elif len(masses) > 1:
        fit_step = float(np.median(np.diff(masses)))

    # --- 1. Amplitude Magnitudes Plot ---
    if magnitudes:
        fig, ax = plt.subplots(figsize=(9, 6), dpi=300)
        colors = plt.cm.tab10(np.linspace(0, 1, len(magnitudes)))

        for idx, k in enumerate(sorted(magnitudes)):
            y_nom = nominal_amps[k]
            label = format_amplitude_label(k)
            is_mass_dep = is_amplitude_mass_dependent(k, config)

            if is_mass_dep:
                # Mass-dependent amplitude -> smooth continuous curve
                m_smooth, y_smooth = smooth_curve(masses, y_nom, num_points=200)
                ax.plot(m_smooth, y_smooth, "-", linewidth=2.0, color=colors[idx], label=label)

                if k in mcmc_samples:
                    lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom)
                    if len(masses) >= 4:
                        m_s, lo_smooth = smooth_curve(masses, lo_68, num_points=200)
                        m_s, hi_smooth = smooth_curve(masses, hi_68, num_points=200)
                        ax.fill_between(m_s, lo_smooth, hi_smooth, color=colors[idx], alpha=0.25)
                    else:
                        ax.fill_between(masses, lo_68, hi_68, color=colors[idx], alpha=0.25)
            else:
                # Mass-independent (piecewise) amplitude -> histogram (stairs)
                edges = compute_bin_edges(masses, default_step=fit_step)
                if hasattr(ax, "stairs"):
                    ax.stairs(y_nom, edges, linewidth=2.0, color=colors[idx], label=label)
                else:
                    ax.step(masses, y_nom, where="mid", linewidth=2.0, color=colors[idx], label=label)

                if k in mcmc_samples:
                    lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom)
                    x_step = np.repeat(edges, 2)[1:-1]
                    lo_step = np.repeat(lo_68, 2)
                    hi_step = np.repeat(hi_68, 2)
                    ax.fill_between(x_step, lo_step, hi_step, color=colors[idx], alpha=0.25)

        ax.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
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
            label = format_amplitude_label(k)
            is_mass_dep = is_amplitude_mass_dependent(k, config)

            if is_mass_dep:
                # Mass-dependent amplitude -> smooth continuous curve
                m_smooth, y_smooth = smooth_curve(masses, y_nom, num_points=200)
                ax.plot(m_smooth, y_smooth, "-", linewidth=2.0, color=colors[idx], label=label)

                if k in mcmc_samples:
                    lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom)
                    if len(masses) >= 4:
                        m_s, lo_smooth = smooth_curve(masses, lo_68, num_points=200)
                        m_s, hi_smooth = smooth_curve(masses, hi_68, num_points=200)
                        ax.fill_between(m_s, lo_smooth, hi_smooth, color=colors[idx], alpha=0.25)
                    else:
                        ax.fill_between(masses, lo_68, hi_68, color=colors[idx], alpha=0.25)
            else:
                # Mass-independent (piecewise) amplitude -> histogram (stairs)
                edges = compute_bin_edges(masses, default_step=fit_step)
                if hasattr(ax, "stairs"):
                    ax.stairs(y_nom, edges, linewidth=2.0, color=colors[idx], label=label)
                else:
                    ax.step(masses, y_nom, where="mid", linewidth=2.0, color=colors[idx], label=label)

                if k in mcmc_samples:
                    lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom)
                    x_step = np.repeat(edges, 2)[1:-1]
                    lo_step = np.repeat(lo_68, 2)
                    hi_step = np.repeat(hi_68, 2)
                    ax.fill_between(x_step, lo_step, hi_step, color=colors[idx], alpha=0.25)

        ax.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
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
        # Resolve target t-bin if applicable
        target_t = getattr(config.dataset, "t_bin", None)
        if not target_t and rundir and rundir.root.name.startswith("t0"):
            target_t = rundir.root.name

        mom_dir_str = config.dataset.get_moments_dir(target_t) if hasattr(config.dataset, "get_moments_dir") else config.dataset.moments_dir
        resolved_dir = Path(mom_dir_str)
        t_name = resolved_dir.name

        search_dirs = [resolved_dir]
        if target_t and (resolved_dir / target_t).is_dir():
            search_dirs.append(resolved_dir / target_t)
        if resolved_dir.parent.is_dir() and resolved_dir.parent != resolved_dir:
            search_dirs.append(resolved_dir.parent)

        cand_names = [
            "UnnormalizedMoments_vs_Mass.root",
            "UnnormalizedMomentsvMass.root",
            "UnnormalizedMoments.root",
            f"UnnormalizedMoments_vs_Mass_{t_name}.root",
            f"UnnormalizedMomentsvMass_{t_name}.root",
            f"UnnormalizedMoments_{t_name}.root",
        ]
        if target_t:
            cand_names.extend([
                f"UnnormalizedMoments_vs_Mass_{target_t}.root",
                f"UnnormalizedMomentsvMass_{target_t}.root",
                f"UnnormalizedMoments_{target_t}.root",
            ])

        for sdir in search_dirs:
            for fname in cand_names:
                cand = sdir / fname
                if cand.is_file():
                    unnorm_path = cand
                    break
            if unnorm_path:
                break
            if sdir.is_dir():
                matches = sorted(sdir.glob("*UnnormalizedMoments*.root"))
                if matches:
                    unnorm_path = matches[0]
                    break

    if not unnorm_path or not unnorm_path.is_file():
        resolved_dir_str = config.dataset.get_moments_dir() if hasattr(config.dataset, "get_moments_dir") else config.dataset.moments_dir
        print(f"[WARNING] Unnormalized moments file not found in {resolved_dir_str}. Skipping moment comparison plots.")
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

    env = get_root_subprocess_env(config)

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
    elif any(m.lower().strip() == "fitted" for m in req_moments):
        selected_moments = [k for k in data_moments.keys() if is_moment_fitted(k, config)]
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

    # Determine fitted mass range from config or model fit points
    fit_mass_min = None
    fit_mass_max = None
    fit_mass_step = 0.04
    if hasattr(config, "dataset") and getattr(config.dataset, "mass_bins", None):
        mb = config.dataset.mass_bins
        fit_mass_min = getattr(mb, "min", None)
        fit_mass_max = getattr(mb, "max", None)
        fit_mass_step = getattr(mb, "step", 0.04) or 0.04

    if fits_list:
        all_masses = []
        for f in fits_list:
            for m_obj in f.get("moments", {}).values():
                all_masses.extend(m_obj.get("mass", []))
        if all_masses:
            fit_masses_unique = sorted(list(set(all_masses)))
            if fit_mass_min is None:
                fit_mass_min = fit_masses_unique[0]
            if fit_mass_max is None:
                fit_mass_max = fit_masses_unique[-1]
            if len(fit_masses_unique) > 1:
                fit_mass_step = float(np.median(np.diff(fit_masses_unique)))

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
                    m_smooth, y_smooth = smooth_curve(mf, yf, num_points=200)
                    if idx == 0:
                        ax.plot(
                            m_smooth, y_smooth, "-", color="red", linewidth=2.0, alpha=0.9,
                            label=f"Rank 1 ($\chi^2 = {fit.get('chi2', 0):.1f}$)", zorder=5
                        )
                    else:
                        lbl = f"Top {len(fits_list)} Lowest $\chi^2$" if idx == 1 else None
                        ax.plot(
                            m_smooth, y_smooth, "-", color="#d62728", linewidth=1.0, alpha=0.35,
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
    # For moments_overview, draw only moments in the fitted range unless explicitly requested otherwise
    if any(m.lower().strip() in ("all", "fitted") for m in req_moments):
        overview_moments = [m for m in selected_moments if is_moment_fitted(m, config)]
        if not overview_moments:
            overview_moments = selected_moments
    else:
        overview_moments = selected_moments

    if len(overview_moments) >= 1:
        K = len(overview_moments)
        ncols = min(4, K)
        nrows = math.ceil(K / ncols)

        fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.2 * nrows), dpi=200)
        axes_flat = np.array(axes).reshape(-1)

        for idx, mom_name in enumerate(overview_moments):
            ax = axes_flat[idx]
            d_val = data_moments[mom_name]
            m_data = d_val.get("mass", [])
            y_data = d_val.get("val", [])
            e_data = d_val.get("err", [])

            # Restrict data points to the fitted mass range
            if fit_mass_min is not None and fit_mass_max is not None:
                tol = fit_mass_step * 0.5
                filtered_pts = [
                    (m, y, e) for m, y, e in zip(m_data, y_data, e_data)
                    if (fit_mass_min - tol) <= m <= (fit_mass_max + tol)
                ]
                if filtered_pts:
                    m_data = [p[0] for p in filtered_pts]
                    y_data = [p[1] for p in filtered_pts]
                    e_data = [p[2] for p in filtered_pts]

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
                        m_smooth, y_smooth = smooth_curve(mf, yf, num_points=200)
                        if f_idx == 0:
                            ax.plot(m_smooth, y_smooth, "-", color="red", linewidth=1.4, alpha=0.9, zorder=5)
                        else:
                            ax.plot(m_smooth, y_smooth, "-", color="#d62728", linewidth=0.8, alpha=0.35, zorder=4)

            # Restrict x-axis limits strictly to the fitted mass range
            if fit_mass_min is not None and fit_mass_max is not None:
                pad = fit_mass_step * 0.6
                ax.set_xlim(fit_mass_min - pad, fit_mass_max + pad)

            ax.set_title(f"${latex_title}$", fontsize=11, fontweight="bold")
            ax.grid(True, linestyle=":", alpha=0.4)
            ax.tick_params(axis="both", which="major", labelsize=8)
            if (idx + ncols >= K) or (idx >= (nrows - 1) * ncols):
                ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=9)

        for idx in range(K, len(axes_flat)):
            axes_flat[idx].set_visible(False)

        fig.suptitle(
            f"{config.project.name} - Unnormalized Moments vs Mass (Top {len(fits_list)} Lowest $\chi^2$ Fits)",
            fontsize=14,
            fontweight="bold",
            y=0.995,
        )
        plt.tight_layout(rect=[0, 0, 1, 0.97])

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


def _plot_single_tbin_worker(
    cfg: FitConfig,
    t_name: str,
    run_dir: Optional[str] = None,
    formats: Tuple[str, ...] = ("png", "pdf"),
) -> Tuple[str, bool, int, str]:
    """
    Worker function to generate diagnostics and publication plots for a single t-bin.
    Returns (t_name, success, num_plots, message).
    """
    try:
        import copy
        import os
        # Set a dedicated MPLCONFIGDIR for each worker to avoid font cache lock contention
        os.environ["MPLCONFIGDIR"] = f"/tmp/mpl_worker_{os.getpid()}"

        from .schema import FitConfig
        from .rundir import get_run_directory

        cfg_t = copy.deepcopy(cfg)
        cfg_t.dataset.t_bin = t_name
        rundir = get_run_directory(cfg_t, run_dir)

        if not rundir.nominal_best_root.is_file():
            legacy_best = rundir.root / "central_best_params.root"
            if legacy_best.is_file():
                rundir.save_nominal_best(legacy_best)
            else:
                return (t_name, False, 0, f"No nominal best fit found at {rundir.nominal_best_root}. Skipping.")

        res = run_all_diagnostics_and_plots(cfg_t, rundir, formats=formats)
        num_plots = len(res.get("plots", []))
        return (t_name, True, num_plots, f"Generated {num_plots} plot files in {rundir.plots_dir}")
    except Exception as e:
        import traceback
        return (t_name, False, 0, f"Error generating plots: {e}\n{traceback.format_exc()}")


def plot_tbins_parallel(
    config: FitConfig,
    tbins: Optional[List[str]] = None,
    run_dir: Optional[str] = None,
    formats: Tuple[str, ...] = ("png", "pdf"),
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Generate diagnostics and publication plots concurrently across all configured t-bins.
    """
    tbins_list = tbins if tbins is not None else config.dataset.get_t_bins()
    if not tbins_list:
        return {"success": False, "total_plots": 0, "tbins": {}}

    half_cpus = get_default_max_workers()
    if max_workers is None:
        cfg_workers = getattr(config.multi_start, "n_tbin_workers", None) if getattr(config, "multi_start", None) else None
        max_workers = min(len(tbins_list), int(cfg_workers)) if cfg_workers else min(len(tbins_list), half_cpus)
    max_workers = max(1, min(max_workers, half_cpus))

    print(f"\n=======================================================")
    print(f" Moments2Amplitudes Publication Plotting Across t-Bins")
    print(f" Project:      {config.project.name}")
    print(f" t-Bins ({len(tbins_list)}): {', '.join(tbins_list)}")
    print(f" Parallel:     {max_workers} worker processes")
    print(f"=======================================================\n")

    results_by_tbin: Dict[str, Dict[str, Any]] = {}
    futures = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for t_name in tbins_list:
            f = executor.submit(
                _plot_single_tbin_worker,
                config,
                t_name,
                run_dir,
                formats,
            )
            futures[f] = t_name

        for f in as_completed(futures):
            t_name = futures[f]
            try:
                t_ret, success, n_plots, msg = f.result()
                results_by_tbin[t_ret] = {
                    "success": success,
                    "num_plots": n_plots,
                    "message": msg,
                }
                status_tag = "[OK]" if success else "[WARNING]"
                print(f"{status_tag} [{t_name}]: {msg}")
            except Exception as e:
                print(f"[ERROR] [{t_name}] Worker crashed: {e}")
                results_by_tbin[t_name] = {
                    "success": False,
                    "num_plots": 0,
                    "message": str(e),
                }

    total_plots = sum(r.get("num_plots", 0) for r in results_by_tbin.values())
    all_success = any(r.get("success", False) for r in results_by_tbin.values())

    print(f"\n=========================================================================================")
    print(f" Parallel Plotting Summary ({len(tbins_list)} t-bins)")
    print(f"=========================================================================================")
    print(f" {'t-bin':<10} {'Status':<16} {'Plots Generated'}")
    print(f" -----------------------------------------------------------------------------------------")
    for t_name in tbins_list:
        r = results_by_tbin.get(t_name, {})
        status_str = "SUCCESS" if r.get("success", False) else "FAILED/SKIPPED"
        n_p = r.get("num_plots", 0)
        print(f" {t_name:<10} {status_str:<16} {n_p}")
    print(f"=========================================================================================\n")
    print(f"[SUCCESS] Generated {total_plots} plot files across {len(tbins_list)} t-bins in parallel.")

    return {
        "success": all_success,
        "total_plots": total_plots,
        "tbins": results_by_tbin,
    }


"""
Stage 4: Plot unnormalized moments vs mass for Linear Algebra and MCMC individually and overlaid.
"""

from __future__ import annotations

import functools
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from .config import MomentFitConfig

print = functools.partial(print, flush=True)


def format_moment_label(moment_name: str) -> str:
    """
    Format moment name (e.g. 'H_0_2_0' or 'H0_2_0' or 'gr_unnormalized_moments_0_2_0_vs_mass')
    into a clean LaTeX string: e.g. Re[H_0(2, 0)], Re[H_1(1, 0)], Im[H_2(1, 1)].
    """
    clean = (
        moment_name.replace("gr_unnormalized_moments_", "")
        .replace("_vs_mass", "")
        .replace("h1_norm_", "")
        .replace("H_", "H")
    )
    parts = clean.split("_")
    if len(parts) >= 3:
        alpha_str = parts[0][1:] if parts[0].startswith("H") else parts[0]
        l_str = parts[1]
        m_str = parts[2]
        if alpha_str in ("0", "1"):
            return rf"$\mathrm{{Re}}[H_{{{alpha_str}}}({l_str}, {m_str})]$"
        elif alpha_str == "2":
            return rf"$\mathrm{{Im}}[H_{{{alpha_str}}}({l_str}, {m_str})]$"
        return rf"$H_{{{alpha_str}}}({l_str}, {m_str})$"
    return moment_name


def parse_moment_tuple(moment_name: str) -> Optional[Tuple[int, int, int]]:
    """Parse (alpha, L, M) integer tuple from a moment name."""
    clean = (
        moment_name.replace("gr_unnormalized_moments_", "")
        .replace("_vs_mass", "")
        .replace("h1_norm_", "")
        .replace("H_", "H")
    )
    parts = clean.split("_")
    if len(parts) >= 3:
        try:
            alpha = int(parts[0][1:] if parts[0].startswith("H") else parts[0])
            l = int(parts[1])
            m = int(parts[2])
            return (alpha, l, m)
        except ValueError:
            return None
    return None


def canonical_moment_name(moment_name: str) -> str:
    """Return canonical 'H_alpha_L_M' name from any variant."""
    tup = parse_moment_tuple(moment_name)
    if tup is not None:
        return f"H_{tup[0]}_{tup[1]}_{tup[2]}"
    return moment_name


def _init_moments_main_env(moments_main_dir: str) -> bool:
    """Set up environment to unpickle Moments.Main MomentResults objects."""
    if moments_main_dir not in sys.path:
        sys.path.insert(0, moments_main_dir)
    try:
        import ROOT
        ROOT.gROOT.SetBatch(True)
        ROOT.gROOT.SetMacroPath(f"{moments_main_dir}:{ROOT.gROOT.GetMacroPath()}")
        import moments.MomentCalculator
        sys.modules["MomentCalculator"] = moments.MomentCalculator
        try:
            from workflow.RootUtilities import loadBasisFunctionsLibrary
            loadBasisFunctionsLibrary()
        except Exception:
            pass
        return True
    except Exception as e:
        print(f"  [Notice] Could not initialize full ROOT Moments.Main environment: {e}")
        return False


def load_linear_algebra_moments(
    pkl_path: Union[str, Path],
    moments_main_dir: str,
) -> Optional[Dict[str, Dict[str, np.ndarray]]]:
    """
    Loads unnormalized physical moments from linear algebra pkl file.

    Returns:
        dict: { 'H_alpha_L_M': {'mass': np.ndarray, 'val': np.ndarray, 'err': np.ndarray} }
        or None if file is missing.
    """
    path = Path(pkl_path)
    if not path.is_file():
        return None

    _init_moments_main_env(moments_main_dir)

    try:
        with open(path, "rb") as f:
            res = pickle.load(f)
    except Exception as e:
        print(f"  [Warning] Failed to unpickle {path}: {e}")
        return None

    moments_data: Dict[str, Dict[str, list]] = {}

    for bin_res in res:
        # Extract mass bin center
        centers = getattr(bin_res, "binCenters", {})
        if not centers:
            continue
        mass = float(list(centers.values())[0])

        indices = getattr(bin_res, "indices", None)
        if indices is None:
            continue

        for qn in indices.qnIndices:
            key = f"H_{qn.momentIndex}_{qn.L}_{qn.M}"
            if key not in moments_data:
                moments_data[key] = {"mass": [], "val": [], "err": []}

            val_obj = bin_res[qn]
            if qn.momentIndex in (0, 1):
                val = float(val_obj.val.real)
                err = float(val_obj.uncertRe)
            else:
                val = float(val_obj.val.imag)
                err = float(val_obj.uncertIm)

            moments_data[key]["mass"].append(mass)
            moments_data[key]["val"].append(val)
            moments_data[key]["err"].append(err)

    # Convert to sorted numpy arrays
    final_data: Dict[str, Dict[str, np.ndarray]] = {}
    for key, d in moments_data.items():
        masses = np.array(d["mass"], dtype=np.float64)
        vals = np.array(d["val"], dtype=np.float64)
        errs = np.array(d["err"], dtype=np.float64)

        sort_idx = np.argsort(masses)
        final_data[key] = {
            "mass": masses[sort_idx],
            "val": vals[sort_idx],
            "err": errs[sort_idx],
        }

    return final_data


def load_mcmc_moments(
    mcmc_tbin_dir: Union[str, Path],
    max_l: int,
) -> Optional[Dict[str, Dict[str, np.ndarray]]]:
    """
    Loads unnormalized physical moments from MCMC outputs.
    Looks for unnormalized_moments_vs_mass*.root or falls back to BruFitResults.MomentsMCMC.

    Returns:
        dict: { 'H_alpha_L_M': {'mass': np.ndarray, 'val': np.ndarray, 'err': np.ndarray} }
        or None if MCMC results are not found.
    """
    tbin_dir = Path(mcmc_tbin_dir)
    if not tbin_dir.is_dir():
        return None

    # Check for pre-generated ROOT files from MonitorMomentFit
    cand_files = [
        tbin_dir / "unnormalized_moments_vs_mass_by_ess.root",
        tbin_dir / "unnormalized_moments_vs_mass.root",
    ]
    root_file = next((f for f in cand_files if f.is_file()), None)

    if root_file is not None:
        try:
            import ROOT
            ROOT.gROOT.SetBatch(True)
            f = ROOT.TFile.Open(str(root_file), "READ")
            if f and not f.IsZombie():
                final_data: Dict[str, Dict[str, np.ndarray]] = {}
                keys = f.GetListOfKeys()
                for key_obj in keys:
                    name = key_obj.GetName()
                    if name.startswith("gr_unnormalized_moments_") and name.endswith("_vs_mass"):
                        tup = parse_moment_tuple(name)
                        if tup is None:
                            continue
                        canon = f"H_{tup[0]}_{tup[1]}_{tup[2]}"
                        gr = f.Get(name)
                        if gr and not gr.IsZombie():
                            n_pts = gr.GetN()
                            masses = np.array([gr.GetPointX(i) for i in range(n_pts)], dtype=np.float64)
                            vals = np.array([gr.GetPointY(i) for i in range(n_pts)], dtype=np.float64)
                            errs = np.array([gr.GetErrorY(i) for i in range(n_pts)], dtype=np.float64)
                            sort_idx = np.argsort(masses)
                            final_data[canon] = {
                                "mass": masses[sort_idx],
                                "val": vals[sort_idx],
                                "err": errs[sort_idx],
                            }
                f.Close()
                if final_data:
                    return final_data
        except Exception as e:
            print(f"  [Notice] Could not read graphs from {root_file}: {e}")

    # Fallback: compute via BruFitResults directly
    try:
        import BruFitResults as BFR
        mcmc = BFR.MomentsMCMC(
            fit_result_dir=str(tbin_dir),
            root_filename="ResultsBruMcmcFixedCovariance.root",
        )
        final_data = {}
        all_moments = list(mcmc.get_moment_name_list())
        if "H_0_0_0" not in all_moments:
            all_moments.append("H_0_0_0")

        for m_name in all_moments:
            try:
                m_masses, m_means, m_errs = mcmc.get_unnormalized_moments_vs_mass(m_name, skip_failed=True)
                if len(m_masses) > 0:
                    canon = canonical_moment_name(m_name)
                    sort_idx = np.argsort(m_masses)
                    final_data[canon] = {
                        "mass": np.array(m_masses[sort_idx], dtype=np.float64),
                        "val": np.array(m_means[sort_idx], dtype=np.float64),
                        "err": np.array(m_errs[sort_idx], dtype=np.float64),
                    }
            except Exception:
                continue

        if final_data:
            return final_data
    except Exception as e:
        print(f"  [Notice] Could not load MCMC results directly via BruFitResults: {e}")

    return None


def plot_single_moment(
    moment_name: str,
    tbin: str,
    max_l: int,
    out_dir: str,
    la_data: Optional[Dict[str, np.ndarray]] = None,
    mcmc_data: Optional[Dict[str, np.ndarray]] = None,
    mode: str = "overlay",  # "linear_algebra", "mcmc", or "overlay"
    formats: Tuple[str, ...] = ("pdf", "png"),
) -> List[Path]:
    """
    Plot unnormalized moment vs mass for linear algebra, MCMC, or overlaid on top of each other.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    saved_paths: List[Path] = []
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=150)

    title_tex = format_moment_label(moment_name)
    has_la = (la_data is not None and len(la_data.get("mass", [])) > 0)
    has_mcmc = (mcmc_data is not None and len(mcmc_data.get("mass", [])) > 0)

    if not has_la and not has_mcmc:
        plt.close(fig)
        return saved_paths

    # Draw Linear Algebra
    if (mode in ("linear_algebra", "overlay")) and has_la:
        ax.errorbar(
            la_data["mass"],
            la_data["val"],
            yerr=la_data["err"],
            fmt="o",
            color="#1f77b4",
            ecolor="#1f77b4",
            elinewidth=1.3,
            capsize=3.0,
            capthick=1.0,
            markersize=5.0,
            label="Linear Algebra (Moments.Main)",
            zorder=3,
        )

    # Draw MCMC
    if (mode in ("mcmc", "overlay")) and has_mcmc:
        ax.errorbar(
            mcmc_data["mass"],
            mcmc_data["val"],
            yerr=mcmc_data["err"],
            fmt="s",
            color="#d62728",
            ecolor="#d62728",
            elinewidth=1.3,
            capsize=3.0,
            capthick=1.0,
            markersize=5.0,
            label="BruFit MCMC (Fixed Covariance)",
            zorder=4,
        )

    # Styling and Labels
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7, zorder=1)
    ax.set_xlabel(r"$m_{\eta\pi^0}\ \left[\mathrm{GeV}/c^2\right]$", fontsize=12)
    ax.set_ylabel("Unnormalized Moment Value", fontsize=12)

    mode_title_map = {
        "linear_algebra": f"{title_tex} — Linear Algebra",
        "mcmc": f"{title_tex} — BruFit MCMC",
        "overlay": f"{title_tex} — Linear Algebra vs MCMC",
    }
    ax.set_title(
        f"{mode_title_map.get(mode, title_tex)} ({tbin}, $L_{{\\max}}={max_l}$)",
        fontsize=13,
        pad=10,
    )

    ax.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=10.5, loc="best")
    ax.grid(True, linestyle=":", alpha=0.6, zorder=0)

    plt.tight_layout()

    prefix_map = {
        "linear_algebra": "moments_linalg",
        "mcmc": "moments_mcmc",
        "overlay": "moments_overlay",
    }
    base_filename = f"{prefix_map.get(mode, 'moments')}_{moment_name}"

    for fmt in formats:
        target = out_path / f"{base_filename}.{fmt}"
        fig.savefig(target, dpi=300 if fmt == "png" else None, bbox_inches="tight")
        saved_paths.append(target)

    plt.close(fig)
    return saved_paths


def plot_overview_grids(
    tbin: str,
    max_l: int,
    out_dir: str,
    la_all: Optional[Dict[str, Dict[str, np.ndarray]]] = None,
    mcmc_all: Optional[Dict[str, Dict[str, np.ndarray]]] = None,
    mode: str = "overlay",
    formats: Tuple[str, ...] = ("pdf", "png"),
) -> List[Path]:
    """
    Generate multi-panel overview grids for H0, H1, and H2 multipoles.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    saved_paths: List[Path] = []
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Gather all available moments
    all_keys = set()
    if la_all:
        all_keys.update(la_all.keys())
    if mcmc_all:
        all_keys.update(mcmc_all.keys())

    if not all_keys:
        return saved_paths

    # Group by alpha (0, 1, 2)
    alpha_groups: Dict[int, List[str]] = {0: [], 1: [], 2: []}
    for k in sorted(all_keys):
        tup = parse_moment_tuple(k)
        if tup is not None and tup[0] in alpha_groups:
            alpha_groups[tup[0]].append(k)

    prefix_map = {
        "linear_algebra": "overview_linalg_grid",
        "mcmc": "overview_mcmc_grid",
        "overlay": "overview_overlay_grid",
    }

    for alpha, moments_list in alpha_groups.items():
        if not moments_list:
            continue

        n_moments = len(moments_list)
        n_cols = 3
        n_rows = int(np.ceil(n_moments / n_cols))

        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(5.0 * n_cols, 3.6 * n_rows),
            dpi=150,
            sharex=True,
        )
        axes_flat = axes.flatten() if n_moments > 1 else [axes]

        alpha_name = {0: "H_0", 1: "H_1", 2: "H_2"}.get(alpha, f"H_{alpha}")

        for idx, m_key in enumerate(moments_list):
            ax = axes_flat[idx]
            la_data = la_all.get(m_key) if la_all else None
            mcmc_data = mcmc_all.get(m_key) if mcmc_all else None

            has_la = (la_data is not None and len(la_data.get("mass", [])) > 0)
            has_mcmc = (mcmc_data is not None and len(mcmc_data.get("mass", [])) > 0)

            if (mode in ("linear_algebra", "overlay")) and has_la:
                ax.errorbar(
                    la_data["mass"], la_data["val"], yerr=la_data["err"],
                    fmt="o", color="#1f77b4", ecolor="#1f77b4", elinewidth=1.0,
                    capsize=2.0, markersize=3.5, label="Linear Alg.", zorder=3,
                )

            if (mode in ("mcmc", "overlay")) and has_mcmc:
                ax.errorbar(
                    mcmc_data["mass"], mcmc_data["val"], yerr=mcmc_data["err"],
                    fmt="s", color="#d62728", ecolor="#d62728", elinewidth=1.0,
                    capsize=2.0, markersize=3.5, label="BruFit MCMC", zorder=4,
                )

            ax.axhline(0, color="gray", linestyle="--", linewidth=0.7, alpha=0.7)
            ax.set_title(format_moment_label(m_key), fontsize=11, pad=4)
            ax.grid(True, linestyle=":", alpha=0.5)

            if idx >= (n_rows - 1) * n_cols:
                ax.set_xlabel(r"$m_{\eta\pi^0}\ [\mathrm{GeV}/c^2]$", fontsize=10)
            if idx % n_cols == 0:
                ax.set_ylabel("Unnorm. Moment", fontsize=10)

        # Hide unused subplots
        for idx in range(n_moments, len(axes_flat)):
            axes_flat[idx].set_visible(False)

        # Add single legend to first subplot or top
        handles, labels = axes_flat[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.98, 0.99), frameon=True, fontsize=10)

        fig.suptitle(
            f"Unnormalized ${alpha_name}$ Multipoles ({mode.replace('_', ' ').title()}) — {tbin}, $L_{{\\max}}={max_l}$",
            fontsize=14,
            y=0.995,
        )
        plt.tight_layout(rect=[0, 0, 1, 0.98])

        base_name = f"{prefix_map.get(mode, 'overview')}_{alpha_name}"
        for fmt in formats:
            target = out_path / f"{base_name}.{fmt}"
            fig.savefig(target, dpi=250 if fmt == "png" else None, bbox_inches="tight")
            saved_paths.append(target)

        plt.close(fig)

    return saved_paths


def run_stage4(
    cfg: MomentFitConfig,
    tbins_filter: Optional[List[str]] = None,
    max_l_filter: Optional[List[int]] = None,
    plot_type: Optional[str] = None,
    moments_filter: Optional[List[str]] = None,
    dry_run: bool = False,
) -> bool:
    """
    Executes Stage 4: Plot unnormalized moments vs mass for Linear Algebra and MCMC.
    Supports individual plots and overlay comparison plots.
    """
    tbins = tbins_filter if tbins_filter is not None else cfg.analysis.tbins
    max_ls = max_l_filter if max_l_filter is not None else cfg.analysis.max_l
    p_cfg = cfg.plotting

    # Determine plot modes to execute
    # CLI override has highest precedence; otherwise use YAML flags
    modes_to_run: List[str] = []
    if plot_type:
        pt = plot_type.lower().strip()
        if pt in ("all", "both"):
            modes_to_run = ["linear_algebra", "mcmc", "overlay"]
        elif pt in ("la", "linalg", "linear_algebra"):
            modes_to_run = ["linear_algebra"]
        elif pt in ("mcmc", "brufit"):
            modes_to_run = ["mcmc"]
        elif pt in ("overlay", "comparison", "compare"):
            modes_to_run = ["overlay"]
        else:
            print(f"❌ Unknown plot-type '{plot_type}'. Expected: all, linear_algebra, mcmc, or overlay.")
            return False
    else:
        if p_cfg.plot_linear_algebra:
            modes_to_run.append("linear_algebra")
        if p_cfg.plot_mcmc:
            modes_to_run.append("mcmc")
        if p_cfg.plot_overlay:
            modes_to_run.append("overlay")

    if not modes_to_run:
        print("ℹ️  No plot types selected in configuration or CLI. Skipping plotting.")
        return True

    # Filter requested moments
    active_moments = moments_filter if moments_filter is not None else p_cfg.moments
    if not active_moments:
        active_moments = ["all"]
    filter_all = any(m.lower().strip() == "all" for m in active_moments)

    print("\n" + "=" * 70)
    print(" STAGE 4: Plot Unnormalized Moments vs Mass")
    print("=" * 70)
    print(f" Target t-bins:      {', '.join(tbins)}")
    print(f" Multipoles Max L:   {max_ls}")
    print(f" Plot Modes:         {', '.join(modes_to_run)}")
    print(f" Moments Filter:     {'All moments' if filter_all else ', '.join(active_moments)}")
    print(f" Overview Grids:     {p_cfg.grid}")
    print(f" Output Base Dir:    {p_cfg.output_dir}")
    print(f" Output Formats:     {p_cfg.formats}")
    print("=" * 70 + "\n")

    if dry_run:
        print("Dry-run requested: skipping file generation.")
        return True

    total_plots_generated = 0

    for tbin in tbins:
        for max_l in max_ls:
            print(f"\n--- Processing Plots for {tbin} (maxL = {max_l}) ---")
            sub_dir = p_cfg.get_output_dir(tbin, max_l)

            # 1. Load Linear Algebra moments
            la_out_dir = cfg.linear_algebra.get_output_dir(tbin, max_l, cfg.analysis.beam_pol_label)
            phys_pkl = os.path.join(la_out_dir, "unnorm_moments_phys.pkl")
            la_data = load_linear_algebra_moments(phys_pkl, cfg.environment.moments_main_dir)
            if la_data is not None:
                print(f"  ✅ Loaded Linear Algebra moments ({len(la_data)} moments) from {phys_pkl}")
            else:
                print(f"  ⚠️  Linear Algebra moments pkl not found: {phys_pkl}")

            # 2. Load MCMC moments
            m_min = cfg.analysis.mass_binning.min
            m_max = cfg.analysis.mass_binning.max
            mcmc_run_dir = cfg.moment_fit.get_run_dir(max_l, m_min, m_max)
            mcmc_tbin_dir = os.path.join(mcmc_run_dir, tbin)
            mcmc_data = load_mcmc_moments(mcmc_tbin_dir, max_l)
            if mcmc_data is not None:
                print(f"  ✅ Loaded BruFit MCMC moments ({len(mcmc_data)} moments) from {mcmc_tbin_dir}")
            else:
                print(f"  ⚠️  BruFit MCMC results not found in: {mcmc_tbin_dir}")

            # Determine union of moments to plot
            avail_moments = set()
            if la_data:
                avail_moments.update(la_data.keys())
            if mcmc_data:
                avail_moments.update(mcmc_data.keys())

            if not avail_moments:
                print(f"  ❌ No moment data available for {tbin} maxL={max_l}. Skipping.")
                continue

            target_moments = sorted(
                avail_moments if filter_all else [m for m in avail_moments if m in active_moments]
            )

            print(f"  Generating plots for {len(target_moments)} moments...")

            for mode in modes_to_run:
                # Check data availability for mode
                if mode == "linear_algebra" and not la_data:
                    print(f"    [Skip] Linear algebra individual plots: data unavailable.")
                    continue
                if mode == "mcmc" and not mcmc_data:
                    print(f"    [Skip] MCMC individual plots: data unavailable.")
                    continue
                if mode == "overlay" and (not la_data or not mcmc_data):
                    print(f"    [Notice] Overlay requires both Linear Algebra and MCMC data. "
                          f"(LA={la_data is not None}, MCMC={mcmc_data is not None}). Skipping overlay.")
                    continue

                mode_dir = os.path.join(sub_dir, mode)
                mode_count = 0

                for m_key in target_moments:
                    la_moment = la_data.get(m_key) if la_data else None
                    mcmc_moment = mcmc_data.get(m_key) if mcmc_data else None

                    saved = plot_single_moment(
                        moment_name=m_key,
                        tbin=tbin,
                        max_l=max_l,
                        out_dir=mode_dir,
                        la_data=la_moment,
                        mcmc_data=mcmc_moment,
                        mode=mode,
                        formats=tuple(p_cfg.formats),
                    )
                    mode_count += len(saved)

                # Generate overview grids if requested
                if p_cfg.grid:
                    grid_saved = plot_overview_grids(
                        tbin=tbin,
                        max_l=max_l,
                        out_dir=mode_dir,
                        la_all=la_data,
                        mcmc_all=mcmc_data,
                        mode=mode,
                        formats=tuple(p_cfg.formats),
                    )
                    mode_count += len(grid_saved)

                total_plots_generated += mode_count
                print(f"    [OK] Mode '{mode}': generated {mode_count} plot files in {mode_dir}")

    print("\n" + "=" * 70)
    print(f"✅ Stage 4 complete! Total plot files generated: {total_plots_generated}")
    print("=" * 70 + "\n")
    return True

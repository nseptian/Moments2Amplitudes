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
import re
import json
import copy
import subprocess
import shutil
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

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
from .cross_section import dynamic_breit_wigner_width


# Curated colorblind-safe publication palette (Tol / Okabe-Ito high-contrast)
PHYSICS_PALETTE: List[str] = [
    "#0072B2",  # Deep Navy / Royal Blue
    "#D55E00",  # Vermilion / Rust
    "#009E73",  # Deep Teal / Pine Green
    "#785EF0",  # Purple / Indigo
    "#E69F00",  # Warm Amber / Ochre
    "#CC79A7",  # Muted Rose
    "#56B4E9",  # Slate / Sky Blue
    "#4B5563",  # Charcoal Slate
    "#882255",  # Wine / Claret
    "#117733",  # Forest Green
]

# Standard linestyles for monochromatic physics publication plots
PHYSICS_LINESTYLES: List[Any] = [
    "--",
    "-.",
    ":",
    (0, (3, 1, 1, 1)),  # densely dashdotdotted
    (0, (5, 2)),        # loosely dashed
    (0, (1, 1)),        # densely dotted
    (0, (3, 5, 1, 5)),  # dashdot
    (0, (3, 1, 3, 1)),  # dashdash
]

# Phase axis limits and standard radian tick marks [-pi, +pi]
PHASE_Y_LIMITS: Tuple[float, float] = (-math.pi, math.pi)
PHASE_Y_TICKS: List[float] = [-math.pi, -math.pi / 2.0, 0.0, math.pi / 2.0, math.pi]
PHASE_Y_TICKLABELS: List[str] = [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$+\pi/2$", r"$+\pi$"]


def configure_phase_axis(
    ax: Any,
    limits: Tuple[float, float] = PHASE_Y_LIMITS,
    ticks: Optional[List[float]] = None,
    ticklabels: Optional[List[str]] = None,
) -> None:
    """
    Configure a phase subplot axis with range [-pi, +pi]
    and clean radian tick marks [-pi, -pi/2, 0, +pi/2, +pi].
    """
    import matplotlib.ticker as ticker
    ax.set_ylim(limits)
    if ticks is None:
        ticks = PHASE_Y_TICKS
    if ticklabels is None:
        ticklabels = PHASE_Y_TICKLABELS
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticklabels)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

PHYSICS_RCPARAMS: Dict[str, Any] = {
    # Inward ticks on all 4 sides with minor ticks
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size": 5.0,
    "ytick.major.size": 5.0,
    "xtick.minor.size": 2.5,
    "ytick.minor.size": 2.5,
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    # Axes lines and fonts
    "axes.linewidth": 1.0,
    "axes.labelsize": 11.5,
    "axes.titlesize": 12.0,
    "axes.labelpad": 4.0,
    "axes.titlepad": 6.0,
    "xtick.labelsize": 10.0,
    "ytick.labelsize": 10.0,
    # Fonts
    "font.size": 11.0,
    "font.family": "sans-serif",
    "mathtext.fontset": "dejavusans",
    # Grid
    "axes.grid": False,
    "grid.linestyle": ":",
    "grid.alpha": 0.45,
    "grid.color": "#B0B0B0",
    # Legend
    "legend.fontsize": 9.0,
    "legend.frameon": True,
    "legend.framealpha": 0.9,
    "legend.edgecolor": "#CCCCCC",
    "legend.fancybox": False,
    # Saving
    "figure.autolayout": False,
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
}


def apply_plot_style(style: Optional[str] = "default") -> None:
    """
    Apply high-energy physics plotting style using mplhep, or clean publication physics style.
    Defaults to 'default' (publication-quality physics style).
    Supports LHCb, LHCb1, LHCb2, ROOT, ATLAS, CMS, ALICE, DUNE, or 'default'.
    Falls back gracefully if mplhep or style is unavailable.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chosen = (style or "default").strip()
    if chosen.lower() in ("default", "none", "matplotlib", "standard", ""):
        plt.style.use("default")
        matplotlib.rcParams.update(PHYSICS_RCPARAMS)
    else:
        try:
            import mplhep as hep
            style_map = {
                "lhcb": hep.style.LHCb,
                "lhcb1": hep.style.LHCb1,
                "lhcb2": hep.style.LHCb2,
                "root": hep.style.ROOT,
                "atlas": hep.style.ATLAS,
                "atlas1": hep.style.ATLAS1,
                "atlas2": hep.style.ATLAS2,
                "cms": hep.style.CMS,
                "alice": hep.style.ALICE,
                "dune": hep.style.DUNE,
            }
            if chosen.lower() in style_map:
                hep.style.use(style_map[chosen.lower()])
            else:
                hep.style.use(chosen)
        except Exception:
            try:
                plt.style.use(chosen)
            except Exception:
                plt.style.use("default")
                matplotlib.rcParams.update(PHYSICS_RCPARAMS)

    # Ensure background grid is suppressed by default across all styles
    matplotlib.rcParams["axes.grid"] = False
    plt.rcParams["axes.grid"] = False


def should_show_grid(config: Optional[Any] = None) -> bool:
    """Return True only if show_grid is explicitly enabled in plotting configuration."""
    if config and getattr(config, "plotting", None):
        return bool(getattr(config.plotting, "show_grid", False))
    return False


def discover_lowest_chi2_fits(rundir: RunDirectory, n_lowest: int = 10) -> List[Tuple[float, Path]]:
    """
    Discover completed fit runs in fit_starts_dir, sort by chi2 ascending,
    and return the top n_lowest (chi2, root_file) tuples.
    Only includes strictly converged fits (status == 0 and is_valid == True).
    """
    valid_fits: List[Tuple[float, Path]] = []
    has_starts_logs = False

    if rundir.fit_starts_dir.is_dir():
        for log_f in rundir.fit_starts_dir.glob("fit_seed_*.log"):
            has_starts_logs = True
            try:
                with open(log_f, "r") as f:
                    for line in f:
                        if "[RESULT]" in line:
                            tokens = dict(part.split("=") for part in line.strip().split() if "=" in part)
                            if "CHI2" in tokens and "VALID" in tokens and "STATUS" in tokens:
                                is_valid = tokens["VALID"] == "1"
                                at_limit = tokens.get("AT_LIMIT") == "1"
                                status = int(tokens["STATUS"])
                                chi2_val = float(tokens["CHI2"])
                                seed_val = int(tokens.get("SEED", 0))

                                if is_valid and not at_limit and status == 0 and chi2_val < 1e8:
                                    cand1 = rundir.fit_starts_dir / f"MassDependentResults_{seed_val}.root"
                                    cand2 = rundir.get_seed_root_file(seed_val)
                                    root_cand = cand1 if cand1.is_file() else (cand2 if cand2.is_file() else None)
                                    if root_cand:
                                        valid_fits.append((chi2_val, root_cand))
            except Exception:
                pass

    valid_fits.sort(key=lambda x: x[0])

    # Ensure nominal_best is included as the primary best fit if available and valid
    if rundir.nominal_best_root.is_file():
        if valid_fits:
            # Replace the best seed file with nominal_best_root to avoid duplicate entries for rank 1
            valid_fits[0] = (valid_fits[0][0], rundir.nominal_best_root)
        elif not has_starts_logs:
            # Only include standalone nominal_best if no failure recorded in summary
            is_nominal_converged = True
            if rundir.summary_txt.is_file():
                try:
                    summary_content = rundir.summary_txt.read_text()
                    if "Minimizer Status:" in summary_content and "Minimizer Status:   0" not in summary_content:
                        is_nominal_converged = False
                    if "Valid Minimum:      No" in summary_content:
                        is_nominal_converged = False
                    if "NOT CONVERGED" in summary_content or "FAILED (" in summary_content or "AT LIMIT" in summary_content:
                        is_nominal_converged = False
                except Exception:
                    pass
            if is_nominal_converged:
                valid_fits.append((0.0, rundir.nominal_best_root))

    return valid_fits[:n_lowest]


def _get_hist_bins(
    n_samples: int,
    user_bins: Optional[int] = None,
    is_overview: bool = False,
) -> int:
    """Compute finer adaptive histogram bins for posterior distributions."""
    if user_bins is not None and user_bins > 0:
        nb = user_bins
    elif is_overview:
        nb = min(60, max(35, int(np.sqrt(n_samples) * 2.8)))
    else:
        nb = min(100, max(50, int(np.sqrt(n_samples) * 4.0)))
    if n_samples < nb:
        nb = max(10, n_samples // 2) if n_samples > 20 else max(2, n_samples)
    return max(2, nb)


def _get_chi2_landscape_bins(
    n_samples: int,
    user_bins: Optional[int] = None,
) -> int:
    """Compute finer adaptive histogram bins for multi-start Chi2 landscape plots."""
    if user_bins is not None and user_bins > 0:
        nb = user_bins
    else:
        nb = min(120, max(60, int(np.sqrt(n_samples) * 3.5)))
    if n_samples < nb:
        nb = max(10, n_samples // 2) if n_samples > 20 else max(2, n_samples)
    return max(2, nb)


def plot_multistart_landscape(
    results: List[Dict[str, Any]],
    output_dir: Path | str,
    project_name: str = "Moments2Amplitudes",
    plateau_delta: float = 10.0,
    formats: Tuple[str, ...] = ("png", "pdf"),
    style: Optional[str] = "default",
    show_grid: bool = False,
    bins: Optional[int] = None,
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

    apply_plot_style(style)

    chi2s = [r["chi2"] for r in results if r.get("is_valid", False) and r.get("status", -1) == 0 and r.get("chi2", 1e9) < 1e6]
    if not chi2s:
        return generated

    min_chi2 = min(chi2s)
    plateau_thresh = min_chi2 + plateau_delta

    fig, ax = plt.subplots(
        figsize=(7.5, 4.8), dpi=300,
        gridspec_kw={"top": 0.92, "bottom": 0.12, "left": 0.11, "right": 0.96}
    )
    n_hist_bins = _get_chi2_landscape_bins(len(chi2s), user_bins=bins)
    counts, hist_bins, patches = ax.hist(
        chi2s, bins=n_hist_bins, color="#B0C4DE", edgecolor="#1A2E40", alpha=0.75,
        label=f"Converged Fits (N={len(chi2s)})"
    )

    ax.axvline(min_chi2, color="#D55E00", linestyle="--", linewidth=1.8, label=f"Global Min $\\chi^2 = {min_chi2:.2f}$")
    ax.axvspan(min_chi2, plateau_thresh, color="#009E73", alpha=0.18, label=f"Plateau $\\leq +{plateau_delta:.1f}$")

    ax.set_xlabel(r"Total $\chi^2$", fontsize=11.5)
    ax.set_ylabel("Number of Starts", fontsize=11.5)
    ax.set_title(f"{project_name} - Multi-Start $\\chi^2$ Landscape (Converged Fits)", fontsize=12.5, pad=8)
    ax.legend(frameon=True, fontsize=9.5)
    if show_grid:
        ax.grid(True, linestyle=":", alpha=0.45)
    else:
        ax.grid(False)

    for fmt in formats:
        p = out_dir / f"multistart_chi2_landscape.{fmt}"
        fig.savefig(p, bbox_inches="tight")
        generated.append(p)

    plt.close(fig)
    return generated


def _format_t_val(v: float) -> str:
    """Format t value preserving 3 decimals for fractional bins (e.g. 0.325) and 2 decimals for others."""
    s3 = f"{v:.3f}"
    if s3.endswith("5") or (not s3.endswith("0")):
        return s3
    return f"{v:.2f}"


def format_tbin_label(tbin: str) -> str:
    """Format t-bin name like 't010020' or 't020032' into '$0.10 < -t < 0.20\\text{ GeV}^2$'."""
    from .cross_section import parse_t_bin_from_string
    parsed = parse_t_bin_from_string(tbin)
    if parsed:
        low, high = parsed
        return f"${_format_t_val(low)} < -t < {_format_t_val(high)}\\text{{ GeV}}^2$"
    return tbin


def plot_multistart_landscape_all_tbins(
    results_by_tbin: Dict[str, Any],
    output_dir: Path | str,
    project_name: str = "Moments2Amplitudes",
    plateau_delta: float = 10.0,
    formats: Tuple[str, ...] = ("png", "pdf"),
    style: Optional[str] = "default",
    show_grid: bool = False,
    bins: Optional[int] = None,
) -> List[Path]:
    """
    Plot combined multi-start Chi2 distribution landscapes across all t-bins side-by-side in one row (1 x N_tbins).
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
        print("[WARNING] matplotlib not found. Skipping multi-start combined plot.")
        return generated

    def _parse_logs_for_chi2s(starts_dir: Path) -> List[float]:
        vals: List[float] = []
        if starts_dir.is_dir():
            for log_f in sorted(starts_dir.glob("fit_seed_*.log")):
                try:
                    for line in log_f.read_text().splitlines():
                        if "[RESULT]" in line:
                            toks = dict(p.split("=") for p in line.strip().split() if "=" in p)
                            if toks.get("VALID") == "1" and toks.get("STATUS") == "0":
                                val = float(toks.get("CHI2", 1e9))
                                if val < 1e6:
                                    vals.append(val)
                except Exception:
                    pass
        return vals

    # Filter t-bins that have valid results
    tbins_data: Dict[str, List[float]] = {}
    for t_name, t_info in results_by_tbin.items():
        chi2s: List[float] = []
        if isinstance(t_info, dict):
            if "chi2s" in t_info and t_info["chi2s"]:
                chi2s = [float(v) for v in t_info["chi2s"] if float(v) < 1e6]
            elif "chi2_list" in t_info and t_info["chi2_list"]:
                chi2s = [float(v) for v in t_info["chi2_list"] if float(v) < 1e6]
            elif "results" in t_info and t_info["results"]:
                chi2s = [
                    r["chi2"] for r in t_info["results"]
                    if isinstance(r, dict) and r.get("is_valid", False) and r.get("status", -1) == 0 and r.get("chi2", 1e9) < 1e6
                ]
            elif "rundir" in t_info and t_info["rundir"]:
                rundir_obj = t_info["rundir"]
                if hasattr(rundir_obj, "fit_starts_dir") and rundir_obj.fit_starts_dir.is_dir():
                    chi2s = _parse_logs_for_chi2s(rundir_obj.fit_starts_dir)
                else:
                    p = Path(rundir_obj)
                    starts_dir = p / "fit" / "starts" if (p / "fit" / "starts").is_dir() else (p if p.is_dir() else None)
                    if starts_dir:
                        chi2s = _parse_logs_for_chi2s(starts_dir)
        elif isinstance(t_info, list):
            if t_info and isinstance(t_info[0], (int, float)):
                chi2s = [float(v) for v in t_info if float(v) < 1e6]
            else:
                chi2s = [
                    r["chi2"] for r in t_info
                    if isinstance(r, dict) and r.get("is_valid", False) and r.get("status", -1) == 0 and r.get("chi2", 1e9) < 1e6
                ]
        elif hasattr(t_info, "fit_starts_dir"):
            chi2s = _parse_logs_for_chi2s(t_info.fit_starts_dir)
        elif isinstance(t_info, (str, Path)):
            p = Path(t_info)
            starts_dir = p / "fit" / "starts" if (p / "fit" / "starts").is_dir() else (p if p.is_dir() else None)
            if starts_dir:
                chi2s = _parse_logs_for_chi2s(starts_dir)

        if chi2s:
            tbins_data[t_name] = chi2s

    if not tbins_data:
        return generated

    # Deduplicate t-bin keys canonicalizing by t_name regex (e.g. t010020)
    canonical_data: Dict[str, List[float]] = {}
    for t_name, chi2_list in tbins_data.items():
        m = re.search(r"t\d{6}", t_name)
        canon_name = m.group(0) if m else t_name
        if canon_name not in canonical_data or len(chi2_list) > len(canonical_data[canon_name]):
            canonical_data[canon_name] = chi2_list
    tbins_data = canonical_data

    def _tbin_sort_key(name: str):
        m = re.search(r"(\d+)", name)
        return int(m.group(1)) if m else name

    sorted_tbins = sorted(tbins_data.keys(), key=_tbin_sort_key)
    tbins_data = {k: tbins_data[k] for k in sorted_tbins}

    apply_plot_style(style)

    n_bins = len(tbins_data)
    nrows = 1
    ncols = n_bins
    figsize = (max(4.0 * n_bins, 7.5), 4.6)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, dpi=300)
    axes_flat = [axes] if n_bins == 1 else (axes.flatten() if hasattr(axes, "flatten") else [axes])

    for idx, (t_name, chi2_list) in enumerate(tbins_data.items()):
        ax = axes_flat[idx]
        min_chi2 = min(chi2_list)
        plateau_thresh = min_chi2 + plateau_delta
        n_plateau = sum(1 for c in chi2_list if c <= plateau_thresh)

        n_hist_bins = _get_chi2_landscape_bins(len(chi2_list), user_bins=bins)
        ax.hist(
            chi2_list, bins=n_hist_bins, color="#B0C4DE", edgecolor="#1A2E40", alpha=0.75,
            label=f"Converged (N={len(chi2_list)})"
        )
        ax.axvline(
            min_chi2, color="#D55E00", linestyle="--", linewidth=1.8,
            label=f"Min $\\chi^2 = {min_chi2:.2f}$"
        )
        ax.axvspan(
            min_chi2, plateau_thresh, color="#009E73", alpha=0.18,
            label=f"Plateau $\\leq +{plateau_delta:.1f}$ ($N={n_plateau}$)"
        )

        t_lbl = format_tbin_label(t_name)
        ax.set_title(t_lbl, fontsize=11.5, pad=6)
        ax.set_xlabel(r"Total $\chi^2$", fontsize=10.5)
        ax.set_ylabel("Number of Starts", fontsize=10.5)
        ax.legend(frameon=True, fontsize=8.5, loc="upper right")
        if show_grid:
            ax.grid(True, linestyle=":", alpha=0.5)
        else:
            ax.grid(False)

    fig.suptitle(f"{project_name} - Multi-Start $\\chi^2$ Landscape Across $t$-Bins", fontsize=13.5, y=1.02)
    fig.tight_layout()

    for fmt in formats:
        p = out_dir / f"multistart_chi2_landscape_all_tbins.{fmt}"
        fig.savefig(p, bbox_inches="tight")
        generated.append(p)

    plt.close(fig)
    return generated


L_SPECTROSCOPIC: Dict[int, str] = {
    0: "S",
    1: "P",
    2: "D",
    3: "F",
    4: "G",
    5: "H",
}


def get_amplitude_symbol_parts(name: str) -> Optional[Tuple[str, str, str, str]]:
    """
    Parse amplitude branch name (e.g. 'a_2_-1', 'b_0_0', 'aphi_1_1') or spectroscopic string
    (e.g. 'D+_2', 'D-_1', '$D^{+}_{2}$') into:
    (l_symbol, refl, m_str, wave_kind)
    where:
      l_symbol: 'S', 'P', 'D', 'F', ...
      refl: '+' for 'a'/'aphi', '-' for 'b'/'bphi'
      m_str: '-1', '0', '1', etc.
      wave_kind: 'amp' or 'phase'
    """
    if not name:
        return None
    parts = name.split("_")
    if len(parts) >= 3:
        prefix = parts[0]
        l_str = parts[1]
        m_str = parts[2]
        if prefix in ("a", "b", "aphi", "bphi"):
            try:
                l_val = int(l_str)
                l_sym = L_SPECTROSCOPIC.get(l_val, f"L_{{{l_val}}}")
                int(m_str)
            except ValueError:
                l_sym = l_str
            refl = "+" if prefix in ("a", "aphi") else "-"
            kind = "phase" if prefix in ("aphi", "bphi") else "amp"
            return (l_sym, refl, m_str, kind)

    # Check spectroscopic format: e.g. D+_2, D-_1, P+_1, D^+_2, $D^{+}_{2}$
    import re
    m = re.match(r"^\$?([SPDFGHspdfgh])(?:\^\{?|\s*)([\+\-])(?:\}?_?\{?|_+)(-?\d+)\}?\$?$", name.strip())
    if m:
        l_sym = m.group(1).upper()
        refl = m.group(2)
        m_str = m.group(3)
        return (l_sym, refl, m_str, "amp")

    return None


def format_amplitude_symbol(name: str) -> str:
    """Format amplitude name into bare LaTeX mathtext (e.g. a_2_-1 -> D^{+}_{-1}, a_0_0 -> S^{+}_{0})."""
    parsed = get_amplitude_symbol_parts(name)
    if parsed:
        l_sym, refl, m_str, kind = parsed
        if kind == "phase":
            return f"\\phi({l_sym}^{{{refl}}}_{{{m_str}}})"
        return f"{l_sym}^{{{refl}}}_{{{m_str}}}"
    clean_name = name.replace("_", r"\_")
    return clean_name


def format_amplitude_label(name: str, plus_minus: bool = False) -> str:
    """Format amplitude name into LaTeX math string (e.g. a_2_-1 -> $D^{+}_{-1}$, aphi_2_1 -> $\\phi(D^{+}_{1})$)."""
    parsed = get_amplitude_symbol_parts(name)
    if parsed:
        l_sym, refl, m_str, kind = parsed
        if kind == "phase":
            prefix = "\\pm\\phi" if plus_minus else "\\phi"
            return f"${prefix}({l_sym}^{{{refl}}}_{{{m_str}}})$"
        return f"${l_sym}^{{{refl}}}_{{{m_str}}}$"
    clean_name = name.replace("_", r"\_")
    return f"${clean_name}$"


def format_amplitude_squared_symbol(name: str) -> str:
    """Format amplitude squared name into bare LaTeX mathtext (e.g. a_2_-1 -> |D^{+}_{-1}|^2)."""
    parsed = get_amplitude_symbol_parts(name)
    if parsed:
        l_sym, refl, m_str, _ = parsed
        return f"|{l_sym}^{{{refl}}}_{{{m_str}}}|^2"
    clean_name = name.replace("_", r"\_")
    return f"|{clean_name}|^2"


def format_amplitude_squared_label(name: str) -> str:
    """Format amplitude squared name into LaTeX math string (e.g. a_2_-1 -> $|D^{+}_{-1}|^2$)."""
    return f"${format_amplitude_squared_symbol(name)}$"


def wrap_angle_pi(angles: Any) -> Any:
    """Wrap angle(s) into [-pi, pi] interval."""
    arr = np.asarray(angles, dtype=float)
    return (arr + np.pi) % (2.0 * np.pi) - np.pi


def split_phase_jumps(m_arr: Any, y_arr: Any, jump_threshold: float = 3.0) -> Tuple[np.ndarray, np.ndarray]:
    """Insert np.nan where adjacent points have phase jumps > jump_threshold (e.g. branch cuts across +/- pi)."""
    m_arr = np.asarray(m_arr, dtype=float)
    y_arr = np.asarray(y_arr, dtype=float)
    if len(y_arr) < 2:
        return m_arr, y_arr
    jumps = np.where(np.abs(np.diff(y_arr)) > jump_threshold)[0]
    if len(jumps) == 0:
        return m_arr, y_arr
    m_out = list(m_arr)
    y_out = list(y_arr)
    for idx in reversed(jumps):
        m_out.insert(idx + 1, np.nan)
        y_out.insert(idx + 1, np.nan)
    return np.array(m_out), np.array(y_out)


def format_relative_phase_symbol(name: str, anchor_name: Optional[str] = None) -> str:
    """Format relative phase bare symbol (e.g. aphi_1_1 vs aphi_2_1 -> \\Delta\\phi(P^{+}_{1} - D^{+}_{1}))."""
    p_name = get_amplitude_symbol_parts(name)
    p_anch = get_amplitude_symbol_parts(anchor_name) if anchor_name else None
    if p_name and p_anch:
        s_name = f"{p_name[0]}^{{{p_name[1]}}}_{{{p_name[2]}}}"
        s_anch = f"{p_anch[0]}^{{{p_anch[1]}}}_{{{p_anch[2]}}}"
        return f"\\Delta\\phi({s_name} - {s_anch})"
    elif p_name:
        s_name = f"{p_name[0]}^{{{p_name[1]}}}_{{{p_name[2]}}}"
        if anchor_name:
            anch_str = anchor_name.replace("_", r"\_")
            if "1320" in anchor_name:
                anch_str = r"a_{2}(1320)"
            elif "1700" in anchor_name:
                anch_str = r"a_{2}(1700)"
            elif "1400" in anchor_name:
                anch_str = r"\pi_{1}(1400)"
            elif "1600" in anchor_name:
                anch_str = r"\pi_{1}(1600)"
            elif "980" in anchor_name:
                anch_str = r"a_{0}(980)"
            return f"\\Delta\\phi({s_name} - {anch_str})"
        return f"\\Delta\\phi({s_name})"
    clean_name = name.replace("_", r"\_")
    if anchor_name:
        clean_anchor = anchor_name.replace("_", r"\_")
        return f"\\Delta\\phi({clean_name} - {clean_anchor})"
    return f"\\Delta\\phi({clean_name})"


def format_relative_phase_label(name: str, anchor_name: Optional[str] = None, plus_minus: bool = False) -> str:
    """Format relative phase into LaTeX math string (e.g. $\\Delta\\phi(P^{+}_{1} - D^{+}_{1})$)."""
    p_name = get_amplitude_symbol_parts(name)
    p_anch = get_amplitude_symbol_parts(anchor_name) if anchor_name else None
    prefix = "\\pm\\Delta\\phi" if plus_minus else "\\Delta\\phi"
    if p_name and p_anch:
        s_name = f"{p_name[0]}^{{{p_name[1]}}}_{{{p_name[2]}}}"
        s_anch = f"{p_anch[0]}^{{{p_anch[1]}}}_{{{p_anch[2]}}}"
        return f"${prefix}({s_name} - {s_anch})$"
    elif p_name:
        s_name = f"{p_name[0]}^{{{p_name[1]}}}_{{{p_name[2]}}}"
        if anchor_name:
            anch_str = anchor_name.replace("_", r"\_")
            if "1320" in anchor_name:
                anch_str = r"a_{2}(1320)"
            elif "1700" in anchor_name:
                anch_str = r"a_{2}(1700)"
            elif "1400" in anchor_name:
                anch_str = r"\pi_{1}(1400)"
            elif "1600" in anchor_name:
                anch_str = r"\pi_{1}(1600)"
            elif "980" in anchor_name:
                anch_str = r"a_{0}(980)"
            return f"${prefix}({s_name} - {anch_str})$"
        return f"${prefix}({s_name})$"
    clean_name = name.replace("_", r"\_")
    if anchor_name:
        clean_anchor = anchor_name.replace("_", r"\_")
        return f"${prefix}({clean_name} - {clean_anchor})$"
    return f"${prefix}({clean_name})$"


def format_reference_partial_wave_label(
    ref_info: Dict[str, Any],
    refl_char: Optional[str] = None,
    phases: Optional[Any] = None,
    include_refl_tags: bool = True,
) -> str:
    """
    Format reference partial-wave(s) into an explicit LaTeX string (e.g. '$D^{+}_{2}$' or '$D^{+}_{2}\\ (+),\\ D^{-}_{2}\\ (-)').
    """
    ref_names = ref_info.get("ref_names", {})
    if refl_char:
        ref_amp = ref_names.get(refl_char, ref_info.get("ref_name", "a_2_2"))
        return format_amplitude_label(ref_amp)

    ref_a = ref_names.get("a", ref_info.get("ref_name", "a_2_2"))
    ref_b = ref_names.get("b", "b_2_2")
    lbl_a = format_amplitude_label(ref_a) if ref_a else None
    lbl_b = format_amplitude_label(ref_b) if ref_b else None

    if phases is not None:
        phases_list = list(phases)
        has_pos = any(k.startswith("aphi_") or k.startswith("a_") for k in phases_list)
        has_neg = any(k.startswith("bphi_") or k.startswith("b_") for k in phases_list)
    else:
        has_pos = bool(ref_a)
        has_neg = bool(ref_b)

    if has_pos and has_neg and lbl_a and lbl_b:
        if include_refl_tags:
            return f"{lbl_a} (+), {lbl_b} (-)"
        return f"{lbl_a}, {lbl_b}"
    elif has_neg and not has_pos and lbl_b:
        return lbl_b
    elif lbl_a:
        return lbl_a
    elif lbl_b:
        return lbl_b
    return "$D^{+}_{2}$"


def resolve_anchor_amplitudes(
    config: Any,
    nominal_amps: Dict[str, List[float]],
    active_phases: List[str],
) -> Dict[str, Optional[str]]:
    """
    Resolve reference/anchor amplitude phase branch names for positive ('a') and negative ('b') reflectivities.
    Returns dict: {'a': pos_anchor_phi, 'b': neg_anchor_phi}.
    
    Precedence:
    1. Explicit config (config.plotting.reference_amplitude or anchor_amplitude)
    2. Auto-detection from config.waves:
       - Find wave / resonance with phase_fixed=True (e.g. D-wave a2(1320))
       - Select the dominant active amplitude in that wave (max integrated intensity)
    3. Auto-detection from active phases:
       - Wave with constant phase across mass
       - Dominant wave in nominal_amps
    """
    anchors: Dict[str, Optional[str]] = {"a": None, "b": None}

    def _normalize_phi_key(item: str) -> Tuple[str, str]:
        s = item.strip()
        if s.startswith("aphi_"):
            return "a", s
        elif s.startswith("bphi_"):
            return "b", s
        elif s.startswith("a_"):
            return "a", "aphi_" + s[2:]
        elif s.startswith("b_"):
            return "b", "bphi_" + s[2:]
        p = get_amplitude_symbol_parts(s)
        if p:
            l_sym, refl, m_str, _ = p
            inv_l = {v: k for k, v in L_SPECTROSCOPIC.items()}
            l_val = inv_l.get(l_sym)
            if l_val is not None:
                refl_char = "a" if refl == "+" else "b"
                return refl_char, f"{refl_char}phi_{l_val}_{m_str}"
        return "a", f"aphi_{s}"

    # 1. Check explicit setting in config
    explicit_val = None
    if hasattr(config, "plotting") and getattr(config, "plotting", None) is not None:
        explicit_val = getattr(config.plotting, "reference_amplitude", None) or getattr(config.plotting, "anchor_amplitude", None)

    if explicit_val:
        items = [s for s in str(explicit_val).replace(",", " ").split() if s]
        for it in items:
            refl, phi_key = _normalize_phi_key(it)
            anchors[refl] = phi_key

    # If positive anchor was given explicitly, try matching negative anchor if not set
    if anchors["a"] and not anchors["b"]:
        matching_b = "bphi_" + anchors["a"][5:]
        if matching_b in active_phases:
            anchors["b"] = matching_b

    # 2. Auto-detection for missing reflectivities
    def _integrated_intensity(phase_key: str) -> float:
        refl_prefix = "a_" if phase_key.startswith("aphi_") else "b_"
        mag_key = refl_prefix + phase_key[5:]
        vals = nominal_amps.get(mag_key, [])
        return sum(v**2 for v in vals) if vals else 0.0

    waves = getattr(config, "waves", []) if config else []

    for refl in ("a", "b"):
        if anchors[refl] and anchors[refl] in active_phases:
            continue

        cands = [p for p in active_phases if p.startswith(f"{refl}phi_")]
        if not cands:
            anchors[refl] = None
            continue

        # Strategy A: Waves with phase_fixed
        fixed_wave_cands: List[str] = []
        for w in waves:
            is_fixed = getattr(w, "phase_fixed", False)
            if not is_fixed:
                for r in getattr(w, "resonances", []):
                    if getattr(r, "phase_fixed", False) or getattr(r, "fixed_phase", None) is not None:
                        is_fixed = True
                        break
            if is_fixed:
                w_l = getattr(w, "l", None)
                for c in cands:
                    if extract_amplitude_l(c) == w_l and c not in fixed_wave_cands:
                        fixed_wave_cands.append(c)

        if fixed_wave_cands:
            anchors[refl] = max(fixed_wave_cands, key=_integrated_intensity)
            continue

        # Strategy B: Phase that is constant across mass (e.g. constant 0.0)
        const_cands = []
        for c in cands:
            vals = nominal_amps.get(c, [])
            if len(vals) > 1 and np.allclose(vals, vals[0], atol=1e-6):
                const_cands.append(c)
        if const_cands:
            anchors[refl] = max(const_cands, key=_integrated_intensity)
            continue

        # Strategy C: Dominant wave overall in this reflectivity
        anchors[refl] = max(cands, key=_integrated_intensity)

    return anchors


def resolve_reference_phase_info(
    config: Any,
    nominal_amps: Dict[str, List[float]],
    active_phases: List[str],
    masses: List[float],
    nominal_params: Optional[Dict[str, float]] = None,
    mcmc_params: Optional[Dict[str, List[float]]] = None,
    mcmc_samples: Optional[Dict[str, Dict[float, List[float]]]] = None,
) -> Dict[str, Any]:
    """
    Resolve the phase reference (the a2(1320) partial wave with the biggest magnitude, e.g. D+_2 / a_2_2)
    for positive ('a') and negative ('b') reflectivities.
    
    Physics logic:
    - Relative phases should always be referenced to the dominant partial wave of the reference resonance
      (specifically a2(1320), where D+_2 has the largest magnitude for positive reflectivity and D-_2 for negative).
    - The actual fitted phase array of this partial wave (nominal_amps[anchor]) is subtracted so that
      the relative phase DeltaPhi(wave - D+_2) represents the true physical phase difference.
    - For the reference wave itself, DeltaPhi == 0 everywhere; it is placed in omit_phases so it is
      not drawn as a redundant zero line, while all other waves (P-waves, S-waves, and other D-wave projections)
      are drawn relative to D+_2 (or D-_2).
    """
    anchors = resolve_anchor_amplitudes(config, nominal_amps, active_phases)
    waves = getattr(config, "waves", []) if config else []

    # 1. Identify reference wave and resonance
    explicit_ref = None
    if hasattr(config, "plotting") and getattr(config, "plotting", None) is not None:
        explicit_ref = getattr(config.plotting, "reference_amplitude", None) or getattr(config.plotting, "anchor_amplitude", None)

    ref_wave = None
    ref_res = None

    if explicit_ref:
        exp_s = str(explicit_ref).strip().lower()
        for w in waves:
            if getattr(w, "name", "").lower() in exp_s or f"l={getattr(w, 'l', -1)}" in exp_s or f"l{getattr(w, 'l', -1)}" in exp_s:
                ref_wave = w
                break
            for r in getattr(w, "resonances", []):
                if getattr(r, "name", "").lower() in exp_s:
                    ref_wave = w
                    ref_res = r
                    break
            if ref_wave:
                break

    if ref_wave is None:
        for w in waves:
            if getattr(w, "phase_fixed", False):
                ref_wave = w
                break
            for r in getattr(w, "resonances", []):
                if getattr(r, "phase_fixed", False) or getattr(r, "fixed_phase", None) is not None:
                    ref_wave = w
                    ref_res = r
                    break
            if ref_wave:
                break

    if ref_wave is None and anchors.get("a"):
        anc_l = extract_amplitude_l(anchors["a"])
        if anc_l is not None:
            for w in waves:
                if getattr(w, "l", None) == anc_l:
                    ref_wave = w
                    break

    # 2. Determine model type and multi-resonance status
    is_two_bw = False
    is_single_bw = False
    ref_wave_l = getattr(ref_wave, "l", None) if ref_wave else None

    if ref_wave:
        model_str = str(getattr(ref_wave, "model", "")).lower()
        n_res = len(getattr(ref_wave, "resonances", []))
        if "two_breit_wigner" in model_str or n_res >= 2:
            is_two_bw = True
        elif "breit_wigner" in model_str or n_res == 1:
            is_single_bw = True

    # 3. Find reference resonance parameters (M1, W1)
    if ref_wave and getattr(ref_wave, "resonances", []):
        if ref_res is None:
            for r in ref_wave.resonances:
                if getattr(r, "phase_fixed", False) or getattr(r, "fixed_phase", None) is not None:
                    ref_res = r
                    break
            if ref_res is None:
                ref_res = ref_wave.resonances[0]

    res_name = getattr(ref_res, "name", "a2_1320") if ref_res else "a2_1320"
    m1_val = getattr(ref_res, "mass", 1.318) if ref_res else 1.318
    w1_val = getattr(ref_res, "width", 0.107) if ref_res else 0.107
    if w1_val is None:
        w1_val = 0.107

    if nominal_params:
        for pk, pv in nominal_params.items():
            pk_lower = pk.lower()
            if ("mass1" in pk_lower or (res_name.lower() in pk_lower and "mass" in pk_lower)) and not pk.endswith("_error"):
                m1_val = pv
                break
        for pk, pv in nominal_params.items():
            pk_lower = pk.lower()
            if ("width1" in pk_lower or (res_name.lower() in pk_lower and "width" in pk_lower)) and not pk.endswith("_error"):
                w1_val = pv
                break

    has_barrier = False
    barrier_l = ref_wave_l if ref_wave_l is not None else 2
    barrier_r = 5.0677
    if ref_wave:
        model_str = str(getattr(ref_wave, "model", "")).lower()
        if "barrier" in model_str or getattr(ref_wave, "barrier", False) or (getattr(ref_wave, "barrier_l", None) is not None):
            has_barrier = True
        if getattr(ref_wave, "barrier_l", None) is not None:
            barrier_l = ref_wave.barrier_l
        if getattr(ref_wave, "barrier_radius", None) is not None:
            barrier_r = ref_wave.barrier_radius

    def _calc_bw_phase(m: float, m0: float, w0: float) -> float:
        if has_barrier:
            gamma = dynamic_breit_wigner_width(m, m0, w0, l=barrier_l, r=barrier_r)
        else:
            gamma = w0
        return math.atan2(m0 * gamma, m0 * m0 - m * m)

    # 4. Resolve partial wave amplitude reference names
    def _phi_to_amp(phi_key: Optional[str]) -> Optional[str]:
        if not phi_key:
            return None
        if phi_key.startswith("aphi_"):
            return "a_" + phi_key[5:]
        if phi_key.startswith("bphi_"):
            return "b_" + phi_key[5:]
        return phi_key

    ref_amp_a = _phi_to_amp(anchors.get("a")) or "a_2_2"
    ref_amp_b = _phi_to_amp(anchors.get("b")) or "b_2_2"
    ref_names = {"a": ref_amp_a, "b": ref_amp_b}
    ref_name = ref_amp_a

    # 5. Compute nominal reference phase from the anchor partial waves
    anc_a = anchors.get("a")
    anc_b = anchors.get("b")

    if anc_a and anc_a in nominal_amps:
        nom_ref_a = list(nominal_amps[anc_a])
    elif ref_wave and (is_two_bw or is_single_bw):
        nom_ref_a = [_calc_bw_phase(m, m1_val, w1_val) for m in masses]
    else:
        nom_ref_a = [0.0] * len(masses)

    if anc_b and anc_b in nominal_amps:
        nom_ref_b = list(nominal_amps[anc_b])
    else:
        nom_ref_b = nom_ref_a

    nominal_ref = {"a": nom_ref_a, "b": nom_ref_b}

    # 6. Compute MCMC reference phase if samples available
    ref_mcmc_dict = None
    if mcmc_params:
        m1_samps = None
        w1_samps = None
        for pk, pvals in mcmc_params.items():
            pk_lower = pk.lower()
            if ("mass1" in pk_lower or (res_name.lower() in pk_lower and "mass" in pk_lower)) and not pk.endswith("_error"):
                m1_samps = pvals
                break
        for pk, pvals in mcmc_params.items():
            pk_lower = pk.lower()
            if ("width1" in pk_lower or (res_name.lower() in pk_lower and "width" in pk_lower)) and not pk.endswith("_error"):
                w1_samps = pvals
                break
        if m1_samps and w1_samps and len(m1_samps) == len(w1_samps):
            ref_mcmc_dict = {}
            for m in masses:
                ref_mcmc_dict[round(m, 5)] = [
                    _calc_bw_phase(m, m0, w0) for m0, w0 in zip(m1_samps, w1_samps)
                ]

    ref_mcmc_a = None
    ref_mcmc_b = None
    if mcmc_samples:
        if anc_a and anc_a in mcmc_samples:
            ref_mcmc_a = mcmc_samples[anc_a]
        if anc_b and anc_b in mcmc_samples:
            ref_mcmc_b = mcmc_samples[anc_b]

    if ref_mcmc_a is None and ref_mcmc_dict is not None:
        ref_mcmc_a = ref_mcmc_dict
    if ref_mcmc_b is None and ref_mcmc_dict is not None:
        ref_mcmc_b = ref_mcmc_dict

    mcmc_ref = {"a": ref_mcmc_a, "b": ref_mcmc_b}

    # 7. Determine omit_phases
    # The reference partial wave itself has identically DeltaPhi == 0; omit it from non-zero relative phase curves
    omit_phases: set[str] = set()
    if is_single_bw:
        if ref_wave_l is not None:
            omit_phases = {p for p in active_phases if extract_amplitude_l(p) == ref_wave_l}
    for refl in ("a", "b"):
        if anchors.get(refl):
            omit_phases.add(anchors[refl])

    return {
        "nominal_ref": nominal_ref,
        "mcmc_ref": mcmc_ref,
        "ref_name": ref_name,
        "ref_names": ref_names,
        "is_two_bw": is_two_bw,
        "is_single_bw": is_single_bw,
        "ref_wave_l": ref_wave_l,
        "omit_phases": omit_phases,
        "anchors": anchors,
    }


def compute_relative_phase_data(
    masses: List[float],
    nominal_phase: List[float],
    anchor_nominal_phase: Optional[List[float]],
    mcmc_phase_dict: Optional[Dict[float, List[float]]] = None,
    anchor_mcmc_dict: Optional[Dict[float, List[float]]] = None,
) -> Tuple[List[float], Optional[Dict[float, List[float]]]]:
    """
    Compute relative phase: DeltaPhi = wrap_angle_pi(Phi - Phi_anchor).
    Propagates MCMC sample-by-sample differences to properly preserve posterior correlations.
    """
    if anchor_nominal_phase is None:
        return list(nominal_phase), mcmc_phase_dict

    y_nom_arr = np.asarray(nominal_phase, dtype=float)
    y_anc_arr = np.asarray(anchor_nominal_phase, dtype=float)
    min_len = min(len(y_nom_arr), len(y_anc_arr))
    y_rel_arr = wrap_angle_pi(y_nom_arr[:min_len] - y_anc_arr[:min_len])
    y_rel = y_rel_arr.tolist()

    mcmc_rel: Optional[Dict[float, List[float]]] = None
    if mcmc_phase_dict is not None:
        mcmc_rel = {}
        anc_nom_lookup = {round(m, 5): v for m, v in zip(masses, y_anc_arr)}

        for m_val, samps in mcmc_phase_dict.items():
            if not samps:
                mcmc_rel[m_val] = []
                continue
            s_arr = np.asarray(samps, dtype=float)
            anc_samps = anchor_mcmc_dict.get(m_val) if anchor_mcmc_dict is not None else None
            if anc_samps and len(anc_samps) > 0:
                anc_s_arr = np.asarray(anc_samps, dtype=float)
                n_eval = min(len(s_arr), len(anc_s_arr))
                diff = s_arr[:n_eval] - anc_s_arr[:n_eval]
            else:
                m_round = round(m_val, 5)
                anc_nom_val = anc_nom_lookup.get(m_round, 0.0)
                diff = s_arr - anc_nom_val
            mcmc_rel[m_val] = wrap_angle_pi(diff).tolist()

    return y_rel, mcmc_rel


def _safe_mathtext_label(label: str, fallback_raw: str = "") -> str:
    """Validate mathtext with matplotlib parser; if invalid, return escaped plain text."""
    try:
        if "$" in label:
            import matplotlib.mathtext as mathtext
            mathtext.MathTextParser("path").parse(label)
        return label
    except Exception:
        fallback = fallback_raw or label
        return fallback.replace("$", "").replace("_", " ")


def format_mcmc_parameter_label(pname: str, short: bool = False) -> str:
    """
    Format internal parameter names (MD_*, MI_*, chi2) into valid, clean LaTeX mathtext labels.
    Guarantees no double-subscripts or invalid characters that would cause mathtext parse errors.
    """
    if pname == "chi2":
        return r"$\chi^2$"

    # 1. Piecewise mass-independent parameters: MI_<mass>_<amp>
    if pname.startswith("MI_"):
        m_parts = pname.split("_")
        if len(m_parts) >= 4:
            mb = m_parts[1]
            amp = "_".join(m_parts[2:])
            try:
                mb_val = f"{float(mb):.2f}"
            except ValueError:
                mb_val = mb
            sym = format_amplitude_symbol(amp)
            if short:
                lbl = f"${sym}$ ({mb_val})"
            else:
                lbl = f"${sym}$ ({mb_val} GeV)"
            return _safe_mathtext_label(lbl, pname)
        return _safe_mathtext_label(pname.replace("_", r"\_"), pname)

    # 2. Shared resonance parameters (Mass, Width, couplings like g_etapi)
    if "Mass" in pname or "Width" in pname:
        m_idx = re.search(r"(Mass|Width)([0-9]*)", pname)
        if m_idx:
            par_type = m_idx.group(1)
            idx = m_idx.group(2)
            sym = "M" if par_type == "Mass" else r"\Gamma"
            idx_sub = f"_{{{idx}}}" if idx else ""
            if short:
                return _safe_mathtext_label(f"${sym}{idx_sub}$", pname)
            clean_res = pname.replace("MD_", "").replace("shared_", "")
            clean_res = re.sub(r"_(Mass|Width)[0-9]*$", "", clean_res).replace("_", r"\_")
            lbl = f"${sym}{idx_sub}$ ({clean_res})" if clean_res else f"${sym}{idx_sub}$"
            return _safe_mathtext_label(lbl, pname)

    if "_g_" in pname or pname.endswith("_g"):
        g_match = re.search(r"g_([a-zA-Z0-9]+)", pname)
        if g_match:
            channel = g_match.group(1).replace("pi", r"\pi").replace("eta", r"\eta")
            return _safe_mathtext_label(f"$g_{{{channel}}}$", pname)

    # 3. Global phase parameters
    if "globalphase" in pname or "phase" in pname:
        refl_match = re.search(r"shared([ab])_", pname)
        refl = "^{(+)}" if refl_match and refl_match.group(1) == "a" else ("^{(-)}" if refl_match and refl_match.group(1) == "b" else "")
        idx_match = re.search(r"globalphase([0-9]+)", pname)
        idx_sub = f"_{{{idx_match.group(1)}}}" if idx_match else ""
        clean_res = pname.replace("MD_", "").replace("shareda_", "").replace("sharedb_", "").replace("shared_", "")
        clean_res = re.sub(r"_globalphase[0-9]*", "", clean_res).replace("_phase", "").replace("_", r"\_")
        if short:
            lbl = f"$\\phi{idx_sub}{refl}$" if not clean_res else f"$\\phi{idx_sub}{refl}({clean_res})$"
        else:
            lbl = f"$\\phi{idx_sub}{refl}({clean_res})$" if clean_res else f"$\\phi{idx_sub}{refl}$"
        return _safe_mathtext_label(lbl, pname)

    # 4. Wave coupling parameters: MD_<refl>_<L>_<M>_<wave>_<param>
    # e.g. MD_a_2_-2_a2_1320_1700_k1, MD_a_2_0_a2_1320_k, MD_b_2_1_k2
    match_wave = re.match(r"^MD_([ab]|aphi|bphi)_([0-9]+)_(-?[0-9]+)(?:_(.*))?$", pname)
    if match_wave:
        refl_pre = match_wave.group(1)
        l_str = match_wave.group(2)
        m_str = match_wave.group(3)
        rest = match_wave.group(4) or ""
        amp_key = f"{refl_pre}_{l_str}_{m_str}"
        sym = format_amplitude_symbol(amp_key)

        if not rest:
            lbl = f"${sym}$"
        elif rest == "k" or rest.endswith("_k"):
            lbl = f"$k({sym})$"
        elif rest == "k1" or rest.endswith("_k1") or rest.endswith("_1"):
            lbl = f"$k_1({sym})$"
        elif rest == "k2" or rest.endswith("_k2") or rest.endswith("_2"):
            lbl = f"$k_2({sym})$"
        elif "re_" in rest or "im_" in rest:
            re_m = re.search(r"(re|im)_([0-9]+)", rest)
            if re_m:
                part_kind = r"\mathrm{Re}" if re_m.group(1) == "re" else r"\mathrm{Im}"
                c_idx = re_m.group(2)
                lbl = f"${part_kind}_{{{c_idx}}}({sym})$"
            else:
                clean_rest = rest.replace("_", r"\_")
                lbl = f"$k_{{{clean_rest}}}({sym})$"
        else:
            clean_rest = rest.replace("_", r"\_")
            lbl = f"$k_{{{clean_rest}}}({sym})$"
        return _safe_mathtext_label(lbl, pname)

    # 5. Generic fallback
    clean_fallback = pname.replace("MD_", "").replace("_", r"\_")
    return _safe_mathtext_label(f"${clean_fallback}$", pname)


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

    included_L = getattr(chi2_cfg, "included_L", None)
    if included_L is not None and len(included_L) > 0:
        if L not in included_L:
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


def sort_amp_key(name: str) -> Tuple[int, int, int, int]:
    """
    Sort key for partial wave amplitude and phase names:
    Orders by:
      1. Reflectivity: '+' ('a') before '-' ('b')
      2. Angular momentum L: 0, 1, 2, ...
      3. Projection M: ..., -2, -1, 0, 1, 2, ...
      4. Kind: amplitude (0) before phase (1)
    """
    parts = get_amplitude_symbol_parts(name)
    if parts:
        l_sym, refl, m_str, wave_kind = parts
        refl_order = 0 if refl == "+" else 1
        phi_order = 1 if wave_kind == "phase" else 0
        l_val = extract_amplitude_l(name)
        l_num = l_val if l_val is not None else 999
        try:
            return (refl_order, l_num, int(m_str), phi_order)
        except ValueError:
            pass
    return (9, 999, 999, 0)


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


def compute_weighted_percentile(arr: np.ndarray, q: float, weights: Optional[np.ndarray] = None) -> float:
    """Compute percentile q in [0, 100] taking into account sample weights if provided."""
    if len(arr) == 0:
        return 0.0
    if weights is None or len(weights) != len(arr) or np.sum(weights) <= 0:
        return float(np.percentile(arr, q))
    sorter = np.argsort(arr)
    s_arr = arr[sorter]
    s_w = weights[sorter]
    cum_w = np.cumsum(s_w)
    tot_w = cum_w[-1]
    if tot_w <= 0:
        return float(np.percentile(arr, q))
    target = (q / 100.0) * tot_w
    return float(np.interp(target, cum_w, s_arr))


def extract_credible_bounds(
    masses: List[float],
    mcmc_dict: Dict[float, List[float]],
    nom_vals: List[float],
    percentiles: Tuple[float, float] = (16.0, 84.0),
    squared: bool = False,
    mcmc_weights: Optional[Dict[float, List[float]]] = None,
) -> Tuple[List[float], List[float]]:
    """Extract lower and upper percentiles for credible interval band across mass bins."""
    import numpy as np
    lo, hi = [], []
    p_lo, p_hi = percentiles
    for idx, m in enumerate(masses):
        m_round = round(m, 4)
        samps = mcmc_dict.get(m_round)
        wts = mcmc_weights.get(m_round) if mcmc_weights else None
        if samps is None:
            for mk, mv in mcmc_dict.items():
                if abs(mk - m_round) < 0.005:
                    samps = mv
                    if mcmc_weights:
                        wts = mcmc_weights.get(mk)
                    break
        if samps:
            vals = np.array([v**2 for v in samps] if squared else samps, dtype=float)
            w_arr = np.array(wts, dtype=float) if wts is not None else None
            lo.append(compute_weighted_percentile(vals, p_lo, w_arr))
            hi.append(compute_weighted_percentile(vals, p_hi, w_arr))
        else:
            ref = nom_vals[idx]**2 if squared else nom_vals[idx]
            lo.append(ref)
            hi.append(ref)
    return lo, hi


def plot_posterior_density_band(
    ax: Any,
    masses: List[float],
    mcmc_dict: Dict[float, List[float]],
    is_mass_dep: bool,
    color: Any = "black",
    squared: bool = False,
    n_levels: int = 20,
    base_alpha: Optional[float] = None,
    max_alpha: float = 0.60,
    fit_step: float = 0.04,
    zorder: int = 3,
    mirrored: bool = False,
    mcmc_weights: Optional[Dict[float, List[float]]] = None,
    linewidth: float = 1.0,
    line_alpha: Optional[float] = None,
    max_traces: int = 500,
    **kwargs: Any,
) -> None:
    """
    Draw individual propagated fit results as transparent lines with visible width.
    Where fit results concentrate in a certain region, overlapping lines naturally
    accumulate opacity (optical density blending), producing an intense modal core
    while outliers and wide spreads appear faint.
    """
    import numpy as np
    from matplotlib.collections import LineCollection
    import matplotlib.colors as mcolors

    if not mcmc_dict:
        return

    has_valid_samples = any(len(samps) > 0 for samps in mcmc_dict.values())
    if not has_valid_samples:
        return

    # 1. Collect aligned sample values across valid mass bins
    valid_masses: List[float] = []
    matrix_vals: List[np.ndarray] = []
    matrix_wts: List[Optional[np.ndarray]] = []

    for m in masses:
        m_round = round(m, 4)
        samps = mcmc_dict.get(m_round)
        wts = mcmc_weights.get(m_round) if mcmc_weights else None
        if samps is None:
            for mk, mv in mcmc_dict.items():
                if abs(mk - m_round) < 0.005:
                    samps = mv
                    if mcmc_weights:
                        wts = mcmc_weights.get(mk)
                    break
        if samps and len(samps) > 0:
            valid_masses.append(m)
            arr = np.array(samps, dtype=float)
            if squared:
                vals = arr**2
            elif mirrored:
                vals = -arr
            else:
                vals = arr
            matrix_vals.append(vals)
            matrix_wts.append(np.array(wts, dtype=float) if wts is not None and len(wts) == len(samps) else None)

    if not matrix_vals or not valid_masses:
        return

    n_samples = min(len(v) for v in matrix_vals)
    if n_samples <= 0:
        return

    # Subsample if sample count is excessively large to avoid inflating PDF/raster size
    if n_samples > max_traces:
        sub_indices = np.linspace(0, n_samples - 1, max_traces, dtype=int)
    else:
        sub_indices = np.arange(n_samples)

    n_traces = len(sub_indices)
    Y = np.column_stack([v[sub_indices] for v in matrix_vals])  # shape: (n_traces, n_masses)

    # 2. Calibrate alpha per trace to achieve natural optical density blending
    if line_alpha is not None:
        eff_alpha = line_alpha
    else:
        eff_alpha = max(0.015, min(0.08, 6.0 / max(1, n_traces)))

    # Per-trace weights if available
    w_vec = None
    for w_cand in matrix_wts:
        if w_cand is not None and len(w_cand) >= n_samples:
            w_vec = w_cand[sub_indices]
            break

    if w_vec is not None and np.sum(w_vec) > 0:
        w_norm = w_vec / np.mean(w_vec)
        alphas = np.clip(eff_alpha * np.sqrt(w_norm), 0.005, 0.50)
    else:
        alphas = np.full(n_traces, eff_alpha)

    # 3. Build line segments
    if is_mass_dep:
        x_arr = np.array(valid_masses, dtype=float)
        sort_idx = np.argsort(x_arr)
        xs = x_arr[sort_idx]
        Ys = Y[:, sort_idx]
        unique_mask = np.concatenate(([True], np.diff(xs) > 1e-7))
        xs = xs[unique_mask]
        Ys = Ys[:, unique_mask]

        if len(xs) >= 4:
            num_points = max(100, min(250, len(xs) * 10))
            x_dense = np.linspace(xs[0], xs[-1], num_points)
            try:
                from scipy.interpolate import CubicSpline
                cs = CubicSpline(xs, Ys, axis=1, bc_type="natural")
                Y_dense = cs(x_dense)
                segments = [np.column_stack([x_dense, Y_dense[s]]) for s in range(n_traces)]
            except Exception:
                segments = [np.column_stack([xs, Ys[s]]) for s in range(n_traces)]
        else:
            segments = [np.column_stack([xs, Ys[s]]) for s in range(n_traces)]
    else:
        edges = compute_bin_edges(valid_masses, default_step=fit_step)
        x_step = np.repeat(edges, 2)[1:-1]
        segments = [np.column_stack([x_step, np.repeat(Y[s], 2)]) for s in range(n_traces)]

    # 4. Construct LineCollection with transparent alphas and visible linewidth
    base_rgb = mcolors.to_rgb(color)
    rgba = np.zeros((n_traces, 4))
    rgba[:, :3] = base_rgb
    rgba[:, 3] = alphas

    lc = LineCollection(segments, colors=rgba, linewidths=linewidth, zorder=zorder)
    ax.add_collection(lc)

    # Update axis data limits
    try:
        ax.update_datalim(lc.get_datalim(ax.transData))
    except Exception:
        pass


def plot_phase_with_ambiguity(
    ax: Any,
    masses: List[float],
    y_nom: List[float],
    is_mass_dep: bool,
    color: Any,
    label: Optional[str] = None,
    label_neg: Optional[str] = None,
    mcmc_dict: Optional[Dict[float, List[float]]] = None,
    fit_step: float = 0.04,
    zorder: int = 4,
    show_ambiguity: bool = True,
    linewidth: float = 2.0,
    alpha_band: float = 0.32,
    base_alpha: Optional[float] = None,
    max_alpha: Optional[float] = None,
    linestyle_pos: str = "-",
    linestyle_neg: str = "--",
    density_band: bool = True,
    mcmc_weights: Optional[Dict[float, List[float]]] = None,
) -> None:
    """
    Plot fitted phase vs mass including both positive and negative branches
    to represent partial wave sign ambiguities (phi -> -phi).
    Draws primary branch (solid) and sign-ambiguous mirrored branch (dashed)
    along with their MCMC credible intervals or continuous posterior density bands.
    """
    import numpy as np

    y_nom_arr = np.array(y_nom, dtype=float)
    is_trivial_zero = bool(np.all(np.abs(y_nom_arr) < 1e-9)) and (
        mcmc_dict is None or all(np.all(np.abs(v) < 1e-9) for v in mcmc_dict.values() if len(v) > 0)
    )

    eff_max_alpha = max_alpha if max_alpha is not None else 0.60
    eff_base_alpha = base_alpha if base_alpha is not None else 0.20

    if is_mass_dep:
        if len(masses) >= 4:
            y_unwrapped = np.unwrap(y_nom_arr)
            m_s, y_s_unwrapped = smooth_curve(masses, y_unwrapped, num_points=200)
            y_s = wrap_angle_pi(y_s_unwrapped)
        else:
            m_s, y_s = np.array(masses, dtype=float), y_nom_arr

        # Split across branch cuts to avoid vertical connecting streaks
        m_plot, y_plot = split_phase_jumps(m_s, y_s)

        # Primary branch
        ax.plot(m_plot, y_plot, linestyle=linestyle_pos, linewidth=linewidth, color=color, label=label, zorder=zorder)

        # Primary MCMC credible band
        if mcmc_dict is not None:
            if density_band:
                plot_posterior_density_band(
                    ax=ax, masses=masses, mcmc_dict=mcmc_dict, is_mass_dep=True,
                    color=color, squared=False, n_levels=20, base_alpha=eff_base_alpha,
                    max_alpha=eff_max_alpha, fit_step=fit_step, zorder=zorder - 1, mirrored=False,
                    mcmc_weights=mcmc_weights
                )
            else:
                lo_68, hi_68 = extract_credible_bounds(masses, mcmc_dict, y_nom, (16.0, 84.0), mcmc_weights=mcmc_weights)
                if len(masses) >= 4:
                    m_s68, lo_s = smooth_curve(masses, np.unwrap(lo_68), num_points=200)
                    m_s68, hi_s = smooth_curve(masses, np.unwrap(hi_68), num_points=200)
                    ax.fill_between(m_s68, wrap_angle_pi(lo_s), wrap_angle_pi(hi_s), color=color, alpha=alpha_band, zorder=zorder - 1)
                else:
                    ax.fill_between(masses, lo_68, hi_68, color=color, alpha=alpha_band, zorder=zorder - 1)

        # Sign-ambiguous branch (phi -> -phi)
        if show_ambiguity and not is_trivial_zero:
            m_plot_neg, y_plot_neg = split_phase_jumps(m_s, -y_s)
            ax.plot(m_plot_neg, y_plot_neg, linestyle=linestyle_neg, linewidth=linewidth, color=color, label=label_neg, zorder=zorder)
            if mcmc_dict is not None:
                if density_band:
                    plot_posterior_density_band(
                        ax=ax, masses=masses, mcmc_dict=mcmc_dict, is_mass_dep=True,
                        color=color, squared=False, n_levels=20, base_alpha=eff_base_alpha,
                        max_alpha=eff_max_alpha, fit_step=fit_step, zorder=zorder - 1, mirrored=True,
                        mcmc_weights=mcmc_weights
                    )
                else:
                    if len(masses) >= 4:
                        ax.fill_between(m_s68, -wrap_angle_pi(hi_s), -wrap_angle_pi(lo_s), color=color, alpha=alpha_band, zorder=zorder - 1)
                    else:
                        ax.fill_between(masses, -hi_68, -lo_68, color=color, alpha=alpha_band, zorder=zorder - 1)

    else:
        edges = compute_bin_edges(masses, default_step=fit_step)

        # Primary branch (stairs or step)
        if hasattr(ax, "stairs"):
            ax.stairs(y_nom, edges, linestyle=linestyle_pos, linewidth=linewidth, color=color, label=label, zorder=zorder)
        else:
            ax.step(masses, y_nom, where="mid", linestyle=linestyle_pos, linewidth=linewidth, color=color, label=label, zorder=zorder)

        if mcmc_dict is not None:
            if density_band:
                plot_posterior_density_band(
                    ax=ax, masses=masses, mcmc_dict=mcmc_dict, is_mass_dep=False,
                    color=color, squared=False, n_levels=20, base_alpha=eff_base_alpha,
                    max_alpha=eff_max_alpha, fit_step=fit_step, zorder=zorder - 1, mirrored=False,
                    mcmc_weights=mcmc_weights
                )
            else:
                lo_68, hi_68 = extract_credible_bounds(masses, mcmc_dict, y_nom, (16.0, 84.0), mcmc_weights=mcmc_weights)
                x_step = np.repeat(edges, 2)[1:-1]
                lo_step = np.repeat(lo_68, 2)
                hi_step = np.repeat(hi_68, 2)
                ax.fill_between(x_step, lo_step, hi_step, color=color, alpha=alpha_band, zorder=zorder - 1)

        # Sign-ambiguous branch (phi -> -phi)
        if show_ambiguity and not is_trivial_zero:
            y_neg = [-v for v in y_nom]
            if hasattr(ax, "stairs"):
                ax.stairs(y_neg, edges, linestyle=linestyle_neg, linewidth=linewidth, color=color, label=label_neg, zorder=zorder)
            else:
                ax.step(masses, y_neg, where="mid", linestyle=linestyle_neg, linewidth=linewidth, color=color, label=label_neg, zorder=zorder)

            if mcmc_dict is not None:
                if density_band:
                    plot_posterior_density_band(
                        ax=ax, masses=masses, mcmc_dict=mcmc_dict, is_mass_dep=False,
                        color=color, squared=False, n_levels=20, base_alpha=eff_base_alpha,
                        max_alpha=eff_max_alpha, fit_step=fit_step, zorder=zorder - 1, mirrored=True,
                        mcmc_weights=mcmc_weights
                    )
                else:
                    ax.fill_between(x_step, -hi_step, -lo_step, color=color, alpha=alpha_band, zorder=zorder - 1)


def resolve_unnormalized_moments_file(
    config: FitConfig,
    rundir: Optional[RunDirectory] = None,
    unnormalized_moments_file: Optional[Path | str] = None,
) -> Optional[Path]:
    """Resolve the unnormalized moments ROOT file from arguments, config, or dataset directory."""
    if unnormalized_moments_file:
        cand = Path(unnormalized_moments_file)
        if cand.is_file():
            return cand
    if getattr(config, "plotting", None) and getattr(config.plotting, "unnormalized_moments_file", None):
        cand = Path(config.plotting.unnormalized_moments_file)
        if cand.is_file():
            return cand
    if getattr(config, "dataset", None) and getattr(config.dataset, "moments_dir", None):
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
                    return cand
            if sdir.is_dir():
                matches = sorted(sdir.glob("*UnnormalizedMoments*.root"))
                if matches:
                    return matches[0]

    return None


def load_h0_moment_data(
    unnorm_path: Optional[Path],
    cache_json: Optional[Path] = None,
) -> Optional[Dict[str, List[float]]]:
    """Load H_0_0_0 (acceptance corrected yield) moment data points and uncertainties."""
    if cache_json and cache_json.is_file():
        try:
            with open(cache_json, "r") as f:
                cdata = json.load(f)
                if "h0_data" in cdata and cdata["h0_data"].get("mass") and cdata["h0_data"].get("val"):
                    return cdata["h0_data"]
                data_moms = cdata.get("data", {})
                if "H_0_0_0" in data_moms and data_moms["H_0_0_0"].get("mass"):
                    return data_moms["H_0_0_0"]
        except Exception:
            pass

    if unnorm_path and unnorm_path.is_file():
        try:
            import uproot
            with uproot.open(unnorm_path) as uf:
                target_key = None
                for k in uf.keys():
                    if "unnormalized_moments_0_0_0" in k or "H_0_0_0" in k:
                        target_key = k
                        break
                if target_key:
                    gr = uf[target_key]
                    x = list(gr.member("fX")) if hasattr(gr, "member") else list(gr.values()[0])
                    y = list(gr.member("fY")) if hasattr(gr, "member") else list(gr.values()[1])
                    ey = list(gr.member("fEY")) if hasattr(gr, "member") else [0.0] * len(y)
                    return {"mass": [float(v) for v in x], "val": [float(v) for v in y], "err": [float(v) for v in ey]}
        except Exception:
            pass

    return None


def _read_nominal_params_from_root(nom_path: Path) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Read mass-dependent nominal parameter values and Minuit errors from a ROOT fit result file."""
    vals: Dict[str, float] = {}
    errs: Dict[str, float] = {}
    if not nom_path or not nom_path.is_file() or nom_path.stat().st_size < 100:
        return vals, errs
    try:
        import uproot
        with uproot.open(nom_path) as f:
            for tree_name in ("mass_dependent_params", "all_params"):
                if tree_name not in f:
                    continue
                ptree = f[tree_name]
                for k in ptree.keys():
                    if k.endswith("_error"):
                        base_k = k[:-6]
                        if tree_name == "all_params" and not base_k.startswith("MD_"):
                            continue
                        if base_k not in errs:
                            try:
                                errs[base_k] = float(ptree[k].array()[0])
                            except Exception:
                                pass
                    else:
                        if tree_name == "all_params" and not k.startswith("MD_"):
                            continue
                        if k not in vals:
                            try:
                                vals[k] = float(ptree[k].array()[0])
                            except Exception:
                                pass
    except Exception:
        pass
    return vals, errs


def load_or_export_amplitude_data(
    config: FitConfig,
    rundir: Optional[RunDirectory] = None,
    nominal_root: Optional[Path | str] = None,
    mcmc_root: Optional[Path | str] = None,
    unnormalized_moments_file: Optional[Path | str] = None,
    cache_dir: Optional[Path | str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Load cached amplitude JSON or export/extract it from nominal ROOT and MCMC ROOT files.
    Tries cached JSON first; if missing or stale, reads trees with uproot or falls back to ExportAmplitudes.C.
    Returns dict containing:
      - 'masses': list of float mass bin centers
      - 'nominal': dict of {amp_name: [nominal_values]}
      - 'mcmc': dict of {amp_name: {mass_val: [sample_values]}}
      - 'mcmc_weights': dict of {amp_name: {mass_val: [sample_weights]}}
      - 'mcmc_params': dict of {param_name: [sample_values]}
      - 'nominal_params': dict of {param_name: nominal_value}
      - 'nominal_param_errors': dict of {param_name: minuit_error}
      - 'h0_data': dict of {'mass': [...], 'val': [...], 'err': [...]} (or None)
    """
    # 1. Resolve ROOT file paths
    nom_path = Path(nominal_root) if nominal_root else (rundir.nominal_best_root if rundir else None)
    if nom_path and not nom_path.is_file() and rundir:
        legacy_best = rundir.root / "central_best_params.root"
        if legacy_best.is_file():
            nom_path = legacy_best
        elif rundir.mcmc_propagated_root.is_file():
            nom_path = rundir.mcmc_propagated_root

    mcmc_path = Path(mcmc_root) if mcmc_root else (rundir.mcmc_propagated_root if (rundir and rundir.mcmc_propagated_root.is_file()) else None)
    if mcmc_path and not mcmc_path.is_file():
        mcmc_path = None

    unnorm_path = resolve_unnormalized_moments_file(config, rundir, unnormalized_moments_file)

    # 2. Check candidate cache locations
    candidate_cache_dirs: List[Path] = []
    if cache_dir:
        candidate_cache_dirs.append(Path(cache_dir))
    if rundir:
        candidate_cache_dirs.extend([
            rundir.root / ".cache",
            rundir.plots_dir / ".cache",
            rundir.plots_dir / "amplitudes" / ".cache",
        ])

    target_cache_dir = candidate_cache_dirs[0] if candidate_cache_dirs else Path(".cache")

    for cdir in candidate_cache_dirs:
        jpath = cdir / "amplitudes_data.json"
        if jpath.is_file() and jpath.stat().st_size > 10:
            is_stale = False
            jmtime = jpath.stat().st_mtime
            if nom_path and nom_path.is_file() and jmtime < (nom_path.stat().st_mtime - 0.01):
                is_stale = True
            if mcmc_path and mcmc_path.is_file() and jmtime < (mcmc_path.stat().st_mtime - 0.01):
                is_stale = True
            if unnorm_path and unnorm_path.is_file() and jmtime < (unnorm_path.stat().st_mtime - 0.01):
                is_stale = True

            if not is_stale:
                try:
                    with open(jpath, "r") as jf:
                        jdata = json.load(jf)
                    raw_masses = jdata.get("masses", [])
                    raw_nominal = jdata.get("nominal", {})
                    raw_mcmc = jdata.get("mcmc", {})
                    raw_params = jdata.get("mcmc_params", {})
                    if mcmc_path and mcmc_path.is_file() and (not raw_mcmc or not raw_params):
                        is_stale = True
                    elif raw_masses and raw_nominal:
                        nom_p: Dict[str, float] = {}
                        nom_errs: Dict[str, float] = {
                            pk: float(pv) for pk, pv in jdata.get("nominal_param_errors", {}).items()
                        }
                        for pk, pv in jdata.get("nominal_params", {}).items():
                            if pk.endswith("_error"):
                                nom_errs.setdefault(pk[:-6], float(pv))
                            else:
                                nom_p[pk] = float(pv)

                        if not nom_errs and nom_path and nom_path.is_file() and nom_path.stat().st_size > 100:
                            root_vals, root_errs = _read_nominal_params_from_root(nom_path)
                            for rk, rv in root_vals.items():
                                nom_p.setdefault(rk, rv)
                            for rk, re_val in root_errs.items():
                                nom_errs.setdefault(rk, re_val)

                        res: Dict[str, Any] = {
                            "masses": [float(m) for m in raw_masses],
                            "nominal": {k: [float(v) for v in vals] for k, vals in raw_nominal.items()},
                            "mcmc": {},
                            "mcmc_weights": {},
                            "mcmc_params": {pk: [float(v) for v in pvals] for pk, pvals in raw_params.items()},
                            "nominal_params": nom_p,
                            "nominal_param_errors": nom_errs,
                            "h0_data": jdata.get("h0_data"),
                        }
                        for k, m_dict in raw_mcmc.items():
                            res["mcmc"][k] = {}
                            for m_str, vals in m_dict.items():
                                try:
                                    res["mcmc"][k][round(float(m_str), 4)] = [float(v) for v in vals]
                                except ValueError:
                                    pass
                        for k, w_dict in jdata.get("mcmc_weights", {}).items():
                            res["mcmc_weights"][k] = {}
                            for m_str, vals in w_dict.items():
                                try:
                                    res["mcmc_weights"][k][round(float(m_str), 4)] = [float(v) for v in vals]
                                except ValueError:
                                    pass
                        if not res["h0_data"] and unnorm_path:
                            cache_comp = (rundir.root / ".cache" / "moment_comparison_data.json") if rundir else None
                            res["h0_data"] = load_h0_moment_data(unnorm_path, cache_comp)
                        return res
                except Exception:
                    pass

    # 3. Extract with uproot or ExportAmplitudes.C
    masses: List[float] = []
    nominal_amps: Dict[str, List[float]] = {}
    mcmc_samples: Dict[str, Dict[float, List[float]]] = {}
    mcmc_weights: Dict[str, Dict[float, List[float]]] = {}
    mcmc_params: Dict[str, List[float]] = {}
    nominal_params: Dict[str, float] = {}
    nominal_param_errors: Dict[str, float] = {}
    h0_data: Optional[Dict[str, List[float]]] = None

    loaded = False
    if nom_path and nom_path.is_file():
        try:
            import uproot
            with uproot.open(nom_path) as f:
                tname = "result" if "result" in f else ("central_result" if "central_result" in f else None)
                if tname:
                    tree = f[tname]
                    masses = list(tree["mass_bin"].array())
                    for k in tree.keys():
                        if k != "mass_bin" and not k.startswith("H_"):
                            nominal_amps[k] = list(tree[k].array())
                    loaded = True
            nominal_params, nominal_param_errors = _read_nominal_params_from_root(nom_path)

            if mcmc_path and mcmc_path.is_file():
                with uproot.open(mcmc_path) as mf:
                    conv_entries = set()
                    has_conv_tree = False
                    for ctree_name in ("mcmc_params", "mcmc_chi2"):
                        if ctree_name in mf:
                            has_conv_tree = True
                            ctree = mf[ctree_name]
                            ckeys = ctree.keys()
                            statuses = list(ctree["Status"].array()) if "Status" in ckeys else None
                            valids = list(ctree["isValid"].array()) if "isValid" in ckeys else None
                            at_limits = list(ctree["atLimit"].array()) if "atLimit" in ckeys else None
                            m_entries = list(ctree["mcmc_entry"].array()) if "mcmc_entry" in ckeys else None
                            n_tot = len(statuses) if statuses is not None else (len(valids) if valids is not None else 0)
                            for idx in range(n_tot):
                                st_ok = (statuses[idx] == 0) if statuses is not None else True
                                v_ok = bool(valids[idx]) if valids is not None else True
                                lim_ok = (not bool(at_limits[idx])) if at_limits is not None else True
                                if st_ok and v_ok and lim_ok:
                                    entry_id = m_entries[idx] if m_entries is not None else idx
                                    conv_entries.add(entry_id)
                            break

                    mtname = "mcmc_result" if "mcmc_result" in mf else ("result" if "result" in mf else None)
                    if mtname:
                        mtree = mf[mtname]
                        mbins = list(mtree["mass_bin"].array())
                        m_result_entries = list(mtree["mcmc_entry"].array()) if "mcmc_entry" in mtree else None
                        m_weights = list(mtree["weight"].array()) if "weight" in mtree else None
                        n_mass_bins = len(masses) if masses else 1
                        for k in nominal_amps.keys():
                            if k in mtree:
                                arr = list(mtree[k].array())
                                mcmc_samples[k] = {}
                                mcmc_weights[k] = {}
                                for i, (m, v) in enumerate(zip(mbins, arr)):
                                    if has_conv_tree:
                                        entry_id = m_result_entries[i] if m_result_entries is not None else (i // n_mass_bins)
                                        if entry_id not in conv_entries:
                                            continue
                                    m_round = round(float(m), 4)
                                    mcmc_samples[k].setdefault(m_round, []).append(float(v))
                                    if m_weights is not None:
                                        mcmc_weights[k].setdefault(m_round, []).append(float(m_weights[i]))
                    if "mcmc_params" in mf:
                        ptree = mf["mcmc_params"]
                        p_entries = list(ptree["mcmc_entry"].array()) if "mcmc_entry" in ptree else None
                        for pk in ptree.keys():
                            if pk not in ("mcmc_entry", "Status", "isValid", "isValidError", "atLimit"):
                                vals = list(ptree[pk].array())
                                mcmc_params[pk] = [
                                    float(v) for i, v in enumerate(vals)
                                    if (not has_conv_tree or (p_entries[i] if p_entries is not None else i) in conv_entries)
                                ]
        except Exception:
            pass

    if not loaded:
        # Fallback to ExportAmplitudes.C macro
        try:
            repo_root = Path(os.environ.get("MOMS2AMPS", str(Path(__file__).resolve().parent.parent.parent)))
            export_macro = repo_root / "macros" / "ExportAmplitudes.C"
            target_cache_dir.mkdir(parents=True, exist_ok=True)
            export_json = target_cache_dir / "amplitudes_data.json"

            if export_macro.is_file() and nom_path and nom_path.is_file():
                root_bin = get_root_binary(config)
                mcmc_arg = str(Path(mcmc_path).resolve()) if (mcmc_path and Path(mcmc_path).is_file()) else ""
                unnorm_arg = str(Path(unnorm_path).resolve()) if (unnorm_path and Path(unnorm_path).is_file()) else ""
                cmd = [
                    root_bin,
                    "-l",
                    "-b",
                    "-q",
                    f'{export_macro.resolve()}("{Path(nom_path).resolve()}", "{mcmc_arg}", "{export_json.resolve()}", "{unnorm_arg}")'
                ]
                env = get_root_subprocess_env(config)
                run_res = subprocess.run(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
                if run_res.returncode == 0 and export_json.is_file():
                    with open(export_json, "r") as jf:
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
                        raw_weights = jdata.get("mcmc_weights", {})
                        for k, w_dict in raw_weights.items():
                            mcmc_weights[k] = {}
                            for m_str, vals in w_dict.items():
                                try:
                                    mcmc_weights[k][round(float(m_str), 4)] = [float(v) for v in vals]
                                except ValueError:
                                    pass
                        raw_params = jdata.get("mcmc_params", {})
                        for pk, pvals in raw_params.items():
                            mcmc_params[pk] = [float(v) for v in pvals]
                        raw_nom_params = jdata.get("nominal_params", {})
                        for pk, pv in raw_nom_params.items():
                            try:
                                if pk.endswith("_error"):
                                    nominal_param_errors[pk[:-6]] = float(pv)
                                else:
                                    nominal_params[pk] = float(pv)
                            except Exception:
                                pass
                        if "h0_data" in jdata and jdata["h0_data"].get("mass") and jdata["h0_data"].get("val"):
                            h0_data = jdata["h0_data"]
                        loaded = True
        except Exception:
            pass

    if loaded and masses and nominal_amps:
        if h0_data is None and unnorm_path:
            cache_comp = (rundir.root / ".cache" / "moment_comparison_data.json") if rundir else None
            h0_data = load_h0_moment_data(unnorm_path, cache_comp)

        # Write to target_cache_dir / "amplitudes_data.json"
        target_cache_dir.mkdir(parents=True, exist_ok=True)
        try:
            cache_payload = {
                "masses": masses,
                "nominal": nominal_amps,
                "mcmc": {k: {str(m): vals for m, vals in m_dict.items()} for k, m_dict in mcmc_samples.items()},
                "mcmc_weights": {k: {str(m): vals for m, vals in w_dict.items()} for k, w_dict in mcmc_weights.items()},
                "mcmc_params": mcmc_params,
                "nominal_params": nominal_params,
                "nominal_param_errors": nominal_param_errors,
                "h0_data": h0_data,
            }
            with open(target_cache_dir / "amplitudes_data.json", "w") as jf:
                json.dump(cache_payload, jf)
        except Exception:
            pass

        return {
            "masses": masses,
            "nominal": nominal_amps,
            "mcmc": mcmc_samples,
            "mcmc_weights": mcmc_weights,
            "mcmc_params": mcmc_params,
            "nominal_params": nominal_params,
            "nominal_param_errors": nominal_param_errors,
            "h0_data": h0_data,
        }

    return None


def plot_amplitudes_and_bands(
    nominal_root: Path | str,
    mcmc_root: Optional[Path | str],
    output_dir: Path | str,
    config: FitConfig,
    formats: Tuple[str, ...] = ("png", "pdf"),
    only_distributions: bool = False,
    only_amplitudes: bool = False,
    style: Optional[str] = None,
    unnormalized_moments_file: Optional[Path | str] = None,
    rundir: Optional[RunDirectory] = None,
) -> List[Path]:
    """
    Plot fitted amplitude magnitudes and phases vs mass with shaded MCMC 68%/95% credible intervals,
    individual wave intensities/phases, wave overview panels, and 1D posterior distributions.
    Overlays H_0_0_0 (acceptance corrected yield) on amplitude intensities when available.
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

    effective_style = style or (getattr(config.plotting, "style", "default") if getattr(config, "plotting", None) else "default")
    apply_plot_style(effective_style)

    # Load amplitude data from cache, uproot, or ExportAmplitudes.C
    amp_data = load_or_export_amplitude_data(
        config=config,
        rundir=rundir,
        nominal_root=nominal_root,
        mcmc_root=mcmc_root,
        unnormalized_moments_file=unnormalized_moments_file,
        cache_dir=out_dir / ".cache",
    )
    if not amp_data or not amp_data.get("masses"):
        print(f"[WARNING] Could not read amplitude branches from {nominal_root}")
        return generated

    masses = amp_data["masses"]
    nominal_amps = amp_data.get("nominal", {})
    mcmc_samples = amp_data.get("mcmc", {})
    mcmc_weights = amp_data.get("mcmc_weights", {})
    mcmc_params = amp_data.get("mcmc_params", {})
    h0_data = amp_data.get("h0_data")

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

    magnitudes = sorted([k for k in nominal_amps.keys() if (k.startswith("a_") or k.startswith("b_")) and is_amp_active(k)], key=sort_amp_key)
    phases = sorted([k for k in nominal_amps.keys() if (k.startswith("aphi_") or k.startswith("bphi_")) and is_amp_active(k)], key=sort_amp_key)
    nominal_params = amp_data.get("nominal_params", {})
    ref_info = resolve_reference_phase_info(
        config=config,
        nominal_amps=nominal_amps,
        active_phases=phases,
        masses=masses,
        nominal_params=nominal_params,
        mcmc_params=mcmc_params,
        mcmc_samples=mcmc_samples,
    )
    anchors = ref_info["anchors"]

    fit_step = 0.04
    if hasattr(config, "dataset") and getattr(config.dataset, "mass_bins", None):
        fit_step = getattr(config.dataset.mass_bins, "step", 0.04) or 0.04
    elif len(masses) > 1:
        fit_step = float(np.median(np.diff(masses)))

    show_phase_ambiguity = getattr(config.plotting, "show_phase_ambiguity", True) if hasattr(config, "plotting") else True
    density_bands = getattr(config.plotting, "posterior_density_bands", True) if hasattr(config, "plotting") else True

    # --- 1. Amplitude Magnitudes & Intensities Plot (Overlaid on H_0_0_0) with Phases in Bottom Panel ---
    if not only_distributions and magnitudes:
        has_phases_panel = bool(phases)
        if has_phases_panel:
            fig, (ax, ax_phase) = plt.subplots(
                2, 1, figsize=(8.5, 7.5), dpi=300, sharex=True,
                gridspec_kw={"height_ratios": [2.2, 1.2], "hspace": 0.08, "top": 0.93, "bottom": 0.09, "left": 0.11, "right": 0.96}
            )
        else:
            fig, ax = plt.subplots(
                figsize=(8.5, 5.5), dpi=300,
                gridspec_kw={"top": 0.92, "bottom": 0.11, "left": 0.11, "right": 0.96}
            )
            ax_phase = None

        # Build consistent color palette mapped to wave identity (L, refl, M)
        all_wave_keys = []
        for k in sorted(magnitudes):
            p = get_amplitude_symbol_parts(k)
            ident = (p[0], p[1], p[2]) if p else k
            if ident not in all_wave_keys:
                all_wave_keys.append(ident)
        for k in sorted(phases):
            p = get_amplitude_symbol_parts(k)
            ident = (p[0], p[1], p[2]) if p else k
            if ident not in all_wave_keys:
                all_wave_keys.append(ident)

        wave_linestyle_map = {
            ident: PHYSICS_LINESTYLES[i % len(PHYSICS_LINESTYLES)]
            for i, ident in enumerate(all_wave_keys)
        }

        def get_color(key: str) -> str:
            return "black"

        def get_linestyle(key: str) -> Any:
            p = get_amplitude_symbol_parts(key)
            ident = (p[0], p[1], p[2]) if p else key
            return wave_linestyle_map.get(ident, "--")

        # 1a. Overlay H_0_0_0 (acceptance corrected yield) if available
        has_h0_overlay = False
        if h0_data and h0_data.get("mass") and h0_data.get("val"):
            m_h0 = h0_data["mass"]
            y_h0 = h0_data["val"]
            e_h0 = h0_data.get("err", [0.0] * len(y_h0))
            ax.errorbar(
                m_h0, y_h0, yerr=e_h0,
                fmt="o", color="black", markersize=3.2,
                capsize=1.2, elinewidth=0.9, alpha=0.85,
                label=r"$H_0(0,0)$ (Acceptance Corrected Yield)",
                zorder=2,
            )
            has_h0_overlay = True

        # 1b. Total fitted intensity \sum |V_\alpha|^2
        tot_nom = [sum(nominal_amps[k][i]**2 for k in magnitudes) for i in range(len(masses))]
        if len(masses) >= 4:
            m_tot_s, y_tot_s = smooth_curve(masses, tot_nom, num_points=200)
            ax.plot(m_tot_s, y_tot_s, "-", color="black", linewidth=2.2, alpha=0.95, label=r"Total $\sum |V|^2$", zorder=5)
        else:
            ax.plot(masses, tot_nom, "-", color="black", linewidth=2.2, alpha=0.95, label=r"Total $\sum |V|^2$", zorder=5)

        # 1c. Individual wave squared intensities |V_\alpha|^2
        for idx, k in enumerate(sorted(magnitudes)):
            y_nom = nominal_amps[k]
            y_sqr = [v**2 for v in y_nom]
            label = format_amplitude_squared_label(k)
            is_mass_dep = is_amplitude_mass_dependent(k, config)
            c = "black"
            ls = get_linestyle(k)

            if is_mass_dep:
                m_smooth, y_smooth = smooth_curve(masses, y_sqr, num_points=200)
                ax.plot(m_smooth, y_smooth, linestyle=ls, linewidth=1.8, color=c, label=label, zorder=4)

                if k in mcmc_samples:
                    if density_bands:
                        plot_posterior_density_band(
                            ax=ax, masses=masses, mcmc_dict=mcmc_samples[k],
                            is_mass_dep=True, color=c, squared=True, n_levels=20,
                            base_alpha=0.18, max_alpha=0.55, fit_step=fit_step, zorder=3,
                            mcmc_weights=mcmc_weights.get(k)
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True, mcmc_weights=mcmc_weights.get(k))
                        if len(masses) >= 4:
                            m_s, lo_smooth = smooth_curve(masses, lo_68, num_points=200)
                            m_s, hi_smooth = smooth_curve(masses, hi_68, num_points=200)
                            ax.fill_between(m_s, lo_smooth, hi_smooth, color=c, alpha=0.30, zorder=3)
                        else:
                            ax.fill_between(masses, lo_68, hi_68, color=c, alpha=0.30, zorder=3)
            else:
                edges = compute_bin_edges(masses, default_step=fit_step)
                if hasattr(ax, "stairs"):
                    ax.stairs(y_sqr, edges, linestyle=ls, linewidth=1.8, color=c, label=label, zorder=4)
                else:
                    ax.step(masses, y_sqr, where="mid", linestyle=ls, linewidth=1.8, color=c, label=label, zorder=4)

                if k in mcmc_samples:
                    if density_bands:
                        plot_posterior_density_band(
                            ax=ax, masses=masses, mcmc_dict=mcmc_samples[k],
                            is_mass_dep=False, color=c, squared=True, n_levels=20,
                            base_alpha=0.18, max_alpha=0.55, fit_step=fit_step, zorder=3,
                            mcmc_weights=mcmc_weights.get(k)
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True, mcmc_weights=mcmc_weights.get(k))
                        x_step = np.repeat(edges, 2)[1:-1]
                        lo_step = np.repeat(lo_68, 2)
                        hi_step = np.repeat(hi_68, 2)
                        ax.fill_between(x_step, lo_step, hi_step, color=c, alpha=0.30, zorder=3)

        ax.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
        ylabel_str = f"Intensity / ({fit_step*1000:.0f} MeV/$c^2$)" if fit_step else "Intensity"
        if has_h0_overlay:
            ax.set_ylabel(ylabel_str, fontsize=11.5)
            ax.set_title(f"{config.project.name} - Fitted Amplitudes vs $H_0(0,0)$", fontsize=12.5, pad=8)
        else:
            ax.set_ylabel(ylabel_str, fontsize=11.5)
            ax.set_title(f"{config.project.name} - Fitted Amplitude Intensities", fontsize=12.5, pad=8)
        ax.legend(frameon=True, fontsize=8.5, ncol=2, loc="best")
        if should_show_grid(config):
            ax.grid(True, linestyle=":", alpha=0.45)

        # 1d. Bottom panel: Phase as function of mass (relative to reference amplitude)
        if ax_phase is not None:
            drew_phase = False
            for idx, k in enumerate(sorted(phases)):
                if k in ref_info["omit_phases"]:
                    continue

                refl_char = "a" if k.startswith("aphi_") else "b"
                anchor_nom = ref_info["nominal_ref"].get(refl_char)
                anchor_mcmc = ref_info["mcmc_ref"].get(refl_char)
                y_rel, mcmc_rel = compute_relative_phase_data(
                    masses, nominal_amps[k], anchor_nom, mcmc_samples.get(k), anchor_mcmc
                )

                ref_amp = ref_info.get("ref_names", {}).get(refl_char, ref_info["ref_name"])
                label = format_relative_phase_label(k, ref_amp, plus_minus=True)
                is_mass_dep = is_amplitude_mass_dependent(k, config)
                c = "black"
                ls = get_linestyle(k)
                neg_ls = ":" if ls == "--" else "--"

                plot_phase_with_ambiguity(
                    ax=ax_phase,
                    masses=masses,
                    y_nom=y_rel,
                    is_mass_dep=is_mass_dep,
                    color=c,
                    label=label,
                    mcmc_dict=mcmc_rel,
                    fit_step=fit_step,
                    show_ambiguity=show_phase_ambiguity,
                    alpha_band=0.30,
                    base_alpha=0.18,
                    max_alpha=0.55,
                    linestyle_pos=ls,
                    linestyle_neg=neg_ls,
                    density_band=density_bands,
                    mcmc_weights=mcmc_weights.get(k),
                )
                drew_phase = True

            configure_phase_axis(ax_phase)
            ax_phase.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
            ax_phase.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11.5)
            ref_text = format_reference_partial_wave_label(ref_info, phases=phases)
            ax_phase.set_ylabel(rf"Relative Phase $\Delta\Phi$ (Ref Partial-Wave: {ref_text}) [rad]", fontsize=10.5)
            if drew_phase:
                ax_phase.legend(
                    frameon=True,
                    fontsize=8,
                    ncol=min(3, max(1, len(phases) - 1)),
                    loc="best",
                    title=f"Ref Partial-Wave: {ref_text}",
                    title_fontsize=8.5,
                )
            if should_show_grid(config):
                ax_phase.grid(True, linestyle=":", alpha=0.45)

            # Suppress x-tick labels and x-axis label on top panel
            ax.tick_params(labelbottom=False)
            ax.set_xlabel("")
        else:
            ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11.5)

        for fmt in formats:
            p = out_dir / f"amplitudes_magnitude_vs_mass.{fmt}"
            fig.savefig(p, bbox_inches="tight")
            generated.append(p)
            p_int = out_dir / f"amplitudes_intensity_vs_mass.{fmt}"
            fig.savefig(p_int, bbox_inches="tight")
            generated.append(p_int)
        plt.close(fig)

    # --- 2. Relative Phases Plot ---
    if not only_distributions and phases:
        fig, ax = plt.subplots(
            figsize=(8.5, 5.5), dpi=300,
            gridspec_kw={"top": 0.92, "bottom": 0.11, "left": 0.11, "right": 0.96}
        )
        drew_phase = False

        for idx, k in enumerate(sorted(phases)):
            if k in ref_info["omit_phases"]:
                continue

            refl_char = "a" if k.startswith("aphi_") else "b"
            anchor_nom = ref_info["nominal_ref"].get(refl_char)
            anchor_mcmc = ref_info["mcmc_ref"].get(refl_char)
            y_rel, mcmc_rel = compute_relative_phase_data(
                masses, nominal_amps[k], anchor_nom, mcmc_samples.get(k), anchor_mcmc
            )

            ref_amp = ref_info.get("ref_names", {}).get(refl_char, ref_info["ref_name"])
            label = format_relative_phase_label(k, ref_amp, plus_minus=True)
            is_mass_dep = is_amplitude_mass_dependent(k, config)
            c = "black"
            ls = get_linestyle(k)
            neg_ls = ":" if ls == "--" else "--"

            plot_phase_with_ambiguity(
                ax=ax,
                masses=masses,
                y_nom=y_rel,
                is_mass_dep=is_mass_dep,
                color=c,
                label=label,
                mcmc_dict=mcmc_rel,
                fit_step=fit_step,
                show_ambiguity=show_phase_ambiguity,
                alpha_band=0.30,
                base_alpha=0.18,
                max_alpha=0.55,
                linestyle_pos=ls,
                linestyle_neg=neg_ls,
                density_band=density_bands,
                mcmc_weights=mcmc_weights.get(k),
            )
            drew_phase = True

        configure_phase_axis(ax)
        ax.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
        ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11.5)
        ref_text = format_reference_partial_wave_label(ref_info, phases=phases)
        ax.set_ylabel(rf"Relative Phase $\Delta\Phi$ (Reference Partial-Wave: {ref_text}) [rad]", fontsize=11.0)
        ax.set_title(f"{config.project.name} - Fitted Relative Phases $\\Delta\\Phi$ (Reference Partial-Wave: {ref_text}, $\\pm$ Sign Ambiguity)", fontsize=11.5, pad=8)
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
        if drew_phase:
            ax.legend(frameon=True, fontsize=9.0, ncol=2, title=f"Reference Partial-Wave: {ref_text}", title_fontsize=9.0)
        if should_show_grid(config):
            ax.grid(True, linestyle=":", alpha=0.45)

        for fmt in formats:
            p = out_dir / f"amplitudes_phase_vs_mass.{fmt}"
            fig.savefig(p, bbox_inches="tight")
            generated.append(p)
        plt.close(fig)

    # Subdirectories for organized outputs
    amp_dir = out_dir / "amplitudes"
    indiv_dir = amp_dir / "individual"
    dist_dir = amp_dir / "distributions"
    amp_dir.mkdir(parents=True, exist_ok=True)
    indiv_dir.mkdir(parents=True, exist_ok=True)
    dist_dir.mkdir(parents=True, exist_ok=True)

    if not only_distributions:
        # --- 3. Individual Amplitude Magnitude & Squared Intensity Plots ---
        for k in sorted(magnitudes):
            y_nom = nominal_amps[k]
            label = format_amplitude_label(k)
            sqr_label = format_amplitude_squared_label(k)
            is_mass_dep = is_amplitude_mass_dependent(k, config)
            c_wave = "black"

            # 3a. Magnitude vs Mass
            fig, ax = plt.subplots(
                figsize=(7, 4.8), dpi=300,
                gridspec_kw={"top": 0.92, "bottom": 0.12, "left": 0.12, "right": 0.96}
            )
            if h0_data and h0_data.get("mass") and h0_data.get("val"):
                m_h0 = h0_data["mass"]
                y_sqrt = [math.sqrt(max(0.0, v)) for v in h0_data["val"]]
                e_sqrt = [e / (2.0 * math.sqrt(max(1e-6, v))) for v, e in zip(h0_data["val"], h0_data.get("err", [0.0] * len(y_sqrt)))]
                ax.errorbar(
                    m_h0, y_sqrt, yerr=e_sqrt,
                    fmt="o", color="black", markersize=2.8, alpha=0.85,
                    capsize=1.2, elinewidth=0.9,
                    label=r"$\sqrt{H_0(0,0)}$ (Total)",
                    zorder=2,
                )
            if is_mass_dep:
                m_s, y_s = smooth_curve(masses, y_nom, num_points=200)
                ax.plot(m_s, y_s, "-", linewidth=2.0, color=c_wave, label=f"Nominal {label}", zorder=4)
                if k in mcmc_samples:
                    if density_bands:
                        plot_posterior_density_band(
                            ax=ax, masses=masses, mcmc_dict=mcmc_samples[k],
                            is_mass_dep=True, color=c_wave, squared=False, n_levels=20,
                            base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3,
                            mcmc_weights=mcmc_weights.get(k)
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), mcmc_weights=mcmc_weights.get(k))
                        lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5), mcmc_weights=mcmc_weights.get(k))
                        m_s95, lo_95s = smooth_curve(masses, lo_95, num_points=200)
                        m_s95, hi_95s = smooth_curve(masses, hi_95, num_points=200)
                        m_s68, lo_68s = smooth_curve(masses, lo_68, num_points=200)
                        m_s68, hi_68s = smooth_curve(masses, hi_68, num_points=200)
                        ax.fill_between(m_s95, lo_95s, hi_95s, color=c_wave, alpha=0.16, label="MCMC 95% CI", zorder=2)
                        ax.fill_between(m_s68, lo_68s, hi_68s, color=c_wave, alpha=0.32, label="MCMC 68% CI", zorder=3)
            else:
                edges = compute_bin_edges(masses, default_step=fit_step)
                if hasattr(ax, "stairs"):
                    ax.stairs(y_nom, edges, linewidth=2.0, color=c_wave, label=f"Nominal {label}", zorder=4)
                else:
                    ax.step(masses, y_nom, where="mid", linewidth=2.0, color=c_wave, label=f"Nominal {label}", zorder=4)
                if k in mcmc_samples:
                    if density_bands:
                        plot_posterior_density_band(
                            ax=ax, masses=masses, mcmc_dict=mcmc_samples[k],
                            is_mass_dep=False, color=c_wave, squared=False, n_levels=20,
                            base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3,
                            mcmc_weights=mcmc_weights.get(k)
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), mcmc_weights=mcmc_weights.get(k))
                        lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5), mcmc_weights=mcmc_weights.get(k))
                        x_step = np.repeat(edges, 2)[1:-1]
                        lo_95s = np.repeat(lo_95, 2)
                        hi_95s = np.repeat(hi_95, 2)
                        lo_68s = np.repeat(lo_68, 2)
                        hi_68s = np.repeat(hi_68, 2)
                        ax.fill_between(x_step, lo_95s, hi_95s, color=c_wave, alpha=0.16, label="MCMC 95% CI", zorder=2)
                        ax.fill_between(x_step, lo_68s, hi_68s, color=c_wave, alpha=0.32, label="MCMC 68% CI", zorder=3)

            ax.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
            ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11)
            ax.set_ylabel(f"Magnitude {label}", fontsize=11)
            ax.set_title(f"{config.project.name} - Amplitude {label}", fontsize=12, pad=8)
            if density_bands and k in mcmc_samples:
                from matplotlib.patches import Patch
                handles, labels = ax.get_legend_handles_labels()
                handles.append(Patch(color=c_wave, alpha=0.35))
                labels.append("MCMC Posterior Density")
                ax.legend(handles=handles, labels=labels, frameon=True, fontsize=9, loc="best")
            else:
                ax.legend(frameon=True, fontsize=9, loc="best")
            if should_show_grid(config):
                ax.grid(True, linestyle=":", alpha=0.45)
            for fmt in formats:
                p_indiv = indiv_dir / f"{k}_magnitude.{fmt}"
                fig.savefig(p_indiv, bbox_inches="tight")
                generated.append(p_indiv)
            plt.close(fig)

            # 3b. Squared Intensity vs Mass
            fig, ax = plt.subplots(
                figsize=(7, 4.8), dpi=300,
                gridspec_kw={"top": 0.92, "bottom": 0.12, "left": 0.12, "right": 0.96}
            )
            y_nom_sqr = [v**2 for v in y_nom]
            if h0_data and h0_data.get("mass") and h0_data.get("val"):
                m_h0 = h0_data["mass"]
                y_h0 = h0_data["val"]
                e_h0 = h0_data.get("err", [0.0] * len(y_h0))
                ax.errorbar(
                    m_h0, y_h0, yerr=e_h0,
                    fmt="o", color="black", markersize=2.8, alpha=0.85,
                    capsize=1.2, elinewidth=0.9,
                    label=r"$H_0(0,0)$ (Total Yield)",
                    zorder=2,
                )
            if is_mass_dep:
                m_s, y_s = smooth_curve(masses, y_nom_sqr, num_points=200)
                ax.plot(m_s, y_s, "-", linewidth=2.0, color=c_wave, label=f"Nominal {sqr_label}", zorder=4)
                if k in mcmc_samples:
                    if density_bands:
                        plot_posterior_density_band(
                            ax=ax, masses=masses, mcmc_dict=mcmc_samples[k],
                            is_mass_dep=True, color=c_wave, squared=True, n_levels=20,
                            base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3,
                            mcmc_weights=mcmc_weights.get(k)
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True, mcmc_weights=mcmc_weights.get(k))
                        lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5), squared=True, mcmc_weights=mcmc_weights.get(k))
                        m_s95, lo_95s = smooth_curve(masses, lo_95, num_points=200)
                        m_s95, hi_95s = smooth_curve(masses, hi_95, num_points=200)
                        m_s68, lo_68s = smooth_curve(masses, lo_68, num_points=200)
                        m_s68, hi_68s = smooth_curve(masses, hi_68, num_points=200)
                        ax.fill_between(m_s95, lo_95s, hi_95s, color=c_wave, alpha=0.16, label="MCMC 95% CI", zorder=2)
                        ax.fill_between(m_s68, lo_68s, hi_68s, color=c_wave, alpha=0.32, label="MCMC 68% CI", zorder=3)
            else:
                edges = compute_bin_edges(masses, default_step=fit_step)
                if hasattr(ax, "stairs"):
                    ax.stairs(y_nom_sqr, edges, linewidth=2.0, color=c_wave, label=f"Nominal {sqr_label}", zorder=4)
                else:
                    ax.step(masses, y_nom_sqr, where="mid", linewidth=2.0, color=c_wave, label=f"Nominal {sqr_label}", zorder=4)
                if k in mcmc_samples:
                    if density_bands:
                        plot_posterior_density_band(
                            ax=ax, masses=masses, mcmc_dict=mcmc_samples[k],
                            is_mass_dep=False, color=c_wave, squared=True, n_levels=20,
                            base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3,
                            mcmc_weights=mcmc_weights.get(k)
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True, mcmc_weights=mcmc_weights.get(k))
                        lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5), squared=True, mcmc_weights=mcmc_weights.get(k))
                        x_step = np.repeat(edges, 2)[1:-1]
                        lo_95s = np.repeat(lo_95, 2)
                        hi_95s = np.repeat(hi_95, 2)
                        lo_68s = np.repeat(lo_68, 2)
                        hi_68s = np.repeat(hi_68, 2)
                        ax.fill_between(x_step, lo_95s, hi_95s, color=c_wave, alpha=0.16, label="MCMC 95% CI", zorder=2)
                        ax.fill_between(x_step, lo_68s, hi_68s, color=c_wave, alpha=0.32, label="MCMC 68% CI", zorder=3)

            ax.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
            ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11)
            ax.set_ylabel(f"Intensity {sqr_label}", fontsize=11)
            ax.set_title(f"{config.project.name} - Intensity {sqr_label}", fontsize=12, pad=8)
            if density_bands and k in mcmc_samples:
                from matplotlib.patches import Patch
                handles, labels = ax.get_legend_handles_labels()
                handles.append(Patch(color=c_wave, alpha=0.35))
                labels.append("MCMC Posterior Density")
                ax.legend(handles=handles, labels=labels, frameon=True, fontsize=9, loc="best")
            else:
                ax.legend(frameon=True, fontsize=9, loc="best")
            if should_show_grid(config):
                ax.grid(True, linestyle=":", alpha=0.45)
            for fmt in formats:
                p_indiv_sqr = indiv_dir / f"{k}_squared.{fmt}"
                fig.savefig(p_indiv_sqr, bbox_inches="tight")
                generated.append(p_indiv_sqr)
            plt.close(fig)

            # 3c. Individual Wave Intensity and Phase (2-panel)
            phase_k = None
            if k.startswith("a_"):
                cand = "aphi_" + k[2:]
                if cand in nominal_amps:
                    phase_k = cand
            elif k.startswith("b_"):
                cand = "bphi_" + k[2:]
                if cand in nominal_amps:
                    phase_k = cand

            refl_char = "a" if k.startswith("a_") else "b"
            is_omitted = bool(phase_k and phase_k in ref_info["omit_phases"])

            # Remove stale files if this phase is omitted (e.g. reference wave whose DeltaPhi == 0)
            if is_omitted:
                for fmt in formats:
                    stale_f = indiv_dir / f"{k}_intensity_and_phase.{fmt}"
                    if stale_f.exists():
                        try:
                            stale_f.unlink()
                        except OSError:
                            pass

            if phase_k and not is_omitted:
                fig, (ax_top, ax_bot) = plt.subplots(
                    2, 1, figsize=(7, 6.5), dpi=300, sharex=True,
                    gridspec_kw={"height_ratios": [2.2, 1.2], "hspace": 0.08, "top": 0.93, "bottom": 0.10, "left": 0.12, "right": 0.96}
                )
                if h0_data and h0_data.get("mass") and h0_data.get("val"):
                    m_h0 = h0_data["mass"]
                    y_h0 = h0_data["val"]
                    e_h0 = h0_data.get("err", [0.0] * len(y_h0))
                    ax_top.errorbar(
                        m_h0, y_h0, yerr=e_h0,
                        fmt="o", color="black", markersize=2.8, alpha=0.85,
                        capsize=1.2, elinewidth=0.9,
                        label=r"$H_0(0,0)$ (Total Yield)",
                        zorder=2,
                    )
                if is_mass_dep:
                    m_s, y_s = smooth_curve(masses, y_nom_sqr, num_points=200)
                    ax_top.plot(m_s, y_s, "-", linewidth=2.0, color=c_wave, label=f"Nominal {sqr_label}", zorder=4)
                    if k in mcmc_samples:
                        if density_bands:
                            plot_posterior_density_band(
                                ax=ax_top, masses=masses, mcmc_dict=mcmc_samples[k],
                                is_mass_dep=True, color=c_wave, squared=True, n_levels=20,
                                base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3,
                                mcmc_weights=mcmc_weights.get(k)
                            )
                        else:
                            lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True, mcmc_weights=mcmc_weights.get(k))
                            lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5), squared=True, mcmc_weights=mcmc_weights.get(k))
                            m_s95, lo_95s = smooth_curve(masses, lo_95, num_points=200)
                            m_s95, hi_95s = smooth_curve(masses, hi_95, num_points=200)
                            m_s68, lo_68s = smooth_curve(masses, lo_68, num_points=200)
                            m_s68, hi_68s = smooth_curve(masses, hi_68, num_points=200)
                            ax_top.fill_between(m_s95, lo_95s, hi_95s, color=c_wave, alpha=0.16, label="MCMC 95% CI", zorder=2)
                            ax_top.fill_between(m_s68, lo_68s, hi_68s, color=c_wave, alpha=0.32, label="MCMC 68% CI", zorder=3)
                else:
                    edges = compute_bin_edges(masses, default_step=fit_step)
                    if hasattr(ax_top, "stairs"):
                        ax_top.stairs(y_nom_sqr, edges, linewidth=2.0, color=c_wave, label=f"Nominal {sqr_label}", zorder=4)
                    else:
                        ax_top.step(masses, y_nom_sqr, where="mid", linewidth=2.0, color=c_wave, label=f"Nominal {sqr_label}", zorder=4)
                    if k in mcmc_samples:
                        if density_bands:
                            plot_posterior_density_band(
                                ax=ax_top, masses=masses, mcmc_dict=mcmc_samples[k],
                                is_mass_dep=False, color=c_wave, squared=True, n_levels=20,
                                base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3,
                                mcmc_weights=mcmc_weights.get(k)
                            )
                        else:
                            lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True, mcmc_weights=mcmc_weights.get(k))
                            lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5), squared=True, mcmc_weights=mcmc_weights.get(k))
                            x_step = np.repeat(edges, 2)[1:-1]
                            lo_95s = np.repeat(lo_95, 2)
                            hi_95s = np.repeat(hi_95, 2)
                            lo_68s = np.repeat(lo_68, 2)
                            hi_68s = np.repeat(hi_68, 2)
                            ax_top.fill_between(x_step, lo_95s, hi_95s, color=c_wave, alpha=0.16, label="MCMC 95% CI", zorder=2)
                            ax_top.fill_between(x_step, lo_68s, hi_68s, color=c_wave, alpha=0.32, label="MCMC 68% CI", zorder=3)

                ax_top.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
                ax_top.set_ylabel(f"Intensity {sqr_label}", fontsize=11)
                ax_top.set_title(f"{config.project.name} - Wave {label}", fontsize=12, pad=8)
                if density_bands and k in mcmc_samples:
                    from matplotlib.patches import Patch
                    handles, labels = ax_top.get_legend_handles_labels()
                    handles.append(Patch(color=c_wave, alpha=0.35))
                    labels.append("MCMC Posterior Density")
                    ax_top.legend(handles=handles, labels=labels, frameon=True, fontsize=8, loc="best")
                else:
                    ax_top.legend(frameon=True, fontsize=8, loc="best")
                if should_show_grid(config):
                    ax_top.grid(True, linestyle=":", alpha=0.45)
                ax_top.tick_params(labelbottom=False)

                anchor_nom = ref_info["nominal_ref"].get(refl_char)
                anchor_mcmc = ref_info["mcmc_ref"].get(refl_char)
                y_rel_phi, mcmc_rel_phi = compute_relative_phase_data(
                    masses, nominal_amps[phase_k], anchor_nom, mcmc_samples.get(phase_k), anchor_mcmc
                )
                ref_amp = ref_info.get("ref_names", {}).get(refl_char, ref_info["ref_name"])
                ref_lbl_refl = format_reference_partial_wave_label(ref_info, refl_char=refl_char)
                phi_label = format_relative_phase_label(phase_k, ref_amp)
                phi_sym = format_relative_phase_symbol(phase_k, ref_amp)
                is_phi_mass_dep = is_amplitude_mass_dependent(phase_k, config)

                plot_phase_with_ambiguity(
                    ax=ax_bot,
                    masses=masses,
                    y_nom=y_rel_phi,
                    is_mass_dep=is_phi_mass_dep,
                    color=c_wave,
                    label=f"Nominal +${phi_sym}$",
                    label_neg=f"Sign Ambiguity -${phi_sym}$",
                    mcmc_dict=mcmc_rel_phi,
                    fit_step=fit_step,
                    show_ambiguity=show_phase_ambiguity,
                    alpha_band=0.32,
                    base_alpha=0.20,
                    max_alpha=0.60,
                    density_band=density_bands,
                    mcmc_weights=mcmc_weights.get(phase_k),
                )

                configure_phase_axis(ax_bot)
                ax_bot.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
                ax_bot.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11)
                ax_bot.set_ylabel(f"Relative Phase {phi_label}\n(Ref Partial-Wave: {ref_lbl_refl}) [rad]", fontsize=10.0)
                ax_bot.legend(frameon=True, fontsize=8, loc="best", title=f"Ref Partial-Wave: {ref_lbl_refl}", title_fontsize=8.0)
                if should_show_grid(config):
                    ax_bot.grid(True, linestyle=":", alpha=0.45)

                for fmt in formats:
                    p_indiv_combo = indiv_dir / f"{k}_intensity_and_phase.{fmt}"
                    fig.savefig(p_indiv_combo, bbox_inches="tight")
                    generated.append(p_indiv_combo)
                plt.close(fig)

        # --- 4. Individual Phase Plots ---
        for k in sorted(phases):
            if k in ref_info["omit_phases"]:
                for fmt in formats:
                    stale_f = indiv_dir / f"{k}_phase.{fmt}"
                    if stale_f.exists():
                        try:
                            stale_f.unlink()
                        except OSError:
                            pass
                continue

            refl_char = "a" if k.startswith("aphi_") else "b"
            ref_amp = ref_info.get("ref_names", {}).get(refl_char, ref_info["ref_name"])
            ref_lbl_refl = format_reference_partial_wave_label(ref_info, refl_char=refl_char)

            anchor_nom = ref_info["nominal_ref"].get(refl_char)
            anchor_mcmc = ref_info["mcmc_ref"].get(refl_char)
            y_rel, mcmc_rel = compute_relative_phase_data(
                masses, nominal_amps[k], anchor_nom, mcmc_samples.get(k), anchor_mcmc
            )

            label = format_relative_phase_label(k, ref_amp)
            sym = format_relative_phase_symbol(k, ref_amp)
            is_mass_dep = is_amplitude_mass_dependent(k, config)
            c_wave = "black"

            fig, ax = plt.subplots(
                figsize=(7, 4.8), dpi=300,
                gridspec_kw={"top": 0.92, "bottom": 0.12, "left": 0.12, "right": 0.96}
            )
            plot_phase_with_ambiguity(
                ax=ax,
                masses=masses,
                y_nom=y_rel,
                is_mass_dep=is_mass_dep,
                color=c_wave,
                label=f"Nominal +${sym}$",
                label_neg=f"Sign Ambiguity -${sym}$",
                mcmc_dict=mcmc_rel,
                fit_step=fit_step,
                show_ambiguity=show_phase_ambiguity,
                alpha_band=0.32,
                base_alpha=0.20,
                max_alpha=0.60,
                density_band=density_bands,
                mcmc_weights=mcmc_weights.get(k),
            )

            configure_phase_axis(ax)
            ax.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
            ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11)
            ax.set_ylabel(f"Relative Phase {label}\n(Reference Partial-Wave: {ref_lbl_refl}) [rad]", fontsize=10.5)
            ax.set_title(f"{config.project.name} - Relative Phase {label} (Reference Partial-Wave: {ref_lbl_refl}, $\\pm$ Sign Ambiguity)", fontsize=11.5, pad=8)
            ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
            ax.legend(frameon=True, fontsize=8.5, title=f"Reference Partial-Wave: {ref_lbl_refl}", title_fontsize=8.5)
            if should_show_grid(config):
                ax.grid(True, linestyle=":", alpha=0.45)
            for fmt in formats:
                p_indiv_phi = indiv_dir / f"{k}_phase.{fmt}"
                fig.savefig(p_indiv_phi, bbox_inches="tight")
                generated.append(p_indiv_phi)
            plt.close(fig)

        # --- 5. Grouped Wave Overview Panels (Sharing Y-Axis for Background H_0_0_0 & Phases Below) ---
        def find_phase_key(amp_key: str) -> Optional[str]:
            if amp_key.startswith("a_"):
                cand = "aphi_" + amp_key[2:]
                if cand in nominal_amps:
                    return cand
            elif amp_key.startswith("b_"):
                cand = "bphi_" + amp_key[2:]
                if cand in nominal_amps:
                    return cand
            return None

        def draw_wave_grid(wave_keys: List[str], grid_title: str, out_filename_base: str):
            if not wave_keys:
                return
            n_w = len(wave_keys)
            ncols = min(5, n_w)
            n_blocks = math.ceil(n_w / ncols)
            has_phases = any(find_phase_key(k) is not None for k in wave_keys)

            fig_height = (3.8 * n_blocks if has_phases else 2.6 * n_blocks) + 0.6
            fig = plt.figure(figsize=(3.2 * ncols, fig_height), dpi=300)

            top_margin = 0.86 if n_blocks == 1 else 0.90
            bottom_margin = 0.11 if n_blocks == 1 else 0.07
            block_hspace = 0.38 if has_phases else 0.28

            gs_outer = fig.add_gridspec(
                n_blocks, 1,
                top=top_margin, bottom=bottom_margin, left=0.07, right=0.98,
                hspace=block_hspace if n_blocks > 1 else 0.0
            )

            ref_int_ax = None
            ref_phi_ax = None

            for b in range(n_blocks):
                if has_phases:
                    gs_inner = gs_outer[b, 0].subgridspec(
                        2, ncols, height_ratios=[2.2, 1.2], hspace=0.06, wspace=0.06
                    )
                else:
                    gs_inner = gs_outer[b, 0].subgridspec(
                        1, ncols, wspace=0.06
                    )

                for c in range(ncols):
                    i = b * ncols + c
                    if i < n_w:
                        ax_int = fig.add_subplot(gs_inner[0, c])
                        ax_phi = fig.add_subplot(gs_inner[1, c]) if has_phases else None
                        k = wave_keys[i]
                        y_nom_sqr = [v**2 for v in nominal_amps[k]]
                        sqr_label = format_amplitude_squared_label(k)
                        is_mass_dep = is_amplitude_mass_dependent(k, config)
                        c_wave = "black"

                        # Share y-axis among all intensity subplots so H_0_0_0 background scale is identical
                        if ref_int_ax is None:
                            ref_int_ax = ax_int
                        else:
                            ax_int.sharey(ref_int_ax)

                        if has_phases and ax_phi is not None:
                            # Share x-axis between intensity and phase
                            ax_int.sharex(ax_phi)
                            # Share y-axis among all phase subplots
                            if ref_phi_ax is None:
                                ref_phi_ax = ax_phi
                            else:
                                ax_phi.sharey(ref_phi_ax)

                        # Overlay H_0(0,0) acceptance corrected yield in background
                        if h0_data and h0_data.get("mass") and h0_data.get("val"):
                            ax_int.errorbar(
                                h0_data["mass"], h0_data["val"], yerr=h0_data.get("err", [0.0] * len(h0_data["val"])),
                                fmt="o", color="black", markersize=2.2, alpha=0.85,
                                capsize=1.0, elinewidth=0.8, zorder=2,
                                label=r"$H_0(0,0)$" if (i == 0) else None,
                            )

                        # Wave intensity curve and credible bands
                        if is_mass_dep:
                            m_s, y_s = smooth_curve(masses, y_nom_sqr, num_points=200)
                            ax_int.plot(m_s, y_s, "-", linewidth=1.8, color="black", zorder=4)
                            if k in mcmc_samples:
                                if density_bands:
                                    plot_posterior_density_band(
                                        ax=ax_int, masses=masses, mcmc_dict=mcmc_samples[k],
                                        is_mass_dep=True, color="black", squared=True, n_levels=20,
                                        base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3,
                                        mcmc_weights=mcmc_weights.get(k)
                                    )
                                else:
                                    lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], nominal_amps[k], (16.0, 84.0), squared=True, mcmc_weights=mcmc_weights.get(k))
                                    m_s68, lo_s = smooth_curve(masses, lo_68, num_points=200)
                                    m_s68, hi_s = smooth_curve(masses, hi_68, num_points=200)
                                    ax_int.fill_between(m_s68, lo_s, hi_s, color="black", alpha=0.30, zorder=3)
                        else:
                            edges = compute_bin_edges(masses, default_step=fit_step)
                            if hasattr(ax_int, "stairs"):
                                ax_int.stairs(y_nom_sqr, edges, linewidth=1.8, color="black", zorder=4)
                            else:
                                ax_int.step(masses, y_nom_sqr, where="mid", linewidth=1.8, color="black", zorder=4)
                            if k in mcmc_samples:
                                if density_bands:
                                    plot_posterior_density_band(
                                        ax=ax_int, masses=masses, mcmc_dict=mcmc_samples[k],
                                        is_mass_dep=False, color="black", squared=True, n_levels=20,
                                        base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3,
                                        mcmc_weights=mcmc_weights.get(k)
                                    )
                                else:
                                    lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], nominal_amps[k], (16.0, 84.0), squared=True, mcmc_weights=mcmc_weights.get(k))
                                    x_step = np.repeat(edges, 2)[1:-1]
                                    lo_s = np.repeat(lo_68, 2)
                                    hi_s = np.repeat(hi_68, 2)
                                    ax_int.fill_between(x_step, lo_s, hi_s, color="black", alpha=0.30, zorder=3)

                        ax_int.set_title(sqr_label, fontsize=10.5, pad=4)
                        ax_int.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
                        if should_show_grid(config):
                            ax_int.grid(True, linestyle=":", alpha=0.45)

                        # Y-axis labelling: only leftmost column shows y-axis ticks and title
                        if c == 0:
                            ylabel_str = f"Intensity / ({fit_step*1000:.0f} MeV/$c^2$)" if fit_step else "Intensity"
                            ax_int.set_ylabel(ylabel_str, fontsize=9.5)
                            ax_int.tick_params(labelleft=True)
                            if has_phases and ax_phi is not None:
                                configure_phase_axis(ax_phi)
                                ref_text = format_reference_partial_wave_label(ref_info, phases=nominal_amps.keys())
                                ax_phi.set_ylabel(rf"$\Delta\Phi$ [rad]" + f"\n(Ref: {ref_text})", fontsize=8.5)
                                ax_phi.tick_params(labelleft=True)
                        else:
                            ax_int.tick_params(labelleft=False)
                            ax_int.set_ylabel("")
                            if has_phases and ax_phi is not None:
                                configure_phase_axis(ax_phi)
                                ax_phi.tick_params(labelleft=False)
                                ax_phi.set_ylabel("")

                        # X-axis labelling: bottom of each block
                        if has_phases and ax_phi is not None:
                            ax_int.tick_params(labelbottom=False)
                            ax_int.set_xlabel("")
                            ax_phi.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=9.5)
                            ax_phi.tick_params(labelbottom=True)

                            # Phase plotting in bottom panel (relative to phase anchor / reference resonance)
                            pk = find_phase_key(k)
                            refl_char = "a" if k.startswith("a_") else "b"
                            is_omitted = bool(pk and pk in ref_info["omit_phases"])
                            ref_lbl_refl = format_reference_partial_wave_label(ref_info, refl_char=refl_char)

                            if is_omitted:
                                ax_phi.plot([masses[0] - fit_step*0.6, masses[-1] + fit_step*0.6], [0, 0], "--", color="gray", linewidth=1.0, alpha=0.7)
                                ax_phi.text(
                                    0.5, 0.5, r"$\Delta\phi \equiv 0$" + f"\n(Reference Partial-Wave: {ref_lbl_refl})",
                                    transform=ax_phi.transAxes, ha="center", va="center",
                                    fontsize=8, color="gray", style="italic",
                                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.7),
                                )
                            elif pk:
                                anchor_nom = ref_info["nominal_ref"].get(refl_char)
                                anchor_mcmc = ref_info["mcmc_ref"].get(refl_char)
                                y_rel_phi, mcmc_rel_phi = compute_relative_phase_data(
                                    masses, nominal_amps[pk], anchor_nom, mcmc_samples.get(pk), anchor_mcmc
                                )
                                is_phi_mass_dep = is_amplitude_mass_dependent(pk, config)
                                plot_phase_with_ambiguity(
                                    ax=ax_phi,
                                    masses=masses,
                                    y_nom=y_rel_phi,
                                    is_mass_dep=is_phi_mass_dep,
                                    color="black",
                                    mcmc_dict=mcmc_rel_phi,
                                    fit_step=fit_step,
                                    show_ambiguity=show_phase_ambiguity,
                                    alpha_band=0.30,
                                    base_alpha=0.20,
                                    max_alpha=0.60,
                                    density_band=density_bands,
                                    mcmc_weights=mcmc_weights.get(pk),
                                )
                            else:
                                ax_phi.plot([masses[0] - fit_step*0.6, masses[-1] + fit_step*0.6], [0, 0], "--", color="gray", linewidth=1.0, alpha=0.7)
                                ax_phi.text(
                                    0.5, 0.5, r"$\phi = 0$ (Fixed)",
                                    transform=ax_phi.transAxes, ha="center", va="center",
                                    fontsize=8, color="gray", style="italic",
                                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.7),
                                )

                            ax_phi.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
                            if should_show_grid(config):
                                ax_phi.grid(True, linestyle=":", alpha=0.45)
                        else:
                            ax_int.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=9.5)
                            ax_int.tick_params(labelbottom=True)
                    else:
                        ax_int = fig.add_subplot(gs_inner[0, c])
                        ax_int.set_visible(False)
                        if has_phases:
                            ax_phi = fig.add_subplot(gs_inner[1, c])
                            ax_phi.set_visible(False)

            if ref_phi_ax is not None:
                configure_phase_axis(ref_phi_ax)

            ref_text_grid = format_reference_partial_wave_label(ref_info, phases=nominal_amps.keys())
            if has_phases:
                fig.suptitle(f"{config.project.name} - {grid_title} (Reference Partial-Wave: {ref_text_grid})", fontsize=11.5, y=0.965 if n_blocks == 1 else 0.970)
            else:
                fig.suptitle(f"{config.project.name} - {grid_title}", fontsize=12.5, y=0.965 if n_blocks == 1 else 0.970)
            for fmt in formats:
                p_grid = amp_dir / f"{out_filename_base}.{fmt}"
                fig.savefig(p_grid, bbox_inches="tight")
                generated.append(p_grid)
                if out_filename_base == "overview_all_waves":
                    p_main = out_dir / f"amplitudes_overview_vs_mass.{fmt}"
                    fig.savefig(p_main, bbox_inches="tight")
                    generated.append(p_main)
            plt.close(fig)

        # Collect and plot overview wave grids for each active L (S, P, D, etc.)
        for l_val in sorted(set(extract_amplitude_l(k) for k in magnitudes if extract_amplitude_l(k) is not None)):
            l_sym = L_SPECTROSCOPIC.get(l_val, f"L{l_val}")
            a_l = [k for k in magnitudes if k.startswith(f"a_{l_val}_")]
            b_l = [k for k in magnitudes if k.startswith(f"b_{l_val}_")]
            if l_val == 0:
                s_all = a_l + b_l
                if s_all:
                    draw_wave_grid(s_all, f"{l_sym}-Wave vs Mass", f"overview_{l_sym}_waves")
            else:
                if a_l:
                    draw_wave_grid(a_l, f"{l_sym}-Waves Positive Reflectivity vs Mass", f"overview_{l_sym}_waves_positive")
                if b_l:
                    draw_wave_grid(b_l, f"{l_sym}-Waves Negative Reflectivity vs Mass", f"overview_{l_sym}_waves_negative")
                if a_l or b_l:
                    draw_wave_grid(a_l + b_l, f"{l_sym}-Waves vs Mass", f"overview_{l_sym}_waves")

        draw_wave_grid(magnitudes, "All Fitted Waves Overview vs Mass", "overview_all_waves")

    # --- 6. 1D MCMC Posterior Distributions ---
    plot_dist_cfg = getattr(config.plotting, "distributions", None)
    if plot_dist_cfg is None:
        plot_dist_cfg = getattr(config.plotting, "plot_distributions", True)
    cfg_dist_bins = getattr(config.plotting, "distribution_bins", None) if getattr(config, "plotting", None) else None
    cfg_overview_bins = getattr(config.plotting, "overview_bins", None) if getattr(config, "plotting", None) else None

    if not only_amplitudes and plot_dist_cfg:
        # 6a. Fitted Parameter Posterior Distributions
        if mcmc_params:
            weights = mcmc_params.get("weight", None)
            w_arr = np.array(weights, dtype=float) if weights is not None else None

            for pname, pvals in mcmc_params.items():
                if pname == "weight":
                    continue
                vals = np.array(pvals, dtype=float)
                if len(vals) < 2 or np.std(vals) < 1e-12:
                    continue
                fig, ax = plt.subplots(
                    figsize=(5.5, 4.2), dpi=300,
                    gridspec_kw={"top": 0.91, "bottom": 0.13, "left": 0.13, "right": 0.95}
                )
                if w_arr is not None and len(w_arr) == len(vals):
                    sorter = np.argsort(vals)
                    s_vals = vals[sorter]
                    s_wts = w_arr[sorter]
                    cum_w = np.cumsum(s_wts)
                    tot_w = cum_w[-1]
                    if tot_w > 0:
                        q16, q50, q84 = np.interp([0.16 * tot_w, 0.50 * tot_w, 0.84 * tot_w], cum_w, s_vals)
                    else:
                        q16, q50, q84 = np.percentile(vals, [16, 50, 84])
                    mean_v = np.average(vals, weights=w_arr)
                    std_v = np.sqrt(np.average((vals - mean_v)**2, weights=w_arr))
                else:
                    q16, q50, q84 = np.percentile(vals, [16, 50, 84])
                    mean_v = np.mean(vals)
                    std_v = np.std(vals)

                n_bins = _get_hist_bins(len(vals), user_bins=cfg_dist_bins, is_overview=False)
                if w_arr is not None and len(w_arr) == len(vals):
                    ax.hist(vals, weights=w_arr, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                else:
                    ax.hist(vals, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                ax.axvline(q50, color="#1A2E40", linestyle="--", linewidth=1.5, label=f"Median: {q50:.4f}")
                ax.axvline(q16, color="#2B4C7E", linestyle=":", linewidth=1.2, label=f"68% CI: [{q16:.4f}, {q84:.4f}]")
                ax.axvline(q84, color="#2B4C7E", linestyle=":", linewidth=1.2)
                ax.axvspan(q16, q84, color="#2B4C7E", alpha=0.15)

                clean_name = format_mcmc_parameter_label(pname, short=False)

                ax.set_title(f"MCMC Posterior Distribution: {clean_name}", fontsize=11, pad=6)
                ax.set_xlabel(f"{clean_name}", fontsize=10)
                ax.set_ylabel("Probability Density", fontsize=10)
                ax.legend(frameon=True, fontsize=8, loc="upper right")
                if should_show_grid(config):
                    ax.grid(True, linestyle=":", alpha=0.4)

                safe_name = re.sub(r"[^a-zA-Z0-9_\-\.]", "_", pname)
                for fmt in formats:
                    p_dist = dist_dir / f"{safe_name}_dist.{fmt}"
                    fig.savefig(p_dist, bbox_inches="tight")
                    generated.append(p_dist)
                plt.close(fig)

            # Overview grid of resonance/key parameters (MD_* and chi2)
            key_params = [p for p in mcmc_params.keys() if p.startswith("MD_") or p == "chi2"]
            if key_params:
                ncols = min(4, len(key_params))
                nrows = math.ceil(len(key_params) / ncols)
                fig, axes = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 3.0 * nrows), dpi=300)
                axes_flat = [axes] if len(key_params) == 1 else (axes.flatten() if hasattr(axes, "flatten") else [axes])

                for i, pname in enumerate(key_params):
                    ax = axes_flat[i]
                    vals = np.array(mcmc_params[pname], dtype=float)
                    if w_arr is not None and len(w_arr) == len(vals):
                        sorter = np.argsort(vals)
                        s_vals = vals[sorter]
                        s_wts = w_arr[sorter]
                        cum_w = np.cumsum(s_wts)
                        tot_w = cum_w[-1]
                        if tot_w > 0:
                            q16, q50, q84 = np.interp([0.16 * tot_w, 0.50 * tot_w, 0.84 * tot_w], cum_w, s_vals)
                        else:
                            q16, q50, q84 = np.percentile(vals, [16, 50, 84])
                    else:
                        q16, q50, q84 = np.percentile(vals, [16, 50, 84])
                    n_bins = _get_hist_bins(len(vals), user_bins=cfg_overview_bins, is_overview=True)
                    if w_arr is not None and len(w_arr) == len(vals):
                        ax.hist(vals, weights=w_arr, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                    else:
                        ax.hist(vals, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                    ax.axvline(q50, color="#1A2E40", linestyle="--", linewidth=1.5)
                    ax.axvline(q16, color="#2B4C7E", linestyle=":", linewidth=1.2)
                    ax.axvline(q84, color="#2B4C7E", linestyle=":", linewidth=1.2)
                    ax.axvspan(q16, q84, color="#2B4C7E", alpha=0.15)

                    short_name = format_mcmc_parameter_label(pname, short=True)

                    ax.set_title(f"{short_name} [{q50:.3f}]", fontsize=10, pad=4)
                    if should_show_grid(config):
                        ax.grid(True, linestyle=":", alpha=0.4)
                    ax.tick_params(axis="both", which="major", labelsize=8)

                for i in range(len(key_params), len(axes_flat)):
                    axes_flat[i].set_visible(False)

                fig.suptitle(f"{config.project.name} - MCMC Resonance Parameters Posterior Distributions", fontsize=12.5, y=0.975)
                plt.tight_layout(rect=[0, 0.02, 1, 0.94])
                for fmt in formats:
                    p_dist_grid = dist_dir / f"distributions_overview.{fmt}"
                    fig.savefig(p_dist_grid, bbox_inches="tight")
                    generated.append(p_dist_grid)
                plt.close(fig)

        # 6b. Amplitude Yield / Integrated Intensity Posterior Distributions
        if mcmc_samples and masses:
            yield_samples: Dict[str, List[float]] = {}
            first_m = round(masses[0], 4)
            for k in sorted(magnitudes):
                m_dict = mcmc_samples.get(k, {})
                n_s = len(m_dict.get(first_m, []))
                if n_s > 0:
                    y_samples = []
                    for s_idx in range(n_s):
                        tot = 0.0
                        for m in masses:
                            m_r = round(m, 4)
                            samps = m_dict.get(m_r, [])
                            if s_idx < len(samps):
                                tot += (samps[s_idx] ** 2) * fit_step
                        y_samples.append(tot)
                    if y_samples and np.std(y_samples) > 1e-12:
                        yield_samples[k] = y_samples

            for k, y_vals in yield_samples.items():
                vals = np.array(y_vals, dtype=float)
                w_samps = mcmc_weights.get(k, {}).get(first_m, None)
                if w_samps is None and "weight" in mcmc_params:
                    w_samps = mcmc_params["weight"]
                w_arr = np.array(w_samps, dtype=float) if (w_samps is not None and len(w_samps) == len(vals)) else None

                fig, ax = plt.subplots(
                    figsize=(5.5, 4.2), dpi=300,
                    gridspec_kw={"top": 0.91, "bottom": 0.13, "left": 0.13, "right": 0.95}
                )
                if w_arr is not None:
                    q16 = compute_weighted_percentile(vals, 16.0, w_arr)
                    q50 = compute_weighted_percentile(vals, 50.0, w_arr)
                    q84 = compute_weighted_percentile(vals, 84.0, w_arr)
                else:
                    q16, q50, q84 = np.percentile(vals, [16, 50, 84])
                nom_yield = sum((v**2) * fit_step for v in nominal_amps.get(k, []))

                n_bins = _get_hist_bins(len(vals), user_bins=cfg_dist_bins, is_overview=False)
                if w_arr is not None:
                    ax.hist(vals, weights=w_arr, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                else:
                    ax.hist(vals, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                ax.axvline(q50, color="#1A2E40", linestyle="--", linewidth=1.5, label=f"Median: {q50:.2f}")
                ax.axvline(q16, color="#2B4C7E", linestyle=":", linewidth=1.2, label=f"68% CI: [{q16:.2f}, {q84:.2f}]")
                ax.axvline(q84, color="#2B4C7E", linestyle=":", linewidth=1.2)
                ax.axvspan(q16, q84, color="#2B4C7E", alpha=0.15)
                if nom_yield > 0:
                    ax.axvline(nom_yield, color="#1A2E40", linestyle="-", linewidth=1.5, label=f"Nominal: {nom_yield:.2f}")

                sqr_label = format_amplitude_squared_label(k)
                sqr_sym = format_amplitude_squared_symbol(k)
                ax.set_title(f"Propagated Amplitude Posterior: $\\int {sqr_sym}\\ dM$", fontsize=11, pad=6)
                ax.set_xlabel(f"Integrated Intensity $\\int {sqr_sym}\\ dM$", fontsize=10)
                ax.set_ylabel("Probability Density", fontsize=10)
                ax.legend(frameon=True, fontsize=8, loc="upper right")
                if should_show_grid(config):
                    ax.grid(True, linestyle=":", alpha=0.4)

                for fmt in formats:
                    p_ydist = dist_dir / f"{k}_yield_dist.{fmt}"
                    fig.savefig(p_ydist, bbox_inches="tight")
                    generated.append(p_ydist)
                plt.close(fig)

            # Overview grid of amplitude yields
            if yield_samples:
                y_keys = list(yield_samples.keys())
                ncols = min(4, len(y_keys))
                nrows = math.ceil(len(y_keys) / ncols)
                fig, axes = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 3.0 * nrows), dpi=300)
                axes_flat = [axes] if len(y_keys) == 1 else (axes.flatten() if hasattr(axes, "flatten") else [axes])

                for i, k in enumerate(y_keys):
                    ax = axes_flat[i]
                    vals = np.array(yield_samples[k], dtype=float)
                    w_samps = mcmc_weights.get(k, {}).get(first_m, None)
                    if w_samps is None and "weight" in mcmc_params:
                        w_samps = mcmc_params["weight"]
                    w_arr = np.array(w_samps, dtype=float) if (w_samps is not None and len(w_samps) == len(vals)) else None

                    if w_arr is not None:
                        q16 = compute_weighted_percentile(vals, 16.0, w_arr)
                        q50 = compute_weighted_percentile(vals, 50.0, w_arr)
                        q84 = compute_weighted_percentile(vals, 84.0, w_arr)
                    else:
                        q16, q50, q84 = np.percentile(vals, [16, 50, 84])
                    n_bins = _get_hist_bins(len(vals), user_bins=cfg_overview_bins, is_overview=True)
                    if w_arr is not None:
                        ax.hist(vals, weights=w_arr, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                    else:
                        ax.hist(vals, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                    ax.axvline(q50, color="#1A2E40", linestyle="--", linewidth=1.5)
                    ax.axvline(q16, color="#2B4C7E", linestyle=":", linewidth=1.2)
                    ax.axvline(q84, color="#2B4C7E", linestyle=":", linewidth=1.2)
                    ax.axvspan(q16, q84, color="#2B4C7E", alpha=0.15)

                    sqr_label = format_amplitude_squared_label(k)
                    ax.set_title(f"{sqr_label} [{q50:.1f}]", fontsize=10, pad=4)
                    if should_show_grid(config):
                        ax.grid(True, linestyle=":", alpha=0.4)
                    ax.tick_params(axis="both", which="major", labelsize=8)

                for i in range(len(y_keys), len(axes_flat)):
                    axes_flat[i].set_visible(False)

                fig.suptitle(f"{config.project.name} - Propagated Amplitude Yields Posterior Distributions", fontsize=12.5, y=0.975)
                plt.tight_layout(rect=[0, 0.02, 1, 0.94])
                for fmt in formats:
                    p_ydist_grid = dist_dir / f"amplitude_yields_overview.{fmt}"
                    fig.savefig(p_ydist_grid, bbox_inches="tight")
                    generated.append(p_ydist_grid)
                plt.close(fig)

    return generated


def load_or_export_moment_comparison_data(
    config: FitConfig,
    rundir: RunDirectory,
    unnormalized_moments_file: Optional[Path | str] = None,
    n_lowest_chi2: Optional[int] = None,
) -> Optional[Dict[str, Any]]:
    """Load cached moment comparison JSON or export it using ExportMomentComparison.C."""
    cache_dir = rundir.root / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    json_path = cache_dir / "moment_comparison_data.json"

    unnorm_path = resolve_unnormalized_moments_file(config, rundir, unnormalized_moments_file)

    # Check if cache exists and is fresh
    if json_path.is_file() and json_path.stat().st_size > 10:
        is_stale = False
        json_mtime = json_path.stat().st_mtime
        if unnorm_path and unnorm_path.is_file():
            if json_mtime < unnorm_path.stat().st_mtime:
                is_stale = True
        if rundir.nominal_best_root.is_file() and json_mtime < rundir.nominal_best_root.stat().st_mtime:
            is_stale = True
        if not is_stale:
            try:
                with open(json_path, "r") as f:
                    data = json.load(f)
                if "data" in data and data["data"]:
                    return data
            except Exception:
                pass

    if not unnorm_path or not unnorm_path.is_file():
        resolved_dir_str = config.dataset.get_moments_dir() if hasattr(config.dataset, "get_moments_dir") else config.dataset.moments_dir
        print(f"[WARNING] Unnormalized moments file not found in {resolved_dir_str}.")
        return None

    # Resolve lowest chi2 fit files
    n_lowest = n_lowest_chi2 or 10
    best_fits = discover_lowest_chi2_fits(rundir, n_lowest=n_lowest)
    if not best_fits and rundir.nominal_best_root.is_file():
        best_fits = [(0.0, rundir.nominal_best_root)]

    if not best_fits:
        print(f"[WARNING] No fit result ROOT files found to overlay in {rundir.fit_dir}.")
        return None

    fit_files = [str(f[1].resolve()) for f in best_fits]
    files_str = ",".join(fit_files)

    repo_root = Path(os.environ.get("MOMS2AMPS", str(Path(__file__).resolve().parent.parent.parent)))
    export_macro = repo_root / "macros" / "ExportMomentComparison.C"
    if not export_macro.is_file():
        print(f"[WARNING] Export macro not found at {export_macro}.")
        return None

    root_bin = get_root_binary(config)
    cmd = [
        root_bin,
        "-l",
        "-b",
        "-q",
        f'{export_macro}("{unnorm_path.resolve()}", "{files_str}", "{json_path.resolve()}")'
    ]
    env = get_root_subprocess_env(config)
    run_res = subprocess.run(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if run_res.returncode != 0 or not json_path.is_file():
        print(f"[WARNING] ExportMomentComparison failed:\n{run_res.stderr}\n{run_res.stdout}")
        return None

    try:
        with open(json_path, "r") as f:
            return json.load(f)
    except Exception as e:
        print(f"[WARNING] Could not parse exported moment comparison JSON: {e}")
        return None


def plot_moments_comparison(
    config: FitConfig,
    rundir: RunDirectory,
    n_lowest_chi2: Optional[int] = None,
    moments: Optional[List[str]] = None,
    unnormalized_moments_file: Optional[Path | str] = None,
    formats: Tuple[str, ...] = ("png", "pdf"),
    style: Optional[str] = None,
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

    effective_style = style or (getattr(config.plotting, "style", "default") if getattr(config, "plotting", None) else "default")
    apply_plot_style(effective_style)

    # Resolve number of lowest chi2 results to draw
    n_lowest = n_lowest_chi2
    if n_lowest is None and getattr(config, "plotting", None):
        n_lowest = getattr(config.plotting, "n_lowest_chi2", None)
    if n_lowest is None and getattr(config, "multi_start", None):
        n_lowest = getattr(config.multi_start, "n_lowest_chi2", None)
    if n_lowest is None or n_lowest <= 0:
        n_lowest = 10

    comp_data = load_or_export_moment_comparison_data(
        config=config,
        rundir=rundir,
        unnormalized_moments_file=unnormalized_moments_file,
        n_lowest_chi2=n_lowest,
    )
    if not comp_data:
        return generated

    data_moments: Dict[str, Any] = comp_data.get("data", {})
    raw_fits: List[Dict[str, Any]] = comp_data.get("fits", [])
    fits_list: List[Dict[str, Any]] = [
        f for f in raw_fits
        if f.get("status", 0) == 0 and f.get("is_valid", True)
    ]

    if not data_moments:
        print(f"[WARNING] No moment data found for {rundir.root}.")
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

        # Restrict data points strictly to the fitted mass range
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

        fig, ax = plt.subplots(
            figsize=(7, 4.8), dpi=300,
            gridspec_kw={"top": 0.92, "bottom": 0.12, "left": 0.12, "right": 0.96}
        )
        ax.errorbar(
            m_data, y_data, yerr=e_data,
            fmt="o", color="black", markersize=3.2,
            capsize=1.2, elinewidth=0.9,
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
                            m_smooth, y_smooth, "-", color="#D55E00", linewidth=1.8, alpha=0.95,
                            label=fr"Rank 1 ($\chi^2 = {fit.get('chi2', 0):.1f}$)", zorder=5
                        )
                    else:
                        lbl = f"Top {len(fits_list)} Converged Fits" if idx == 1 else None
                        ax.plot(
                            m_smooth, y_smooth, "-", color="#D55E00", linewidth=0.9, alpha=0.30,
                            label=lbl, zorder=4
                        )

        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
        if fit_mass_min is not None and fit_mass_max is not None:
            pad = fit_mass_step * 0.6
            ax.set_xlim(fit_mass_min - pad, fit_mass_max + pad)

        ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11.5)
        ax.set_ylabel("Unnormalized Moment", fontsize=11.5)
        ax.set_title(f"{config.project.name} - ${latex_title}$", fontsize=12.5, pad=8)
        ax.legend(frameon=True, fontsize=9.5, loc="best")
        if should_show_grid(config):
            ax.grid(True, linestyle=":", alpha=0.45)

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

        fig, axes = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 2.8 * nrows), dpi=200)
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
                fmt="o", color="black", markersize=2.0,
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
                            ax.plot(m_smooth, y_smooth, "-", color="#D55E00", linewidth=1.4, alpha=0.95, zorder=5)
                        else:
                            ax.plot(m_smooth, y_smooth, "-", color="#D55E00", linewidth=0.8, alpha=0.30, zorder=4)

            # Restrict x-axis limits strictly to the fitted mass range
            if fit_mass_min is not None and fit_mass_max is not None:
                pad = fit_mass_step * 0.6
                ax.set_xlim(fit_mass_min - pad, fit_mass_max + pad)

            ax.set_title(f"${latex_title}$", fontsize=10.5, pad=4)
            if should_show_grid(config):
                ax.grid(True, linestyle=":", alpha=0.4)
            ax.tick_params(axis="both", which="major", labelsize=8)
            if (idx + ncols >= K) or (idx >= (nrows - 1) * ncols):
                ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=9)

        for idx in range(K, len(axes_flat)):
            axes_flat[idx].set_visible(False)

        fig.suptitle(
            f"{config.project.name} - Unnormalized Moments vs Mass (Top {len(fits_list)} Converged Fits)",
            fontsize=13,
            y=0.975,
        )
        plt.tight_layout(rect=[0, 0.02, 1, 0.94])

        for fmt in formats:
            p_grid = rundir.plots_dir / f"moments_overview.{fmt}"
            fig.savefig(p_grid, bbox_inches="tight")
            generated.append(p_grid)
        plt.close(fig)

    print(f"[OK] Generated {len(generated)} moment comparison plot files in {rundir.plots_dir}")
    return generated


def compute_moment_chi2(
    data_moment: Dict[str, Any],
    fit_moment: Optional[Dict[str, Any]],
    fit_mass_min: Optional[float] = None,
    fit_mass_max: Optional[float] = None,
    fit_mass_step: float = 0.04,
) -> Optional[Tuple[float, int]]:
    """
    Compute chi-square between experimental data points and model curve for an individual moment.
    Restricted to the fitted mass range [fit_mass_min, fit_mass_max].
    Returns (chi2, num_points) or None if data or fit is unavailable.
    """
    if not data_moment or not fit_moment:
        return None

    import numpy as np

    m_data = np.asarray(data_moment.get("mass", []))
    y_data = np.asarray(data_moment.get("val", []))
    e_data = np.asarray(data_moment.get("err", []))

    if len(m_data) == 0 or len(y_data) == 0 or len(e_data) == 0:
        return None

    # Restrict to fitted mass range
    if fit_mass_min is not None and fit_mass_max is not None:
        tol = fit_mass_step * 0.5
        mask = (m_data >= (fit_mass_min - tol)) & (m_data <= (fit_mass_max + tol))
        m_data = m_data[mask]
        y_data = y_data[mask]
        e_data = e_data[mask]

    # Filter valid points with positive error and finite values
    valid = (e_data > 0) & np.isfinite(y_data) & np.isfinite(e_data)
    m_data = m_data[valid]
    y_data = y_data[valid]
    e_data = e_data[valid]

    if len(m_data) == 0:
        return None

    mf = np.asarray(fit_moment.get("mass", []))
    yf = np.asarray(fit_moment.get("val", []))
    if len(mf) == 0 or len(yf) == 0:
        return None

    # Sort fit arrays if needed
    sort_idx = np.argsort(mf)
    mf = mf[sort_idx]
    yf = yf[sort_idx]

    # Interpolate model to data mass points
    y_model = np.interp(m_data, mf, yf)
    residuals = (y_data - y_model) / e_data
    chi2 = float(np.sum(residuals ** 2))
    return chi2, len(m_data)


def _draw_overview_grid_cells(
    axes: Any,
    moments_subset: List[str],
    sorted_tbins: List[str],
    t_data: Dict[str, Any],
    fit_mass_min: Optional[float],
    fit_mass_max: Optional[float],
    fit_mass_step: float,
    draw_fits: bool,
    show_grid: bool = False,
) -> None:
    """Populate grid cells for multi-t-bin moments overview matrix."""
    import numpy as np
    n_r = len(moments_subset)
    n_c = len(sorted_tbins)
    pad = fit_mass_step * 0.6 if (fit_mass_min is not None and fit_mass_max is not None) else 0.024

    for r_idx, mom_name in enumerate(moments_subset):
        parts = mom_name.split("_")
        if len(parts) >= 4:
            latex_mom = f"H_{{{parts[1]}}}({parts[2]},{parts[3]})"
        else:
            latex_mom = mom_name

        for c_idx, t_name in enumerate(sorted_tbins):
            ax = axes[r_idx, c_idx]
            t_info = t_data[t_name]
            d_val = t_info["data"].get(mom_name)

            f_mom_best = None
            if draw_fits and t_info.get("best_fit"):
                f_mom_best = t_info["best_fit"].get("moments", {}).get(mom_name)

            if d_val:
                m_pts = d_val.get("mass", [])
                y_pts = d_val.get("val", [])
                e_pts = d_val.get("err", [])

                if fit_mass_min is not None and fit_mass_max is not None:
                    tol = fit_mass_step * 0.5
                    filtered = [
                        (m, y, e) for m, y, e in zip(m_pts, y_pts, e_pts)
                        if (fit_mass_min - tol) <= m <= (fit_mass_max + tol)
                    ]
                    if filtered:
                        m_pts = [p[0] for p in filtered]
                        y_pts = [p[1] for p in filtered]
                        e_pts = [p[2] for p in filtered]

                ax.errorbar(
                    m_pts, y_pts, yerr=e_pts,
                    fmt="o", color="black", markersize=2.0,
                    capsize=0.8, elinewidth=0.6, zorder=3,
                )

            if draw_fits and f_mom_best:
                mf = f_mom_best.get("mass", [])
                yf = f_mom_best.get("val", [])
                if len(mf) > 3:
                    ms, ys = smooth_curve(mf, yf, num_points=200)
                    ax.plot(ms, ys, "-", color="#D55E00", linewidth=1.3, zorder=4)

            # Compute moment chi2 (and show global chi2 explicitly for H_0_0_0)
            if d_val and f_mom_best:
                chi2_info = compute_moment_chi2(
                    d_val,
                    f_mom_best,
                    fit_mass_min=fit_mass_min,
                    fit_mass_max=fit_mass_max,
                    fit_mass_step=fit_mass_step,
                )
                c2_glob = t_info["best_fit"].get("chi2", None) if t_info.get("best_fit") else None
                if mom_name == "H_0_0_0" and c2_glob is not None:
                    txt = f"Global $\\chi^2={c2_glob:.1f}$"
                    if chi2_info is not None:
                        c2_mom, _ = chi2_info
                        txt += f"\n$\\chi^2_{{\\mathrm{{mom}}}}={c2_mom:.1f}$"
                    ax.text(
                        0.96, 0.90, txt,
                        transform=ax.transAxes, ha="right", va="top",
                        fontsize=7.0,
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85, edgecolor="none"),
                        zorder=5,
                    )
                elif chi2_info is not None:
                    c2_mom, _ = chi2_info
                    ax.text(
                        0.96, 0.90, f"$\\chi^2={c2_mom:.1f}$",
                        transform=ax.transAxes, ha="right", va="top",
                        fontsize=7.5,
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85, edgecolor="none"),
                        zorder=5,
                    )

            ax.axhline(0, color="gray", linestyle="--", linewidth=0.6, alpha=0.5, zorder=1)
            if fit_mass_min is not None and fit_mass_max is not None:
                ax.set_xlim(fit_mass_min - pad, fit_mass_max + pad)

            if show_grid:
                ax.grid(True, linestyle=":", alpha=0.35)
            ax.tick_params(axis="both", which="major", labelsize=7.5)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)

            if r_idx == 0:
                ax.set_title(t_info["label"], fontsize=9.5, pad=4)
            if c_idx == 0:
                ax.set_ylabel(f"${latex_mom}$", fontsize=10.5)
            if r_idx == n_r - 1:
                ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=8.5)


def plot_moments_across_tbins(
    config: FitConfig,
    tbins: Optional[List[str]] = None,
    run_dir: Optional[str] = None,
    output_dir: Optional[Path | str] = None,
    moments: Optional[List[str]] = None,
    formats: Tuple[str, ...] = ("png", "pdf"),
    style: Optional[str] = None,
    draw_best_fits: Optional[bool] = None,
    sharey: bool = False,
) -> List[Path]:
    """
    Plot each moment as a function of mass across all t-bins side-by-side (N_tbins x 1 subplots).
    Enables inspecting how angular moments evolve across Mandelstam -t transfers.
    Draws subplots from lowest t-bin (left) to highest t-bin (right) with data points
    and nominal model fit curves for each t-bin.

    Generates:
      - Individual moment N_tbins x 1 plots: <out_dir>/moments_vs_t/<moment_name>.<fmt>
      - Multi-panel overview grid of key moments: <out_dir>/moments_vs_t/moments_overview_all_tbins.<fmt>
    """
    generated: List[Path] = []

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("[WARNING] matplotlib not found. Skipping multi-t-bin moment comparison plots.")
        return generated

    tbins_list = tbins if tbins is not None else config.dataset.get_t_bins()
    if not tbins_list or len(tbins_list) < 2:
        return generated

    effective_style = style or (getattr(config.plotting, "style", "default") if getattr(config, "plotting", None) else "default")
    apply_plot_style(effective_style)

    from .cross_section import parse_t_bin_from_string
    from .rundir import get_run_directory

    # Resolve output directory
    if output_dir:
        out_dir = Path(output_dir)
    else:
        base_dir = None
        if getattr(config, "multi_start", None) and getattr(config.multi_start, "output_dir", None):
            cand = Path(config.multi_start.output_dir)
            if cand.is_dir():
                base_dir = cand
        if base_dir is None and run_dir:
            cand = Path(run_dir)
            if cand.is_dir():
                base_dir = cand
        if base_dir is None:
            cfg_first = copy.deepcopy(config)
            cfg_first.dataset.t_bin = tbins_list[0]
            first_rdir = get_run_directory(cfg_first, run_dir)
            base_dir = first_rdir.root.parent

        out_dir = base_dir / "plots" / "moments_vs_t"

    indiv_dir = out_dir / "individual"
    out_dir.mkdir(parents=True, exist_ok=True)
    indiv_dir.mkdir(parents=True, exist_ok=True)

    # Load moment data for each t-bin
    t_data: Dict[str, Dict[str, Any]] = {}
    for t_name in tbins_list:
        cfg_t = copy.deepcopy(config)
        cfg_t.dataset.t_bin = t_name
        t_rundir = get_run_directory(cfg_t, run_dir)
        n_chi2 = getattr(config.plotting, "n_lowest_chi2", 1) if getattr(config, "plotting", None) else 1
        comp = load_or_export_moment_comparison_data(cfg_t, t_rundir, n_lowest_chi2=n_chi2)
        if not comp:
            continue

        parsed_t = parse_t_bin_from_string(t_name)
        if parsed_t:
            t_min, t_max = parsed_t
            t_lbl = f"${_format_t_val(t_min)} < -t < {_format_t_val(t_max)}\\ \\mathrm{{GeV}}^2$"
        else:
            t_min, t_max = 0.0, 0.0
            t_lbl = t_name

        raw_fits = comp.get("fits", [])
        valid_fits = [
            f for f in raw_fits
            if f.get("status", 0) == 0 and f.get("is_valid", True)
        ]
        best_fit = valid_fits[0] if valid_fits else None

        t_data[t_name] = {
            "name": t_name,
            "t_min": t_min,
            "t_max": t_max,
            "label": t_lbl,
            "data": comp.get("data", {}),
            "best_fit": best_fit,
            "fits": valid_fits,
            "rundir": t_rundir,
        }

    sorted_tbins = [t for t in tbins_list if t in t_data]
    sorted_tbins.sort(key=lambda t: (t_data[t]["t_min"] if t_data[t]["t_min"] is not None else 0.0))

    if not sorted_tbins:
        print("[WARNING] No moment comparison data available across t-bins.")
        return generated

    # Filter requested moments
    req_moments = moments if moments is not None else getattr(config.plotting, "moments", ["all"])
    if not req_moments:
        req_moments = ["all"]

    all_moment_keys = set()
    for t in sorted_tbins:
        all_moment_keys.update(t_data[t]["data"].keys())

    selected_moments: List[str] = []
    if any(m.lower().strip() == "all" for m in req_moments):
        selected_moments = list(all_moment_keys)
    elif any(m.lower().strip() == "fitted" for m in req_moments):
        selected_moments = [k for k in all_moment_keys if is_moment_fitted(k, config)]
    else:
        for m_req in req_moments:
            m_str = m_req.strip()
            m_lower = m_str.lower()
            if m_lower.startswith("l="):
                try:
                    target_l = int(m_lower.replace("l=", ""))
                    for k in all_moment_keys:
                        parts = k.split("_")
                        if len(parts) >= 4 and int(parts[2]) == target_l:
                            if k not in selected_moments:
                                selected_moments.append(k)
                except ValueError:
                    pass
            elif m_str in all_moment_keys:
                if m_str not in selected_moments:
                    selected_moments.append(m_str)
            else:
                if not m_str.startswith("H_"):
                    alt = f"H_{m_str}"
                    if alt in all_moment_keys and alt not in selected_moments:
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
        return generated

    draw_fits = draw_best_fits if draw_best_fits is not None else getattr(config.plotting, "draw_best_fits", True)

    # Fitted mass range from config
    fit_mass_min = None
    fit_mass_max = None
    fit_mass_step = 0.04
    if hasattr(config, "dataset") and getattr(config.dataset, "mass_bins", None):
        mb = config.dataset.mass_bins
        fit_mass_min = getattr(mb, "min", None)
        fit_mass_max = getattr(mb, "max", None)
        fit_mass_step = getattr(mb, "step", 0.04) or 0.04

    if fit_mass_min is None or fit_mass_max is None:
        all_fit_masses = []
        for t_info in t_data.values():
            if t_info.get("best_fit"):
                for m_obj in t_info["best_fit"].get("moments", {}).values():
                    all_fit_masses.extend(m_obj.get("mass", []))
        if all_fit_masses:
            unique_m = sorted(list(set(all_fit_masses)))
            if fit_mass_min is None:
                fit_mass_min = unique_m[0]
            if fit_mass_max is None:
                fit_mass_max = unique_m[-1]
            if len(unique_m) > 1:
                fit_mass_step = float(np.median(np.diff(unique_m)))

    palette = PHYSICS_PALETTE
    markers = ["o", "s", "^", "D", "v", "<", ">", "p", "*", "X"]

    # 1. Individual Moment Plots across all t-bins (side-by-side from lowest to highest -t)
    n_cols = len(sorted_tbins)
    fig_width = max(3.8 * n_cols, 8.0)
    fig_height = 4.2

    for mom_name in selected_moments:
        parts = mom_name.split("_")
        if len(parts) >= 4:
            alpha, l_val, m_val = parts[1], parts[2], parts[3]
            latex_title = f"H_{{{alpha}}}({l_val},{m_val})"
        else:
            latex_title = mom_name

        fig, axes = plt.subplots(
            1, n_cols,
            figsize=(fig_width, fig_height),
            sharey=sharey,
            dpi=300,
        )
        axes_list = np.atleast_1d(axes)

        for idx, t_name in enumerate(sorted_tbins):
            ax = axes_list[idx]
            t_info = t_data[t_name]
            t_lbl = t_info.get("label", t_info.get("name", t_name))

            d_val = t_info["data"].get(mom_name)
            if d_val:
                m_pts = d_val.get("mass", [])
                y_pts = d_val.get("val", [])
                e_pts = d_val.get("err", [])

                if fit_mass_min is not None and fit_mass_max is not None:
                    tol = fit_mass_step * 0.5
                    filtered = [
                        (m, y, e) for m, y, e in zip(m_pts, y_pts, e_pts)
                        if (fit_mass_min - tol) <= m <= (fit_mass_max + tol)
                    ]
                    if filtered:
                        m_pts = [p[0] for p in filtered]
                        y_pts = [p[1] for p in filtered]
                        e_pts = [p[2] for p in filtered]

                ax.errorbar(
                    m_pts, y_pts, yerr=e_pts,
                    fmt="o", color="black", markersize=3.2,
                    capsize=1.2, elinewidth=0.85,
                    label="Data" if idx == 0 else None,
                    zorder=3,
                )

            f_mom_best = None
            if draw_fits:
                fits_list = t_info.get("fits", [])
                if not fits_list and t_info.get("best_fit"):
                    fits_list = [t_info["best_fit"]]
                for f_idx, fit_obj in enumerate(fits_list):
                    f_mom = fit_obj.get("moments", {}).get(mom_name)
                    if f_mom:
                        if f_idx == 0:
                            f_mom_best = f_mom
                        mf = f_mom.get("mass", [])
                        yf = f_mom.get("val", [])
                        if len(mf) > 3:
                            ms, ys = smooth_curve(mf, yf, num_points=200)
                            if f_idx == 0:
                                ax.plot(
                                    ms, ys, "-", color="#D55E00", linewidth=1.8, alpha=0.95,
                                    label="Rank 1 Fit" if idx == 0 else None, zorder=5
                                )
                            else:
                                lbl = f"Top {len(fits_list)} Fits" if (idx == 0 and f_idx == 1) else None
                                ax.plot(
                                    ms, ys, "-", color="#D55E00", linewidth=0.9, alpha=0.25,
                                    label=lbl, zorder=4
                                )

            ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6, zorder=1)
            if fit_mass_min is not None and fit_mass_max is not None:
                pad = fit_mass_step * 0.6
                ax.set_xlim(fit_mass_min - pad, fit_mass_max + pad)

            # Compute chi2 for this individual moment (and show global chi2 explicitly for H_0_0_0)
            chi2_str = ""
            c2_glob = t_info["best_fit"].get("chi2", None) if t_info.get("best_fit") else None
            chi2_info = None
            if d_val and f_mom_best:
                chi2_info = compute_moment_chi2(
                    d_val,
                    f_mom_best,
                    fit_mass_min=fit_mass_min,
                    fit_mass_max=fit_mass_max,
                    fit_mass_step=fit_mass_step,
                )

            if mom_name == "H_0_0_0" and c2_glob is not None:
                if chi2_info is not None:
                    c2_mom, _ = chi2_info
                    chi2_str = f"\n(Global $\\chi^2 = {c2_glob:.1f}$, $\\chi^2_{{\\mathrm{{mom}}}} = {c2_mom:.1f}$)"
                else:
                    chi2_str = f"\n(Global $\\chi^2 = {c2_glob:.1f}$)"
            elif chi2_info is not None:
                c2_mom, _ = chi2_info
                chi2_str = f"\n($\\chi^2 = {c2_mom:.1f}$)"

            ax.set_title(f"{t_lbl}{chi2_str}", fontsize=9.8, pad=5)
            if should_show_grid(config):
                ax.grid(True, linestyle=":", alpha=0.45)
            ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=10)
            ax.tick_params(axis="both", which="major", labelsize=9)

            if idx == 0:
                ax.set_ylabel(f"${latex_title}$", fontsize=11.5)
                ax.legend(frameon=True, fontsize=8.5, loc="best")
            elif not sharey:
                ax.set_ylabel(f"${latex_title}$", fontsize=10, alpha=0.85)

            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)

        fig.suptitle(f"{config.project.name} - ${latex_title}$ Across $-t$ Bins", fontsize=13, y=0.98)
        fig.tight_layout(rect=[0, 0.02, 1, 0.93])

        for fmt in formats:
            p_indiv = indiv_dir / f"{mom_name}_vs_t.{fmt}"
            p_top = out_dir / f"{mom_name}.{fmt}"
            fig.savefig(p_indiv, bbox_inches="tight")
            fig.savefig(p_top, bbox_inches="tight")
            generated.append(p_indiv)
            generated.append(p_top)
        plt.close(fig)

    # 2. Multi-panel Overview Grid (Matrix: Rows = Moments, Columns = -t Bins)
    overview_moments = [m for m in selected_moments if is_moment_fitted(m, config)]
    if not overview_moments:
        overview_moments = selected_moments

    K = len(overview_moments)
    n_c = len(sorted_tbins)

    if K >= 1 and n_c >= 1:
        # Multi-page PDF for clear, readable letter/slide viewing
        if "pdf" in formats:
            from matplotlib.backends.backend_pdf import PdfPages
            p_pdf = out_dir / "moments_overview_all_tbins.pdf"
            moms_per_page = 5
            total_pages = math.ceil(K / moms_per_page)
            with PdfPages(p_pdf) as pdf:
                for page_idx in range(0, K, moms_per_page):
                    page_moms = overview_moments[page_idx : page_idx + moms_per_page]
                    curr_nr = len(page_moms)
                    fig_p, axes_p = plt.subplots(
                        curr_nr, n_c,
                        figsize=(3.2 * n_c, 2.2 * curr_nr + 0.6),
                        sharex=True, dpi=150,
                    )
                    axes_p = np.atleast_2d(axes_p)
                    fig_p.subplots_adjust(hspace=0.28, wspace=0.32, left=0.08, right=0.98, top=0.92, bottom=0.08)

                    _draw_overview_grid_cells(
                        axes_p, page_moms, sorted_tbins, t_data,
                        fit_mass_min, fit_mass_max, fit_mass_step, draw_fits,
                        show_grid=should_show_grid(config)
                    )
                    page_num = (page_idx // moms_per_page) + 1
                    fig_p.suptitle(
                        f"{config.project.name} - Moments vs Mass Across $-t$ Bins (Page {page_num}/{total_pages})",
                        fontsize=12, y=0.98
                    )
                    pdf.savefig(fig_p)
                    plt.close(fig_p)
            generated.append(p_pdf)

        # Full PNG matrix grid
        if "png" in formats:
            p_png = out_dir / "moments_overview_all_tbins.png"
            fig_full, axes_full = plt.subplots(
                K, n_c,
                figsize=(3.2 * n_c, max(2.1 * K + 0.8, 6.0)),
                sharex=True, dpi=130,
            )
            axes_full = np.atleast_2d(axes_full)
            fig_full.subplots_adjust(hspace=0.28, wspace=0.32, left=0.08, right=0.98, top=0.985, bottom=0.02)

            _draw_overview_grid_cells(
                axes_full, overview_moments, sorted_tbins, t_data,
                fit_mass_min, fit_mass_max, fit_mass_step, draw_fits,
                show_grid=should_show_grid(config)
            )
            fig_full.suptitle(
                f"{config.project.name} - Moments vs Mass Across All $-t$ Bins",
                fontsize=13, y=0.998
            )
            fig_full.savefig(p_png, bbox_inches="tight")
            plt.close(fig_full)
            generated.append(p_png)

        # Copy overview files to each individual t-bin's plots directory
        for fmt in formats:
            src_f = out_dir / f"moments_overview_all_tbins.{fmt}"
            if src_f.is_file():
                for t_name in sorted_tbins:
                    t_pdir = t_data[t_name]["rundir"].plots_dir
                    if t_pdir.is_dir():
                        dst_p = t_pdir / f"moments_overview_all_tbins.{fmt}"
                        shutil.copyfile(src_f, dst_p)
                        generated.append(dst_p)

    print(f"[OK] Generated {len(generated)} multi-t-bin moment comparison plot files in {out_dir}")
    return generated


def _format_resonance_display_name(res_name: str) -> str:
    """Format resonance name like 'a2_1320' or 'pi1_1400' into LaTeX math (e.g. $a_{2}(1320)$)."""
    s = str(res_name).strip()
    m = re.match(r"^([a-zA-Z]+?)([0-9]*)_([0-9]+)$", s)
    if m:
        base, j_sub, mass_tag = m.group(1), m.group(2), m.group(3)
        greek_map = {
            "pi": r"\pi",
            "rho": r"\rho",
            "eta": r"\eta",
            "omega": r"\omega",
            "phi": r"\phi",
        }
        base_tex = greek_map.get(base.lower(), base)
        if j_sub:
            lbl = f"${base_tex}_{{{j_sub}}}({mass_tag})$"
        else:
            lbl = f"${base_tex}({mass_tag})$"
        return _safe_mathtext_label(lbl, s)
    clean = s.replace("_", r"\_")
    return _safe_mathtext_label(f"${clean}$", s)


def _derive_two_bw_wave_name(w: Any) -> str:
    """Mirror driver's TwoBreitWigner canonical wave name without circular imports."""
    w_name = getattr(w, "name", "") or ""
    if w_name and "_" in w_name and not w_name.lower().endswith("wave"):
        return w_name
    resonances = getattr(w, "resonances", []) or []
    if len(resonances) >= 2:
        r1, r2 = resonances[0].name, resonances[1].name
        parts1 = r1.split("_", 1)
        parts2 = r2.split("_", 1)
        if len(parts1) == 2 and len(parts2) == 2 and parts1[0] == parts2[0]:
            return f"{r1}_{parts2[1]}"
        return f"{r1}_{r2}"
    return w_name or "two_bw"


def _is_unused_poly_param(pname: str, config: Optional[FitConfig]) -> bool:
    """Return True for polynomial coefficients when all matching waves have poly_order < 0."""
    is_poly = (
        pname.endswith("m_threshold")
        or pname.endswith("m_expansion")
        or "_re_" in pname
        or "_im_" in pname
    )
    if not is_poly:
        return False
    waves = getattr(config, "waves", []) or [] if config else []
    md_waves = [
        w for w in waves
        if str(getattr(getattr(w, "type", ""), "value", getattr(w, "type", ""))).lower() in ("mass_dependent", "model")
    ]
    for w in md_waves:
        if getattr(w, "poly_order", -1) < 0:
            for r in getattr(w, "resonances", []) or []:
                if getattr(r, "name", "") in pname:
                    return True
            if len(md_waves) == 1:
                return True
    return False


def _extract_t_axis_coords(
    t_data: Dict[str, Dict[str, Any]],
    sorted_tbins: List[str],
) -> Tuple[List[float], List[float], bool]:
    """
    Compute (t_centers, t_half_widths, has_numeric_t) from sorted t-bins.
    Falls back to ordinal indices [0..N-1] with tick labels if any bin lacks numeric bounds.
    """
    centers: List[float] = []
    half_widths: List[float] = []
    has_numeric = True
    for t_name in sorted_tbins:
        info = t_data[t_name]
        t_min = float(info.get("t_min", 0.0) or 0.0)
        t_max = float(info.get("t_max", 0.0) or 0.0)
        if t_max > t_min:
            centers.append(0.5 * (t_min + t_max))
            half_widths.append(0.5 * (t_max - t_min))
        else:
            has_numeric = False
            break
    if not has_numeric:
        centers = [float(i) for i in range(len(sorted_tbins))]
        half_widths = [0.0] * len(sorted_tbins)
    return centers, half_widths, has_numeric


def _build_resonance_param_metadata(
    config: Optional[FitConfig],
    pname: str,
) -> Dict[str, Any]:
    """
    Match an internal MD_* parameter name against configured resonances to determine:
      - res_key, res_display
      - role ('mass', 'width', 'coupling_g', 'phase', 'coupling', 'poly', 'other')
      - label (clean LaTeX mathtext), unit, sort_key, init_val, limits
    """
    waves = getattr(config, "waves", []) or [] if config else []

    matched_res = None
    matched_res_idx: Optional[int] = None
    matched_wave_l: Optional[int] = None

    for w in waves:
        w_type = str(getattr(getattr(w, "type", ""), "value", getattr(w, "type", ""))).lower()
        if w_type not in ("mass_dependent", "model"):
            continue
        w_model = str(getattr(getattr(w, "model", ""), "value", getattr(w, "model", ""))).lower()
        resonances = getattr(w, "resonances", []) or []
        if "two_breit_wigner" in w_model and len(resonances) >= 2:
            wave_name = _derive_two_bw_wave_name(w)
            if wave_name in pname or any(r.name in pname for r in resonances[:2]):
                matched_wave_l = getattr(w, "l", None)
                if re.search(r"(?:Mass|Width|k|globalphase|phase)1(?:$|_)", pname):
                    matched_res = resonances[0]
                    matched_res_idx = 1
                elif re.search(r"(?:Mass|Width|k|globalphase|phase)2(?:$|_)", pname) or pname.endswith("_globalphase"):
                    matched_res = resonances[1]
                    matched_res_idx = 2
                else:
                    for idx_r, r in enumerate(resonances[:2], start=1):
                        if r.name in pname:
                            matched_res = r
                            matched_res_idx = idx_r
                            break
                if matched_res is not None:
                    break
        else:
            for r in resonances:
                if r.name and r.name in pname:
                    matched_res = r
                    matched_wave_l = getattr(w, "l", None)
                    break
            if matched_res is not None:
                break

    res_key = matched_res.name if matched_res is not None else "resonance_other"
    res_disp = _format_resonance_display_name(res_key) if matched_res is not None else "Other"

    # Strip wrapping $ from res_disp for compound expressions
    res_tex_inner = res_disp[1:-1] if (res_disp.startswith("$") and res_disp.endswith("$")) else res_disp

    role = "other"
    unit = ""
    init_val: Optional[float] = None
    limits: Optional[List[float]] = None
    label = format_mcmc_parameter_label(pname, short=False)
    role_order = (6, 9, 999, pname)

    m_mass = re.search(r"_Mass([0-9]*)$", pname)
    m_width = re.search(r"_Width([0-9]*)$", pname)
    if m_mass:
        role = "mass"
        unit = r"$\mathrm{GeV}/c^2$"
        idx_s = m_mass.group(1)
        idx_sub = f"_{{{idx_s}}}" if idx_s else ""
        if matched_res is not None:
            init_val = getattr(matched_res, "mass", None)
            limits = getattr(matched_res, "mass_limits", None)
            label = _safe_mathtext_label(f"$M({res_tex_inner})$", pname)
        else:
            label = _safe_mathtext_label(f"$M{idx_sub}$", pname)
        role_order = (0, 0, int(idx_s) if idx_s else 0, pname)
    elif m_width:
        role = "width"
        unit = r"$\mathrm{GeV}/c^2$"
        idx_s = m_width.group(1)
        idx_sub = f"_{{{idx_s}}}" if idx_s else ""
        if matched_res is not None:
            init_val = getattr(matched_res, "width", None)
            limits = getattr(matched_res, "width_limits", None)
            label = _safe_mathtext_label(rf"$\Gamma({res_tex_inner})$", pname)
        else:
            label = _safe_mathtext_label(rf"$\Gamma{idx_sub}$", pname)
        role_order = (1, 0, int(idx_s) if idx_s else 0, pname)
    elif "_g_" in pname or pname.endswith("_g"):
        role = "coupling_g"
        g_match = re.search(r"g_([a-zA-Z0-9]+)$", pname)
        ch = g_match.group(1) if g_match else ""
        ch_low = ch.lower()
        if matched_res is not None:
            if ch_low in ("s", "etapi"):
                init_val = getattr(matched_res, "g_etapi", None) or getattr(matched_res, "g_s", None)
                limits = getattr(matched_res, "g_etapi_limits", None) or getattr(matched_res, "g_s_limits", None)
            elif ch_low == "kk":
                init_val = getattr(matched_res, "g_KK", None)
                limits = getattr(matched_res, "g_KK_limits", None)
            elif ch_low == "etaprimepi":
                init_val = getattr(matched_res, "g_etaprimepi", None)
                limits = getattr(matched_res, "g_etaprimepi_limits", None)
        ch_tex = ch.replace("prime", "'").replace("pi", r"\pi").replace("eta", r"\eta")
        if matched_res is not None:
            label = _safe_mathtext_label(f"$g_{{{ch_tex}}}({res_tex_inner})$", pname)
        else:
            label = _safe_mathtext_label(f"$g_{{{ch_tex}}}$", pname)
        role_order = (2, 0, 0, pname)
    elif "globalphase" in pname or pname.endswith("_phase") or re.search(r"_phase[0-9]+$", pname):
        role = "phase"
        unit = "rad"
        if matched_res is not None:
            limits = getattr(matched_res, "phase_limits", None)
        refl_match = re.search(r"shared([ab])_", pname)
        refl_sym = "^{(+)}" if refl_match and refl_match.group(1) == "a" else ("^{(-)}" if refl_match and refl_match.group(1) == "b" else "")
        refl_rank = 0 if (refl_match and refl_match.group(1) == "a") else (1 if (refl_match and refl_match.group(1) == "b") else 2)
        if matched_res is not None:
            label = _safe_mathtext_label(f"$\\phi{refl_sym}({res_tex_inner})$", pname)
        else:
            label = format_mcmc_parameter_label(pname, short=False)
        role_order = (3, refl_rank, matched_res_idx or 0, pname)
    else:
        m_wave = re.match(r"^MD_([ab])_([0-9]+)_(-?[0-9]+)(?:_(.*))?$", pname)
        if m_wave:
            refl_c, l_s, m_s, rest = m_wave.group(1), m_wave.group(2), m_wave.group(3), (m_wave.group(4) or "")
            amp_sym = format_amplitude_symbol(f"{refl_c}_{l_s}_{m_s}")
            refl_rank = 0 if refl_c == "a" else 1
            if rest == "k" or rest.endswith("_k") or re.search(r"_k[0-9]+$", rest):
                role = "coupling"
                k_m = re.search(r"k([0-9]*)$", rest)
                k_sub = f"_{{{k_m.group(1)}}}" if (k_m and k_m.group(1)) else ""
                if matched_res is not None:
                    label = _safe_mathtext_label(f"$k{k_sub}({amp_sym};\\, {res_tex_inner})$", pname)
                else:
                    label = _safe_mathtext_label(f"$k{k_sub}({amp_sym})$", pname)
                role_order = (4, refl_rank, int(l_s) * 100 + int(m_s), pname)
            elif "re_" in rest or "im_" in rest:
                role = "poly"
                role_order = (5, refl_rank, int(l_s) * 100 + int(m_s), pname)

    return {
        "name": pname,
        "res_key": res_key,
        "res_display": res_disp,
        "wave_l": matched_wave_l,
        "res_idx": matched_res_idx,
        "role": role,
        "unit": unit,
        "label": label,
        "sort_key": role_order,
        "init_val": float(init_val) if init_val is not None else None,
        "limits": limits,
    }


def _plot_single_param_vs_t(
    ax: Any,
    entry: Dict[str, Any],
    show_grid: bool = False,
) -> None:
    """Render one fitted resonance parameter as a function of Mandelstam -t onto axis `ax`."""
    import numpy as np

    t_c = entry["t_centers"]
    t_hw = entry["t_half_widths"]
    y_nom = entry["nom_vals"]
    y_err = entry["nom_errs"]
    has_mcmc = entry["has_mcmc"]
    y_med = entry["mcmc_meds"]
    y_lo = entry["mcmc_lo_errs"]
    y_hi = entry["mcmc_hi_errs"]
    has_numeric_t = entry["has_numeric_t"]
    t_labels = entry["t_labels"]

    # Track plotted y-extent carefully so tight errorbars never get crushed
    all_y_vals = list(y_nom)
    for yv, ye in zip(y_nom, y_err):
        all_y_vals.extend([yv - ye, yv + ye])

    # Offset slightly when both Minuit and MCMC are overlaid cleanly on the same t-bin
    minuit_x = [
        (x - min(0.010, 0.15 * hw if hw > 0 else 0.010)) if (has_mcmc and hm is not None and has_numeric_t) else x
        for x, hw, hm in zip(t_c, t_hw, y_med)
    ]

    ax.plot(t_c, y_nom, ":", color="#666666", linewidth=0.9, alpha=0.55, zorder=2)
    ax.errorbar(
        minuit_x,
        y_nom,
        yerr=y_err if any(abs(e) > 0 for e in y_err) else None,
        xerr=t_hw if has_numeric_t else None,
        fmt="o",
        markersize=5.2,
        color="#1A1A1A",
        ecolor="#1A1A1A",
        elinewidth=1.35,
        capsize=2.8,
        label="Nominal (Minuit)",
        zorder=5,
    )

    if has_mcmc:
        mx, my, my_lo, my_hi, mhw = [], [], [], [], []
        for x, hw, med, lo, hi in zip(t_c, t_hw, y_med, y_lo, y_hi):
            if med is not None:
                dx = min(0.010, 0.15 * hw if hw > 0 else 0.010) if has_numeric_t else 0.12
                mx.append(x + dx)
                my.append(med)
                my_lo.append(lo or 0.0)
                my_hi.append(hi or 0.0)
                mhw.append(hw if has_numeric_t else 0.0)
                all_y_vals.extend([med - (lo or 0.0), med + (hi or 0.0)])
        if mx:
            ax.errorbar(
                mx,
                my,
                yerr=[my_lo, my_hi],
                xerr=mhw if has_numeric_t else None,
                fmt="D",
                markersize=4.6,
                color="#0072B2",
                ecolor="#0072B2",
                elinewidth=1.35,
                capsize=2.8,
                alpha=0.90,
                label="MCMC Median (68% CI)",
                zorder=6,
            )

    # Config / PDG initial value reference line (for shape parameters)
    init_val = entry.get("init_val")
    if init_val is not None and entry["role"] in ("mass", "width", "coupling_g"):
        ax.axhline(
            init_val,
            color="#009E73",
            linestyle="--",
            linewidth=1.15,
            alpha=0.80,
            label=f"Init / PDG ({init_val:.4g})",
            zorder=1,
        )
        all_y_vals.append(init_val)

    # Show configured parameter limit line ONLY if a fitted point lies near/on that limit
    limits = entry.get("limits")
    if limits and len(limits) == 2:
        lim_lo, lim_hi = float(limits[0]), float(limits[1])
        span = max(lim_hi - lim_lo, 1e-6)
        near_lo = any(abs(v - lim_lo) <= 0.08 * span for v in y_nom)
        near_hi = any(abs(lim_hi - v) <= 0.08 * span for v in y_nom)
        if near_lo:
            ax.axhline(lim_lo, color="#D55E00", linestyle=":", linewidth=1.3, alpha=0.85, label=f"Lower Limit ({lim_lo:.3g})", zorder=1)
            all_y_vals.append(lim_lo)
        if near_hi:
            ax.axhline(lim_hi, color="#D55E00", linestyle="-.", linewidth=1.3, alpha=0.85, label=f"Upper Limit ({lim_hi:.3g})", zorder=1)
            all_y_vals.append(lim_hi)

    if entry["role"] == "phase":
        ax.axhline(0.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5, zorder=1)

    if all_y_vals:
        y_min_v = float(np.min(all_y_vals))
        y_max_v = float(np.max(all_y_vals))
        y_span = max(y_max_v - y_min_v, max(abs(y_max_v), 1.0) * 0.05)
        ax.set_ylim(y_min_v - 0.14 * y_span, y_max_v + 0.28 * y_span)

    if has_numeric_t and t_c:
        t_left = min(c - hw for c, hw in zip(t_c, t_hw))
        t_right = max(c + hw for c, hw in zip(t_c, t_hw))
        t_span = max(t_right - t_left, 0.1)
        ax.set_xlim(t_left - 0.04 * t_span, t_right + 0.04 * t_span)
        ax.set_xlabel(r"$-t\ [\mathrm{GeV}^2]$", fontsize=10.5)
    elif t_c:
        ax.set_xticks(t_c)
        ax.set_xticklabels(t_labels, fontsize=8.5)
        ax.set_xlabel(r"$-t$ Bin", fontsize=10.5)

    y_unit = f" [{entry['unit']}]" if entry.get("unit") else ""
    ax.set_ylabel(f"{entry['label']}{y_unit}", fontsize=10.5)
    if show_grid:
        ax.grid(True, linestyle=":", alpha=0.45)


def plot_resonance_parameters_vs_t(
    config: FitConfig,
    t_data: Dict[str, Dict[str, Any]],
    sorted_tbins: List[str],
    output_dir: Path | str,
    formats: Tuple[str, ...] = ("png", "pdf"),
    style: Optional[str] = None,
) -> List[Path]:
    """
    Plot each fitted resonance's parameter values (Mass, Width, Flatté couplings,
    global phases, and wave couplings) as a function of Mandelstam -t.

    Saves:
      - <output_dir>/resonance_parameters_vs_t.<fmt> (overview of shape & phase parameters)
      - <output_dir>/resonance_<res_name>_params_vs_t.<fmt> (per-resonance parameter grid)
      - <output_dir>/parameters_vs_t/<safe_pname>_vs_t.<fmt> (individual parameter plots)
    """
    generated: List[Path] = []
    if len(sorted_tbins) < 2:
        return generated

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return generated

    if style:
        apply_plot_style(style)

    out_dir = Path(output_dir)
    param_dir = out_dir / "parameters_vs_t"
    out_dir.mkdir(parents=True, exist_ok=True)
    param_dir.mkdir(parents=True, exist_ok=True)

    # Collect all MD_* parameter names seen in nominal_params or mcmc_params
    all_md_names = set()
    for t_name in sorted_tbins:
        t_info = t_data[t_name]
        for pk in t_info.get("nominal_params", {}).keys():
            if pk.startswith("MD_") and not pk.endswith("_error"):
                all_md_names.add(pk)
        for pk in t_info.get("mcmc_params", {}).keys():
            if pk.startswith("MD_") and not pk.endswith("_error"):
                all_md_names.add(pk)

    if not all_md_names:
        return generated

    t_centers_all, t_hw_all, has_numeric_t = _extract_t_axis_coords(t_data, sorted_tbins)
    t_labels_all = [t_data[t].get("name", t) for t in sorted_tbins]
    show_grid = should_show_grid(config)

    param_entries: List[Dict[str, Any]] = []
    for pname in sorted(all_md_names):
        if _is_unused_poly_param(pname, config):
            continue

        t_c_list: List[float] = []
        t_hw_list: List[float] = []
        t_lbl_list: List[str] = []
        nom_vals: List[float] = []
        nom_errs: List[float] = []
        mcmc_meds: List[Optional[float]] = []
        mcmc_lo_errs: List[Optional[float]] = []
        mcmc_hi_errs: List[Optional[float]] = []
        has_mcmc_spread = False

        for idx, t_name in enumerate(sorted_tbins):
            t_info = t_data[t_name]
            nom_p = t_info.get("nominal_params", {})
            nom_e = t_info.get("nominal_param_errors", {})
            mc_p = t_info.get("mcmc_params", {})

            val = nom_p.get(pname)
            err = nom_e.get(pname, 0.0)
            mc_arr_raw = mc_p.get(pname)

            if val is None and not mc_arr_raw:
                continue

            # Skip inactive Flatté dummy couplings (fixed at -99.0 with zero error)
            if val is not None and val <= -90.0 and abs(err) < 1e-12:
                continue

            med_v: Optional[float] = None
            lo_e: Optional[float] = None
            hi_e: Optional[float] = None
            if mc_arr_raw and len(mc_arr_raw) >= 2:
                arr = np.asarray(mc_arr_raw, dtype=float)
                if np.std(arr) > 1e-12:
                    has_mcmc_spread = True
                    w_raw = mc_p.get("weight")
                    w_arr = np.asarray(w_raw, dtype=float) if (w_raw is not None and len(w_raw) == len(arr)) else None
                    q16 = compute_weighted_percentile(arr, 16.0, w_arr)
                    q50 = compute_weighted_percentile(arr, 50.0, w_arr)
                    q84 = compute_weighted_percentile(arr, 84.0, w_arr)
                    med_v = q50
                    lo_e = max(0.0, q50 - q16)
                    hi_e = max(0.0, q84 - q50)

            eff_val = float(val) if val is not None else float(med_v if med_v is not None else np.mean(mc_arr_raw))
            t_c_list.append(t_centers_all[idx])
            t_hw_list.append(t_hw_all[idx])
            t_lbl_list.append(t_labels_all[idx])
            nom_vals.append(eff_val)
            nom_errs.append(abs(float(err)))
            mcmc_meds.append(med_v)
            mcmc_lo_errs.append(lo_e)
            mcmc_hi_errs.append(hi_e)

        if not nom_vals:
            continue

        has_minuit_err = any(e > 1e-12 for e in nom_errs)
        has_variation = (max(nom_vals) - min(nom_vals) > 1e-9) if len(nom_vals) > 1 else False
        # Skip strictly fixed / non-fitted constants (e.g. fixed reference phase 0.0 +- 0.0)
        if not (has_minuit_err or has_mcmc_spread or has_variation):
            continue

        meta = _build_resonance_param_metadata(config, pname)
        meta.update({
            "t_centers": t_c_list,
            "t_half_widths": t_hw_list,
            "t_labels": t_lbl_list,
            "has_numeric_t": has_numeric_t,
            "nom_vals": nom_vals,
            "nom_errs": nom_errs,
            "has_mcmc": has_mcmc_spread,
            "mcmc_meds": mcmc_meds,
            "mcmc_lo_errs": mcmc_lo_errs,
            "mcmc_hi_errs": mcmc_hi_errs,
        })
        param_entries.append(meta)

    if not param_entries:
        return generated

    # 1. Individual plot per fitted parameter vs -t
    proj_name = config.project.name if config.project else "PartialWaveFit"
    for entry in param_entries:
        fig, ax = plt.subplots(
            figsize=(6.8, 4.4), dpi=300,
            gridspec_kw={"top": 0.90, "bottom": 0.14, "left": 0.14, "right": 0.96},
        )
        _plot_single_param_vs_t(ax, entry, show_grid=show_grid)
        ax.set_title(f"{proj_name} - {entry['label']} vs $-t$", fontsize=11.5, pad=8)
        ax.legend(frameon=True, fontsize=8.5, loc="best")
        safe_p = re.sub(r"[^a-zA-Z0-9_\-\.]", "_", entry["name"])
        for fmt in formats:
            p_out = param_dir / f"{safe_p}_vs_t.{fmt}"
            fig.savefig(p_out, bbox_inches="tight")
            generated.append(p_out)
        plt.close(fig)

    # Helper to render a multi-panel grid of parameter entries
    def _save_param_grid(entries: List[Dict[str, Any]], title: str, out_base: str) -> None:
        if not entries:
            return
        ordered = sorted(entries, key=lambda e: e["sort_key"])
        n_p = len(ordered)
        ncols = min(3, n_p)
        nrows = math.ceil(n_p / ncols)
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(4.6 * ncols, 3.4 * nrows + 0.5),
            dpi=300,
        )
        axes_flat = [axes] if n_p == 1 else list(np.asarray(axes).ravel())
        for idx, entry in enumerate(ordered):
            ax = axes_flat[idx]
            _plot_single_param_vs_t(ax, entry, show_grid=show_grid)
            ax.set_title(entry["label"], fontsize=11.0, pad=5)
            ax.legend(frameon=True, fontsize=7.5, loc="best")
        for idx in range(n_p, len(axes_flat)):
            axes_flat[idx].set_visible(False)

        fig.suptitle(title, fontsize=12.5, y=0.985)
        fig.tight_layout(rect=[0, 0.01, 1, 0.955])
        for fmt in formats:
            p_out = out_dir / f"{out_base}.{fmt}"
            fig.savefig(p_out, bbox_inches="tight")
            generated.append(p_out)
        plt.close(fig)

    # 2. Combined overview of all fitted resonance shape & global phase parameters vs -t
    shape_phase_entries = [
        e for e in param_entries if e["role"] in ("mass", "width", "coupling_g", "phase")
    ]
    if not shape_phase_entries:
        shape_phase_entries = param_entries
    _save_param_grid(
        shape_phase_entries,
        f"{proj_name} - Fitted Resonance Shape & Phase Parameters vs $-t$",
        "resonance_parameters_vs_t",
    )

    # 3. Per-resonance complete parameter evolution grid (shape, phase, and wave couplings k)
    by_resonance: Dict[str, List[Dict[str, Any]]] = {}
    for entry in param_entries:
        by_resonance.setdefault(entry["res_key"], []).append(entry)

    for res_key, r_entries in by_resonance.items():
        if res_key == "resonance_other":
            continue
        res_disp = r_entries[0]["res_display"]
        safe_rk = re.sub(r"[^a-zA-Z0-9_\-\.]", "_", res_key)
        _save_param_grid(
            r_entries,
            f"{proj_name} - Fitted {res_disp} Parameters vs $-t$",
            f"resonance_{safe_rk}_params_vs_t",
        )

    return generated


def plot_spd_total_intensity_vs_t(
    config: FitConfig,
    t_data: Dict[str, Dict[str, Any]],
    sorted_tbins: List[str],
    magnitudes: List[str],
    out_dir: Path | str,
    fit_mass_step: float = 0.04,
    formats: Tuple[str, ...] = ("png", "pdf"),
    style: Optional[str] = None,
) -> List[Path]:
    """
    Plot total and reflectivity-separated magnitude squared for S, P, D (and any
    active higher-L) waves across Mandelstam -t:
      1. <out_dir>/spd_waves_total_reflectivity_vs_t.<fmt>:
         Top panel: Integrated sum |V|^2 vs -t for S, P, D, Total (epsilon = +-1 summed) and H_0(0,0).
         Bottom panel: Wave intensity fractions (%) vs -t.
      2. <out_dir>/spd_waves_vs_mass_all_tbins.<fmt>:
         1 x N_tbins mass spectrum comparing Total, S, P, D (epsilon = +-1 summed).
      3. <out_dir>/spd_waves_refl_separated_vs_mass_all_tbins.<fmt>:
         2 x N_tbins mass spectrum separated by reflectivity: top row natural parity
         (epsilon = +1: S+, P+, D+), bottom row unnatural parity (epsilon = -1: S-, P-, D-),
         overlaid on H_0(0,0) data and Total fit curve with shared y-axis.
      4. <out_dir>/spd_waves_refl_separated_vs_t.<fmt>:
         Integrated reflectivity-separated wave yields (S+, P+, D+ vs S-, P-, D-)
         and fractional shares of total intensity vs -t.
    """
    generated: List[Path] = []
    if len(sorted_tbins) < 2 or not magnitudes:
        return generated

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return generated

    if style:
        apply_plot_style(style)

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    show_grid = should_show_grid(config)

    l_vals = sorted(set(extract_amplitude_l(k) for k in magnitudes if extract_amplitude_l(k) is not None))
    if not l_vals:
        return generated

    t_centers, t_hw, has_numeric_t = _extract_t_axis_coords(t_data, sorted_tbins)
    t_labels = [t_data[t].get("name", t) for t in sorted_tbins]

    l_colors = {
        0: "#0072B2",  # S-wave: Royal Blue
        1: "#D55E00",  # P-wave: Vermilion
        2: "#009E73",  # D-wave: Pine Green
        3: "#785EF0",
        4: "#E69F00",
    }
    l_markers = {
        0: "s",
        1: "^",
        2: "D",
        3: "P",
        4: "X",
    }
    l_linestyles = {
        0: "--",
        1: "-.",
        2: ":",
        3: (0, (3, 1, 1, 1)),
        4: (0, (5, 2)),
    }

    def _compute_integrated_wave_stats(wave_subset: List[str], t_info: Dict[str, Any]) -> Dict[str, Optional[float]]:
        masses = t_info.get("masses", [])
        nom_dict = t_info.get("nominal", {})
        nom_total = sum(
            sum(nom_dict[w][i] ** 2 for w in wave_subset if w in nom_dict and i < len(nom_dict[w]))
            for i in range(len(masses))
        )

        mcmc_dict = t_info.get("mcmc", {})
        avail_w = [w for w in wave_subset if w in mcmc_dict and mcmc_dict[w]]
        med_v: Optional[float] = None
        lo_e: Optional[float] = None
        hi_e: Optional[float] = None
        sample_totals: Optional[List[float]] = None
        sample_weights: Optional[List[float]] = None

        if avail_w and masses:
            per_mass_samps: List[List[np.ndarray]] = []
            w_arr_list: List[Optional[np.ndarray]] = []
            mc_weights = t_info.get("mcmc_weights", {})
            for m in masses:
                m_r = round(float(m), 4)
                ws_at_m: List[np.ndarray] = []
                w_at_m: Optional[np.ndarray] = None
                for w in avail_w:
                    s = mcmc_dict[w].get(m_r)
                    wt = mc_weights.get(w, {}).get(m_r)
                    if s is None:
                        for mk, mv in mcmc_dict[w].items():
                            if abs(mk - m_r) < 0.005:
                                s = mv
                                wt = mc_weights.get(w, {}).get(mk)
                                break
                    if s:
                        ws_at_m.append(np.asarray(s, dtype=float))
                        if wt is not None and w_at_m is None:
                            w_at_m = np.asarray(wt, dtype=float)
                if ws_at_m:
                    per_mass_samps.append(ws_at_m)
                    w_arr_list.append(w_at_m)

            if per_mass_samps:
                n_samps = min(min(len(arr) for arr in row) for row in per_mass_samps)
                if n_samps > 0:
                    tot_vec = np.zeros(n_samps, dtype=float)
                    for row in per_mass_samps:
                        for arr in row:
                            tot_vec += arr[:n_samps] ** 2
                    sample_totals = tot_vec.tolist()
                    w_ref = next((wa[:n_samps] for wa in w_arr_list if wa is not None and len(wa) >= n_samps), None)
                    if w_ref is not None:
                        sample_weights = w_ref.tolist()
                    q16 = compute_weighted_percentile(tot_vec, 16.0, w_ref)
                    q50 = compute_weighted_percentile(tot_vec, 50.0, w_ref)
                    q84 = compute_weighted_percentile(tot_vec, 84.0, w_ref)
                    med_v = q50
                    lo_e = max(0.0, q50 - q16)
                    hi_e = max(0.0, q84 - q50)

        return {
            "nom": float(nom_total),
            "med": med_v,
            "lo": lo_e,
            "hi": hi_e,
            "samples": sample_totals,
            "weights": sample_weights,
        }

    # Gather integrated quantities across all t-bins
    waves_by_l = {l_v: [k for k in magnitudes if extract_amplitude_l(k) == l_v] for l_v in l_vals}
    tot_stats: List[Dict[str, Any]] = []
    l_stats: Dict[int, List[Dict[str, Any]]] = {l_v: [] for l_v in l_vals}
    h0_tot_vals: List[Optional[float]] = []
    h0_tot_errs: List[Optional[float]] = []

    for t_name in sorted_tbins:
        t_info = t_data[t_name]
        masses = t_info.get("masses", [])
        tot_s = _compute_integrated_wave_stats(magnitudes, t_info)
        tot_stats.append(tot_s)
        for l_v in l_vals:
            l_stats[l_v].append(_compute_integrated_wave_stats(waves_by_l[l_v], t_info))

        # Sum H_0(0,0) across fitted mass window
        h0_d = t_info.get("h0_data")
        if h0_d and h0_d.get("mass") and h0_d.get("val") and masses:
            m_min_w = min(masses) - 0.5 * fit_mass_step
            m_max_w = max(masses) + 0.5 * fit_mass_step
            h0_sum = 0.0
            h0_var = 0.0
            n_h0 = 0
            errs_list = h0_d.get("err", [0.0] * len(h0_d["val"]))
            for mh, yh, eh in zip(h0_d["mass"], h0_d["val"], errs_list):
                if m_min_w <= mh <= m_max_w:
                    h0_sum += float(yh)
                    h0_var += float(eh) ** 2
                    n_h0 += 1
            if n_h0 > 0:
                h0_tot_vals.append(h0_sum)
                h0_tot_errs.append(math.sqrt(h0_var))
            else:
                h0_tot_vals.append(None)
                h0_tot_errs.append(None)
        else:
            h0_tot_vals.append(None)
            h0_tot_errs.append(None)

    # --- Figure 1: 1D Integrated Total Magnitude Squared (S, P, D in Total Reflectivities) + Wave Fractions vs -t ---
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1,
        figsize=(9.0, 7.2),
        dpi=300,
        sharex=True,
        gridspec_kw={"height_ratios": [1.8, 1.0], "hspace": 0.08, "top": 0.92, "bottom": 0.10, "left": 0.12, "right": 0.96},
    )

    # Overlay H0(0,0) data sum if present across t-bins
    valid_h0_idx = [i for i, hv in enumerate(h0_tot_vals) if hv is not None]
    if valid_h0_idx:
        hx = [t_centers[i] for i in valid_h0_idx]
        hhw = [t_hw[i] for i in valid_h0_idx]
        hy = [h0_tot_vals[i] for i in valid_h0_idx]
        he = [h0_tot_errs[i] or 0.0 for i in valid_h0_idx]
        ax_top.errorbar(
            hx, hy, yerr=he, xerr=hhw if has_numeric_t else None,
            fmt="o", markersize=6.0, mfc="none", mec="#1A1A1A", ecolor="#555555",
            elinewidth=1.2, capsize=2.8, label=r"Data Yield $\sum H_0(0,0)$", zorder=3,
        )

    # Total intensity curve across all waves & reflectivities
    y_tot_nom = [s["nom"] for s in tot_stats]
    has_mcmc_tot = any(s["med"] is not None for s in tot_stats)
    tot_yerr = (
        [[s["lo"] or 0.0 for s in tot_stats], [s["hi"] or 0.0 for s in tot_stats]]
        if has_mcmc_tot else None
    )
    ax_top.plot(t_centers, y_tot_nom, "-", color="#1A1A1A", linewidth=1.8, zorder=4)
    ax_top.errorbar(
        t_centers, y_tot_nom, yerr=tot_yerr, xerr=t_hw if has_numeric_t else None,
        fmt="o", markersize=5.2, color="#1A1A1A", ecolor="#1A1A1A", elinewidth=1.3,
        capsize=2.8, label=r"Total Fit $\sum |V|^2$ ($\epsilon = \pm 1$)", zorder=5,
    )

    # S, P, D wave curves in total reflectivities
    for i_l, l_v in enumerate(l_vals):
        l_sym = L_SPECTROSCOPIC.get(l_v, f"L{l_v}")
        c_l = l_colors.get(l_v, PHYSICS_PALETTE[(i_l + 3) % len(PHYSICS_PALETTE)])
        mk_l = l_markers.get(l_v, "o")
        ls_l = l_linestyles.get(l_v, "-")

        y_l_nom = [s["nom"] for s in l_stats[l_v]]
        has_mc_l = any(s["med"] is not None for s in l_stats[l_v])
        y_l_err = (
            [[s["lo"] or 0.0 for s in l_stats[l_v]], [s["hi"] or 0.0 for s in l_stats[l_v]]]
            if has_mc_l else None
        )
        if l_v == 0:
            lbl_l = r"$S$-wave ($|S^+_0|^2 + |S^-_0|^2$)"
            short_lbl = r"$S$-wave"
        else:
            lbl_l = rf"${l_sym}$-wave ($\sum_{{m,\pm}} |{l_sym}^{{\pm}}_{{m}}|^2$)"
            short_lbl = rf"${l_sym}$-wave"

        dx_l = (i_l - 0.5 * (len(l_vals) - 1)) * 0.004 if (has_mcmc_tot and has_numeric_t) else 0.0
        xl = [xc + dx_l for xc in t_centers]

        ax_top.plot(xl, y_l_nom, linestyle=ls_l, color=c_l, linewidth=1.8, zorder=5 + i_l)
        ax_top.errorbar(
            xl, y_l_nom, yerr=y_l_err, xerr=t_hw if has_numeric_t else None,
            fmt=mk_l, markersize=5.5, color=c_l, ecolor=c_l, elinewidth=1.3,
            capsize=2.6, label=lbl_l, zorder=6 + i_l,
        )

        # Wave fraction (%) on bottom panel
        frac_nom = [
            (100.0 * yl / yt) if yt > 0 else 0.0
            for yl, yt in zip(y_l_nom, y_tot_nom)
        ]
        frac_lo: List[float] = []
        frac_hi: List[float] = []
        has_frac_err = False
        for idx_t in range(len(sorted_tbins)):
            sl = l_stats[l_v][idx_t]["samples"]
            st = tot_stats[idx_t]["samples"]
            wt = tot_stats[idx_t]["weights"]
            if sl and st:
                n_s = min(len(sl), len(st))
                f_arr = np.array(
                    [100.0 * sl[j] / st[j] if st[j] > 0 else 0.0 for j in range(n_s)],
                    dtype=float,
                )
                w_arr = np.asarray(wt[:n_s], dtype=float) if (wt is not None and len(wt) >= n_s) else None
                if np.std(f_arr) > 1e-9:
                    has_frac_err = True
                    q16 = compute_weighted_percentile(f_arr, 16.0, w_arr)
                    q50 = compute_weighted_percentile(f_arr, 50.0, w_arr)
                    q84 = compute_weighted_percentile(f_arr, 84.0, w_arr)
                    frac_lo.append(max(0.0, q50 - q16))
                    frac_hi.append(max(0.0, q84 - q50))
                else:
                    frac_lo.append(0.0)
                    frac_hi.append(0.0)
            else:
                frac_lo.append(0.0)
                frac_hi.append(0.0)

        ax_bot.plot(xl, frac_nom, linestyle=ls_l, color=c_l, linewidth=1.6, zorder=4 + i_l)
        ax_bot.errorbar(
            xl, frac_nom,
            yerr=[frac_lo, frac_hi] if has_frac_err else None,
            xerr=t_hw if has_numeric_t else None,
            fmt=mk_l, markersize=5.0, color=c_l, ecolor=c_l, elinewidth=1.2,
            capsize=2.4, label=short_lbl, zorder=5 + i_l,
        )

    ax_top.set_ylabel(r"Integrated Magnitude Squared $\sum |V|^2$", fontsize=11.0)
    ax_top.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)
    ax_top.set_ylim(0, None)
    ax_top.legend(frameon=True, fontsize=8.8, loc="best")
    proj_name = config.project.name if config.project else "PartialWaveFit"
    ax_top.set_title(
        f"{proj_name} - Total Magnitude Squared by Wave ($\\epsilon = \\pm 1$ Summed) vs $-t$",
        fontsize=12.0, pad=8,
    )
    if show_grid:
        ax_top.grid(True, linestyle=":", alpha=0.45)

    ax_bot.set_ylabel("Wave Fraction [%]", fontsize=11.0)
    ax_bot.set_ylim(-3, 108)
    ax_bot.set_yticks([0, 20, 40, 60, 80, 100])
    if has_numeric_t and t_centers:
        t_left = min(c - hw for c, hw in zip(t_centers, t_hw))
        t_right = max(c + hw for c, hw in zip(t_centers, t_hw))
        t_span = max(t_right - t_left, 0.1)
        ax_bot.set_xlim(t_left - 0.04 * t_span, t_right + 0.04 * t_span)
        ax_bot.set_xlabel(r"$-t\ [\mathrm{GeV}^2]$", fontsize=11.5)
    elif t_centers:
        ax_bot.set_xticks(t_centers)
        ax_bot.set_xticklabels(t_labels, fontsize=9.0)
        ax_bot.set_xlabel(r"$-t$ Bin", fontsize=11.5)
    ax_bot.legend(frameon=True, fontsize=8.5, ncol=min(4, len(l_vals)), loc="best")
    if show_grid:
        ax_bot.grid(True, linestyle=":", alpha=0.45)

    for fmt in formats:
        p_spd = out_path / f"spd_waves_total_reflectivity_vs_t.{fmt}"
        fig.savefig(p_spd, bbox_inches="tight")
        generated.append(p_spd)
    plt.close(fig)

    # --- Figure 2 & 3 shared helpers ---
    def _compute_mass_sum_and_band(
        wave_subset: List[str],
        t_info: Dict[str, Any],
    ) -> Tuple[List[float], Optional[List[float]], Optional[List[float]]]:
        masses = t_info.get("masses", [])
        nom_d = t_info.get("nominal", {})
        y_nom = [
            sum(nom_d[w][i] ** 2 for w in wave_subset if w in nom_d and i < len(nom_d[w]))
            for i in range(len(masses))
        ]
        mcmc_dict = t_info.get("mcmc", {})
        mc_weights = t_info.get("mcmc_weights", {})
        avail_w = [w for w in wave_subset if w in mcmc_dict and mcmc_dict[w]]
        if not (avail_w and masses):
            return y_nom, None, None

        lo_band: List[float] = []
        hi_band: List[float] = []
        has_band = False
        for i, m in enumerate(masses):
            m_r = round(float(m), 4)
            ws_at_m: List[np.ndarray] = []
            w_at_m: Optional[np.ndarray] = None
            for w in avail_w:
                s = mcmc_dict[w].get(m_r)
                wt = mc_weights.get(w, {}).get(m_r)
                if s is None:
                    for mk, mv in mcmc_dict[w].items():
                        if abs(mk - m_r) < 0.005:
                            s = mv
                            wt = mc_weights.get(w, {}).get(mk)
                            break
                if s:
                    ws_at_m.append(np.asarray(s, dtype=float))
                    if wt is not None and w_at_m is None:
                        w_at_m = np.asarray(wt, dtype=float)
            if ws_at_m:
                n_s = min(len(arr) for arr in ws_at_m)
                if n_s >= 2:
                    tot_m_vec = np.zeros(n_s, dtype=float)
                    for arr in ws_at_m:
                        tot_m_vec += arr[:n_s] ** 2
                    w_ref = w_at_m[:n_s] if (w_at_m is not None and len(w_at_m) >= n_s) else None
                    q16 = compute_weighted_percentile(tot_m_vec, 16.0, w_ref)
                    q84 = compute_weighted_percentile(tot_m_vec, 84.0, w_ref)
                    if q84 - q16 > 1e-9:
                        has_band = True
                    lo_band.append(max(0.0, q16))
                    hi_band.append(max(0.0, q84))
                    continue
            v_fallback = y_nom[i] if i < len(y_nom) else 0.0
            lo_band.append(v_fallback)
            hi_band.append(v_fallback)

        if not has_band:
            return y_nom, None, None
        return y_nom, lo_band, hi_band

    def _draw_spd_mass_curve(
        ax: Any,
        masses: List[float],
        y_vals: List[float],
        lo_band: Optional[List[float]],
        hi_band: Optional[List[float]],
        color: str,
        linestyle: Any,
        label: Optional[str],
        is_md: bool,
        zorder: int,
        linewidth: float = 1.8,
    ) -> None:
        if not masses or not y_vals:
            return
        if is_md and len(masses) >= 4:
            ms, ys = smooth_curve(masses, y_vals, num_points=200)
            ax.plot(ms, ys, linestyle=linestyle, color=color, linewidth=linewidth, label=label, zorder=zorder)
            if lo_band is not None and hi_band is not None:
                _, y_lo_s = smooth_curve(masses, lo_band, num_points=200)
                _, y_hi_s = smooth_curve(masses, hi_band, num_points=200)
                ax.fill_between(ms, np.maximum(0.0, y_lo_s), y_hi_s, color=color, alpha=0.18, zorder=zorder - 1)
        else:
            edges = compute_bin_edges(masses, default_step=fit_mass_step)
            if hasattr(ax, "stairs"):
                ax.stairs(y_vals, edges, linestyle=linestyle, color=color, linewidth=linewidth, label=label, zorder=zorder)
            else:
                ax.step(masses, y_vals, where="mid", linestyle=linestyle, color=color, linewidth=linewidth, label=label, zorder=zorder)
            if lo_band is not None and hi_band is not None:
                ax.fill_between(edges[:-1], lo_band, hi_band, step="post", color=color, alpha=0.18, zorder=zorder - 1)

    # --- Figure 2: Clean 1 x N_tbins Mass-Spectrum Comparison of S, P, D (Total Reflectivities) Across -t ---
    n_cols = len(sorted_tbins)
    fig_width = max(3.8 * n_cols, 8.0)
    fig_height = 4.2
    sharey_cfg = getattr(config.plotting, "amplitudes_vs_t_sharey", False) if getattr(config, "plotting", None) else False

    fig2, axes2 = plt.subplots(
        1, n_cols, figsize=(fig_width, fig_height), sharey=sharey_cfg, dpi=300,
    )
    axes2_list = list(np.atleast_1d(axes2))

    for idx, t_name in enumerate(sorted_tbins):
        ax = axes2_list[idx]
        t_info = t_data[t_name]
        t_lbl = t_info.get("label", t_info.get("name", t_name))
        masses = t_info.get("masses", [])
        if not masses:
            continue

        if t_info.get("h0_data"):
            h0_d = t_info["h0_data"]
            if h0_d.get("mass") and h0_d.get("val"):
                ax.errorbar(
                    h0_d["mass"], h0_d["val"], yerr=h0_d.get("err", [0.0] * len(h0_d["val"])),
                    fmt="o", color="black", markersize=2.5, alpha=0.80,
                    capsize=1.0, elinewidth=0.8,
                    label=r"$H_0(0,0)$ (Total)" if idx == 0 else None, zorder=2,
                )

        tot_m, tot_lo_m, tot_hi_m = _compute_mass_sum_and_band(magnitudes, t_info)
        _draw_spd_mass_curve(
            ax, masses, tot_m, tot_lo_m, tot_hi_m,
            color="#1A1A1A", linestyle="-",
            label=r"Total $\sum |V|^2$" if idx == 0 else None,
            is_md=(len(masses) >= 4), zorder=4, linewidth=2.0,
        )

        for i_l, l_v in enumerate(l_vals):
            l_sym = L_SPECTROSCOPIC.get(l_v, f"L{l_v}")
            c_l = l_colors.get(l_v, PHYSICS_PALETTE[(i_l + 3) % len(PHYSICS_PALETTE)])
            ls_l = l_linestyles.get(l_v, "--")
            w_sub = waves_by_l[l_v]
            y_l_m, lo_l_m, hi_l_m = _compute_mass_sum_and_band(w_sub, t_info)
            lbl_m = (
                r"Total $S$-wave ($|S|^2_{\pm}$)" if l_v == 0
                else rf"Total ${l_sym}$-wave ($\sum |{l_sym}^\pm|^2$)"
            )
            is_md = all(is_amplitude_mass_dependent(w, config) for w in w_sub)
            _draw_spd_mass_curve(
                ax, masses, y_l_m, lo_l_m, hi_l_m,
                color=c_l, linestyle=ls_l,
                label=lbl_m if idx == 0 else None,
                is_md=is_md, zorder=5 + i_l,
            )

        ax.set_xlim(masses[0] - fit_mass_step * 0.6, masses[-1] + fit_mass_step * 0.6)
        ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=10.5)
        ax.set_title(t_lbl, fontsize=11, pad=6)
        if idx == 0:
            ax.set_ylabel("Intensity $|V|^2$ (Total Reflectivity)", fontsize=10.5)
        if show_grid:
            ax.grid(True, linestyle=":", alpha=0.45)

    handles2, labels2 = axes2_list[0].get_legend_handles_labels()
    if handles2:
        fig2.legend(
            handles2, labels2, loc="upper center", bbox_to_anchor=(0.5, 0.03),
            ncol=min(5, len(handles2)), frameon=True, fontsize=8.8,
        )
    fig2.suptitle(f"{proj_name} - Total Magnitude Squared by Wave ($\\epsilon = \\pm 1$) vs $-t$", fontsize=12.0, y=0.98)
    fig2.tight_layout(rect=[0, 0.06, 1, 0.94])
    for fmt in formats:
        p_spd_m = out_path / f"spd_waves_vs_mass_all_tbins.{fmt}"
        fig2.savefig(p_spd_m, bbox_inches="tight")
        generated.append(p_spd_m)
    plt.close(fig2)

    # --- Figure 3: 2 x N_tbins Reflectivity-Separated Wave Sum vs Mass Across -t (Natural vs Unnatural Parity) ---
    refl_rows: List[Tuple[str, str, str]] = []
    if any(k.startswith("a_") for k in magnitudes):
        refl_rows.append(("a_", "+", r"Positive Reflectivity ($\epsilon = +1$)"))
    if any(k.startswith("b_") for k in magnitudes):
        refl_rows.append(("b_", "-", r"Negative Reflectivity ($\epsilon = -1$)"))

    if refl_rows:
        refl_wave_colors = {
            0: "#009E73",  # S-wave: Teal Green
            1: "#0072B2",  # P-wave: Royal Blue
            2: "#CC79A7",  # D-wave: Magenta/Rose
            3: "#D55E00",
            4: "#E69F00",
        }
        refl_wave_linestyles = {
            0: "-",
            1: "-",
            2: "-",
            3: "--",
            4: "-.",
        }

        n_r = len(refl_rows)
        fig3, axes3 = plt.subplots(
            n_r, n_cols,
            figsize=(max(3.6 * n_cols + 1.4, 8.6), 3.6 * n_r + 0.7),
            sharex=True, sharey=True, dpi=300,
        )
        axes3_arr = np.atleast_2d(axes3)
        if axes3_arr.shape != (n_r, n_cols):
            axes3_arr = axes3_arr.reshape(n_r, n_cols)

        global_ymax = 0.0
        for t_name in sorted_tbins:
            t_info = t_data[t_name]
            masses = t_info.get("masses", [])
            if masses:
                tot_vals, _, _ = _compute_mass_sum_and_band(magnitudes, t_info)
                if tot_vals:
                    global_ymax = max(global_ymax, max(tot_vals))
            h0_d = t_info.get("h0_data")
            if h0_d and h0_d.get("val"):
                global_ymax = max(global_ymax, max(float(v) for v in h0_d["val"]))
        y_top_lim = global_ymax * 1.15 if global_ymax > 0 else 1.0

        for r_idx, (prefix, refl_sign, refl_desc) in enumerate(refl_rows):
            for c_idx, t_name in enumerate(sorted_tbins):
                ax = axes3_arr[r_idx, c_idx]
                t_info = t_data[t_name]
                t_lbl = t_info.get("label", t_info.get("name", t_name))
                masses = t_info.get("masses", [])
                if not masses:
                    continue

                is_first_col = (c_idx == 0)
                if t_info.get("h0_data"):
                    h0_d = t_info["h0_data"]
                    if h0_d.get("mass") and h0_d.get("val"):
                        ax.errorbar(
                            h0_d["mass"], h0_d["val"], yerr=h0_d.get("err", [0.0] * len(h0_d["val"])),
                            fmt="o", color="#1A1A1A", markersize=2.5, alpha=0.85,
                            capsize=1.0, elinewidth=0.8,
                            label=r"$H_0(0,0)$ Data" if is_first_col else None, zorder=2,
                        )

                tot_m, tot_lo_m, tot_hi_m = _compute_mass_sum_and_band(magnitudes, t_info)
                _draw_spd_mass_curve(
                    ax, masses, tot_m, tot_lo_m, tot_hi_m,
                    color="#D62728", linestyle="-",
                    label=r"Total $\sum |V|^2$" if is_first_col else None,
                    is_md=(len(masses) >= 4), zorder=4, linewidth=1.8,
                )

                for i_l, l_v in enumerate(l_vals):
                    w_sub_r = [w for w in waves_by_l[l_v] if w.startswith(prefix)]
                    if not w_sub_r:
                        continue
                    l_sym = L_SPECTROSCOPIC.get(l_v, f"L{l_v}")
                    c_r = refl_wave_colors.get(l_v, PHYSICS_PALETTE[(i_l + 3) % len(PHYSICS_PALETTE)])
                    ls_r = refl_wave_linestyles.get(l_v, "-")
                    y_r, lo_r, hi_r = _compute_mass_sum_and_band(w_sub_r, t_info)
                    lbl_r = (
                        rf"${l_sym}^{{{refl_sign}}}$ ($|{l_sym}^{{{refl_sign}}}_0|^2$)" if l_v == 0
                        else rf"${l_sym}^{{{refl_sign}}}$ ($\sum_m |{l_sym}^{{{refl_sign}}}_m|^2$)"
                    )
                    is_md = all(is_amplitude_mass_dependent(w, config) for w in w_sub_r)
                    _draw_spd_mass_curve(
                        ax, masses, y_r, lo_r, hi_r,
                        color=c_r, linestyle=ls_r,
                        label=lbl_r if is_first_col else None,
                        is_md=is_md, zorder=5 + i_l, linewidth=1.8,
                    )

                ax.set_xlim(masses[0] - fit_mass_step * 0.6, masses[-1] + fit_mass_step * 0.6)
                ax.set_ylim(0.0, y_top_lim)
                if r_idx == 0:
                    ax.set_title(t_lbl, fontsize=11.0, pad=6)
                if r_idx == n_r - 1:
                    ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=10.5)
                if c_idx == 0:
                    ax.set_ylabel(f"Intensity $|V|^2$\n{refl_desc}", fontsize=10.0)
                    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)
                if show_grid:
                    ax.grid(True, linestyle=":", alpha=0.45)

            row_handles, row_labels = axes3_arr[r_idx, 0].get_legend_handles_labels()
            if row_handles:
                axes3_arr[r_idx, -1].legend(
                    row_handles, row_labels,
                    loc="center left", bbox_to_anchor=(1.02, 0.5),
                    frameon=True, fontsize=8.8,
                )

        fig3.suptitle(
            f"{proj_name} - Reflectivity-Separated Total Wave Intensity vs Mass Across $-t$",
            fontsize=12.5, y=0.985,
        )
        fig3.tight_layout(rect=[0, 0.01, 0.99, 0.95], w_pad=0.6, h_pad=0.8)
        for fmt in formats:
            p_refl_m = out_path / f"spd_waves_refl_separated_vs_mass_all_tbins.{fmt}"
            fig3.savefig(p_refl_m, bbox_inches="tight")
            generated.append(p_refl_m)
        plt.close(fig3)

        # --- Figure 4: Integrated Reflectivity-Separated Wave Intensities (S+, P+, D+ vs S-, P-, D-) vs -t ---
        n_r = len(refl_rows)
        fig4, axes4 = plt.subplots(
            2, n_r,
            figsize=(6.2 * n_r + 0.6, 7.2),
            dpi=300,
            sharex="col", sharey="row",
            gridspec_kw={"height_ratios": [1.8, 1.0], "hspace": 0.09, "wspace": 0.08, "top": 0.90, "bottom": 0.10, "left": 0.10, "right": 0.97},
        )
        axes4_arr = np.atleast_2d(axes4)
        if axes4_arr.shape != (2, n_r):
            axes4_arr = axes4_arr.reshape(2, n_r)

        for r_idx, (prefix, refl_sign, refl_desc) in enumerate(refl_rows):
            ax_t = axes4_arr[0, r_idx]
            ax_b = axes4_arr[1, r_idx]

            refl_waves_all = [k for k in magnitudes if k.startswith(prefix)]
            refl_tot_stats = [_compute_integrated_wave_stats(refl_waves_all, t_data[t]) for t in sorted_tbins]
            y_rtot = [s["nom"] for s in refl_tot_stats]
            has_mc_rt = any(s["med"] is not None for s in refl_tot_stats)
            y_rtot_err = (
                [[s["lo"] or 0.0 for s in refl_tot_stats], [s["hi"] or 0.0 for s in refl_tot_stats]]
                if has_mc_rt else None
            )

            # Subtle reference of total fit intensity summed over both reflectivities
            ax_t.plot(
                t_centers, y_tot_nom, "--", color="#888888", linewidth=1.2, alpha=0.80,
                label=r"Total Fit $\sum |V|^2$ ($\epsilon = \pm 1$)", zorder=3,
            )
            ax_t.plot(t_centers, y_rtot, "-", color="#1A1A1A", linewidth=1.8, zorder=4)
            ax_t.errorbar(
                t_centers, y_rtot, yerr=y_rtot_err, xerr=t_hw if has_numeric_t else None,
                fmt="o", markersize=5.0, color="#1A1A1A", ecolor="#1A1A1A", elinewidth=1.2,
                capsize=2.6, label=rf"Total $\sum |V^{{{refl_sign}}}|^2$", zorder=5,
            )

            for i_l, l_v in enumerate(l_vals):
                w_sub_r = [w for w in waves_by_l[l_v] if w.startswith(prefix)]
                if not w_sub_r:
                    continue
                l_sym = L_SPECTROSCOPIC.get(l_v, f"L{l_v}")
                c_l = refl_wave_colors.get(l_v, PHYSICS_PALETTE[(i_l + 3) % len(PHYSICS_PALETTE)])
                mk_l = l_markers.get(l_v, "o")
                ls_l = l_linestyles.get(l_v, "--")

                stats_rl = [_compute_integrated_wave_stats(w_sub_r, t_data[t]) for t in sorted_tbins]
                y_rl = [s["nom"] for s in stats_rl]
                has_mc_rl = any(s["med"] is not None for s in stats_rl)
                y_rl_err = (
                    [[s["lo"] or 0.0 for s in stats_rl], [s["hi"] or 0.0 for s in stats_rl]]
                    if has_mc_rl else None
                )
                lbl_rl = (
                    rf"${l_sym}^{{{refl_sign}}}$ ($|{l_sym}^{{{refl_sign}}}_0|^2$)" if l_v == 0
                    else rf"${l_sym}^{{{refl_sign}}}$ ($\sum_m |{l_sym}^{{{refl_sign}}}_m|^2$)"
                )
                short_rl = rf"${l_sym}^{{{refl_sign}}}$"
                dx_l = (i_l - 0.5 * (len(l_vals) - 1)) * 0.004 if (has_mcmc_tot and has_numeric_t) else 0.0
                xl = [xc + dx_l for xc in t_centers]

                ax_t.plot(xl, y_rl, linestyle=ls_l, color=c_l, linewidth=1.8, zorder=5 + i_l)
                ax_t.errorbar(
                    xl, y_rl, yerr=y_rl_err, xerr=t_hw if has_numeric_t else None,
                    fmt=mk_l, markersize=5.2, color=c_l, ecolor=c_l, elinewidth=1.2,
                    capsize=2.5, label=lbl_rl, zorder=6 + i_l,
                )

                # Fraction of total (both-reflectivity) intensity vs -t
                frac_r_nom = [(100.0 * yl / yt) if yt > 0 else 0.0 for yl, yt in zip(y_rl, y_tot_nom)]
                frac_r_lo: List[float] = []
                frac_r_hi: List[float] = []
                has_fr_err = False
                for idx_t in range(len(sorted_tbins)):
                    sl = stats_rl[idx_t]["samples"]
                    st = tot_stats[idx_t]["samples"]
                    wt = tot_stats[idx_t]["weights"]
                    if sl and st:
                        n_s = min(len(sl), len(st))
                        f_arr = np.array(
                            [100.0 * sl[j] / st[j] if st[j] > 0 else 0.0 for j in range(n_s)],
                            dtype=float,
                        )
                        w_arr = np.asarray(wt[:n_s], dtype=float) if (wt is not None and len(wt) >= n_s) else None
                        if np.std(f_arr) > 1e-9:
                            has_fr_err = True
                            q16 = compute_weighted_percentile(f_arr, 16.0, w_arr)
                            q50 = compute_weighted_percentile(f_arr, 50.0, w_arr)
                            q84 = compute_weighted_percentile(f_arr, 84.0, w_arr)
                            frac_r_lo.append(max(0.0, q50 - q16))
                            frac_r_hi.append(max(0.0, q84 - q50))
                        else:
                            frac_r_lo.append(0.0)
                            frac_r_hi.append(0.0)
                    else:
                        frac_r_lo.append(0.0)
                        frac_r_hi.append(0.0)

                ax_b.plot(xl, frac_r_nom, linestyle=ls_l, color=c_l, linewidth=1.6, zorder=4 + i_l)
                ax_b.errorbar(
                    xl, frac_r_nom,
                    yerr=[frac_r_lo, frac_r_hi] if has_fr_err else None,
                    xerr=t_hw if has_numeric_t else None,
                    fmt=mk_l, markersize=4.8, color=c_l, ecolor=c_l, elinewidth=1.1,
                    capsize=2.3, label=short_rl, zorder=5 + i_l,
                )

            ax_t.set_title(refl_desc, fontsize=11.0, pad=7)
            ax_t.set_ylim(0, None)
            ax_t.legend(frameon=True, fontsize=8.4, loc="best")
            if r_idx == 0:
                ax_t.set_ylabel(r"Integrated Magnitude Squared $\sum |V|^2$", fontsize=10.5)
                ax_t.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)
                ax_b.set_ylabel("Share of Total Fit [%]", fontsize=10.5)
            if show_grid:
                ax_t.grid(True, linestyle=":", alpha=0.45)
                ax_b.grid(True, linestyle=":", alpha=0.45)

            ax_b.set_ylim(-3, 108)
            ax_b.set_yticks([0, 20, 40, 60, 80, 100])
            ax_b.legend(frameon=True, fontsize=8.2, ncol=min(4, len(l_vals)), loc="best")
            if has_numeric_t and t_centers:
                t_left = min(c - hw for c, hw in zip(t_centers, t_hw))
                t_right = max(c + hw for c, hw in zip(t_centers, t_hw))
                t_span = max(t_right - t_left, 0.1)
                ax_b.set_xlim(t_left - 0.04 * t_span, t_right + 0.04 * t_span)
                ax_b.set_xlabel(r"$-t\ [\mathrm{GeV}^2]$", fontsize=11.0)
            elif t_centers:
                ax_b.set_xticks(t_centers)
                ax_b.set_xticklabels(t_labels, fontsize=9.0)
                ax_b.set_xlabel(r"$-t$ Bin", fontsize=11.0)

        fig4.suptitle(
            f"{proj_name} - Reflectivity-Separated Total Magnitude Squared by Wave vs $-t$",
            fontsize=12.0, y=0.975,
        )
        for fmt in formats:
            p_refl_t = out_path / f"spd_waves_refl_separated_vs_t.{fmt}"
            fig4.savefig(p_refl_t, bbox_inches="tight")
            generated.append(p_refl_t)
        plt.close(fig4)

    return generated


def plot_amplitudes_across_tbins(
    config: FitConfig,
    tbins: Optional[List[str]] = None,
    run_dir: Optional[str] = None,
    output_dir: Optional[Path | str] = None,
    formats: Tuple[str, ...] = ("png", "pdf"),
    style: Optional[str] = None,
    sharey: bool = False,
) -> List[Path]:
    """
    Plot fitted partial wave amplitudes, phases, and total intensity across all t-bins
    side-by-side (1 x N_tbins subplots) ordered from lowest to highest Mandelstam -t.

    Generates:
      - Individual wave intensity 1 x N_tbins plots:
          <out_dir>/individual/<wave>_intensity_vs_t.<fmt>
      - Individual wave phase 1 x N_tbins plots:
          <out_dir>/individual/<wave>_phase_vs_t.<fmt>
      - Summed wave totals 1 x N_tbins plots:
          <out_dir>/total_intensity_vs_t.<fmt> (with H(0,0,0) overlay)
          <out_dir>/D_waves_total_intensity_vs_t.<fmt>
          <out_dir>/S_wave_total_intensity_vs_t.<fmt>
      - All-wave overview comparison 1 x N_tbins plot:
          <out_dir>/amplitudes_overview_all_tbins.<fmt>
    """
    generated: List[Path] = []

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("[WARNING] matplotlib not found. Skipping multi-t-bin amplitude comparison plots.")
        return generated

    tbins_list = tbins if tbins is not None else config.dataset.get_t_bins()
    if not tbins_list or len(tbins_list) < 2:
        return generated

    effective_style = style or (getattr(config.plotting, "style", "default") if getattr(config, "plotting", None) else "default")
    apply_plot_style(effective_style)

    from .cross_section import parse_t_bin_from_string
    from .rundir import get_run_directory

    # Resolve output directory
    if output_dir:
        out_dir = Path(output_dir)
    else:
        base_dir = None
        if getattr(config, "multi_start", None) and getattr(config.multi_start, "output_dir", None):
            cand = Path(config.multi_start.output_dir)
            if cand.is_dir():
                base_dir = cand
        if base_dir is None and run_dir:
            cand = Path(run_dir)
            if cand.is_dir():
                base_dir = cand
        if base_dir is None:
            cfg_first = copy.deepcopy(config)
            cfg_first.dataset.t_bin = tbins_list[0]
            first_rdir = get_run_directory(cfg_first, run_dir)
            base_dir = first_rdir.root.parent

        out_dir = base_dir / "plots" / "amplitudes_vs_t"

    indiv_dir = out_dir / "individual"
    out_dir.mkdir(parents=True, exist_ok=True)
    indiv_dir.mkdir(parents=True, exist_ok=True)

    # Load amplitude data for each t-bin
    t_data: Dict[str, Dict[str, Any]] = {}
    for t_name in tbins_list:
        cfg_t = copy.deepcopy(config)
        cfg_t.dataset.t_bin = t_name
        t_rundir = get_run_directory(cfg_t, run_dir)
        amp_data = load_or_export_amplitude_data(cfg_t, t_rundir)
        if not amp_data or not amp_data.get("masses"):
            continue

        parsed_t = parse_t_bin_from_string(t_name)
        if parsed_t:
            t_min, t_max = parsed_t
            t_lbl = f"${_format_t_val(t_min)} < -t < {_format_t_val(t_max)}\\ \\mathrm{{GeV}}^2$"
        else:
            t_min, t_max = 0.0, 0.0
            t_lbl = t_name

        t_data[t_name] = {
            "name": t_name,
            "t_min": t_min,
            "t_max": t_max,
            "t_center": 0.5 * (t_min + t_max),
            "label": t_lbl,
            "masses": amp_data["masses"],
            "nominal": amp_data.get("nominal", {}),
            "nominal_params": amp_data.get("nominal_params", {}),
            "nominal_param_errors": amp_data.get("nominal_param_errors", {}),
            "mcmc": amp_data.get("mcmc", {}),
            "mcmc_weights": amp_data.get("mcmc_weights", {}),
            "mcmc_params": amp_data.get("mcmc_params", {}),
            "h0_data": amp_data.get("h0_data"),
            "rundir": t_rundir,
        }

    sorted_tbins = [t for t in tbins_list if t in t_data]
    sorted_tbins.sort(key=lambda t: (t_data[t]["t_min"] if t_data[t]["t_min"] is not None else 0.0))

    if not sorted_tbins:
        print("[WARNING] No amplitude data available across t-bins.")
        return generated

    active_waves_l = [w.l for w in getattr(config, "waves", [])] if getattr(config, "waves", None) else []

    def is_amp_active_any(k: str) -> bool:
        for t in sorted_tbins:
            y_vals = t_data[t]["nominal"].get(k, [])
            if any(abs(v) > 1e-9 for v in y_vals):
                return True
        l_val = extract_amplitude_l(k)
        if l_val is not None and l_val in active_waves_l:
            return True
        return False

    all_mag_keys = set()
    all_phi_keys = set()
    for t in sorted_tbins:
        for k in t_data[t]["nominal"].keys():
            if (k.startswith("a_") or k.startswith("b_")) and not (k.startswith("aphi_") or k.startswith("bphi_")):
                all_mag_keys.add(k)
            elif k.startswith("aphi_") or k.startswith("bphi_"):
                all_phi_keys.add(k)

    magnitudes = sorted([k for k in all_mag_keys if is_amp_active_any(k)], key=sort_amp_key)
    phases = sorted([k for k in all_phi_keys if is_amp_active_any(k)], key=sort_amp_key)

    n_cols = len(sorted_tbins)
    fig_width = max(3.8 * n_cols, 8.0)
    fig_height = 4.2

    # Mass limits
    fit_mass_min = None
    fit_mass_max = None
    fit_mass_step = 0.04
    if hasattr(config, "dataset") and getattr(config.dataset, "mass_bins", None):
        mb = config.dataset.mass_bins
        fit_mass_min = getattr(mb, "min", None)
        fit_mass_max = getattr(mb, "max", None)
        fit_mass_step = getattr(mb, "step", 0.04) or 0.04

    ref_masses = t_data[sorted_tbins[0]]["masses"]
    if fit_mass_min is None and ref_masses:
        fit_mass_min = ref_masses[0]
    if fit_mass_max is None and ref_masses:
        fit_mass_max = ref_masses[-1]
    if len(ref_masses) > 1:
        fit_mass_step = float(np.median(np.diff(ref_masses)))

    show_phase_ambiguity = getattr(config.plotting, "show_phase_ambiguity", True) if hasattr(config, "plotting") else True
    density_bands = getattr(config.plotting, "posterior_density_bands", True) if hasattr(config, "plotting") else True

    # --- 1. Individual Partial Wave Intensity vs -t ---
    for k in magnitudes:
        sqr_label = format_amplitude_squared_label(k)
        label = format_amplitude_label(k)
        is_mass_dep = is_amplitude_mass_dependent(k, config)

        fig, axes = plt.subplots(
            1, n_cols,
            figsize=(fig_width, fig_height),
            sharey=sharey,
            dpi=300,
        )
        axes_list = np.atleast_1d(axes)

        for idx, t_name in enumerate(sorted_tbins):
            ax = axes_list[idx]
            t_info = t_data[t_name]
            t_lbl = t_info.get("label", t_info.get("name", t_name))
            masses = t_info["masses"]
            y_nom = t_info["nominal"].get(k, [])
            y_sqr = [v**2 for v in y_nom]
            mcmc_dict = t_info["mcmc"].get(k)
            mcmc_w = t_info["mcmc_weights"].get(k)

            if is_mass_dep:
                if len(masses) >= 4:
                    m_s, y_s = smooth_curve(masses, y_sqr, num_points=200)
                    ax.plot(m_s, y_s, "-", linewidth=2.0, color="black", label=f"Nominal {sqr_label}", zorder=4)
                else:
                    ax.plot(masses, y_sqr, "-", linewidth=2.0, color="black", label=f"Nominal {sqr_label}", zorder=4)

                if mcmc_dict:
                    if density_bands:
                        plot_posterior_density_band(
                            ax=ax, masses=masses, mcmc_dict=mcmc_dict,
                            is_mass_dep=True, color="black", squared=True, n_levels=20,
                            base_alpha=0.22, max_alpha=0.62, fit_step=fit_mass_step, zorder=3,
                            mcmc_weights=mcmc_w,
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_dict, y_nom, (16.0, 84.0), squared=True, mcmc_weights=mcmc_w)
                        lo_95, hi_95 = extract_credible_bounds(masses, mcmc_dict, y_nom, (2.5, 97.5), squared=True, mcmc_weights=mcmc_w)
                        if len(masses) >= 4:
                            m_s, lo_95s = smooth_curve(masses, lo_95, num_points=200)
                            m_s, hi_95s = smooth_curve(masses, hi_95, num_points=200)
                            m_s, lo_68s = smooth_curve(masses, lo_68, num_points=200)
                            m_s, hi_68s = smooth_curve(masses, hi_68, num_points=200)
                            ax.fill_between(m_s, lo_95s, hi_95s, color="black", alpha=0.16, label="MCMC 95% CI", zorder=2)
                            ax.fill_between(m_s, lo_68s, hi_68s, color="black", alpha=0.32, label="MCMC 68% CI", zorder=3)
                        else:
                            ax.fill_between(masses, lo_95, hi_95, color="black", alpha=0.16, label="MCMC 95% CI", zorder=2)
                            ax.fill_between(masses, lo_68, hi_68, color="black", alpha=0.32, label="MCMC 68% CI", zorder=3)
            else:
                edges = compute_bin_edges(masses, default_step=fit_mass_step)
                if hasattr(ax, "stairs"):
                    ax.stairs(y_sqr, edges, linewidth=2.0, color="black", label=f"Nominal {sqr_label}", zorder=4)
                else:
                    ax.step(masses, y_sqr, where="mid", linewidth=2.0, color="black", label=f"Nominal {sqr_label}", zorder=4)
                if mcmc_dict:
                    if density_bands:
                        plot_posterior_density_band(
                            ax=ax, masses=masses, mcmc_dict=mcmc_dict,
                            is_mass_dep=False, color="black", squared=True, n_levels=20,
                            base_alpha=0.22, max_alpha=0.62, fit_step=fit_mass_step, zorder=3,
                            mcmc_weights=mcmc_w,
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_dict, y_nom, (16.0, 84.0), squared=True, mcmc_weights=mcmc_w)
                        lo_95, hi_95 = extract_credible_bounds(masses, mcmc_dict, y_nom, (2.5, 97.5), squared=True, mcmc_weights=mcmc_w)
                        x_step = np.repeat(edges, 2)[1:-1]
                        lo_95s = np.repeat(lo_95, 2)
                        hi_95s = np.repeat(hi_95, 2)
                        lo_68s = np.repeat(lo_68, 2)
                        hi_68s = np.repeat(hi_68, 2)
                        ax.fill_between(x_step, lo_95s, hi_95s, color="black", alpha=0.16, label="MCMC 95% CI", zorder=2)
                        ax.fill_between(x_step, lo_68s, hi_68s, color="black", alpha=0.32, label="MCMC 68% CI", zorder=3)

            ax.set_xlim(masses[0] - fit_mass_step * 0.6, masses[-1] + fit_mass_step * 0.6)
            ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=10.5)
            ax.set_title(t_lbl, fontsize=11, pad=6)
            if idx == 0:
                ax.set_ylabel(f"Intensity {sqr_label}", fontsize=11)
                ax.legend(frameon=True, fontsize=8, loc="best")
            if should_show_grid(config):
                ax.grid(True, linestyle=":", alpha=0.45)

        fig.suptitle(f"{config.project.name} - Intensity {sqr_label} vs $-t$", fontsize=12, y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        for fmt in formats:
            p = indiv_dir / f"{k}_intensity_vs_t.{fmt}"
            fig.savefig(p, bbox_inches="tight")
            generated.append(p)
        plt.close(fig)

    # Resolve reference phase info for each t-bin
    ref_info_by_tbin: Dict[str, Dict[str, Any]] = {}
    for t_name in sorted_tbins:
        t_info = t_data[t_name]
        ref_info_by_tbin[t_name] = resolve_reference_phase_info(
            config=config,
            nominal_amps=t_info.get("nominal", {}),
            active_phases=phases,
            masses=t_info.get("masses", []),
            nominal_params=t_info.get("nominal_params"),
            mcmc_params=t_info.get("mcmc_params"),
            mcmc_samples=t_info.get("mcmc"),
        )

    # --- 2. Individual Partial Wave Relative Phase vs -t ---
    for phase_k in phases:
        refl_char = "a" if phase_k.startswith("aphi_") else "b"
        first_t = sorted_tbins[0]
        first_ref_info = ref_info_by_tbin[first_t]

        # Reference amplitude phase is excluded from relative phase vs -t plots
        if phase_k in first_ref_info["omit_phases"]:
            for fmt in formats:
                stale_f = indiv_dir / f"{phase_k}_phase_vs_t.{fmt}"
                if stale_f.exists():
                    try:
                        stale_f.unlink()
                    except OSError:
                        pass
            continue

        ref_amp = first_ref_info.get("ref_names", {}).get(refl_char, first_ref_info["ref_name"])
        ref_lbl_str = format_reference_partial_wave_label(first_ref_info, refl_char=refl_char)
        phi_label = format_relative_phase_label(phase_k, ref_amp)
        phi_sym = format_relative_phase_symbol(phase_k, ref_amp)
        is_phi_mass_dep = is_amplitude_mass_dependent(phase_k, config)

        fig, axes = plt.subplots(
            1, n_cols,
            figsize=(fig_width, fig_height),
            sharey=True,
            dpi=300,
        )
        axes_list = np.atleast_1d(axes)

        for idx, t_name in enumerate(sorted_tbins):
            ax = axes_list[idx]
            t_info = t_data[t_name]
            t_ref_info = ref_info_by_tbin[t_name]
            t_lbl = t_info.get("label", t_info.get("name", t_name))
            masses = t_info["masses"]
            anchor_nom = t_ref_info["nominal_ref"].get(refl_char)
            anchor_mcmc = t_ref_info["mcmc_ref"].get(refl_char)
            y_rel, mcmc_rel = compute_relative_phase_data(
                masses, t_info["nominal"].get(phase_k, []), anchor_nom,
                t_info["mcmc"].get(phase_k), anchor_mcmc
            )
            mcmc_w = t_info["mcmc_weights"].get(phase_k)

            plot_phase_with_ambiguity(
                ax=ax,
                masses=masses,
                y_nom=y_rel,
                is_mass_dep=is_phi_mass_dep,
                color="black",
                label=f"Nominal +${phi_sym}$",
                label_neg=f"Sign Ambiguity -${phi_sym}$",
                mcmc_dict=mcmc_rel,
                fit_step=fit_mass_step,
                show_ambiguity=show_phase_ambiguity,
                alpha_band=0.32,
                base_alpha=0.20,
                max_alpha=0.60,
                density_band=density_bands,
                mcmc_weights=mcmc_w,
            )

            configure_phase_axis(ax)
            ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
            ax.set_xlim(masses[0] - fit_mass_step * 0.6, masses[-1] + fit_mass_step * 0.6)
            ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=10.5)
            ax.set_title(t_lbl, fontsize=11, pad=6)
            if idx == 0:
                ax.set_ylabel(f"Relative Phase {phi_label}\n(Reference Partial-Wave: {ref_lbl_str}) [rad]", fontsize=10.0)
                ax.legend(frameon=True, fontsize=7.5, loc="best", title=f"Ref Partial-Wave: {ref_lbl_str}", title_fontsize=7.5)
            if should_show_grid(config):
                ax.grid(True, linestyle=":", alpha=0.45)

        fig.suptitle(f"{config.project.name} - Relative Phase {phi_label} vs $-t$ (Reference Partial-Wave: {ref_lbl_str})", fontsize=11.5, y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        for fmt in formats:
            p = indiv_dir / f"{phase_k}_phase_vs_t.{fmt}"
            fig.savefig(p, bbox_inches="tight")
            generated.append(p)
        plt.close(fig)

    # --- 3. Summed Wave Totals vs -t ---
    def _compute_summed_mcmc(wave_subset: List[str], t_info: Dict[str, Any], masses: List[float]) -> Optional[Dict[float, List[float]]]:
        mcmc_dict = t_info.get("mcmc", {})
        avail_waves = [w for w in wave_subset if w in mcmc_dict and mcmc_dict[w]]
        if not avail_waves:
            return None
        res_mcmc: Dict[float, List[float]] = {}
        for m in masses:
            m_round = round(m, 4)
            wave_samps = []
            for w in avail_waves:
                s = mcmc_dict[w].get(m_round)
                if s is None:
                    for mk, mv in mcmc_dict[w].items():
                        if abs(mk - m_round) < 0.005:
                            s = mv
                            break
                if s:
                    wave_samps.append(s)
            if wave_samps:
                min_len = min(len(ws) for ws in wave_samps)
                if min_len > 0:
                    summed = [sum(ws[s_idx]**2 for ws in wave_samps) for s_idx in range(min_len)]
                    res_mcmc[m_round] = summed
        return res_mcmc if res_mcmc else None

    def _plot_summed_intensity(
        wave_subset: List[str],
        title_name: str,
        out_name: str,
        overlay_h0: bool = False,
    ):
        if not wave_subset:
            return
        fig, axes = plt.subplots(1, n_cols, figsize=(fig_width, fig_height), sharey=sharey, dpi=300)
        axes_list = np.atleast_1d(axes)

        is_mass_dep = all(is_amplitude_mass_dependent(w, config) for w in wave_subset)

        for idx, t_name in enumerate(sorted_tbins):
            ax = axes_list[idx]
            t_info = t_data[t_name]
            t_lbl = t_info.get("label", t_info.get("name", t_name))
            masses = t_info["masses"]

            # Overlay H0(0,0) if requested and available
            if overlay_h0 and t_info.get("h0_data"):
                h0_d = t_info["h0_data"]
                if h0_d.get("mass") and h0_d.get("val"):
                    m_h0 = h0_d["mass"]
                    y_h0 = h0_d["val"]
                    e_h0 = h0_d.get("err", [0.0] * len(y_h0))
                    ax.errorbar(
                        m_h0, y_h0, yerr=e_h0,
                        fmt="o", color="black", markersize=3.0, alpha=0.85,
                        capsize=1.2, elinewidth=0.9,
                        label=r"$H_0(0,0)$ (Total Yield)",
                        zorder=2,
                    )

            # Sum of nominal intensities
            y_nom_sum = [
                sum(t_info["nominal"][w][i]**2 for w in wave_subset if w in t_info["nominal"] and i < len(t_info["nominal"][w]))
                for i in range(len(masses))
            ]

            if is_mass_dep:
                if len(masses) >= 4:
                    m_s, y_s = smooth_curve(masses, y_nom_sum, num_points=200)
                    ax.plot(m_s, y_s, "-", linewidth=2.2, color="black", label=f"Nominal {title_name}", zorder=4)
                else:
                    ax.plot(masses, y_nom_sum, "-", linewidth=2.2, color="black", label=f"Nominal {title_name}", zorder=4)
            else:
                edges = compute_bin_edges(masses, default_step=fit_mass_step)
                if hasattr(ax, "stairs"):
                    ax.stairs(y_nom_sum, edges, linewidth=2.2, color="black", label=f"Nominal {title_name}", zorder=4)
                else:
                    ax.step(masses, y_nom_sum, where="mid", linewidth=2.2, color="black", label=f"Nominal {title_name}", zorder=4)

            # MCMC band
            summed_mcmc = _compute_summed_mcmc(wave_subset, t_info, masses)
            if summed_mcmc:
                if density_bands:
                    plot_posterior_density_band(
                        ax=ax, masses=masses, mcmc_dict=summed_mcmc,
                        is_mass_dep=is_mass_dep, color="black", squared=False, n_levels=20,
                        base_alpha=0.22, max_alpha=0.62, fit_step=fit_mass_step, zorder=3,
                    )
                else:
                    lo_68, hi_68 = extract_credible_bounds(masses, summed_mcmc, y_nom_sum, (16.0, 84.0), squared=False)
                    lo_95, hi_95 = extract_credible_bounds(masses, summed_mcmc, y_nom_sum, (2.5, 97.5), squared=False)
                    if is_mass_dep and len(masses) >= 4:
                        m_s, lo_95s = smooth_curve(masses, lo_95, num_points=200)
                        m_s, hi_95s = smooth_curve(masses, hi_95, num_points=200)
                        m_s, lo_68s = smooth_curve(masses, lo_68, num_points=200)
                        m_s, hi_68s = smooth_curve(masses, hi_68, num_points=200)
                        ax.fill_between(m_s, lo_95s, hi_95s, color="black", alpha=0.16, label="MCMC 95% CI", zorder=2)
                        ax.fill_between(m_s, lo_68s, hi_68s, color="black", alpha=0.32, label="MCMC 68% CI", zorder=3)
                    else:
                        edges = compute_bin_edges(masses, default_step=fit_mass_step)
                        x_step = np.repeat(edges, 2)[1:-1]
                        lo_95s = np.repeat(lo_95, 2)
                        hi_95s = np.repeat(hi_95, 2)
                        lo_68s = np.repeat(lo_68, 2)
                        hi_68s = np.repeat(hi_68, 2)
                        ax.fill_between(x_step, lo_95s, hi_95s, color="black", alpha=0.16, label="MCMC 95% CI", zorder=2)
                        ax.fill_between(x_step, lo_68s, hi_68s, color="black", alpha=0.32, label="MCMC 68% CI", zorder=3)

            ax.set_xlim(masses[0] - fit_mass_step * 0.6, masses[-1] + fit_mass_step * 0.6)
            ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=10.5)
            ax.set_title(t_lbl, fontsize=11, pad=6)
            if idx == 0:
                ax.set_ylabel(f"Intensity {title_name}", fontsize=11)
                ax.legend(frameon=True, fontsize=8, loc="best")
            if should_show_grid(config):
                ax.grid(True, linestyle=":", alpha=0.45)

        fig.suptitle(f"{config.project.name} - {title_name} vs $-t$", fontsize=12, y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        for fmt in formats:
            p_out = out_dir / f"{out_name}.{fmt}"
            fig.savefig(p_out, bbox_inches="tight")
            generated.append(p_out)
        plt.close(fig)

    # 3a. Total Intensity
    _plot_summed_intensity(magnitudes, r"Total $\sum |V|^2$", "total_intensity_vs_t", overlay_h0=True)

    # 3b. Summed intensities for each active L (S, P, D, etc.)
    for l_val in sorted(set(extract_amplitude_l(k) for k in magnitudes if extract_amplitude_l(k) is not None)):
        l_sym = L_SPECTROSCOPIC.get(l_val, f"L{l_val}")
        l_waves = [k for k in magnitudes if extract_amplitude_l(k) == l_val]
        if l_waves:
            if l_val == 0:
                _plot_summed_intensity(l_waves, rf"Total ${l_sym}$-Wave $|{l_sym}|^2$", f"{l_sym}_wave_total_intensity_vs_t")
            else:
                _plot_summed_intensity(l_waves, rf"Total ${l_sym}$-Wave $\sum |{l_sym}|^2$", f"{l_sym}_waves_total_intensity_vs_t")

    # --- 4. Overview Comparison Plot Across All t-bins ---
    fig, axes = plt.subplots(1, n_cols, figsize=(fig_width, fig_height), sharey=sharey, dpi=300)
    axes_list = np.atleast_1d(axes)

    all_wave_keys = []
    for k in magnitudes:
        p = get_amplitude_symbol_parts(k)
        ident = (p[0], p[1], p[2]) if p else k
        if ident not in all_wave_keys:
            all_wave_keys.append(ident)

    wave_linestyle_map = {
        ident: PHYSICS_LINESTYLES[i % len(PHYSICS_LINESTYLES)]
        for i, ident in enumerate(all_wave_keys)
    }

    for idx, t_name in enumerate(sorted_tbins):
        ax = axes_list[idx]
        t_info = t_data[t_name]
        t_lbl = t_info.get("label", t_info.get("name", t_name))
        masses = t_info["masses"]

        # 1. Overlay H0(0,0)
        if t_info.get("h0_data"):
            h0_d = t_info["h0_data"]
            if h0_d.get("mass") and h0_d.get("val"):
                m_h0 = h0_d["mass"]
                y_h0 = h0_d["val"]
                e_h0 = h0_d.get("err", [0.0] * len(y_h0))
                ax.errorbar(
                    m_h0, y_h0, yerr=e_h0,
                    fmt="o", color="black", markersize=2.6, alpha=0.8,
                    capsize=1.0, elinewidth=0.8,
                    label=r"$H_0(0,0)$ (Total)",
                    zorder=2,
                )

        # 2. Total intensity
        tot_nom = [
            sum(t_info["nominal"][k][i]**2 for k in magnitudes if k in t_info["nominal"] and i < len(t_info["nominal"][k]))
            for i in range(len(masses))
        ]
        if len(masses) >= 4:
            m_s, y_s = smooth_curve(masses, tot_nom, num_points=200)
            ax.plot(m_s, y_s, "-", linewidth=2.2, color="black", alpha=0.95, label=r"Total $\sum |V|^2$", zorder=5)
        else:
            ax.plot(masses, tot_nom, "-", linewidth=2.2, color="black", alpha=0.95, label=r"Total $\sum |V|^2$", zorder=5)

        # 3. Individual waves
        for k in magnitudes:
            if k not in t_info["nominal"]:
                continue
            y_w = [v**2 for v in t_info["nominal"][k]]
            p = get_amplitude_symbol_parts(k)
            ident = (p[0], p[1], p[2]) if p else k
            ls = wave_linestyle_map.get(ident, "--")
            sqr_label = format_amplitude_squared_label(k)
            is_mass_dep = is_amplitude_mass_dependent(k, config)
            if is_mass_dep and len(masses) >= 4:
                m_w_s, y_w_s = smooth_curve(masses, y_w, num_points=200)
                ax.plot(m_w_s, y_w_s, linestyle=ls, linewidth=1.6, color="black", label=sqr_label, zorder=4)
            else:
                edges = compute_bin_edges(masses, default_step=fit_mass_step)
                if hasattr(ax, "stairs"):
                    ax.stairs(y_w, edges, linestyle=ls, linewidth=1.6, color="black", label=sqr_label, zorder=4)
                else:
                    ax.step(masses, y_w, where="mid", linestyle=ls, linewidth=1.6, color="black", label=sqr_label, zorder=4)

        ax.set_xlim(masses[0] - fit_mass_step * 0.6, masses[-1] + fit_mass_step * 0.6)
        ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=10.5)
        ax.set_title(t_lbl, fontsize=11, pad=6)
        if idx == 0:
            ax.set_ylabel("Intensity $|V|^2$", fontsize=11)
        if should_show_grid(config):
            ax.grid(True, linestyle=":", alpha=0.45)

    handles_ov, labels_ov = axes_list[0].get_legend_handles_labels()
    if handles_ov:
        if len(handles_ov) > 6:
            fig.legend(
                handles_ov, labels_ov,
                loc="upper center", bbox_to_anchor=(0.5, 0.02),
                ncol=min(7, len(handles_ov)), frameon=True, fontsize=7.8,
            )
            layout_rect = [0, 0.06, 1, 0.94]
        else:
            axes_list[0].legend(frameon=True, fontsize=7.5, loc="best")
            layout_rect = [0, 0, 1, 0.95]
    else:
        layout_rect = [0, 0, 1, 0.95]

    fig.suptitle(f"{config.project.name} - Partial Wave Intensities Overview vs $-t$", fontsize=12, y=0.98)
    fig.tight_layout(rect=layout_rect)
    for fmt in formats:
        p_overview = out_dir / f"amplitudes_overview_all_tbins.{fmt}"
        fig.savefig(p_overview, bbox_inches="tight")
        generated.append(p_overview)

        # Copy to each t-bin's plots directory
        for t_name in sorted_tbins:
            t_pdir = t_data[t_name]["rundir"].plots_dir
            if t_pdir.is_dir():
                dst_p = t_pdir / f"amplitudes_overview_all_tbins.{fmt}"
                shutil.copyfile(p_overview, dst_p)
                generated.append(dst_p)
    plt.close(fig)

    # --- 5. Total Magnitude Squared for S, P, D Waves in Total Reflectivities vs -t ---
    try:
        spd_plots = plot_spd_total_intensity_vs_t(
            config=config,
            t_data=t_data,
            sorted_tbins=sorted_tbins,
            magnitudes=magnitudes,
            out_dir=out_dir,
            fit_mass_step=fit_mass_step,
            formats=formats,
            style=None,
        )
        generated.extend(spd_plots)
    except Exception as e:
        print(f"[WARNING] Could not generate S/P/D total reflectivity vs -t plots: {e}")

    # --- 6. Fitted Resonance Parameter Values vs -t ---
    try:
        res_param_plots = plot_resonance_parameters_vs_t(
            config=config,
            t_data=t_data,
            sorted_tbins=sorted_tbins,
            output_dir=out_dir,
            formats=formats,
            style=None,
        )
        generated.extend(res_param_plots)
    except Exception as e:
        print(f"[WARNING] Could not generate resonance parameter vs -t plots: {e}")

    print(f"[OK] Generated {len(generated)} multi-t-bin amplitude comparison plot files in {out_dir}")
    return generated


def run_all_diagnostics_and_plots(
    config: FitConfig,
    rundir: RunDirectory,
    formats: Tuple[str, ...] = ("png", "pdf"),
    only_propagated: bool = False,
    only_distributions: bool = False,
    only_amplitudes: bool = False,
    style: Optional[str] = None,
) -> Dict[str, Any]:
    """Run full diagnostic report generation and publication plotting."""
    from .summary import FitSummary

    effective_style = style or (getattr(config.plotting, "style", "default") if getattr(config, "plotting", None) else "default")

    # Fast path: only generate propagated amplitude/distribution plots
    if only_propagated or only_distributions or only_amplitudes:
        plot_files = plot_amplitudes_and_bands(
            rundir.nominal_best_root,
            rundir.mcmc_propagated_root,
            rundir.plots_dir,
            config,
            formats=formats,
            only_distributions=only_distributions,
            only_amplitudes=only_amplitudes,
            style=effective_style,
            rundir=rundir,
            unnormalized_moments_file=getattr(config.plotting, "unnormalized_moments_file", None),
        )
        return {
            "summary": None,
            "reports": [],
            "plots": plot_files,
        }

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
        only_distributions=only_distributions,
        only_amplitudes=only_amplitudes,
        style=effective_style,
        rundir=rundir,
        unnormalized_moments_file=getattr(config.plotting, "unnormalized_moments_file", None),
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
                                is_valid = tokens["VALID"] == "1"
                                status = int(tokens.get("STATUS", -1))
                                if is_valid and status == 0:
                                    ms_results.append({
                                        "seed": int(tokens.get("SEED", 0)),
                                        "chi2": float(tokens["CHI2"]),
                                        "is_valid": True,
                                        "status": 0,
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
                style=effective_style,
                show_grid=should_show_grid(config),
                bins=getattr(config.plotting, "chi2_landscape_bins", None) if getattr(config, "plotting", None) else None,
            )
            plot_files.extend(ms_plots)

    # Generate moment comparison plots if enabled
    if getattr(config, "plotting", None) and config.plotting.enabled:
        mom_plots = plot_moments_comparison(
            config,
            rundir,
            formats=formats,
            style=effective_style,
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
                style=effective_style,
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
    only_propagated: bool = False,
    only_distributions: bool = False,
    only_amplitudes: bool = False,
    style: Optional[str] = None,
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
            elif rundir.mcmc_propagated_root.is_file():
                rundir.save_nominal_best(rundir.mcmc_propagated_root)
            else:
                return (t_name, False, 0, f"No nominal best fit found at {rundir.nominal_best_root}. Skipping.")

        res = run_all_diagnostics_and_plots(
            cfg_t,
            rundir,
            formats=formats,
            only_propagated=only_propagated,
            only_distributions=only_distributions,
            only_amplitudes=only_amplitudes,
            style=style,
        )
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
    only_propagated: bool = False,
    only_distributions: bool = False,
    only_amplitudes: bool = False,
    style: Optional[str] = None,
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

    effective_style = style or (getattr(config.plotting, "style", "default") if getattr(config, "plotting", None) else "default")

    mode_label = "MCMC Posterior Distributions" if only_distributions else (
        "Propagated Amplitude Plots" if only_propagated else "Publication Plots"
    )
    print(f"\n=======================================================")
    print(f" Moments2Amplitudes {mode_label} Across t-Bins")
    print(f" Project:      {config.project.name}")
    print(f" t-Bins ({len(tbins_list)}): {', '.join(tbins_list)}")
    if effective_style and effective_style.lower() not in ("default", "none", "matplotlib", "standard"):
        print(f" Style:        {effective_style} (mplhep)")
    else:
        print(f" Style:        default (standard matplotlib)")
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
                only_propagated,
                only_distributions,
                only_amplitudes,
                effective_style,
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

    # Generate multi-t-bin moment evolution plots across all t-bins
    plot_mom_t = getattr(config.plotting, "moments_vs_t", True) if getattr(config, "plotting", None) else True
    if plot_mom_t and len(tbins_list) > 1 and not (only_distributions or only_amplitudes):
        print(f"\n[INFO] Generating cross-t-bin moment evolution plots across {len(tbins_list)} t-bins...")
        try:
            t_mom_plots = plot_moments_across_tbins(
                config,
                tbins=tbins_list,
                run_dir=run_dir,
                formats=formats,
                style=effective_style,
                sharey=getattr(config.plotting, "moments_vs_t_sharey", False) if getattr(config, "plotting", None) else False,
            )
            total_plots += len(t_mom_plots)
        except Exception as e:
            print(f"[WARNING] Could not generate multi-t-bin moment plots: {e}")

    # Generate multi-t-bin amplitude evolution plots across all t-bins
    plot_amp_t = getattr(config.plotting, "amplitudes_vs_t", True) if getattr(config, "plotting", None) else True
    if plot_amp_t and len(tbins_list) > 1 and not only_distributions:
        print(f"\n[INFO] Generating cross-t-bin amplitude evolution plots across {len(tbins_list)} t-bins...")
        try:
            t_amp_plots = plot_amplitudes_across_tbins(
                config,
                tbins=tbins_list,
                run_dir=run_dir,
                formats=formats,
                style=effective_style,
                sharey=getattr(config.plotting, "amplitudes_vs_t_sharey", False) if getattr(config, "plotting", None) else False,
            )
            total_plots += len(t_amp_plots)
        except Exception as e:
            print(f"[WARNING] Could not generate multi-t-bin amplitude plots: {e}")

    return {
        "success": all_success,
        "total_plots": total_plots,
        "tbins": results_by_tbin,
    }


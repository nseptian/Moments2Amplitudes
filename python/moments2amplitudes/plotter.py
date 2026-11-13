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
        return

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
                except Exception:
                    pass
            if is_nominal_converged:
                valid_fits.append((0.0, rundir.nominal_best_root))

    return valid_fits[:n_lowest]


def plot_multistart_landscape(
    results: List[Dict[str, Any]],
    output_dir: Path | str,
    project_name: str = "Moments2Amplitudes",
    plateau_delta: float = 10.0,
    formats: Tuple[str, ...] = ("png", "pdf"),
    style: Optional[str] = "default",
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
    counts, bins, patches = ax.hist(
        chi2s, bins=40, color="#B0C4DE", edgecolor="#1A2E40", alpha=0.75,
        label=f"Converged Fits (N={len(chi2s)})"
    )

    ax.axvline(min_chi2, color="#D55E00", linestyle="--", linewidth=1.8, label=f"Global Min $\\chi^2 = {min_chi2:.2f}$")
    ax.axvspan(min_chi2, plateau_thresh, color="#009E73", alpha=0.18, label=f"Plateau $\\leq +{plateau_delta:.1f}$")

    ax.set_xlabel(r"Total $\chi^2$", fontsize=11.5)
    ax.set_ylabel("Number of Starts", fontsize=11.5)
    ax.set_title(f"{project_name} - Multi-Start $\\chi^2$ Landscape (Converged Fits)", fontsize=12.5, pad=8)
    ax.legend(frameon=True, fontsize=9.5)
    ax.grid(True, linestyle=":", alpha=0.45)

    for fmt in formats:
        p = out_dir / f"multistart_chi2_landscape.{fmt}"
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
    Parse amplitude branch name (e.g. 'a_2_-1', 'b_0_0', 'aphi_1_1') into:
    (l_symbol, refl, m_str, wave_kind)
    where:
      l_symbol: 'S', 'P', 'D', 'F', ...
      refl: '+' for 'a'/'aphi', '-' for 'b'/'bphi'
      m_str: '-1', '0', '1', etc.
      wave_kind: 'amp' or 'phase'
    """
    parts = name.split("_")
    if len(parts) >= 3:
        prefix = parts[0]
        l_str = parts[1]
        m_str = "_".join(parts[2:])
        if prefix in ("a", "b", "aphi", "bphi"):
            try:
                l_val = int(l_str)
                l_sym = L_SPECTROSCOPIC.get(l_val, f"L_{{{l_val}}}")
            except ValueError:
                l_sym = l_str
            refl = "+" if prefix in ("a", "aphi") else "-"
            kind = "phase" if prefix in ("aphi", "bphi") else "amp"
            return (l_sym, refl, m_str, kind)
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
    percentiles: Tuple[float, float] = (16.0, 84.0),
    squared: bool = False,
) -> Tuple[List[float], List[float]]:
    """Extract lower and upper percentiles for credible interval band across mass bins."""
    import numpy as np
    lo, hi = [], []
    p_lo, p_hi = percentiles
    for idx, m in enumerate(masses):
        m_round = round(m, 4)
        samps = mcmc_dict.get(m_round)
        if samps is None:
            for mk, mv in mcmc_dict.items():
                if abs(mk - m_round) < 0.005:
                    samps = mv
                    break
        if samps:
            vals = [v**2 for v in samps] if squared else samps
            lo.append(float(np.percentile(vals, p_lo)))
            hi.append(float(np.percentile(vals, p_hi)))
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
) -> None:
    """
    Draw a continuous credible uncertainty band directly from the empirical
    propagated posterior distribution, where the shading intensity/density
    reflects the local probability density p(y | mass).

    A visible base envelope is drawn over the credible region (e.g. 95% CI)
    followed by stacked quantile shells centered on the median so that the band
    is clearly visible on paper and optical density gracefully climbs
    towards the posterior mode/median.
    """
    import numpy as np

    if not mcmc_dict:
        return

    has_valid_samples = any(len(samps) > 0 for samps in mcmc_dict.values())
    if not has_valid_samples:
        return

    # Calibrate alpha parameters: ensure a visibly prominent base envelope
    if base_alpha is None:
        eff_base_alpha = max(0.18, max_alpha * 0.35)
    else:
        eff_base_alpha = max(0.10, min(base_alpha, max_alpha * 0.70))

    coverages = np.linspace(0.95, 0.05, n_levels)
    if n_levels > 1 and max_alpha > eff_base_alpha:
        rem = (1.0 - max_alpha) / (1.0 - eff_base_alpha)
        d_alpha = max(0.01, 1.0 - rem**(1.0 / (n_levels - 1)))
    else:
        d_alpha = max(0.01, max_alpha / float(max(1, n_levels)))

    for idx, cov in enumerate(coverages):
        p_lo = (1.0 - cov) / 2.0 * 100.0
        p_hi = (1.0 + cov) / 2.0 * 100.0
        lo_list: List[float] = []
        hi_list: List[float] = []

        for m in masses:
            m_round = round(m, 4)
            samps = mcmc_dict.get(m_round)
            if samps is None:
                for mk, mv in mcmc_dict.items():
                    if abs(mk - m_round) < 0.005:
                        samps = mv
                        break
            if samps:
                arr = np.array(samps, dtype=float)
                if squared:
                    vals = arr**2
                elif mirrored:
                    vals = -arr
                else:
                    vals = arr
                lo_list.append(float(np.percentile(vals, p_lo)))
                hi_list.append(float(np.percentile(vals, p_hi)))
            else:
                lo_list.append(0.0)
                hi_list.append(0.0)

        layer_alpha = eff_base_alpha if idx == 0 else d_alpha

        if is_mass_dep:
            if len(masses) >= 4:
                m_s, lo_s = smooth_curve(masses, lo_list, num_points=200)
                m_s, hi_s = smooth_curve(masses, hi_list, num_points=200)
                ax.fill_between(m_s, lo_s, hi_s, color=color, alpha=layer_alpha, linewidth=0, zorder=zorder)
            else:
                ax.fill_between(masses, lo_list, hi_list, color=color, alpha=layer_alpha, linewidth=0, zorder=zorder)
        else:
            edges = compute_bin_edges(masses, default_step=fit_step)
            x_step = np.repeat(edges, 2)[1:-1]
            lo_step = np.repeat(lo_list, 2)
            hi_step = np.repeat(hi_list, 2)
            ax.fill_between(x_step, lo_step, hi_step, color=color, alpha=layer_alpha, linewidth=0, zorder=zorder)


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
            m_s, y_s = smooth_curve(masses, y_nom, num_points=200)
        else:
            m_s, y_s = np.array(masses, dtype=float), y_nom_arr

        # Primary branch
        ax.plot(m_s, y_s, linestyle=linestyle_pos, linewidth=linewidth, color=color, label=label, zorder=zorder)

        # Primary MCMC credible band
        if mcmc_dict is not None:
            if density_band:
                plot_posterior_density_band(
                    ax=ax, masses=masses, mcmc_dict=mcmc_dict, is_mass_dep=True,
                    color=color, squared=False, n_levels=20, base_alpha=eff_base_alpha,
                    max_alpha=eff_max_alpha, fit_step=fit_step, zorder=zorder - 1, mirrored=False
                )
            else:
                lo_68, hi_68 = extract_credible_bounds(masses, mcmc_dict, y_nom, (16.0, 84.0))
                if len(masses) >= 4:
                    m_s68, lo_s = smooth_curve(masses, lo_68, num_points=200)
                    m_s68, hi_s = smooth_curve(masses, hi_68, num_points=200)
                    ax.fill_between(m_s68, lo_s, hi_s, color=color, alpha=alpha_band, zorder=zorder - 1)
                else:
                    ax.fill_between(masses, lo_68, hi_68, color=color, alpha=alpha_band, zorder=zorder - 1)

        # Sign-ambiguous branch (phi -> -phi)
        if show_ambiguity and not is_trivial_zero:
            ax.plot(m_s, -y_s, linestyle=linestyle_neg, linewidth=linewidth, color=color, label=label_neg, zorder=zorder)
            if mcmc_dict is not None:
                if density_band:
                    plot_posterior_density_band(
                        ax=ax, masses=masses, mcmc_dict=mcmc_dict, is_mass_dep=True,
                        color=color, squared=False, n_levels=20, base_alpha=eff_base_alpha,
                        max_alpha=eff_max_alpha, fit_step=fit_step, zorder=zorder - 1, mirrored=True
                    )
                else:
                    if len(masses) >= 4:
                        ax.fill_between(m_s68, -hi_s, -lo_s, color=color, alpha=alpha_band, zorder=zorder - 1)
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
                    max_alpha=eff_max_alpha, fit_step=fit_step, zorder=zorder - 1, mirrored=False
                )
            else:
                lo_68, hi_68 = extract_credible_bounds(masses, mcmc_dict, y_nom, (16.0, 84.0))
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
                        max_alpha=eff_max_alpha, fit_step=fit_step, zorder=zorder - 1, mirrored=True
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

    # Resolve unnormalized moments file for H_0_0_0 reference
    unnorm_path = resolve_unnormalized_moments_file(config, rundir, unnormalized_moments_file)
    h0_data: Optional[Dict[str, List[float]]] = None

    # Extract nominal result tree
    masses: List[float] = []
    nominal_amps: Dict[str, List[float]] = {}
    mcmc_samples: Dict[str, Dict[float, List[float]]] = {}
    mcmc_params: Dict[str, List[float]] = {}

    # Read with uproot or PyROOT
    loaded = False
    try:
        import uproot
        with uproot.open(nominal_root) as f:
            tname = "result" if "result" in f else ("central_result" if "central_result" in f else None)
            if tname:
                tree = f[tname]
                masses = list(tree["mass_bin"].array())
                for k in tree.keys():
                    if k != "mass_bin" and not k.startswith("H_"):
                        nominal_amps[k] = list(tree[k].array())
                loaded = True

        if mcmc_root and Path(mcmc_root).is_file():
            with uproot.open(mcmc_root) as mf:
                # Identify converged MCMC entries (Status == 0 and isValid == True)
                conv_entries = set()
                has_conv_tree = False
                for ctree_name in ("mcmc_params", "mcmc_chi2"):
                    if ctree_name in mf:
                        has_conv_tree = True
                        ctree = mf[ctree_name]
                        ckeys = ctree.keys()
                        statuses = list(ctree["Status"].array()) if "Status" in ckeys else None
                        valids = list(ctree["isValid"].array()) if "isValid" in ckeys else None
                        m_entries = list(ctree["mcmc_entry"].array()) if "mcmc_entry" in ckeys else None
                        n_tot = len(statuses) if statuses is not None else (len(valids) if valids is not None else 0)
                        for idx in range(n_tot):
                            st_ok = (statuses[idx] == 0) if statuses is not None else True
                            v_ok = bool(valids[idx]) if valids is not None else True
                            if st_ok and v_ok:
                                entry_id = m_entries[idx] if m_entries is not None else idx
                                conv_entries.add(entry_id)
                        break

                mtname = "mcmc_result" if "mcmc_result" in mf else ("result" if "result" in mf else None)
                if mtname:
                    mtree = mf[mtname]
                    mbins = list(mtree["mass_bin"].array())
                    m_result_entries = list(mtree["mcmc_entry"].array()) if "mcmc_entry" in mtree else None
                    n_mass_bins = len(masses) if masses else 1
                    for k in nominal_amps.keys():
                        if k in mtree:
                            arr = list(mtree[k].array())
                            mcmc_samples[k] = {}
                            for i, (m, v) in enumerate(zip(mbins, arr)):
                                if has_conv_tree:
                                    entry_id = m_result_entries[i] if m_result_entries is not None else (i // n_mass_bins)
                                    if entry_id not in conv_entries:
                                        continue
                                mcmc_samples[k].setdefault(round(float(m), 4), []).append(float(v))
                if "mcmc_params" in mf:
                    ptree = mf["mcmc_params"]
                    p_entries = list(ptree["mcmc_entry"].array()) if "mcmc_entry" in ptree else None
                    for pk in ptree.keys():
                        if pk not in ("mcmc_entry", "Status", "isValid", "isValidError"):
                            vals = list(ptree[pk].array())
                            mcmc_params[pk] = [
                                float(v) for i, v in enumerate(vals)
                                if (not has_conv_tree or (p_entries[i] if p_entries is not None else i) in conv_entries)
                            ]
    except Exception:
        pass

    if not loaded:
        # Fallback: export to JSON using ROOT macro ExportAmplitudes.C
        try:
            repo_root = Path(os.environ.get("MOMS2AMPS", str(Path(__file__).resolve().parent.parent.parent)))
            export_macro = repo_root / "macros" / "ExportAmplitudes.C"
            cache_dir = out_dir / ".cache"
            cache_dir.mkdir(parents=True, exist_ok=True)
            json_path = cache_dir / "amplitudes_data.json"

            # Check if cached amplitudes JSON already exists and is fresh
            if json_path.is_file() and json_path.stat().st_size > 10:
                is_stale = False
                json_mtime = json_path.stat().st_mtime
                if Path(nominal_root).is_file() and json_mtime < Path(nominal_root).stat().st_mtime:
                    is_stale = True
                if mcmc_root and Path(mcmc_root).is_file():
                    if json_mtime < Path(mcmc_root).stat().st_mtime:
                        is_stale = True
                if unnorm_path and Path(unnorm_path).is_file():
                    if json_mtime < Path(unnorm_path).stat().st_mtime:
                        is_stale = True

                if not is_stale:
                    try:
                        with open(json_path, "r") as jf:
                            jdata = json.load(jf)
                            raw_masses = jdata.get("masses", [])
                            raw_nominal = jdata.get("nominal", {})
                            raw_mcmc = jdata.get("mcmc", {})
                            raw_params = jdata.get("mcmc_params", {})
                            # If MCMC root file was provided and non-empty, but cached JSON has no MCMC samples, cache is stale!
                            if mcmc_root and Path(mcmc_root).is_file() and (not raw_mcmc or not raw_params):
                                is_stale = True
                            elif raw_masses and raw_nominal:
                                masses = [float(m) for m in raw_masses]
                                nominal_amps = {k: [float(v) for v in vals] for k, vals in raw_nominal.items()}
                                for k, m_dict in raw_mcmc.items():
                                    mcmc_samples[k] = {}
                                    for m_str, vals in m_dict.items():
                                        try:
                                            mcmc_samples[k][round(float(m_str), 4)] = [float(v) for v in vals]
                                        except ValueError:
                                            pass
                                for pk, pvals in raw_params.items():
                                    mcmc_params[pk] = [float(v) for v in pvals]
                                if "h0_data" in jdata and jdata["h0_data"].get("mass") and jdata["h0_data"].get("val"):
                                    h0_data = jdata["h0_data"]
                                loaded = True
                    except Exception:
                        loaded = False

            if not loaded and export_macro.is_file() and Path(nominal_root).is_file():
                root_bin = get_root_binary(config)
                mcmc_arg = str(Path(mcmc_root).resolve()) if (mcmc_root and Path(mcmc_root).is_file()) else ""
                unnorm_arg = str(unnorm_path.resolve()) if (unnorm_path and unnorm_path.is_file()) else ""
                cmd = [
                    root_bin,
                    "-l",
                    "-b",
                    "-q",
                    f'{export_macro.resolve()}("{Path(nominal_root).resolve()}", "{mcmc_arg}", "{json_path.resolve()}", "{unnorm_arg}")'
                ]
                env = get_root_subprocess_env(config)
                run_res = subprocess.run(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
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
                            raw_params = jdata.get("mcmc_params", {})
                            for pk, pvals in raw_params.items():
                                mcmc_params[pk] = [float(v) for v in pvals]
                            if "h0_data" in jdata and jdata["h0_data"].get("mass") and jdata["h0_data"].get("val"):
                                h0_data = jdata["h0_data"]
                            loaded = True
                else:
                    if run_res.returncode != 0:
                        print(f"[WARNING] ExportAmplitudes macro execution failed:\n{run_res.stderr}\n{run_res.stdout}")
        except Exception as e:
            print(f"[WARNING] Error running ExportAmplitudes fallback: {e}")

    if not loaded or not masses:
        print(f"[WARNING] Could not read amplitude branches from {nominal_root}")
        return generated

    if h0_data is None:
        cache_comp_json = (rundir.root / ".cache" / "moment_comparison_data.json") if rundir else (out_dir.parent / ".cache" / "moment_comparison_data.json")
        h0_data = load_h0_moment_data(unnorm_path, cache_comp_json)

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
                            base_alpha=0.18, max_alpha=0.55, fit_step=fit_step, zorder=3
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True)
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
                            base_alpha=0.18, max_alpha=0.55, fit_step=fit_step, zorder=3
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True)
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
        ax.grid(True, linestyle=":", alpha=0.45)

        # 1d. Bottom panel: Phase as function of mass
        if ax_phase is not None:
            for idx, k in enumerate(sorted(phases)):
                y_nom = nominal_amps[k]
                label = format_amplitude_label(k, plus_minus=True)
                is_mass_dep = is_amplitude_mass_dependent(k, config)
                c = "black"
                ls = get_linestyle(k)
                neg_ls = ":" if ls == "--" else "--"

                plot_phase_with_ambiguity(
                    ax=ax_phase,
                    masses=masses,
                    y_nom=y_nom,
                    is_mass_dep=is_mass_dep,
                    color=c,
                    label=label,
                    mcmc_dict=mcmc_samples.get(k),
                    fit_step=fit_step,
                    show_ambiguity=show_phase_ambiguity,
                    alpha_band=0.30,
                    base_alpha=0.18,
                    max_alpha=0.55,
                    linestyle_pos=ls,
                    linestyle_neg=neg_ls,
                    density_band=density_bands,
                )

            configure_phase_axis(ax_phase)
            ax_phase.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
            ax_phase.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11.5)
            ax_phase.set_ylabel("Phase [rad]", fontsize=11.5)
            ax_phase.legend(frameon=True, fontsize=8, ncol=min(3, max(1, len(phases))), loc="best")
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

        for idx, k in enumerate(sorted(phases)):
            y_nom = nominal_amps[k]
            label = format_amplitude_label(k, plus_minus=True)
            is_mass_dep = is_amplitude_mass_dependent(k, config)
            c = "black"
            ls = get_linestyle(k)
            neg_ls = ":" if ls == "--" else "--"

            plot_phase_with_ambiguity(
                ax=ax,
                masses=masses,
                y_nom=y_nom,
                is_mass_dep=is_mass_dep,
                color=c,
                label=label,
                mcmc_dict=mcmc_samples.get(k),
                fit_step=fit_step,
                show_ambiguity=show_phase_ambiguity,
                alpha_band=0.30,
                base_alpha=0.18,
                max_alpha=0.55,
                linestyle_pos=ls,
                linestyle_neg=neg_ls,
                density_band=density_bands,
            )

        configure_phase_axis(ax)
        ax.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
        ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11.5)
        ax.set_ylabel("Phase [rad]", fontsize=11.5)
        ax.set_title(f"{config.project.name} - Fitted Wave Phases ($\\pm$ Sign Ambiguity)", fontsize=12.5, pad=8)
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.legend(frameon=True, fontsize=9.5, ncol=2)
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
                            base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0))
                        lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5))
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
                            base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0))
                        lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5))
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
                            base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True)
                        lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5), squared=True)
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
                            base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3
                        )
                    else:
                        lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True)
                        lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5), squared=True)
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
                ax.legend(handles=handles, labels=labels, frameon=True, fontsize=9, loc="best")
            else:
                ax.legend(frameon=True, fontsize=9, loc="best")
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

            if phase_k:
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
                                base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3
                            )
                        else:
                            lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True)
                            lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5), squared=True)
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
                                base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3
                            )
                        else:
                            lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (16.0, 84.0), squared=True)
                            lo_95, hi_95 = extract_credible_bounds(masses, mcmc_samples[k], y_nom, (2.5, 97.5), squared=True)
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
                ax_top.grid(True, linestyle=":", alpha=0.45)
                ax_top.tick_params(labelbottom=False)

                y_nom_phi = nominal_amps[phase_k]
                phi_label = format_amplitude_label(phase_k)
                phi_sym = format_amplitude_symbol(phase_k)
                is_phi_mass_dep = is_amplitude_mass_dependent(phase_k, config)

                plot_phase_with_ambiguity(
                    ax=ax_bot,
                    masses=masses,
                    y_nom=y_nom_phi,
                    is_mass_dep=is_phi_mass_dep,
                    color=c_wave,
                    label=f"Nominal +{phi_sym}",
                    label_neg=f"Sign Ambiguity -{phi_sym}",
                    mcmc_dict=mcmc_samples.get(phase_k),
                    fit_step=fit_step,
                    show_ambiguity=show_phase_ambiguity,
                    alpha_band=0.32,
                    base_alpha=0.20,
                    max_alpha=0.60,
                    density_band=density_bands,
                )

                configure_phase_axis(ax_bot)
                ax_bot.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
                ax_bot.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11)
                ax_bot.set_ylabel(f"Phase {phi_label} [rad]", fontsize=10.5)
                ax_bot.legend(frameon=True, fontsize=8, loc="best")
                ax_bot.grid(True, linestyle=":", alpha=0.45)

                for fmt in formats:
                    p_indiv_combo = indiv_dir / f"{k}_intensity_and_phase.{fmt}"
                    fig.savefig(p_indiv_combo, bbox_inches="tight")
                    generated.append(p_indiv_combo)
                plt.close(fig)

        # --- 4. Individual Phase Plots ---
        for k in sorted(phases):
            y_nom = nominal_amps[k]
            label = format_amplitude_label(k)
            sym = format_amplitude_symbol(k)
            is_mass_dep = is_amplitude_mass_dependent(k, config)
            c_wave = "black"

            fig, ax = plt.subplots(
                figsize=(7, 4.8), dpi=300,
                gridspec_kw={"top": 0.92, "bottom": 0.12, "left": 0.12, "right": 0.96}
            )
            plot_phase_with_ambiguity(
                ax=ax,
                masses=masses,
                y_nom=y_nom,
                is_mass_dep=is_mass_dep,
                color=c_wave,
                label=f"Nominal +{sym}",
                label_neg=f"Sign Ambiguity -{sym}",
                mcmc_dict=mcmc_samples.get(k),
                fit_step=fit_step,
                show_ambiguity=show_phase_ambiguity,
                alpha_band=0.32,
                base_alpha=0.20,
                max_alpha=0.60,
                density_band=density_bands,
            )

            configure_phase_axis(ax)
            ax.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
            ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=11)
            ax.set_ylabel(f"Phase {label} [rad]", fontsize=11)
            ax.set_title(f"{config.project.name} - Phase {label} ($\\pm$ Sign Ambiguity)", fontsize=12, pad=8)
            ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
            ax.legend(frameon=True, fontsize=9)
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
                                        base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3
                                    )
                                else:
                                    lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], nominal_amps[k], (16.0, 84.0), squared=True)
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
                                        base_alpha=0.22, max_alpha=0.62, fit_step=fit_step, zorder=3
                                    )
                                else:
                                    lo_68, hi_68 = extract_credible_bounds(masses, mcmc_samples[k], nominal_amps[k], (16.0, 84.0), squared=True)
                                    x_step = np.repeat(edges, 2)[1:-1]
                                    lo_s = np.repeat(lo_68, 2)
                                    hi_s = np.repeat(hi_68, 2)
                                    ax_int.fill_between(x_step, lo_s, hi_s, color="black", alpha=0.30, zorder=3)

                        ax_int.set_title(sqr_label, fontsize=10.5, pad=4)
                        ax_int.set_xlim(masses[0] - fit_step * 0.6, masses[-1] + fit_step * 0.6)
                        ax_int.grid(True, linestyle=":", alpha=0.45)

                        # Y-axis labelling: only leftmost column shows y-axis ticks and title
                        if c == 0:
                            ylabel_str = f"Intensity / ({fit_step*1000:.0f} MeV/$c^2$)" if fit_step else "Intensity"
                            ax_int.set_ylabel(ylabel_str, fontsize=9.5)
                            ax_int.tick_params(labelleft=True)
                            if has_phases and ax_phi is not None:
                                configure_phase_axis(ax_phi)
                                ax_phi.set_ylabel("Phase [rad]", fontsize=9.5)
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

                            # Phase plotting in bottom panel
                            pk = find_phase_key(k)
                            if pk:
                                y_nom_phi = nominal_amps[pk]
                                is_phi_mass_dep = is_amplitude_mass_dependent(pk, config)
                                plot_phase_with_ambiguity(
                                    ax=ax_phi,
                                    masses=masses,
                                    y_nom=y_nom_phi,
                                    is_mass_dep=is_phi_mass_dep,
                                    color="black",
                                    mcmc_dict=mcmc_samples.get(pk),
                                    fit_step=fit_step,
                                    show_ambiguity=show_phase_ambiguity,
                                    alpha_band=0.30,
                                    base_alpha=0.20,
                                    max_alpha=0.60,
                                    density_band=density_bands,
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

        a_2_waves = [k for k in magnitudes if k.startswith("a_2_")]
        b_2_waves = [k for k in magnitudes if k.startswith("b_2_")]
        s_waves = [k for k in magnitudes if k.startswith("a_0_") or k.startswith("b_0_")]

        draw_wave_grid(a_2_waves, "D-Waves Positive Reflectivity vs Mass", "overview_D_waves_positive")
        draw_wave_grid(b_2_waves, "D-Waves Negative Reflectivity vs Mass", "overview_D_waves_negative")
        draw_wave_grid(s_waves, "S-Wave vs Mass", "overview_S_waves")
        draw_wave_grid(magnitudes, "All Fitted Waves Overview vs Mass", "overview_all_waves")

    # --- 6. 1D MCMC Posterior Distributions ---
    plot_dist_cfg = getattr(config.plotting, "distributions", None)
    if plot_dist_cfg is None:
        plot_dist_cfg = getattr(config.plotting, "plot_distributions", True)

    if not only_amplitudes and plot_dist_cfg:
        # 6a. Fitted Parameter Posterior Distributions
        if mcmc_params:
            for pname, pvals in mcmc_params.items():
                vals = np.array(pvals, dtype=float)
                if len(vals) < 2 or np.std(vals) < 1e-12:
                    continue
                fig, ax = plt.subplots(
                    figsize=(5.5, 4.2), dpi=300,
                    gridspec_kw={"top": 0.91, "bottom": 0.13, "left": 0.13, "right": 0.95}
                )
                q16, q50, q84 = np.percentile(vals, [16, 50, 84])
                mean_v = np.mean(vals)
                std_v = np.std(vals)

                n_bins = min(50, max(15, int(np.sqrt(len(vals)))))
                ax.hist(vals, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                ax.axvline(q50, color="#1A2E40", linestyle="--", linewidth=1.5, label=f"Median: {q50:.4f}")
                ax.axvline(q16, color="#2B4C7E", linestyle=":", linewidth=1.2, label=f"68% CI: [{q16:.4f}, {q84:.4f}]")
                ax.axvline(q84, color="#2B4C7E", linestyle=":", linewidth=1.2)
                ax.axvspan(q16, q84, color="#2B4C7E", alpha=0.15)

                clean_name = pname
                if pname.startswith("MI_"):
                    m_parts = pname.split("_")
                    if len(m_parts) >= 4:
                        mb = m_parts[1]
                        amp = "_".join(m_parts[2:])
                        clean_name = f"{format_amplitude_label(amp)} @ {mb} GeV"
                    else:
                        clean_name = pname.replace("_", "\\_")
                elif pname.startswith("MD_"):
                    parts = pname[3:].split("_")
                    if len(parts) >= 3 and parts[0] in ("a", "b", "aphi", "bphi"):
                        amp_key = f"{parts[0]}_{parts[1]}_{parts[2]}"
                        sym = format_amplitude_symbol(amp_key)
                        suffix = "_".join(parts[3:]) if len(parts) > 3 else ""
                        clean_name = f"$k({sym})$" if suffix == "k" else f"${sym}$"
                    else:
                        clean_name = pname.replace("_", "\\_")
                else:
                    clean_name = pname.replace("_", "\\_")

                ax.set_title(f"MCMC Posterior Distribution: {clean_name}", fontsize=11, pad=6)
                ax.set_xlabel(f"{clean_name}", fontsize=10)
                ax.set_ylabel("Probability Density", fontsize=10)
                ax.legend(frameon=True, fontsize=8, loc="upper right")
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
                    q16, q50, q84 = np.percentile(vals, [16, 50, 84])
                    n_bins = min(40, max(12, int(np.sqrt(len(vals)))))
                    ax.hist(vals, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                    ax.axvline(q50, color="#1A2E40", linestyle="--", linewidth=1.5)
                    ax.axvline(q16, color="#2B4C7E", linestyle=":", linewidth=1.2)
                    ax.axvline(q84, color="#2B4C7E", linestyle=":", linewidth=1.2)
                    ax.axvspan(q16, q84, color="#2B4C7E", alpha=0.15)

                    raw_short = pname.replace("MD_", "").replace("_a2_1320", "")
                    p_parts = raw_short.split("_")
                    if len(p_parts) >= 3 and p_parts[0] in ("a", "b", "aphi", "bphi"):
                        amp_key = f"{p_parts[0]}_{p_parts[1]}_{p_parts[2]}"
                        sym = format_amplitude_symbol(amp_key)
                        suffix = "_".join(p_parts[3:]) if len(p_parts) > 3 else ""
                        if suffix == "k":
                            short_name = f"$k({sym})$"
                        elif suffix:
                            short_name = f"${sym}_{{{suffix}}}$"
                        else:
                            short_name = f"${sym}$"
                    else:
                        short_name = raw_short.replace("_", "\\_")

                    ax.set_title(f"{short_name} [{q50:.3f}]", fontsize=10, pad=4)
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
                fig, ax = plt.subplots(
                    figsize=(5.5, 4.2), dpi=300,
                    gridspec_kw={"top": 0.91, "bottom": 0.13, "left": 0.13, "right": 0.95}
                )
                q16, q50, q84 = np.percentile(vals, [16, 50, 84])
                nom_yield = sum((v**2) * fit_step for v in nominal_amps.get(k, []))

                n_bins = min(50, max(15, int(np.sqrt(len(vals)))))
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
                    q16, q50, q84 = np.percentile(vals, [16, 50, 84])
                    n_bins = min(40, max(12, int(np.sqrt(len(vals)))))
                    ax.hist(vals, bins=n_bins, color="#B0C4DE", edgecolor="#2B4C7E", alpha=0.7, density=True)
                    ax.axvline(q50, color="#1A2E40", linestyle="--", linewidth=1.5)
                    ax.axvline(q16, color="#2B4C7E", linestyle=":", linewidth=1.2)
                    ax.axvline(q84, color="#2B4C7E", linestyle=":", linewidth=1.2)
                    ax.axvspan(q16, q84, color="#2B4C7E", alpha=0.15)

                    sqr_label = format_amplitude_squared_label(k)
                    ax.set_title(f"{sqr_label} [{q50:.1f}]", fontsize=10, pad=4)
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
                            label=f"Rank 1 ($\chi^2 = {fit.get('chi2', 0):.1f}$)", zorder=5
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
            t_lbl = f"${t_min:.2f} < -t < {t_max:.2f}\\ \\mathrm{{GeV}}^2$"
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
            t_lbl = t_info["label"]

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

            if draw_fits:
                fits_list = t_info.get("fits", [])
                if not fits_list and t_info.get("best_fit"):
                    fits_list = [t_info["best_fit"]]
                for f_idx, fit_obj in enumerate(fits_list):
                    f_mom = fit_obj.get("moments", {}).get(mom_name)
                    if f_mom:
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

            chi2_str = ""
            if t_info.get("best_fit"):
                c2 = t_info["best_fit"].get("chi2", None)
                if c2 is not None:
                    chi2_str = f"\n($\\chi^2 = {c2:.1f}$)"

            ax.set_title(f"{t_lbl}{chi2_str}", fontsize=10, pad=5)
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

    # 2. Multi-panel Overview Grid
    overview_moments = [m for m in selected_moments if is_moment_fitted(m, config)]
    if not overview_moments:
        overview_moments = selected_moments

    K = len(overview_moments)
    if K >= 1:
        ncols = min(4, K)
        nrows = math.ceil(K / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 3.0 * nrows), dpi=250)
        axes_flat = np.array(axes).reshape(-1)

        for idx, mom_name in enumerate(overview_moments):
            ax = axes_flat[idx]
            parts = mom_name.split("_")
            if len(parts) >= 4:
                alpha, l_val, m_val = parts[1], parts[2], parts[3]
                latex_title = f"H_{{{alpha}}}({l_val},{m_val})"
            else:
                latex_title = mom_name

            for t_idx, t_name in enumerate(sorted_tbins):
                t_info = t_data[t_name]
                c = palette[t_idx % len(palette)]
                m_sym = markers[t_idx % len(markers)]
                lbl = t_info["label"]

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
                        fmt=m_sym, color=c, markersize=2.2,
                        capsize=1.0, elinewidth=0.7,
                        label=lbl if idx == 0 else None,
                        alpha=0.9, zorder=4 + t_idx
                    )

                if draw_fits and t_info["best_fit"]:
                    f_mom = t_info["best_fit"].get("moments", {}).get(mom_name)
                    if f_mom:
                        mf = f_mom.get("mass", [])
                        yf = f_mom.get("val", [])
                        if len(mf) > 3:
                            ms, ys = smooth_curve(mf, yf, num_points=200)
                            ax.plot(ms, ys, "-", color=c, linewidth=1.4, alpha=0.85, zorder=3 + t_idx)

            ax.axhline(0, color="gray", linestyle="--", linewidth=0.7, alpha=0.5)
            if fit_mass_min is not None and fit_mass_max is not None:
                pad = fit_mass_step * 0.6
                ax.set_xlim(fit_mass_min - pad, fit_mass_max + pad)

            ax.set_title(f"${latex_title}$", fontsize=10.5, pad=4)
            ax.grid(True, linestyle=":", alpha=0.4)
            ax.tick_params(axis="both", which="major", labelsize=8)
            if (idx + ncols >= K) or (idx >= (nrows - 1) * ncols):
                ax.set_xlabel(r"$M(\eta\pi^0)\ [\mathrm{GeV}/c^2]$", fontsize=9)

        for idx in range(K, len(axes_flat)):
            axes_flat[idx].set_visible(False)

        # Unified figure legend across the top
        handles, labels = axes_flat[0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles, labels,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.995),
                ncol=min(5, len(sorted_tbins)),
                frameon=True,
                fontsize=9
            )

        fig.suptitle(
            f"{config.project.name} - Moments vs Mass Across All $-t$ Bins",
            fontsize=12.5,
            y=0.965 if handles else 0.985,
        )
        plt.tight_layout(rect=[0, 0.02, 1, 0.93 if handles else 0.96])

        for fmt in formats:
            p_overview = out_dir / f"moments_overview_all_tbins.{fmt}"
            fig.savefig(p_overview, bbox_inches="tight")
            generated.append(p_overview)

            # Also copy overview plot to each individual t-bin's plots directory
            for t_name in sorted_tbins:
                t_pdir = t_data[t_name]["rundir"].plots_dir
                if t_pdir.is_dir():
                    dst_p = t_pdir / f"moments_overview_all_tbins.{fmt}"
                    shutil.copyfile(p_overview, dst_p)
                    generated.append(dst_p)

        plt.close(fig)

    print(f"[OK] Generated {len(generated)} multi-t-bin moment comparison plot files in {out_dir}")
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

    return {
        "success": all_success,
        "total_plots": total_plots,
        "tbins": results_by_tbin,
    }


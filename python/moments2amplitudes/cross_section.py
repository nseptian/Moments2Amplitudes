"""
Reflectivity-Separated Cross-Section Evaluator for Resonances in Moments2Amplitudes.

Computes:
1. Wave-by-wave and reflectivity-separated acceptance-corrected event yields:
   Y_w = (1 / Delta_m) * int_{m_min}^{m_max} |A_w(m)|^2 dm
   Y(+) = sum_{a in PosRefl} Y_a
   Y(-) = sum_{b in NegRefl} Y_b
   Y(tot) = Y(+) + Y(-)
2. Differential cross sections:
   dsigma / dt = Y / (den_factor * Luminosity * Delta_t * BranchingRatio)
3. Full parameter error propagation from Minuit / MCMC.
4. Publication-ready cross-section vs -t figures comparing with published Phase-1 PWA and theory bands.
"""

from __future__ import annotations
import math
import os
import re
import json
import subprocess
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Tuple, Any

# Ensure matplotlib has a writable config directory
os.environ.setdefault("MPLCONFIGDIR", "/tmp")

from .schema import FitConfig, CrossSectionConfig
from .rundir import RunDirectory
from .environment import get_root_binary, get_root_subprocess_env


# Published Phase-1 PWA cross-section reference points for a2(1320)
PHASE1_A2_TBIN_CENTERS = [0.15, 0.2625, 0.4125, 0.625, 0.875]
PHASE1_A2_POS_REFL_VALS = [0.187, 0.090, 0.037, 0.038, 0.032]
PHASE1_A2_POS_REFL_STAT = [0.021, 0.012, 0.007, 0.004, 0.005]
PHASE1_A2_POS_REFL_SYST = [0.038, 0.028, 0.013, 0.005, 0.013]
PHASE1_A2_NEG_REFL_VALS = [0.071, 0.026, 0.022, 0.011, 0.005]
PHASE1_A2_NEG_REFL_STAT = [0.013, 0.009, 0.005, 0.003, 0.002]
PHASE1_A2_NEG_REFL_SYST = [0.029, 0.017, 0.010, 0.004, 0.005]
PHASE1_A2_TOTAL_VALS = [0.258, 0.116, 0.059, 0.049, 0.037]
PHASE1_A2_TOTAL_STAT = [0.024, 0.015, 0.008, 0.006, 0.007]
PHASE1_A2_TOTAL_SYST = [0.048, 0.033, 0.016, 0.007, 0.014]


@dataclass
class WaveResonanceYield:
    wave_name: str
    reflectivity: int          # +1 for PosRefl (a_*), -1 for NegRefl (b_*)
    coupling: float
    coupling_error: float
    yield_val: float
    yield_err: float
    coupling_mcmc_err: float = 0.0
    yield_mcmc_err: float = 0.0
    yield_mcmc_err_lo: float = 0.0
    yield_mcmc_err_hi: float = 0.0


@dataclass
class CrossSectionResult:
    resonance: str
    mass: float
    mass_error: float
    width: float
    width_error: float
    t_min: float
    t_max: float
    t_center: float
    delta_t: float
    luminosity: float          # in pb^-1
    branching_ratio: float
    den_factor: float
    pos_refl_yield: float
    pos_refl_yield_err: float
    neg_refl_yield: float
    neg_refl_yield_err: float
    total_yield: float
    total_yield_err: float
    cs_pos: float              # dsigma/dt in mub / GeV^2
    cs_pos_err: float
    cs_neg: float              # dsigma/dt in mub / GeV^2
    cs_neg_err: float
    cs_total: float            # dsigma/dt in mub / GeV^2
    cs_total_err: float
    pos_refl_fraction: float   # Y(+) / Y(tot)
    wave_yields: Dict[str, WaveResonanceYield] = field(default_factory=dict)
    run_dir: str = ""

    # MCMC posterior uncertainty propagation fields
    has_mcmc: bool = False
    mcmc_samples_count: int = 0
    mcmc_file: Optional[str] = None
    pos_refl_yield_mcmc_err: float = 0.0
    pos_refl_yield_mcmc_err_lo: float = 0.0
    pos_refl_yield_mcmc_err_hi: float = 0.0
    neg_refl_yield_mcmc_err: float = 0.0
    neg_refl_yield_mcmc_err_lo: float = 0.0
    neg_refl_yield_mcmc_err_hi: float = 0.0
    total_yield_mcmc_err: float = 0.0
    total_yield_mcmc_err_lo: float = 0.0
    total_yield_mcmc_err_hi: float = 0.0
    cs_pos_mcmc_err: float = 0.0
    cs_pos_mcmc_err_lo: float = 0.0
    cs_pos_mcmc_err_hi: float = 0.0
    cs_neg_mcmc_err: float = 0.0
    cs_neg_mcmc_err_lo: float = 0.0
    cs_neg_mcmc_err_hi: float = 0.0
    cs_total_mcmc_err: float = 0.0
    cs_total_mcmc_err_lo: float = 0.0
    cs_total_mcmc_err_hi: float = 0.0
    pos_refl_fraction_mcmc_err: float = 0.0
    cs_pos_mcmc_mean: float = 0.0
    cs_pos_mcmc_median: float = 0.0
    cs_neg_mcmc_mean: float = 0.0
    cs_neg_mcmc_median: float = 0.0
    cs_total_mcmc_mean: float = 0.0
    cs_total_mcmc_median: float = 0.0
    pos_refl_fraction_mcmc_mean: float = 0.0
    pos_refl_fraction_mcmc_median: float = 0.0
    has_kinematic: bool = False

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        return d


def parse_t_bin_from_string(text: str) -> Optional[Tuple[float, float]]:
    """
    Parse (t_min, t_max) in GeV^2 from canonical strings or directory names:
    - t010020 -> (0.10, 0.20)
    - t020032 -> (0.20, 0.325)
    - t032050 -> (0.325, 0.50)
    - t050075 -> (0.50, 0.75)
    - t075100 -> (0.75, 1.00)
    - 0.1:0.2 or 0.1_0.2 -> (0.1, 0.2)
    """
    if not text:
        return None

    # Known canonical GlueX -t bins
    canon_map = {
        "t010020": (0.10, 0.20),
        "t020032": (0.20, 0.325),
        "t032050": (0.325, 0.50),
        "t050075": (0.50, 0.75),
        "t075100": (0.75, 1.00),
    }
    for k, v in canon_map.items():
        if k in text.lower():
            return v

    # Regex for t(\d+)_(\d+) or t(\d{3})(\d{3})
    m_six = re.search(r't(\d{3})(\d{3})', text)
    if m_six:
        v1 = int(m_six.group(1)) / 100.0 if int(m_six.group(1)) < 100 else int(m_six.group(1)) / 1000.0
        v2 = int(m_six.group(2)) / 100.0 if int(m_six.group(2)) < 100 else int(m_six.group(2)) / 1000.0
        return (v1, v2)

    # Regex for float:float or float_float
    m_pair = re.search(r'(\d+\.\d+)[:_\-](\d+\.\d+)', text)
    if m_pair:
        return (float(m_pair.group(1)), float(m_pair.group(2)))

    return None


def compute_bw_integral(
    coupling: float,
    mass: float,
    width: float,
    m_min: float = 0.8,
    m_max: float = 2.0,
    bin_width: float = 0.04,
    has_kinematic: bool = False,
    m_d1: float = 0.547853,
    m_d2: float = 0.13957,
) -> float:
    """
    Evaluate Breit-Wigner line-shape integral normalized by mass bin width:
    Y = (1 / bin_width) * int_{m_min}^{m_max} [ rho(m) * k^2 / ((m^2 - M^2)^2 + (M*Gamma)^2) ] dm
    """
    if mass <= 0 or width <= 0 or coupling == 0:
        return 0.0

    m_th = m_d1 + m_d2
    m_diff = abs(m_d1 - m_d2)

    def phase_space(m: float) -> float:
        if m <= m_th:
            return 0.0
        q_sq = max(0.0, (m * m - m_th * m_th) * (m * m - m_diff * m_diff))
        q = math.sqrt(q_sq) / (2.0 * m)
        return 2.0 * q / m

    try:
        from scipy.integrate import quad
        def bw(m: float) -> float:
            factor = phase_space(m) if has_kinematic else 1.0
            den = ((m * m - mass * mass) ** 2) + ((mass * width) ** 2)
            return factor * (coupling * coupling) / den
        val, _ = quad(bw, m_min, m_max, limit=100)
        return val / bin_width
    except Exception:
        # High-precision Simpson integration fallback (N=5000)
        import numpy as np
        m_pts = np.linspace(m_min, m_max, 5001)
        den = ((m_pts * m_pts - mass * mass) ** 2) + ((mass * width) ** 2)
        factors = np.array([phase_space(float(m)) for m in m_pts]) if has_kinematic else 1.0
        y_pts = factors * (coupling * coupling) / den
        integral = np.trapz(y_pts, m_pts)
        return float(integral / bin_width)


def extract_resonance_parameters_from_root(
    root_file: Path | str,
    resonance_name: str,
    config: Optional[FitConfig] = None,
) -> Dict[str, Any]:
    """
    Extract coupling constants, mass, and width for a resonance from ROOT nominal_best.root.
    Decoupled from PyROOT ABI via macros/ExportResonanceCouplings.C.
    """
    root_p = Path(root_file).resolve()
    if not root_p.is_file():
        raise FileNotFoundError(f"Fit ROOT file not found: {root_p}")

    repo_root = Path(os.environ.get("MOMS2AMPS", str(Path(__file__).resolve().parent.parent.parent)))
    export_macro = repo_root / "macros" / "ExportResonanceCouplings.C"
    if not export_macro.is_file():
        raise FileNotFoundError(f"Export macro not found at: {export_macro}")

    json_tmp = root_p.parent / f".export_resonance_{resonance_name}.json"
    root_bin = get_root_binary(config)

    cmd = [
        root_bin,
        "-l",
        "-b",
        "-q",
        f'{export_macro}("{root_p}", "{resonance_name}", "{json_tmp}")'
    ]

    env = get_root_subprocess_env(config)

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if proc.returncode != 0 or not json_tmp.is_file():
        raise RuntimeError(f"ExportResonanceCouplings failed (code {proc.returncode}):\n{proc.stderr}\n{proc.stdout}")

    try:
        with open(json_tmp, "r") as f:
            data = json.load(f)
    finally:
        if json_tmp.is_file():
            json_tmp.unlink()

    return data


def find_mcmc_root_file(
    config: Optional[FitConfig] = None,
    rundir: Optional[RunDirectory] = None,
    explicit_file: Optional[Path | str] = None,
) -> Optional[Path]:
    """
    Locate an MCMC propagated ROOT file from explicit argument, config, rundir,
    or standard working directory locations.
    """
    if explicit_file:
        p = Path(explicit_file).resolve()
        if p.is_file():
            return p

    if config and hasattr(config, "cross_section") and config.cross_section:
        if config.cross_section.mcmc_file:
            p = Path(config.cross_section.mcmc_file).resolve()
            if p.is_file():
                return p

    candidates: List[Path] = []
    if rundir:
        candidates.append(rundir.mcmc_propagated_root)
        candidates.append(rundir.root / "mcmc_propagation_spiecewise_results.root")
        candidates.append(rundir.root / "mcmc_propagation_results.root")
        candidates.append(rundir.root / "mcmc" / "mcmc_worker_0.root")
        if config and config.mcmc_propagation and config.mcmc_propagation.output_file:
            candidates.append(rundir.root / config.mcmc_propagation.output_file)
    else:
        if config and config.mcmc_propagation and config.mcmc_propagation.output_file:
            candidates.append(Path(config.mcmc_propagation.output_file).resolve())
        # Current working directory fallbacks only when no explicit rundir is specified
        candidates.append(Path("mcmc_propagation_spiecewise_results.root").resolve())
        candidates.append(Path("mcmc_propagation_results.root").resolve())

    for c in candidates:
        if c and c.is_file():
            return c

    return None


def extract_mcmc_resonance_parameters_from_root(
    mcmc_root_file: Path | str,
    resonance_name: str,
    config: Optional[FitConfig] = None,
    tree_name: str = "mcmc_params",
) -> Dict[str, Any]:
    """
    Extract MCMC sample couplings, mass, width, and fit status for a resonance from ROOT.
    Decoupled from PyROOT ABI via macros/ExportMCMCResonanceCouplings.C.
    """
    mcmc_p = Path(mcmc_root_file).resolve()
    if not mcmc_p.is_file():
        raise FileNotFoundError(f"MCMC ROOT file not found: {mcmc_p}")

    repo_root = Path(os.environ.get("MOMS2AMPS", str(Path(__file__).resolve().parent.parent.parent)))
    export_macro = repo_root / "macros" / "ExportMCMCResonanceCouplings.C"
    if not export_macro.is_file():
        raise FileNotFoundError(f"Export macro not found at: {export_macro}")

    json_tmp = mcmc_p.parent / f".export_mcmc_resonance_{resonance_name}_{os.getpid()}.json"
    root_bin = get_root_binary(config)

    cmd = [
        root_bin,
        "-l",
        "-b",
        "-q",
        f'{export_macro}("{mcmc_p}", "{resonance_name}", "{json_tmp}", "{tree_name}")'
    ]

    env = get_root_subprocess_env(config)

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if proc.returncode != 0 or not json_tmp.is_file():
        raise RuntimeError(f"ExportMCMCResonanceCouplings failed (code {proc.returncode}):\n{proc.stderr}\n{proc.stdout}")

    try:
        with open(json_tmp, "r") as f:
            data = json.load(f)
    finally:
        if json_tmp.is_file():
            json_tmp.unlink()

    return data


def propagate_mcmc_cross_section(
    res: CrossSectionResult,
    mcmc_data: Dict[str, Any],
    mass_limits: Optional[Tuple[float, float]] = None,
    bin_width: Optional[float] = None,
    filter_valid: bool = True,
    mcmc_file_path: Optional[str] = None,
    has_kinematic: bool = False,
) -> CrossSectionResult:
    """
    Propagate MCMC sample couplings into reflectivity-separated yields and cross section.
    Calculates sample distributions, standard deviations, and 68% credible intervals.
    """
    import numpy as np

    waves = mcmc_data.get("waves", [])
    refls = mcmc_data.get("reflectivities", [])
    couplings = mcmc_data.get("couplings", {})
    n_entries = mcmc_data.get("entries_count", 0)

    if n_entries == 0 or not waves or not couplings:
        return res

    # Filter valid entries
    if filter_valid and "is_valid" in mcmc_data:
        val_arr = np.array(mcmc_data["is_valid"], dtype=bool)
        status_arr = np.array(mcmc_data.get("status", [0] * n_entries), dtype=int)
        mask = (val_arr == True) & (status_arr == 0)
        if not np.any(mask):
            mask = (val_arr == True)
        if not np.any(mask):
            mask = np.ones(n_entries, dtype=bool)
    else:
        mask = np.ones(n_entries, dtype=bool)

    valid_indices = np.where(mask)[0]
    n_valid = len(valid_indices)
    if n_valid == 0:
        return res

    m_lims = mass_limits or (0.8, 2.0)
    b_w = bin_width if bin_width is not None and bin_width > 0 else 0.04

    # Resolve mass and width per entry
    masses = np.array(mcmc_data.get("mass", []), dtype=float)
    widths = np.array(mcmc_data.get("width", []), dtype=float)

    if len(masses) != n_entries or np.all(masses <= 0):
        masses = np.full(n_entries, res.mass)
    if len(widths) != n_entries or np.all(widths <= 0):
        widths = np.full(n_entries, res.width)

    val_masses = masses[valid_indices]
    val_widths = widths[valid_indices]

    weights = np.array(mcmc_data.get("weights", []), dtype=float)
    has_weights = (len(weights) == n_entries and np.sum(weights) > 0)
    w_valid = weights[valid_indices] if has_weights else None

    # Evaluate Breit-Wigner integral factor I0 per entry
    is_const_mw = np.allclose(val_masses, val_masses[0]) and np.allclose(val_widths, val_widths[0])
    if is_const_mw:
        i0 = compute_bw_integral(1.0, float(val_masses[0]), float(val_widths[0]), m_lims[0], m_lims[1], b_w, has_kinematic=has_kinematic)
        i0_arr = np.full(n_valid, i0)
    else:
        i0_arr = np.array([
            compute_bw_integral(1.0, float(val_masses[i]), float(val_widths[i]), m_lims[0], m_lims[1], b_w, has_kinematic=has_kinematic)
            for i in range(n_valid)
        ])

    def calc_std(arr: np.ndarray) -> float:
        if w_valid is None:
            return float(np.std(arr))
        mean = np.average(arr, weights=w_valid)
        return float(np.sqrt(np.average((arr - mean)**2, weights=w_valid)))

    def calc_stats(arr: np.ndarray) -> Tuple[float, float, float, float]:
        if w_valid is None:
            med = float(np.median(arr))
            mean = float(np.mean(arr))
            p16 = float(np.percentile(arr, 16))
            p84 = float(np.percentile(arr, 84))
        else:
            mean = float(np.average(arr, weights=w_valid))
            sorter = np.argsort(arr)
            s_arr = arr[sorter]
            s_w = w_valid[sorter]
            cum_w = np.cumsum(s_w)
            tot_w = cum_w[-1]
            if tot_w > 0:
                p16, med, p84 = np.interp([0.16 * tot_w, 0.50 * tot_w, 0.84 * tot_w], cum_w, s_arr)
                p16, med, p84 = float(p16), float(med), float(p84)
            else:
                med = float(np.median(arr))
                p16 = float(np.percentile(arr, 16))
                p84 = float(np.percentile(arr, 84))
        return mean, med, max(0.0, med - p16), max(0.0, p84 - med)

    # Calculate wave yields for each MCMC sample
    y_pos_arr = np.zeros(n_valid)
    y_neg_arr = np.zeros(n_valid)

    for w_name, refl in zip(waves, refls):
        if w_name not in couplings:
            continue
        k_vals = np.array(couplings[w_name], dtype=float)[valid_indices]
        y_w_arr = (k_vals ** 2) * i0_arr

        # Wave statistics
        k_mcmc_err = calc_std(k_vals)
        y_mcmc_err = calc_std(y_w_arr)
        _, y_med, y_lo, y_hi = calc_stats(y_w_arr)

        if w_name in res.wave_yields:
            w_yield = res.wave_yields[w_name]
            w_yield.coupling_mcmc_err = k_mcmc_err
            w_yield.yield_mcmc_err = y_mcmc_err
            w_yield.yield_mcmc_err_lo = y_lo
            w_yield.yield_mcmc_err_hi = y_hi

        if refl == 1:
            y_pos_arr += y_w_arr
        else:
            y_neg_arr += y_w_arr

    y_tot_arr = y_pos_arr + y_neg_arr

    # Cross sections per sample
    denom = res.den_factor * res.luminosity * res.delta_t * res.branching_ratio
    if denom <= 0:
        denom = 1.0

    cs_pos_arr = y_pos_arr / denom
    cs_neg_arr = y_neg_arr / denom
    cs_tot_arr = y_tot_arr / denom
    pos_frac_arr = np.where(y_tot_arr > 0, y_pos_arr / y_tot_arr, 0.5)

    y_pos_mean, y_pos_med, y_pos_lo, y_pos_hi = calc_stats(y_pos_arr)
    y_neg_mean, y_neg_med, y_neg_lo, y_neg_hi = calc_stats(y_neg_arr)
    y_tot_mean, y_tot_med, y_tot_lo, y_tot_hi = calc_stats(y_tot_arr)

    cs_pos_mean, cs_pos_med, cs_pos_lo, cs_pos_hi = calc_stats(cs_pos_arr)
    cs_neg_mean, cs_neg_med, cs_neg_lo, cs_neg_hi = calc_stats(cs_neg_arr)
    cs_tot_mean, cs_tot_med, cs_tot_lo, cs_tot_hi = calc_stats(cs_tot_arr)

    pf_mean, pf_med, pf_lo, pf_hi = calc_stats(pos_frac_arr)

    res.has_mcmc = True
    res.mcmc_samples_count = n_valid
    res.mcmc_file = str(mcmc_file_path) if mcmc_file_path else None

    res.pos_refl_yield_mcmc_err = calc_std(y_pos_arr)
    res.pos_refl_yield_mcmc_err_lo = y_pos_lo
    res.pos_refl_yield_mcmc_err_hi = y_pos_hi

    res.neg_refl_yield_mcmc_err = calc_std(y_neg_arr)
    res.neg_refl_yield_mcmc_err_lo = y_neg_lo
    res.neg_refl_yield_mcmc_err_hi = y_neg_hi

    res.total_yield_mcmc_err = calc_std(y_tot_arr)
    res.total_yield_mcmc_err_lo = y_tot_lo
    res.total_yield_mcmc_err_hi = y_tot_hi

    res.cs_pos_mcmc_err = calc_std(cs_pos_arr)
    res.cs_pos_mcmc_err_lo = cs_pos_lo
    res.cs_pos_mcmc_err_hi = cs_pos_hi

    res.cs_neg_mcmc_err = calc_std(cs_neg_arr)
    res.cs_neg_mcmc_err_lo = cs_neg_lo
    res.cs_neg_mcmc_err_hi = cs_neg_hi

    res.cs_total_mcmc_err = calc_std(cs_tot_arr)
    res.cs_total_mcmc_err_lo = cs_tot_lo
    res.cs_total_mcmc_err_hi = cs_tot_hi

    res.pos_refl_fraction_mcmc_err = calc_std(pos_frac_arr)

    res.cs_pos_mcmc_mean = cs_pos_mean
    res.cs_pos_mcmc_median = cs_pos_med
    res.cs_neg_mcmc_mean = cs_neg_mean
    res.cs_neg_mcmc_median = cs_neg_med
    res.cs_total_mcmc_mean = cs_tot_mean
    res.cs_total_mcmc_median = cs_tot_med

    res.pos_refl_fraction_mcmc_mean = pf_mean
    res.pos_refl_fraction_mcmc_median = pf_med

    return res


def calculate_resonance_cross_section(
    config: FitConfig,
    rundir: RunDirectory,
    resonance_name: Optional[str] = None,
    t_min: Optional[float] = None,
    t_max: Optional[float] = None,
    luminosity: Optional[float] = None,
    branching_ratio: Optional[float] = None,
    mass_limits: Optional[Tuple[float, float]] = None,
    mcmc_root_file: Optional[Path | str] = None,
    propagate_mcmc: Optional[bool] = None,
) -> CrossSectionResult:
    """
    Calculate reflectivity-separated resonance yields and differential cross section dsigma/dt.
    Propagates Minuit parabolic parameter errors and MCMC posterior sample uncertainties.
    """
    cs_cfg = getattr(config, "cross_section", CrossSectionConfig())
    target_res = resonance_name or cs_cfg.resonance
    lumi = luminosity if luminosity is not None else cs_cfg.luminosity
    br = branching_ratio if branching_ratio is not None else cs_cfg.branching_ratio
    m_lims = mass_limits if mass_limits is not None else cs_cfg.mass_limits

    # Resolve mass bin width
    bin_width = cs_cfg.mass_bin_width
    if bin_width is None:
        bin_width = config.dataset.mass_bins.step
    if bin_width is None or bin_width <= 0:
        bin_width = 0.04

    # Resolve -t binning: prioritize explicit t_min/t_max, then auto-detection from rundir/dataset, then config fallback
    t_lo = t_min
    t_hi = t_max

    if t_lo is None or t_hi is None:
        parsed = parse_t_bin_from_string(str(rundir.root.name))
        if not parsed and rundir.root.parent:
            parsed = parse_t_bin_from_string(str(rundir.root.parent.name))
        if not parsed:
            parsed = parse_t_bin_from_string(str(rundir.root))
        if not parsed and getattr(config.dataset, "t_bin", None):
            parsed = parse_t_bin_from_string(config.dataset.t_bin)
        if not parsed:
            resolved_dir = config.dataset.get_moments_dir() if hasattr(config.dataset, "get_moments_dir") else config.dataset.moments_dir
            parsed = parse_t_bin_from_string(resolved_dir)
        if not parsed:
            parsed = parse_t_bin_from_string(config.dataset.var_bin_prefix)

        if parsed:
            t_lo, t_hi = parsed
        else:
            t_lo = cs_cfg.t_min if cs_cfg.t_min is not None else 0.1
            t_hi = cs_cfg.t_max if cs_cfg.t_max is not None else 0.2

    delta_t = t_hi - t_lo
    t_center = (t_lo + t_hi) / 2.0

    # Denominator factor for luminosity units (default pb^-1 -> 1e6 to convert to mub^-1)
    lumi_unit = cs_cfg.luminosity_unit.lower()
    if "pb" in lumi_unit:
        den_factor = 1.0e6
    elif "nb" in lumi_unit:
        den_factor = 1.0e3
    else:
        den_factor = 1.0

    # Locate nominal best root file
    root_target = rundir.nominal_best_root
    if not root_target.is_file():
        legacy_best = rundir.root / "central_best_params.root"
        if legacy_best.is_file():
            root_target = legacy_best
        else:
            raise FileNotFoundError(f"No nominal best fit file found at {root_target}")

    # Extract resonance couplings from ROOT
    res_data = extract_resonance_parameters_from_root(root_target, target_res, config)
    mass = res_data.get("mass", 0.0)
    mass_err = res_data.get("mass_error", 0.0)
    width = res_data.get("width", 0.0)
    width_err = res_data.get("width_error", 0.0)

    # Fallback mass/width from config if zero in tree
    if mass <= 0.0 or width <= 0.0:
        for w in config.waves:
            for r in w.resonances:
                if r.name == target_res:
                    if mass <= 0.0: mass = r.mass
                    if width <= 0.0: width = r.width

    couplings_dict: Dict[str, Dict[str, Any]] = res_data.get("couplings", {})
    if not couplings_dict:
        raise ValueError(f"No couplings found for resonance '{target_res}' in {root_target}")

    # Detect whether kinematic factor sqrt(rho_etapi) is used in the model
    has_kinematic = False
    if getattr(cs_cfg, "has_kinematic", None) is not None:
        has_kinematic = bool(cs_cfg.has_kinematic)
    elif config and getattr(config, "waves", None):
        for w in config.waves:
            matches_res = False
            if getattr(w, "resonances", None):
                for r in w.resonances:
                    if r.name == target_res:
                        matches_res = True
            elif target_res in getattr(w, "name", ""):
                matches_res = True

            if matches_res and getattr(w, "formula", None) and ("rho_etapi" in w.formula or "rho" in w.formula):
                has_kinematic = True
                break

    # Compute yields per wave
    wave_yields: Dict[str, WaveResonanceYield] = {}
    pos_yield = 0.0
    pos_var = 0.0
    neg_yield = 0.0
    neg_var = 0.0

    for w_name, c_info in couplings_dict.items():
        k = float(c_info.get("k", 0.0))
        k_err = float(c_info.get("k_error", 0.0))
        refl = int(c_info.get("reflectivity", 1 if w_name.startswith("a_") else -1))

        # Y_w = Integral(BW) / bin_width
        y_val = compute_bw_integral(k, mass, width, m_lims[0], m_lims[1], bin_width, has_kinematic=has_kinematic)
        # dY / dk = 2 * Y / k  ==>  sigma(Y) = 2 * (Y / |k|) * sigma(k)
        if abs(k) > 1e-9:
            y_err = abs(2.0 * (y_val / k) * k_err)
        else:
            y_err = 0.0

        w_yield = WaveResonanceYield(
            wave_name=w_name,
            reflectivity=refl,
            coupling=k,
            coupling_error=k_err,
            yield_val=y_val,
            yield_err=y_err,
        )
        wave_yields[w_name] = w_yield

        if refl == 1:
            pos_yield += y_val
            pos_var += (y_err * y_err)
        else:
            neg_yield += y_val
            neg_var += (y_err * y_err)

    pos_yield_err = math.sqrt(pos_var)
    neg_yield_err = math.sqrt(neg_var)
    total_yield = pos_yield + neg_yield
    total_yield_err = math.sqrt(pos_var + neg_var)

    # Cross section: dsigma/dt = Y / (den_factor * lumi * delta_t * BR)
    denom = den_factor * lumi * delta_t * br
    if denom <= 0:
        raise ValueError(f"Invalid cross-section denominator: {denom} (lumi={lumi}, delta_t={delta_t}, BR={br})")

    cs_pos = pos_yield / denom
    cs_pos_err = pos_yield_err / denom
    cs_neg = neg_yield / denom
    cs_neg_err = neg_yield_err / denom
    cs_total = total_yield / denom
    cs_total_err = total_yield_err / denom

    pos_frac = (pos_yield / total_yield) if total_yield > 0 else 0.5

    res = CrossSectionResult(
        resonance=target_res,
        mass=mass,
        mass_error=mass_err,
        width=width,
        width_error=width_err,
        t_min=t_lo,
        t_max=t_hi,
        t_center=t_center,
        delta_t=delta_t,
        luminosity=lumi,
        branching_ratio=br,
        den_factor=den_factor,
        pos_refl_yield=pos_yield,
        pos_refl_yield_err=pos_yield_err,
        neg_refl_yield=neg_yield,
        neg_refl_yield_err=neg_yield_err,
        total_yield=total_yield,
        total_yield_err=total_yield_err,
        cs_pos=cs_pos,
        cs_pos_err=cs_pos_err,
        cs_neg=cs_neg,
        cs_neg_err=cs_neg_err,
        cs_total=cs_total,
        cs_total_err=cs_total_err,
        pos_refl_fraction=pos_frac,
        wave_yields=wave_yields,
        run_dir=str(rundir.root),
        has_kinematic=has_kinematic,
    )

    # Propagate MCMC uncertainties if requested
    do_mcmc = propagate_mcmc if propagate_mcmc is not None else getattr(cs_cfg, "propagate_mcmc", True)
    if do_mcmc:
        mcmc_target = find_mcmc_root_file(config=config, rundir=rundir, explicit_file=mcmc_root_file)
        if mcmc_target:
            try:
                mcmc_data = extract_mcmc_resonance_parameters_from_root(mcmc_target, target_res, config)
                res = propagate_mcmc_cross_section(
                    res,
                    mcmc_data,
                    mass_limits=m_lims,
                    bin_width=bin_width,
                    mcmc_file_path=str(mcmc_target),
                    has_kinematic=has_kinematic,
                )
                print(f"[OK] Propagated MCMC statistical uncertainty from: {mcmc_target} ({res.mcmc_samples_count} converged samples)")
            except Exception as e:
                print(f"[WARNING] Could not propagate MCMC cross-section uncertainty from {mcmc_target}: {e}")

    return res


def format_cross_section_table(res: CrossSectionResult) -> str:
    """Format a clean, informative ASCII table of the cross-section calculation."""
    lines = [
        "========================================================================================================================",
        f" Reflectivity-Separated Cross-Section: {res.resonance.upper()}",
        "========================================================================================================================",
        f" Resonance Mass:       {res.mass:.4f} ± {res.mass_error:.4f} GeV",
        f" Resonance Width:      {res.width:.4f} ± {res.width_error:.4f} GeV",
        f" Kinematic Factor:     {'sqrt(rho_etapi) enabled' if res.has_kinematic else 'None (Standard Breit-Wigner)'}",
        f" -t Binning:           {res.t_min:.3f} < -t < {res.t_max:.3f} GeV^2  (Delta t = {res.delta_t:.3f}, Center = {res.t_center:.4f})",
        f" Integrated Lumi:      {res.luminosity:.2f} pb^-1  (factor: {res.den_factor:.0e} to mub^-1)",
        f" Branching Fraction:   {res.branching_ratio:.4e}",
    ]

    if res.has_mcmc:
        mcmc_src = Path(res.mcmc_file).name if res.mcmc_file else "MCMC Tree"
        lines.append(f" MCMC Propagation:     Enabled ({res.mcmc_samples_count} samples from {mcmc_src})")
    else:
        lines.append(" MCMC Propagation:     Disabled / Not available (Minuit parabolic error propagation)")

    lines.extend([
        "------------------------------------------------------------------------------------------------------------------------",
        " Wave Breakdown & Acceptance-Corrected Yields:",
    ])

    if res.has_mcmc:
        lines.append(f" {'Wave':<10} | {'Refl':^6} | {'Coupling k':^12} | {'Minuit Err':^10} | {'MCMC Err':^10} | {'Yield':^13} | {'Minuit Err':^10} | {'MCMC Err':^10}")
        lines.append("------------------------------------------------------------------------------------------------------------------------")
        for w_name, w_info in sorted(res.wave_yields.items(), key=lambda x: (x[1].reflectivity, x[0])):
            refl_str = "+1" if w_info.reflectivity == 1 else "-1"
            lines.append(
                f" {w_name:<10} | {refl_str:^6} | {w_info.coupling:12.4f} | {w_info.coupling_error:10.4f} | {w_info.coupling_mcmc_err:10.4f} | "
                f"{w_info.yield_val:13.1f} | {w_info.yield_err:10.1f} | {w_info.yield_mcmc_err:10.1f}"
            )
    else:
        lines.append(f" {'Wave':<10} | {'Refl':^6} | {'Coupling k':^14} | {'Error':^12} | {'Yield':^14} | {'Yield Err':^12}")
        lines.append("------------------------------------------------------------------------------------------------------------------------")
        for w_name, w_info in sorted(res.wave_yields.items(), key=lambda x: (x[1].reflectivity, x[0])):
            refl_str = "+1" if w_info.reflectivity == 1 else "-1"
            lines.append(
                f" {w_name:<10} | {refl_str:^6} | {w_info.coupling:14.4f} | {w_info.coupling_error:12.4f} | "
                f"{w_info.yield_val:14.1f} | {w_info.yield_err:12.1f}"
            )

    lines.append("------------------------------------------------------------------------------------------------------------------------")
    lines.append(" Reflectivity-Separated Cross-Section Results:")

    if res.has_mcmc:
        lines.extend([
            f"  * Positive Reflectivity (+): Yield = {res.pos_refl_yield:12.1f} ± {res.pos_refl_yield_err:8.1f} (Minuit) ± {res.pos_refl_yield_mcmc_err:8.1f} (MCMC)",
            f"                               dsigma/dt = {res.cs_pos:.5f} ± {res.cs_pos_err:.5f} (Minuit) ± {res.cs_pos_mcmc_err:.5f} (MCMC) mub / GeV^2",
            f"                               [MCMC 68% CI: {res.cs_pos_mcmc_median:.5f} (-{res.cs_pos_mcmc_err_lo:.5f}, +{res.cs_pos_mcmc_err_hi:.5f})]",
            f"  * Negative Reflectivity (-): Yield = {res.neg_refl_yield:12.1f} ± {res.neg_refl_yield_err:8.1f} (Minuit) ± {res.neg_refl_yield_mcmc_err:8.1f} (MCMC)",
            f"                               dsigma/dt = {res.cs_neg:.5f} ± {res.cs_neg_err:.5f} (Minuit) ± {res.cs_neg_mcmc_err:.5f} (MCMC) mub / GeV^2",
            f"                               [MCMC 68% CI: {res.cs_neg_mcmc_median:.5f} (-{res.cs_neg_mcmc_err_lo:.5f}, +{res.cs_neg_mcmc_err_hi:.5f})]",
            f"  * Total Cross Section:       Yield = {res.total_yield:12.1f} ± {res.total_yield_err:8.1f} (Minuit) ± {res.total_yield_mcmc_err:8.1f} (MCMC)",
            f"                               dsigma/dt = {res.cs_total:.5f} ± {res.cs_total_err:.5f} (Minuit) ± {res.cs_total_mcmc_err:.5f} (MCMC) mub / GeV^2",
            f"                               [MCMC 68% CI: {res.cs_total_mcmc_median:.5f} (-{res.cs_total_mcmc_err_lo:.5f}, +{res.cs_total_mcmc_err_hi:.5f})]",
            f"  * Positive Reflectivity Fraction: {res.pos_refl_fraction * 100:.2f} % (Nominal) | {res.pos_refl_fraction_mcmc_mean * 100:.2f} ± {res.pos_refl_fraction_mcmc_err * 100:.2f} % (MCMC)",
        ])
    else:
        lines.extend([
            f"  * Positive Reflectivity (+): Yield = {res.pos_refl_yield:12.1f} ± {res.pos_refl_yield_err:10.1f}",
            f"                               dsigma/dt = {res.cs_pos:.5f} ± {res.cs_pos_err:.5f} mub / GeV^2",
            f"  * Negative Reflectivity (-): Yield = {res.neg_refl_yield:12.1f} ± {res.neg_refl_yield_err:10.1f}",
            f"                               dsigma/dt = {res.cs_neg:.5f} ± {res.cs_neg_err:.5f} mub / GeV^2",
            f"  * Total Cross Section:       Yield = {res.total_yield:12.1f} ± {res.total_yield_err:10.1f}",
            f"                               dsigma/dt = {res.cs_total:.5f} ± {res.cs_total_err:.5f} mub / GeV^2",
            f"  * Positive Reflectivity Fraction: {res.pos_refl_fraction * 100:.2f} %",
        ])

    lines.append("========================================================================================================================")
    return "\n".join(lines)


def save_cross_section_reports(res: CrossSectionResult, output_dir: Path | str) -> Tuple[Path, Path, Path]:
    """Save ASCII, Markdown, and JSON cross-section summaries."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    txt_p = out_dir / f"cross_section_{res.resonance}.txt"
    md_p = out_dir / f"cross_section_{res.resonance}.md"
    json_p = out_dir / f"cross_section_{res.resonance}.json"

    # TXT
    txt_content = format_cross_section_table(res)
    txt_p.write_text(txt_content)

    # MD
    md_lines = [
        f"# Cross-Section Report: {res.resonance}",
        "",
        f"- **Resonance:** `{res.resonance}` (Mass = {res.mass:.4f} GeV, Width = {res.width:.4f} GeV)",
        f"- **Kinematics:** ${res.t_min:.3f} < -t < {res.t_max:.3f}$ GeV$^2$ (Center = {res.t_center:.4f} GeV$^2$, $\\Delta t = {res.delta_t:.3f}$)",
        f"- **Luminosity:** ${res.luminosity:.2f}\\ \\text{{pb}}^{{-1}}$, **Total BR:** `{res.branching_ratio:.4e}`",
    ]

    if res.has_mcmc:
        md_lines.append(f"- **MCMC Propagation:** `{res.mcmc_samples_count}` converged posterior fits from `{Path(res.mcmc_file).name if res.mcmc_file else 'MCMC Tree'}`")

    md_lines.extend([
        "",
        "## Cross-Section Results",
        "",
    ])

    if res.has_mcmc:
        md_lines.extend([
            "| Quantity | Yield (Events) | Minuit Err | MCMC Stat Err | $d\\sigma/dt$ ($\\mu\\text{b}/\\text{GeV}^2$) | Minuit Err | MCMC Stat Err | MCMC 68% CI | Fraction |",
            "|---|---|---|---|---|---|---|---|---|",
            f"| **Positive Reflectivity (+)** | {res.pos_refl_yield:.1f} | ±{res.pos_refl_yield_err:.1f} | ±{res.pos_refl_yield_mcmc_err:.1f} | **{res.cs_pos:.5f}** | ±{res.cs_pos_err:.5f} | **±{res.cs_pos_mcmc_err:.5f}** | [-{res.cs_pos_mcmc_err_lo:.5f}, +{res.cs_pos_mcmc_err_hi:.5f}] | {res.pos_refl_fraction*100:.1f}% |",
            f"| **Negative Reflectivity (-)** | {res.neg_refl_yield:.1f} | ±{res.neg_refl_yield_err:.1f} | ±{res.neg_refl_yield_mcmc_err:.1f} | **{res.cs_neg:.5f}** | ±{res.cs_neg_err:.5f} | **±{res.cs_neg_mcmc_err:.5f}** | [-{res.cs_neg_mcmc_err_lo:.5f}, +{res.cs_neg_mcmc_err_hi:.5f}] | {(1.0-res.pos_refl_fraction)*100:.1f}% |",
            f"| **Total Cross Section** | {res.total_yield:.1f} | ±{res.total_yield_err:.1f} | ±{res.total_yield_mcmc_err:.1f} | **{res.cs_total:.5f}** | ±{res.cs_total_err:.5f} | **±{res.cs_total_mcmc_err:.5f}** | [-{res.cs_total_mcmc_err_lo:.5f}, +{res.cs_total_mcmc_err_hi:.5f}] | 100.0% |",
            "",
            "## Wave Decomposition",
            "",
            "| Wave | Reflectivity | Coupling $k$ | Minuit Err | MCMC Err | Yield | Minuit Err | MCMC Err |",
            "|---|---|---|---|---|---|---|---|",
        ])
        for w_name, w_info in sorted(res.wave_yields.items(), key=lambda x: (x[1].reflectivity, x[0])):
            refl_str = "$+1$" if w_info.reflectivity == 1 else "$-1$"
            md_lines.append(f"| `{w_name}` | {refl_str} | {w_info.coupling:.4f} | {w_info.coupling_error:.4f} | {w_info.coupling_mcmc_err:.4f} | {w_info.yield_val:.1f} | ±{w_info.yield_err:.1f} | ±{w_info.yield_mcmc_err:.1f} |")
    else:
        md_lines.extend([
            "| Quantity | Yield (Events) | $d\\sigma/dt$ ($\\mu\\text{b}/\\text{GeV}^2$) | Fraction |",
            "|---|---|---|---|",
            f"| **Positive Reflectivity (+)** | {res.pos_refl_yield:.1f} ± {res.pos_refl_yield_err:.1f} | **{res.cs_pos:.5f} ± {res.cs_pos_err:.5f}** | {res.pos_refl_fraction*100:.1f}% |",
            f"| **Negative Reflectivity (-)** | {res.neg_refl_yield:.1f} ± {res.neg_refl_yield_err:.1f} | **{res.cs_neg:.5f} ± {res.cs_neg_err:.5f}** | {(1.0-res.pos_refl_fraction)*100:.1f}% |",
            f"| **Total Cross Section** | {res.total_yield:.1f} ± {res.total_yield_err:.1f} | **{res.cs_total:.5f} ± {res.cs_total_err:.5f}** | 100.0% |",
            "",
            "## Wave Decomposition",
            "",
            "| Wave | Reflectivity | Coupling $k$ | Error | Yield |",
            "|---|---|---|---|---|",
        ])
        for w_name, w_info in sorted(res.wave_yields.items(), key=lambda x: (x[1].reflectivity, x[0])):
            refl_str = "$+1$" if w_info.reflectivity == 1 else "$-1$"
            md_lines.append(f"| `{w_name}` | {refl_str} | {w_info.coupling:.4f} | {w_info.coupling_error:.4f} | {w_info.yield_val:.1f} ± {w_info.yield_err:.1f} |")

    md_p.write_text("\n".join(md_lines))

    # JSON
    with open(json_p, "w") as f:
        json.dump(res.to_dict(), f, indent=2)

    return txt_p, md_p, json_p


def plot_cross_section_vs_t(
    results: List[CrossSectionResult],
    output_dir: Path | str,
    formats: Tuple[str, ...] = ("png", "pdf"),
    phase1_comparison: bool = True,
    theory_band_file: Optional[str] = None,
    style: Optional[str] = "default",
) -> List[Path]:
    """
    Plot differential cross-section dsigma/dt vs -t separated by reflectivity.
    Overlays published Phase-1 PWA benchmarks and theoretical model bands.
    Plots MCMC propagated posterior statistical uncertainties when available.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    generated: List[Path] = []

    if not results:
        return generated

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.patches import Rectangle
    except ImportError:
        print("[WARNING] matplotlib not found. Skipping cross-section plot.")
        return generated

    try:
        from .plotter import apply_plot_style
        apply_plot_style(style)
    except Exception:
        pass

    # Sort results by -t center
    sorted_res = sorted(results, key=lambda x: x.t_center)
    t_centers = [r.t_center for r in sorted_res]
    cs_pos = [r.cs_pos for r in sorted_res]
    cs_neg = [r.cs_neg for r in sorted_res]
    cs_tot = [r.cs_total for r in sorted_res]

    # Use MCMC propagated uncertainty if available, fallback to Minuit
    has_mcmc_any = any(r.has_mcmc for r in sorted_res)
    cs_pos_err = [(r.cs_pos_mcmc_err if r.has_mcmc else r.cs_pos_err) for r in sorted_res]
    cs_neg_err = [(r.cs_neg_mcmc_err if r.has_mcmc else r.cs_neg_err) for r in sorted_res]
    cs_tot_err = [(r.cs_total_mcmc_err if r.has_mcmc else r.cs_total_err) for r in sorted_res]

    fig, ax = plt.subplots(figsize=(12, 6.5), dpi=300)

    # 1. Theoretical prediction band
    th_file = theory_band_file
    if not th_file or not Path(th_file).is_file():
        candidates = [
            Path(__file__).resolve().parent / "data" / "theory_pred" / "Bands_TMD_A2_cs.txt",
            Path(__file__).resolve().parent.parent.parent / "data" / "theory_pred" / "Bands_TMD_A2_cs.txt",
            Path("/d/home/septian/EtaPi0Analysis/run_legacy/theory_pred/Bands_TMD_A2_cs.txt"),
        ]
        for cand in candidates:
            if cand.is_file():
                th_file = str(cand)
                break

    th_data = None
    if th_file and Path(th_file).is_file():
        try:
            th_data = np.loadtxt(th_file, skiprows=1)
        except Exception as e:
            print(f"[WARNING] Could not load theory band from {th_file}: {e}")

    # 2. Overlay published Phase-1 benchmark values if requested and resonance is a2(1320)
    has_phase1 = phase1_comparison and any("a2" in r.resonance.lower() for r in sorted_res)
    if has_phase1:
        x_offset = 0.01
        p1_x = [x - x_offset for x in PHASE1_A2_TBIN_CENTERS]
        box_half_w = 0.007

        # Statistical uncertainties as shaded transparent boxes
        for i, x in enumerate(p1_x):
            rect_pos = Rectangle(
                (x - box_half_w, PHASE1_A2_POS_REFL_VALS[i] - PHASE1_A2_POS_REFL_STAT[i]),
                2 * box_half_w, 2 * PHASE1_A2_POS_REFL_STAT[i],
                facecolor="red", edgecolor="red", alpha=0.20, linewidth=1, zorder=2
            )
            ax.add_patch(rect_pos)

            rect_neg = Rectangle(
                (x - box_half_w, PHASE1_A2_NEG_REFL_VALS[i] - PHASE1_A2_NEG_REFL_STAT[i]),
                2 * box_half_w, 2 * PHASE1_A2_NEG_REFL_STAT[i],
                facecolor="blue", edgecolor="blue", alpha=0.20, linewidth=1, zorder=2
            )
            ax.add_patch(rect_neg)

            rect_tot = Rectangle(
                (x - box_half_w, PHASE1_A2_TOTAL_VALS[i] - PHASE1_A2_TOTAL_STAT[i]),
                2 * box_half_w, 2 * PHASE1_A2_TOTAL_STAT[i],
                facecolor="black", edgecolor="black", alpha=0.15, linewidth=1, zorder=2
            )
            ax.add_patch(rect_tot)

        # Published PWA central values (all circles)
        ax.scatter(p1_x, PHASE1_A2_TOTAL_VALS, marker="o", color="black", label="Total (Published PWA)", zorder=3)
        ax.scatter(p1_x, PHASE1_A2_POS_REFL_VALS, marker="o", color="red", label="Positive Reflectivity (Published PWA)", zorder=3)
        ax.scatter(p1_x, PHASE1_A2_NEG_REFL_VALS, marker="o", color="blue", label="Negative Reflectivity (Published PWA)", zorder=3)

        # Published PWA systematic uncertainties (open circles)
        ax.errorbar(
            p1_x, PHASE1_A2_TOTAL_VALS, yerr=PHASE1_A2_TOTAL_SYST,
            fmt="o", mfc="none", mec="black", ecolor="black", capsize=3, linestyle="none", alpha=0.8, zorder=3
        )
        ax.errorbar(
            p1_x, PHASE1_A2_POS_REFL_VALS, yerr=PHASE1_A2_POS_REFL_SYST,
            fmt="o", mfc="none", mec="red", ecolor="red", capsize=3, linestyle="none", alpha=0.8, zorder=3
        )
        ax.errorbar(
            p1_x, PHASE1_A2_NEG_REFL_VALS, yerr=PHASE1_A2_NEG_REFL_SYST,
            fmt="o", mfc="none", mec="blue", ecolor="blue", capsize=3, linestyle="none", alpha=0.8, zorder=3
        )

    # 3. Current Extraction Points (inverted triangles with error bars)
    curr_x = [x + (0.01 if has_phase1 else 0.0) for x in t_centers]

    ax.errorbar(
        curr_x, cs_tot, yerr=cs_tot_err,
        fmt="v", markersize=7, color="black", ecolor="black", capsize=3, elinewidth=1.6, alpha=0.85,
        label="Total (Moment Analysis)", zorder=5
    )
    ax.errorbar(
        curr_x, cs_pos, yerr=cs_pos_err,
        fmt="v", markersize=7, color="red", ecolor="red", capsize=3, elinewidth=1.6, alpha=0.85,
        label="Positive Reflectivity (Moment Analysis)", zorder=5
    )
    ax.errorbar(
        curr_x, cs_neg, yerr=cs_neg_err,
        fmt="v", markersize=7, color="blue", ecolor="blue", capsize=3, elinewidth=1.6, alpha=0.85,
        label="Negative Reflectivity (Moment Analysis)", zorder=5
    )

    # 4. Theoretical Model Band (drawn with label matching notebook)
    if th_data is not None:
        t_th = th_data[:, 0]
        cs_min = th_data[:, 1]
        cs_max = th_data[:, 2] if th_data.shape[1] >= 3 else th_data[:, 1]
        ax.fill_between(
            t_th, cs_min, cs_max,
            color="black", alpha=0.3,
            label="Theoretical Model [Phys. Rev. D 102, 014003]",
            zorder=1
        )

    ax.set_xlabel(r"$-t$ ($\mathrm{GeV}^2$)", fontsize=16)
    ax.set_ylabel(r"$\frac{d\sigma}{dt}\ (\frac{\mu\mathrm{b}}{\mathrm{GeV}^2})$", fontsize=16)
    max_t_extent = max(t_centers + (PHASE1_A2_TBIN_CENTERS if has_phase1 else []))
    ax.set_xlim(0.08, max_t_extent + 0.12)
    y_max = max(cs_tot + (PHASE1_A2_TOTAL_VALS if has_phase1 else [0.3])) * 1.35
    ax.set_ylim(0, y_max)
    ax.legend(frameon=True, fontsize=11, loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.5)

    base_name = f"{sorted_res[0].resonance}_CrossSection_vs_t"
    for fmt in formats:
        p = out_dir / f"{base_name}.{fmt}"
        fig.savefig(p, bbox_inches="tight")
        generated.append(p)
    plt.close(fig)

    print(f"[OK] Saved cross-section vs -t plot to: {out_dir / base_name}.png")
    return generated
    ax.set_title(f"Differential Cross Section vs $-t$ for {res_label}", fontsize=14, fontweight="bold")
    ax.set_xlim(0.08, max(t_centers) + 0.15)
    ax.set_ylim(0, max(cs_tot + ([0.3] if not cs_tot else [])) * 1.35)
    ax.legend(frameon=True, fontsize=10, loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.5)

    base_name = f"{sorted_res[0].resonance}_CrossSection_vs_t"
    for fmt in formats:
        p = out_dir / f"{base_name}.{fmt}"
        fig.savefig(p, bbox_inches="tight")
        generated.append(p)
    plt.close(fig)

    print(f"[OK] Saved cross-section vs -t plot to: {out_dir / base_name}.png")
    return generated


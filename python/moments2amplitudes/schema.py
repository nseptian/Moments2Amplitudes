"""
Schema definitions for Moments2Amplitudes fit configurations.
Uses standard Python dataclasses for lightweight, dependency-free validation.
"""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from enum import Enum
from pathlib import Path
from typing import List, Dict, Optional, Union, Any
import numpy as np
from .environment import get_default_max_workers


class WaveType(str, Enum):
    PIECEWISE = "piecewise"              # Mass-independent per bin
    MASS_DEPENDENT = "mass_dependent"    # Mass-dependent analytic model


class ModelType(str, Enum):
    POLYNOMIAL = "polynomial"
    FLATTE = "flatte"
    BREIT_WIGNER = "breit_wigner"
    TWO_BREIT_WIGNER = "two_breit_wigner"
    FLATTE_BREIT_WIGNER = "flatte_breit_wigner"
    BARRIER_BREIT_WIGNER = "barrier_breit_wigner"
    BARRIER_TWO_BREIT_WIGNER = "barrier_two_breit_wigner"
    CUSTOM = "custom"


class InitStrategy(str, Enum):
    DEFAULT = "default"
    FIXED_ZERO = "fixed_zero"
    FIXED_VALUE = "fixed_value"
    RANDOM = "random"
    FROM_FILE = "from_file"


@dataclass
class EnvironmentConfig:
    brufit_dir: Optional[str] = None
    root_dir: Optional[str] = None
    extra_root_includes: List[str] = field(default_factory=list)


@dataclass
class MassBinsConfig:
    min: Optional[float] = None
    max: Optional[float] = None
    step: Optional[float] = None
    bins: Optional[List[float]] = None

    def get_bin_centers(self) -> List[float]:
        if self.bins is not None and len(self.bins) > 0:
            return [round(b, 6) for b in self.bins]
        if self.min is not None and self.max is not None and self.step is not None:
            # Generate bin centers inclusive of max within floating point epsilon
            n_bins = int(round((self.max - self.min) / self.step)) + 1
            vals = [round(self.min + i * self.step, 6) for i in range(n_bins)]
            return vals
        raise ValueError("MassBinsConfig must specify either (min, max, step) or an explicit list of 'bins'")


@dataclass
class DatasetConfig:
    moments_dir: str
    t_bin: Optional[str] = None
    t_bins: List[str] = field(default_factory=list)
    results_filename: str = "ResultsBruMcmcFixedCovariance.root"
    bin_prefix: str = "Mpi0eta"
    var_bin_prefix: str = "t010030"
    mass_bins: MassBinsConfig = field(default_factory=MassBinsConfig)
    use_covariance: bool = True

    def get_t_bins(self) -> List[str]:
        """Return the list of configured t-bins. If t_bins is empty, falls back to [t_bin] if set."""
        if self.t_bins:
            return list(self.t_bins)
        tbins_alias = getattr(self, "tbins", None)
        if tbins_alias:
            return list(tbins_alias)
        if self.t_bin:
            return [self.t_bin]
        return []

    def get_moments_dir(self, target_t_bin: Optional[str] = None) -> str:
        """
        Return the resolved directory containing the moment mass-bin folders.
        If target_t_bin (or self.t_bin) is configured and not already the terminal
        directory in moments_dir, it resolves moments_dir / t_bin.
        """
        t_choice = target_t_bin if target_t_bin is not None else self.t_bin
        if not t_choice:
            return self.moments_dir

        t_clean = str(t_choice).strip()
        if not t_clean:
            return self.moments_dir

        t_path = Path(t_clean)
        if t_path.is_absolute():
            return str(t_path)

        base_path = Path(self.moments_dir)
        # Avoid duplicating if moments_dir already ends with t_clean (e.g. ".../t010020")
        if base_path.name == t_clean:
            return str(base_path)

        return str(base_path / t_clean)


@dataclass
class Chi2Config:
    included_L: List[int] = field(default_factory=lambda: [0, 2, 4])
    max_l: int = 4
    include_all: bool = False
    no_use: List[str] = field(default_factory=lambda: ["H_3"])

    def __init__(
        self,
        included_L: Optional[List[int]] = None,
        max_l: int = 4,
        include_all: bool = False,
        no_use: Optional[List[str]] = None,
        **kwargs,
    ):
        self.included_L = included_L if included_L is not None else [0, 2, 4]
        self.max_l = max_l
        self.include_all = include_all

        if no_use is not None:
            self.no_use = no_use
        elif 3 in self.included_L:
            self.no_use = []
        else:
            self.no_use = ["H_3"]


@dataclass
class ResonanceConfig:
    name: str
    mass: float = 1.318
    width: Optional[float] = None
    mass_fixed: bool = False
    width_fixed: bool = False
    mass_limits: Optional[List[float]] = None
    width_limits: Optional[List[float]] = None
    g_s: Optional[float] = None
    g_s_fixed: bool = False
    g_s_limits: Optional[List[float]] = None
    g_etapi: Optional[float] = None
    g_etapi_fixed: bool = False
    g_etapi_limits: Optional[List[float]] = None
    g_KK: Optional[float] = None
    g_KK_fixed: bool = False
    g_KK_limits: Optional[List[float]] = None
    g_etaprimepi: Optional[float] = None
    g_etaprimepi_fixed: bool = False
    g_etaprimepi_limits: Optional[List[float]] = None
    gaussian_constraints: Dict[str, List[float]] = field(default_factory=dict)  # param -> [mean, sigma]
    phase: Optional[float] = None
    phase_fixed: bool = False
    fixed_phase: Optional[float] = None
    phase_limits: Optional[List[float]] = None
    shared_phase: Optional[Union[bool, str, int]] = None
    fixed_zero: bool = False
    k_fixed: bool = False

    def __post_init__(self):
        if self.g_etapi is None and self.g_s is not None:
            self.g_etapi = self.g_s
            self.g_etapi_fixed = self.g_s_fixed
            self.g_etapi_limits = self.g_s_limits
        elif self.g_s is None and self.g_etapi is not None:
            self.g_s = self.g_etapi
            self.g_s_fixed = self.g_etapi_fixed
            self.g_s_limits = self.g_etapi_limits


@dataclass
class ParameterInitRule:
    strategy: InitStrategy = InitStrategy.DEFAULT
    min_val: float = 0.0
    max_val: float = 350.0
    value: float = 0.0
    is_fixed: bool = False
    spread: Optional[float] = None


@dataclass
class WaveInitConfig:
    strategy: InitStrategy = InitStrategy.DEFAULT
    fixed_zero: bool = False
    file_path: Optional[str] = None
    option: int = 0
    shared_phase: Union[bool, str, int] = False
    magnitude: ParameterInitRule = field(
        default_factory=lambda: ParameterInitRule(strategy=InitStrategy.RANDOM, min_val=0.0, max_val=350.0)
    )
    coupling: ParameterInitRule = field(
        default_factory=lambda: ParameterInitRule(strategy=InitStrategy.DEFAULT, min_val=-80.0, max_val=80.0)
    )
    phase: ParameterInitRule = field(
        default_factory=lambda: ParameterInitRule(strategy=InitStrategy.FIXED_VALUE, value=0.0, is_fixed=False)
    )
    poly_scale: float = 250.0
    poly_growth: float = 3.8
    spread: Optional[float] = None
    magnitude_spread: Optional[float] = None
    coupling_spread: Optional[float] = None
    phase_spread: Optional[float] = None
    abs_floor: float = 10.0
    coupling_floor: float = 5.0
    phase_limits_specified: bool = False


@dataclass
class WaveConfig:
    l: int
    name: str = ""
    type: WaveType = WaveType.MASS_DEPENDENT
    model: ModelType = ModelType.BREIT_WIGNER
    formula: Optional[str] = None
    parameters: List[str] = field(default_factory=list)
    poly_order: int = -1
    resonances: List[ResonanceConfig] = field(default_factory=list)
    init: WaveInitConfig = field(default_factory=WaveInitConfig)
    phase: Optional[float] = None
    phase_fixed: bool = False
    shared_phase: Optional[Union[bool, str, int]] = None
    barrier: bool = False
    barrier_l: Optional[int] = None
    barrier_radius: float = 5.0677

    def __post_init__(self):
        if not self.name:
            if self.l == 0:
                self.name = "S-wave"
            elif self.l == 1:
                self.name = "P-wave"
            elif self.l == 2:
                self.name = "D-wave"
            else:
                self.name = f"L{self.l}-wave"

    @property
    def option(self) -> int:
        """BruFit parameter sharing option (0, 1, or 2), automatically derived from resonances."""
        return self.init.option


@dataclass
class MinimizerConfig:
    algorithm: str = "Minuit2"
    strategy: int = 1
    tolerance: float = 1.0e-3
    max_calls: int = 1000000
    max_iterations: int = 100000
    print_level: int = 0
    reject_at_limit: bool = True
    limit_tolerance: float = 1.0e-5


def parse_seed_range(val: Any) -> Optional[Tuple[int, int]]:
    """
    Parse a seed range specification into an inclusive (start, end) integer tuple.
    Supported formats:
      - [5001, 10000] or (5001, 10000)
      - "5001:10000"
      - "5001-10000"
      - "5001..10000"
      - "5001, 10000"
      - "5001 10000"
    """
    if val is None:
        return None
    if isinstance(val, (list, tuple)):
        if len(val) != 2:
            raise ValueError(f"Seed range list/tuple must have exactly 2 elements, got {val}")
        start, end = int(val[0]), int(val[1])
        if start < 1 or end < 1:
            raise ValueError(f"Seed values must be positive integers >= 1, got ({start}, {end})")
        if start > end:
            raise ValueError(f"Seed range start ({start}) must be <= end ({end})")
        return (start, end)
    if isinstance(val, str):
        s = val.strip()
        for sep in [":", "..", ",", "-"]:
            if sep in s:
                parts = [p.strip() for p in s.split(sep, 1)]
                if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                    start, end = int(parts[0]), int(parts[1])
                    if start < 1 or end < 1:
                        raise ValueError(f"Seed values must be positive integers >= 1, got ({start}, {end})")
                    if start > end:
                        raise ValueError(f"Seed range start ({start}) must be <= end ({end})")
                    return (start, end)
        if " " in s:
            parts = s.split()
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                start, end = int(parts[0]), int(parts[1])
                if start < 1 or end < 1:
                    raise ValueError(f"Seed values must be positive integers >= 1, got ({start}, {end})")
                if start > end:
                    raise ValueError(f"Seed range start ({start}) must be <= end ({end})")
                return (start, end)
    raise ValueError(f"Invalid seed range format: '{val}'. Expected '[start, end]' or 'start:end'.")


@dataclass
class MultiStartConfig:
    enabled: bool = True
    n_starts: int = 100
    n_threads: Optional[int] = None
    n_tbin_workers: Optional[int] = None
    output_dir: str = "logs"
    plateau_delta_chi2: float = 0.05
    seed_offset: int = 1
    seed_start: Optional[int] = None
    seed_end: Optional[int] = None
    seed_range: Optional[Union[List[int], Tuple[int, int], str]] = None
    skip_existing: bool = False
    n_lowest_chi2: int = 10

    def get_seed_range(self, n_starts_override: Optional[int] = None) -> Tuple[int, int]:
        """Return resolved (start_seed, end_seed) tuple (inclusive)."""
        if self.seed_range is not None:
            return parse_seed_range(self.seed_range)
        if self.seed_start is not None and self.seed_end is not None:
            if self.seed_start < 1 or self.seed_end < 1:
                raise ValueError(f"Seed values must be positive integers >= 1, got ({self.seed_start}, {self.seed_end})")
            if self.seed_start > self.seed_end:
                raise ValueError(f"seed_start ({self.seed_start}) must be <= seed_end ({self.seed_end})")
            return (self.seed_start, self.seed_end)
        start = self.seed_start if self.seed_start is not None else self.seed_offset
        count = n_starts_override if n_starts_override is not None else self.n_starts
        return (start, start + count - 1)

    def get_seeds(self, n_starts_override: Optional[int] = None) -> List[int]:
        """Return the list of seeds to run."""
        s_start, s_end = self.get_seed_range(n_starts_override=n_starts_override)
        return list(range(s_start, s_end + 1))


@dataclass
class PlottingConfig:
    enabled: bool = True
    n_lowest_chi2: int = 10
    draw_best_fits: bool = True
    unnormalized_moments_file: Optional[str] = None
    reference_amplitude: Optional[str] = None
    anchor_amplitude: Optional[str] = None
    moments: List[str] = field(default_factory=lambda: ["all"])
    formats: List[str] = field(default_factory=lambda: ["png", "pdf"])
    style: str = "default"
    show_phase_ambiguity: bool = True
    posterior_density_bands: bool = True
    distributions: bool = True
    moments_vs_t: bool = True
    moments_vs_t_sharey: bool = False
    amplitudes_vs_t: bool = True
    amplitudes_vs_t_sharey: bool = False
    show_grid: bool = False
    distribution_bins: int = 50
    overview_bins: int = 35
    chi2_landscape_bins: int = 60


@dataclass
class CrossSectionConfig:
    enabled: bool = True
    resonance: str = "a2_1320"
    branching_ratio: float = 0.0565
    luminosity: float = 257.4
    luminosity_unit: str = "pb^-1"
    t_min: Optional[float] = None
    t_max: Optional[float] = None
    delta_t: Optional[float] = None
    mass_limits: Tuple[float, float] = (0.8, 2.0)
    mass_bin_width: Optional[float] = None
    theory_band_file: Optional[str] = None
    phase1_comparison: bool = True
    mcmc_file: Optional[str] = None
    propagate_mcmc: bool = True
    has_kinematic: Optional[bool] = None
    has_barrier: Optional[bool] = None
    barrier_l: int = 2
    barrier_radius: float = 5.0677


@dataclass
class MCMCPropagationConfig:
    enabled: bool = False
    mcmc_tree_path: str = ""
    burn_in: int = 0
    target_entries: int = 200
    stride: Optional[int] = None
    n_workers: int = field(default_factory=get_default_max_workers)
    n_tbin_workers: Optional[int] = None
    output_file: str = "mcmc_propagation_results.root"
    log_dir: str = "logs_mcmc"

    @property
    def max_entries(self) -> int:
        return self.target_entries

    @max_entries.setter
    def max_entries(self, val: int):
        self.target_entries = val


@dataclass
class ProjectConfig:
    name: str = "PartialWaveFit"
    description: str = ""
    tag: str = ""


@dataclass
class FitConfig:
    project: ProjectConfig = field(default_factory=ProjectConfig)
    environment: EnvironmentConfig = field(default_factory=EnvironmentConfig)
    dataset: DatasetConfig = field(default_factory=lambda: DatasetConfig(moments_dir=""))
    chi2: Chi2Config = field(default_factory=Chi2Config)
    waves: List[WaveConfig] = field(default_factory=list)
    minimizer: MinimizerConfig = field(default_factory=MinimizerConfig)
    multi_start: MultiStartConfig = field(default_factory=MultiStartConfig)
    mcmc_propagation: MCMCPropagationConfig = field(default_factory=MCMCPropagationConfig)
    plotting: PlottingConfig = field(default_factory=PlottingConfig)
    cross_section: CrossSectionConfig = field(default_factory=CrossSectionConfig)
    shared_phases_by_reflectivity: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def get_mass_bins(self) -> List[float]:
        return self.dataset.mass_bins.get_bin_centers()

    def get_wave_for_l(self, l: int) -> Optional[WaveConfig]:
        for w in self.waves:
            if w.l == l:
                return w
        return None


def get_two_bw_wave_name(w: WaveConfig) -> str:
    """Derive canonical wave name for a TwoBreitWigner wave (e.g. 'a2_1320_1700')."""
    if w.name and "_" in w.name and not w.name.lower().endswith("wave"):
        return w.name
    if len(w.resonances) >= 2:
        r1, r2 = w.resonances[0].name, w.resonances[1].name
        parts1 = r1.split("_", 1)
        parts2 = r2.split("_", 1)
        if len(parts1) == 2 and len(parts2) == 2 and parts1[0] == parts2[0]:
            return f"{r1}_{parts2[1]}"
        return f"{r1}_{r2}"
    elif len(w.resonances) == 1:
        return w.resonances[0].name
    return f"wave_{w.l}"


def resolve_gaussian_constraint_param(
    par_key: str,
    res_name: str,
    wave_name: Optional[str] = None,
    res_idx: Optional[int] = None,
    gs_key: str = "g_etapi",
    model: Optional[ModelType] = None,
) -> str:
    """Resolve user-friendly constraint keys (e.g. 'mass', 'width') to internal MD parameter names."""
    k_raw = str(par_key).strip()

    if res_idx is not None and wave_name:
        if k_raw == f"MD_shared_{res_name}_Mass":
            return f"MD_shared_{wave_name}_Mass{res_idx}"
        if k_raw == f"MD_shared_{res_name}_Width":
            return f"MD_shared_{wave_name}_Width{res_idx}"
        if k_raw.startswith("MD_") or k_raw.startswith("MI_"):
            return k_raw

        k_low = k_raw.lower()
        alias_map = {
            "mass": f"MD_shared_{wave_name}_Mass{res_idx}",
            "m": f"MD_shared_{wave_name}_Mass{res_idx}",
            f"mass{res_idx}": f"MD_shared_{wave_name}_Mass{res_idx}",
            f"m{res_idx}": f"MD_shared_{wave_name}_Mass{res_idx}",
            "mass1": f"MD_shared_{wave_name}_Mass1",
            "m1": f"MD_shared_{wave_name}_Mass1",
            "mass2": f"MD_shared_{wave_name}_Mass2",
            "m2": f"MD_shared_{wave_name}_Mass2",
            "width": f"MD_shared_{wave_name}_Width{res_idx}",
            "w": f"MD_shared_{wave_name}_Width{res_idx}",
            f"width{res_idx}": f"MD_shared_{wave_name}_Width{res_idx}",
            f"w{res_idx}": f"MD_shared_{wave_name}_Width{res_idx}",
            "width1": f"MD_shared_{wave_name}_Width1",
            "w1": f"MD_shared_{wave_name}_Width1",
            "width2": f"MD_shared_{wave_name}_Width2",
            "w2": f"MD_shared_{wave_name}_Width2",
        }
        return alias_map.get(k_low, k_raw)

    if k_raw.startswith("MD_") or k_raw.startswith("MI_"):
        return k_raw

    k_low = k_raw.lower()
    if model == ModelType.FLATTE_BREIT_WIGNER:
        alias_map = {
            "mass_flatte": f"MD_shared_{res_name}_MassFlatte",
            "massflatte": f"MD_shared_{res_name}_MassFlatte",
            "m_flatte": f"MD_shared_{res_name}_MassFlatte",
            "mass_bw": f"MD_shared_{res_name}_MassBW",
            "massbw": f"MD_shared_{res_name}_MassBW",
            "m_bw": f"MD_shared_{res_name}_MassBW",
            "width_bw": f"MD_shared_{res_name}_WidthBW",
            "widthbw": f"MD_shared_{res_name}_WidthBW",
            "w_bw": f"MD_shared_{res_name}_WidthBW",
            "mass": f"MD_shared_{res_name}_MassFlatte",
            "m": f"MD_shared_{res_name}_MassFlatte",
            "width": f"MD_shared_{res_name}_WidthBW",
            "w": f"MD_shared_{res_name}_WidthBW",
            "g_s": f"MD_shared_{res_name}_{gs_key}",
            "g_etapi": f"MD_shared_{res_name}_{gs_key}",
            "g_kk": f"MD_shared_{res_name}_g_KK",
            "g_etaprimepi": f"MD_shared_{res_name}_g_etaprimepi",
        }
        return alias_map.get(k_low, k_raw)

    alias_map = {
        "mass": f"MD_shared_{res_name}_Mass",
        "m": f"MD_shared_{res_name}_Mass",
        "width": f"MD_shared_{res_name}_Width",
        "w": f"MD_shared_{res_name}_Width",
        "g_s": f"MD_shared_{res_name}_{gs_key}",
        "g_etapi": f"MD_shared_{res_name}_{gs_key}",
        "g_kk": f"MD_shared_{res_name}_g_KK",
        "g_etaprimepi": f"MD_shared_{res_name}_g_etaprimepi",
    }
    return alias_map.get(k_low, k_raw)

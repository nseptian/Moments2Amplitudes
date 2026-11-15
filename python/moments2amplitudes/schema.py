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
    included_l: List[int] = field(default_factory=lambda: [0, 2, 4])
    excluded_l: List[int] = field(default_factory=list)
    max_l: int = 4
    include_all: bool = False
    no_use: List[str] = field(default_factory=lambda: ["H_3"])


@dataclass
class ResonanceConfig:
    name: str
    mass: float = 1.318
    width: float = 0.105
    mass_fixed: bool = False
    width_fixed: bool = False
    mass_limits: Optional[List[float]] = None
    width_limits: Optional[List[float]] = None
    gaussian_constraints: Dict[str, List[float]] = field(default_factory=dict)  # param -> [mean, sigma]
    phase_reference: bool = False  # Fixed to zero as phase reference
    phase_fixed: bool = False
    fixed_phase: Optional[float] = None
    phase_limits: Optional[List[float]] = None


@dataclass
class ParameterInitRule:
    strategy: InitStrategy = InitStrategy.DEFAULT
    min_val: float = 0.0
    max_val: float = 200.0
    value: float = 0.0
    is_fixed: bool = False


@dataclass
class WaveInitConfig:
    strategy: InitStrategy = InitStrategy.DEFAULT
    fixed_zero: bool = False
    file_path: Optional[str] = None
    option: int = 0
    shared_phase: Union[bool, str, int] = False
    magnitude: ParameterInitRule = field(default_factory=lambda: ParameterInitRule(strategy=InitStrategy.RANDOM, min_val=0.0, max_val=200.0))
    phase: ParameterInitRule = field(default_factory=lambda: ParameterInitRule(strategy=InitStrategy.FIXED_VALUE, value=0.0, is_fixed=False))


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
    phase_reference: Union[bool, str] = False
    shared_phase: Optional[Union[bool, str, int]] = None

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


@dataclass
class MinimizerConfig:
    algorithm: str = "Minuit2"
    strategy: int = 1
    tolerance: float = 1.0e-3
    max_calls: int = 1000000
    max_iterations: int = 100000
    print_level: int = 0


@dataclass
class MultiStartConfig:
    enabled: bool = True
    n_starts: int = 100
    n_threads: int = field(default_factory=get_default_max_workers)
    n_tbin_workers: Optional[int] = None
    output_dir: str = "logs"
    plateau_delta_chi2: float = 0.05
    seed_offset: int = 1
    n_lowest_chi2: int = 10


@dataclass
class PlottingConfig:
    enabled: bool = True
    n_lowest_chi2: int = 10
    draw_best_fits: bool = True
    unnormalized_moments_file: Optional[str] = None
    moments: List[str] = field(default_factory=lambda: ["all"])
    formats: List[str] = field(default_factory=lambda: ["png", "pdf"])
    style: str = "default"
    show_phase_ambiguity: bool = True
    posterior_density_bands: bool = True
    distributions: bool = True
    moments_vs_t: bool = True
    moments_vs_t_sharey: bool = False
    show_grid: bool = False


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


@dataclass
class MCMCPropagationConfig:
    enabled: bool = False
    mcmc_tree_path: str = ""
    burn_in: int = 0
    max_entries: int = 200
    stride: int = 1
    n_workers: int = field(default_factory=get_default_max_workers)
    n_tbin_workers: Optional[int] = None
    output_file: str = "mcmc_propagation_results.root"
    log_dir: str = "logs_mcmc"


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
    phase_reference: Optional[Union[str, List[str]]] = None
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

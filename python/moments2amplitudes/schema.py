"""
Schema definitions for Moments2Amplitudes fit configurations.
Uses standard Python dataclasses for lightweight, dependency-free validation.
"""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import List, Dict, Optional, Union, Any
import numpy as np


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
    results_filename: str = "ResultsBruMcmcFixedCovariance.root"
    bin_prefix: str = "Mpi0eta"
    var_bin_prefix: str = "t010030"
    mass_bins: MassBinsConfig = field(default_factory=MassBinsConfig)
    use_covariance: bool = True


@dataclass
class Chi2Config:
    included_l: List[int] = field(default_factory=lambda: [0, 2, 4])
    excluded_l: List[int] = field(default_factory=list)
    max_l: int = 4
    include_all: bool = False


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
    n_threads: int = 8
    output_dir: str = "logs"
    plateau_delta_chi2: float = 0.05
    seed_offset: int = 1


@dataclass
class MCMCPropagationConfig:
    enabled: bool = False
    mcmc_tree_path: str = ""
    burn_in: int = 0
    max_entries: int = 200
    stride: int = 1
    n_workers: int = 8
    output_file: str = "mcmc_propagation_results.root"
    log_dir: str = "logs_mcmc"


@dataclass
class ProjectConfig:
    name: str = "PartialWaveFit"
    description: str = ""


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

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def get_mass_bins(self) -> List[float]:
        return self.dataset.mass_bins.get_bin_centers()

    def get_wave_for_l(self, l: int) -> Optional[WaveConfig]:
        for w in self.waves:
            if w.l == l:
                return w
        return None

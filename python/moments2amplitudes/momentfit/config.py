"""
Configuration dataclasses, loading, and validation for the Moments pipeline.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import yaml


@dataclass
class EnvironmentConfig:
    """Environment paths and binary locations."""
    pixi_manifest: str = "/d/home/septian/Moments2Amplitudes/pixi.toml"
    moments_main_dir: str = "/d/grid17/septian/Moments-dev/Moments.Main"
    mom2amp_dir: str = "/d/home/septian/Moments2Amplitudes"
    etapi0_dir: str = "/d/home/septian/EtaPi0Analysis"
    brufit_dir: str = "/d/grid17/septian/brufit_dev/brufit_negative"
    root_sh: str = "/d/grid17/septian/root_v6-36-08_cuda_install/bin/thisroot.sh"
    python_bin: str = "/d/home/septian/Moments2Amplitudes/.pixi/envs/default/bin/python"


@dataclass
class MassBinningConfig:
    """Mass binning definition."""
    variable: str = "Mpi0eta"
    label: str = "#it{m}_{#it{#eta}#it{#pi}^{0}}"
    unit: str = "GeV/#it{c}^{2}"
    min: float = 0.8
    max: float = 2.0
    n_bins: int = 30
    digits: int = 3

    @property
    def bin_width(self) -> float:
        return (self.max - self.min) / self.n_bins

    def get_bin_range(self, bin_idx: int) -> Tuple[float, float]:
        low = self.min + bin_idx * self.bin_width
        high = self.min + (bin_idx + 1) * self.bin_width
        return (low, high)

    def get_bin_center(self, bin_idx: int) -> float:
        low, high = self.get_bin_range(bin_idx)
        return (low + high) / 2.0

    def get_bin_filter(self, bin_idx: int) -> str:
        low, high = self.get_bin_range(bin_idx)
        return f"({self.variable} >= {low:.6f}) && ({self.variable} < {high:.6f})"

    def get_bin_label(self, bin_idx: int) -> str:
        center = self.get_bin_center(bin_idx)
        return f"{self.variable}_{center:.{self.digits}f}"


@dataclass
class AnalysisConfig:
    """Core analysis settings."""
    subsystem: str = "EtaPi0"
    frame: str = "GJ"  # "GJ" (Gottfried-Jackson) or "HF" (Helicity)
    data_period: str = "merged"
    beam_pol_label: str = "All"
    polarization: str = "Pol"  # column name or float
    tbins: List[str] = field(default_factory=lambda: [
        "t010020", "t020032", "t032050", "t050075", "t075100"
    ])
    max_l: List[int] = field(default_factory=lambda: [4, 6, 8])
    mass_binning: MassBinningConfig = field(default_factory=MassBinningConfig)


@dataclass
class InputDataConfig:
    """Input datasets, tree and column specifications."""
    base_dir: str = "/d/home/septian/Moments2Amplitudes/data/nominal"
    tree_name: str = "kin"
    weight_branch: str = "weightASBS"
    data_file_pattern: str = "{data_dir}/{tbin}/data_flat_All.root"
    ps_acc_file_pattern: str = "{data_dir}/{tbin}/phaseSpace_acc_flat_All.root"
    ps_gen_file_pattern: str = "{data_dir}/{tbin}/phaseSpace_gen_flat_All.root"
    columns: Dict[str, str] = field(default_factory=lambda: {
        "theta": "std::acos(cosTheta_eta_gj)",
        "phi": "phi_eta_gj",
        "Phi": "Phi",
    })
    additional_cuts: List[str] = field(default_factory=list)

    def resolve_path(self, pattern: str, tbin: str) -> str:
        return pattern.format(data_dir=self.base_dir, tbin=tbin)


@dataclass
class LinearAlgebraConfig:
    """Stage 1: Linear algebra moments calculation settings."""
    output_base_dir: str = "/d/home/septian/Moments2Amplitudes/Moments.out/nominal"
    method: str = "LIN_ALG_BG_SUBTR_NEG_WEIGHTS"
    weight_branch: str = "weightASBS"
    normalize_moments: bool = False
    bootstrap_samples: int = 0
    force_integral_matrix: bool = True
    calc_acc_ps_moments: bool = False
    openmp_threads: int = 4
    max_parallel_workers: int = 4

    def get_output_dir(self, tbin: str, max_l: int, beam_pol_label: str = "All") -> str:
        return f"{self.output_base_dir}/{tbin}/{beam_pol_label}.maxL_{max_l}"


@dataclass
class ExportCovarianceConfig:
    """Stage 2: Covariance and initial moment value export settings."""
    scale_physics_moments: float = 2.0
    min_eigen_floor: float = 1.0e-4
    yield_variance_floor: float = 10000.0
    add_cholesky_jitter: bool = True
    macro_path: str = "/d/home/septian/Moments2Amplitudes/macros/ExportCovarianceC.C"
    plot_covariance: bool = False


@dataclass
class MCMCConfig:
    """MCMC parameters for MomentFit."""
    init_scale: float = 1.0
    step: float = 1.0
    iterations: int = 200000
    burnin: int = 2000
    apply_mc_variance: bool = False
    allow_negative_pdf: bool = True
    use_tuning: bool = True
    gibbs_block_size: int = 0
    gibbs_block_size_by_max_l: Dict[int, int] = field(default_factory=dict)
    cold_warmup_acc: float = 0.15
    min_acceptance: float = 0.05
    max_acceptance: float = 0.15
    target_acceptance: Optional[float] = None
    tune_mode: str = "MappedRhat"
    rhat_window: int = 1000
    rhat_target: float = 1.05
    ess_fraction: float = 0.02
    rhat_max_retries: int = 5
    bailout_threshold: int = 10000
    bailout_min_acceptance: float = 0.0001
    mc_integration_bins: int = 50
    fixed_covariance: bool = True
    use_covariance: bool = True
    set_init_value: bool = True
    use_gj_frame: bool = True
    num_cpu: int = 1
    seed: int = 111

    def get_gibbs_block_size(self, max_l: int) -> int:
        """Resolve Gibbs block size for a given max_l, falling back to global gibbs_block_size."""
        if self.gibbs_block_size_by_max_l and max_l in self.gibbs_block_size_by_max_l:
            return int(self.gibbs_block_size_by_max_l[max_l])
        return int(self.gibbs_block_size)

    def get_target_acceptance(self) -> float:
        """Resolve target acceptance rate, defaulting to midpoint of [min_acceptance, max_acceptance]."""
        if self.target_acceptance is not None:
            return float(self.target_acceptance)
        return 0.5 * (float(self.min_acceptance) + float(self.max_acceptance))


@dataclass
class DiagnosticsConfig:
    """Diagnostics and auto-retry settings for MomentFit."""
    ess_threshold: float = 20.0
    ess_metric: str = "mean"  # "mean" or "min"
    min_ess_threshold: Optional[float] = 10.0  # Safety floor on worst-moment ESS
    auto_retry: bool = False
    max_retries: int = 2
    max_retry_iters: int = 200000
    burnin_fraction: float = 0.10
    low_acceptance_threshold: float = 0.03
    gibbs_low_acceptance_threshold: float = 0.12
    high_acceptance_threshold: float = 0.35
    gibbs_activation_acceptance: float = 0.05
    recondition_covariance: bool = False
    recondition_shrinkage: float = 0.2
    recondition_min_eigen_floor: float = 1.0e-4


@dataclass
class BruFitDataConfig:
    """Stage 3: BruFit MCMC input data paths configuration."""
    base_dir: str = "/d/home/septian/Moments2Amplitudes/data/data/nominal_tbins"
    data: str = "{data_dir}/data_merged_{tbin}.root"
    mcrecon: str = "{data_dir}/reconmc_merged_{tbin}.root"
    mcgen: str = "{data_dir}/genmc_merged_{tbin}.root"
    weight: str = "{data_dir}/weights_data_merged_{tbin}.root"
    tree_name: str = "kin"
    weight_species: str = "Acc"
    weight_name: str = "weight"
    use_generated_mc: bool = True

    def resolve_path(self, pattern: str, tbin: str) -> str:
        """Resolve a file path pattern for a given t-bin, expanding placeholders."""
        path = pattern.replace("{data_dir}", self.base_dir).replace("{base_dir}", self.base_dir)
        path = path.replace("{tbin}", tbin)
        if not os.path.isabs(path):
            path = os.path.join(self.base_dir, path)
        return os.path.normpath(path)

    def get_data_path(self, tbin: str) -> str:
        return self.resolve_path(self.data, tbin)

    def get_mcrecon_path(self, tbin: str) -> str:
        return self.resolve_path(self.mcrecon, tbin)

    def get_mcgen_path(self, tbin: str) -> str:
        return self.resolve_path(self.mcgen, tbin)

    def get_weight_path(self, tbin: str) -> str:
        return self.resolve_path(self.weight, tbin)

    # Backward-compatible property aliases
    @property
    def data_file_pattern(self) -> str:
        return self.data

    @property
    def mcrecon_file_pattern(self) -> str:
        return self.mcrecon

    @property
    def mcgen_file_pattern(self) -> str:
        return self.mcgen

    @property
    def weight_file_pattern(self) -> str:
        return self.weight


# Backward-compatible alias
Stage3DataConfig = BruFitDataConfig


@dataclass
class BruFitStageConfig:
    """Stage 3: MomentFit with BruFit MCMC settings."""
    macro_path: str = "/d/home/septian/Moments2Amplitudes/macros/MomentFit.C"
    run_script: str = "/d/home/septian/Moments2Amplitudes/macros/RunMomentFit.sh"
    monitor_script: str = "/d/home/septian/Moments2Amplitudes/macros/MonitorMomentFit.py"
    collector_script: str = "/d/home/septian/Moments2Amplitudes/macros/MomentFitCollector.py"
    acceptance_macro: str = "/d/home/septian/Moments2Amplitudes/macros/CalculateAcceptance_m080200.C"
    skip_acceptance: bool = False
    force_acceptance: bool = False
    output_dir_template: str = (
        "/d/home/septian/Moments2Amplitudes/run/"
        "MomentFit_m{tag_min:03d}{tag_max:03d}_FixedCov_Lmax{maxL}_InitValueScale{init_scale:.1f}_step{step:.1f}/"
    )
    data: BruFitDataConfig = field(default_factory=BruFitDataConfig)
    mcmc: MCMCConfig = field(default_factory=MCMCConfig)
    diagnostics: DiagnosticsConfig = field(default_factory=DiagnosticsConfig)
    max_parallel: int = 2
    overwrite: bool = False

    def get_run_dir(self, max_l: int, m_min: float, m_max: float) -> str:
        tag_min = int(round(m_min * 100))
        tag_max = int(round(m_max * 100))
        return self.output_dir_template.format(
            tag_min=tag_min,
            tag_max=tag_max,
            maxL=max_l,
            init_scale=self.mcmc.init_scale,
            step=self.mcmc.step,
        )


@dataclass
class PlottingConfig:
    """Plotting settings for unnormalized moments vs mass."""
    enabled: bool = True
    plot_linear_algebra: bool = True  # Plot linear algebra unnormalized moments vs mass
    plot_mcmc: bool = True            # Plot BruFit MCMC unnormalized moments vs mass
    plot_overlay: bool = True         # Overlay linear algebra and MCMC on top of each other
    moments: List[str] = field(default_factory=lambda: ["all"])
    output_dir: str = "/d/home/septian/Moments2Amplitudes/plots/momentfit"
    formats: List[str] = field(default_factory=lambda: ["pdf", "png"])
    grid: bool = True                 # Multi-panel overview grid plots
    style: str = "default"

    def get_output_dir(self, tbin: str, max_l: int) -> str:
        return f"{self.output_dir}/{tbin}/maxL_{max_l}"


# Backward-compatible aliases
BruFitConfig = BruFitStageConfig
Stage3MomentFitConfig = BruFitStageConfig


@dataclass
class MomentFitConfig:
    """Top-level MomentFit workflow configuration."""
    config_path: str = ""
    stages: List[str] = field(default_factory=lambda: [
        "linear_algebra",
        "export_covariance",
        "moment_fit",
    ])
    environment: EnvironmentConfig = field(default_factory=EnvironmentConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    input_data: InputDataConfig = field(default_factory=InputDataConfig)
    linear_algebra: LinearAlgebraConfig = field(default_factory=LinearAlgebraConfig)
    export_covariance: ExportCovarianceConfig = field(default_factory=ExportCovarianceConfig)
    moment_fit: BruFitStageConfig = field(default_factory=BruFitStageConfig)
    plotting: PlottingConfig = field(default_factory=PlottingConfig)


# Backward-compatible aliases
PipelineConfig = MomentFitConfig
MomentFitWorkflowConfig = MomentFitConfig


def _parse_brufit_data_config(data: Dict[str, Any]) -> BruFitDataConfig:
    """Parse BruFit data config dictionary handling aliases."""
    if not isinstance(data, dict):
        return BruFitDataConfig()
    kwargs = {}
    key_mapping = {
        "base_dir": ["base_dir", "data_dir", "dir", "input_data_dir"],
        "data": ["data", "data_file", "data_file_pattern", "data_path"],
        "mcrecon": ["mcrecon", "mcrecon_file", "mcrecon_file_pattern", "reconmc", "reconmc_file", "simulated", "simulated_file"],
        "mcgen": ["mcgen", "mcgen_file", "mcgen_file_pattern", "genmc", "genmc_file", "generated", "generated_file"],
        "weight": ["weight", "weight_file", "weight_file_pattern", "weights", "weights_file", "weights_file_pattern"],
        "tree_name": ["tree_name", "treename", "tree"],
        "weight_species": ["weight_species", "species"],
        "weight_name": ["weight_name", "weight_branch", "branch_name"],
        "use_generated_mc": ["use_generated_mc", "acceptance_correction"],
    }
    for target_key, aliases in key_mapping.items():
        for alias in aliases:
            if alias in data:
                kwargs[target_key] = data[alias]
                break
    return BruFitDataConfig(**kwargs)


def _parse_dict_to_dataclass(cls: Any, data: Dict[str, Any]) -> Any:
    """Recursively parse dictionary into dataclass instances."""
    if not isinstance(data, dict):
        return data
    kwargs = {}
    for f in cls.__dataclass_fields__.values():
        name = f.name
        if name in data:
            val = data[name]
            field_type = f.type
            if name == "gibbs_block_size_by_max_l":
                if val is None:
                    kwargs[name] = {}
                elif isinstance(val, dict):
                    kwargs[name] = {int(k): int(v) for k, v in val.items()}
                else:
                    raise ValueError(f"Expected dict for '{name}', got {type(val).__name__}")
            # Handle nested dataclasses
            elif hasattr(field_type, "__dataclass_fields__"):
                kwargs[name] = _parse_dict_to_dataclass(field_type, val)
            elif isinstance(val, dict) and hasattr(f.default_factory, "__self__"):
                # fallback for default factory classes
                kwargs[name] = val
            else:
                kwargs[name] = val
    return cls(**kwargs)


def load_momentfit_config(config_path: Union[str, Path]) -> MomentFitConfig:
    """Load and parse MomentFit YAML configuration."""
    path = Path(config_path)
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}

    cfg = MomentFitConfig()
    cfg.config_path = str(path.resolve())

    if "stages" in data:
        cfg.stages = list(data["stages"])
    if "environment" in data:
        cfg.environment = _parse_dict_to_dataclass(EnvironmentConfig, data["environment"])
    if "analysis" in data:
        analysis_data = dict(data["analysis"])
        if "mass_binning" in analysis_data:
            mass_bin_data = analysis_data.pop("mass_binning")
            cfg.analysis.mass_binning = _parse_dict_to_dataclass(MassBinningConfig, mass_bin_data)
        for k, v in analysis_data.items():
            if hasattr(cfg.analysis, k):
                setattr(cfg.analysis, k, v)
    if "input_data" in data:
        cfg.input_data = _parse_dict_to_dataclass(InputDataConfig, data["input_data"])
    if "linear_algebra" in data:
        cfg.linear_algebra = _parse_dict_to_dataclass(LinearAlgebraConfig, data["linear_algebra"])
    if "export_covariance" in data:
        cfg.export_covariance = _parse_dict_to_dataclass(ExportCovarianceConfig, data["export_covariance"])
    if "moment_fit" in data:
        mf_data = dict(data["moment_fit"])
        if "data" in mf_data and isinstance(mf_data["data"], dict):
            cfg.moment_fit.data = _parse_brufit_data_config(mf_data.pop("data"))
        elif "input_data" in mf_data and isinstance(mf_data["input_data"], dict):
            cfg.moment_fit.data = _parse_brufit_data_config(mf_data.pop("input_data"))
        elif "paths" in mf_data and isinstance(mf_data["paths"], dict):
            cfg.moment_fit.data = _parse_brufit_data_config(mf_data.pop("paths"))
        elif any(k in mf_data for k in ("mcrecon", "mcgen", "weight")):
            direct_data = {}
            for k in ("base_dir", "data_dir", "data", "mcrecon", "mcgen", "weight", "tree_name", "weight_species", "weight_name", "use_generated_mc"):
                if k in mf_data:
                    direct_data[k] = mf_data.pop(k)
            cfg.moment_fit.data = _parse_brufit_data_config(direct_data)
        if "mcmc" in mf_data:
            cfg.moment_fit.mcmc = _parse_dict_to_dataclass(MCMCConfig, mf_data.pop("mcmc"))
        if "diagnostics" in mf_data:
            cfg.moment_fit.diagnostics = _parse_dict_to_dataclass(DiagnosticsConfig, mf_data.pop("diagnostics"))
        for k, v in mf_data.items():
            if hasattr(cfg.moment_fit, k):
                setattr(cfg.moment_fit, k, v)
    if "plotting" in data:
        cfg.plotting = _parse_dict_to_dataclass(PlottingConfig, data["plotting"])

    return cfg


def validate_momentfit_config(cfg: MomentFitConfig, check_files: bool = True) -> List[str]:
    """Validate MomentFit configuration and check existence of critical directories."""
    errors: List[str] = []

    if check_files:
        # Check Moments.Main directory
        if not os.path.isdir(cfg.environment.moments_main_dir):
            errors.append(f"Moments.Main directory does not exist: {cfg.environment.moments_main_dir}")

        # Check basisFunctions.C in Moments.Main
        bf_path = os.path.join(cfg.environment.moments_main_dir, "cpp", "basisFunctions.C")
        if not os.path.isfile(bf_path):
            errors.append(f"basisFunctions.C not found in Moments.Main: {bf_path}")

        # Check data directory
        if not os.path.isdir(cfg.input_data.base_dir):
            errors.append(f"Input data directory does not exist: {cfg.input_data.base_dir}")

    # Check t-bins
    if not cfg.analysis.tbins:
        errors.append("No t-bins defined in analysis configuration.")
    if check_files:
        for tbin in cfg.analysis.tbins:
            data_f = cfg.input_data.resolve_path(cfg.input_data.data_file_pattern, tbin)
            if not os.path.isfile(data_f):
                errors.append(f"Real data file not found for t-bin '{tbin}': {data_f}")

    # Check maxL
    if not cfg.analysis.max_l:
        errors.append("No max_l values defined in analysis configuration.")
    for l in cfg.analysis.max_l:
        if l <= 0 or l % 2 != 0:
            errors.append(f"max_l must be positive even integer, got {l}")

    # Check export covariance macro
    if check_files and not os.path.isfile(cfg.export_covariance.macro_path):
        errors.append(f"ExportCovarianceC.C macro not found: {cfg.export_covariance.macro_path}")

    # Check MomentFit script and macro
    if "moment_fit" in cfg.stages:
        if check_files and not os.path.isfile(cfg.moment_fit.run_script):
            errors.append(f"MomentFit run script not found: {cfg.moment_fit.run_script}")
        if check_files and not os.path.isfile(cfg.moment_fit.macro_path):
            errors.append(f"MomentFit.C macro not found: {cfg.moment_fit.macro_path}")
        if check_files and cfg.moment_fit.acceptance_macro and not os.path.isfile(cfg.moment_fit.acceptance_macro):
            errors.append(f"Acceptance macro not found: {cfg.moment_fit.acceptance_macro}")

        # Check MomentFit data base_dir if configured
        if check_files and cfg.moment_fit.data.base_dir:
            base = cfg.moment_fit.data.base_dir
            if not os.path.exists(base) and not any(
                os.path.islink(p) for p in [base, os.path.dirname(base), os.path.dirname(os.path.dirname(base))]
            ):
                errors.append(f"MomentFit data base_dir does not exist: {base}")

        # Check MomentFit MCMC parameter consistency
        mcmc_cfg = cfg.moment_fit.mcmc
        if not (0.0 < mcmc_cfg.min_acceptance < mcmc_cfg.max_acceptance <= 1.0):
            errors.append(
                f"mcmc.min_acceptance ({mcmc_cfg.min_acceptance}) and mcmc.max_acceptance ({mcmc_cfg.max_acceptance}) must satisfy 0 < min < max <= 1."
            )
        if mcmc_cfg.target_acceptance is not None and not (mcmc_cfg.min_acceptance <= mcmc_cfg.target_acceptance <= mcmc_cfg.max_acceptance):
            errors.append(
                f"mcmc.target_acceptance ({mcmc_cfg.target_acceptance}) must lie within [{mcmc_cfg.min_acceptance}, {mcmc_cfg.max_acceptance}]."
            )
        if mcmc_cfg.tune_mode not in ("MappedRhat", "Acceptance", "kMappedRhat", "kAcceptance"):
            errors.append(f"mcmc.tune_mode must be 'MappedRhat' or 'Acceptance', got '{mcmc_cfg.tune_mode}'.")
        if mcmc_cfg.gibbs_block_size < 0 or any(v < 0 for v in mcmc_cfg.gibbs_block_size_by_max_l.values()):
            errors.append("mcmc.gibbs_block_size and all mcmc.gibbs_block_size_by_max_l values must be >= 0.")
        if mcmc_cfg.iterations <= 0 or mcmc_cfg.burnin < 0 or mcmc_cfg.burnin >= mcmc_cfg.iterations:
            errors.append(f"mcmc.burnin ({mcmc_cfg.burnin}) must be >= 0 and < mcmc.iterations ({mcmc_cfg.iterations}).")
        if mcmc_cfg.mc_integration_bins <= 0:
            errors.append(f"mcmc.mc_integration_bins must be > 0, got {mcmc_cfg.mc_integration_bins}.")
        if mcmc_cfg.num_cpu < 1:
            errors.append(f"mcmc.num_cpu must be >= 1, got {mcmc_cfg.num_cpu}.")
        if mcmc_cfg.bailout_threshold < 0 or not (0.0 <= mcmc_cfg.bailout_min_acceptance <= 1.0):
            errors.append("mcmc.bailout_threshold must be >= 0 and bailout_min_acceptance in [0, 1].")

        # Check MomentFit diagnostics thresholds
        diag_cfg = cfg.moment_fit.diagnostics
        if diag_cfg.ess_metric not in ("mean", "min"):
            errors.append(f"diagnostics.ess_metric must be 'mean' or 'min', got '{diag_cfg.ess_metric}'.")
        if diag_cfg.max_retry_iters <= 0:
            errors.append(f"diagnostics.max_retry_iters must be > 0, got {diag_cfg.max_retry_iters}.")
        if not (0.0 <= diag_cfg.burnin_fraction < 0.5):
            errors.append(f"diagnostics.burnin_fraction must be in [0.0, 0.5), got {diag_cfg.burnin_fraction}.")
        if not (0.0 < diag_cfg.low_acceptance_threshold < diag_cfg.high_acceptance_threshold <= 1.0):
            errors.append("diagnostics.low_acceptance_threshold must be in (0, high_acceptance_threshold) and <= 1.")
        if not (0.0 < diag_cfg.gibbs_low_acceptance_threshold < diag_cfg.high_acceptance_threshold <= 1.0):
            errors.append("diagnostics.gibbs_low_acceptance_threshold must be in (0, high_acceptance_threshold) and <= 1.")
        if not (0.0 <= diag_cfg.gibbs_activation_acceptance <= 1.0):
            errors.append(f"diagnostics.gibbs_activation_acceptance must be in [0, 1], got {diag_cfg.gibbs_activation_acceptance}.")
        if not (0.0 <= diag_cfg.recondition_shrinkage <= 1.0):
            errors.append(f"diagnostics.recondition_shrinkage must be in [0, 1], got {diag_cfg.recondition_shrinkage}.")
        if diag_cfg.recondition_min_eigen_floor <= 0.0:
            errors.append(f"diagnostics.recondition_min_eigen_floor must be > 0, got {diag_cfg.recondition_min_eigen_floor}.")

    return errors


# Backward-compatible function aliases
load_pipeline_config = load_momentfit_config
validate_pipeline_config = validate_momentfit_config

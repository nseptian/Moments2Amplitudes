"""
Configuration loader, parser, and validator for Moments2Amplitudes.
Supports YAML files, environment variable substitution (${VAR} and ${VAR:-default}),
and semantic validation.
"""

from __future__ import annotations
import os
import re
import yaml
import json
from pathlib import Path
from typing import Dict, Any, List, Optional

from .schema import (
    FitConfig,
    ProjectConfig,
    EnvironmentConfig,
    DatasetConfig,
    MassBinsConfig,
    Chi2Config,
    WaveConfig,
    WaveType,
    ModelType,
    InitStrategy,
    WaveInitConfig,
    ParameterInitRule,
    ResonanceConfig,
    MinimizerConfig,
    MultiStartConfig,
    MCMCPropagationConfig,
)


def expand_env_vars(val: Any) -> Any:
    """Recursively expand environment variables like ${VAR} or ${VAR:-default} in strings."""
    if isinstance(val, str):
        # Match pattern ${VAR:-default} or ${VAR}
        pattern = re.compile(r"\$\{([A-Za-z0-9_]+)(?::-([^}]*))?\}")

        def repl(match):
            var_name = match.group(1)
            default_val = match.group(2) if match.group(2) is not None else ""
            return os.environ.get(var_name, default_val)

        res = pattern.sub(repl, val)
        # Also expand standard $VAR if not caught
        res = os.path.expandvars(res)
        return res
    elif isinstance(val, dict):
        return {k: expand_env_vars(v) for k, v in val.items()}
    elif isinstance(val, list):
        return [expand_env_vars(item) for item in val]
    return val


def parse_wave_config(wave_dict: Dict[str, Any]) -> WaveConfig:
    l_val = wave_dict.get("l", 0)
    name = wave_dict.get("name", "")
    type_str = wave_dict.get("type", "mass_dependent").lower()
    wave_type = WaveType.PIECEWISE if "piecewise" in type_str or "mass_independent" in type_str else WaveType.MASS_DEPENDENT

    model_str = wave_dict.get("model", "breit_wigner").lower()
    if "custom" in model_str:
        model_type = ModelType.CUSTOM
    elif "two" in model_str:
        model_type = ModelType.TWO_BREIT_WIGNER
    elif "flatte" in model_str and "breit" in model_str:
        model_type = ModelType.FLATTE_BREIT_WIGNER
    elif "flatte" in model_str:
        model_type = ModelType.FLATTE
    elif "poly" in model_str:
        model_type = ModelType.POLYNOMIAL
    else:
        model_type = ModelType.BREIT_WIGNER

    formula = wave_dict.get("formula")
    parameters = wave_dict.get("parameters", [])

    # Assign default explicit formulas if not specified
    if not formula:
        if model_type == ModelType.BREIT_WIGNER:
            formula = "k / (M^2 - s - i * M * width)"
            if not parameters:
                parameters = ["k", "Mass", "Width"]
        elif model_type == ModelType.TWO_BREIT_WIGNER:
            formula = "k1 / (M1^2 - s - i * M1 * w1) + (k2 / (M2^2 - s - i * M2 * w2)) * exp(i * phase)"
            if not parameters:
                parameters = ["k1", "Mass1", "Width1", "k2", "Mass2", "Width2"]
        elif model_type == ModelType.FLATTE:
            formula = "k / (M0^2 - s - i * (g_etapi^2 * rho_etapi + g_KK^2 * rho_KK + g_etaprimepi^2 * rho_etaprimepi))"
            if not parameters:
                parameters = ["k", "Mass", "g_etapi", "g_KK", "g_etaprimepi"]
        elif model_type == ModelType.FLATTE_BREIT_WIGNER:
            formula = "k_fl / (M0^2 - s - i * (g_etapi^2 * rho_etapi + g_KK^2 * rho_KK)) + k_bw / (M_bw^2 - s - i * M_bw * width_bw)"
            if not parameters:
                parameters = ["k_fl", "MassFlatte", "g_etapi", "g_KK", "g_etaprimepi", "k_bw", "MassBW", "WidthBW"]
        elif model_type == ModelType.POLYNOMIAL:
            formula = "sum_{n=0}^{order} (re_n + i*im_n) * z^n"
    elif not parameters:
        # If formula was provided without parameters, use standard defaults if applicable
        if model_type == ModelType.BREIT_WIGNER:
            parameters = ["k", "Mass", "Width"]
        elif model_type == ModelType.FLATTE:
            parameters = ["k", "Mass", "g_etapi", "g_KK", "g_etaprimepi"]

    poly_order = wave_dict.get("poly_order", -1)

    # Resonances
    resonances: List[ResonanceConfig] = []
    for res_dict in wave_dict.get("resonances", []):
        resonances.append(
            ResonanceConfig(
                name=res_dict.get("name", "resonance"),
                mass=float(res_dict.get("mass", 1.318)),
                width=float(res_dict.get("width", 0.105)),
                mass_fixed=bool(res_dict.get("mass_fixed", False)),
                width_fixed=bool(res_dict.get("width_fixed", False)),
                mass_limits=res_dict.get("mass_limits"),
                width_limits=res_dict.get("width_limits"),
                gaussian_constraints=res_dict.get("gaussian_constraints", {}),
            )
        )

    # Init config
    init_dict = wave_dict.get("init", {})
    strat_str = init_dict.get("strategy", "default").lower()
    fixed_zero = bool(init_dict.get("fixed_zero", False))
    if fixed_zero:
        init_strat = InitStrategy.FIXED_ZERO
    elif strat_str == "from_file":
        init_strat = InitStrategy.FROM_FILE
    elif strat_str == "fixed_zero":
        init_strat = InitStrategy.FIXED_ZERO
    elif strat_str == "fixed":
        init_strat = InitStrategy.FIXED_VALUE
    elif strat_str == "random":
        init_strat = InitStrategy.RANDOM
    else:
        init_strat = InitStrategy.DEFAULT

    mag_dict = init_dict.get("magnitude", {})
    mag_rule = ParameterInitRule(
        strategy=InitStrategy.RANDOM if mag_dict.get("strategy", "random") == "random" else InitStrategy.FIXED_VALUE,
        min_val=float(mag_dict.get("min", 0.0)),
        max_val=float(mag_dict.get("max", 200.0)),
        value=float(mag_dict.get("value", 0.0)),
        is_fixed=bool(mag_dict.get("fixed", False)),
    )

    phase_dict = init_dict.get("phase", {})
    phase_strat = phase_dict.get("strategy", "fixed").lower()
    phase_rule = ParameterInitRule(
        strategy=InitStrategy.FIXED_VALUE if phase_strat == "fixed" else InitStrategy.RANDOM,
        min_val=float(phase_dict.get("min", -3.14159265)),
        max_val=float(phase_dict.get("max", 3.14159265)),
        value=float(phase_dict.get("value", 0.0)),
        is_fixed=bool(phase_dict.get("fixed", False)),
    )

    init_cfg = WaveInitConfig(
        strategy=init_strat,
        fixed_zero=fixed_zero,
        file_path=init_dict.get("file_path"),
        option=int(init_dict.get("option", 0)),
        magnitude=mag_rule,
        phase=phase_rule,
    )

    return WaveConfig(
        l=l_val,
        name=name,
        type=wave_type,
        model=model_type,
        formula=formula,
        parameters=parameters,
        poly_order=poly_order,
        resonances=resonances,
        init=init_cfg,
    )


def load_config(file_path: Union[str, Path]) -> FitConfig:
    """Load and parse a YAML or JSON configuration file with environment variable expansion."""
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {file_path}")

    with open(path, "r") as f:
        raw_data = yaml.safe_load(f)

    if not isinstance(raw_data, dict):
        raise ValueError(f"Invalid configuration file: root must be a mapping, got {type(raw_data)}")

    data = expand_env_vars(raw_data)

    # Project
    proj_dict = data.get("project", {})
    project = ProjectConfig(
        name=proj_dict.get("name", path.stem),
        description=proj_dict.get("description", ""),
    )

    # Environment
    env_dict = data.get("environment", {})
    brufit_dir = env_dict.get("brufit_dir") or os.environ.get("BRUFIT", "/d/grid17/septian/brufit_dev/brufit")
    root_dir = env_dict.get("root_dir") or os.environ.get("ROOTSYS")
    environment = EnvironmentConfig(
        brufit_dir=brufit_dir,
        root_dir=root_dir,
        extra_root_includes=env_dict.get("extra_root_includes", []),
    )

    # Dataset
    ds_dict = data.get("dataset", {})
    mb_dict = ds_dict.get("mass_bins", {})
    mass_bins = MassBinsConfig(
        min=float(mb_dict["min"]) if "min" in mb_dict and mb_dict["min"] is not None else None,
        max=float(mb_dict["max"]) if "max" in mb_dict and mb_dict["max"] is not None else None,
        step=float(mb_dict["step"]) if "step" in mb_dict and mb_dict["step"] is not None else None,
        bins=[float(b) for b in mb_dict["bins"]] if "bins" in mb_dict and mb_dict["bins"] is not None else None,
    )

    dataset = DatasetConfig(
        moments_dir=ds_dict.get("moments_dir", ""),
        results_filename=ds_dict.get("results_filename", "ResultsBruMcmcFixedCovariance.root"),
        bin_prefix=ds_dict.get("bin_prefix", "Mpi0eta"),
        var_bin_prefix=ds_dict.get("var_bin_prefix", "t010030"),
        mass_bins=mass_bins,
        use_covariance=bool(ds_dict.get("use_covariance", True)),
    )

    # Chi2
    chi2_dict = data.get("chi2", {})
    chi2 = Chi2Config(
        included_l=[int(x) for x in chi2_dict.get("included_l", [0, 2, 4])],
        excluded_l=[int(x) for x in chi2_dict.get("excluded_l", [])],
        max_l=int(chi2_dict.get("max_l", 4)),
        include_all=bool(chi2_dict.get("include_all", False)),
    )

    # Waves
    waves: List[WaveConfig] = []
    for w in data.get("waves", []):
        waves.append(parse_wave_config(w))

    # Minimizer
    min_dict = data.get("minimizer", {})
    minimizer = MinimizerConfig(
        algorithm=min_dict.get("algorithm", "Minuit2"),
        strategy=int(min_dict.get("strategy", 1)),
        tolerance=float(min_dict.get("tolerance", 1.0e-3)),
        max_calls=int(min_dict.get("max_calls", 1000000)),
        max_iterations=int(min_dict.get("max_iterations", 100000)),
        print_level=int(min_dict.get("print_level", 0)),
    )

    # MultiStart
    ms_dict = data.get("multi_start", {})
    multi_start = MultiStartConfig(
        enabled=bool(ms_dict.get("enabled", True)),
        n_starts=int(ms_dict.get("n_starts", 100)),
        n_threads=int(ms_dict.get("n_threads", 8)),
        output_dir=ms_dict.get("output_dir", "logs"),
        plateau_delta_chi2=float(ms_dict.get("plateau_delta_chi2", 0.05)),
        seed_offset=int(ms_dict.get("seed_offset", 1)),
    )

    # MCMCPropagation
    mcmc_dict = data.get("mcmc_propagation", {})
    mcmc_propagation = MCMCPropagationConfig(
        enabled=bool(mcmc_dict.get("enabled", False)),
        mcmc_tree_path=mcmc_dict.get("mcmc_tree_path", ""),
        burn_in=int(mcmc_dict.get("burn_in", 0)),
        max_entries=int(mcmc_dict.get("max_entries", 200)),
        stride=int(mcmc_dict.get("stride", 1)),
        n_workers=int(mcmc_dict.get("n_workers", 8)),
        output_file=mcmc_dict.get("output_file", "mcmc_propagation_results.root"),
        log_dir=mcmc_dict.get("log_dir", "logs_mcmc"),
    )

    return FitConfig(
        project=project,
        environment=environment,
        dataset=dataset,
        chi2=chi2,
        waves=waves,
        minimizer=minimizer,
        multi_start=multi_start,
        mcmc_propagation=mcmc_propagation,
    )


def save_config(config: FitConfig, file_path: Union[str, Path]) -> None:
    """Save FitConfig to YAML or JSON."""
    path = Path(file_path)
    data = config.to_dict()
    with open(path, "w") as f:
        if path.suffix.lower() == ".json":
            json.dump(data, f, indent=2)
        else:
            yaml.dump(data, f, sort_keys=False, default_flow_style=False)


def validate_config(config: FitConfig) -> List[str]:
    """
    Validate configuration and return a list of warning or error messages.
    Empty list indicates configuration is valid.
    """
    errors: List[str] = []

    # Dataset validation
    if not config.dataset.moments_dir:
        errors.append("dataset.moments_dir must not be empty.")

    try:
        bins = config.get_mass_bins()
        if len(bins) == 0:
            errors.append("dataset.mass_bins resolved to 0 mass bins.")
    except Exception as e:
        errors.append(f"Failed to generate mass bins: {e}")

    # Waves validation
    if len(config.waves) == 0:
        errors.append("No waves defined in configuration.")

    l_set = set()
    for wave in config.waves:
        if wave.l in l_set:
            errors.append(f"Duplicate wave configuration for L={wave.l}.")
        l_set.add(wave.l)

        if wave.type == WaveType.MASS_DEPENDENT:
            if wave.model in (ModelType.BREIT_WIGNER, ModelType.TWO_BREIT_WIGNER, ModelType.FLATTE_BREIT_WIGNER):
                if len(wave.resonances) == 0:
                    errors.append(f"Wave L={wave.l} uses model {wave.model.value} but defines no resonances.")

    # MCMC validation
    if config.mcmc_propagation.enabled:
        if not config.mcmc_propagation.mcmc_tree_path:
            errors.append("mcmc_propagation is enabled but mcmc_tree_path is not specified.")

    return errors


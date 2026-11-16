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
    PlottingConfig,
    CrossSectionConfig,
)
from .environment import get_default_max_workers


TBIN_INDEX_MAP = {
    0: "t010020",
    1: "t020032",
    2: "t032050",
    3: "t050075",
    4: "t075100",
}
STANDARD_TBINS = ["t010020", "t020032", "t032050", "t050075", "t075100"]


def parse_t_bins_input(val: Any) -> List[str]:
    """
    Parse t-bins from various user input forms:
    - "all" / ["all"] -> all 5 standard t-bins
    - [0, 1, 2] -> ["t010020", "t020032", "t032050"]
    - "t010020,t020032" or ["t010020", "t020032"]
    """
    if val is None:
        return []

    if isinstance(val, int):
        return [TBIN_INDEX_MAP.get(val, f"t{val}")]

    if isinstance(val, str):
        val_str = val.strip()
        if not val_str:
            return []
        if val_str.lower() == "all":
            return list(STANDARD_TBINS)
        parts = [p.strip() for p in val_str.replace(",", " ").split() if p.strip()]
        out: List[str] = []
        for p in parts:
            if p.lower() == "all":
                return list(STANDARD_TBINS)
            if p.isdigit() and int(p) in TBIN_INDEX_MAP:
                out.append(TBIN_INDEX_MAP[int(p)])
            else:
                out.append(p)
        return out

    if isinstance(val, (list, tuple)):
        out_list: List[str] = []
        for item in val:
            sub = parse_t_bins_input(item)
            for s in sub:
                if s not in out_list:
                    out_list.append(s)
        return out_list

    return [str(val).strip()]


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
        phase_raw = res_dict.get("phase", None)
        phase_val = None
        phase_fixed = None
        phase_limits = res_dict.get("phase_limits", None)

        if isinstance(phase_raw, dict):
            phase_val = float(phase_raw.get("value", 0.0))
            phase_fixed = bool(phase_raw.get("fixed", phase_raw.get("is_fixed", True)))
            if "limits" in phase_raw:
                phase_limits = phase_raw["limits"]
        elif isinstance(phase_raw, (int, float)):
            phase_val = float(phase_raw)
            if "phase_fixed" in res_dict:
                phase_fixed = bool(res_dict["phase_fixed"])
            elif "fixed" in res_dict:
                phase_fixed = bool(res_dict["fixed"])
            else:
                # Direct scalar e.g. phase: 0.0 -> fixed to zero as reference
                phase_fixed = True
        elif "phase_fixed" in res_dict:
            phase_fixed = bool(res_dict["phase_fixed"])
            phase_val = float(res_dict.get("fixed_phase", res_dict.get("phase", 0.0)))
        elif "fixed_phase" in res_dict and res_dict["fixed_phase"] is not None:
            phase_val = float(res_dict["fixed_phase"])
            phase_fixed = True

        if phase_fixed is None:
            phase_fixed = False

        res_fixed_val = phase_val if phase_fixed else None

        res_shared_phase_raw = None
        for key in ("shared_phase", "shared_phase_by_reflectivity", "shared_phase_per_reflectivity", "shared_phases_by_reflectivity"):
            if key in res_dict:
                res_shared_phase_raw = res_dict[key]
                break

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
                phase=phase_val,
                phase_fixed=phase_fixed,
                fixed_phase=res_fixed_val,
                phase_limits=phase_limits,
                shared_phase=res_shared_phase_raw,
            )
        )

    # Wave-level general phase
    wave_phase_raw = wave_dict.get("phase", None)
    wave_phase_fixed = bool(wave_dict.get("phase_fixed", False))
    wave_phase_val = None

    if isinstance(wave_phase_raw, dict):
        wave_phase_fixed = bool(wave_phase_raw.get("fixed", wave_phase_raw.get("is_fixed", True)))
        wave_phase_val = float(wave_phase_raw.get("value", 0.0))
    elif isinstance(wave_phase_raw, (int, float)):
        wave_phase_val = float(wave_phase_raw)
        if "phase_fixed" in wave_dict:
            wave_phase_fixed = bool(wave_dict["phase_fixed"])
        else:
            wave_phase_fixed = True

    if wave_phase_val is not None and wave_phase_fixed:
        for r in resonances:
            if r.phase is None:
                r.phase = wave_phase_val
                r.phase_fixed = True
                r.fixed_phase = wave_phase_val

    # Init config
    init_dict = wave_dict.get("init", {})

    # Wave-level shared phase (optional wave-level default)
    wave_shared_phase_raw = None
    for key in ("shared_phase", "shared_phase_by_reflectivity", "shared_phase_per_reflectivity", "shared_phases_by_reflectivity"):
        if key in wave_dict:
            wave_shared_phase_raw = wave_dict[key]
            break
        if key in init_dict:  # legacy fallback
            wave_shared_phase_raw = init_dict[key]
            break

    # If set at wave level, propagate as default to resonances that didn't specify it
    if wave_shared_phase_raw is not None:
        for r in resonances:
            if r.shared_phase is None:
                r.shared_phase = wave_shared_phase_raw

    # Wave-level shared phase summary value
    shared_phase_summary = wave_shared_phase_raw
    if shared_phase_summary is None:
        for r in resonances:
            if r.shared_phase is not None:
                shared_phase_summary = r.shared_phase
                break

    def _parse_shared_phase_opt(val) -> int:
        if isinstance(val, bool):
            return 2 if val else 0
        elif isinstance(val, str):
            s = val.strip().lower()
            if s in ("true", "1", "yes", "reflectivity", "by_reflectivity", "per_reflectivity", "shared_reflectivity"):
                return 2
            elif s in ("l", "by_l", "same_l", "same-l"):
                return 1
            elif s in ("false", "0", "no", "none", "disabled"):
                return 0
        elif isinstance(val, (int, float)):
            return int(val)
        return 0

    option = int(init_dict["option"]) if "option" in init_dict else None
    if option is None:
        res_opts = [_parse_shared_phase_opt(r.shared_phase) for r in resonances if r.shared_phase is not None]
        if res_opts:
            option = max(res_opts)
        elif wave_shared_phase_raw is not None:
            option = _parse_shared_phase_opt(wave_shared_phase_raw)
        elif (any(r.phase_fixed for r in resonances) or (wave_phase_val is not None and wave_phase_fixed)):
            # Default to shared phase per reflectivity when a phase is fixed
            option = 2
        else:
            option = 0

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

    phase_dict = init_dict.get("phase", wave_dict.get("phase", {}))
    if isinstance(phase_dict, (int, float)):
        phase_rule = ParameterInitRule(
            strategy=InitStrategy.FIXED_VALUE,
            value=float(phase_dict),
            is_fixed=bool(init_dict.get("phase_fixed", wave_dict.get("phase_fixed", True))),
        )
    elif isinstance(phase_dict, dict):
        phase_strat = phase_dict.get("strategy", "fixed").lower()
        phase_rule = ParameterInitRule(
            strategy=InitStrategy.FIXED_VALUE if phase_strat == "fixed" else InitStrategy.RANDOM,
            min_val=float(phase_dict.get("min", -3.14159265)),
            max_val=float(phase_dict.get("max", 3.14159265)),
            value=float(phase_dict.get("value", 0.0)),
            is_fixed=bool(phase_dict.get("fixed", init_dict.get("phase_fixed", wave_dict.get("phase_fixed", False)))),
        )
    else:
        phase_rule = ParameterInitRule(
            strategy=InitStrategy.FIXED_VALUE,
            value=0.0,
            is_fixed=False,
        )

    init_cfg = WaveInitConfig(
        strategy=init_strat,
        fixed_zero=fixed_zero,
        file_path=init_dict.get("file_path"),
        option=option,
        shared_phase=shared_phase_summary if shared_phase_summary is not None else False,
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
        phase=wave_phase_val,
        phase_fixed=wave_phase_fixed,
        shared_phase=shared_phase_summary,
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
        tag=str(proj_dict.get("tag", data.get("tag", ""))),
    )

    # Environment (optional, defaults to environment variables or system paths)
    env_dict = data.get("environment") or {}
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

    # Resolve t-bins (single or list, configurable via t_bins, tbins, t_bin, tbin, t_bin_dir, tbins_dir)
    raw_tbins = None
    for k in ("t_bins", "tbins"):
        if k in ds_dict and ds_dict[k] is not None:
            raw_tbins = ds_dict[k]
            break

    raw_single = None
    for k in ("t_bin", "tbin", "t_bin_dir", "tbins_dir"):
        if k in ds_dict and ds_dict[k] is not None:
            raw_single = ds_dict[k]
            break

    t_bins: List[str] = []
    if raw_tbins is not None:
        t_bins = parse_t_bins_input(raw_tbins)
    elif raw_single is not None:
        t_bins = parse_t_bins_input(raw_single)

    t_bin = t_bins[0] if t_bins else None

    dataset = DatasetConfig(
        moments_dir=ds_dict.get("moments_dir", ""),
        t_bin=t_bin,
        t_bins=t_bins,
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
    n_lowest_ms = ms_dict.get("n_lowest_chi2", 10)
    n_tbin_workers = (
        int(ms_dict["n_tbin_workers"])
        if "n_tbin_workers" in ms_dict and ms_dict["n_tbin_workers"] is not None
        else None
    )
    multi_start = MultiStartConfig(
        enabled=bool(ms_dict.get("enabled", True)),
        n_starts=int(ms_dict.get("n_starts", 100)),
        n_threads=int(ms_dict["n_threads"]) if "n_threads" in ms_dict and ms_dict["n_threads"] is not None else get_default_max_workers(),
        n_tbin_workers=n_tbin_workers,
        output_dir=ms_dict.get("output_dir", "logs"),
        plateau_delta_chi2=float(ms_dict.get("plateau_delta_chi2", 0.05)),
        seed_offset=int(ms_dict.get("seed_offset", 1)),
        n_lowest_chi2=int(n_lowest_ms),
    )

    # Plotting
    plot_dict = data.get("plotting", {})
    n_lowest_plot = plot_dict.get("n_lowest_chi2", n_lowest_ms)
    plotting = PlottingConfig(
        enabled=bool(plot_dict.get("enabled", True)),
        n_lowest_chi2=int(n_lowest_plot),
        draw_best_fits=bool(plot_dict.get("draw_best_fits", True)),
        unnormalized_moments_file=plot_dict.get("unnormalized_moments_file"),
        moments=plot_dict.get("moments", ["all"]),
        formats=plot_dict.get("formats", ["png", "pdf"]),
        style=str(plot_dict.get("style", "default")),
        show_phase_ambiguity=bool(plot_dict.get("show_phase_ambiguity", True)),
        posterior_density_bands=bool(plot_dict.get("posterior_density_bands", True)),
    )

    # MCMCPropagation
    mcmc_dict = data.get("mcmc_propagation", {})
    n_tbin_workers_mcmc = (
        int(mcmc_dict["n_tbin_workers"])
        if "n_tbin_workers" in mcmc_dict and mcmc_dict["n_tbin_workers"] is not None
        else None
    )
    mcmc_propagation = MCMCPropagationConfig(
        enabled=bool(mcmc_dict.get("enabled", False)),
        mcmc_tree_path=mcmc_dict.get("mcmc_tree_path", ""),
        burn_in=int(mcmc_dict.get("burn_in", 0)),
        max_entries=int(mcmc_dict.get("max_entries", 200)),
        stride=int(mcmc_dict.get("stride", 1)),
        n_workers=int(mcmc_dict["n_workers"]) if "n_workers" in mcmc_dict and mcmc_dict["n_workers"] is not None else get_default_max_workers(),
        n_tbin_workers=n_tbin_workers_mcmc,
        output_file=mcmc_dict.get("output_file", "mcmc_propagation_results.root"),
        log_dir=mcmc_dict.get("log_dir", "logs_mcmc"),
    )

    # CrossSection
    cs_dict = data.get("cross_section", {})
    mass_lims = cs_dict.get("mass_limits", [0.8, 2.0])
    cross_section = CrossSectionConfig(
        enabled=bool(cs_dict.get("enabled", True)),
        resonance=cs_dict.get("resonance", "a2_1320"),
        branching_ratio=float(cs_dict.get("branching_ratio", 0.0565)),
        luminosity=float(cs_dict.get("luminosity", 257.4)),
        luminosity_unit=cs_dict.get("luminosity_unit", "pb^-1"),
        t_min=float(cs_dict["t_min"]) if "t_min" in cs_dict and cs_dict["t_min"] is not None else None,
        t_max=float(cs_dict["t_max"]) if "t_max" in cs_dict and cs_dict["t_max"] is not None else None,
        delta_t=float(cs_dict["delta_t"]) if "delta_t" in cs_dict and cs_dict["delta_t"] is not None else None,
        mass_limits=(float(mass_lims[0]), float(mass_lims[1])) if len(mass_lims) == 2 else (0.8, 2.0),
        mass_bin_width=float(cs_dict["mass_bin_width"]) if "mass_bin_width" in cs_dict and cs_dict["mass_bin_width"] is not None else None,
        theory_band_file=cs_dict.get("theory_band_file"),
        phase1_comparison=bool(cs_dict.get("phase1_comparison", True)),
        mcmc_file=cs_dict.get("mcmc_file"),
        propagate_mcmc=bool(cs_dict.get("propagate_mcmc", True)),
    )

    # Top-level shared phase options
    shared_phases_by_refl = bool(
        data.get("shared_phases_by_reflectivity", data.get("shared_phase_per_reflectivity", False))
    )

    if shared_phases_by_refl:
        for w in waves:
            if w.type == WaveType.MASS_DEPENDENT and w.init.option == 0:
                w.init.option = 2
                w.init.shared_phase = True
                w.shared_phase = True
                for r in w.resonances:
                    if r.shared_phase is None:
                        r.shared_phase = True

    return FitConfig(
        project=project,
        environment=environment,
        dataset=dataset,
        chi2=chi2,
        waves=waves,
        minimizer=minimizer,
        multi_start=multi_start,
        mcmc_propagation=mcmc_propagation,
        plotting=plotting,
        cross_section=cross_section,
        shared_phases_by_reflectivity=shared_phases_by_refl,
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
        if config.mcmc_propagation.burn_in < 0:
            errors.append("mcmc_propagation burn_in must be non-negative.")
        if config.mcmc_propagation.stride < 1:
            errors.append("mcmc_propagation stride must be at least 1.")

    return errors


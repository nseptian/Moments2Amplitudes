"""
Execution driver and ROOT macro generator for Moments2Amplitudes fits.
Translates declarative YAML FitConfig into C++ engine execution, orchestrates
parallel multi-start fits, and propagates MCMC posterior uncertainties.
"""

from __future__ import annotations
import os
import sys
import copy
import time
import subprocess
import threading
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
import numpy as np

import re
from .schema import (
    FitConfig,
    WaveType,
    ModelType,
    InitStrategy,
)
from .rundir import RunDirectory, get_run_directory, resolve_init_file_for_tbin
from .plotter import plot_multistart_landscape, plot_amplitudes_and_bands, run_all_diagnostics_and_plots
from .summary import FitSummary
from .environment import (
    get_root_binary,
    get_root_subprocess_env,
    resolve_rootsys,
    resolve_brufit,
    get_repo_root,
    get_default_max_workers,
)
from .compiler import ensure_aclic_updated


def resolve_mcmc_tree_path(raw_path: str, t_bin: Optional[str] = None) -> str:
    """Resolve MCMC tree file path for a specific t-bin."""
    if not raw_path or not t_bin:
        return raw_path
    if "{tbin}" in raw_path or "{t_bin}" in raw_path or "{t-bin}" in raw_path:
        return raw_path.replace("{t_bin}", t_bin).replace("{tbin}", t_bin).replace("{t-bin}", t_bin)
    # If the path explicitly embeds a t-bin name like /t010020/, replace it with /<t_bin>/
    subbed = re.sub(r"/t\d{6}/", f"/{t_bin}/", raw_path)
    if Path(subbed).is_file():
        return subbed
    return raw_path


def _translate_subexpr(expr: str) -> str:
    res = expr.strip()
    # Power replacement: X^2 -> (X * X)
    res = re.sub(r'([A-Za-z0-9_]+)\^2', r'(\1 * \1)', res)
    res = re.sub(r'([A-Za-z0-9_]+)\^([0-9]+)', r'std::pow(\1, \2)', res)

    # Breakup momentum / phase space factors
    res = res.replace("rho_etapi", "m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.547853, 0.13957)")
    res = res.replace("rho_KK", "m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.493677, 0.493677)")
    res = res.replace("rho_etaprimepi", "m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.95778, 0.13957)")

    # Complex imaginary unit replacements
    res = re.sub(r'\bi\s*\*\s*', 'std::complex<double>(0.0, 1.0) * ', res)
    res = re.sub(r'\bi\b', 'std::complex<double>(0.0, 1.0)', res)
    res = re.sub(r'\bexp\(', 'std::exp(', res)
    res = re.sub(r'\bsqrt\(', 'std::sqrt(', res)
    return res


def translate_formula_to_cpp(formula: str, parameters: List[str]) -> Tuple[List[str], str]:
    """
    Translate a mathematical formula into C++ statements and return expression.
    Supports optional kinematic factor prefixes, e.g. sqrt(rho_etapi) * k / (M^2 - s - i*M*width).
    """
    param_lines = ["const double s = mass * mass;"]
    for i, p in enumerate(parameters):
        param_lines.append(f"const double {p} = (p.size() > {i}) ? p[{i}] : 0.0;")

    # Check if formula matches Two Breit-Wigner form (with optional prefix, e.g. sqrt(rho_etapi) * (...))
    two_bw_regex = re.compile(
        r'^\s*(?:(?P<prefix>.+?)\s*\*\s*)?(?:\(\s*)?'
        r'(?:\(\s*)?(?P<k1>[A-Za-z0-9_]+)\s*/\s*\(\s*(?P<m1>[A-Za-z0-9_]+)(?:\^2|\s*\*\s*(?P=m1))\s*-\s*s\s*-\s*i\s*\*\s*(?:(?P=m1)\s*\*\s*(?P<w1>[A-Za-z0-9_]+)|(?P<w1_b>[A-Za-z0-9_]+)\s*\*\s*(?P=m1))\s*\)(?:\s*\))?'
        r'(?:\s*\*\s*exp\(\s*i\s*\*\s*(?P<phase1>[A-Za-z0-9_]+)\s*\))?\s*\+\s*'
        r'(?:\(\s*)?(?P<k2>[A-Za-z0-9_]+)\s*/\s*\(\s*(?P<m2>[A-Za-z0-9_]+)(?:\^2|\s*\*\s*(?P=m2))\s*-\s*s\s*-\s*i\s*\*\s*(?:(?P=m2)\s*\*\s*(?P<w2>[A-Za-z0-9_]+)|(?P<w2_b>[A-Za-z0-9_]+)\s*\*\s*(?P=m2))\s*\)(?:\s*\))?'
        r'(?:\s*\*\s*exp\(\s*i\s*\*\s*(?P<phase2>[A-Za-z0-9_]+)\s*\))?'
        r'(?:\s*\))?\s*$'
    )
    two_bw_m = two_bw_regex.match(formula)
    if two_bw_m:
        prefix = two_bw_m.group("prefix")
        k1 = two_bw_m.group("k1")
        m1 = two_bw_m.group("m1")
        w1 = two_bw_m.group("w1") or two_bw_m.group("w1_b")
        k2 = two_bw_m.group("k2")
        m2_var = two_bw_m.group("m2")
        w2 = two_bw_m.group("w2") or two_bw_m.group("w2_b")
        phase1_var = two_bw_m.group("phase1")
        phase2_var = two_bw_m.group("phase2")
        if phase1_var and phase2_var and "phase1" in parameters and "phase2" in parameters:
            # When both phase1 and phase2 are listed in parameters, support both 8-param custom
            # callers and 7-param TwoBreitWigner callers where phase1 is reference (fixed to 0)
            # and phase2 is the relative phase passed as p[6].
            param_lines = [
                "const double s = mass * mass;",
                f"const double {k1} = (p.size() > 0) ? p[0] : 0.0;",
                f"const double {m1} = (p.size() > 1) ? p[1] : 0.0;",
                f"const double {w1} = (p.size() > 2) ? p[2] : 0.0;",
                f"const double {phase1_var} = (p.size() > 7) ? p[3] : 0.0; // Reference phase (fixed to 0)",
                f"const double {k2} = (p.size() > 7) ? p[4] : ((p.size() > 3) ? p[3] : 0.0);",
                f"const double {m2_var} = (p.size() > 7) ? p[5] : ((p.size() > 4) ? p[4] : 0.0);",
                f"const double {w2} = (p.size() > 7) ? p[6] : ((p.size() > 5) ? p[5] : 0.0);",
                f"const double {phase2_var} = (p.size() > 7) ? p[7] : ((p.size() > 6) ? p[6] : 0.0); // Relative phase",
            ]

        if phase1_var:
            bw1_str = f"({k1} / std::complex<double>(({m1} * {m1}) - s, -{m1} * {w1})) * std::polar(1.0, {phase1_var})"
        else:
            bw1_str = f"{k1} / std::complex<double>(({m1} * {m1}) - s, -{m1} * {w1})"
        if phase2_var:
            bw2_str = f"({k2} / std::complex<double>(({m2_var} * {m2_var}) - s, -{m2_var} * {w2})) * std::polar(1.0, {phase2_var})"
        else:
            bw2_str = f"{k2} / std::complex<double>(({m2_var} * {m2_var}) - s, -{m2_var} * {w2})"
        core_expr = f"{bw1_str} + {bw2_str}"
        if prefix:
            trans_prefix = _translate_subexpr(prefix)
            return param_lines, f"{trans_prefix} * ({core_expr})"
        return param_lines, core_expr

    # Check if formula matches canonical Breit-Wigner form (with optional prefix, e.g. sqrt(rho_etapi) * k / ...)
    bw_regex = re.compile(
        r'^\s*(?:(?P<prefix>.+?)\s*\*\s*)?(?:\(\s*)?'
        r'(?P<k>[A-Za-z0-9_]+)\s*/\s*\(\s*(?P<m>[A-Za-z0-9_]+)(?:\^2|\s*\*\s*(?P=m))\s*-\s*s\s*-\s*i\s*\*\s*(?:(?P=m)\s*\*\s*(?P<w>[A-Za-z0-9_]+)|(?P<w_b>[A-Za-z0-9_]+)\s*\*\s*(?P=m))\s*\)'
        r'(?:\s*\))?\s*$'
    )
    bw_m = bw_regex.match(formula)
    if bw_m:
        prefix = bw_m.group("prefix")
        k_var = bw_m.group("k")
        m_var = bw_m.group("m")
        w_var = bw_m.group("w") or bw_m.group("w_b")
        bw_str = f"{k_var} / std::complex<double>(({m_var} * {m_var}) - s, -{m_var} * {w_var})"
        if prefix:
            trans_prefix = _translate_subexpr(prefix)
            return param_lines, f"{trans_prefix} * ({bw_str})"
        return param_lines, bw_str

    # Check if formula matches canonical Flatte form: k / (M0^2 - s - i * (g_etapi^2 * rho_etapi + ...))
    if "flatte" in formula.lower() or ("g_etapi" in formula and "rho" in formula):
        param_lines.extend([
            "const double g_etapi_sqr = g_etapi * g_etapi;",
            "const double g_KK_sqr = (g_KK < 0.0) ? (g_etapi_sqr * 1.05) : (g_KK * g_KK);",
            "const double g_etaprimepi_sqr = (g_etaprimepi < 0.0) ? (g_etapi_sqr * 0.772) : (g_etaprimepi * g_etaprimepi);",
            "const double rho_etapi = m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.547853, 0.13957);",
            "const double rho_KK = m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.493677, 0.493677);",
            "const double rho_etaprimepi = m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.95778, 0.13957);",
            "const double gamma_total = g_etapi_sqr * rho_etapi + g_KK_sqr * rho_KK + g_etaprimepi_sqr * rho_etaprimepi;",
        ])
        mass_param = "Mass" if "Mass" in parameters else ("M0" if "M0" in parameters else "M")
        k_param = "k" if "k" in parameters else "k_fl"
        flatte_str = f"{k_param} / std::complex<double>(({mass_param} * {mass_param}) - s, -gamma_total)"
        prefix_m = re.match(r'^\s*(?P<prefix>.+?)\s*\*\s*(?:\(\s*)?[A-Za-z0-9_]+\s*/\s*\(', formula)
        if prefix_m:
            trans_prefix = _translate_subexpr(prefix_m.group("prefix"))
            return param_lines, f"{trans_prefix} * ({flatte_str})"
        return param_lines, flatte_str

    return param_lines, _translate_subexpr(formula)


def _get_two_bw_wave_name(w: WaveConfig) -> str:
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


def generate_root_macro(config: FitConfig, output_macro_path: Path) -> Path:
    """
    Generate a self-contained, robust C++ ROOT macro that initializes and runs
    MassDependentFitter according to the FitConfig specifications.
    """
    output_macro_path = Path(output_macro_path).resolve()
    output_macro_path.parent.mkdir(parents=True, exist_ok=True)

    mass_bins = config.get_mass_bins()
    n_bins = len(mass_bins)
    mass_bins_str = ", ".join(f"{b:.6f}" for b in mass_bins)

    # Moments config
    inc_l_str = ", ".join(str(l) for l in config.chi2.included_l)
    exc_l_str = ", ".join(str(l) for l in config.chi2.excluded_l)

    # Classify waves into mass-dependent and piecewise
    md_waves = [w for w in config.waves if w.type == WaveType.MASS_DEPENDENT]
    mi_waves = [w for w in config.waves if w.type == WaveType.PIECEWISE]

    mi_l_list = [w.l for w in mi_waves]
    mi_l_str = ", ".join(str(l) for l in mi_l_list)

    max_wave_l = max((w.l for w in config.waves), default=2)

    # Determine project directory and brufit path
    repo_root = os.environ.get("MOMS2AMPS", str(Path(__file__).resolve().parent.parent.parent))
    brufit_path = resolve_brufit(config)

    # Explicit formula functions at file scope to avoid Cling lambda-mangler issues
    formula_functions: List[str] = []
    formula_registrations: List[str] = []

    for w in md_waves:
        if not w.formula:
            continue
        form_esc = w.formula.replace('"', '\\"')
        param_lines, expr = translate_formula_to_cpp(w.formula, w.parameters)
        param_code = "\n    ".join(param_lines)

        if w.model == ModelType.TWO_BREIT_WIGNER:
            wave_name = _get_two_bw_wave_name(w)
            safe_name = f"{w.l}_{wave_name.replace('-', '_')}"
            func_name = f"eval_formula_{safe_name}"

            formula_functions.extend([
                f"// Explicit formula for L={w.l} {wave_name}: {w.formula}",
                f"std::complex<double> {func_name}(double mass, const std::vector<double>& p) {{",
                f"    {param_code}",
                f"    return {expr};",
                "}",
                "",
            ])
            formula_registrations.extend([
                f'    // Explicit formula configured for L={w.l} {wave_name}: {w.formula}',
                f'    fitter.SetWaveCustomFormula({w.l}, "{wave_name}", {func_name}, "{form_esc}");',
            ])
        else:
            target_res = w.resonances if w.resonances else [None]
            for res in target_res:
                r_name = res.name if res else f"wave_{w.l}"
                safe_name = f"{w.l}_{r_name.replace('-', '_')}"
                func_name = f"eval_formula_{safe_name}"

                formula_functions.extend([
                    f"// Explicit formula for L={w.l} {r_name}: {w.formula}",
                    f"std::complex<double> {func_name}(double mass, const std::vector<double>& p) {{",
                    f"    {param_code}",
                    f"    return {expr};",
                    "}",
                    "",
                ])
                formula_registrations.extend([
                    f'    // Explicit formula configured for L={w.l} {r_name}: {w.formula}',
                    f'    fitter.SetWaveCustomFormula({w.l}, "{r_name}", {func_name}, "{form_esc}");',
                ])

    lines: List[str] = [
        "// Auto-generated ROOT driver macro by Moments2Amplitudes",
        f"// Project: {config.project.name}",
        "#include <iostream>",
        "#include <vector>",
        "#include <map>",
        "#include <TString.h>",
        "#include <TSystem.h>",
        "#include <TFile.h>",
        "#include <TTree.h>",
        "#include <TStopwatch.h>",
        "",
        '#include "Setup.h"',
        '#include "ConfigureAmpsNoValues.C"',
        '#include "MomentHelper.h"',
        '#include "MassDependentMoments.h"',
        '#include "MassDependentFitter.h"',
        "",
        "using namespace m2pw;",
        "using namespace std;",
        "",
        "void EnsureLoaded() {",
        '    if (TClass::GetClass("m2pw::MassDependentFitter")) return;',
        f'    gSystem->Setenv("MOMS2AMPS", "{repo_root}");',
        f'    gSystem->Setenv("BRUFIT", "{brufit_path}");',
        f'    TString loadScript = Form("{repo_root}/macros/Load.C");',
        '    if (gSystem->AccessPathName(loadScript) == 0) {',
        '        gROOT->ProcessLine(Form(".x %s", loadScript.Data()));',
        "    }",
        "}",
        "",
    ]

    lines.extend(formula_functions)

    no_use_items = getattr(config.chi2, "no_use", None) or ["H_3"]
    no_use_str = ", ".join(f'"{x}"' for x in no_use_items)

    lines.extend([
        "void RunWorker(const char* mode = \"single\", int seed = 1, int mcmcStart = 0, int mcmcEnd = -1, const char* outDir = \".\", const char* bestParamsFile = \"\") {",
        "    EnsureLoaded();",
        "    TStopwatch timer;",
        "    timer.Start();",
        "",
        "    // Setup mass bins",
        f"    const std::vector<double> massBins = {{{mass_bins_str}}};",
        f"    const double lMax = {config.chi2.max_l};",
        "    const double noise = 0.0;",
        f"    const std::vector<TString> noUse = {{{no_use_str}}};",
        "",
        "    // Setup MassDependentMoments",
        f'    const TString momentsDir = "{config.dataset.get_moments_dir()}";',
        f'    const TString resultsFile = "{config.dataset.results_filename}";',
        f'    const TString binPrefix = "{config.dataset.bin_prefix}";',
        f'    const TString varBinPrefix = "{config.dataset.var_bin_prefix}";',
        "",
        "    std::cout << \"[Worker] Loading dataset moments across \" << massBins.size() << \" mass bins...\" << std::endl;",
        "    MassDependentMoments moms(massBins, momentsDir, resultsFile, binPrefix, varBinPrefix);",
        "",
        "    // Configure mass-dependent models and moments selection",
        "    MassDependentFitter::MassDependenceConfig massDepConfig;",
        "    massDepConfig.massDependentWaves.clear();",
        f"    massDepConfig.massIndependentL = {{{mi_l_str}}};",
    ])

    # Mass-dependent wave models
    for w in md_waves:
        if w.model == ModelType.CUSTOM:
            m_type = "MassDependentFitter::MassDependenceConfig::ModelType::CustomFormula"
        elif w.model == ModelType.TWO_BREIT_WIGNER:
            m_type = "MassDependentFitter::MassDependenceConfig::ModelType::TwoBreitWignerPlusPolynomial"
        elif w.model == ModelType.FLATTE_BREIT_WIGNER:
            m_type = "MassDependentFitter::MassDependenceConfig::ModelType::FlatteBreitWignerPlusPolynomial"
        elif w.model == ModelType.FLATTE:
            m_type = "MassDependentFitter::MassDependenceConfig::ModelType::FlattePlusPolynomial"
        elif w.model == ModelType.POLYNOMIAL:
            m_type = "MassDependentFitter::MassDependenceConfig::ModelType::Polynomial"
        else:
            m_type = "MassDependentFitter::MassDependenceConfig::ModelType::BreitWignerPlusPolynomial"

        form_str = (w.formula or "").replace('"', '\\"')
        param_list_str = ", ".join(f'"{p}"' for p in w.parameters) if w.parameters else ""

        if w.model == ModelType.TWO_BREIT_WIGNER:
            wave_name = _get_two_bw_wave_name(w)
            lines.append(f'    massDepConfig.massDependentWaves[{w.l}].push_back({{"{wave_name}", {m_type}, {w.poly_order}, "{form_str}", {{{param_list_str}}}}});')
        else:
            target_res = w.resonances if w.resonances else [None]
            for res in target_res:
                r_name = res.name if res else f"wave_{w.l}"
                lines.append(f'    massDepConfig.massDependentWaves[{w.l}].push_back({{"{r_name}", {m_type}, {w.poly_order}, "{form_str}", {{{param_list_str}}}}});')

    # Moments selection
    lines.extend([
        "",
        "    MassDependentFitter::MomentsConfig momentsConfig;",
        f"    momentsConfig.includeAll = {str(config.chi2.include_all).lower()};",
        f"    momentsConfig.includedL = {{{inc_l_str}}};",
        f"    momentsConfig.excludedL = {{{exc_l_str}}};",
        "",
        "    // Create fitter instance with configured partial wave setup",
        f"    auto& fitSetup = ConfigureAmpsNoValues({max_wave_l}, {max_wave_l}, 2);",
        "    MassDependentFitter fitter(fitSetup, massBins, lMax, noise, noUse, massDepConfig, momentsConfig);",
        f"    fitter.SetUseCovarianceMatrix({str(config.dataset.use_covariance).lower()});",
        f"    fitter.SetMinimizerPrintLevel({config.minimizer.print_level});",
        f"    fitter.SetMaxFunctionCalls({config.minimizer.max_calls});",
        f"    fitter.SetMaxIterations({config.minimizer.max_iterations});",
        f"    fitter.SetMinimizerTolerance({config.minimizer.tolerance});",
        f"    fitter.SetMinimizerStrategy({config.minimizer.strategy});",
        "    fitter.SetEquationValues(moms);",
    ])

    # Assign explicit formulas to wave models
    lines.extend(formula_registrations)

    lines.extend([
        "",
        "    // Initialize ParameterManager",
        "    MassDependentFitter::ParameterManager paramManager(fitter);",
    ])

    # Wave initialization rules for piecewise waves
    for w in mi_waves:
        if w.init.fixed_zero:
            lines.append(f"    paramManager.SetWaveInitRule({w.l}, MassDependentFitter::WaveInitRule(MassDependentFitter::WaveInitRule::Strategy::FixedZero, 0.0, 0.0, 0.0, true));")
        elif w.init.strategy == InitStrategy.FIXED_VALUE:
            lines.append(f"    paramManager.SetWaveInitRule({w.l}, MassDependentFitter::WaveInitRule(MassDependentFitter::WaveInitRule::Strategy::FixedValue, 0.0, 0.0, {w.init.magnitude.value}, true));")
        elif w.init.strategy == InitStrategy.RANDOM or w.init.magnitude.strategy == InitStrategy.RANDOM:
            lines.append(f"    paramManager.SetWaveInitRule({w.l}, MassDependentFitter::WaveInitRule(MassDependentFitter::WaveInitRule::Strategy::RandomPositiveMagnitudeZeroPhase, {w.init.magnitude.min_val}, {w.init.magnitude.max_val}, 0.0, false));")

    # Add piecewise parameters
    mi_from_file = [w for w in mi_waves if w.init.strategy == InitStrategy.FROM_FILE and w.init.file_path]
    mi_random = [w for w in mi_waves if not (w.init.strategy == InitStrategy.FROM_FILE and w.init.file_path)]

    for w in mi_from_file:
        init_file = resolve_init_file_for_tbin(w.init.file_path, config.dataset.t_bin, repo_root)
        lines.append(f'    paramManager.AddMassIndependentParameter({{{w.l}}}, "{init_file}");')

    if mi_random:
        rand_l_str = ", ".join(str(w.l) for w in mi_random)
        lines.append(f"    paramManager.AddMassIndependentParameter({{{rand_l_str}}}, seed);")

    # Fix phase for piecewise waves if phase is fixed
    for w in mi_waves:
        if w.init.phase.is_fixed or w.phase_fixed:
            val = w.init.phase.value if w.init.phase.is_fixed else (w.phase if w.phase is not None else 0.0)
            lines.extend([
                f"    // Fixed phase for piecewise wave L={w.l}",
                "    for (const double mb : massBins) {",
                f'        TString aPhaseName = Form("MI_%1.6f_aphi_{w.l}_0", mb);',
                f'        TString bPhaseName = Form("MI_%1.6f_bphi_{w.l}_0", mb);',
                "        if (paramManager.HasParameter(aPhaseName)) {",
                f"            paramManager.SetParameterInitialValue(aPhaseName, {val});",
                "            paramManager.FixParameter(aPhaseName);",
                "        }",
                "        if (paramManager.HasParameter(bPhaseName)) {",
                f"            paramManager.SetParameterInitialValue(bPhaseName, {val});",
                "            paramManager.FixParameter(bPhaseName);",
                "        }",
                "    }",
            ])

    # Add mass-dependent parameters
    for w in md_waves:
        if w.init.strategy == InitStrategy.FROM_FILE and w.init.file_path:
            init_file = resolve_init_file_for_tbin(w.init.file_path, config.dataset.t_bin, repo_root)
            lines.append(f'    paramManager.AddMassDependentParameter({{{w.l}}}, "{init_file}", {w.init.option});')
        else:
            lines.append(f"    paramManager.AddMassDependentParameter({{{w.l}}}, seed, {w.init.option});")

        # Set resonance parameters defaults/limits/fixes
        if w.model == ModelType.TWO_BREIT_WIGNER and len(w.resonances) >= 2:
            res1 = w.resonances[0]
            res2 = w.resonances[1]
            wave_name = _get_two_bw_wave_name(w)

            # Resonance 1 (e.g. a2_1320)
            res1_mass_name = f"MD_shared_{wave_name}_Mass1"
            res1_width_name = f"MD_shared_{wave_name}_Width1"
            lines.append(f'    paramManager.SetParameterInitialValue("{res1_mass_name}", {res1.mass});')
            lines.append(f'    paramManager.SetParameterInitialValue("{res1_width_name}", {res1.width});')
            if res1.mass_fixed:
                lines.append(f'    paramManager.FixParameter("{res1_mass_name}");')
            elif res1.mass_limits and len(res1.mass_limits) == 2:
                lines.append(f'    paramManager.SetParameterLimits("{res1_mass_name}", {res1.mass_limits[0]}, {res1.mass_limits[1]});')

            if res1.width_fixed:
                lines.append(f'    paramManager.FixParameter("{res1_width_name}");')
            elif res1.width_limits and len(res1.width_limits) == 2:
                lines.append(f'    paramManager.SetParameterLimits("{res1_width_name}", {res1.width_limits[0]}, {res1.width_limits[1]});')

            for par_key, constr in res1.gaussian_constraints.items():
                if len(constr) == 2:
                    lines.append(f'    paramManager.SetGaussianConstraint("{par_key}", {constr[0]}, {constr[1]});')

            # Resonance 2 (e.g. a2_1700)
            res2_mass_name = f"MD_shared_{wave_name}_Mass2"
            res2_width_name = f"MD_shared_{wave_name}_Width2"
            lines.append(f'    paramManager.SetParameterInitialValue("{res2_mass_name}", {res2.mass});')
            lines.append(f'    paramManager.SetParameterInitialValue("{res2_width_name}", {res2.width});')
            if res2.mass_fixed:
                lines.append(f'    paramManager.FixParameter("{res2_mass_name}");')
            elif res2.mass_limits and len(res2.mass_limits) == 2:
                lines.append(f'    paramManager.SetParameterLimits("{res2_mass_name}", {res2.mass_limits[0]}, {res2.mass_limits[1]});')

            if res2.width_fixed:
                lines.append(f'    paramManager.FixParameter("{res2_width_name}");')
            elif res2.width_limits and len(res2.width_limits) == 2:
                lines.append(f'    paramManager.SetParameterLimits("{res2_width_name}", {res2.width_limits[0]}, {res2.width_limits[1]});')

            for par_key, constr in res2.gaussian_constraints.items():
                if len(constr) == 2:
                    lines.append(f'    paramManager.SetGaussianConstraint("{par_key}", {constr[0]}, {constr[1]});')

            # Phase handling for TwoBreitWigner:
            # Resonance 1 is the phase reference (phase fixed to 0 in built-in formula).
            # Resonance 2 carries the relative phase (represented by wave's global phase).
            phase_res = res2 if not res2.phase_fixed else (res1 if not res1.phase_fixed else res2)
            res_shared_opt = 2
            if phase_res.shared_phase is not None:
                if isinstance(phase_res.shared_phase, bool):
                    res_shared_opt = 2 if phase_res.shared_phase else 0
                elif isinstance(phase_res.shared_phase, str):
                    s = phase_res.shared_phase.strip().lower()
                    if s in ("true", "1", "yes", "reflectivity", "by_reflectivity", "per_reflectivity"):
                        res_shared_opt = 2
                    elif s in ("l", "by_l", "same_l", "same-l"):
                        res_shared_opt = 1
                    else:
                        res_shared_opt = 0
                elif isinstance(phase_res.shared_phase, (int, float)):
                    res_shared_opt = int(phase_res.shared_phase)
            elif w.init.option > 0:
                res_shared_opt = w.init.option

            pnames_list: List[str] = []
            if res_shared_opt == 2:
                pnames_list.extend([
                    f'"MD_shareda_{wave_name}_globalphase"',
                    f'"MD_sharedb_{wave_name}_globalphase"',
                    f'"MD_shareda_{res2.name}_globalphase"',
                    f'"MD_sharedb_{res2.name}_globalphase"',
                    f'"MD_shareda_{res1.name}_globalphase"',
                    f'"MD_sharedb_{res1.name}_globalphase"',
                ])
            elif res_shared_opt == 1:
                pnames_list.extend([
                    f'"MD_{wave_name}_globalphase"',
                    f'"MD_{res2.name}_globalphase"',
                    f'"MD_{res1.name}_globalphase"',
                ])
            pnames_list.extend([
                f'"MD_{wave_name}_phase"',
                f'"MD_{res2.name}_phase"',
                f'"MD_{res1.name}_phase"',
            ])
            cand_pnames_str = ",\n        ".join(pnames_list)

            if phase_res.phase_fixed:
                fix_val = phase_res.fixed_phase if phase_res.fixed_phase is not None else (phase_res.phase if phase_res.phase is not None else 0.0)
                lines.extend([
                    f'    // Fixed relative phase for {wave_name} (fixed to {fix_val})',
                    f'    for (const char* pName : {{',
                    f'        {cand_pnames_str}',
                    f'    }}) {{',
                    f'        if (paramManager.HasParameter(pName)) {{',
                    f'            paramManager.SetParameterInitialValue(pName, {fix_val});',
                    f'            paramManager.FixParameter(pName);',
                    f'        }}',
                    f'    }}',
                ])
            elif phase_res.phase is not None:
                init_val = phase_res.phase
                lines.extend([
                    f'    // Initial relative phase for {wave_name}',
                    f'    for (const char* pName : {{',
                    f'        {cand_pnames_str}',
                    f'    }}) {{',
                    f'        if (paramManager.HasParameter(pName)) {{',
                    f'            paramManager.SetParameterInitialValue(pName, {init_val});',
                    f'        }}',
                    f'    }}',
                ])

            if phase_res.phase_limits and len(phase_res.phase_limits) == 2 and not phase_res.phase_fixed:
                lines.extend([
                    f'    // Set relative phase limits for {wave_name}',
                    f'    for (const char* pName : {{',
                    f'        {cand_pnames_str}',
                    f'    }}) {{',
                    f'        if (paramManager.HasParameter(pName)) {{',
                    f'            paramManager.SetParameterLimits(pName, {phase_res.phase_limits[0]}, {phase_res.phase_limits[1]});',
                    f'        }}',
                    f'    }}',
                ])
        else:
            for res in w.resonances:
                lines.append(f'    paramManager.ConfigureBreitWignerDefaultsForWave("{res.name}", {res.mass}, {res.width});')
                res_mass_name = f"MD_shared_{res.name}_Mass"
                res_width_name = f"MD_shared_{res.name}_Width"
                if res.mass_fixed:
                    lines.append(f'    paramManager.FixParameter("{res_mass_name}");')
                elif res.mass_limits and len(res.mass_limits) == 2:
                    lines.append(f'    paramManager.SetParameterLimits("{res_mass_name}", {res.mass_limits[0]}, {res.mass_limits[1]});')

                if res.width_fixed:
                    lines.append(f'    paramManager.FixParameter("{res_width_name}");')
                elif res.width_limits and len(res.width_limits) == 2:
                    lines.append(f'    paramManager.SetParameterLimits("{res_width_name}", {res.width_limits[0]}, {res.width_limits[1]});')

                # Gaussian constraints
                for par_key, constr in res.gaussian_constraints.items():
                    if len(constr) == 2:
                        lines.append(f'    paramManager.SetGaussianConstraint("{par_key}", {constr[0]}, {constr[1]});')

                # Phase handling for resonance
                res_shared_opt = 2
                if res.shared_phase is not None:
                    if isinstance(res.shared_phase, bool):
                        res_shared_opt = 2 if res.shared_phase else 0
                    elif isinstance(res.shared_phase, str):
                        s = res.shared_phase.strip().lower()
                        if s in ("true", "1", "yes", "reflectivity", "by_reflectivity", "per_reflectivity"):
                            res_shared_opt = 2
                        elif s in ("l", "by_l", "same_l", "same-l"):
                            res_shared_opt = 1
                        else:
                            res_shared_opt = 0
                    elif isinstance(res.shared_phase, (int, float)):
                        res_shared_opt = int(res.shared_phase)
                elif w.init.option > 0:
                    res_shared_opt = w.init.option

                pnames_list: List[str] = []
                if res_shared_opt == 2:
                    pnames_list.extend([
                        f'"MD_shareda_{res.name}_globalphase"',
                        f'"MD_sharedb_{res.name}_globalphase"',
                        f'"MD_shareda_{res.name}_globalphase_flatte"',
                        f'"MD_shareda_{res.name}_globalphase_bw"',
                        f'"MD_sharedb_{res.name}_globalphase_flatte"',
                        f'"MD_sharedb_{res.name}_globalphase_bw"',
                        f'"MD_shareda_{res.name}_globalphase1"',
                        f'"MD_shareda_{res.name}_globalphase2"',
                        f'"MD_sharedb_{res.name}_globalphase1"',
                        f'"MD_sharedb_{res.name}_globalphase2"',
                    ])
                elif res_shared_opt == 1:
                    pnames_list.extend([
                        f'"MD_{res.name}_globalphase"',
                        f'"MD_{res.name}_globalphase_flatte"',
                        f'"MD_{res.name}_globalphase_bw"',
                        f'"MD_{res.name}_globalphase1"',
                        f'"MD_{res.name}_globalphase2"',
                    ])
                pnames_list.extend([
                    f'"MD_{res.name}_phase"',
                    f'"MD_{res.name}_phase1"',
                    f'"MD_{res.name}_phase2"',
                ])
                cand_pnames_str = ",\n        ".join(pnames_list)

                if res.phase_fixed:
                    fix_val = res.fixed_phase if res.fixed_phase is not None else (res.phase if res.phase is not None else 0.0)
                    lines.extend([
                        f'    // Fixed phase for resonance {res.name} (phases fixed to {fix_val})',
                        f'    for (const char* pName : {{',
                        f'        {cand_pnames_str}',
                        f'    }}) {{',
                        f'        if (paramManager.HasParameter(pName)) {{',
                        f'            paramManager.SetParameterInitialValue(pName, {fix_val});',
                        f'            paramManager.FixParameter(pName);',
                        f'        }}',
                        f'    }}',
                    ])
                elif res.phase is not None:
                    init_val = res.phase
                    lines.extend([
                        f'    // Initial phase for resonance {res.name}',
                        f'    for (const char* pName : {{',
                        f'        {cand_pnames_str}',
                        f'    }}) {{',
                        f'        if (paramManager.HasParameter(pName)) {{',
                        f'            paramManager.SetParameterInitialValue(pName, {init_val});',
                        f'        }}',
                        f'    }}',
                    ])

                if res.phase_limits and len(res.phase_limits) == 2 and not res.phase_fixed:
                    lines.extend([
                        f'    // Set phase limits for resonance {res.name}',
                        f'    for (const char* pName : {{',
                        f'        {cand_pnames_str}',
                        f'    }}) {{',
                        f'        if (paramManager.HasParameter(pName)) {{',
                        f'            paramManager.SetParameterLimits(pName, {res.phase_limits[0]}, {res.phase_limits[1]});',
                        f'        }}',
                        f'    }}',
                    ])

    # If bestParamsFile is provided in MCMC mode, load best parameters into paramManager
    lines.extend([
        "",
        '    if (TString(mode) == "mcmc" && TString(bestParamsFile).Length() > 0) {',
        '        TFile bestF(bestParamsFile, "READ");',
        "        if (bestF.IsOpen() && !bestF.IsZombie()) {",
        '            TTree* bTree = dynamic_cast<TTree*>(bestF.Get("result"));',
        "            if (bTree && bTree->GetEntries() > 0) {",
        "                bTree->GetEntry(0);",
        "                for (const auto& pName : paramManager.GetParameterNames()) {",
        "                    TLeaf* leaf = bTree->GetLeaf(pName);",
        "                    if (leaf) {",
        "                        paramManager.SetParameterInitialValue(pName, leaf->GetValue());",
        "                    }",
        "                }",
        "            }",
        "            bestF.Close();",
        "        }",
        "    }",
        "",
        '    if (TString(mode) == "single") {',
        "        std::cout << \"[Worker] Running minimization with seed=\" << seed << std::endl;",
        "        fitter.MinimizeChi2(paramManager);",
        "",
        "        const double chi2 = fitter.GetLastChi2();",
        "        const bool isValid = fitter.IsMinimizerValid();",
        "        const int status = fitter.GetMinimizerStatus();",
        "",
        '        TString outFile = Form("%s/MassDependentResults_%d.root", outDir, seed);',
        "        fitter.MakeResultTree(paramManager, outFile);",
        "",
        '        std::cout << "[RESULT] SEED=" << seed',
        '                  << " CHI2=" << Form("%.4f", chi2)',
        '                  << " VALID=" << (isValid ? 1 : 0)',
        '                  << " STATUS=" << status',
        '                  << " TIME=" << Form("%.2f", timer.RealTime())',
        '                  << " OUT=" << outFile << std::endl;',
        '    } else if (TString(mode) == "mcmc") {',
        "        std::cout << \"[Worker] Running MCMC propagation entries [\" << mcmcStart << \", \" << mcmcEnd << \")...\" << std::endl;",
        '        TString outFile = Form("%s/mcmc_worker_%d.root", outDir, seed);',
        f"        fitter.FitMCMCTree(paramManager, moms, outFile, mcmcStart, mcmcEnd, {config.mcmc_propagation.stride}, outDir);",
        '        std::cout << "[RESULT_MCMC] WORKER=" << seed',
        '                  << " ENTRIES=" << (mcmcEnd - mcmcStart)',
        '                  << " TIME=" << Form("%.2f", timer.RealTime())',
        '                  << " OUT=" << outFile << std::endl;',
        "    }",
        "}",
        "",
        f"void {output_macro_path.stem}(const char* mode = \"single\", int seed = 1, int mcmcStart = 0, int mcmcEnd = -1, const char* outDir = \".\", const char* bestParamsFile = \"\") {{",
        "    RunWorker(mode, seed, mcmcStart, mcmcEnd, outDir, bestParamsFile);",
        "}",
    ])

    if output_macro_path.stem != "worker":
        lines.extend([
            "",
            "void worker(const char* mode = \"single\", int seed = 1, int mcmcStart = 0, int mcmcEnd = -1, const char* outDir = \".\", const char* bestParamsFile = \"\") {",
            "    RunWorker(mode, seed, mcmcStart, mcmcEnd, outDir, bestParamsFile);",
            "}",
        ])

    with open(output_macro_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return output_macro_path


def execute_root_command(
    macro_path: Path,
    args: Tuple,
    log_path: Path,
    config: Optional[FitConfig] = None,
    timeout: Optional[int] = None,
) -> Tuple[int, str]:
    """Execute a single ROOT invocation in batch mode and capture output using ROOT 6.36."""
    arg_strs = []
    for a in args:
        if isinstance(a, str):
            arg_strs.append(f'"{a}"')
        else:
            arg_strs.append(str(a))
    call_expr = f'{macro_path.name}({", ".join(arg_strs)})'

    repo_root = Path(os.environ.get("MOMS2AMPS", str(Path(__file__).resolve().parent.parent.parent)))
    load_script = repo_root / "macros" / "Load.C"

    root_bin = get_root_binary(config)

    cmd = [
        root_bin,
        "-l",
        "-b",
        "-q",
    ]
    if load_script.is_file() and macro_path.resolve() != load_script.resolve():
        cmd.append(str(load_script))
    cmd.append(f"{macro_path}{call_expr[len(macro_path.name):]}")

    log_path.parent.mkdir(parents=True, exist_ok=True)
    env = get_root_subprocess_env(config)

    ret_code = -1
    log_content = ""
    try:
        with open(log_path, "w") as log_f:
            proc = subprocess.run(
                cmd,
                stdin=subprocess.DEVNULL,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                cwd=str(macro_path.parent),
                text=True,
                env=env,
                timeout=timeout,
            )
            ret_code = proc.returncode
    except subprocess.TimeoutExpired:
        with open(log_path, "a") as log_f:
            log_f.write(f"\n[ERROR] Process timed out after {timeout} seconds\n")
        log_content = f"[TIMEOUT] Execution timed out after {timeout} seconds"
        return -1, log_content

    # Read log
    if log_path.is_file():
        with open(log_path, "r") as f:
            log_content = f.read()

    return ret_code, log_content


def run_multi_start(
    config: FitConfig,
    n_starts: Optional[int] = None,
    n_threads: Optional[int] = None,
    output_dir: Optional[str] = None,
    config_source_path: Optional[Path | str] = None,
    tag: str = "",
) -> Dict[str, Any]:
    """
    Run parallel multi-start randomized fits to find the central global minimum.
    Generates summary plots, identifies the plateau, and saves central_best_params.root.
    """
    rundir = get_run_directory(config, output_dir)
    rundir.ensure_directories()

    # Ensure C++ ACLiC shared libraries are compiled and up-to-date before parallel workers launch
    ensure_aclic_updated(config)

    if config_source_path:
        rundir.save_config_snapshot(config_source_path)

    # Propagate t-bin and resolve initial parameters files
    repo_root = get_repo_root()
    for w in config.waves:
        if w.init.strategy == InitStrategy.FROM_FILE and w.init.file_path:
            init_file = resolve_init_file_for_tbin(w.init.file_path, config.dataset.t_bin, repo_root)
            w.init.file_path = str(init_file)

    half_cpus = get_default_max_workers()
    starts = n_starts if n_starts is not None else config.multi_start.n_starts
    threads = n_threads if n_threads is not None else (config.multi_start.n_threads or half_cpus)
    if n_threads is None:
        threads = min(threads, half_cpus)

    macro_path = rundir.worker_macro
    generate_root_macro(config, macro_path)

    header_tag = f" {tag}" if tag else ""
    print(f"\n=======================================================")
    print(f" Moments2Amplitudes Multi-Start Fit{header_tag}")
    print(f" Project:      {config.project.name}")
    if config.dataset.t_bin:
        print(f" t-Bin:        {config.dataset.t_bin}")
    print(f" Total Trials: {starts}")
    print(f" Threads:      {threads}")
    print(f" Run Dir:      {rundir.root}")
    print(f"=======================================================\n")

    seeds = [config.multi_start.seed_offset + i for i in range(starts)]
    results: List[Dict[str, Any]] = []

    # Parallel execution using concurrent.futures
    import concurrent.futures

    def run_seed_task(seed_val: int) -> Dict[str, Any]:
        log_f = rundir.get_seed_log_file(seed_val)
        seed_timeout = getattr(config.minimizer, "timeout", None) if config and hasattr(config, "minimizer") else None
        if seed_timeout is None:
            seed_timeout = 600
        ret_code, output = execute_root_command(
            macro_path,
            ("single", seed_val, 0, -1, str(rundir.fit_starts_dir), ""),
            log_f,
            config,
            timeout=seed_timeout,
        )

        # Parse [RESULT] line
        chi2 = 1e9
        is_valid = False
        status = -1
        out_root = str(rundir.fit_starts_dir / f"MassDependentResults_{seed_val}.root")

        for line in output.splitlines():
            if "[RESULT]" in line and f"SEED={seed_val}" in line:
                for token in line.split():
                    if token.startswith("CHI2="):
                        try:
                            chi2 = float(token.split("=")[1])
                        except ValueError:
                            pass
                    elif token.startswith("VALID="):
                        is_valid = token.split("=")[1] == "1"
                    elif token.startswith("STATUS="):
                        try:
                            status = int(token.split("=")[1])
                        except ValueError:
                            pass
                    elif token.startswith("OUT="):
                        out_root = token.split("=")[1]

        return {
            "seed": seed_val,
            "chi2": chi2,
            "is_valid": is_valid,
            "status": status,
            "ret_code": ret_code,
            "root_file": out_root,
            "log_file": str(log_f),
        }

    start_time = time.time()
    completed = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        future_map = {executor.submit(run_seed_task, s): s for s in seeds}
        for future in concurrent.futures.as_completed(future_map):
            res = future.result()
            results.append(res)
            completed += 1
            status_tag = "CONVERGED" if res["is_valid"] and res["status"] == 0 else "FAIL"
            prefix = f"{tag} " if tag else ""
            print(
                f"{prefix}[{completed:3d}/{starts:3d}] Seed {res['seed']:4d} | "
                f"Chi2: {res['chi2']:10.2f} | Status: {status_tag} ({res['status']})"
            )

    elapsed = time.time() - start_time
    valid_results = [r for r in results if r["is_valid"] and r["status"] == 0 and r["chi2"] < 1e8]
    valid_results.sort(key=lambda x: x["chi2"])

    summary_tag = f" {tag}" if tag else ""
    print("\n-------------------------------------------------------")
    print(f" Multi-Start Summary{summary_tag}: {len(valid_results)}/{starts} converged in {elapsed:.1f}s")

    if not valid_results:
        print(f"[ERROR]{summary_tag} No fits converged cleanly. Check logs in:", rundir.fit_starts_dir)
        for r in results:
            if r.get("log_file") and Path(r["log_file"]).is_file():
                try:
                    with open(r["log_file"], "r") as lf:
                        log_text = lf.read()
                    for line in log_text.splitlines():
                        if "*** Break ***" in line or "illegal instruction" in line or "segmentation violation" in line or "bus error" in line or "fatal error" in line:
                            print(f"[ERROR]{summary_tag} Process crashed with fatal error: {line.strip()}")
                            break
                    break
                except Exception:
                    pass
        return {"converged": False, "results": results, "rundir": rundir}

    min_chi2 = valid_results[0]["chi2"]
    best_seed = valid_results[0]["seed"]
    best_file = valid_results[0]["root_file"]

    plateau_thresh = min_chi2 + config.multi_start.plateau_delta_chi2
    plateau_fits = [r for r in valid_results if r["chi2"] <= plateau_thresh]

    print(f" Best Minimum Chi2:    {min_chi2:.4f} (Seed {best_seed})")
    print(f" Plateau Count (<= {plateau_thresh:.4f}): {len(plateau_fits)}/{len(valid_results)}")
    print("-------------------------------------------------------\n")

    # Save nominal best root file
    best_target = rundir.save_nominal_best(best_file)
    # Also write central_best_params.root in root dir for backward compatibility
    legacy_best = rundir.root / "central_best_params.root"
    if Path(best_file).is_file():
        import shutil
        shutil.copyfile(best_file, legacy_best)
        print(f"[OK]{summary_tag} Saved nominal best parameters to: {best_target}")

    # Plot multi-start landscape using plotter
    try:
        plot_multistart_landscape(
            valid_results,
            output_dir=rundir.plots_dir,
            project_name=config.project.name,
            plateau_delta=config.multi_start.plateau_delta_chi2,
        )
    except Exception as e:
        print(f"[WARNING]{summary_tag} Could not generate multi-start plot: {e}")

    # Generate summary report tables
    summary = FitSummary(config, rundir)
    summary.chi2 = min_chi2
    summary.status = valid_results[0]["status"]
    summary.is_valid = valid_results[0]["is_valid"]
    summary.best_seed = best_seed
    summary.total_starts = starts
    summary.converged_starts = len(valid_results)
    summary.plateau_count = len(plateau_fits)
    summary.elapsed_seconds = elapsed

    summary.load_from_root(best_target)
    if not summary.parameters and valid_results[0].get("log_file"):
        summary.load_from_log(valid_results[0]["log_file"])
    txt_rep, md_rep = summary.save_reports()
    print(f"\n{summary.generate_ascii_table()}\n")

    return {
        "converged": True,
        "best_seed": best_seed,
        "best_chi2": min_chi2,
        "plateau_count": len(plateau_fits),
        "best_file": str(best_target),
        "summary_file": str(txt_rep),
        "markdown_report": str(md_rep),
        "valid_count": len(valid_results),
        "rundir": rundir,
    }


def run_multi_start_tbins(
    config: FitConfig,
    tbins: Optional[List[str]] = None,
    n_starts: Optional[int] = None,
    n_threads: Optional[int] = None,
    n_tbin_workers: Optional[int] = None,
    output_dir: Optional[str] = None,
    config_source_path: Optional[Path | str] = None,
) -> Dict[str, Any]:
    """
    Run multi-start fits across one or more t-bins.
    When multiple t-bins are provided, executes them in parallel across concurrent worker threads.
    """
    import copy
    import concurrent.futures

    tbins_list = tbins if tbins is not None else config.dataset.get_t_bins()

    # Single t-bin fallback
    if len(tbins_list) <= 1:
        single_cfg = copy.deepcopy(config)
        if tbins_list:
            single_cfg.dataset.t_bin = tbins_list[0]
            if output_dir:
                resolved_out = output_dir.replace("{t_bin}", tbins_list[0]).replace("{tbin}", tbins_list[0])
            else:
                resolved_out = None
        else:
            resolved_out = output_dir
        return run_multi_start(
            single_cfg,
            n_starts=n_starts,
            n_threads=n_threads,
            output_dir=resolved_out,
            config_source_path=config_source_path,
        )

    # Ensure C++ ACLiC shared libraries are compiled before multi-tbin workers start
    ensure_aclic_updated(config)

    # Multi-t-bin parallel execution
    half_cpus = get_default_max_workers()
    max_tbin_workers = (
        n_tbin_workers
        or getattr(config.multi_start, "n_tbin_workers", None)
        or min(len(tbins_list), half_cpus)
    )
    max_tbin_workers = max(1, min(len(tbins_list), max_tbin_workers, half_cpus))

    starts = n_starts if n_starts is not None else config.multi_start.n_starts
    threads = n_threads if n_threads is not None else (config.multi_start.n_threads or half_cpus)
    if n_threads is None:
        threads = min(threads, max(1, half_cpus // max_tbin_workers))

    print(f"\n=======================================================")
    print(f" Moments2Amplitudes Parallel t-Bins Execution")
    print(f" Project:           {config.project.name}")
    print(f" t-Bins ({len(tbins_list)}):     {', '.join(tbins_list)}")
    print(f" Starts per t-bin:  {starts}")
    print(f" Threads per fit:   {threads}")
    print(f" Concurrent t-Bins: {max_tbin_workers}")
    print(f"=======================================================\n")

    def run_single_tbin(t_name: str) -> Dict[str, Any]:
        cfg_t = copy.deepcopy(config)
        cfg_t.dataset.t_bin = t_name

        # Propagate t_bin to wave initialization paths
        repo_root = get_repo_root()
        for w in cfg_t.waves:
            if w.init.strategy == InitStrategy.FROM_FILE and w.init.file_path:
                resolved_f = resolve_init_file_for_tbin(w.init.file_path, t_name, repo_root)
                w.init.file_path = str(resolved_f)

        if cfg_t.mcmc_propagation.mcmc_tree_path:
            cfg_t.mcmc_propagation.mcmc_tree_path = resolve_mcmc_tree_path(
                cfg_t.mcmc_propagation.mcmc_tree_path, t_name
            )

        if output_dir:
            if "{t_bin}" in output_dir or "{tbin}" in output_dir:
                t_out = output_dir.replace("{t_bin}", t_name).replace("{tbin}", t_name)
            else:
                t_out = str(Path(output_dir) / t_name)
        elif cfg_t.multi_start.output_dir and cfg_t.multi_start.output_dir not in ("logs", "run/central_fits", ""):
            raw_out = cfg_t.multi_start.output_dir
            if "{t_bin}" in raw_out or "{tbin}" in raw_out:
                t_out = raw_out.replace("{t_bin}", t_name).replace("{tbin}", t_name)
            else:
                t_out = str(Path(raw_out) / t_name)
        else:
            t_out = str(Path("runs") / f"{cfg_t.project.name}_{t_name}")

        res = run_multi_start(
            cfg_t,
            n_starts=n_starts,
            n_threads=threads,
            output_dir=t_out,
            config_source_path=config_source_path,
            tag=f"[{t_name}]",
        )
        res["t_bin"] = t_name
        return res

    results_by_tbin: Dict[str, Dict[str, Any]] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_tbin_workers) as executor:
        future_to_tbin = {executor.submit(run_single_tbin, t): t for t in tbins_list}
        for future in concurrent.futures.as_completed(future_to_tbin):
            t_name = future_to_tbin[future]
            try:
                res = future.result()
                results_by_tbin[t_name] = res
            except Exception as e:
                print(f"[ERROR] Failed fit for t-bin {t_name}: {e}")
                results_by_tbin[t_name] = {"t_bin": t_name, "converged": False, "error": str(e)}

    # Consolidated summary table across t-bins
    print(f"\n=========================================================================================")
    print(f" Parallel t-Bins Execution Summary ({len(tbins_list)} t-bins)")
    print(f"=========================================================================================")
    print(f" {'t-bin':<10} {'Chi2':<12} {'Best Seed':<12} {'Plateau':<14} {'Status':<12} {'Run Directory'}")
    print(f" -----------------------------------------------------------------------------------------")
    all_converged = True
    for t_name in tbins_list:
        res = results_by_tbin.get(t_name, {})
        conv = res.get("converged", False)
        if not conv:
            all_converged = False
        status_str = "CONVERGED" if conv else "FAILED"
        rundir_obj = res.get("rundir")
        r_dir_str = str(rundir_obj.root) if rundir_obj else "-"
        best_seed_str = str(res.get("best_seed", "-"))
        chi2_val = res.get("best_chi2")
        chi2_str = f"{chi2_val:.4f}" if chi2_val is not None else "-"
        plateau_count = res.get("plateau_count", "-")
        valid_count = res.get("valid_count", "-")
        plat_str = f"{plateau_count}/{valid_count}" if plateau_count != "-" else "-"
        print(f" {t_name:<10} {chi2_str:<12} {best_seed_str:<12} {plat_str:<14} {status_str:<12} {r_dir_str}")
    print(f"=========================================================================================\n")

    return {
        "converged": all_converged,
        "tbins": results_by_tbin,
    }


# Process-level lock to ensure thread-safe plotting when multiple t-bins complete propagation concurrently
_mcmc_plot_lock = threading.Lock()


def run_mcmc_propagation(
    config: FitConfig,
    best_params_file: Optional[str] = None,
    n_workers: Optional[int] = None,
    max_entries: Optional[int] = None,
    output_file: Optional[str] = None,
    run_dir: Optional[str] = None,
    tag: str = "",
) -> Dict[str, Any]:
    """
    Run parallel MCMC posterior uncertainty propagation across workers and merge with hadd.
    """
    rundir = get_run_directory(config, run_dir)
    rundir.ensure_directories()

    # Ensure C++ ACLiC shared libraries are compiled before MCMC workers start
    ensure_aclic_updated(config)

    half_cpus = get_default_max_workers()
    workers = n_workers if n_workers is not None else (config.mcmc_propagation.n_workers or half_cpus)
    if n_workers is None:
        workers = min(workers, half_cpus)
    total_entries = max_entries if max_entries is not None else config.mcmc_propagation.max_entries
    burn_in = config.mcmc_propagation.burn_in
    out_file = Path(output_file or rundir.mcmc_propagated_root).resolve()

    # Determine best parameters file
    repo_root = get_repo_root()
    if not best_params_file:
        if rundir.nominal_best_root.is_file():
            best_params_file = str(rundir.nominal_best_root)
        elif (rundir.root / "central_best_params.root").is_file():
            best_params_file = str(rundir.root / "central_best_params.root")
        elif config.multi_start.output_dir:
            cand_p = resolve_init_file_for_tbin(config.multi_start.output_dir, config.dataset.t_bin, repo_root)
            if cand_p.is_file():
                best_params_file = str(cand_p)
        if not best_params_file:
            raise FileNotFoundError(f"Best parameters file not found in {rundir.root}. Run multi-start fit first or pass --best-params.")
    else:
        resolved_best = resolve_init_file_for_tbin(best_params_file, config.dataset.t_bin, repo_root)
        if resolved_best.is_file():
            best_params_file = str(resolved_best)

    macro_path = rundir.worker_macro
    generate_root_macro(config, macro_path)

    tag_str = f"{tag} " if tag else ""
    print(f"\n=======================================================")
    print(f" {tag_str}Moments2Amplitudes MCMC Propagation")
    print(f" Project:      {config.project.name}")
    print(f" Total Entries:{total_entries} (burn-in: {burn_in})")
    print(f" Workers:      {workers}")
    print(f" Best Params:  {best_params_file}")
    print(f" Output File:  {out_file}")
    print(f"=======================================================\n")

    # Compute entry chunks per worker
    chunk_size = total_entries // workers
    worker_tasks: List[Tuple[int, int, int]] = []  # (worker_id, start_entry, end_entry)
    for w in range(workers):
        start = burn_in + w * chunk_size
        end = burn_in + total_entries if w == workers - 1 else start + chunk_size
        worker_tasks.append((w, start, end))

    worker_outputs: List[str] = []
    import concurrent.futures

    def run_worker_task(w_info: Tuple[int, int, int]) -> Dict[str, Any]:
        wid, wstart, wend = w_info
        log_f = rundir.get_mcmc_worker_log_file(wid)
        ret_code, output = execute_root_command(
            macro_path,
            ("mcmc", wid, wstart, wend, str(rundir.mcmc_workers_dir), str(best_params_file)),
            log_f,
            config,
        )
        out_root = str(rundir.mcmc_workers_dir / f"mcmc_worker_{wid}.root")
        return {"worker_id": wid, "ret_code": ret_code, "output_file": out_root, "entries": wend - wstart}

    start_time = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(run_worker_task, task) for task in worker_tasks]
        for f in concurrent.futures.as_completed(futures):
            res = f.result()
            print(f"{tag}[Worker {res['worker_id']}] Completed {res['entries']} entries (exit {res['ret_code']})")
            if Path(res["output_file"]).is_file():
                worker_outputs.append(res["output_file"])

    elapsed = time.time() - start_time
    print(f"\n{tag}[OK] All MCMC workers finished in {elapsed:.1f}s.")

    # Merge worker outputs using hadd
    if len(worker_outputs) == 0:
        raise RuntimeError("No MCMC worker root files were generated.")

    print(f"{tag}[Merge] Merging {len(worker_outputs)} root files with hadd -> {out_file}...")
    rootsys_val = resolve_rootsys(config)
    hadd_bin = str(Path(rootsys_val) / "bin" / "hadd") if (Path(rootsys_val) / "bin" / "hadd").is_file() else "hadd"
    hadd_cmd = [hadd_bin, "-f", str(out_file)] + worker_outputs
    hadd_env = get_root_subprocess_env(config)
    hadd_res = subprocess.run(hadd_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=hadd_env)
    if hadd_res.returncode != 0:
        print("[WARNING] hadd output:\n", hadd_res.stderr)
        raise RuntimeError(f"hadd failed with exit code {hadd_res.returncode}")

    print(f"{tag}[SUCCESS] Successfully created merged MCMC propagation results: {out_file}")

    # Generate full diagnostics and plots with MCMC error bands
    print(f"{tag}[Plot] Generating amplitude diagnostic plots and MCMC posterior bands...")
    try:
        with _mcmc_plot_lock:
            diag_res = run_all_diagnostics_and_plots(config, rundir)
        print(f"{tag}[OK] Generated amplitude plots and updated diagnostics in: {rundir.plots_dir}")
    except Exception as e:
        print(f"{tag}[WARNING] Diagnostic plot generation encountered: {e}")

    return {"output_file": str(out_file), "merged_workers": len(worker_outputs), "rundir": rundir}


def run_mcmc_propagation_tbins(
    config: FitConfig,
    tbins: Optional[List[str]] = None,
    best_params_file: Optional[str] = None,
    n_workers: Optional[int] = None,
    max_entries: Optional[int] = None,
    output_file: Optional[str] = None,
    run_dir: Optional[str] = None,
    n_tbin_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Run parallel MCMC posterior uncertainty propagation concurrently across multiple t-bins.
    """
    import concurrent.futures

    tbins_list = tbins if tbins is not None else config.dataset.get_t_bins()
    if len(tbins_list) <= 1:
        single_cfg = copy.deepcopy(config)
        if tbins_list:
            single_cfg.dataset.t_bin = tbins_list[0]
        return run_mcmc_propagation(
            single_cfg,
            best_params_file=best_params_file,
            n_workers=n_workers,
            max_entries=max_entries,
            output_file=output_file,
            run_dir=run_dir,
        )

    # Ensure C++ ACLiC shared libraries are compiled before MCMC workers start
    ensure_aclic_updated(config)

    half_cpus = get_default_max_workers()
    max_tbin_workers = (
        n_tbin_workers
        or getattr(config.mcmc_propagation, "n_tbin_workers", None)
        or getattr(config.multi_start, "n_tbin_workers", None)
        or min(len(tbins_list), half_cpus)
    )
    max_tbin_workers = max(1, min(len(tbins_list), max_tbin_workers, half_cpus))

    effective_workers = n_workers
    if effective_workers is None:
        cfg_workers = getattr(config.mcmc_propagation, "n_workers", None) or half_cpus
        effective_workers = min(cfg_workers, max(1, half_cpus // max_tbin_workers))

    print(f"\n=======================================================")
    print(f" Moments2Amplitudes Parallel MCMC Propagation Across t-Bins")
    print(f" Project:           {config.project.name}")
    print(f" t-Bins ({len(tbins_list)}):     {', '.join(tbins_list)}")
    print(f" Concurrent t-Bins: {max_tbin_workers}")
    print(f"=======================================================\n")

    repo_root = get_repo_root()

    def run_single_tbin_mcmc(t_name: str) -> Dict[str, Any]:
        cfg_t = copy.deepcopy(config)
        cfg_t.dataset.t_bin = t_name
        if cfg_t.mcmc_propagation.mcmc_tree_path:
            cfg_t.mcmc_propagation.mcmc_tree_path = resolve_mcmc_tree_path(
                cfg_t.mcmc_propagation.mcmc_tree_path, t_name
            )

        t_rundir = get_run_directory(cfg_t, run_dir)
        t_best_params = None
        if best_params_file:
            t_best_params = str(resolve_init_file_for_tbin(best_params_file, t_name, repo_root))
        elif t_rundir.nominal_best_root.is_file():
            t_best_params = str(t_rundir.nominal_best_root)
        elif (t_rundir.root / "central_best_params.root").is_file():
            t_best_params = str(t_rundir.root / "central_best_params.root")

        t_outfile = None
        if output_file:
            t_outfile = output_file.replace("{t_bin}", t_name).replace("{tbin}", t_name)

        res = run_mcmc_propagation(
            cfg_t,
            best_params_file=t_best_params,
            n_workers=effective_workers,
            max_entries=max_entries,
            output_file=t_outfile,
            run_dir=str(t_rundir.root),
            tag=f"[{t_name}]",
        )
        res["t_bin"] = t_name
        return res

    results_by_tbin: Dict[str, Dict[str, Any]] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_tbin_workers) as executor:
        future_to_tbin = {executor.submit(run_single_tbin_mcmc, t): t for t in tbins_list}
        for future in concurrent.futures.as_completed(future_to_tbin):
            t_name = future_to_tbin[future]
            try:
                res = future.result()
                results_by_tbin[t_name] = res
            except Exception as e:
                print(f"[ERROR] [{t_name}] MCMC propagation failed: {e}")
                results_by_tbin[t_name] = {"t_bin": t_name, "error": str(e)}

    # Summary across t-bins
    print(f"\n=========================================================================================")
    print(f" Parallel MCMC Propagation Summary ({len(tbins_list)} t-bins)")
    print(f"=========================================================================================")
    print(f" {'t-bin':<10} {'Status':<12} {'MCMC Output File'}")
    print(f" -----------------------------------------------------------------------------------------")
    all_success = True
    for t_name in tbins_list:
        res = results_by_tbin.get(t_name, {})
        succ = "output_file" in res
        if not succ:
            all_success = False
        status_str = "SUCCESS" if succ else "FAILED"
        out_f = res.get("output_file", res.get("error", "-"))
        print(f" {t_name:<10} {status_str:<12} {out_f}")
    print(f"=========================================================================================\n")

    return {
        "success": all_success,
        "tbins": results_by_tbin,
    }



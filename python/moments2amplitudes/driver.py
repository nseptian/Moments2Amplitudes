"""
Execution driver and ROOT macro generator for Moments2Amplitudes fits.
Translates declarative YAML FitConfig into C++ engine execution, orchestrates
parallel multi-start fits, and propagates MCMC posterior uncertainties.
"""

from __future__ import annotations
import os
import sys
import time
import subprocess
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


def translate_formula_to_cpp(formula: str, parameters: List[str]) -> Tuple[List[str], str]:
    """
    Translate a mathematical formula into C++ statements and return expression.
    """
    param_lines = ["const double s = mass * mass;"]
    for i, p in enumerate(parameters):
        param_lines.append(f"const double {p} = (p.size() > {i}) ? p[{i}] : 0.0;")

    # Check if formula matches Two Breit-Wigner form: k1 / (M1^2 - s - i*M1*w1) + (k2 / (M2^2 - s - i*M2*w2)) * exp(i*phase)
    two_bw_regex = re.compile(
        r'([A-Za-z0-9_]+)\s*/\s*\(\s*([A-Za-z0-9_]+)(?:\^2|\s*\*\s*\2)\s*-\s*s\s*-\s*i\s*\*\s*\2\s*\*\s*([A-Za-z0-9_]+)\s*\)\s*\+\s*'
        r'(?:\(\s*)?([A-Za-z0-9_]+)\s*/\s*\(\s*([A-Za-z0-9_]+)(?:\^2|\s*\*\s*\5)\s*-\s*s\s*-\s*i\s*\*\s*\5\s*\*\s*([A-Za-z0-9_]+)\s*\)(?:\s*\))?'
        r'(?:\s*\*\s*exp\(\s*i\s*\*\s*([A-Za-z0-9_]+)\s*\))?'
    )
    two_bw_m = two_bw_regex.search(formula)
    if two_bw_m:
        k1, m1, w1 = two_bw_m.group(1), two_bw_m.group(2), two_bw_m.group(3)
        k2, m2_var, w2 = two_bw_m.group(4), two_bw_m.group(5), two_bw_m.group(6)
        phase_var = two_bw_m.group(7)
        bw1_str = f"{k1} / std::complex<double>(({m1} * {m1}) - s, -{m1} * {w1})"
        if phase_var:
            bw2_str = f"({k2} / std::complex<double>(({m2_var} * {m2_var}) - s, -{m2_var} * {w2})) * std::polar(1.0, {phase_var})"
        else:
            bw2_str = f"{k2} / std::complex<double>(({m2_var} * {m2_var}) - s, -{m2_var} * {w2})"
        return param_lines, f"{bw1_str} + {bw2_str}"

    # Check if formula matches canonical Breit-Wigner form: k / (M^2 - s - i * M * width)
    bw_regex = re.compile(
        r'([A-Za-z0-9_]+)\s*/\s*\(\s*([A-Za-z0-9_]+)(?:\^2|\s*\*\s*\2)\s*-\s*s\s*-\s*i\s*\*\s*\2\s*\*\s*([A-Za-z0-9_]+)\s*\)'
    )
    bw_m = bw_regex.search(formula)
    if bw_m:
        k_var, m_var, w_var = bw_m.group(1), bw_m.group(2), bw_m.group(3)
        return param_lines, f"{k_var} / std::complex<double>(({m_var} * {m_var}) - s, -{m_var} * {w_var})"

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
        return param_lines, f"{k_param} / std::complex<double>(({mass_param} * {mass_param}) - s, -gamma_total)"

    expr = formula
    # Power replacement: X^2 -> (X * X)
    expr = re.sub(r'([A-Za-z0-9_]+)\^2', r'(\1 * \1)', expr)
    expr = re.sub(r'([A-Za-z0-9_]+)\^([0-9]+)', r'std::pow(\1, \2)', expr)

    # Breakup momentum / phase space factors
    expr = expr.replace("rho_etapi", "m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.547853, 0.13957)")
    expr = expr.replace("rho_KK", "m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.493677, 0.493677)")
    expr = expr.replace("rho_etaprimepi", "m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.95778, 0.13957)")

    # Complex imaginary unit replacements
    expr = re.sub(r'\bi\s*\*\s*', 'std::complex<double>(0.0, 1.0) * ', expr)
    expr = re.sub(r'\bi\b', 'std::complex<double>(0.0, 1.0)', expr)
    expr = re.sub(r'\bexp\(', 'std::exp(', expr)
    expr = re.sub(r'\bsqrt\(', 'std::sqrt(', expr)

    return param_lines, expr


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

    # Determine project directory and brufit path
    repo_root = Path(__file__).resolve().parent.parent.parent
    brufit_path = config.environment.brufit_dir or os.environ.get("BRUFIT", "/d/grid17/septian/brufit_dev/brufit")

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
        '#include "MomentHelper.h"',
        '#include "MassDependentMoments.h"',
        '#include "MassDependentFitter.h"',
        "",
        "using namespace m2pw;",
        "using namespace std;",
        "",
        "void EnsureLoaded() {",
        f'    gSystem->Setenv("MOMS2AMPS", "{repo_root}");',
        f'    gSystem->Setenv("BRUFIT", "{brufit_path}");',
        f'    TString loadScript = Form("{repo_root}/macros/Load.C");',
        '    if (gSystem->AccessPathName(loadScript) == 0) {',
        '        gROOT->ProcessLine(Form(".x %s", loadScript.Data()));',
        "    }",
        "}",
        "",
        "void RunWorker(const char* mode = \"single\", int seed = 1, int mcmcStart = 0, int mcmcEnd = -1, const char* outDir = \".\", const char* bestParamsFile = \"\") {",
        "    EnsureLoaded();",
        "    TStopwatch timer;",
        "    timer.Start();",
        "",
        "    // Setup mass bins",
        f"    const std::vector<double> massBins = {{{mass_bins_str}}};",
        f"    const double lMax = {config.chi2.max_l};",
        "    const double noise = 0.0;",
        "    const std::vector<TString> noUse = {};",
        "",
        "    // Setup MassDependentMoments",
        f'    const TString momentsDir = "{config.dataset.moments_dir}";',
        f'    const TString resultsFile = "{config.dataset.results_filename}";',
        f'    const TString binPrefix = "{config.dataset.bin_prefix}";',
        f'    const TString varBinPrefix = "{config.dataset.var_bin_prefix}";',
        "",
        "    std::cout << \"[Worker] Loading dataset moments across \" << massBins.size() << \" mass bins...\" << std::endl;",
        "    MassDependentMoments moms(massBins, momentsDir, resultsFile, binPrefix, varBinPrefix);",
        "",
        "    // Configure mass-dependent models and moments selection",
        "    MassDependentFitter::MassDependenceConfig massDepConfig;",
        f"    massDepConfig.massIndependentL = {{{mi_l_str}}};",
    ]

    # Mass-dependent wave models
    for w in md_waves:
        for res in w.resonances:
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
            lines.append(f'    massDepConfig.massDependentWaves[{w.l}].push_back({{"{res.name}", {m_type}, {w.poly_order}, "{form_str}", {{{param_list_str}}}}});')

    # Moments selection
    lines.extend([
        "",
        "    MassDependentFitter::MomentsConfig momentsConfig;",
        f"    momentsConfig.includeAll = {str(config.chi2.include_all).lower()};",
        f"    momentsConfig.includedL = {{{inc_l_str}}};",
        f"    momentsConfig.excludedL = {{{exc_l_str}}};",
        "",
        "    // Create fitter instance",
        "    Setup fitSetup;",
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
    for w in md_waves:
        for res in w.resonances:
            if w.formula:
                param_lines, expr = translate_formula_to_cpp(w.formula, w.parameters)
                param_code = "\n        ".join(param_lines)
                form_esc = w.formula.replace('"', '\\"')
                lines.extend([
                    f'    // Explicit formula configured for L={w.l} {res.name}: {w.formula}',
                    f'    fitter.SetWaveCustomFormula({w.l}, "{res.name}", [](double mass, const std::vector<double>& p) -> std::complex<double> {{',
                    f'        {param_code}',
                    f'        return {expr};',
                    f'    }}, "{form_esc}");',
                ])

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
    if mi_l_list:
        lines.append(f"    paramManager.AddMassIndependentParameter({{{mi_l_str}}}, seed);")

    # Add mass-dependent parameters
    for w in md_waves:
        if w.init.strategy == InitStrategy.FROM_FILE and w.init.file_path:
            lines.append(f'    paramManager.AddMassDependentParameter({{{w.l}}}, "{w.init.file_path}", {w.init.option});')
        else:
            lines.append(f"    paramManager.AddMassDependentParameter({{{w.l}}}, seed, {w.init.option});")

        # Set resonance parameters defaults/limits/fixes
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
        "void generated_driver(const char* mode = \"single\", int seed = 1, int mcmcStart = 0, int mcmcEnd = -1, const char* outDir = \".\", const char* bestParamsFile = \"\") {",
        "    RunWorker(mode, seed, mcmcStart, mcmcEnd, outDir, bestParamsFile);",
        "}",
    ])

    with open(output_macro_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return output_macro_path


def execute_root_command(macro_path: Path, args: Tuple, log_path: Path) -> Tuple[int, str]:
    """Execute a single ROOT invocation in batch mode and capture output."""
    arg_strs = []
    for a in args:
        if isinstance(a, str):
            arg_strs.append(f'"{a}"')
        else:
            arg_strs.append(str(a))
    call_expr = f'{macro_path.name}({", ".join(arg_strs)})'

    cmd = [
        "root",
        "-l",
        "-b",
        "-q",
        f"{macro_path}{call_expr[len(macro_path.name):]}",
    ]

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log_f:
        proc = subprocess.run(
            cmd,
            stdout=log_f,
            stderr=subprocess.STDOUT,
            cwd=str(macro_path.parent),
            text=True,
        )

    # Read log
    log_content = ""
    if log_path.is_file():
        with open(log_path, "r") as f:
            log_content = f.read()

    return proc.returncode, log_content


def run_multi_start(
    config: FitConfig,
    n_starts: Optional[int] = None,
    n_threads: Optional[int] = None,
    output_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run parallel multi-start randomized fits to find the central global minimum.
    Generates summary plots, identifies the plateau, and saves central_best_params.root.
    """
    starts = n_starts if n_starts is not None else config.multi_start.n_starts
    threads = n_threads if n_threads is not None else config.multi_start.n_threads
    out_dir = Path(output_dir or config.multi_start.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    central_fits_dir = out_dir / "central_fits"
    central_fits_dir.mkdir(parents=True, exist_ok=True)

    macro_path = out_dir / "generated_driver.C"
    generate_root_macro(config, macro_path)

    print(f"\n=======================================================")
    print(f" Moments2Amplitudes Multi-Start Fit")
    print(f" Project:      {config.project.name}")
    print(f" Total Trials: {starts}")
    print(f" Threads:      {threads}")
    print(f" Output Dir:   {out_dir}")
    print(f"=======================================================\n")

    seeds = [config.multi_start.seed_offset + i for i in range(starts)]
    results: List[Dict[str, Any]] = []

    # Parallel execution using concurrent.futures
    import concurrent.futures

    def run_seed_task(seed_val: int) -> Dict[str, Any]:
        log_f = central_fits_dir / f"fit_seed_{seed_val}.log"
        ret_code, output = execute_root_command(
            macro_path,
            ("single", seed_val, 0, -1, str(central_fits_dir), ""),
            log_f,
        )

        # Parse [RESULT] line
        chi2 = 1e9
        is_valid = False
        status = -1
        out_root = str(central_fits_dir / f"MassDependentResults_{seed_val}.root")

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
            print(
                f"[{completed:3d}/{starts:3d}] Seed {res['seed']:4d} | "
                f"Chi2: {res['chi2']:10.2f} | Status: {status_tag} ({res['status']})"
            )

    elapsed = time.time() - start_time
    valid_results = [r for r in results if r["is_valid"] and r["status"] == 0 and r["chi2"] < 1e8]
    valid_results.sort(key=lambda x: x["chi2"])

    print("\n-------------------------------------------------------")
    print(f" Multi-Start Summary: {len(valid_results)}/{starts} converged in {elapsed:.1f}s")

    if not valid_results:
        print("[ERROR] No fits converged cleanly. Check logs in:", central_fits_dir)
        return {"converged": False, "results": results}

    min_chi2 = valid_results[0]["chi2"]
    best_seed = valid_results[0]["seed"]
    best_file = valid_results[0]["root_file"]

    plateau_thresh = min_chi2 + config.multi_start.plateau_delta_chi2
    plateau_fits = [r for r in valid_results if r["chi2"] <= plateau_thresh]

    print(f" Best Minimum Chi2:    {min_chi2:.4f} (Seed {best_seed})")
    print(f" Plateau Count (<= {plateau_thresh:.4f}): {len(plateau_fits)}/{len(valid_results)}")
    print("-------------------------------------------------------\n")

    # Copy best params
    best_target = out_dir / "central_best_params.root"
    if Path(best_file).is_file():
        import shutil
        shutil.copyfile(best_file, best_target)
        print(f"[OK] Saved central best parameters to: {best_target}")

    # Write summary text file
    summary_path = out_dir / "central_fits_summary.txt"
    with open(summary_path, "w") as sf:
        sf.write(f"Project: {config.project.name}\n")
        sf.write(f"Total Trials: {starts}\n")
        sf.write(f"Converged: {len(valid_results)}\n")
        sf.write(f"Best Chi2: {min_chi2:.4f} (Seed {best_seed})\n")
        sf.write(f"Plateau Delta Chi2: {config.multi_start.plateau_delta_chi2}\n")
        sf.write(f"Plateau Fits Count: {len(plateau_fits)}\n\n")
        sf.write("Rank | Seed | Chi2 | Valid | Status | RootFile\n")
        for i, r in enumerate(valid_results):
            sf.write(f"{i+1:4d} | {r['seed']:4d} | {r['chi2']:10.4f} | {r['is_valid']} | {r['status']:2d} | {r['root_file']}\n")

    # Plot Chi2 distribution using matplotlib
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        chi2_vals = [r["chi2"] for r in valid_results]
        fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
        ax.hist(chi2_vals, bins=30, color="#2b5c8f", edgecolor="black", alpha=0.8)
        ax.axvline(min_chi2, color="red", linestyle="--", linewidth=2, label=f"Best $\\chi^2 = {min_chi2:.2f}$")
        ax.axvline(plateau_thresh, color="orange", linestyle=":", linewidth=2, label=f"Plateau threshold ($\\Delta \\chi^2 \\le {config.multi_start.plateau_delta_chi2}$)")
        ax.set_title(f"{config.project.name} - $\\chi^2$ Distribution ({len(valid_results)} converged)", fontsize=13)
        ax.set_xlabel("$\\chi^2$", fontsize=11)
        ax.set_ylabel("Count", fontsize=11)
        ax.grid(True, linestyle="--", alpha=0.5)
        ax.legend(fontsize=10)

        dist_png = out_dir / "central_chi2_distribution.png"
        dist_pdf = out_dir / "central_chi2_distribution.pdf"
        fig.tight_layout()
        fig.savefig(dist_png)
        fig.savefig(dist_pdf)
        plt.close(fig)
        print(f"[OK] Generated Chi2 distribution plot: {dist_png}")
    except Exception as e:
        print(f"[WARNING] Could not generate matplotlib plot: {e}")

    return {
        "converged": True,
        "best_seed": best_seed,
        "best_chi2": min_chi2,
        "plateau_count": len(plateau_fits),
        "best_file": str(best_target),
        "summary_file": str(summary_path),
        "valid_count": len(valid_results),
    }


def run_mcmc_propagation(
    config: FitConfig,
    best_params_file: Optional[str] = None,
    n_workers: Optional[int] = None,
    max_entries: Optional[int] = None,
    output_file: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run parallel MCMC posterior uncertainty propagation across workers and merge with hadd.
    """
    workers = n_workers if n_workers is not None else config.mcmc_propagation.n_workers
    total_entries = max_entries if max_entries is not None else config.mcmc_propagation.max_entries
    burn_in = config.mcmc_propagation.burn_in
    out_file = Path(output_file or config.mcmc_propagation.output_file).resolve()
    log_dir = Path(config.mcmc_propagation.log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)

    # Determine best parameters file
    if not best_params_file:
        candidate = Path(config.multi_start.output_dir) / "central_best_params.root"
        if candidate.is_file():
            best_params_file = str(candidate)
        else:
            raise FileNotFoundError("Best parameters file not found. Run multi-start fit first or pass --best-params-file.")

    macro_path = log_dir / "generated_driver.C"
    generate_root_macro(config, macro_path)

    print(f"\n=======================================================")
    print(f" Moments2Amplitudes MCMC Propagation")
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
        log_f = log_dir / f"mcmc_worker_{wid}.log"
        ret_code, output = execute_root_command(
            macro_path,
            ("mcmc", wid, wstart, wend, str(log_dir), str(best_params_file)),
            log_f,
        )
        out_root = str(log_dir / f"mcmc_worker_{wid}.root")
        return {"worker_id": wid, "ret_code": ret_code, "output_file": out_root, "entries": wend - wstart}

    start_time = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(run_worker_task, task) for task in worker_tasks]
        for f in concurrent.futures.as_completed(futures):
            res = f.result()
            print(f"[Worker {res['worker_id']}] Completed {res['entries']} entries (exit {res['ret_code']})")
            if Path(res["output_file"]).is_file():
                worker_outputs.append(res["output_file"])

    elapsed = time.time() - start_time
    print(f"\n[OK] All MCMC workers finished in {elapsed:.1f}s.")

    # Merge worker outputs using hadd
    if len(worker_outputs) == 0:
        raise RuntimeError("No MCMC worker root files were generated.")

    print(f"[Merge] Merging {len(worker_outputs)} root files with hadd -> {out_file}...")
    hadd_cmd = ["hadd", "-f", str(out_file)] + worker_outputs
    hadd_res = subprocess.run(hadd_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if hadd_res.returncode != 0:
        print("[WARNING] hadd output:\n", hadd_res.stderr)
        raise RuntimeError(f"hadd failed with exit code {hadd_res.returncode}")

    print(f"[SUCCESS] Successfully created merged MCMC propagation results: {out_file}")
    return {"output_file": str(out_file), "merged_workers": len(worker_outputs)}


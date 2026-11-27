"""
Stage 1: Linear algebra moments calculation using Moments.Main.
"""

from __future__ import annotations

import functools
import multiprocessing
import os
import shutil
import subprocess
import sys
import time
from copy import deepcopy
from pathlib import Path
from typing import Any, List, Optional, Tuple

import threadpoolctl
from wurlitzer import STDOUT, pipes

from .config import PipelineConfig

# Ensure immediate print flushing
print = functools.partial(print, flush=True)


def _get_py39_env(mom2amp_dir: str) -> dict:
    """Constructs a clean isolated environment for Python 3.9 + ROOT 6.28."""
    repo_root = Path(mom2amp_dir)
    py39_env_dir = repo_root / ".pixi" / "envs" / "py39"
    if not (py39_env_dir / "bin" / "python").is_file():
        py39_env_dir = Path("/d/grid17/septian/Moments-dev/Moments.Nizar/.pixi/envs/default")

    clean_pypath = [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p and not any(k in p for k in ("root_v6", "cuda_install", "cuda_build", "brufit"))]
    clean_ldpath = [p for p in os.environ.get("LD_LIBRARY_PATH", "").split(os.pathsep) if p and not any(k in p for k in ("root_v6", "cuda_install", "cuda_build", "brufit"))]

    env = os.environ.copy()
    for k in list(env.keys()):
        if k.startswith("ROOT"):
            env.pop(k, None)
    env["PATH"] = f"{py39_env_dir}/bin:/usr/bin:/bin"
    env["LD_LIBRARY_PATH"] = f"{py39_env_dir}/lib"
    env["PYTHONPATH"] = f"/d/grid17/septian/Moments-dev/Moments.Main:{repo_root}/python"
    env["PYTHONNOUSERSITE"] = "1"
    env["MPLCONFIGDIR"] = "/tmp"
    return env


def _find_py39_runner(mom2amp_dir: str) -> Tuple[List[str], dict]:
    """
    Locates the Python 3.9 / ROOT 6.28 runner command prefix for linear algebra (Moments.Main).
    """
    repo_root = Path(mom2amp_dir)
    env = _get_py39_env(mom2amp_dir)

    # 1. Direct py39 environment binary if present
    py39_bin = repo_root / ".pixi" / "envs" / "py39" / "bin" / "python"
    if py39_bin.is_file():
        return [str(py39_bin)], env

    # 2. Moments.Nizar fallback
    nizar_bin = Path("/d/grid17/septian/Moments-dev/Moments.Nizar/.pixi/envs/default/bin/python")
    if nizar_bin.is_file():
        return [str(nizar_bin)], env

    # 3. Pixi run
    pixi_bin = shutil.which("pixi") or os.path.expanduser("~/.pixi/bin/pixi")
    manifest_path = repo_root / "pixi.toml"
    if os.path.isfile(pixi_bin) and manifest_path.is_file():
        return [pixi_bin, "run", "--manifest-path", str(manifest_path), "-e", "py39", "python"], env

    raise RuntimeError(
        "Could not locate a valid Python 3.9 + ROOT 6.28 environment for Moments.Main linear algebra.\n"
        "Please ensure '.pixi/envs/py39' or '/d/grid17/septian/Moments-dev/Moments.Nizar' is available."
    )


def _init_root_and_moments_main(moments_main_dir: str) -> None:
    """Initialize ROOT environment and load Moments.Main basis functions library."""
    if moments_main_dir not in sys.path:
        sys.path.insert(0, moments_main_dir)

    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetMacroPath(f"{moments_main_dir}:{ROOT.gROOT.GetMacroPath()}")

    # Ensure pickle compatibility between 'moments.MomentCalculator' and 'MomentCalculator'
    import moments.MomentCalculator
    sys.modules["MomentCalculator"] = moments.MomentCalculator

    from workflow.RootUtilities import loadBasisFunctionsLibrary
    loadBasisFunctionsLibrary()


def run_linear_algebra_worker(
    cfg: PipelineConfig,
    tbin: str,
    max_l: int,
    log_to_file: bool = True,
) -> int:
    """
    Executes moment analysis for a single (t-bin, maxL) combination in a dedicated worker process.
    """
    moments_main_dir = cfg.environment.moments_main_dir
    _init_root_and_moments_main(moments_main_dir)

    import ROOT
    from moments.MomentCalculator import (
        DataSet,
        KinematicBinningVariable,
        MomentCalculator,
        MomentCalculatorsKinematicBinning,
        MomentIndices,
    )
    from workflow.PlottingUtilities import HistAxisBinning
    from workflow import Utilities

    out_dir = cfg.linear_algebra.get_output_dir(tbin, max_l, cfg.analysis.beam_pol_label)
    os.makedirs(out_dir, exist_ok=True)

    log_file_path = os.path.join(out_dir, "calc_moments.log")

    def _execute():
        print(f"=== Starting Linear Algebra Moments Calculation ===")
        print(f"  t-bin:        {tbin}")
        print(f"  Max L:        {max_l}")
        print(f"  Output Dir:   {out_dir}")
        print(f"  Moments.Main: {moments_main_dir}")

        timer = Utilities.Timer()
        timer.start("Total calculation time")

        # Setup thread limit for worker process
        thread_controller = threadpoolctl.ThreadpoolController()
        with thread_controller.limit(limits=cfg.linear_algebra.openmp_threads):
            print(f"OpenMP thread limits: {cfg.linear_algebra.openmp_threads}")

            # 1. Resolve file paths
            data_file = cfg.input_data.resolve_path(cfg.input_data.data_file_pattern, tbin)
            acc_file = cfg.input_data.resolve_path(cfg.input_data.ps_acc_file_pattern, tbin)
            gen_file = cfg.input_data.resolve_path(cfg.input_data.ps_gen_file_pattern, tbin)

            tree_name = cfg.input_data.tree_name

            print(f"Loading data from '{data_file}' (tree: '{tree_name}')")
            df_data = ROOT.RDataFrame(tree_name, data_file)

            print(f"Loading accepted MC from '{acc_file}' (tree: '{tree_name}')")
            df_acc = ROOT.RDataFrame(tree_name, acc_file)

            print(f"Loading generated MC from '{gen_file}' (tree: '{tree_name}')")
            df_gen = ROOT.RDataFrame(tree_name, gen_file)

            # 2. Configure weight branch and required columns for Moments.Main
            weight_branch = getattr(cfg.linear_algebra, "weight_branch", None) or getattr(cfg.input_data, "weight_branch", "weightASBS")
            if not df_data.HasColumn(weight_branch):
                raise ValueError(
                    f"Weight column '{weight_branch}' not found in real-data tree '{tree_name}' of file '{data_file}'. "
                    f"Available columns: {list(df_data.GetColumnNames())}"
                )
            print(f"Using weight branch '{weight_branch}' for linear algebra moment calculation")

            for col_name, col_formula in cfg.input_data.columns.items():
                if not df_data.HasColumn(col_name):
                    print(f"Defining real-data column '{col_name}' = '{col_formula}'")
                    df_data = df_data.Define(col_name, col_formula)
                if not df_acc.HasColumn(col_name):
                    print(f"Defining acc-MC column '{col_name}' = '{col_formula}'")
                    df_acc = df_acc.Define(col_name, col_formula)

            # 3. Apply additional cuts if any
            for cut in cfg.input_data.additional_cuts:
                print(f"Applying cut: {cut}")
                df_data = df_data.Filter(cut)
                df_acc = df_acc.Filter(cut)

            # 4. Configure MomentIndices and Mass binning
            is_polarized = (cfg.analysis.polarization is not None)
            indices_meas = MomentIndices(maxL=max_l, polarized=is_polarized)
            indices_phys = MomentIndices(maxL=max_l, polarized=is_polarized)

            bin_cfg = cfg.analysis.mass_binning
            var_mass = KinematicBinningVariable(
                name=bin_cfg.variable,
                label=bin_cfg.label,
                unit=bin_cfg.unit,
                nmbDigits=bin_cfg.digits,
            )

            moment_calculators = MomentCalculatorsKinematicBinning([])

            print(f"Preparing {bin_cfg.n_bins} mass bins in range [{bin_cfg.min}, {bin_cfg.max}] {bin_cfg.unit}:")

            for bin_idx in range(bin_cfg.n_bins):
                bin_range = bin_cfg.get_bin_range(bin_idx)
                bin_center = bin_cfg.get_bin_center(bin_idx)
                bin_filter = bin_cfg.get_bin_filter(bin_idx)

                df_data_bin = df_data.Filter(bin_filter)
                df_acc_bin = df_acc.Filter(bin_filter)
                df_gen_bin = df_gen.Filter(bin_filter)

                n_data = df_data_bin.Count().GetValue()
                n_acc = df_acc_bin.Count().GetValue()
                n_gen = df_gen_bin.Count().GetValue()

                eff = (n_acc / n_gen) if n_gen > 0 else 0.0
                print(
                    f"  Bin [{bin_idx + 1:2d}/{bin_cfg.n_bins:2d}] center={bin_center:.3f} {bin_cfg.unit} "
                    f"range=[{bin_range[0]:.3f}, {bin_range[1]:.3f}] | Events: data={n_data}, acc={n_acc}, gen={n_gen} (eff={eff:.3f})"
                )

                dataset = DataSet(
                    data=df_data_bin,
                    phaseSpaceData=df_acc_bin,
                    nmbGenEvents=n_gen if n_gen > 0 else n_acc,
                    polarization=cfg.analysis.polarization,
                    weightColumn=weight_branch,
                )

                calc = MomentCalculator(
                    indicesMeas=indices_meas,
                    indicesPhys=indices_phys,
                    dataSet=dataset,
                    binCenters={var_mass: bin_center},
                    integralFileBaseName=f"{out_dir}/integralMatrix",
                    flipSignOfWeights=False,
                )
                moment_calculators.append(calc)

            # 5. Acceptance Integral Matrix Calculation
            with timer.timeThis(f"Integral matrix calculation for {len(moment_calculators)} bins"):
                print(f"Calculating acceptance integral matrices (force={cfg.linear_algebra.force_integral_matrix})...")
                moment_calculators.calculateIntegralMatrices(
                    forceCalculation=cfg.linear_algebra.force_integral_matrix
                )

            # 6. Physical & Measured Moment Calculation
            with timer.timeThis(f"Moments calculation for {len(moment_calculators)} bins"):
                print(f"Calculating moments (normalize={cfg.linear_algebra.normalize_moments}, bootstrap={cfg.linear_algebra.bootstrap_samples})...")
                moment_calculators.calculateMoments(
                    normalize=cfg.linear_algebra.normalize_moments,
                    nmbBootstrapSamples=cfg.linear_algebra.bootstrap_samples,
                )

            # 7. Save Pickled Moment Results
            prefix = "norm" if cfg.linear_algebra.normalize_moments else "unnorm"
            phys_pkl = os.path.join(out_dir, f"{prefix}_moments_phys.pkl")
            meas_pkl = os.path.join(out_dir, f"{prefix}_moments_meas.pkl")

            moment_calculators.momentResultsPhys.savePickle(phys_pkl)
            moment_calculators.momentResultsMeas.savePickle(meas_pkl)

            timer.stop("Total calculation time")
            print(f"Calculation complete! Summary:\n{timer.summary}")
            print(f"Saved physical moments: {phys_pkl}")
            print(f"Saved measured moments: {meas_pkl}")

    if log_to_file:
        print(f"  Logging detailed output to: {log_file_path}")
        with open(log_file_path, "w") as log_file, pipes(stdout=log_file, stderr=STDOUT):
            _execute()
    else:
        _execute()

    return 0


def run_stage1(
    cfg: PipelineConfig,
    tbins_filter: Optional[List[str]] = None,
    max_l_filter: Optional[List[int]] = None,
    dry_run: bool = False,
) -> bool:
    """
    Executes Stage 1 across all configured t-bins and max_l multipoles.
    """
    tbins = tbins_filter if tbins_filter is not None else cfg.analysis.tbins
    max_ls = max_l_filter if max_l_filter is not None else cfg.analysis.max_l

    tasks: List[Tuple[str, int]] = [(tbin, l) for tbin in tbins for l in max_ls]

    print("\n" + "=" * 70)
    print(" STAGE 1: Linear Algebra Moments Calculation (Moments.Main)")
    print("=" * 70)
    print(f" Target t-bins:      {', '.join(tbins)}")
    print(f" Multipoles Max L:   {max_ls}")
    print(f" Total Tasks:        {len(tasks)}")
    print(f" Output Base Dir:    {cfg.linear_algebra.output_base_dir}")
    print(f" Max Parallel Jobs:  {cfg.linear_algebra.max_parallel_workers}")
    print("=" * 70 + "\n")

    if dry_run:
        print("[DRY-RUN] Tasks to execute:")
        for tbin, l in tasks:
            out = cfg.linear_algebra.get_output_dir(tbin, l, cfg.analysis.beam_pol_label)
            print(f"  * tbin={tbin}, maxL={l} -> {out}")
        return True

    is_py312 = sys.version_info >= (3, 10)
    cmd_prefix = None
    py39_env = None
    if is_py312:
        cmd_prefix, py39_env = _find_py39_runner(cfg.environment.mom2amp_dir)
        print(f"  Stage 1 Runner: {' '.join(cmd_prefix)} (Python 3.9 / ROOT 6.28)")

    worker_context = multiprocessing.get_context("fork") if not is_py312 else None
    active_processes: List[Tuple[Any, str, int]] = []
    concurrency_limit = max(1, cfg.linear_algebra.max_parallel_workers)

    failed_tasks: List[Tuple[str, int, int]] = []
    completed_count = 0

    task_queue = list(tasks)

    while task_queue or active_processes:
        # Launch new workers up to concurrency limit
        while task_queue and len(active_processes) < concurrency_limit:
            tbin, l = task_queue.pop(0)
            print(f"[LAUNCH] Worker for tbin={tbin}, maxL={l}")
            if is_py312:
                cmd = cmd_prefix + [
                    "-m", "moments2amplitudes.momentfit.stage1_linear_algebra",
                    "--config", cfg.config_path,
                    "--tbin", tbin,
                    "--max-l", str(l),
                ]
                proc = subprocess.Popen(cmd, env=py39_env)
            else:
                proc = worker_context.Process(
                    target=run_linear_algebra_worker,
                    args=(cfg, tbin, l, True),
                )
                proc.start()
            active_processes.append((proc, tbin, l))

        # Check running workers
        time.sleep(0.5)
        still_running = []
        for proc, tbin, l in active_processes:
            is_alive = proc.poll() is None if is_py312 else proc.is_alive()
            if not is_alive:
                if not is_py312:
                    proc.join()
                    code = proc.exitcode
                else:
                    code = proc.returncode
                completed_count += 1
                if code == 0:
                    print(f"[FINISHED] [{completed_count}/{len(tasks)}] ✅ Task tbin={tbin}, maxL={l} completed successfully.")
                else:
                    print(f"[ERROR] [{completed_count}/{len(tasks)}] ❌ Task tbin={tbin}, maxL={l} failed with exit code {code}!")
                    failed_tasks.append((tbin, l, code))
            else:
                still_running.append((proc, tbin, l))
        active_processes = still_running

    print("\n" + "-" * 70)
    if failed_tasks:
        print(f"❌ Stage 1 completed with {len(failed_tasks)} failed task(s):")
        for tbin, l, code in failed_tasks:
            print(f"  - tbin={tbin}, maxL={l} (exit code: {code})")
        return False
    else:
        print(f"✅ Stage 1 completed successfully! All {len(tasks)} tasks finished.")
        return True


def main() -> None:
    import argparse
    from .config import load_pipeline_config

    parser = argparse.ArgumentParser(description="Stage 1 Linear Algebra Worker (Python 3.9 / ROOT 6.28)")
    parser.add_argument("--config", "-c", required=True, help="Path to YAML configuration")
    parser.add_argument("--tbin", "-t", required=True, help="Target t-bin")
    parser.add_argument("--max-l", "-l", type=int, required=True, help="Target max L")
    parser.add_argument("--no-log", action="store_true", help="Log to stdout instead of log file")
    args = parser.parse_args()

    cfg = load_pipeline_config(args.config)
    ret = run_linear_algebra_worker(cfg, args.tbin, args.max_l, log_to_file=not args.no_log)
    sys.exit(ret)


if __name__ == "__main__":
    main()

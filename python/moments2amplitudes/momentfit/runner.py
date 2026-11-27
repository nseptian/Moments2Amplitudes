"""
MomentFit runner and CLI entry point for Moments analysis workflow.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional

from .config import (
    MomentFitConfig,
    PipelineConfig,
    load_momentfit_config,
    load_pipeline_config,
    validate_momentfit_config,
    validate_pipeline_config,
)


def find_script(script_path: Optional[str], default_rel_path: str, repo_root: Path) -> Optional[str]:
    """Find executable script path by config, relative repo path, or fallback."""
    if script_path and os.path.isfile(script_path):
        return script_path
    local_path = repo_root / default_rel_path
    if local_path.is_file():
        return str(local_path)
    alt_path = Path("/d/home/septian/EtaPi0Analysis") / default_rel_path
    if alt_path.is_file():
        return str(alt_path)
    return None


def run_monitor(
    cfg: Optional[MomentFitConfig] = None,
    run_dir: Optional[str] = None,
    tbins_filter: Optional[List[str]] = None,
    max_l_filter: Optional[List[int]] = None,
    watch: bool = False,
    interval: float = 10.0,
    progress_only: bool = False,
    cancel: bool = False,
    force: bool = False,
    plot: bool = False,
    skip_trace: bool = False,
    skip_posterior: bool = False,
    skip_moments: bool = False,
    downsample_trace: Optional[int] = None,
    replot_from_root: bool = False,
    ess_threshold: Optional[float] = None,
    ess_metric: Optional[str] = None,
    root_filename: Optional[str] = None,
    dry_run: bool = False,
) -> bool:
    """
    Launches MonitorMomentFit.py across target run directories.
    """
    repo_root = Path(__file__).resolve().parent.parent.parent.parent
    configured_script = getattr(cfg.moment_fit, "monitor_script", None) if cfg else None
    script_path = find_script(configured_script, "macros/MonitorMomentFit.py", repo_root)
    if not script_path:
        print("❌ Error: MonitorMomentFit.py script not found.")
        return False

    target_entries: List[tuple[Optional[int], str]] = []
    if run_dir:
        explicit_l = max_l_filter[0] if max_l_filter and len(max_l_filter) == 1 else None
        target_entries.append((explicit_l, run_dir))
    elif cfg:
        max_ls = max_l_filter if max_l_filter is not None else cfg.analysis.max_l
        m_min = cfg.analysis.mass_binning.min
        m_max = cfg.analysis.mass_binning.max
        for l in max_ls:
            target_entries.append((l, cfg.moment_fit.get_run_dir(l, m_min, m_max)))
    else:
        print("❌ Error: Either a configuration file or --run-dir is required to monitor.")
        return False

    env = os.environ.copy()
    if cfg and getattr(cfg.environment, "brufit_dir", None) and "BRUFIT" not in env:
        env["BRUFIT"] = cfg.environment.brufit_dir

    root_sh = getattr(cfg.environment, "root_sh", None) if cfg else None
    need_source_root = bool(root_sh and os.path.isfile(root_sh) and not env.get("ROOTSYS"))

    all_ok = True
    for l_val, rdir in target_entries:
        cmd = [sys.executable, script_path, "--run-dir", rdir]
        if tbins_filter:
            cmd.extend(["--tbin", tbins_filter[0] if len(tbins_filter) == 1 else " ".join(tbins_filter)])
        else:
            cmd.extend(["--tbin", "all"])
        if l_val is not None:
            cmd.extend(["--max-l", str(l_val)])
        if cfg is not None:
            cmd.extend([
                "--m-min", str(cfg.analysis.mass_binning.min),
                "--m-max", str(cfg.analysis.mass_binning.max),
                "--n-bins", str(cfg.analysis.mass_binning.n_bins),
                "--burn-in", str(cfg.moment_fit.mcmc.burnin),
            ])
        if watch:
            cmd.append("--watch")
            cmd.extend(["--interval", str(interval)])
        if progress_only:
            cmd.append("--progress-only")
        if cancel:
            cmd.append("--cancel")
        if force:
            cmd.append("--force")
        if plot:
            cmd.append("--plot")
        if skip_trace:
            cmd.append("--skip-trace")
        if skip_posterior:
            cmd.append("--skip-posterior")
        if skip_moments:
            cmd.append("--skip-moments")
        if downsample_trace is not None:
            cmd.extend(["--downsample-trace", str(downsample_trace)])
        if replot_from_root:
            cmd.append("--replot-from-root")
        if ess_threshold is not None:
            cmd.extend(["--ess-threshold", str(ess_threshold)])
        elif cfg and hasattr(cfg.moment_fit, "diagnostics"):
            cmd.extend(["--ess-threshold", str(cfg.moment_fit.diagnostics.ess_threshold)])
        if ess_metric:
            cmd.extend(["--ess-metric", str(ess_metric)])
        elif cfg and hasattr(cfg.moment_fit, "diagnostics"):
            cmd.extend(["--ess-metric", str(cfg.moment_fit.diagnostics.ess_metric)])
        if root_filename:
            cmd.extend(["--root-filename", root_filename])

        print(f"\n[MONITOR] Launching MonitorMomentFit for run directory:\n  {rdir}")
        if dry_run:
            print(f"  [DRY-RUN] Command: {' '.join(cmd)}")
            continue

        exec_cmd = cmd
        if need_source_root:
            exec_cmd = ["bash", "-c", f"source {shlex.quote(root_sh)} && exec \"$@\"", "_", *cmd]

        try:
            ret = subprocess.run(exec_cmd, env=env)
            if ret.returncode != 0:
                all_ok = False
        except KeyboardInterrupt:
            print("\n[INFO] Monitoring stopped by user.")
            return True

    return all_ok


def run_collector(
    cfg: Optional[MomentFitConfig] = None,
    run_dir: Optional[str] = None,
    tbins_filter: Optional[List[str]] = None,
    max_l_filter: Optional[List[int]] = None,
    auto_retry: bool = False,
    no_exec: bool = False,
    non_interactive: bool = True,
    ess_threshold: Optional[float] = None,
    ess_metric: Optional[str] = None,
    min_ess_threshold: Optional[float] = None,
    max_parallel: Optional[int] = None,
    macro_path: Optional[str] = None,
    dry_run: bool = False,
) -> bool:
    """
    Launches MomentFitCollector.py across target run directories.
    """
    repo_root = Path(__file__).resolve().parent.parent.parent.parent
    configured_script = getattr(cfg.moment_fit, "collector_script", None) if cfg else None
    script_path = find_script(configured_script, "macros/MomentFitCollector.py", repo_root)
    if not script_path:
        print("❌ Error: MomentFitCollector.py script not found.")
        return False

    from .stage3_moment_fit import build_mcmc_and_diagnostics_flags
    import re

    target_entries: List[Tuple[Optional[int], str]] = []
    if run_dir:
        m_l = re.search(r"Lmax(\d+)", run_dir)
        target_entries.append((int(m_l.group(1)) if m_l else None, run_dir))
    elif cfg:
        max_ls = max_l_filter if max_l_filter is not None else cfg.analysis.max_l
        m_min = cfg.analysis.mass_binning.min
        m_max = cfg.analysis.mass_binning.max
        for l in max_ls:
            target_entries.append((l, cfg.moment_fit.get_run_dir(l, m_min, m_max)))
    else:
        print("❌ Error: Either a configuration file or --run-dir is required to collect diagnostics.")
        return False

    all_ok = True
    for dir_max_l, rdir in target_entries:
        cmd = [sys.executable, script_path, "--run-dir", rdir]
        if dir_max_l is not None:
            cmd.extend(["--max-l", str(dir_max_l)])
        if tbins_filter:
            cmd.extend(["--tbin", tbins_filter[0] if len(tbins_filter) == 1 else " ".join(tbins_filter)])
        else:
            cmd.extend(["--tbin", "all"])
        if auto_retry:
            cmd.append("--auto-retry")
        if no_exec:
            cmd.append("--no-exec")
        if non_interactive:
            cmd.append("--non-interactive")
        if ess_threshold is not None:
            cmd.extend(["--ess-threshold", str(ess_threshold)])
        elif cfg and hasattr(cfg.moment_fit, "diagnostics"):
            cmd.extend(["--ess-threshold", str(cfg.moment_fit.diagnostics.ess_threshold)])
        if ess_metric:
            cmd.extend(["--ess-metric", str(ess_metric)])
        elif cfg and hasattr(cfg.moment_fit, "diagnostics"):
            cmd.extend(["--ess-metric", str(cfg.moment_fit.diagnostics.ess_metric)])
        if min_ess_threshold is not None:
            cmd.extend(["--min-ess-threshold", str(min_ess_threshold)])
        elif cfg and hasattr(cfg.moment_fit, "diagnostics"):
            cfg_min_ess = getattr(cfg.moment_fit.diagnostics, "min_ess_threshold", None)
            if cfg_min_ess is not None:
                cmd.extend(["--min-ess-threshold", str(cfg_min_ess)])
        if max_parallel is not None:
            cmd.extend(["--max-parallel", str(max_parallel)])
        elif cfg:
            cmd.extend(["--max-parallel", str(cfg.moment_fit.max_parallel)])
        if macro_path:
            cmd.extend(["--macro-path", macro_path])
        elif cfg and hasattr(cfg.moment_fit, "macro_path"):
            cmd.extend(["--macro-path", cfg.moment_fit.macro_path])
        if cfg and hasattr(cfg, "linear_algebra") and getattr(cfg.linear_algebra, "output_base_dir", None):
            cmd.extend(["--cov-dir", cfg.linear_algebra.output_base_dir])
        if cfg and hasattr(cfg.moment_fit, "data"):
            data_cfg = cfg.moment_fit.data
            if getattr(data_cfg, "base_dir", None):
                cmd.extend(["--data-dir", data_cfg.base_dir])
            if getattr(data_cfg, "data", None):
                cmd.extend(["--data-pattern", data_cfg.data])
            if getattr(data_cfg, "mcrecon", None):
                cmd.extend(["--mcrecon-pattern", data_cfg.mcrecon])
            if getattr(data_cfg, "mcgen", None):
                cmd.extend(["--mcgen-pattern", data_cfg.mcgen])
            if getattr(data_cfg, "weight", None):
                cmd.extend(["--weight-pattern", data_cfg.weight])
        if cfg and hasattr(cfg.moment_fit, "mcmc"):
            mcmc_cfg = cfg.moment_fit.mcmc
            cmd.extend([
                "--init-scale", str(mcmc_cfg.init_scale),
                "--step", str(mcmc_cfg.step),
                "--burn-in", str(mcmc_cfg.burnin),
                "--iters", str(mcmc_cfg.iterations),
            ])
            cmd.extend(build_mcmc_and_diagnostics_flags(cfg, max_l=dir_max_l))

        print(f"\n[COLLECTOR] Launching MomentFitCollector for run directory:\n  {rdir}")
        if dry_run:
            print(f"  [DRY-RUN] Command: {' '.join(cmd)}")
            continue

        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            all_ok = False

    return all_ok


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="End-to-End MomentFit: Linear Algebra -> Export Covariance -> BruFit MCMC -> Diagnostics & Plotting",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "config",
        nargs="?",
        default=None,
        help="Path to MomentFit YAML configuration file (e.g. configs/momentfit_nominal.yaml).",
    )
    parser.add_argument(
        "--config", "-c",
        dest="config_opt",
        type=str,
        default=None,
        help="Path to YAML configuration file override.",
    )
    parser.add_argument(
        "--stage", "-s",
        nargs="+",
        choices=[
            "1", "2", "3", "4",
            "linear_algebra", "export_covariance", "moment_fit",
            "plot", "plotting",
            "monitor", "status",
            "collect", "collector", "diagnostics",
            "all",
        ],
        default=None,
        help="Specify which stage(s) to run (default: runs stages specified in YAML).",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Run Stage 4 (plot unnormalized moments vs mass).",
    )
    parser.add_argument(
        "--plot-type",
        choices=["all", "linear_algebra", "mcmc", "overlay"],
        default=None,
        help="Which unnormalized moment plots to generate: linear_algebra, mcmc, overlay, or all.",
    )
    parser.add_argument(
        "--plot-moments",
        nargs="+",
        default=None,
        help="Specific moment names to plot (e.g. 'H_0_0_0' 'H_0_2_0' or 'all').",
    )
    parser.add_argument(
        "--monitor", "--status",
        dest="monitor",
        action="store_true",
        help="Monitor BruFit MomentFit running progress, view status tables, and generate diagnostics.",
    )
    parser.add_argument(
        "--collect", "--diagnose",
        dest="collect",
        action="store_true",
        help="Diagnose BruFit MomentFit results across bins and generate targeted retries for unsuccessful bins.",
    )
    parser.add_argument(
        "--watch", "-w",
        action="store_true",
        help="Continuously monitor and refresh progress until interrupted (Ctrl+C).",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=10.0,
        help="Refresh interval in seconds when --watch is enabled (default: 10.0s).",
    )
    parser.add_argument(
        "--progress-only", "--status-only",
        dest="progress_only",
        action="store_true",
        help="Display status and progress tables only; do not generate PDF plots.",
    )
    parser.add_argument(
        "--cancel", "--stop",
        action="store_true",
        help="Cancel running MomentFit jobs for the specified t-bin(s) or run directory.",
    )
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Force terminate immediately using SIGKILL instead of graceful SIGTERM (with --cancel).",
    )
    parser.add_argument(
        "--auto-retry",
        action="store_true",
        help="Automatically execute targeted retries generated by the collector.",
    )
    parser.add_argument(
        "--no-exec",
        action="store_true",
        help="Do not execute retry script from within Python (collector script generation only).",
    )
    parser.add_argument(
        "--run-dir", "-d",
        type=str,
        default=None,
        help="Override or specify the base MomentFit run directory directly.",
    )
    parser.add_argument(
        "--ess-threshold",
        type=float,
        default=None,
        help="Effective Sample Size (ESS) threshold for fit evaluation (default: from YAML or 20.0).",
    )
    parser.add_argument(
        "--ess-metric",
        type=str,
        choices=["mean", "min"],
        default=None,
        help="ESS metric to evaluate against threshold: 'mean' or 'min' (default: from YAML or 'mean').",
    )
    parser.add_argument(
        "--min-ess-threshold",
        type=float,
        default=None,
        help="Safety floor on worst-moment ESS even when --ess-metric=mean (default: from YAML or 10.0).",
    )
    parser.add_argument(
        "--plot-traces",
        action="store_true",
        help="Generate MCMC trace plots across all moments via MonitorMomentFit.",
    )
    parser.add_argument(
        "--plot-posteriors",
        action="store_true",
        help="Generate 1D posterior distributions across all moments via MonitorMomentFit.",
    )
    parser.add_argument(
        "--downsample-trace",
        type=int,
        nargs="?",
        const=300,
        default=None,
        help="Downsample MCMC trace plots to N points per moment (default: 300).",
    )
    parser.add_argument(
        "--replot-from-root",
        action="store_true",
        help="Re-render all PDF plots directly from saved ROOT files without re-evaluating MCMC trees.",
    )
    parser.add_argument(
        "--plot-cov", "--plot-covariance",
        dest="plot_cov",
        action="store_true",
        help="Plot comparison of MCMC proposal covariance matrix before and after alignment (SPD regularization).",
    )
    parser.add_argument(
        "--plot-cov-bins",
        nargs="+",
        default=None,
        help="Filter specific mass bin labels or centers to plot covariance for (default: all).",
    )
    parser.add_argument(
        "--tbin", "-t",
        nargs="+",
        default=None,
        help="Filter specific t-bin(s) to run (e.g. 't010020' or '0').",
    )
    parser.add_argument(
        "--max-l", "-l",
        nargs="+",
        type=int,
        default=None,
        help="Filter specific maxL multipole(s) to run (e.g. 4 or 4 6 8).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview tasks and command lines without executing fits.",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate configuration file and prerequisite paths, then exit.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output directories without prompting.",
    )
    parser.add_argument(
        "--force-integral",
        action="store_true",
        help="Force recalculation of acceptance integral matrices in Stage 1.",
    )
    args = parser.parse_args(argv)

    # Infer stages from action flags if --stage was not explicitly passed
    if args.monitor or args.watch or args.cancel or args.progress_only:
        if args.stage is None:
            args.stage = ["monitor"]
        elif "monitor" not in args.stage:
            args.stage.append("monitor")

    if args.collect or args.auto_retry:
        if args.stage is None:
            args.stage = ["collect"]
        elif "collect" not in args.stage:
            args.stage.append("collect")

    if args.plot or args.plot_traces or args.plot_posteriors:
        if args.stage is None:
            args.stage = ["plot"]
        elif "plot" not in args.stage and "4" not in args.stage:
            args.stage.append("plot")

    if args.plot_cov:
        if args.stage is None:
            args.stage = ["plot_covariance"]
        elif "plot_covariance" not in args.stage:
            args.stage.append("plot_covariance")

    cfg_file = args.config_opt or args.config
    if not cfg_file:
        # Allow running monitoring, collector, or covariance plotting directly if --run-dir is provided
        if args.run_dir and (args.monitor or args.watch or args.collect or args.cancel or args.plot_cov or (args.stage and any(s in ("monitor", "collect", "status", "diagnostics", "plot_covariance", "cov") for s in args.stage))):
            args.config = None
            return args

        parser.error(
            "A configuration file is required.\n"
            "Usage:\n"
            "  pixi run momentfit <config.yaml> [options]\n"
            "  python run_momentfit.py <config.yaml> [options]\n"
            "Example:\n"
            "  pixi run momentfit configs/momentfit_nominal.yaml"
        )
    args.config = cfg_file
    return args


def normalize_stage_names(stages: List[str]) -> List[str]:
    """Map numeric and alias stage names to canonical names."""
    mapping = {
        "1": "linear_algebra",
        "linear_algebra": "linear_algebra",
        "2": "export_covariance",
        "export_covariance": "export_covariance",
        "3": "moment_fit",
        "moment_fit": "moment_fit",
        "4": "plot",
        "plot": "plot",
        "plotting": "plot",
        "monitor": "monitor",
        "status": "monitor",
        "collect": "collect",
        "collector": "collect",
        "diagnostics": "collect",
        "plot_covariance": "plot_covariance",
        "plot_cov": "plot_covariance",
        "cov": "plot_covariance",
        "covariance": "plot_covariance",
    }
    canonical = []
    for s in stages:
        s_lower = str(s).lower()
        if s_lower == "all":
            return ["linear_algebra", "export_covariance", "moment_fit", "plot"]
        if s_lower in mapping:
            canonical.append(mapping[s_lower])
        else:
            raise ValueError(f"Unknown stage name: {s}")
    deduped = list(dict.fromkeys(canonical))  # preserve order without duplicates
    # "collect" and "monitor" operate on Stage 3 (moment_fit) outputs; when passed alongside
    # Stage 3 / moment_fit, skip running the fitting process and execute monitoring/collection directly.
    if "moment_fit" in deduped and any(s in ("collect", "monitor") for s in deduped):
        deduped = [s for s in deduped if s != "moment_fit"]
    return deduped


def normalize_tbin_filter(tbins_in: Optional[List[str]]) -> Optional[List[str]]:
    """Normalize t-bin input (supports indices 0..4 or names t010020..)."""
    if tbins_in is None:
        return None
    index_map = {
        "0": "t010020",
        "1": "t020032",
        "2": "t032050",
        "3": "t050075",
        "4": "t075100",
    }
    normalized = []
    for item in tbins_in:
        for sub in item.replace(",", " ").split():
            if sub.lower() == "all":
                return None
            if sub in index_map:
                normalized.append(index_map[sub])
            else:
                normalized.append(sub)
    return normalized


def run_momentfit(
    config_path: Optional[str] = None,
    stages: Optional[List[str]] = None,
    tbins_filter: Optional[List[str]] = None,
    max_l_filter: Optional[List[int]] = None,
    dry_run: bool = False,
    validate_only: bool = False,
    overwrite: bool = False,
    force_integral: bool = False,
    plot_type: Optional[str] = None,
    moments_filter: Optional[List[str]] = None,
    run_dir: Optional[str] = None,
    watch: bool = False,
    interval: float = 10.0,
    progress_only: bool = False,
    cancel: bool = False,
    force: bool = False,
    auto_retry: bool = False,
    no_exec: bool = False,
    ess_threshold: Optional[float] = None,
    ess_metric: Optional[str] = None,
    min_ess_threshold: Optional[float] = None,
    plot_traces: bool = False,
    plot_posteriors: bool = False,
    downsample_trace: Optional[int] = None,
    replot_from_root: bool = False,
    plot_cov: bool = False,
    plot_cov_bins: Optional[List[str]] = None,
) -> int:
    """
    Main MomentFit orchestration routine.
    """
    total_start_time = time.time()
    print("=" * 80)
    print(" MOMENTFIT ANALYSIS WORKFLOW")
    if config_path:
        print(f" Config file: {config_path}")
    if run_dir:
        print(f" Run directory: {run_dir}")
    print("=" * 80)

    cfg: Optional[MomentFitConfig] = None
    if config_path:
        try:
            cfg = load_momentfit_config(config_path)
        except Exception as e:
            print(f"❌ Failed to load configuration: {e}")
            return 1

        # Apply CLI overrides
        if overwrite:
            cfg.moment_fit.overwrite = True
        if force_integral:
            cfg.linear_algebra.force_integral_matrix = True
    elif not run_dir:
        print("❌ Error: Either a configuration file or --run-dir is required.")
        return 1

    tbins_norm = normalize_tbin_filter(tbins_filter)

    # Determine which stages to run
    if stages:
        active_stages = normalize_stage_names(stages)
    elif cfg:
        active_stages = normalize_stage_names(cfg.stages)
    else:
        active_stages = ["monitor"]

    # Validate configuration if running core pipeline stages (1, 2, 3)
    if cfg and not set(active_stages).issubset({"monitor", "collect", "plot", "plot_covariance"}):
        errors = validate_momentfit_config(cfg)
        if errors:
            print("\n[VALIDATION FAILED] Please address the following issues:")
            for err in errors:
                print(f"  ❌ {err}")
            return 1
        else:
            print("✅ Configuration is valid!")
    elif cfg and validate_only:
        errors = validate_momentfit_config(cfg)
        if errors:
            print("\n[VALIDATION FAILED] Please address the following issues:")
            for err in errors:
                print(f"  ❌ {err}")
            return 1
        print("✅ Configuration is valid!")
        return 0

    if validate_only:
        print("\nValidation mode requested: exiting with success.")
        return 0

    target_tbins_str = tbins_norm if tbins_norm else (cfg.analysis.tbins if cfg else "all")
    target_max_l_str = max_l_filter if max_l_filter else (cfg.analysis.max_l if cfg else "auto")
    print(f"\nPipeline Execution Sequence: {' -> '.join(active_stages)}")
    print(f"Target t-bins:             {target_tbins_str}")
    print(f"Multipoles Max L:          {target_max_l_str}")
    print(f"Dry-run mode:              {dry_run}\n")

    stage_timings = {}
    stage_statuses = {}

    from .stage1_linear_algebra import run_stage1
    from .stage2_export_covariance import run_stage2
    from .stage3_moment_fit import run_stage3
    from .stage4_plot_moments import run_stage4

    # --------------------------------------------------------------------------
    # STAGE 1: Linear Algebra Moments Calculation
    # --------------------------------------------------------------------------
    if "linear_algebra" in active_stages:
        if cfg is None:
            print("❌ Error: Stage 1 requires a configuration YAML file.")
            return 1
        t0 = time.time()
        success = run_stage1(
            cfg=cfg,
            tbins_filter=tbins_norm,
            max_l_filter=max_l_filter,
            dry_run=dry_run,
        )
        stage_timings["linear_algebra"] = time.time() - t0
        stage_statuses["linear_algebra"] = "SUCCESS" if success else "FAILED"
        if not success and not dry_run:
            print("❌ Pipeline halted due to error in Stage 1.")
            return 1

    # --------------------------------------------------------------------------
    # STAGE 2: Export Covariance & Format for BruFit
    # --------------------------------------------------------------------------
    if "export_covariance" in active_stages:
        if cfg is None:
            print("❌ Error: Stage 2 requires a configuration YAML file.")
            return 1
        t0 = time.time()
        success = run_stage2(
            cfg=cfg,
            tbins_filter=tbins_norm,
            max_l_filter=max_l_filter,
            dry_run=dry_run,
            plot_cov=plot_cov,
            plot_cov_bins=plot_cov_bins,
        )
        stage_timings["export_covariance"] = time.time() - t0
        stage_statuses["export_covariance"] = "SUCCESS" if success else "FAILED"
        if not success and not dry_run:
            print("❌ Pipeline halted due to error in Stage 2.")
            return 1

    # --------------------------------------------------------------------------
    # STAGE 3: MomentFit with BruFit MCMC
    # --------------------------------------------------------------------------
    if "moment_fit" in active_stages:
        if cfg is None:
            print("❌ Error: Stage 3 requires a configuration YAML file.")
            return 1
        t0 = time.time()
        success = run_stage3(
            cfg=cfg,
            tbins_filter=tbins_norm,
            max_l_filter=max_l_filter,
            dry_run=dry_run,
        )
        stage_timings["moment_fit"] = time.time() - t0
        stage_statuses["moment_fit"] = "SUCCESS" if success else "FAILED"
        if not success and not dry_run:
            print("❌ Workflow halted due to error in Stage 3.")
            return 1

    # --------------------------------------------------------------------------
    # STAGE: Monitoring & Progress Inspection (MonitorMomentFit.py)
    # --------------------------------------------------------------------------
    if "monitor" in active_stages:
        t0 = time.time()
        success = run_monitor(
            cfg=cfg,
            run_dir=run_dir,
            tbins_filter=tbins_norm,
            max_l_filter=max_l_filter,
            watch=watch,
            interval=interval,
            progress_only=progress_only,
            cancel=cancel,
            force=force,
            plot=plot_traces or plot_posteriors,
            skip_trace=not plot_traces if (plot_traces or plot_posteriors) else False,
            skip_posterior=not plot_posteriors if (plot_traces or plot_posteriors) else False,
            downsample_trace=downsample_trace,
            replot_from_root=replot_from_root,
            ess_threshold=ess_threshold,
            ess_metric=ess_metric,
            min_ess_threshold=min_ess_threshold,
            dry_run=dry_run,
        )
        stage_timings["monitor"] = time.time() - t0
        stage_statuses["monitor"] = "SUCCESS" if success else "FAILED"
        if not success and not dry_run:
            print("❌ Workflow halted due to error in Monitor.")
            return 1

    # --------------------------------------------------------------------------
    # STAGE: Collector & Diagnostics / Targeted Retries (MomentFitCollector.py)
    # --------------------------------------------------------------------------
    if "collect" in active_stages:
        t0 = time.time()
        success = run_collector(
            cfg=cfg,
            run_dir=run_dir,
            tbins_filter=tbins_norm,
            max_l_filter=max_l_filter,
            auto_retry=auto_retry,
            no_exec=no_exec,
            ess_threshold=ess_threshold,
            ess_metric=ess_metric,
            min_ess_threshold=min_ess_threshold,
            dry_run=dry_run,
        )
        stage_timings["collect"] = time.time() - t0
        stage_statuses["collect"] = "SUCCESS" if success else "FAILED"
        if not success and not dry_run:
            print("❌ Workflow halted due to error in Collector.")
            return 1

    # --------------------------------------------------------------------------
    # STAGE 4: Plot Unnormalized Moments vs Mass
    # --------------------------------------------------------------------------
    if "plot" in active_stages:
        t0 = time.time()
        success = True
        if cfg is not None:
            success = run_stage4(
                cfg=cfg,
                tbins_filter=tbins_norm,
                max_l_filter=max_l_filter,
                plot_type=plot_type,
                moments_filter=moments_filter,
                dry_run=dry_run,
            )
        if plot_traces or plot_posteriors or replot_from_root:
            m_ok = run_monitor(
                cfg=cfg,
                run_dir=run_dir,
                tbins_filter=tbins_norm,
                max_l_filter=max_l_filter,
                plot=True,
                skip_moments=True,
                skip_trace=not plot_traces if (plot_traces or plot_posteriors) else False,
                skip_posterior=not plot_posteriors if (plot_traces or plot_posteriors) else False,
                downsample_trace=downsample_trace,
                replot_from_root=replot_from_root,
                ess_threshold=ess_threshold,
                ess_metric=ess_metric,
                dry_run=dry_run,
            )
            success = success and m_ok

        stage_timings["plot"] = time.time() - t0
        stage_statuses["plot"] = "SUCCESS" if success else "FAILED"
        if not success and not dry_run:
            print("❌ Workflow halted due to error in Stage 4.")
            return 1

    # --------------------------------------------------------------------------
    # STAGE: Plot Covariance Matrices (Before and After Alignment)
    # --------------------------------------------------------------------------
    if "plot_covariance" in active_stages:
        t0 = time.time()
        from .plot_covariance import plot_covariance_from_config, plot_covariance_for_directory
        if run_dir:
            plots = plot_covariance_for_directory(
                cov_dir=run_dir,
                mass_bins_filter=plot_cov_bins,
                tbin=tbins_norm[0] if tbins_norm else "",
                max_l=max_l_filter[0] if max_l_filter else "",
            )
            success = len(plots) > 0
        elif cfg is not None:
            res = plot_covariance_from_config(
                cfg=cfg,
                tbins_filter=tbins_norm,
                max_l_filter=max_l_filter,
                mass_bins_filter=plot_cov_bins,
            )
            success = any(len(p) > 0 for p in res.values())
        else:
            print("❌ Error: Either a configuration file or --run-dir is required to plot covariance matrices.")
            success = False

        stage_timings["plot_covariance"] = time.time() - t0
        stage_statuses["plot_covariance"] = "SUCCESS" if success else "FAILED"
        if not success and not dry_run:
            print("❌ Workflow halted due to error in Covariance Plotting.")
            return 1

    # --------------------------------------------------------------------------
    # Workflow Summary
    # --------------------------------------------------------------------------
    total_time = time.time() - total_start_time
    print("\n" + "=" * 80)
    print(" MOMENTFIT WORKFLOW EXECUTION SUMMARY")
    print("=" * 80)
    print(f" {'Stage':<25} | {'Status':<10} | {'Duration':<12}")
    print("-" * 55)
    for stage in active_stages:
        st = stage_statuses.get(stage, "SKIPPED")
        dur = f"{stage_timings.get(stage, 0.0):.1f} s"
        icon = "✅" if st == "SUCCESS" else ("❌" if st == "FAILED" else "⏭️")
        print(f" {stage:<25} | {icon} {st:<7} | {dur:<12}")
    print("-" * 55)
    print(f" Total Wall Time: {total_time:.1f} s")
    print("=" * 80 + "\n")

    return 0


# Backward-compatible function alias
run_pipeline = run_momentfit


def main():
    args = parse_args()
    ret = run_momentfit(
        config_path=args.config,
        stages=args.stage,
        tbins_filter=args.tbin,
        max_l_filter=args.max_l,
        dry_run=args.dry_run,
        validate_only=args.validate,
        overwrite=args.overwrite,
        force_integral=args.force_integral,
        plot_type=args.plot_type,
        moments_filter=args.plot_moments,
        run_dir=args.run_dir,
        watch=args.watch,
        interval=args.interval,
        progress_only=args.progress_only,
        cancel=args.cancel,
        force=args.force,
        auto_retry=args.auto_retry,
        no_exec=args.no_exec,
        ess_threshold=args.ess_threshold,
        ess_metric=args.ess_metric,
        min_ess_threshold=args.min_ess_threshold,
        plot_traces=args.plot_traces,
        plot_posteriors=args.plot_posteriors,
        downsample_trace=args.downsample_trace,
        replot_from_root=args.replot_from_root,
        plot_cov=args.plot_cov,
        plot_cov_bins=args.plot_cov_bins,
    )
    sys.exit(ret)


if __name__ == "__main__":
    main()

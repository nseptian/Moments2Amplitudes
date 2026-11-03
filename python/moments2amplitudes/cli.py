"""
Command-line interface (CLI) for Moments2Amplitudes (mom2amp).
"""

from __future__ import annotations
import sys
import os
import argparse
from pathlib import Path

from .config import load_config, validate_config
from .driver import run_multi_start, run_mcmc_propagation
from .rundir import get_run_directory
from .summary import FitSummary
from .plotter import run_all_diagnostics_and_plots


def cmd_validate(args: argparse.Namespace) -> int:
    """Validate a YAML configuration file."""
    config_path = Path(args.config)
    print(f"\nValidating configuration: {config_path}")
    try:
        cfg = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1

    errors = validate_config(cfg)
    if errors:
        print("\n[VALIDATION FAILED] Found issues:")
        for err in errors:
            print(f"  - {err}")
        return 1

    print("\n[OK] Configuration is valid!")
    print(f"  Project:        {cfg.project.name}")
    print(f"  Moments Dir:    {cfg.dataset.moments_dir}")
    print(f"  Mass Bins:      {len(cfg.get_mass_bins())} bins ({cfg.get_mass_bins()[0]:.3f} to {cfg.get_mass_bins()[-1]:.3f} GeV)")
    print(f"  Use Covariance: {cfg.dataset.use_covariance}")
    print(f"  Chi2 Included L:{cfg.chi2.included_l}")
    print(f"  Waves Defined:  {len(cfg.waves)}")
    for w in cfg.waves:
        res_names = [r.name for r in w.resonances] if w.resonances else "none"
        print(f"    * L={w.l} ({w.name}): type={w.type.value}, model={w.model.value}, resonances={res_names}")
    print(f"  Multi-start:    {cfg.multi_start.n_starts} trials on {cfg.multi_start.n_threads} threads")
    print(f"  MCMC Prop:      enabled={cfg.mcmc_propagation.enabled} ({cfg.mcmc_propagation.max_entries} entries on {cfg.mcmc_propagation.n_workers} workers)")
    print(f"  BruFit Path:    {cfg.environment.brufit_dir}")
    return 0


def cmd_fit(args: argparse.Namespace) -> int:
    """Run parallel multi-start randomized fit to determine global minimum."""
    config_path = Path(args.config)
    try:
        cfg = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1

    errors = validate_config(cfg)
    if errors:
        print(f"[ERROR] Invalid configuration: {errors}")
        return 1

    res = run_multi_start(
        cfg,
        n_starts=args.n_starts,
        n_threads=args.n_threads,
        output_dir=args.output_dir,
    )

    return 0 if res.get("converged", False) else 1


def cmd_propagate(args: argparse.Namespace) -> int:
    """Run parallel MCMC posterior uncertainty propagation."""
    config_path = Path(args.config)
    try:
        cfg = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1

    cfg.mcmc_propagation.enabled = True
    if args.mcmc_tree:
        cfg.mcmc_propagation.mcmc_tree_path = args.mcmc_tree

    errors = validate_config(cfg)
    if errors:
        print(f"[ERROR] Invalid configuration: {errors}")
        return 1

    res = run_mcmc_propagation(
        cfg,
        best_params_file=args.best_params,
        n_workers=args.n_workers,
        max_entries=args.max_entries,
        output_file=args.output_file,
    )

    return 0 if "output_file" in res else 1


def cmd_run(args: argparse.Namespace) -> int:
    """Run full pipeline: multi-start central fit -> MCMC uncertainty propagation."""
    config_path = Path(args.config)
    try:
        cfg = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1

    errors = validate_config(cfg)
    if errors:
        print(f"[ERROR] Invalid configuration: {errors}")
        return 1

    # Stage 1: Multi-start fit
    print("\n>>> Stage 1: Running Multi-Start Global Minimum Fit...")
    res_fit = run_multi_start(
        cfg,
        n_starts=args.n_starts,
        n_threads=args.n_threads,
        output_dir=args.output_dir,
    )

    if not res_fit.get("converged", False):
        print("[ERROR] Stage 1 failed to converge.")
        return 1

    # Stage 2: MCMC propagation (if enabled or requested)
    if cfg.mcmc_propagation.enabled or args.propagate:
        print("\n>>> Stage 2: Running MCMC Posterior Uncertainty Propagation...")
        best_file = res_fit.get("best_file")
        res_mcmc = run_mcmc_propagation(
            cfg,
            best_params_file=best_file,
            n_workers=args.n_workers,
            max_entries=args.max_entries,
            output_file=args.output_file,
        )
        if "output_file" not in res_mcmc:
            return 1

    print("\n[SUCCESS] Pipeline completed successfully!")
    return 0


def cmd_plot(args: argparse.Namespace) -> int:
    """Generate publication amplitude and phase plots with MCMC error bands."""
    config_path = Path(args.config)
    try:
        cfg = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1

    run_dir = getattr(args, "run_dir", None) or getattr(args, "output_dir", None)
    rundir = get_run_directory(cfg, run_dir)
    print(f"\nGenerating plots from run directory: {rundir.root}")
    
    if not rundir.nominal_best_root.is_file():
        legacy_best = rundir.root / "central_best_params.root"
        if legacy_best.is_file():
            rundir.save_nominal_best(legacy_best)
        else:
            print(f"[ERROR] No nominal best fit found at {rundir.nominal_best_root}")
            return 1

    formats = tuple(args.formats.split(",")) if hasattr(args, "formats") and args.formats else ("png", "pdf")
    res = run_all_diagnostics_and_plots(cfg, rundir, formats=formats)
    print(f"[SUCCESS] Generated {len(res['plots'])} plot files in: {rundir.plots_dir}")
    return 0


def cmd_summary(args: argparse.Namespace) -> int:
    """Display fit summary table and save reports."""
    config_path = Path(args.config)
    try:
        cfg = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1

    run_dir = getattr(args, "run_dir", None) or getattr(args, "output_dir", None)
    rundir = get_run_directory(cfg, run_dir)

    target_root = rundir.nominal_best_root
    if not target_root.is_file():
        legacy_best = rundir.root / "central_best_params.root"
        if legacy_best.is_file():
            target_root = legacy_best
        else:
            print(f"[ERROR] No fit results found in {rundir.root}")
            return 1

    mcmc_target = rundir.mcmc_propagated_root if rundir.mcmc_propagated_root.is_file() else None

    summary = FitSummary(cfg, rundir)
    summary.load_from_root(target_root, mcmc_target)
    if not summary.parameters:
        logs = sorted(rundir.fit_starts_dir.glob("fit_seed_*.log"))
        if logs:
            summary.load_from_log(logs[0])
    txt_path, md_path = summary.save_reports()
    print(f"\n{summary.generate_ascii_table()}\n")
    print(f"[OK] Saved reports:")
    print(f"  - Text:     {txt_path}")
    print(f"  - Markdown: {md_path}")
    return 0


def main():
    parser = argparse.ArgumentParser(
        prog="mom2amp",
        description="Moments2Amplitudes - Mass-dependent partial wave analysis framework",
    )
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Validate
    p_val = subparsers.add_parser("validate", help="Validate a YAML fit configuration")
    p_val.add_argument("config", help="Path to config.yaml")

    # Fit
    p_fit = subparsers.add_parser("fit", help="Run multi-start global minimum fit")
    p_fit.add_argument("config", help="Path to config.yaml")
    p_fit.add_argument("--n-starts", type=int, default=None, help="Number of randomized starts")
    p_fit.add_argument("--n-threads", type=int, default=None, help="Parallel worker threads")
    p_fit.add_argument("--output-dir", "--run-dir", dest="output_dir", type=str, default=None, help="Run / output directory")

    # Propagate
    p_prop = subparsers.add_parser("propagate", help="Run MCMC posterior uncertainty propagation")
    p_prop.add_argument("config", help="Path to config.yaml")
    p_prop.add_argument("--best-params", type=str, default=None, help="Path to central best parameters ROOT file")
    p_prop.add_argument("--mcmc-tree", type=str, default=None, help="Path to MCMC tree file")
    p_prop.add_argument("--n-workers", type=int, default=None, help="Number of parallel workers")
    p_prop.add_argument("--max-entries", type=int, default=None, help="Max posterior entries to propagate")
    p_prop.add_argument("--output-file", type=str, default=None, help="Output merged ROOT file")
    p_prop.add_argument("--run-dir", type=str, default=None, help="Run directory")

    # Plot
    p_plot = subparsers.add_parser("plot", help="Generate publication amplitude plots and error bands")
    p_plot.add_argument("config", help="Path to config.yaml")
    p_plot.add_argument("--run-dir", "--output-dir", dest="run_dir", type=str, default=None, help="Run directory")
    p_plot.add_argument("--formats", type=str, default="png,pdf", help="Image output formats (comma-separated, default: png,pdf)")

    # Summary
    p_sum = subparsers.add_parser("summary", help="Display fit parameter summary and save report tables")
    p_sum.add_argument("config", help="Path to config.yaml")
    p_sum.add_argument("--run-dir", "--output-dir", dest="run_dir", type=str, default=None, help="Run directory")

    # Run (All in one)
    p_run = subparsers.add_parser("run", help="Run full pipeline (multi-start fit + MCMC propagation)")
    p_run.add_argument("config", help="Path to config.yaml")
    p_run.add_argument("--n-starts", type=int, default=None, help="Number of randomized starts")
    p_run.add_argument("--n-threads", type=int, default=None, help="Parallel worker threads")
    p_run.add_argument("--output-dir", "--run-dir", dest="output_dir", type=str, default=None, help="Run directory")
    p_run.add_argument("--propagate", action="store_true", help="Force MCMC propagation even if disabled in config")
    p_run.add_argument("--n-workers", type=int, default=None, help="Number of MCMC workers")
    p_run.add_argument("--max-entries", type=int, default=None, help="Max posterior entries to propagate")
    p_run.add_argument("--output-file", type=str, default=None, help="Output merged ROOT file")

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(0)

    if args.command == "validate":
        sys.exit(cmd_validate(args))
    elif args.command == "fit":
        sys.exit(cmd_fit(args))
    elif args.command == "propagate":
        sys.exit(cmd_propagate(args))
    elif args.command == "plot":
        sys.exit(cmd_plot(args))
    elif args.command == "summary":
        sys.exit(cmd_summary(args))
    elif args.command == "run":
        sys.exit(cmd_run(args))
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()


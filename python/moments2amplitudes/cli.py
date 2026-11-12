"""
Command-line interface (CLI) for Moments2Amplitudes (mom2amp).
"""

from __future__ import annotations
import sys
import os
import copy
import argparse
from pathlib import Path

from .config import load_config, validate_config, parse_t_bins_input
from .driver import run_multi_start, run_multi_start_tbins, run_mcmc_propagation, run_mcmc_propagation_tbins
from .rundir import get_run_directory, RunDirectory
from .summary import FitSummary
from .plotter import run_all_diagnostics_and_plots, plot_tbins_parallel


def apply_cli_overrides(cfg: FitConfig, args: argparse.Namespace) -> None:
    """Apply common CLI overrides such as t-bins to configuration."""
    tbins_arg = getattr(args, "tbins", None)
    tbin_arg = getattr(args, "t_bin", None)

    if tbins_arg is not None:
        parsed = parse_t_bins_input(tbins_arg)
        if parsed:
            cfg.dataset.t_bins = parsed
            cfg.dataset.t_bin = parsed[0]
    elif tbin_arg is not None:
        parsed = parse_t_bins_input(tbin_arg)
        if parsed:
            cfg.dataset.t_bin = parsed[0]
            cfg.dataset.t_bins = parsed


def cmd_validate(args: argparse.Namespace) -> int:
    """Validate a YAML configuration file."""
    config_path = Path(args.config)
    print(f"\nValidating configuration: {config_path}")
    try:
        cfg = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1

    apply_cli_overrides(cfg, args)

    errors = validate_config(cfg)
    if errors:
        print("\n[VALIDATION FAILED] Found issues:")
        for err in errors:
            print(f"  - {err}")
        return 1

    print("\n[OK] Configuration is valid!")
    print(f"  Project:        {cfg.project.name}")
    print(f"  Moments Dir:    {cfg.dataset.moments_dir}")
    tbins = cfg.dataset.get_t_bins()
    if len(tbins) > 1:
        print(f"  t-Bins ({len(tbins)}):     {', '.join(tbins)}")
        for t in tbins:
            print(f"    * {t}: {cfg.dataset.get_moments_dir(t)}")
    elif len(tbins) == 1:
        print(f"  t-bin:          {tbins[0]} (resolved dir: {cfg.dataset.get_moments_dir()})")
    print(f"  Mass Bins:      {len(cfg.get_mass_bins())} bins ({cfg.get_mass_bins()[0]:.3f} to {cfg.get_mass_bins()[-1]:.3f} GeV)")
    print(f"  Use Covariance: {cfg.dataset.use_covariance}")
    print(f"  Chi2 Included L:{cfg.chi2.included_l}")
    print(f"  Waves Defined:  {len(cfg.waves)}")
    for w in cfg.waves:
        res_names = [r.name for r in w.resonances] if w.resonances else "none"
        print(f"    * L={w.l} ({w.name}): type={w.type.value}, model={w.model.value}, resonances={res_names}")
    print(f"  Multi-start:    {cfg.multi_start.n_starts} trials on {cfg.multi_start.n_threads} threads")
    print(f"  MCMC Prop:      enabled={cfg.mcmc_propagation.enabled} ({cfg.mcmc_propagation.max_entries} entries on {cfg.mcmc_propagation.n_workers} workers)")
    brufit_disp = (
        (cfg.environment.brufit_dir if cfg.environment and cfg.environment.brufit_dir else None)
        or os.environ.get("BRUFIT", "not set (using system default)")
    )
    print(f"  BruFit Path:    {brufit_disp}")
    return 0


def cmd_fit(args: argparse.Namespace) -> int:
    """Run parallel multi-start randomized fit to determine global minimum."""
    config_path = Path(args.config)
    try:
        cfg = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1

    apply_cli_overrides(cfg, args)

    errors = validate_config(cfg)
    if errors:
        print(f"[ERROR] Invalid configuration: {errors}")
        return 1

    tbins = cfg.dataset.get_t_bins()
    n_tbin_workers = getattr(args, "tbin_workers", None)

    res = run_multi_start_tbins(
        cfg,
        tbins=tbins,
        n_starts=args.n_starts,
        n_threads=args.n_threads,
        n_tbin_workers=n_tbin_workers,
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

    apply_cli_overrides(cfg, args)

    cfg.mcmc_propagation.enabled = True
    if args.mcmc_tree:
        cfg.mcmc_propagation.mcmc_tree_path = args.mcmc_tree

    errors = validate_config(cfg)
    if errors:
        print(f"[ERROR] Invalid configuration: {errors}")
        return 1

    tbins = cfg.dataset.get_t_bins()
    n_tbin_workers = getattr(args, "tbin_workers", None)
    res = run_mcmc_propagation_tbins(
        cfg,
        tbins=tbins,
        best_params_file=args.best_params,
        n_workers=args.n_workers,
        n_tbin_workers=n_tbin_workers,
        max_entries=args.max_entries,
        output_file=args.output_file,
        run_dir=getattr(args, "run_dir", None),
    )

    success = res.get("success", False) or "output_file" in res
    if success and getattr(args, "draw_plots", False):
        print("\n>>> Automatically drawing propagated amplitude plots and distributions...")
        args.only_propagated = True
        cmd_plot(args)

    return 0 if success else 1


def cmd_run(args: argparse.Namespace) -> int:
    """Run full pipeline: multi-start central fit -> MCMC uncertainty propagation."""
    config_path = Path(args.config)
    try:
        cfg = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1

    apply_cli_overrides(cfg, args)

    errors = validate_config(cfg)
    if errors:
        print(f"[ERROR] Invalid configuration: {errors}")
        return 1

    # Stage 1: Multi-start fit
    print("\n>>> Stage 1: Running Multi-Start Global Minimum Fit...")
    tbins = cfg.dataset.get_t_bins()
    n_tbin_workers = getattr(args, "tbin_workers", None)
    res_fit = run_multi_start_tbins(
        cfg,
        tbins=tbins,
        n_starts=args.n_starts,
        n_threads=args.n_threads,
        n_tbin_workers=n_tbin_workers,
        output_dir=args.output_dir,
    )

    if not res_fit.get("converged", False):
        print("[ERROR] Stage 1 failed to converge.")
        return 1

    # Stage 2: MCMC propagation (if enabled or requested)
    if cfg.mcmc_propagation.enabled or args.propagate:
        print("\n>>> Stage 2: Running MCMC Posterior Uncertainty Propagation...")
        res_mcmc = run_mcmc_propagation_tbins(
            cfg,
            tbins=tbins,
            best_params_file=args.best_params,
            n_workers=args.n_workers,
            n_tbin_workers=n_tbin_workers,
            max_entries=args.max_entries,
            output_file=args.output_file,
            run_dir=getattr(args, "output_dir", None) or getattr(args, "run_dir", None),
        )
        if not (res_mcmc.get("success", False) or "output_file" in res_mcmc):
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

    apply_cli_overrides(cfg, args)

    if hasattr(args, "n_lowest_chi2") and args.n_lowest_chi2 is not None:
        cfg.plotting.n_lowest_chi2 = args.n_lowest_chi2
    if hasattr(args, "unnormalized_moments_file") and args.unnormalized_moments_file:
        cfg.plotting.unnormalized_moments_file = args.unnormalized_moments_file
    if hasattr(args, "moments") and args.moments:
        cfg.plotting.moments = [m.strip() for m in args.moments.split(",")]

    formats = tuple(args.formats.split(",")) if hasattr(args, "formats") and args.formats else ("png", "pdf")
    tbins = cfg.dataset.get_t_bins()

    only_prop = getattr(args, "only_propagated", False)
    only_dist = getattr(args, "only_distributions", False)
    only_amp = getattr(args, "only_amplitudes", False)

    if len(tbins) > 1 and not getattr(args, "t_bin", None):
        max_workers = getattr(args, "n_workers", None) or getattr(cfg.multi_start, "n_tbin_workers", None)
        res = plot_tbins_parallel(
            cfg,
            tbins=tbins,
            run_dir=getattr(args, "run_dir", None) or getattr(args, "output_dir", None),
            formats=formats,
            max_workers=max_workers,
            only_propagated=only_prop,
            only_distributions=only_dist,
            only_amplitudes=only_amp,
        )
        return 0 if (res.get("success", False) or res.get("total_plots", 0) > 0) else 1

    run_dir = getattr(args, "run_dir", None) or getattr(args, "output_dir", None)
    rundir = get_run_directory(cfg, run_dir)
    print(f"\nGenerating plots from run directory: {rundir.root}")
    
    if not rundir.nominal_best_root.is_file():
        legacy_best = rundir.root / "central_best_params.root"
        if legacy_best.is_file():
            rundir.save_nominal_best(legacy_best)
        elif rundir.mcmc_propagated_root.is_file():
            rundir.save_nominal_best(rundir.mcmc_propagated_root)
        else:
            print(f"[ERROR] No nominal best fit found at {rundir.nominal_best_root}")
            return 1

    res = run_all_diagnostics_and_plots(
        cfg,
        rundir,
        formats=formats,
        only_propagated=only_prop,
        only_distributions=only_dist,
        only_amplitudes=only_amp,
    )
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

    apply_cli_overrides(cfg, args)

    tbins = cfg.dataset.get_t_bins()
    if len(tbins) > 1 and not getattr(args, "t_bin", None):
        print(f"\nGenerating fit summaries across {len(tbins)} t-bins...")
        for t_name in tbins:
            cfg_t = copy.deepcopy(cfg)
            cfg_t.dataset.t_bin = t_name
            rundir = get_run_directory(cfg_t, getattr(args, "run_dir", None))
            print(f"\n=======================================================")
            print(f" Fit Summary: t-bin [{t_name}] ({rundir.root})")
            print(f"=======================================================")
            target_root = rundir.nominal_best_root
            if not target_root.is_file():
                legacy_best = rundir.root / "central_best_params.root"
                if legacy_best.is_file():
                    target_root = legacy_best
                else:
                    print(f"[WARNING] No fit results found in {rundir.root}. Skipping.")
                    continue
            mcmc_target = rundir.mcmc_propagated_root if rundir.mcmc_propagated_root.is_file() else None
            summary = FitSummary(cfg_t, rundir)
            summary.load_from_root(target_root, mcmc_target)
            if not summary.parameters:
                logs = sorted(rundir.fit_starts_dir.glob("fit_seed_*.log"))
                if logs:
                    summary.load_from_log(logs[0])
            txt_path, md_path = summary.save_reports()
            print(f"\n{summary.generate_ascii_table()}\n")
        return 0

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


def cmd_cross_section(args: argparse.Namespace) -> int:
    """Calculate reflectivity-separated cross section for a fitted resonance."""
    config_path = Path(args.config)
    try:
        cfg = load_config(config_path)
    except Exception as e:
        print(f"[ERROR] Failed to load configuration: {e}")
        return 1

    apply_cli_overrides(cfg, args)

    # Overrides
    if hasattr(args, "resonance") and args.resonance:
        cfg.cross_section.resonance = args.resonance
    if hasattr(args, "luminosity") and args.luminosity is not None:
        cfg.cross_section.luminosity = args.luminosity
    if hasattr(args, "branching_ratio") and args.branching_ratio is not None:
        cfg.cross_section.branching_ratio = args.branching_ratio
    if hasattr(args, "t_min") and args.t_min is not None:
        cfg.cross_section.t_min = args.t_min
    if hasattr(args, "t_max") and args.t_max is not None:
        cfg.cross_section.t_max = args.t_max
    if hasattr(args, "delta_t") and args.delta_t is not None:
        cfg.cross_section.delta_t = args.delta_t

    from .cross_section import (
        calculate_resonance_cross_section,
        format_cross_section_table,
        save_cross_section_reports,
        plot_cross_section_vs_t,
        CrossSectionResult,
    )

    mcmc_file = getattr(args, "mcmc_file", None)
    propagate_mcmc = not getattr(args, "no_mcmc", False)

    results: List[CrossSectionResult] = []

    target_runs = args.runs if hasattr(args, "runs") and args.runs else None
    if not target_runs and len(cfg.dataset.get_t_bins()) > 1:
        target_runs = []
        for t_name in cfg.dataset.get_t_bins():
            cfg_t = copy.deepcopy(cfg)
            cfg_t.dataset.t_bin = t_name
            rdir = get_run_directory(cfg_t, getattr(args, "run_dir", None))
            target_runs.append(str(rdir.root))

    if target_runs:
        print(f"\nCalculating cross-section across {len(target_runs)} run directories...")
        for r_dir in target_runs:
            rundir = RunDirectory(r_dir)
            try:
                res = calculate_resonance_cross_section(
                    cfg,
                    rundir,
                    resonance_name=cfg.cross_section.resonance,
                    luminosity=cfg.cross_section.luminosity,
                    branching_ratio=cfg.cross_section.branching_ratio,
                    mcmc_root_file=mcmc_file,
                    propagate_mcmc=propagate_mcmc,
                )
                results.append(res)
                err_str = f"± {res.cs_total_mcmc_err:.4f} (MCMC)" if res.has_mcmc else f"± {res.cs_total_err:.4f} (Minuit)"
                print(f"[{res.t_min:.3f} < -t < {res.t_max:.3f}] dsigma/dt: Total={res.cs_total:.4f} {err_str}, Pos={res.cs_pos:.4f}, Neg={res.cs_neg:.4f} mub/GeV^2")
            except Exception as e:
                print(f"[WARNING] Could not process run dir {r_dir}: {e}")

        if not results:
            print("[ERROR] No cross-section results evaluated.")
            return 1

        primary_rundir = RunDirectory(target_runs[0])
        plot_cross_section_vs_t(
            results,
            primary_rundir.plots_dir,
            phase1_comparison=not getattr(args, "no_phase1", False),
            theory_band_file=getattr(args, "theory_band", None),
        )
        return 0

    run_dir = getattr(args, "run_dir", None) or getattr(args, "output_dir", None)
    rundir = get_run_directory(cfg, run_dir)
    print(f"\nCalculating cross section from run directory: {rundir.root}")

    try:
        res = calculate_resonance_cross_section(
            cfg,
            rundir,
            resonance_name=cfg.cross_section.resonance,
            t_min=cfg.cross_section.t_min,
            t_max=cfg.cross_section.t_max,
            luminosity=cfg.cross_section.luminosity,
            branching_ratio=cfg.cross_section.branching_ratio,
            mcmc_root_file=mcmc_file,
            propagate_mcmc=propagate_mcmc,
        )
    except Exception as e:
        print(f"[ERROR] Failed to calculate cross section: {e}")
        return 1

    print(f"\n{format_cross_section_table(res)}\n")
    txt_p, md_p, json_p = save_cross_section_reports(res, rundir.reports_dir)
    print(f"[OK] Saved reports:")
    print(f"  - Text:     {txt_p}")
    print(f"  - Markdown: {md_p}")
    print(f"  - JSON:     {json_p}")

    plot_cross_section_vs_t(
        [res],
        rundir.plots_dir,
        phase1_comparison=not getattr(args, "no_phase1", False),
        theory_band_file=getattr(args, "theory_band", None),
    )
    return 0


def cmd_compile(args: argparse.Namespace) -> int:
    """Check, clean, and recompile ACLiC C++ shared libraries."""
    from .compiler import is_aclic_updated, clean_aclic_artifacts, compile_aclic_libraries, ensure_aclic_updated

    if getattr(args, "clean", False):
        removed = clean_aclic_artifacts()
        print(f"[ACLiC] Cleaned {len(removed)} build artifacts.")

    if getattr(args, "check", False):
        up_to_date, reason = is_aclic_updated()
        status = "UP-TO-DATE" if up_to_date else "OUTDATED"
        print(f"[ACLiC] Status: {status} ({reason})")
        return 0 if up_to_date else 1

    force = getattr(args, "force", False) or getattr(args, "clean", False)
    try:
        ensure_aclic_updated(force=force)
        return 0
    except Exception as e:
        print(f"[ERROR] Compilation failed: {e}")
        return 1


def main():
    parser = argparse.ArgumentParser(
        prog="mom2amp",
        description="Moments2Amplitudes - Mass-dependent partial wave analysis framework",
    )
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Validate
    p_val = subparsers.add_parser("validate", help="Validate a YAML fit configuration")
    p_val.add_argument("config", help="Path to config.yaml")
    p_val.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_val.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to evaluate (e.g. 'all', 't010020 t020032', '0 1 2')")

    # Fit
    p_fit = subparsers.add_parser("fit", help="Run multi-start global minimum fit")
    p_fit.add_argument("config", help="Path to config.yaml")
    p_fit.add_argument("--n-starts", type=int, default=None, help="Number of randomized starts")
    p_fit.add_argument("--n-threads", type=int, default=None, help="Parallel worker threads per fit")
    p_fit.add_argument("--tbin-workers", type=int, default=None, help="Max concurrent t-bins running in parallel")
    p_fit.add_argument("--output-dir", "--run-dir", dest="output_dir", type=str, default=None, help="Run / output directory (supports {t_bin})")
    p_fit.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_fit.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to run in parallel (e.g. 'all', 't010020 t020032', '0 1 2')")

    # Propagate
    p_prop = subparsers.add_parser("propagate", help="Run MCMC posterior uncertainty propagation")
    p_prop.add_argument("config", help="Path to config.yaml")
    p_prop.add_argument("--best-params", type=str, default=None, help="Path to central best parameters ROOT file")
    p_prop.add_argument("--mcmc-tree", type=str, default=None, help="Path to MCMC tree file")
    p_prop.add_argument("--n-workers", type=int, default=None, help="Number of parallel workers per t-bin")
    p_prop.add_argument("--tbin-workers", type=int, default=None, help="Max concurrent t-bins running in parallel")
    p_prop.add_argument("--max-entries", type=int, default=None, help="Max posterior entries to propagate")
    p_prop.add_argument("--output-file", type=str, default=None, help="Output merged ROOT file")
    p_prop.add_argument("--run-dir", type=str, default=None, help="Run directory")
    p_prop.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_prop.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to evaluate (e.g. 'all', 't010020 t020032', '0 1 2')")
    p_prop.add_argument("--plot", "--draw", dest="draw_plots", action="store_true", help="Automatically draw propagated amplitude distributions and credible bands after propagation completes")

    # Plot
    p_plot = subparsers.add_parser("plot", help="Generate publication amplitude plots and error bands")
    p_plot.add_argument("config", help="Path to config.yaml")
    p_plot.add_argument("--run-dir", "--output-dir", dest="run_dir", type=str, default=None, help="Run directory")
    p_plot.add_argument("--formats", type=str, default="png,pdf", help="Image output formats (comma-separated, default: png,pdf)")
    p_plot.add_argument("--n-lowest-chi2", type=int, default=None, help="Number of lowest chi2 fit curves to draw (default: from config or 10)")
    p_plot.add_argument("--unnormalized-moments-file", type=str, default=None, help="Path to UnnormalizedMoments_vs_Mass.root")
    p_plot.add_argument("--moments", type=str, default=None, help="Moments to plot (e.g. 'all', 'L=2', 'H_0_0_0,H_0_2_0')")
    p_plot.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_plot.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to evaluate (e.g. 'all', 't010020 t020032', '0 1 2')")
    p_plot.add_argument("--n-workers", "--tbin-workers", dest="n_workers", type=int, default=None, help="Number of parallel worker processes across t-bins")
    p_plot.add_argument("--only-propagated", "--propagated-only", dest="only_propagated", action="store_true", help="Only draw plots derived from propagated MCMC results (credible bands, individual amplitudes, and posterior distributions)")
    p_plot.add_argument("--only-distributions", "--distributions-only", "--only-distribution", "--distribution-only", dest="only_distributions", action="store_true", help="Only draw 1D posterior distributions of propagated MCMC parameters and yields")
    p_plot.add_argument("--only-amplitudes", "--amplitudes-only", dest="only_amplitudes", action="store_true", help="Only draw amplitude curves and credible bands (individual amplitudes and wave overviews)")

    # Summary
    p_sum = subparsers.add_parser("summary", help="Display fit parameter summary and save report tables")
    p_sum.add_argument("config", help="Path to config.yaml")
    p_sum.add_argument("--run-dir", "--output-dir", dest="run_dir", type=str, default=None, help="Run directory")
    p_sum.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_sum.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to evaluate (e.g. 'all', 't010020 t020032', '0 1 2')")

    # Cross-Section
    p_cs = subparsers.add_parser("cross-section", aliases=["cs"], help="Calculate reflectivity-separated cross-section for a fitted resonance")
    p_cs.add_argument("config", help="Path to config.yaml")
    p_cs.add_argument("--resonance", type=str, default=None, help="Target resonance name (default: from config or a2_1320)")
    p_cs.add_argument("--luminosity", type=float, default=None, help="Integrated luminosity in pb^-1")
    p_cs.add_argument("--br", "--branching-ratio", dest="branching_ratio", type=float, default=None, help="Total branching ratio")
    p_cs.add_argument("--t-min", type=float, default=None, help="Minimum -t (GeV^2)")
    p_cs.add_argument("--t-max", type=float, default=None, help="Maximum -t (GeV^2)")
    p_cs.add_argument("--delta-t", type=float, default=None, help="Bin width delta t (GeV^2)")
    p_cs.add_argument("--run-dir", "--output-dir", dest="run_dir", type=str, default=None, help="Run directory")
    p_cs.add_argument("--runs", nargs="*", default=None, help="Multiple run directories to calculate cross section across t-bins")
    p_cs.add_argument("--theory-band", type=str, default=None, help="Path to theory prediction band file")
    p_cs.add_argument("--no-phase1", action="store_true", help="Disable overlay of published Phase-1 values")
    p_cs.add_argument("--mcmc-file", type=str, default=None, help="Path to MCMC propagated ROOT file for posterior statistical uncertainty")
    p_cs.add_argument("--no-mcmc", action="store_true", help="Disable MCMC uncertainty propagation even if MCMC file is found")
    p_cs.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_cs.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to evaluate (e.g. 'all', 't010020 t020032', '0 1 2')")

    # Run (All in one)
    p_run = subparsers.add_parser("run", help="Run full pipeline (multi-start fit + MCMC propagation)")
    p_run.add_argument("config", help="Path to config.yaml")
    p_run.add_argument("--n-starts", type=int, default=None, help="Number of randomized starts")
    p_run.add_argument("--n-threads", type=int, default=None, help="Parallel worker threads")
    p_run.add_argument("--tbin-workers", type=int, default=None, help="Max concurrent t-bins running in parallel")
    p_run.add_argument("--output-dir", "--run-dir", dest="output_dir", type=str, default=None, help="Run directory")
    p_run.add_argument("--propagate", action="store_true", help="Force MCMC propagation even if disabled in config")
    p_run.add_argument("--n-workers", type=int, default=None, help="Number of MCMC workers")
    p_run.add_argument("--max-entries", type=int, default=None, help="Max posterior entries to propagate")
    p_run.add_argument("--output-file", type=str, default=None, help="Output merged ROOT file")
    p_run.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_run.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to run in parallel (e.g. 'all', 't010020 t020032', '0 1 2')")

    # Compile
    p_comp = subparsers.add_parser("compile", help="Check, clean, and recompile ACLiC C++ shared libraries")
    p_comp.add_argument("--clean", action="store_true", help="Clean stale ACLiC build artifacts before compiling")
    p_comp.add_argument("--force", "-f", action="store_true", help="Force unconditional recompilation")
    p_comp.add_argument("--check", action="store_true", help="Only check if libraries are up-to-date without compiling")

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(0)

    if args.command == "compile":
        sys.exit(cmd_compile(args))
    elif args.command == "validate":
        sys.exit(cmd_validate(args))
    elif args.command == "fit":
        sys.exit(cmd_fit(args))
    elif args.command == "propagate":
        sys.exit(cmd_propagate(args))
    elif args.command == "plot":
        sys.exit(cmd_plot(args))
    elif args.command == "summary":
        sys.exit(cmd_summary(args))
    elif args.command in ("cross-section", "cs"):
        sys.exit(cmd_cross_section(args))
    elif args.command == "run":
        sys.exit(cmd_run(args))
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()


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
from .plotter import run_all_diagnostics_and_plots, plot_tbins_parallel, plot_moments_across_tbins, plot_amplitudes_across_tbins, plot_multistart_landscape_all_tbins
from .environment import get_default_max_workers
from .schema import FitConfig, parse_seed_range


def apply_cli_overrides(cfg: FitConfig, args: argparse.Namespace) -> None:
    """Apply common CLI overrides such as t-bins and multi-start seed range to configuration."""
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

    # Multi-start CLI overrides
    if getattr(cfg, "multi_start", None) is not None:
        n_starts_arg = getattr(args, "n_starts", None) or getattr(args, "trials", None)
        if n_starts_arg is not None:
            cfg.multi_start.n_starts = int(n_starts_arg)

        seed_range_arg = getattr(args, "seed_range", None)
        if seed_range_arg is not None:
            cfg.multi_start.seed_range = seed_range_arg

        seed_start_arg = getattr(args, "seed_start", None)
        if seed_start_arg is not None:
            cfg.multi_start.seed_start = int(seed_start_arg)

        seed_end_arg = getattr(args, "seed_end", None)
        if seed_end_arg is not None:
            cfg.multi_start.seed_end = int(seed_end_arg)

        skip_existing_arg = getattr(args, "skip_existing", None)
        if skip_existing_arg is not None:
            cfg.multi_start.skip_existing = bool(skip_existing_arg)

    style_arg = getattr(args, "style", None)
    if style_arg is not None and getattr(cfg, "plotting", None) is not None:
        cfg.plotting.style = style_arg

    dist_bins_arg = getattr(args, "distribution_bins", None)
    if dist_bins_arg is not None and getattr(cfg, "plotting", None) is not None:
        cfg.plotting.distribution_bins = int(dist_bins_arg)

    chi2_bins_arg = getattr(args, "chi2_bins", None)
    if chi2_bins_arg is not None and getattr(cfg, "plotting", None) is not None:
        cfg.plotting.chi2_landscape_bins = int(chi2_bins_arg)

    ref_amp_arg = getattr(args, "reference_amplitude", None) or getattr(args, "anchor_amplitude", None)
    if ref_amp_arg is not None and getattr(cfg, "plotting", None) is not None:
        cfg.plotting.reference_amplitude = str(ref_amp_arg).strip()
        cfg.plotting.anchor_amplitude = str(ref_amp_arg).strip()

    stride_arg = getattr(args, "stride", None)
    if stride_arg is not None and getattr(cfg, "mcmc_propagation", None) is not None:
        cfg.mcmc_propagation.stride = stride_arg


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
    print(f"  Chi2 Included L:{cfg.chi2.included_L}")
    print(f"  Waves Defined:  {len(cfg.waves)}")
    for w in cfg.waves:
        res_names = [r.name for r in w.resonances] if w.resonances else "none"
        print(f"    * L={w.l} ({w.name}): type={w.type.value}, model={w.model.value}, resonances={res_names}")
    threads_display = (
        cfg.multi_start.n_threads
        if cfg.multi_start.n_threads is not None
        else get_default_max_workers()
    )
    seeds = cfg.multi_start.get_seeds()
    s_start, s_end = cfg.multi_start.get_seed_range()
    skip_str = ", skip_existing=True" if cfg.multi_start.skip_existing else ""
    print(f"  Multi-start:    {len(seeds)} trials (seeds {s_start}..{s_end}{skip_str}) on {threads_display} threads")
    print(f"  MCMC Prop:      enabled={cfg.mcmc_propagation.enabled} ({cfg.mcmc_propagation.target_entries} entries on {cfg.mcmc_propagation.n_workers} workers)")
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
        n_starts=getattr(args, "n_starts", None) or getattr(args, "trials", None),
        seed_start=getattr(args, "seed_start", None),
        seed_end=getattr(args, "seed_end", None),
        seed_range=getattr(args, "seed_range", None),
        skip_existing=getattr(args, "skip_existing", None),
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
    req_entries = getattr(args, "target_entries", None) if getattr(args, "target_entries", None) is not None else getattr(args, "max_entries", None)
    res = run_mcmc_propagation_tbins(
        cfg,
        tbins=tbins,
        best_params_file=args.best_params,
        n_workers=args.n_workers,
        n_tbin_workers=n_tbin_workers,
        target_entries=req_entries,
        stride=getattr(args, "stride", None),
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
        n_starts=getattr(args, "n_starts", None) or getattr(args, "trials", None),
        seed_start=getattr(args, "seed_start", None),
        seed_end=getattr(args, "seed_end", None),
        seed_range=getattr(args, "seed_range", None),
        skip_existing=getattr(args, "skip_existing", None),
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
            target_entries=getattr(args, "target_entries", None) if getattr(args, "target_entries", None) is not None else getattr(args, "max_entries", None),
            stride=getattr(args, "stride", None),
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

    if getattr(args, "moments_vs_t", False):
        run_dir = getattr(args, "run_dir", None) or getattr(args, "output_dir", None)
        res_plots = plot_moments_across_tbins(
            cfg,
            tbins=tbins,
            run_dir=run_dir,
            formats=formats,
            style=getattr(args, "style", None),
            sharey=getattr(args, "moments_vs_t_sharey", False),
        )
        print(f"[SUCCESS] Generated {len(res_plots)} multi-t-bin moment comparison plot files.")
        return 0 if res_plots else 1

    if getattr(args, "amplitudes_vs_t", False):
        run_dir = getattr(args, "run_dir", None) or getattr(args, "output_dir", None)
        res_plots = plot_amplitudes_across_tbins(
            cfg,
            tbins=tbins,
            run_dir=run_dir,
            formats=formats,
            style=getattr(args, "style", None),
            sharey=getattr(args, "amplitudes_vs_t_sharey", False),
        )
        print(f"[SUCCESS] Generated {len(res_plots)} multi-t-bin amplitude comparison plot files.")
        return 0 if res_plots else 1

    if getattr(args, "chi2_landscape", False) or getattr(args, "multistart_landscape", False):
        run_dir = getattr(args, "run_dir", None) or getattr(args, "output_dir", None)
        parent_out = Path(run_dir) if run_dir else Path(cfg.multi_start.output_dir if cfg.multi_start.output_dir else "runs")
        plots_out = parent_out / "plots" if (parent_out / "plots").is_dir() else parent_out
        results_by_tbin = {}
        for t_name in tbins:
            cfg_t = copy.deepcopy(cfg)
            cfg_t.dataset.t_bin = t_name
            rdir = get_run_directory(cfg_t, run_dir)
            results_by_tbin[t_name] = rdir
        show_grid = getattr(cfg.plotting, "show_grid", False) if getattr(cfg, "plotting", None) else False
        chi2_bins = getattr(args, "chi2_bins", None) or (getattr(cfg.plotting, "chi2_landscape_bins", None) if getattr(cfg, "plotting", None) else None)
        res_plots = plot_multistart_landscape_all_tbins(
            results_by_tbin,
            output_dir=plots_out,
            project_name=cfg.project.name,
            plateau_delta=cfg.multi_start.plateau_delta_chi2,
            formats=formats,
            style=getattr(args, "style", None),
            show_grid=show_grid,
            bins=chi2_bins,
        )
        if res_plots:
            print(f"[SUCCESS] Generated multi-start Chi2 landscape plot across {len(tbins)} t-bins: {res_plots[0]}")
            return 0
        else:
            print(f"[WARNING] No multi-start fit starts found to plot in {parent_out}")
            return 1

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

        # Automatically generate multi-start Chi2 landscape across all t-bins if starts exist and not filtered
        if not (only_prop or only_dist or only_amp):
            run_dir = getattr(args, "run_dir", None) or getattr(args, "output_dir", None)
            parent_out = Path(run_dir) if run_dir else Path(cfg.multi_start.output_dir if cfg.multi_start.output_dir else "runs")
            plots_out = parent_out / "plots" if (parent_out / "plots").is_dir() else parent_out
            results_by_tbin = {}
            for t_name in tbins:
                cfg_t = copy.deepcopy(cfg)
                cfg_t.dataset.t_bin = t_name
                rdir = get_run_directory(cfg_t, run_dir)
                results_by_tbin[t_name] = rdir
            show_grid = getattr(cfg.plotting, "show_grid", False) if getattr(cfg, "plotting", None) else False
            chi2_bins = getattr(args, "chi2_bins", None) or (getattr(cfg.plotting, "chi2_landscape_bins", None) if getattr(cfg, "plotting", None) else None)
            res_plots = plot_multistart_landscape_all_tbins(
                results_by_tbin,
                output_dir=plots_out,
                project_name=cfg.project.name,
                plateau_delta=cfg.multi_start.plateau_delta_chi2,
                formats=formats,
                style=getattr(args, "style", None),
                show_grid=show_grid,
                bins=chi2_bins,
            )
            if res_plots:
                print(f"[OK] Saved multi-start Chi2 landscape across all {len(tbins)} t-bins to: {plots_out / 'multistart_chi2_landscape_all_tbins.png'}")

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
                    summary.load_from_log(logs[0], mcmc_file=mcmc_target)
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
            summary.load_from_log(logs[0], mcmc_file=mcmc_target)
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
        parse_t_bin_from_string,
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
            cfg_t = copy.deepcopy(cfg)
            parsed_t = parse_t_bin_from_string(rundir.root.name)
            if parsed_t:
                cfg_t.dataset.t_bin = rundir.root.name
                t_lo_val, t_hi_val = parsed_t
            else:
                t_lo_val, t_hi_val = None, None
            try:
                res = calculate_resonance_cross_section(
                    cfg_t,
                    rundir,
                    resonance_name=cfg.cross_section.resonance,
                    t_min=t_lo_val,
                    t_max=t_hi_val,
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

        first_p = Path(target_runs[0])
        project_plots_dir = first_p.parent / "plots" if (first_p.parent.is_dir() and first_p.name.startswith("t")) else RunDirectory(target_runs[0]).plots_dir
        project_plots_dir.mkdir(parents=True, exist_ok=True)
        plot_cross_section_vs_t(
            results,
            project_plots_dir,
            phase1_comparison=not getattr(args, "no_phase1", False),
            theory_band_file=getattr(args, "theory_band", None),
            style=getattr(args, "style", None) or getattr(cfg.plotting, "style", "default"),
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
        style=getattr(args, "style", None) or getattr(cfg.plotting, "style", "default"),
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
    up_to_date, reason = is_aclic_updated()
    if up_to_date and not force:
        print(f"[ACLiC] Status: UP-TO-DATE ({reason})")
        print("[ACLiC] All shared libraries in src/ are already built. Use 'pixi run compile --force' (or -f) to recompile.")
        return 0

    try:
        ensure_aclic_updated(force=force)
        return 0
    except Exception as e:
        print(f"[ERROR] Compilation failed: {e}")
        return 1


def cmd_momentfit(args: argparse.Namespace) -> int:
    """Run the MomentFit workflow (linear algebra -> export covariance -> BruFit MCMC -> diagnostics & plotting)."""
    from .momentfit.runner import run_momentfit
    config_path = getattr(args, "config_opt", None) or getattr(args, "config", None)
    run_dir = getattr(args, "run_dir", None)
    stages = args.stage

    if not config_path and not run_dir:
        print("[ERROR] A configuration file or --run-dir is required.\nUsage: mom2amp momentfit <config.yaml> [options]")
        return 1

    if getattr(args, "monitor", False) or getattr(args, "watch", False) or getattr(args, "cancel", False) or getattr(args, "progress_only", False):
        if stages is None:
            stages = ["monitor"]
        elif "monitor" not in stages:
            stages.append("monitor")

    if getattr(args, "calculate_acceptance", None) is not None:
        if stages is None:
            stages = ["calculate_acceptance"]
        elif "calculate_acceptance" not in stages:
            stages.append("calculate_acceptance")

    if getattr(args, "retry", False):
        if stages is None:
            stages = ["retry"]
        elif "retry" not in stages:
            stages.append("retry")
    elif getattr(args, "collect", False) or getattr(args, "auto_retry", False):
        if stages is None:
            stages = ["collect"]
        elif "collect" not in stages:
            stages.append("collect")

    if getattr(args, "plot", False) or getattr(args, "plot_traces", False) or getattr(args, "plot_posteriors", False):
        if stages is None:
            stages = ["plot"]
        elif "plot" not in stages and "4" not in stages:
            stages.append("plot")

    if getattr(args, "plot_cov", False) or getattr(args, "plot_covariance", False):
        if stages is None:
            stages = ["plot_covariance"]
        elif "plot_covariance" not in stages:
            stages.append("plot_covariance")

    return run_momentfit(
        config_path=config_path,
        stages=stages,
        tbins_filter=args.tbin,
        max_l_filter=args.max_l,
        dry_run=args.dry_run,
        validate_only=args.validate,
        overwrite=args.overwrite,
        force_integral=args.force_integral,
        plot_type=getattr(args, "plot_type", None),
        moments_filter=getattr(args, "plot_moments", None),
        run_dir=run_dir,
        watch=getattr(args, "watch", False),
        interval=getattr(args, "interval", 10.0),
        progress_only=getattr(args, "progress_only", False),
        cancel=getattr(args, "cancel", False),
        force=getattr(args, "force", False),
        auto_retry=getattr(args, "auto_retry", False),
        no_exec=getattr(args, "no_exec", False),
        max_retries=getattr(args, "max_retries", None),
        max_parallel=getattr(args, "max_parallel", None),
        ess_threshold=getattr(args, "ess_threshold", None),
        ess_metric=getattr(args, "ess_metric", None),
        min_ess_threshold=getattr(args, "min_ess_threshold", None),
        plot_traces=getattr(args, "plot_traces", False),
        plot_posteriors=getattr(args, "plot_posteriors", False),
        downsample_trace=getattr(args, "downsample_trace", None),
        replot_from_root=getattr(args, "replot_from_root", False),
        plot_cov=getattr(args, "plot_cov", False) or getattr(args, "plot_covariance", False),
        plot_cov_bins=getattr(args, "plot_cov_bins", None),
        calculate_acceptance=getattr(args, "calculate_acceptance", None),
        skip_acceptance=getattr(args, "skip_acceptance", False),
        force_acceptance=getattr(args, "force_acceptance", False),
    )


def main():
    parser = argparse.ArgumentParser(
        prog="mom2amp",
        description="Moments2Amplitudes - Mass-dependent partial wave analysis framework",
    )
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Validate
    p_val = subparsers.add_parser("validate", help="Validate a YAML fit configuration")
    p_val.add_argument("config", help="Path to config.yaml")
    p_val.add_argument("--trials", "--n-starts", dest="n_starts", type=int, default=None, help="Number of randomized starts (trials)")
    p_val.add_argument("--seed-range", dest="seed_range", type=str, default=None, help="Seed range start:end (e.g. 5001:10000 or 5001-10000)")
    p_val.add_argument("--seed-start", dest="seed_start", type=int, default=None, help="Starting seed (1-based, inclusive)")
    p_val.add_argument("--seed-end", dest="seed_end", type=int, default=None, help="Ending seed (inclusive)")
    p_val.add_argument("--skip-existing", dest="skip_existing", action="store_true", default=None, help="Skip seeds that have already completed in the output directory")
    p_val.add_argument("--no-skip-existing", dest="skip_existing", action="store_false", help="Force re-running existing seeds")
    p_val.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_val.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to evaluate (e.g. 'all', 't010020 t020032', '0 1 2')")

    # Fit
    p_fit = subparsers.add_parser("fit", help="Run multi-start global minimum fit")
    p_fit.add_argument("config", help="Path to config.yaml")
    p_fit.add_argument("--n-starts", "--trials", dest="n_starts", type=int, default=None, help="Number of randomized starts (or trials)")
    p_fit.add_argument("--seed-range", dest="seed_range", type=str, default=None, help="Seed range start:end (e.g. 5001:10000 or 5001-10000)")
    p_fit.add_argument("--seed-start", dest="seed_start", type=int, default=None, help="Starting seed (1-based, inclusive)")
    p_fit.add_argument("--seed-end", dest="seed_end", type=int, default=None, help="Ending seed (inclusive)")
    p_fit.add_argument("--skip-existing", dest="skip_existing", action="store_true", default=None, help="Skip seeds whose fit results already exist")
    p_fit.add_argument("--no-skip-existing", dest="skip_existing", action="store_false", help="Force re-running existing seeds")
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
    p_prop.add_argument("--target-entries", "--max-entries", dest="target_entries", type=int, default=None, help="Target number of posterior entries to propagate")
    p_prop.add_argument("--stride", type=int, default=None, help="MCMC thinning stride override (default: automatically calculated to span chain for target_entries)")
    p_prop.add_argument("--output-file", type=str, default=None, help="Output merged ROOT file")
    p_prop.add_argument("--run-dir", type=str, default=None, help="Run directory")
    p_prop.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_prop.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to evaluate (e.g. 'all', 't010020 t020032', '0 1 2')")
    p_prop.add_argument("--plot", "--draw", dest="draw_plots", action="store_true", help="Automatically draw propagated amplitude distributions and credible bands after propagation completes")
    p_prop.add_argument("--style", "--plot-style", dest="style", type=str, default=None, help="Plot style for --plot (default: from config or 'default'; options: default, LHCb, ROOT, ATLAS, CMS, ALICE)")
    p_prop.add_argument("--reference-amplitude", "--anchor-amplitude", "--anchor", "--ref-amp", dest="reference_amplitude", type=str, default=None, help="Reference / anchor amplitude name for relative phase plotting (e.g. 'a_2_1', 'a_2_2', 'a_0_0')")

    # Plot
    p_plot = subparsers.add_parser("plot", help="Generate publication amplitude plots and error bands")
    p_plot.add_argument("config", help="Path to config.yaml")
    p_plot.add_argument("--run-dir", "--output-dir", dest="run_dir", type=str, default=None, help="Run directory")
    p_plot.add_argument("--formats", type=str, default="png,pdf", help="Image output formats (comma-separated, default: png,pdf)")
    p_plot.add_argument("--style", "--plot-style", dest="style", type=str, default=None, help="Plot style (default: from config or 'default'; options: default, LHCb, ROOT, ATLAS, CMS, ALICE)")
    p_plot.add_argument("--n-lowest-chi2", type=int, default=None, help="Number of lowest chi2 fit curves to draw (default: from config or 10)")
    p_plot.add_argument("--unnormalized-moments-file", type=str, default=None, help="Path to UnnormalizedMoments_vs_Mass.root")
    p_plot.add_argument("--moments", type=str, default=None, help="Moments to plot (e.g. 'all', 'L=2', 'H_0_0_0,H_0_2_0')")
    p_plot.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_plot.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to evaluate (e.g. 'all', 't010020 t020032', '0 1 2')")
    p_plot.add_argument("--n-workers", "--tbin-workers", dest="n_workers", type=int, default=None, help="Number of parallel worker processes across t-bins")
    p_plot.add_argument("--only-propagated", "--propagated-only", dest="only_propagated", action="store_true", help="Only draw plots derived from propagated MCMC results (credible bands, individual amplitudes, and posterior distributions)")
    p_plot.add_argument("--only-distributions", "--distributions-only", "--only-distribution", "--distribution-only", dest="only_distributions", action="store_true", help="Only draw 1D posterior distributions of propagated MCMC parameters and yields")
    p_plot.add_argument("--only-amplitudes", "--amplitudes-only", dest="only_amplitudes", action="store_true", help="Only draw amplitude curves and credible bands (individual amplitudes and wave overviews)")
    p_plot.add_argument("--chi2-landscape", "--multistart-landscape", dest="chi2_landscape", action="store_true", help="Plot combined multi-start Chi2 distribution landscapes across all t-bins side-by-side in one row")
    p_plot.add_argument("--moments-vs-t", dest="moments_vs_t", action="store_true", help="Plot each moment vs mass across all t-bins as a side-by-side plot (N_tbins x 1) to inspect evolution with -t")
    p_plot.add_argument("--moments-vs-t-sharey", dest="moments_vs_t_sharey", action="store_true", help="Share y-axis scale across all t-bins in moments-vs-t plots")
    p_plot.add_argument("--amplitudes-vs-t", dest="amplitudes_vs_t", action="store_true", help="Plot each partial wave amplitude and total intensity vs mass across all t-bins as a side-by-side plot (1 x N_tbins) to inspect evolution with -t")
    p_plot.add_argument("--amplitudes-vs-t-sharey", dest="amplitudes_vs_t_sharey", action="store_true", help="Share y-axis scale across all t-bins in amplitudes-vs-t plots")
    p_plot.add_argument("--reference-amplitude", "--anchor-amplitude", "--anchor", "--ref-amp", dest="reference_amplitude", type=str, default=None, help="Reference / anchor amplitude name for relative phase plotting (e.g. 'a_2_1', 'a_2_2', 'a_0_0')")
    p_plot.add_argument("--distribution-bins", "--dist-bins", dest="distribution_bins", type=int, default=None, help="Number of histogram bins for 1D MCMC parameter and yield posterior distribution plots (default: 50)")
    p_plot.add_argument("--chi2-bins", dest="chi2_bins", type=int, default=None, help="Number of histogram bins for multi-start Chi2 landscape plots (default: 60)")

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
    p_cs.add_argument("--style", "--plot-style", dest="style", type=str, default=None, help="Plot style for cross section (default: from config or 'default'; options: default, LHCb, ROOT, ATLAS, CMS, ALICE)")

    # Run (All in one)
    p_run = subparsers.add_parser("run", help="Run full pipeline (multi-start fit + MCMC propagation)")
    p_run.add_argument("config", help="Path to config.yaml")
    p_run.add_argument("--n-starts", "--trials", dest="n_starts", type=int, default=None, help="Number of randomized starts (or trials)")
    p_run.add_argument("--seed-range", dest="seed_range", type=str, default=None, help="Seed range start:end (e.g. 5001:10000 or 5001-10000)")
    p_run.add_argument("--seed-start", dest="seed_start", type=int, default=None, help="Starting seed (1-based, inclusive)")
    p_run.add_argument("--seed-end", dest="seed_end", type=int, default=None, help="Ending seed (inclusive)")
    p_run.add_argument("--skip-existing", dest="skip_existing", action="store_true", default=None, help="Skip seeds whose fit results already exist")
    p_run.add_argument("--no-skip-existing", dest="skip_existing", action="store_false", help="Force re-running existing seeds")
    p_run.add_argument("--n-threads", type=int, default=None, help="Parallel worker threads")
    p_run.add_argument("--tbin-workers", type=int, default=None, help="Max concurrent t-bins running in parallel")
    p_run.add_argument("--output-dir", "--run-dir", dest="output_dir", type=str, default=None, help="Run directory")
    p_run.add_argument("--propagate", action="store_true", help="Force MCMC propagation even if disabled in config")
    p_run.add_argument("--n-workers", type=int, default=None, help="Number of MCMC workers")
    p_run.add_argument("--target-entries", "--max-entries", dest="target_entries", type=int, default=None, help="Target number of posterior entries to propagate")
    p_run.add_argument("--stride", type=int, default=None, help="MCMC thinning stride override (default: automatically calculated to span chain for target_entries)")
    p_run.add_argument("--output-file", type=str, default=None, help="Output merged ROOT file")
    p_run.add_argument("--t-bin", "--tbin", dest="t_bin", type=str, default=None, help="Single t-bin directory override (e.g. t010020, t020032)")
    p_run.add_argument("--tbins", "--t-bins", dest="tbins", nargs="*", default=None, help="t-bins to run in parallel (e.g. 'all', 't010020 t020032', '0 1 2')")

    # Compile
    p_comp = subparsers.add_parser("compile", help="Check, clean, and recompile ACLiC C++ shared libraries")
    p_comp.add_argument("--clean", action="store_true", help="Clean stale ACLiC build artifacts before compiling")
    p_comp.add_argument("--force", "-f", action="store_true", help="Force unconditional recompilation")
    p_comp.add_argument("--check", action="store_true", help="Only check if libraries are up-to-date without compiling")

    # MomentFit Pipeline
    p_mfit = subparsers.add_parser(
        "momentfit",
        aliases=["pipeline"],
        help="Run end-to-end moments pipeline (Linear Algebra -> Export Covariance -> MomentFit MCMC -> Diagnostics & Plotting)",
    )
    p_mfit.add_argument("config", nargs="?", default=None, help="Path to MomentFit YAML config (e.g. configs/momentfit_nominal.yaml)")
    p_mfit.add_argument("--config", "-c", dest="config_opt", default=None, help="Path to MomentFit YAML configuration file override")
    p_mfit.add_argument(
        "--stage", "-s",
        nargs="+",
        choices=[
            "1", "2", "3", "4",
            "linear_algebra", "export_covariance", "moment_fit",
            "plot", "plotting",
            "monitor", "status",
            "calculate_acceptance", "calc_acceptance", "acceptance",
            "collect", "collector", "diagnostics",
            "retry", "refit", "retries",
            "plot_covariance", "plot_cov", "cov",
            "all",
        ],
        default=None,
        help="Specify which stage(s) to run (default: from YAML)",
    )
    p_mfit.add_argument("--monitor", "--status", dest="monitor", action="store_true", help="Monitor BruFit MomentFit running progress, view status tables, and generate diagnostics")
    p_mfit.add_argument("--collect", "--diagnose", dest="collect", action="store_true", help="Diagnose BruFit MomentFit results across bins and generate targeted retries for unsuccessful bins")
    p_mfit.add_argument("--retry", dest="retry", action="store_true", help="Collect diagnostics and execute targeted single-bin retry rounds (RunMomentFit_retry.sh)")
    p_mfit.add_argument("--watch", "-w", action="store_true", help="Continuously monitor and refresh progress until interrupted (Ctrl+C)")
    p_mfit.add_argument("--interval", type=float, default=10.0, help="Refresh interval in seconds when --watch is enabled (default: 10.0s)")
    p_mfit.add_argument("--progress-only", "--status-only", dest="progress_only", action="store_true", help="Display status and progress tables only; do not generate PDF plots")
    p_mfit.add_argument("--cancel", "--stop", action="store_true", help="Cancel running MomentFit jobs for the specified t-bin(s) or run directory")
    p_mfit.add_argument("--force", "-f", action="store_true", help="Force terminate immediately using SIGKILL instead of graceful SIGTERM (with --cancel)")
    p_mfit.add_argument("--auto-retry", action="store_true", help="Automatically execute targeted retries generated by the collector")
    p_mfit.add_argument("--no-exec", action="store_true", help="Do not execute retry script from within Python (collector script generation only)")
    p_mfit.add_argument("--max-retries", type=int, default=None, help="Maximum retry rounds to execute when --retry / --stage retry is used")
    p_mfit.add_argument("--max-parallel", "-j", type=int, default=None, help="Override maximum concurrent single-bin retry jobs")
    p_mfit.add_argument("--run-dir", "-d", type=str, default=None, help="Override or specify the base MomentFit run directory directly")
    p_mfit.add_argument("--ess-threshold", type=float, default=None, help="Effective Sample Size (ESS) threshold for fit evaluation (default: from YAML or 20.0)")
    p_mfit.add_argument("--ess-metric", type=str, choices=["mean", "min"], default=None, help="ESS metric to evaluate against threshold: 'mean' or 'min' (default: from YAML or 'mean')")
    p_mfit.add_argument("--min-ess-threshold", type=float, default=None, help="Safety floor on worst-moment ESS even when --ess-metric=mean")
    p_mfit.add_argument("--plot", action="store_true", help="Run Stage 4 (plot unnormalized moments vs mass)")
    p_mfit.add_argument("--plot-type", choices=["all", "linear_algebra", "mcmc", "overlay"], default=None, help="Which unnormalized moment plots to generate: linear_algebra, mcmc, overlay, or all")
    p_mfit.add_argument("--plot-moments", nargs="+", default=None, help="Specific moment names to plot (e.g. 'H_0_0_0' 'H_0_2_0' or 'all')")
    p_mfit.add_argument("--plot-traces", action="store_true", help="Generate MCMC trace plots across all moments via MonitorMomentFit")
    p_mfit.add_argument("--plot-posteriors", action="store_true", help="Generate 1D posterior distributions across all moments via MonitorMomentFit")
    p_mfit.add_argument("--downsample-trace", type=int, nargs="?", const=300, default=None, help="Downsample MCMC trace plots to N points per moment (default: 300)")
    p_mfit.add_argument("--replot-from-root", action="store_true", help="Re-render all PDF plots directly from saved ROOT files without re-evaluating MCMC trees")
    p_mfit.add_argument("--plot-cov", "--plot-covariance", dest="plot_cov", action="store_true", help="Plot comparison of MCMC proposal covariance matrix before and after alignment (SPD regularization)")
    p_mfit.add_argument("--plot-cov-bins", nargs="+", default=None, help="Filter specific mass bin labels or centers to plot covariance for (default: all)")
    p_mfit.add_argument("--tbin", "-t", nargs="+", default=None, help="Filter specific t-bin(s) to run (e.g. 't010020' or '0')")
    p_mfit.add_argument("--max-l", "-l", nargs="+", type=int, default=None, help="Filter specific maxL multipole(s) to run (e.g. 4 or 4 6 8)")
    p_mfit.add_argument("--dry-run", action="store_true", help="Preview tasks and command lines without executing fits")
    p_mfit.add_argument("--validate", action="store_true", help="Validate configuration file and prerequisite paths, then exit")
    p_mfit.add_argument("--overwrite", action="store_true", help="Overwrite existing output directories without prompting")
    p_mfit.add_argument("--force-integral", action="store_true", help="Force recalculation of acceptance integral matrices in Stage 1")
    p_mfit.add_argument(
        "--calculate-acceptance",
        dest="calculate_acceptance",
        nargs="?",
        const=True,
        default=None,
        metavar="MACRO_PATH",
        help="Run acceptance calculation directly (optionally overriding the CalculateAcceptance macro path)",
    )
    p_mfit.add_argument("--skip-acceptance", action="store_true", help="Skip acceptance calculation after monitoring/fitting")
    p_mfit.add_argument("--force-acceptance", action="store_true", help="Force acceptance calculation even when not all bins meet convergence")

    raw_argv = sys.argv[1:]
    if len(raw_argv) >= 2 and raw_argv[0] in ("momentfit", "pipeline"):
        from .momentfit.runner import preprocess_momentfit_argv
        raw_argv = [raw_argv[0]] + preprocess_momentfit_argv(raw_argv[1:])

    args = parser.parse_args(raw_argv)
    if not args.command:
        parser.print_help()
        sys.exit(0)

    if args.command == "compile":
        sys.exit(cmd_compile(args))
    elif args.command in ("momentfit", "pipeline"):
        sys.exit(cmd_momentfit(args))
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


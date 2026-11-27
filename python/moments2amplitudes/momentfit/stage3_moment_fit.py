"""
Stage 3: MomentFit with BruFit MCMC using fixed covariance from linear algebra.
"""

from __future__ import annotations

import functools
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional, Tuple

from .config import PipelineConfig

print = functools.partial(print, flush=True)


def build_mcmc_and_diagnostics_flags(
    cfg: PipelineConfig,
    max_l: Optional[int] = None,
) -> List[str]:
    """Build CLI flags for all configurable MCMC, frame, covariance, and diagnostics parameters."""
    mcmc_cfg = cfg.moment_fit.mcmc
    diag_cfg = cfg.moment_fit.diagnostics
    resolved_gibbs = (
        mcmc_cfg.get_gibbs_block_size(max_l)
        if max_l is not None
        else int(mcmc_cfg.gibbs_block_size)
    )
    resolved_target_acc = mcmc_cfg.get_target_acceptance()

    flags = [
        "--gibbs-block-size", str(resolved_gibbs),
        "--cold-warmup-acc", str(mcmc_cfg.cold_warmup_acc),
        "--min-acceptance", str(mcmc_cfg.min_acceptance),
        "--max-acceptance", str(mcmc_cfg.max_acceptance),
        "--target-acceptance", str(resolved_target_acc),
        "--tune-mode", str(mcmc_cfg.tune_mode),
        "--rhat-window", str(mcmc_cfg.rhat_window),
        "--rhat-target", str(mcmc_cfg.rhat_target),
        "--ess-fraction", str(mcmc_cfg.ess_fraction),
        "--rhat-max-retries", str(mcmc_cfg.rhat_max_retries),
        "--bailout-threshold", str(mcmc_cfg.bailout_threshold),
        "--bailout-min-acceptance", str(mcmc_cfg.bailout_min_acceptance),
        "--mc-integration-bins", str(mcmc_cfg.mc_integration_bins),
        "--num-cpu", str(mcmc_cfg.num_cpu),
        "--seed", str(mcmc_cfg.seed),
        "--mc-variance" if mcmc_cfg.apply_mc_variance else "--no-mc-variance",
        "--negative-pdf" if mcmc_cfg.allow_negative_pdf else "--no-negative-pdf",
        "--tuning" if mcmc_cfg.use_tuning else "--no-tuning",
        "--fixed-covariance" if mcmc_cfg.fixed_covariance else "--no-fixed-covariance",
        "--use-covariance" if mcmc_cfg.use_covariance else "--no-use-covariance",
        "--set-init-value" if mcmc_cfg.set_init_value else "--no-set-init-value",
        "--gj-frame" if mcmc_cfg.use_gj_frame else "--hel-frame",
        "--max-retry-iters", str(diag_cfg.max_retry_iters),
        "--burnin-fraction", str(diag_cfg.burnin_fraction),
        "--low-acceptance-threshold", str(diag_cfg.low_acceptance_threshold),
        "--gibbs-low-acceptance-threshold", str(diag_cfg.gibbs_low_acceptance_threshold),
        "--high-acceptance-threshold", str(diag_cfg.high_acceptance_threshold),
        "--gibbs-activation-acceptance", str(diag_cfg.gibbs_activation_acceptance),
        "--recondition-covariance" if diag_cfg.recondition_covariance else "--no-recondition-covariance",
        "--recondition-shrinkage", str(diag_cfg.recondition_shrinkage),
        "--recondition-min-eigen-floor", str(diag_cfg.recondition_min_eigen_floor),
    ]
    return flags


TBIN_INDEX_MAP = {
    "t010020": 0,
    "t020032": 1,
    "t032050": 2,
    "t050075": 3,
    "t075100": 4,
}


def build_moment_fit_command(
    cfg: PipelineConfig,
    max_l: int,
    tbin_arg: str,
) -> Tuple[List[str], str]:
    """
    Constructs the command-line invocation for the configured MomentFit run script (RunMomentFit.sh).
    """
    m_min = cfg.analysis.mass_binning.min
    m_max = cfg.analysis.mass_binning.max
    n_bins = cfg.analysis.mass_binning.n_bins
    run_dir = cfg.moment_fit.get_run_dir(max_l, m_min, m_max)

    cmd = [
        "bash",
        cfg.moment_fit.run_script,
        "--run-dir", run_dir,
        "--cov-dir", cfg.linear_algebra.output_base_dir,
        "--tbin", tbin_arg,
        "--max-l", str(max_l),
        "--m-min", str(m_min),
        "--m-max", str(m_max),
        "--n-bins", str(n_bins),
        "--init-scale", str(cfg.moment_fit.mcmc.init_scale),
        "--step", str(cfg.moment_fit.mcmc.step),
        "--iters", str(cfg.moment_fit.mcmc.iterations),
        "--burn-in", str(cfg.moment_fit.mcmc.burnin),
        "--ess-threshold", str(cfg.moment_fit.diagnostics.ess_threshold),
        "--ess-metric", str(cfg.moment_fit.diagnostics.ess_metric),
        "--min-ess-threshold", str(
            cfg.moment_fit.diagnostics.min_ess_threshold
            if getattr(cfg.moment_fit.diagnostics, "min_ess_threshold", None) is not None
            else 0
        ),
        "--max-parallel", str(cfg.moment_fit.max_parallel),
    ]
    cmd.extend(build_mcmc_and_diagnostics_flags(cfg, max_l=max_l))

    if cfg.moment_fit.diagnostics.auto_retry:
        cmd.append("--auto-retry")
        cmd.extend(["--max-retries", str(cfg.moment_fit.diagnostics.max_retries)])

    if cfg.moment_fit.overwrite:
        cmd.append("--overwrite")

    if getattr(cfg.moment_fit, "skip_acceptance", False):
        cmd.append("--skip-acceptance")
    if getattr(cfg.moment_fit, "force_acceptance", False):
        cmd.append("--force-acceptance")

    # Pass configured input data directory and patterns if available
    if hasattr(cfg.moment_fit, "data"):
        data_cfg = cfg.moment_fit.data
        if data_cfg.base_dir:
            cmd.extend(["--data-dir", data_cfg.base_dir])
        if data_cfg.data:
            cmd.extend(["--data-pattern", data_cfg.data])
        if data_cfg.mcrecon:
            cmd.extend(["--mcrecon-pattern", data_cfg.mcrecon])
        if data_cfg.mcgen:
            cmd.extend(["--mcgen-pattern", data_cfg.mcgen])
        if data_cfg.weight:
            cmd.extend(["--weight-pattern", data_cfg.weight])

    return cmd, run_dir


def check_stage3_data_files(
    cfg: PipelineConfig,
    tbins: List[str],
) -> List[str]:
    """
    Checks that required input data, MC, and weight files exist for Stage 3.
    """
    missing = []
    if not hasattr(cfg.moment_fit, "data"):
        return missing
    data_cfg = cfg.moment_fit.data
    for tbin in tbins:
        data_f = data_cfg.get_data_path(tbin)
        recon_f = data_cfg.get_mcrecon_path(tbin)
        weight_f = data_cfg.get_weight_path(tbin)
        if not os.path.isfile(data_f):
            missing.append(data_f)
        if not os.path.isfile(recon_f):
            missing.append(recon_f)
        if not os.path.isfile(weight_f):
            missing.append(weight_f)
        if data_cfg.use_generated_mc:
            gen_f = data_cfg.get_mcgen_path(tbin)
            if not os.path.isfile(gen_f):
                missing.append(gen_f)
    return missing


def check_stage2_prerequisites(
    cfg: PipelineConfig,
    tbins: List[str],
    max_l: int,
) -> List[str]:
    """
    Verifies that Stage 2 exported covariance and initial value files exist.
    """
    missing = []
    bin_cfg = cfg.analysis.mass_binning
    for tbin in tbins:
        out_dir = cfg.linear_algebra.get_output_dir(tbin, max_l, cfg.analysis.beam_pol_label)
        for bin_idx in range(bin_cfg.n_bins):
            center = bin_cfg.get_bin_center(bin_idx)
            lbl = f"{bin_cfg.variable}_{center:.{bin_cfg.digits}f}"
            cov_f = os.path.join(out_dir, f"BruCovMatrix_{lbl}.root")
            init_f = os.path.join(out_dir, f"covariance_matrix_{lbl}.root")
            if not os.path.isfile(cov_f):
                missing.append(cov_f)
            if not os.path.isfile(init_f):
                missing.append(init_f)
    return missing


def run_stage3(
    cfg: PipelineConfig,
    tbins_filter: Optional[List[str]] = None,
    max_l_filter: Optional[List[int]] = None,
    dry_run: bool = False,
) -> bool:
    """
    Executes Stage 3: MomentFit MCMC with BruFit across target t-bins and max_l multipoles.
    """
    tbins = tbins_filter if tbins_filter is not None else cfg.analysis.tbins
    max_ls = max_l_filter if max_l_filter is not None else cfg.analysis.max_l

    print("\n" + "=" * 70)
    print(" STAGE 3: MomentFit with BruFit MCMC (Fixed Linear Algebra Covariance)")
    print("=" * 70)
    print(f" Target t-bins:      {', '.join(tbins)}")
    print(f" Multipoles Max L:   {max_ls}")
    print(f" MCMC Iterations:    {cfg.moment_fit.mcmc.iterations} (burn-in: {cfg.moment_fit.mcmc.burnin})")
    print(f" Init Scale / Step:  {cfg.moment_fit.mcmc.init_scale} / {cfg.moment_fit.mcmc.step}")
    print(f" Tuning / Neg PDF:   {cfg.moment_fit.mcmc.use_tuning} ({cfg.moment_fit.mcmc.tune_mode}) / {cfg.moment_fit.mcmc.allow_negative_pdf}")
    print(f" Fixed Covariance:   {cfg.moment_fit.mcmc.fixed_covariance} (Recondition retries: {cfg.moment_fit.diagnostics.recondition_covariance})")
    print(f" Auto-Retry Loop:    {cfg.moment_fit.diagnostics.auto_retry} (ESS >= {cfg.moment_fit.diagnostics.ess_threshold}, max_iters={cfg.moment_fit.diagnostics.max_retry_iters})")
    print(f" Script:             {cfg.moment_fit.run_script}")
    print("=" * 70 + "\n")

    # Map t-bins to argument format for RunMomentFit.sh
    if set(tbins) == set(TBIN_INDEX_MAP.keys()):
        tbin_arg = "all"
    elif len(tbins) == 1 and tbins[0] in TBIN_INDEX_MAP:
        tbin_arg = str(TBIN_INDEX_MAP[tbins[0]])
    else:
        # Pass mapped numeric indices joined by space or names
        indices = [str(TBIN_INDEX_MAP[t]) for t in tbins if t in TBIN_INDEX_MAP]
        tbin_arg = " ".join(indices) if indices else "all"

    all_success = True

    for l in max_ls:
        cmd, run_dir = build_moment_fit_command(cfg, l, tbin_arg)

        print(f"\n[EXECUTION PLAN] MomentFit for Max L = {l}:")
        print(f"  Run Directory: {run_dir}")
        print(f"  Command:       {' '.join(cmd)}")

        # Verify Stage 2 prerequisites
        if not dry_run:
            if cfg.moment_fit.mcmc.fixed_covariance or cfg.moment_fit.mcmc.set_init_value:
                missing = check_stage2_prerequisites(cfg, tbins, l)
                if missing:
                    print(f"❌ Error: Missing {len(missing)} covariance prerequisite files for Max L = {l}:")
                    for m in missing[:5]:
                        print(f"  - {m}")
                    if len(missing) > 5:
                        print(f"  ... and {len(missing) - 5} more.")
                    print("  Please run Stage 1 and Stage 2 first!")
                    all_success = False
                    continue

            # Verify Stage 3 input data files
            missing_data = check_stage3_data_files(cfg, tbins)
            if missing_data:
                print(f"⚠️  Warning: {len(missing_data)} input data/MC file(s) for Stage 3 were not found on disk:")
                for m in missing_data[:5]:
                    print(f"  - {m}")
                if len(missing_data) > 5:
                    print(f"  ... and {len(missing_data) - 5} more.")

        if dry_run:
            print("[DRY-RUN] Skipping actual fit execution.")
            continue

        run_script_name = os.path.basename(cfg.moment_fit.run_script)
        print(f"\n--- Launching {run_script_name} for Max L = {l} ---")
        start_time = time.time()

        # Run with environment from set_ana_env.sh
        env = os.environ.copy()
        if "BRUFIT" not in env:
            env["BRUFIT"] = cfg.environment.brufit_dir
        if getattr(cfg.moment_fit, "calculate_acceptance", ""):
            env["MOM2AMP_CALCULATE_ACCEPTANCE"] = str(cfg.moment_fit.calculate_acceptance)

        ret = subprocess.run(cmd, env=env)
        elapsed = time.time() - start_time

        if ret.returncode != 0:
            print(f"❌ MomentFit failed for Max L = {l} (exit code: {ret.returncode}, elapsed: {elapsed:.1f}s)")
            all_success = False
        else:
            print(f"✅ MomentFit completed successfully for Max L = {l} (elapsed: {elapsed:.1f}s)")

    print("\n" + "-" * 70)
    if all_success:
        print("✅ Stage 3 completed successfully for all multipoles!")
    else:
        print("❌ Stage 3 encountered failures.")
    return all_success

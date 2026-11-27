"""
Unit tests for MomentFit workflow (Linear Algebra -> Covariance Export -> BruFit MomentFit).
"""

from __future__ import annotations

import os
import unittest
from pathlib import Path

from moments2amplitudes.momentfit.config import (
    MomentFitConfig,
    load_momentfit_config,
    validate_momentfit_config,
)
from moments2amplitudes.momentfit.runner import normalize_stage_names, normalize_tbin_filter
from moments2amplitudes.momentfit.stage3_moment_fit import (
    TBIN_INDEX_MAP,
    build_moment_fit_command,
    check_stage2_prerequisites,
)


class TestMomentFit(unittest.TestCase):
    def setUp(self):
        self.config_path = Path("configs/momentfit_nominal.yaml")

    def test_load_and_validate_nominal_momentfit_config(self):
        self.assertTrue(self.config_path.is_file(), f"Config file not found: {self.config_path}")
        cfg = load_momentfit_config(self.config_path)
        self.assertIsInstance(cfg, MomentFitConfig)

        # Check stages
        self.assertEqual(cfg.stages, ["linear_algebra", "export_covariance", "moment_fit", "plot"])

        # Check analysis settings
        self.assertEqual(cfg.analysis.subsystem, "EtaPi0")
        self.assertEqual(cfg.analysis.frame, "GJ")
        self.assertIn("t010020", cfg.analysis.tbins)
        self.assertIn(4, cfg.analysis.max_l)

        # Check mass binning
        self.assertEqual(cfg.analysis.mass_binning.n_bins, 30)
        self.assertAlmostEqual(cfg.analysis.mass_binning.min, 0.8)
        self.assertAlmostEqual(cfg.analysis.mass_binning.max, 2.0)
        self.assertAlmostEqual(cfg.analysis.mass_binning.bin_width, 0.04)

        # Check weight branch configuration
        self.assertEqual(cfg.input_data.weight_branch, "weightASBS")
        self.assertEqual(cfg.linear_algebra.weight_branch, "weightASBS")
        self.assertNotIn("eventWeight", cfg.input_data.columns)

        # Check moment_fit input data configuration
        data_cfg = cfg.moment_fit.data
        self.assertTrue(data_cfg.base_dir.endswith("nominal_tbins"))
        self.assertIn("{tbin}", data_cfg.data)
        self.assertIn("{tbin}", data_cfg.mcrecon)
        self.assertIn("{tbin}", data_cfg.mcgen)
        self.assertIn("{tbin}", data_cfg.weight)
        self.assertEqual(data_cfg.tree_name, "kin")
        self.assertEqual(data_cfg.weight_species, "Acc")
        self.assertEqual(data_cfg.weight_name, "weight")
        self.assertTrue(data_cfg.use_generated_mc)

        # Check path resolution for tbin
        data_p = data_cfg.get_data_path("t010020")
        mcrecon_p = data_cfg.get_mcrecon_path("t010020")
        mcgen_p = data_cfg.get_mcgen_path("t010020")
        weight_p = data_cfg.get_weight_path("t010020")
        self.assertTrue(data_p.endswith("data_merged_t010020.root"))
        self.assertTrue(mcrecon_p.endswith("reconmc_merged_t010020.root"))
        self.assertTrue(mcgen_p.endswith("genmc_merged_t010020.root"))
        self.assertTrue(weight_p.endswith("weights_data_merged_t010020.root"))

        # Check plotting configuration
        self.assertTrue(cfg.plotting.enabled)
        self.assertTrue(cfg.plotting.plot_linear_algebra)
        self.assertTrue(cfg.plotting.plot_mcmc)
        self.assertTrue(cfg.plotting.plot_overlay)
        self.assertEqual(cfg.plotting.moments, ["all"])

        # Check validation (syntax and configuration values)
        errors = validate_momentfit_config(cfg, check_files=False)
        self.assertEqual(errors, [], f"Validation failed with errors: {errors}")

    def test_stage_normalization(self):
        self.assertEqual(
            normalize_stage_names(["1", "2", "3", "4"]),
            ["linear_algebra", "export_covariance", "moment_fit", "plot"],
        )
        self.assertEqual(
            normalize_stage_names(["linear_algebra", "plot"]),
            ["linear_algebra", "plot"],
        )
        self.assertEqual(
            normalize_stage_names(["monitor", "collect"]),
            ["monitor", "collect"],
        )
        self.assertEqual(
            normalize_stage_names(["status", "diagnostics"]),
            ["monitor", "collect"],
        )
        self.assertEqual(
            normalize_stage_names(["all"]),
            ["linear_algebra", "export_covariance", "moment_fit", "plot"],
        )
        self.assertEqual(normalize_stage_names(["3", "collect"]), ["collect"])
        self.assertEqual(normalize_stage_names(["moment_fit", "collect"]), ["collect"])
        self.assertEqual(normalize_stage_names(["3", "monitor"]), ["monitor"])
        self.assertEqual(normalize_stage_names(["3", "monitor", "collect"]), ["monitor", "collect"])
        with self.assertRaises(ValueError):
            normalize_stage_names(["invalid_stage"])

    def test_stage_3_collect_skips_fitting(self):
        from unittest.mock import patch
        from moments2amplitudes.momentfit.runner import parse_args, run_momentfit

        args = parse_args(["configs/momentfit_test.yaml", "--stage", "3", "collect"])
        self.assertEqual(args.stage, ["3", "collect"])
        self.assertEqual(normalize_stage_names(args.stage), ["collect"])

        with patch("moments2amplitudes.momentfit.stage3_moment_fit.run_stage3") as mock_stage3, \
             patch("moments2amplitudes.momentfit.runner.run_collector", return_value=True) as mock_collect:
            rc = run_momentfit(config_path=args.config, stages=args.stage, dry_run=True)

        self.assertEqual(rc, 0)
        mock_stage3.assert_not_called()
        mock_collect.assert_called_once()

    def test_tbin_filter_normalization(self):
        self.assertIsNone(normalize_tbin_filter(None))
        self.assertIsNone(normalize_tbin_filter(["all"]))
        self.assertEqual(normalize_tbin_filter(["t010020", "t020032"]), ["t010020", "t020032"])
        self.assertEqual(normalize_tbin_filter(["0", "1"]), ["t010020", "t020032"])

    def test_tbin_index_mapping(self):
        self.assertEqual(TBIN_INDEX_MAP["t010020"], 0)
        self.assertEqual(TBIN_INDEX_MAP["t020032"], 1)
        self.assertEqual(TBIN_INDEX_MAP["t032050"], 2)
        self.assertEqual(TBIN_INDEX_MAP["t050075"], 3)
        self.assertEqual(TBIN_INDEX_MAP["t075100"], 4)

    def test_build_moment_fit_command(self):
        cfg = load_momentfit_config(self.config_path)
        cmd, run_dir = build_moment_fit_command(cfg, max_l=4, tbin_arg="0")

        self.assertEqual(cmd[0], "bash")
        self.assertTrue(cmd[1].endswith("RunMomentFit.sh"))
        self.assertIn("--cov-dir", cmd)
        self.assertIn("--max-l", cmd)
        self.assertIn("4", cmd)
        self.assertIn("--tbin", cmd)
        self.assertIn("0", cmd)
        self.assertIn("--data-dir", cmd)
        self.assertIn("--data-pattern", cmd)
        self.assertIn("--mcrecon-pattern", cmd)
        self.assertIn("--mcgen-pattern", cmd)
        self.assertIn("--weight-pattern", cmd)
        self.assertTrue("MomentFit_m080200_FixedCov_Lmax4" in run_dir)

    def test_brufit_data_config_parsing_and_aliases(self):
        from moments2amplitudes.momentfit.config import _parse_brufit_data_config, BruFitDataConfig

        # Test dictionary with alternate aliases
        custom_dict = {
            "data_dir": "/custom/path/data",
            "data_file": "{data_dir}/custom_data_{tbin}.root",
            "reconmc": "{data_dir}/custom_recon_{tbin}.root",
            "genmc": "{data_dir}/custom_gen_{tbin}.root",
            "weights": "{data_dir}/custom_weight_{tbin}.root",
            "tree": "custom_kin",
            "species": "Sig",
            "weight_name": "w",
            "use_generated_mc": False,
        }
        data_cfg = _parse_brufit_data_config(custom_dict)
        self.assertEqual(data_cfg.base_dir, "/custom/path/data")
        self.assertEqual(data_cfg.data, "{data_dir}/custom_data_{tbin}.root")
        self.assertEqual(data_cfg.mcrecon, "{data_dir}/custom_recon_{tbin}.root")
        self.assertEqual(data_cfg.mcgen, "{data_dir}/custom_gen_{tbin}.root")
        self.assertEqual(data_cfg.weight, "{data_dir}/custom_weight_{tbin}.root")
        self.assertEqual(data_cfg.tree_name, "custom_kin")
        self.assertEqual(data_cfg.weight_species, "Sig")
        self.assertEqual(data_cfg.weight_name, "w")
        self.assertFalse(data_cfg.use_generated_mc)

        # Test path resolution
        self.assertEqual(
            data_cfg.get_data_path("t020032"),
            "/custom/path/data/custom_data_t020032.root",
        )
        self.assertEqual(
            data_cfg.get_mcrecon_path("t020032"),
            "/custom/path/data/custom_recon_t020032.root",
        )

    def test_parse_args_requires_config(self):
        from moments2amplitudes.momentfit.runner import parse_args
        # Should exit with error if no config and no run-dir provided
        with self.assertRaises(SystemExit):
            parse_args([])

        # With positional config
        args = parse_args(["configs/momentfit_nominal.yaml"])
        self.assertEqual(args.config, "configs/momentfit_nominal.yaml")

        # With -c flag
        args_flag = parse_args(["-c", "configs/momentfit_nominal.yaml"])
        self.assertEqual(args_flag.config, "configs/momentfit_nominal.yaml")

        # Monitoring with --run-dir does not require config file
        args_mon = parse_args(["--run-dir", "/tmp/run_dir", "--monitor"])
        self.assertIsNone(args_mon.config)
        self.assertEqual(args_mon.run_dir, "/tmp/run_dir")
        self.assertTrue(args_mon.monitor)
        self.assertEqual(args_mon.stage, ["monitor"])

        # Collecting with --auto-retry and config
        args_col = parse_args(["configs/momentfit_nominal.yaml", "--collect", "--auto-retry"])
        self.assertEqual(args_col.config, "configs/momentfit_nominal.yaml")
        self.assertTrue(args_col.collect)
        self.assertTrue(args_col.auto_retry)
        self.assertEqual(args_col.stage, ["collect"])

    def test_plotting_helpers_and_single_moment_overlay(self):
        import tempfile
        import numpy as np
        from moments2amplitudes.momentfit.stage4_plot_moments import (
            parse_moment_tuple,
            canonical_moment_name,
            format_moment_label,
            plot_single_moment,
        )

        self.assertEqual(parse_moment_tuple("H_0_2_0"), (0, 2, 0))
        self.assertEqual(parse_moment_tuple("H2_1_1"), (2, 1, 1))
        self.assertEqual(canonical_moment_name("gr_unnormalized_moments_0_2_0_vs_mass"), "H_0_2_0")
        self.assertIn("Re", format_moment_label("H_0_2_0"))
        self.assertIn("Im", format_moment_label("H_2_1_1"))

        # Test plot_single_moment with synthetic LA and MCMC data
        with tempfile.TemporaryDirectory() as tmpdir:
            la_data = {
                "mass": np.array([0.82, 0.86, 0.90]),
                "val": np.array([1000.0, 1500.0, 2000.0]),
                "err": np.array([50.0, 60.0, 70.0]),
            }
            mcmc_data = {
                "mass": np.array([0.82, 0.86, 0.90]),
                "val": np.array([980.0, 1520.0, 1990.0]),
                "err": np.array([45.0, 55.0, 65.0]),
            }

            # 1. Linear algebra individually
            la_plots = plot_single_moment(
                moment_name="H_0_0_0",
                tbin="t010020",
                max_l=4,
                out_dir=os.path.join(tmpdir, "la"),
                la_data=la_data,
                mcmc_data=None,
                mode="linear_algebra",
                formats=("png", "pdf"),
            )
            self.assertEqual(len(la_plots), 2)
            for p in la_plots:
                self.assertTrue(p.is_file(), f"Expected plot file not created: {p}")

            # 2. MCMC individually
            mcmc_plots = plot_single_moment(
                moment_name="H_0_0_0",
                tbin="t010020",
                max_l=4,
                out_dir=os.path.join(tmpdir, "mcmc"),
                la_data=None,
                mcmc_data=mcmc_data,
                mode="mcmc",
                formats=("png", "pdf"),
            )
            self.assertEqual(len(mcmc_plots), 2)
            for p in mcmc_plots:
                self.assertTrue(p.is_file(), f"Expected plot file not created: {p}")

            # 3. Overlay (on top of each other)
            overlay_plots = plot_single_moment(
                moment_name="H_0_0_0",
                tbin="t010020",
                max_l=4,
                out_dir=os.path.join(tmpdir, "overlay"),
                la_data=la_data,
                mcmc_data=mcmc_data,
                mode="overlay",
                formats=("png", "pdf"),
            )
            self.assertEqual(len(overlay_plots), 2)
            for p in overlay_plots:
                self.assertTrue(p.is_file(), f"Expected plot file not created: {p}")

    def test_parse_args_plot_cov(self):
        from moments2amplitudes.momentfit.runner import parse_args, normalize_stage_names

        args = parse_args(["configs/momentfit_nominal.yaml", "--plot-cov", "--plot-cov-bins", "1.300", "1.340"])
        self.assertTrue(args.plot_cov)
        self.assertEqual(args.plot_cov_bins, ["1.300", "1.340"])
        self.assertEqual(args.stage, ["plot_covariance"])

        self.assertEqual(
            normalize_stage_names(["plot_cov", "cov", "plot_covariance"]),
            ["plot_covariance"],
        )

    def test_covariance_plotting_and_spd_diagnostics(self):
        import tempfile
        from moments2amplitudes.momentfit.plot_covariance import (
            cov_to_correlation,
            check_positive_definite,
            plot_covariance_comparison_single_bin,
            plot_covariance_for_directory,
        )

        # Check existing nominal covariance file pair if available
        cov_dir = "Moments.out/nominal/t010020/All.maxL_8"
        raw_p = os.path.join(cov_dir, "covariance_matrix_Mpi0eta_1.300.root")
        aligned_p = os.path.join(cov_dir, "BruCovMatrix_Mpi0eta_1.300.root")

        if os.path.isfile(raw_p) and os.path.isfile(aligned_p):
            with tempfile.TemporaryDirectory() as tmpdir:
                res = plot_covariance_comparison_single_bin(
                    raw_path=raw_p,
                    aligned_path=aligned_p,
                    output_dir=tmpdir,
                    mass_bin_label="Mpi0eta_1.300",
                    tbin="t010020",
                    max_l=8,
                    formats=["png", "pdf"],
                )
                self.assertIsNotNone(res)
                self.assertFalse(res["raw_is_spd"], "Raw covariance was expected to be indefinite")
                self.assertTrue(res["aligned_is_spd"], "Aligned covariance must be strictly SPD")
                self.assertGreater(res["aligned_min_eigenval"], 0.0)
                self.assertEqual(len(res["saved_plots"]), 2)
                for p in res["saved_plots"]:
                    self.assertTrue(os.path.isfile(p), f"Plot file not created: {p}")

                # Test plot_covariance_for_directory with filter
                dir_plots = plot_covariance_for_directory(
                    cov_dir=cov_dir,
                    output_dir=tmpdir,
                    mass_bins_filter=["1.300"],
                    tbin="t010020",
                    max_l=8,
                    formats=["png"],
                )
                self.assertGreaterEqual(len(dir_plots), 1)

    def test_vectorized_chain_diagnostics_and_collector_cache(self):
        import importlib.util
        import tempfile
        import numpy as np
        import uproot
        from unittest.mock import patch
        import BruFitResults as BFR

        spec = importlib.util.spec_from_file_location("MomentFitCollector", "macros/MomentFitCollector.py")
        mfc = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mfc)

        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            tbin_dir = run_dir / "t010020"
            bin0_dir = tbin_dir / "Mpi0eta0.820000_"
            bin0_dir.mkdir(parents=True)
            root_path = bin0_dir / "ResultsBruMcmcFixedCovariance.root"

            rng = np.random.default_rng(123)
            n_entries = 1200
            weights = rng.integers(1, 5, size=n_entries).astype(np.float64)
            h_0_0_0 = np.cumsum(rng.normal(0.0, 0.2, size=n_entries)).astype(np.float64)
            h_0_2_0 = rng.normal(1.5, 0.3, size=n_entries).astype(np.float64)
            h_0_2_0_err = rng.uniform(0.01, 0.05, size=n_entries).astype(np.float64)

            with uproot.recreate(root_path) as uf:
                uf["MCMCTree"] = {
                    "weight": weights,
                    "H_0_0_0": h_0_0_0,
                    "H_0_2_0": h_0_2_0,
                    "H_0_2_0_err": h_0_2_0_err,
                }
            self.assertGreater(root_path.stat().st_size, 10240)

            # Compare vectorized uproot diagnostics against reference scalar ACF formula
            bfr_diag = BFR.compute_chain_diagnostics_from_root(str(root_path), burn_in=200)
            mfc_diag = mfc.compute_bin_diagnostics_from_root(str(root_path), burn_in=200)
            self.assertIsNotNone(bfr_diag)
            self.assertIsNotNone(mfc_diag)
            self.assertEqual(set(bfr_diag["per_moment"].keys()), {"H_0_0_0", "H_0_2_0"})

            int_w = np.maximum(1, np.round(weights)).astype(int)
            for name, raw_arr in [("H_0_0_0", h_0_0_0), ("H_0_2_0", h_0_2_0)]:
                d = np.repeat(raw_arr, int_w)[200:]
                acf = BFR._compute_autocorrelation_fft(d, max_lag=len(d) // 2)
                tau_ref = 0.5
                for lag in range(1, len(acf)):
                    if acf[lag] <= 0:
                        break
                    tau_ref += acf[lag]
                ess_ref = min(float(len(d)), len(d) / (2.0 * tau_ref))
                self.assertAlmostEqual(bfr_diag["per_moment"][name]["tau"], tau_ref, places=10)
                self.assertAlmostEqual(bfr_diag["per_moment"][name]["ess"], ess_ref, places=10)
                self.assertAlmostEqual(mfc_diag["per_moment"][name]["tau"], tau_ref, places=10)
                self.assertAlmostEqual(mfc_diag["per_moment"][name]["ess"], ess_ref, places=10)

            # Evaluate tbin populates .collector_diag_cache.json
            evals1 = mfc.evaluate_tbin(
                tbin_idx=0, tbin_name="t010020", tbin_dir=str(tbin_dir), run_dir=str(run_dir),
                m_min=0.8, m_max=2.0, n_bins=30, ess_threshold=20.0, burn_in=200,
            )
            self.assertEqual(len(evals1), 30)
            self.assertEqual(evals1[0]["status"], "SUCCESS")
            self.assertAlmostEqual(evals1[0]["mean_ess"], mfc_diag["ess_mean"], places=10)
            cache_file = tbin_dir / ".collector_diag_cache.json"
            self.assertTrue(cache_file.is_file())

            # Second call with same ROOT file hits cache without invoking _inspect_root_bin_task
            with patch.object(mfc, "_inspect_root_bin_task", side_effect=AssertionError("Should hit cache")) as mock_insp:
                evals2 = mfc.evaluate_tbin(
                    tbin_idx=0, tbin_name="t010020", tbin_dir=str(tbin_dir), run_dir=str(run_dir),
                    m_min=0.8, m_max=2.0, n_bins=30, ess_threshold=100000.0, burn_in=200,
                )
                mock_insp.assert_not_called()
            # Threshold recomputed dynamically from cached per_moment_ess
            self.assertFalse(evals2[0]["is_successful"])
            self.assertEqual(evals2[0]["n_low_ess"], 2)

            # Single-job fast path evaluates only requested bin_idx
            evals_single = mfc.evaluate_tbin(
                tbin_idx=0, tbin_name="t010020", tbin_dir=str(tbin_dir), run_dir=str(run_dir),
                m_min=0.8, m_max=2.0, n_bins=30, ess_threshold=20.0, burn_in=200, target_bin_idx=0,
            )
            self.assertEqual(len(evals_single), 1)
            self.assertEqual(evals_single[0]["bin_idx"], 0)

            # End-to-end CLI subprocess run exercising ProcessPoolExecutor (> 2 uncached bins)
            import subprocess, sys, json
            for b_i, b_center in [(1, 0.86), (2, 0.90), (3, 0.94)]:
                bdir = tbin_dir / f"Mpi0eta{b_center:.6f}_"
                bdir.mkdir(parents=True)
                with uproot.recreate(bdir / "ResultsBruMcmcFixedCovariance.root") as uf:
                    uf["MCMCTree"] = {
                        "weight": weights,
                        "H_0_0_0": h_0_0_0,
                        "H_0_2_0": h_0_2_0,
                    }
            summary_json = run_dir / "summary_4bins.json"
            proc = subprocess.run(
                [
                    sys.executable, "macros/MomentFitCollector.py",
                    "--run-dir", str(run_dir), "--tbin", "0",
                    "--m-min", "0.8", "--m-max", "0.96", "--n-bins", "4",
                    "--burn-in", "200", "--ess-threshold", "1.0",
                    "--non-interactive", "--no-exec",
                    "--summary-json", str(summary_json),
                ],
                capture_output=True, text=True, check=True,
            )
            self.assertEqual(proc.returncode, 0)
            summary = json.loads(summary_json.read_text())
            self.assertEqual(summary["total_bins"], 4)
            self.assertEqual(summary["successful_count"], 4)

    def test_mcmc_and_diagnostics_config_yaml_and_validation(self):
        import tempfile
        import yaml

        for yaml_path in ["configs/momentfit_nominal.yaml", "configs/momentfit_test.yaml"]:
            cfg = load_momentfit_config(yaml_path)
            self.assertTrue(cfg.moment_fit.mcmc.fixed_covariance)
            self.assertTrue(cfg.moment_fit.mcmc.use_covariance)
            self.assertTrue(cfg.moment_fit.mcmc.set_init_value)
            self.assertTrue(cfg.moment_fit.mcmc.use_gj_frame)
            self.assertEqual(cfg.moment_fit.mcmc.min_acceptance, 0.05)
            self.assertEqual(cfg.moment_fit.mcmc.max_acceptance, 0.15)
            self.assertEqual(cfg.moment_fit.mcmc.get_target_acceptance(), 0.10)
            self.assertEqual(cfg.moment_fit.mcmc.tune_mode, "MappedRhat")
            self.assertEqual(cfg.moment_fit.mcmc.mc_integration_bins, 50)
            self.assertEqual(cfg.moment_fit.diagnostics.max_retry_iters, 200000)
            self.assertAlmostEqual(cfg.moment_fit.diagnostics.burnin_fraction, 0.10)
            self.assertAlmostEqual(cfg.moment_fit.diagnostics.low_acceptance_threshold, 0.03)
            self.assertAlmostEqual(cfg.moment_fit.diagnostics.gibbs_low_acceptance_threshold, 0.12)
            self.assertAlmostEqual(cfg.moment_fit.diagnostics.high_acceptance_threshold, 0.35)
            self.assertFalse(cfg.moment_fit.diagnostics.recondition_covariance)
            self.assertEqual(validate_momentfit_config(cfg, check_files=False), [])

        with tempfile.TemporaryDirectory() as tmpdir:
            custom_yaml = Path(tmpdir) / "custom.yaml"
            raw = yaml.safe_load(Path("configs/momentfit_test.yaml").read_text())
            raw["moment_fit"]["mcmc"]["gibbs_block_size"] = 0
            raw["moment_fit"]["mcmc"]["gibbs_block_size_by_max_l"] = {"6": 12, 8: 16}
            raw["moment_fit"]["mcmc"]["fixed_covariance"] = False
            raw["moment_fit"]["diagnostics"]["recondition_covariance"] = True
            raw["moment_fit"]["diagnostics"]["recondition_shrinkage"] = 0.25
            custom_yaml.write_text(yaml.safe_dump(raw))

            cfg_custom = load_momentfit_config(custom_yaml)
            self.assertEqual(cfg_custom.moment_fit.mcmc.get_gibbs_block_size(4), 0)
            self.assertEqual(cfg_custom.moment_fit.mcmc.get_gibbs_block_size(6), 12)
            self.assertEqual(cfg_custom.moment_fit.mcmc.get_gibbs_block_size(8), 16)
            self.assertFalse(cfg_custom.moment_fit.mcmc.fixed_covariance)
            self.assertTrue(cfg_custom.moment_fit.diagnostics.recondition_covariance)
            self.assertEqual(validate_momentfit_config(cfg_custom, check_files=False), [])

            # Invalid bounds caught by validation
            cfg_bad = load_momentfit_config(custom_yaml)
            cfg_bad.moment_fit.mcmc.min_acceptance = 0.3
            cfg_bad.moment_fit.mcmc.max_acceptance = 0.2
            cfg_bad.moment_fit.diagnostics.burnin_fraction = 0.6
            cfg_bad.moment_fit.diagnostics.recondition_shrinkage = 1.5
            errs = validate_momentfit_config(cfg_bad, check_files=False)
            self.assertGreaterEqual(len(errs), 3)

    def test_stage3_and_collector_forwarding_with_gibbs_overrides(self):
        from unittest.mock import patch, MagicMock
        from moments2amplitudes.momentfit.runner import run_collector

        cfg = load_momentfit_config("configs/momentfit_test.yaml")
        cfg.moment_fit.mcmc.gibbs_block_size = 0
        cfg.moment_fit.mcmc.gibbs_block_size_by_max_l = {8: 16}
        cfg.moment_fit.mcmc.fixed_covariance = False
        cfg.moment_fit.diagnostics.recondition_covariance = True

        cmd_l4, _ = build_moment_fit_command(cfg, max_l=4, tbin_arg="0")
        cmd_l8, _ = build_moment_fit_command(cfg, max_l=8, tbin_arg="0")
        self.assertEqual(cmd_l4[cmd_l4.index("--gibbs-block-size") + 1], "0")
        self.assertEqual(cmd_l8[cmd_l8.index("--gibbs-block-size") + 1], "16")
        self.assertIn("--no-fixed-covariance", cmd_l8)
        self.assertIn("--recondition-covariance", cmd_l8)
        self.assertIn("--burnin-fraction", cmd_l8)
        self.assertIn("--low-acceptance-threshold", cmd_l8)

        with patch("moments2amplitudes.momentfit.runner.subprocess.run", return_value=MagicMock(returncode=0)) as mock_run:
            ok = run_collector(cfg=cfg, max_l_filter=[8], no_exec=True, dry_run=False)
        self.assertTrue(ok)
        mock_run.assert_called_once()
        col_cmd = mock_run.call_args[0][0]
        self.assertEqual(col_cmd[col_cmd.index("--gibbs-block-size") + 1], "16")
        self.assertEqual(col_cmd[col_cmd.index("--max-retry-iters") + 1], "200000")
        self.assertIn("--no-fixed-covariance", col_cmd)
        self.assertIn("--recondition-covariance", col_cmd)

    def test_aligned_acceptance_dynamic_burnin_and_reconditioning(self):
        import importlib.util
        import tempfile
        from unittest.mock import patch
        import numpy as np
        import uproot

        spec = importlib.util.spec_from_file_location("moment_fit_collector", "macros/MomentFitCollector.py")
        mfc = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mfc)

        # Full-vector acceptance (6%) stays inside [0.03, 0.35]: keep calibrated step, scale iters
        _, step_fv, iters_fv, _ = mfc.recommend_low_ess_params(
            bin_scale=1.0, bin_step=1.0, bin_iters=50000, burn_in=2000,
            eval_ess=10.0, ess_threshold=20.0, mean_tau=3000.0, bin_acceptance=0.06,
            bin_gibbs=0, low_acceptance_threshold=0.03, gibbs_low_acceptance_threshold=0.12,
        )
        self.assertEqual(step_fv, 1.0)
        self.assertGreater(iters_fv, 50000)

        # Genuinely stalled full-vector (2% < 3%): shrink step
        _, step_low, _, _ = mfc.recommend_low_ess_params(
            bin_scale=1.0, bin_step=1.0, bin_iters=50000, burn_in=2000,
            eval_ess=5.0, ess_threshold=20.0, mean_tau=5000.0, bin_acceptance=0.02,
            bin_gibbs=0, low_acceptance_threshold=0.03, gibbs_low_acceptance_threshold=0.12,
        )
        self.assertLess(step_low, 1.0)

        # Blocked Gibbs (8% < 12%): shrink step
        _, step_gibbs, _, _ = mfc.recommend_low_ess_params(
            bin_scale=1.0, bin_step=1.0, bin_iters=50000, burn_in=2000,
            eval_ess=8.0, ess_threshold=20.0, mean_tau=4000.0, bin_acceptance=0.08,
            bin_gibbs=12, low_acceptance_threshold=0.03, gibbs_low_acceptance_threshold=0.12,
        )
        self.assertLess(step_gibbs, 1.0)

        # Dynamic burn-in scaling
        self.assertEqual(mfc.compute_retry_burnin(2000, 200000, 0.10), 20000)
        self.assertEqual(mfc.compute_retry_burnin(2000, 10000, 0.10), 2000)
        self.assertEqual(mfc.compute_retry_burnin(2000, 200000, 0.0), 2000)

        # Shrinkage + SPD projection with correlated prior & rank-deficient empirical sample
        rng = np.random.default_rng(7)
        base = rng.normal(0.0, 1.0, size=(300, 1))
        samples = np.hstack([base, base * 2.0, rng.normal(1000.0, 250.0, size=(300, 1))])
        weights = rng.integers(1, 4, size=300)
        names = ["H_0_0_0", "H_0_2_0", "Yld_Moments"]
        prior = np.array([[1.0, 0.2, 0.0], [0.2, 4.0, 0.0], [0.0, 0.0, 62500.0]])
        spd = mfc._compute_shrinkage_spd_covariance(
            samples, weights, names, prior_cov=prior, shrinkage=0.25, min_eigen_floor=1e-3,
        )
        self.assertEqual(spd.shape, (3, 3))
        self.assertAlmostEqual(spd[0, 2], 0.0)
        self.assertAlmostEqual(spd[2, 0], 0.0)
        np.linalg.cholesky(spd)

        # End-to-end recondition_covariance_for_bin + retry script emission (mocked ROOT matrix I/O)
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir) / "run"
            tbin_dir = run_dir / "t010020" / "Mpi0eta0.820000_"
            src_cov_dir = Path(tmpdir) / "cov" / "t010020" / "All.maxL_4"
            tbin_dir.mkdir(parents=True)
            src_cov_dir.mkdir(parents=True)
            root_path = tbin_dir / "ResultsBruMcmcFixedCovariance.root"
            with uproot.recreate(root_path) as uf:
                uf["MCMCTree"] = {
                    "weight": np.full(600, 2, dtype=np.float64),
                    "H_0_0_0": rng.normal(1.0, 0.4, size=600),
                    "H_0_2_0": rng.normal(0.5, 0.3, size=600),
                    "Yld_Moments": rng.normal(1000.0, 50.0, size=600),
                }
            (src_cov_dir / "BruCovMatrix_Mpi0eta_0.820.root").write_bytes(b"stage2")
            (src_cov_dir / "covariance_matrix_Mpi0eta_0.820.root").write_bytes(b"init2")

            bin_rec = {
                "tbin_idx": 0,
                "tbin_name": "t010020",
                "bin_idx": 0,
                "mass_center": 0.82,
                "bin_dir_path": str(tbin_dir),
                "root_file": str(root_path),
                "failure_type": "LOW_ESS",
                "rec_scale": 1.0,
                "rec_step": 1.0,
                "rec_iters": 150000,
                "rec_burnin": 15000,
                "rec_gibbs": 0,
            }

            saved_paths = []

            def fake_save(out_path, cov_array, param_names, matrix_name="MyCovMatrix"):
                Path(out_path).parent.mkdir(parents=True, exist_ok=True)
                Path(out_path).write_text("mock_spd")
                saved_paths.append((out_path, cov_array.copy(), list(param_names)))
                return True

            with patch.object(mfc, "_load_bru_cov_matrix_root", return_value=(prior, names)), \
                 patch.object(mfc, "_save_bru_cov_matrix_root", side_effect=fake_save):
                recond_dir = mfc.recondition_covariance_for_bin(
                    bin_record=bin_rec,
                    run_dir=str(run_dir),
                    src_cov_dir=str(Path(tmpdir) / "cov"),
                    max_l=4,
                    burn_in=200,
                    shrinkage=0.2,
                )

            self.assertIsNotNone(recond_dir)
            self.assertIn("t010020", recond_dir)
            self.assertEqual(len(saved_paths), 1)
            self.assertTrue((Path(recond_dir) / "covariance_matrix_Mpi0eta_0.820.root").is_file())
            self.assertEqual((src_cov_dir / "BruCovMatrix_Mpi0eta_0.820.root").read_bytes(), b"stage2")

            bin_rec["rec_cov_dir"] = recond_dir
            retry_path = mfc.generate_retry_script(
                run_dir=str(run_dir),
                unsuccessful_bins=[bin_rec],
                macro_path="/d/home/septian/Moments2Amplitudes/macros/MomentFit.C",
                max_l=4,
                burn_in=2000,
                cov_dir=str(Path(tmpdir) / "cov"),
                fixed_covariance=True,
                recondition_covariance=True,
            )
            script_txt = Path(retry_path).read_text()
            self.assertIn("MOM2AMP_FIXED_COV=\"true\"", script_txt)
            self.assertIn(recond_dir, script_txt)
            self.assertIn(", 150000, 15000,", script_txt)

    def test_momentfit_retry_cli_and_bounded_loop(self):
        import json
        import tempfile
        from unittest.mock import patch, MagicMock
        from moments2amplitudes.momentfit.runner import parse_args, run_momentfit, run_retry

        # Supported CLI forms for retry
        args_flag = parse_args(["configs/momentfit_test.yaml", "--retry", "--max-l", "8", "--max-retries", "3", "-j", "4"])
        self.assertTrue(args_flag.retry)
        self.assertEqual(args_flag.stage, ["retry"])
        self.assertEqual(args_flag.max_retries, 3)
        self.assertEqual(args_flag.max_parallel, 4)

        args_sub = parse_args(["retry", "configs/momentfit_test.yaml", "--max-l", "8", "--tbin", "0"])
        self.assertEqual(args_sub.config, "configs/momentfit_test.yaml")
        self.assertTrue(args_sub.retry)
        self.assertEqual(args_sub.stage, ["retry"])

        args_trail = parse_args(["configs/momentfit_test.yaml", "retry", "--max-l", "8"])
        self.assertEqual(args_trail.config, "configs/momentfit_test.yaml")
        self.assertEqual(args_trail.stage, ["retry"])

        args_rundir = parse_args(["retry", "--run-dir", "/tmp/some_run_Lmax8"])
        self.assertIsNone(args_rundir.config)
        self.assertEqual(args_rundir.run_dir, "/tmp/some_run_Lmax8")
        self.assertEqual(args_rundir.stage, ["retry"])

        self.assertEqual(normalize_stage_names(["3", "retry"]), ["retry"])
        self.assertEqual(normalize_stage_names(["refit"]), ["retry"])

        with patch("moments2amplitudes.momentfit.stage3_moment_fit.run_stage3") as mock_s3, \
             patch("moments2amplitudes.momentfit.runner.run_retry", return_value=True) as mock_retry:
            rc = run_momentfit(
                config_path="configs/momentfit_test.yaml",
                stages=["retry"],
                max_retries=2,
                max_parallel=4,
                dry_run=True,
            )
        self.assertEqual(rc, 0)
        mock_s3.assert_not_called()
        mock_retry.assert_called_once()
        self.assertEqual(mock_retry.call_args.kwargs["max_retries"], 2)
        self.assertEqual(mock_retry.call_args.kwargs["max_parallel"], 4)

        # Verify bounded retry loop stops once unsuccessful_count drops to 0
        with tempfile.TemporaryDirectory() as tmpdir:
            rdir = Path(tmpdir) / "MomentFit_m104160_FixedCov_Lmax8_InitValueScale1.0_step1.0"
            rdir.mkdir(parents=True)
            summary_path = rdir / "unsuccessful_bins_summary.json"
            retry_script = rdir / "RunMomentFit_retry.sh"
            retry_script.write_text("#!/bin/bash\nexit 0\n")

            calls = []

            def fake_subproc(cmd, **kwargs):
                calls.append((list(cmd), dict(kwargs.get("env", {}))))
                if "MomentFitCollector.py" in " ".join(cmd):
                    col_count = sum(1 for c, _ in calls if "MomentFitCollector.py" in " ".join(c))
                    summary_path.write_text(json.dumps({"unsuccessful_count": 3 if col_count == 1 else 0}))
                return MagicMock(returncode=0)

            cfg = load_momentfit_config("configs/momentfit_test.yaml")
            with patch("moments2amplitudes.momentfit.runner.subprocess.run", side_effect=fake_subproc):
                ok = run_retry(
                    cfg=cfg,
                    run_dir=str(rdir),
                    tbins_filter=["t010020"],
                    max_retries=3,
                    max_parallel=4,
                    dry_run=False,
                )
            self.assertTrue(ok)
            collector_calls = [c for c, _ in calls if "MomentFitCollector.py" in " ".join(c)]
            bash_calls = [(c, env) for c, env in calls if c[0] == "bash"]
            self.assertEqual(len(collector_calls), 2)
            self.assertEqual(len(bash_calls), 1)
            self.assertIn("--non-interactive", bash_calls[0][0])
            self.assertIn("--no-exec", bash_calls[0][0])
            self.assertEqual(bash_calls[0][1].get("MAX_PARALLEL"), "4")

        # Verify cli.py mom2amp momentfit parser supports retry options
        import sys
        from moments2amplitudes import cli as mom_cli
        with patch("moments2amplitudes.momentfit.runner.run_momentfit", return_value=0) as mock_run_mf, \
             patch.object(sys, "argv", ["mom2amp", "momentfit", "retry", "configs/momentfit_test.yaml", "--max-retries", "2", "-j", "3", "--dry-run"]):
            with self.assertRaises(SystemExit) as cm:
                mom_cli.main()
            self.assertEqual(cm.exception.code, 0)
            mock_run_mf.assert_called_once()
            self.assertEqual(mock_run_mf.call_args.kwargs["stages"], ["retry"])
            self.assertEqual(mock_run_mf.call_args.kwargs["max_retries"], 2)
            self.assertEqual(mock_run_mf.call_args.kwargs["max_parallel"], 3)


class TestGeneralScriptAndAcceptanceRouting(unittest.TestCase):
    """Verify legacy _m080200.sh cleanup, RunMomentFit.sh --config routing, and acceptance CLI options."""

    @staticmethod
    def _repo_root() -> Path:
        return Path(__file__).resolve().parents[1]

    def _write_sandbox_safe_config(
        self,
        tmpdir: str,
        calculate_acceptance: str | None = None,
        use_legacy_alias: bool = False,
    ) -> tuple[Path, Path]:
        import yaml

        root = self._repo_root()
        tmp = Path(tmpdir)
        mm_cpp = tmp / "moments_main" / "cpp"
        mm_cpp.mkdir(parents=True, exist_ok=True)
        (mm_cpp / "basisFunctions.C").write_text("// stub\n", encoding="utf-8")

        data_dir = tmp / "input_data" / "t010020"
        data_dir.mkdir(parents=True, exist_ok=True)
        (data_dir / "data_flat_All.root").write_text("", encoding="utf-8")

        mfit_data_dir = tmp / "mfit_data"
        mfit_data_dir.mkdir(parents=True, exist_ok=True)

        run_dir = tmp / "run_out"
        run_dir.mkdir(parents=True, exist_ok=True)

        cfg_data = yaml.safe_load((root / "configs" / "momentfit_test.yaml").read_text(encoding="utf-8"))
        cfg_data["environment"]["moments_main_dir"] = str(tmp / "moments_main")
        cfg_data["analysis"]["tbins"] = ["t010020"]
        cfg_data["analysis"]["max_l"] = [4]
        cfg_data["input_data"]["base_dir"] = str(tmp / "input_data")
        cfg_data["linear_algebra"]["output_base_dir"] = str(tmp / "moments_out")
        cfg_data["moment_fit"]["data"]["base_dir"] = str(mfit_data_dir)
        cfg_data["moment_fit"]["output_dir_template"] = str(
            tmp / "run_m{tag_min:03d}{tag_max:03d}_L{maxL}_s{init_scale:.1f}_st{step:.1f}/"
        )
        if use_legacy_alias:
            cfg_data["moment_fit"].pop("calculate_acceptance", None)
            cfg_data["moment_fit"]["acceptance_macro"] = (
                calculate_acceptance or str(root / "macros" / "CalculateAcceptance_m080200.C")
            )
        elif calculate_acceptance is not None:
            cfg_data["moment_fit"]["calculate_acceptance"] = calculate_acceptance

        cfg_path = tmp / "momentfit_sandbox.yaml"
        cfg_path.write_text(yaml.safe_dump(cfg_data), encoding="utf-8")
        return cfg_path, run_dir

    def test_legacy_m080200_shell_script_removed(self):
        root = self._repo_root()
        self.assertFalse((root / "macros" / "legacy" / "RunMomentFit_m080200.sh").exists())

        for rel_path in [
            "macros/RunMomentFit.sh",
            "macros/MonitorMomentFit.py",
            "macros/MomentFit.C",
            "python/moments2amplitudes/momentfit/stage3_moment_fit.py",
            "python/moments2amplitudes/momentfit/runner.py",
            "python/moments2amplitudes/cli.py",
            "docs/MOMENTFIT.md",
        ]:
            content = (root / rel_path).read_text(encoding="utf-8")
            self.assertNotIn("RunMomentFit_m080200", content, f"Unexpected RunMomentFit_m080200 in {rel_path}")
            self.assertNotIn("MomentFit_m080200(", content, f"Unexpected MomentFit_m080200( in {rel_path}")
            self.assertNotIn("--acceptance-macro", content, f"Unexpected --acceptance-macro in {rel_path}")

    def test_acceptance_options_validate_stage3_and_monitor(self):
        import contextlib
        import io
        import tempfile
        from moments2amplitudes.momentfit.config import load_momentfit_config, validate_momentfit_config
        from moments2amplitudes.momentfit.stage3_moment_fit import build_moment_fit_command
        from moments2amplitudes.momentfit.runner import run_calculate_acceptance, run_monitor, run_momentfit

        with tempfile.TemporaryDirectory() as tmpdir:
            cfg_path, rdir = self._write_sandbox_safe_config(tmpdir)
            cfg = load_momentfit_config(str(cfg_path))
            self.assertTrue(cfg.moment_fit.calculate_acceptance.endswith("CalculateAcceptance_m080200.C"))
            self.assertEqual(cfg.moment_fit.acceptance_macro, cfg.moment_fit.calculate_acceptance)
            self.assertFalse(cfg.moment_fit.skip_acceptance)
            self.assertFalse(cfg.moment_fit.force_acceptance)
            self.assertEqual(validate_momentfit_config(cfg), [])

            # Legacy YAML key `acceptance_macro` still works as alias
            alias_cfg_path, _ = self._write_sandbox_safe_config(tmpdir, use_legacy_alias=True)
            alias_cfg = load_momentfit_config(str(alias_cfg_path))
            self.assertTrue(alias_cfg.moment_fit.calculate_acceptance.endswith("CalculateAcceptance_m080200.C"))

            # Stage 3 command forwarding
            cfg.moment_fit.skip_acceptance = True
            cfg.moment_fit.force_acceptance = True
            cmd, _ = build_moment_fit_command(cfg, max_l=4, tbin_arg="0")
            self.assertNotIn("--acceptance-macro", cmd)
            self.assertIn("--skip-acceptance", cmd)
            self.assertIn("--force-acceptance", cmd)

            # Invalid calculate_acceptance macro path caught by validator
            cfg.moment_fit.calculate_acceptance = str(Path(tmpdir) / "nonexistent_acc.C")
            errors = validate_momentfit_config(cfg)
            self.assertTrue(any("Acceptance macro not found" in err for err in errors), errors)

            # run_calculate_acceptance and run_momentfit forwarding (both boolean and explicit path)
            buf_acc = io.StringIO()
            with contextlib.redirect_stdout(buf_acc):
                ok_acc = run_calculate_acceptance(
                    cfg=load_momentfit_config(str(cfg_path)),
                    run_dir=str(rdir),
                    calculate_acceptance="/tmp/custom_acceptance.C",
                    skip_acceptance=True,
                    force_acceptance=True,
                    dry_run=True,
                )
            self.assertTrue(ok_acc)
            out_acc = buf_acc.getvalue()
            self.assertIn("--calculate-acceptance /tmp/custom_acceptance.C", out_acc)
            self.assertIn("--skip-acceptance", out_acc)
            self.assertIn("--force-acceptance", out_acc)
            self.assertIn(f"--data-dir {cfg.moment_fit.data.base_dir}", out_acc)

            buf_stage = io.StringIO()
            with contextlib.redirect_stdout(buf_stage):
                rc_acc = run_momentfit(
                    config_path=str(cfg_path),
                    stages=["calculate_acceptance"],
                    run_dir=str(rdir),
                    calculate_acceptance=True,
                    dry_run=True,
                )
            self.assertEqual(rc_acc, 0)
            self.assertIn("CalculateAcceptance_m080200.C", buf_stage.getvalue())

            buf_mon = io.StringIO()
            with contextlib.redirect_stdout(buf_mon):
                rc_mon = run_momentfit(
                    config_path=str(cfg_path),
                    stages=["monitor"],
                    run_dir=str(rdir),
                    min_ess_threshold=5.0,
                    skip_acceptance=True,
                    dry_run=True,
                )
            self.assertEqual(rc_mon, 0)
            self.assertIn("--skip-acceptance", buf_mon.getvalue())

    def test_run_momentfit_sh_config_routing_and_cli_flags(self):
        import os
        import subprocess
        import sys
        import tempfile
        from unittest.mock import patch
        from moments2amplitudes import cli as mom_cli

        root = self._repo_root()
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg_path, rdir = self._write_sandbox_safe_config(tmpdir)
            env = os.environ.copy()
            env["PYTHONPATH"] = f"{root / 'python'}:{env.get('PYTHONPATH', '')}".rstrip(":")

            for flag in ("--config", "-c"):
                res = subprocess.run(
                    [
                        "bash",
                        str(root / "macros" / "RunMomentFit.sh"),
                        flag,
                        str(cfg_path),
                        "--max-l",
                        "4",
                        "--tbin",
                        "0",
                        "--skip-acceptance",
                        "--force-acceptance",
                        "--dry-run",
                    ],
                    cwd=str(root),
                    env=env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(res.returncode, 0, f"{flag} failed:\n{res.stdout}\n{res.stderr}")
                self.assertIn("--skip-acceptance", res.stdout)
                self.assertIn("--force-acceptance", res.stdout)
                self.assertIn("[DRY-RUN]", res.stdout)
                self.assertNotIn("RuntimeWarning", res.stderr)

            # RunMomentFit.sh --config with --calculate-acceptance [path]
            res_acc = subprocess.run(
                [
                    "bash",
                    str(root / "macros" / "RunMomentFit.sh"),
                    "-c",
                    str(cfg_path),
                    "--run-dir",
                    str(rdir),
                    "--calculate-acceptance",
                    "/tmp/custom_acc.C",
                    "--dry-run",
                ],
                cwd=str(root),
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(res_acc.returncode, 0, f"calculate-acceptance failed:\n{res_acc.stdout}\n{res_acc.stderr}")
            self.assertIn("--calculate-acceptance /tmp/custom_acc.C", res_acc.stdout)

        # CLI tests: --calculate-acceptance with path, without path, and before config.yaml
        with patch("moments2amplitudes.momentfit.runner.run_momentfit", return_value=0) as mock_run_mf, \
             patch.object(
                 sys,
                 "argv",
                 [
                     "mom2amp",
                     "momentfit",
                     "configs/momentfit_test.yaml",
                     "--calculate-acceptance",
                     "/tmp/acc_override.C",
                     "--skip-acceptance",
                     "--force-acceptance",
                     "--dry-run",
                 ],
             ):
            with self.assertRaises(SystemExit) as cm:
                mom_cli.main()
            self.assertEqual(cm.exception.code, 0)
            self.assertEqual(mock_run_mf.call_args.kwargs["stages"], ["calculate_acceptance"])
            self.assertEqual(mock_run_mf.call_args.kwargs["calculate_acceptance"], "/tmp/acc_override.C")
            self.assertTrue(mock_run_mf.call_args.kwargs["skip_acceptance"])
            self.assertTrue(mock_run_mf.call_args.kwargs["force_acceptance"])

        with patch("moments2amplitudes.momentfit.runner.run_momentfit", return_value=0) as mock_run_mf, \
             patch.object(
                 sys,
                 "argv",
                 [
                     "mom2amp",
                     "momentfit",
                     "--calculate-acceptance",
                     "configs/momentfit_test.yaml",
                     "--dry-run",
                 ],
             ):
            with self.assertRaises(SystemExit) as cm:
                mom_cli.main()
            self.assertEqual(cm.exception.code, 0)
            self.assertEqual(mock_run_mf.call_args.kwargs["config_path"], "configs/momentfit_test.yaml")
            self.assertEqual(mock_run_mf.call_args.kwargs["stages"], ["calculate_acceptance"])
            self.assertIs(mock_run_mf.call_args.kwargs["calculate_acceptance"], True)


if __name__ == "__main__":
    unittest.main()

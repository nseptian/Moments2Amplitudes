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
        with self.assertRaises(ValueError):
            normalize_stage_names(["invalid_stage"])

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


if __name__ == "__main__":
    unittest.main()

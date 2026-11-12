"""
Tests for publication plotter engine and moment comparison logic.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock
from moments2amplitudes.schema import FitConfig, PlottingConfig, Chi2Config, WaveConfig, WaveType
from moments2amplitudes.rundir import RunDirectory
from moments2amplitudes.plotter import (
    discover_lowest_chi2_fits,
    smooth_curve,
    is_moment_fitted,
    extract_amplitude_l,
    is_amplitude_mass_dependent,
    compute_bin_edges,
    plot_tbins_parallel,
)


class TestPlotter(unittest.TestCase):
    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())
        self.rundir = RunDirectory(self.test_dir).ensure_directories()

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_discover_lowest_chi2_fits(self):
        # Create dummy log files
        for seed, chi2, status, valid in [
            (1, 100.5, 0, 1),
            (2, 85.2, 0, 1),
            (3, 120.0, 1, 0),
            (4, 75.1, 0, 1),
            (5, 95.0, 0, 1),
        ]:
            log_f = self.rundir.get_seed_log_file(seed)
            log_f.write_text(f"[RESULT] SEED={seed} CHI2={chi2} VALID={valid} STATUS={status}\n")
            # Create dummy root file
            root_f = self.rundir.fit_starts_dir / f"MassDependentResults_{seed}.root"
            root_f.touch()

        fits = discover_lowest_chi2_fits(self.rundir, n_lowest=3)
        self.assertEqual(len(fits), 3)
        self.assertAlmostEqual(fits[0][0], 75.1)
        self.assertAlmostEqual(fits[1][0], 85.2)
        self.assertAlmostEqual(fits[2][0], 95.0)

    def test_smooth_curve(self):
        # 1. Four or more points should produce smooth curve of length num_points
        xs = [1.06, 1.10, 1.14, 1.18, 1.22, 1.26]
        ys = [10.0, 20.0, 50.0, 45.0, 25.0, 12.0]
        xd, yd = smooth_curve(xs, ys, num_points=150)
        self.assertEqual(len(xd), 150)
        self.assertEqual(len(yd), 150)
        self.assertAlmostEqual(xd[0], 1.06)
        self.assertAlmostEqual(xd[-1], 1.26)

        # 2. Fewer than four points should fall back to original array
        xs_short = [1.0, 2.0]
        ys_short = [3.0, 4.0]
        xd_s, yd_s = smooth_curve(xs_short, ys_short)
        self.assertEqual(len(xd_s), 2)
        self.assertEqual(len(yd_s), 2)

    def test_is_moment_fitted(self):
        cfg = FitConfig()
        cfg.chi2 = Chi2Config(
            included_l=[4],
            excluded_l=[],
            max_l=4,
            no_use=["H_3"],
        )

        # L=4 moments should be fitted
        self.assertTrue(is_moment_fitted("H_0_4_0", cfg))
        self.assertTrue(is_moment_fitted("H_1_4_2", cfg))
        self.assertTrue(is_moment_fitted("H_2_4_4", cfg))

        # L != 4 should not be fitted
        self.assertFalse(is_moment_fitted("H_0_0_0", cfg))
        self.assertFalse(is_moment_fitted("H_0_2_0", cfg))
        self.assertFalse(is_moment_fitted("H_0_3_0", cfg))

        # Test with Stage 2 config (L in [0, 2, 4])
        cfg2 = FitConfig()
        cfg2.chi2 = Chi2Config(
            included_l=[0, 2, 4],
            excluded_l=[],
            max_l=4,
            no_use=["H_3"],
        )
        self.assertTrue(is_moment_fitted("H_0_0_0", cfg2))
        self.assertTrue(is_moment_fitted("H_0_2_0", cfg2))
        self.assertTrue(is_moment_fitted("H_0_4_0", cfg2))
        self.assertFalse(is_moment_fitted("H_0_1_0", cfg2))
        self.assertFalse(is_moment_fitted("H_0_3_0", cfg2))

        # Test with excluded_l
        cfg3 = FitConfig()
        cfg3.chi2 = Chi2Config(
            included_l=[],
            excluded_l=[2],
            max_l=4,
            no_use=[],
        )
        self.assertTrue(is_moment_fitted("H_0_0_0", cfg3))
        self.assertFalse(is_moment_fitted("H_0_2_0", cfg3))
        self.assertTrue(is_moment_fitted("H_0_4_0", cfg3))

    def test_amplitude_wave_type_classification(self):
        # Setup config with Piecewise S-wave (l=0) and Mass-Dependent D-wave (l=2)
        cfg = FitConfig()
        cfg.waves = [
            WaveConfig(l=0, name="S-wave", type=WaveType.PIECEWISE),
            WaveConfig(l=2, name="D-wave", type=WaveType.MASS_DEPENDENT),
        ]

        # Check L extraction
        self.assertEqual(extract_amplitude_l("a_0_0"), 0)
        self.assertEqual(extract_amplitude_l("b_0_0"), 0)
        self.assertEqual(extract_amplitude_l("aphi_0_0"), 0)
        self.assertEqual(extract_amplitude_l("a_2_-1"), 2)
        self.assertEqual(extract_amplitude_l("bphi_2_2"), 2)
        self.assertIsNone(extract_amplitude_l("mass_bin"))
        self.assertIsNone(extract_amplitude_l("H_0_0_0"))

        # Check mass-dependent vs mass-independent classification
        # S-wave (l=0) is piecewise -> mass-independent (False)
        self.assertFalse(is_amplitude_mass_dependent("a_0_0", cfg))
        self.assertFalse(is_amplitude_mass_dependent("b_0_0", cfg))
        self.assertFalse(is_amplitude_mass_dependent("aphi_0_0", cfg))

        # D-wave (l=2) is mass_dependent -> True
        self.assertTrue(is_amplitude_mass_dependent("a_2_0", cfg))
        self.assertTrue(is_amplitude_mass_dependent("a_2_1", cfg))
        self.assertTrue(is_amplitude_mass_dependent("b_2_-1", cfg))
        self.assertTrue(is_amplitude_mass_dependent("aphi_2_2", cfg))

    def test_compute_bin_edges(self):
        masses = [1.06, 1.10, 1.14]
        edges = compute_bin_edges(masses, default_step=0.04)
        self.assertEqual(len(edges), 4)
        self.assertAlmostEqual(edges[0], 1.04)
        self.assertAlmostEqual(edges[1], 1.08)
        self.assertAlmostEqual(edges[2], 1.12)
        self.assertAlmostEqual(edges[3], 1.16)

    def test_format_amplitude_squared_label(self):
        from moments2amplitudes.plotter import format_amplitude_squared_label
        self.assertEqual(format_amplitude_squared_label("a_2_1"), "$|a_{2,1}|^2$")
        self.assertEqual(format_amplitude_squared_label("b_2_-1"), "$|b_{2,-1}|^2$")
        self.assertEqual(format_amplitude_squared_label("a_0_0"), "$|a_{0,0}|^2$")
        self.assertEqual(format_amplitude_squared_label("other_var"), "$|other\\_var|^2$")

    @patch("moments2amplitudes.plotter.as_completed")
    @patch("moments2amplitudes.plotter.ProcessPoolExecutor")
    def test_plot_tbins_parallel(self, mock_executor_cls, mock_as_completed):
        mock_executor = MagicMock()
        mock_executor_cls.return_value.__enter__.return_value = mock_executor
        
        futures_list = []
        for t in ["t010020", "t020032"]:
            f = MagicMock()
            f.result.return_value = (t, True, 10, f"Generated 10 plots for {t}")
            futures_list.append(f)

        mock_executor.submit.side_effect = futures_list
        mock_as_completed.return_value = futures_list

        cfg = FitConfig()
        cfg.dataset.t_bins = ["t010020", "t020032"]

        res = plot_tbins_parallel(cfg, max_workers=2)
        self.assertTrue(res["success"])
        self.assertEqual(res["total_plots"], 20)
        self.assertEqual(len(res["tbins"]), 2)
        self.assertEqual(mock_executor.submit.call_count, 2)

    @patch("moments2amplitudes.plotter.plot_amplitudes_and_bands")
    def test_run_all_diagnostics_propagated_options(self, mock_plot_amps):
        from moments2amplitudes.plotter import run_all_diagnostics_and_plots
        cfg = FitConfig()
        mock_plot_amps.return_value = [self.rundir.plots_dir / "test.png"]

        # 1. only_distributions=True
        res = run_all_diagnostics_and_plots(cfg, self.rundir, only_distributions=True)
        self.assertEqual(len(res["plots"]), 1)
        self.assertIsNone(res["summary"])
        mock_plot_amps.assert_called_with(
            self.rundir.nominal_best_root,
            self.rundir.mcmc_propagated_root,
            self.rundir.plots_dir,
            cfg,
            formats=("png", "pdf"),
            only_distributions=True,
            only_amplitudes=False,
        )

        # 2. only_propagated=True
        res_prop = run_all_diagnostics_and_plots(cfg, self.rundir, only_propagated=True)
        self.assertEqual(len(res_prop["plots"]), 1)
        self.assertIsNone(res_prop["summary"])


if __name__ == "__main__":
    unittest.main()


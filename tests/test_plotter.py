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
    plot_phase_with_ambiguity,
    plot_moments_across_tbins,
    plot_amplitudes_across_tbins,
    plot_resonance_parameters_vs_t,
    plot_spd_total_intensity_vs_t,
    load_or_export_amplitude_data,
    plot_multistart_landscape_all_tbins,
    _get_hist_bins,
    _get_chi2_landscape_bins,
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

        # Place nominal_best_root
        self.rundir.nominal_best_root.touch()

        fits = discover_lowest_chi2_fits(self.rundir, n_lowest=3)
        self.assertEqual(len(fits), 3)
        # Seed 4 was 75.1 -> replaced by nominal_best_root without duplicating seed 4
        self.assertAlmostEqual(fits[0][0], 75.1)
        self.assertEqual(fits[0][1].resolve(), self.rundir.nominal_best_root.resolve())
        # Next is seed 2
        self.assertAlmostEqual(fits[1][0], 85.2)
        # Next is seed 5
        self.assertAlmostEqual(fits[2][0], 95.0)

    def test_discover_lowest_chi2_fits_all_failed(self):
        # When all multi-starts fail, discover_lowest_chi2_fits returns empty
        for seed in [1, 2, 3]:
            log_f = self.rundir.get_seed_log_file(seed)
            log_f.write_text(f"[RESULT] SEED={seed} CHI2=999.0 VALID=0 STATUS=3\n")
            root_f = self.rundir.fit_starts_dir / f"MassDependentResults_{seed}.root"
            root_f.touch()

        self.rundir.nominal_best_root.touch()
        fits = discover_lowest_chi2_fits(self.rundir, n_lowest=5)
        self.assertEqual(len(fits), 0)

    def test_discover_lowest_chi2_fits_rejects_at_limit(self):
        for seed, chi2, status, valid, at_lim in [
            (1, 50.0, 0, 1, 1),   # Lower chi2, but hit parameter limit -> rejected
            (2, 80.0, 0, 0, 1),   # Rejected by both valid=0 and AT_LIMIT=1
            (3, 95.0, 0, 1, 0),   # Valid interior fit
        ]:
            log_f = self.rundir.get_seed_log_file(seed)
            log_f.write_text(f"[RESULT] SEED={seed} CHI2={chi2} VALID={valid} AT_LIMIT={at_lim} STATUS={status}\n")
            (self.rundir.fit_starts_dir / f"MassDependentResults_{seed}.root").touch()

        fits = discover_lowest_chi2_fits(self.rundir, n_lowest=5)
        self.assertEqual(len(fits), 1)
        self.assertAlmostEqual(fits[0][0], 95.0)

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
            included_L=[4],
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
            included_L=[0, 2, 4],
            max_l=4,
            no_use=["H_3"],
        )
        self.assertTrue(is_moment_fitted("H_0_0_0", cfg2))
        self.assertTrue(is_moment_fitted("H_0_2_0", cfg2))
        self.assertTrue(is_moment_fitted("H_0_4_0", cfg2))
        self.assertFalse(is_moment_fitted("H_0_1_0", cfg2))
        self.assertFalse(is_moment_fitted("H_0_3_0", cfg2))

        # Test with Stage 3 SPD config (L in [0, 1, 2, 3, 4])
        cfg3 = FitConfig()
        cfg3.chi2 = Chi2Config(
            included_L=[0, 1, 2, 3, 4],
            max_l=4,
        )
        self.assertTrue(is_moment_fitted("H_0_0_0", cfg3))
        self.assertTrue(is_moment_fitted("H_0_1_0", cfg3))
        self.assertTrue(is_moment_fitted("H_0_2_0", cfg3))
        self.assertTrue(is_moment_fitted("H_0_3_0", cfg3))
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
        from moments2amplitudes.plotter import (
            format_amplitude_label,
            format_amplitude_symbol,
            format_amplitude_squared_label,
            format_amplitude_squared_symbol,
        )
        # Magnitude labels: [l]^refl_m (l=0,1,2,3 -> S,P,D,F; refl=+ for a, - for b)
        self.assertEqual(format_amplitude_label("a_0_0"), "$S^{+}_{0}$")
        self.assertEqual(format_amplitude_label("b_0_0"), "$S^{-}_{0}$")
        self.assertEqual(format_amplitude_label("a_1_0"), "$P^{+}_{0}$")
        self.assertEqual(format_amplitude_label("b_1_-1"), "$P^{-}_{-1}$")
        self.assertEqual(format_amplitude_label("a_2_1"), "$D^{+}_{1}$")
        self.assertEqual(format_amplitude_label("b_2_-1"), "$D^{-}_{-1}$")
        self.assertEqual(format_amplitude_label("a_3_2"), "$F^{+}_{2}$")
        self.assertEqual(format_amplitude_label("aphi_2_1"), "$\\phi(D^{+}_{1})$")
        self.assertEqual(format_amplitude_label("bphi_2_-1"), "$\\phi(D^{-}_{-1})$")

        # Squared labels: |[l]^refl_m|^2
        self.assertEqual(format_amplitude_squared_label("a_0_0"), "$|S^{+}_{0}|^2$")
        self.assertEqual(format_amplitude_squared_label("b_0_0"), "$|S^{-}_{0}|^2$")
        self.assertEqual(format_amplitude_squared_label("a_2_1"), "$|D^{+}_{1}|^2$")
        self.assertEqual(format_amplitude_squared_label("b_2_-1"), "$|D^{-}_{-1}|^2$")
        self.assertEqual(format_amplitude_squared_label("other_var"), "$|other\\_var|^2$")

        # Bare symbols (without surrounding $)
        self.assertEqual(format_amplitude_symbol("a_2_1"), "D^{+}_{1}")
        self.assertEqual(format_amplitude_squared_symbol("a_2_1"), "|D^{+}_{1}|^2")

    def test_format_mcmc_parameter_label(self):
        import matplotlib.mathtext as mathtext
        from moments2amplitudes.plotter import format_mcmc_parameter_label

        parser = mathtext.MathTextParser("path")
        test_cases = [
            ("chi2", r"$\chi^2$", r"$\chi^2$"),
            ("MI_1.300000_a_2_0", r"$D^{+}_{0}$ (1.30 GeV)", r"$D^{+}_{0}$ (1.30)"),
            ("MD_shared_a2_1320_Mass1", r"$M_{1}$ (a2\_1320)", r"$M_{1}$"),
            ("MD_shared_a2_1320_Width1", r"$\Gamma_{1}$ (a2\_1320)", r"$\Gamma_{1}$"),
            ("MD_shared_a2_1700_Mass2", r"$M_{2}$ (a2\_1700)", r"$M_{2}$"),
            ("MD_shared_a2_1700_Width2", r"$\Gamma_{2}$ (a2\_1700)", r"$\Gamma_{2}$"),
            ("MD_shared_a2_1320_Mass", r"$M$ (a2\_1320)", r"$M$"),
            ("MD_shared_a2_1320_Width", r"$\Gamma$ (a2\_1320)", r"$\Gamma$"),
            ("MD_shareda_a2_1700_globalphase", r"$\phi^{(+)}(a2\_1700)$", r"$\phi^{(+)}(a2\_1700)$"),
            ("MD_a_2_-2_1700_k1", r"$k_1(D^{+}_{-2})$", r"$k_1(D^{+}_{-2})$"),
            ("MD_a_2_-2_a2_1320_1700_k1", r"$k_1(D^{+}_{-2})$", r"$k_1(D^{+}_{-2})$"),
            ("MD_a_2_-2_a2_1320_1700_k2", r"$k_2(D^{+}_{-2})$", r"$k_2(D^{+}_{-2})$"),
            ("MD_a_2_0_a2_1320_k", r"$k(D^{+}_{0})$", r"$k(D^{+}_{0})$"),
            ("MD_b_2_1_a2_1700_k", r"$k(D^{-}_{1})$", r"$k(D^{-}_{1})$"),
            ("MD_shared_f0_980_g_etapi", r"$g_{\eta\pi}$", r"$g_{\eta\pi}$"),
        ]

        for raw_name, exp_long, exp_short in test_cases:
            long_lbl = format_mcmc_parameter_label(raw_name, short=False)
            short_lbl = format_mcmc_parameter_label(raw_name, short=True)
            self.assertEqual(long_lbl, exp_long, f"Mismatch for long label: {raw_name}")
            self.assertEqual(short_lbl, exp_short, f"Mismatch for short label: {raw_name}")
            # Ensure both can be parsed by matplotlib mathtext without double subscript or syntax errors
            try:
                parser.parse(long_lbl)
                parser.parse(short_lbl)
            except Exception as e:
                self.fail(f"MathTextParser failed on {raw_name} -> long: {long_lbl}, short: {short_lbl}: {e}")

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
            style="default",
            rundir=self.rundir,
            unnormalized_moments_file=None,
        )

        # 2. only_propagated=True
        res_prop = run_all_diagnostics_and_plots(cfg, self.rundir, only_propagated=True)
        self.assertEqual(len(res_prop["plots"]), 1)
        self.assertIsNone(res_prop["summary"])

    def test_apply_plot_style_options(self):
        from moments2amplitudes.plotter import apply_plot_style
        # Should cleanly apply without raising exceptions
        for st in ["LHCb", "LHCb2", "ROOT", "ATLAS", "CMS", "ALICE", "default"]:
            apply_plot_style(st)

    def test_resolve_unnormalized_moments_and_h0_overlay(self):
        from moments2amplitudes.plotter import resolve_unnormalized_moments_file, load_h0_moment_data
        import json

        cfg = FitConfig()
        # Direct override
        dummy_file = self.test_dir / "UnnormalizedMoments_vs_Mass.root"
        dummy_file.touch()
        self.assertEqual(resolve_unnormalized_moments_file(cfg, self.rundir, str(dummy_file)), dummy_file)

        # Config plotting path
        cfg.plotting.unnormalized_moments_file = str(dummy_file)
        self.assertEqual(resolve_unnormalized_moments_file(cfg, self.rundir), dummy_file)

        # Load H0 moment data from cached JSON
        cache_json = self.rundir.root / ".cache" / "moment_comparison_data.json"
        cache_json.parent.mkdir(parents=True, exist_ok=True)
        sample_json = {
            "h0_data": {
                "mass": [1.06, 1.10, 1.14],
                "val": [1000.0, 2000.0, 3000.0],
                "err": [50.0, 60.0, 70.0],
            }
        }
        with open(cache_json, "w") as f:
            json.dump(sample_json, f)

        loaded_h0 = load_h0_moment_data(dummy_file, cache_json)
        self.assertIsNotNone(loaded_h0)
        self.assertEqual(loaded_h0["mass"], [1.06, 1.10, 1.14])
        self.assertEqual(loaded_h0["val"], [1000.0, 2000.0, 3000.0])

    def test_plot_amplitudes_and_bands_2panel(self):
        from moments2amplitudes.plotter import plot_amplitudes_and_bands
        import json

        cfg = FitConfig()
        out_dir = self.test_dir / "plots_2panel"
        out_dir.mkdir(parents=True, exist_ok=True)
        cache_dir = out_dir / ".cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        json_path = cache_dir / "amplitudes_data.json"

        data = {
            "masses": [1.1, 1.2, 1.3, 1.4],
            "nominal": {
                "a_2_1": [10.0, 20.0, 30.0, 40.0],
                "aphi_2_1": [0.1, 0.2, 0.3, 0.4],
                "a_0_0": [5.0, 8.0, 12.0, 15.0],
                "aphi_0_0": [0.5, 0.6, 0.7, 0.8],
            },
            "mcmc": {
                "a_2_1": {
                    "1.1": [9.0, 10.0, 11.0],
                    "1.2": [19.0, 20.0, 21.0],
                    "1.3": [29.0, 30.0, 31.0],
                    "1.4": [39.0, 40.0, 41.0],
                },
                "aphi_2_1": {
                    "1.1": [0.08, 0.1, 0.12],
                    "1.2": [0.18, 0.2, 0.22],
                    "1.3": [0.28, 0.3, 0.32],
                    "1.4": [0.38, 0.4, 0.42],
                },
            },
            "h0_data": {
                "mass": [1.1, 1.2, 1.3, 1.4],
                "val": [100.0, 400.0, 900.0, 1600.0],
                "err": [10.0, 20.0, 30.0, 40.0],
            },
        }
        with open(json_path, "w") as f:
            json.dump(data, f)

        dummy_nom = self.test_dir / "dummy_nom.root"
        dummy_nom.touch()

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            plots = plot_amplitudes_and_bands(dummy_nom, None, out_dir, cfg, formats=("png",))

        plot_names = {p.name for p in plots}
        self.assertIn("amplitudes_magnitude_vs_mass.png", plot_names)
        self.assertIn("amplitudes_intensity_vs_mass.png", plot_names)
        self.assertIn("amplitudes_phase_vs_mass.png", plot_names)
        # Reference amplitude phase is excluded from phase plots
        self.assertNotIn("a_2_1_intensity_and_phase.png", plot_names)
        self.assertIn("a_0_0_intensity_and_phase.png", plot_names)
        self.assertIn("overview_D_waves_positive.png", plot_names)
        self.assertIn("overview_S_waves.png", plot_names)
        self.assertIn("overview_all_waves.png", plot_names)
        self.assertIn("amplitudes_overview_vs_mass.png", plot_names)

    def test_plot_phase_with_ambiguity(self):
        import numpy as np
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        masses = [1.1, 1.2, 1.3, 1.4, 1.5]
        y_nom = [0.2, 0.5, 1.1, 1.6, 2.0]
        mcmc_dict = {round(m, 4): [v - 0.1, v + 0.1] for m, v in zip(masses, y_nom)}

        # 1. Test discrete 68% band (density_band=False)
        plot_phase_with_ambiguity(
            ax=ax,
            masses=masses,
            y_nom=y_nom,
            is_mass_dep=True,
            color="black",
            label="$\\pm\\phi$",
            mcmc_dict=mcmc_dict,
            show_ambiguity=True,
            density_band=False,
        )

        # Verify lines: primary solid and sign-ambiguity dashed
        lines = ax.get_lines()
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[0].get_linestyle(), "-")
        self.assertEqual(lines[1].get_linestyle(), "--")
        # Check y values: second line is negated
        np.testing.assert_allclose(lines[1].get_ydata(), -lines[0].get_ydata())

        # Verify credible interval fills: two PolyCollections (positive and negative)
        collections = ax.collections
        self.assertEqual(len(collections), 2)
        plt.close(fig)

        # 2. Test continuous posterior density band (density_band=True) via individual transparent traces
        fig2, ax2 = plt.subplots()
        plot_phase_with_ambiguity(
            ax=ax2,
            masses=masses,
            y_nom=y_nom,
            is_mass_dep=True,
            color="black",
            label="$\\pm\\phi$",
            mcmc_dict=mcmc_dict,
            show_ambiguity=True,
            density_band=True,
        )
        from matplotlib.collections import LineCollection
        self.assertEqual(len(ax2.collections), 2)
        self.assertTrue(all(isinstance(c, LineCollection) for c in ax2.collections))
        plt.close(fig2)

    def test_physics_linestyles_and_monochrome_colors(self):
        from moments2amplitudes.plotter import PHYSICS_LINESTYLES, plot_amplitudes_and_bands
        import json
        import matplotlib.pyplot as plt

        self.assertGreaterEqual(len(PHYSICS_LINESTYLES), 4)

        cfg = FitConfig()
        out_dir = self.test_dir / "plots_mono"
        out_dir.mkdir(parents=True, exist_ok=True)
        cache_dir = out_dir / ".cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        json_path = cache_dir / "amplitudes_data.json"

        data = {
            "masses": [1.1, 1.2, 1.3, 1.4],
            "nominal": {
                "a_2_1": [10.0, 20.0, 30.0, 40.0],
                "aphi_2_1": [0.1, 0.2, 0.3, 0.4],
            },
            "mcmc": {
                "a_2_1": {
                    "1.1": [9.0, 10.0, 11.0],
                    "1.2": [19.0, 20.0, 21.0],
                    "1.3": [29.0, 30.0, 31.0],
                    "1.4": [39.0, 40.0, 41.0],
                },
            },
            "h0_data": {
                "mass": [1.1, 1.2, 1.3, 1.4],
                "val": [100.0, 400.0, 900.0, 1600.0],
                "err": [10.0, 20.0, 30.0, 40.0],
            },
        }
        with open(json_path, "w") as f:
            json.dump(data, f)

        dummy_nom = self.test_dir / "dummy_nom_mono.root"
        dummy_nom.touch()

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            plots = plot_amplitudes_and_bands(dummy_nom, None, out_dir, cfg, formats=("png",))

        self.assertTrue(len(plots) > 0)

    def test_configure_phase_axis_and_limits(self):
        import math
        import matplotlib.pyplot as plt
        from moments2amplitudes.plotter import configure_phase_axis, PHASE_Y_LIMITS, PHASE_Y_TICKS, PHASE_Y_TICKLABELS

        fig, ax = plt.subplots()
        configure_phase_axis(ax)

        y_min, y_max = ax.get_ylim()
        self.assertAlmostEqual(y_min, -math.pi)
        self.assertAlmostEqual(y_max, math.pi)

        ticks = ax.get_yticks()
        self.assertEqual(len(ticks), 5)
        for t_act, t_exp in zip(ticks, PHASE_Y_TICKS):
            self.assertAlmostEqual(t_act, t_exp)

        labels = [t.get_text() for t in ax.get_yticklabels()]
        self.assertEqual(labels, PHASE_Y_TICKLABELS)
        plt.close(fig)

    def test_plot_moments_across_tbins(self):
        import json
        out_dir = self.test_dir / "plots" / "moments_vs_t"
        cfg = FitConfig()
        cfg.project.name = "TestProject"
        cfg.dataset.tbins = ["t010020", "t020032"]

        # Single t-bin should return empty
        self.assertEqual(plot_moments_across_tbins(cfg, tbins=["t010020"]), [])

        # Create dummy moment comparison data in rundir caches
        for t in ["t010020", "t020032"]:
            t_dir = self.test_dir / t
            cache_dir = t_dir / ".cache"
            cache_dir.mkdir(parents=True, exist_ok=True)
            json_file = cache_dir / "moment_comparison_data.json"
            dummy_data = {
                "data": {
                    "H_0_0_0": {
                        "mass": [1.1, 1.2, 1.3],
                        "val": [100.0, 150.0, 120.0],
                        "err": [10.0, 15.0, 12.0],
                    },
                    "H_0_2_0": {
                        "mass": [1.1, 1.2, 1.3],
                        "val": [20.0, 30.0, 25.0],
                        "err": [5.0, 6.0, 5.5],
                    },
                },
                "fits": [
                    {
                        "chi2": 10.0,
                        "status": 0,
                        "is_valid": True,
                        "moments": {
                            "H_0_0_0": {"mass": [1.1, 1.2, 1.3], "val": [98.0, 148.0, 122.0]},
                            "H_0_2_0": {"mass": [1.1, 1.2, 1.3], "val": [19.0, 31.0, 24.0]},
                        },
                    }
                ],
            }
            with open(json_file, "w") as f:
                json.dump(dummy_data, f)

        plots = plot_moments_across_tbins(
            cfg,
            tbins=["t010020", "t020032"],
            run_dir=str(self.test_dir),
            output_dir=out_dir,
            formats=("png",),
        )

        self.assertTrue(len(plots) > 0)
        self.assertTrue((out_dir / "H_0_0_0.png").is_file())
        self.assertTrue((out_dir / "H_0_2_0.png").is_file())
        self.assertTrue((out_dir / "moments_overview_all_tbins.png").is_file())

    def test_plot_multistart_landscape_all_tbins(self):
        out_dir = self.test_dir / "ms_plots"
        dummy_results = {
            "t010020": {
                "results": [
                    {"chi2": 834.5, "status": 0, "is_valid": True},
                    {"chi2": 840.1, "status": 0, "is_valid": True},
                    {"chi2": 850.3, "status": 0, "is_valid": True},
                ]
            },
            "t020032": {
                "results": [
                    {"chi2": 910.2, "status": 0, "is_valid": True},
                    {"chi2": 915.0, "status": 0, "is_valid": True},
                    {"chi2": 925.8, "status": 0, "is_valid": True},
                ]
            },
        }

        generated = plot_multistart_landscape_all_tbins(
            dummy_results,
            output_dir=out_dir,
            project_name="TestProject",
            plateau_delta=10.0,
            formats=("png",),
        )

        self.assertEqual(len(generated), 1)
        self.assertTrue((out_dir / "multistart_chi2_landscape_all_tbins.png").is_file())

    def test_finer_binning_helpers(self):
        # 1D distributions
        # Low N
        self.assertEqual(_get_hist_bins(3), 3)
        # Moderate N (typical MCMC posterior ~150-200)
        self.assertGreaterEqual(_get_hist_bins(150), 50)
        self.assertLessEqual(_get_hist_bins(150), 100)
        # Overview grid (more compact)
        self.assertGreaterEqual(_get_hist_bins(150, is_overview=True), 35)
        self.assertLessEqual(_get_hist_bins(150, is_overview=True), 60)
        # Explicit user bins
        self.assertEqual(_get_hist_bins(200, user_bins=75), 75)

        # Multi-start Chi2 landscape
        # Typical 100 starts
        self.assertGreaterEqual(_get_chi2_landscape_bins(100), 60)
        # Large 1000 starts
        self.assertGreaterEqual(_get_chi2_landscape_bins(1000), 100)
        # Explicit user bins
        self.assertEqual(_get_chi2_landscape_bins(500, user_bins=80), 80)

    def test_plot_multistart_landscape_custom_bins(self):
        out_dir = self.test_dir / "ms_custom_bins"
        dummy_results = {
            "t010020": [800.0 + i * 2.0 for i in range(50)],
        }
        generated = plot_multistart_landscape_all_tbins(
            dummy_results,
            output_dir=out_dir,
            project_name="TestProject",
            bins=45,
            formats=("png",),
        )
        self.assertEqual(len(generated), 1)
        self.assertTrue((out_dir / "multistart_chi2_landscape_all_tbins.png").is_file())

    def test_plot_amplitudes_across_tbins(self):
        import json
        out_dir = self.test_dir / "plots" / "amplitudes_vs_t"
        cfg = FitConfig()
        cfg.project.name = "TestProject"
        cfg.dataset.tbins = ["t010020", "t020032"]

        # Single t-bin should return empty
        self.assertEqual(plot_amplitudes_across_tbins(cfg, tbins=["t010020"]), [])

        # Create dummy amplitude data in rundir caches for two t-bins
        for t in ["t010020", "t020032"]:
            t_dir = self.test_dir / t
            cache_dir = t_dir / ".cache"
            cache_dir.mkdir(parents=True, exist_ok=True)
            p_dir = t_dir / "plots"
            p_dir.mkdir(parents=True, exist_ok=True)
            json_file = cache_dir / "amplitudes_data.json"
            dummy_data = {
                "masses": [1.1, 1.2, 1.3, 1.4],
                "nominal": {
                    "a_2_1": [10.0, 20.0, 30.0, 40.0],
                    "aphi_2_1": [0.1, 0.2, 0.3, 0.4],
                    "a_0_0": [5.0, 8.0, 12.0, 15.0],
                    "aphi_0_0": [0.5, 0.6, 0.7, 0.8],
                },
                "mcmc": {
                    "a_2_1": {
                        "1.1": [9.0, 10.0, 11.0],
                        "1.2": [19.0, 20.0, 21.0],
                        "1.3": [29.0, 30.0, 31.0],
                        "1.4": [39.0, 40.0, 41.0],
                    },
                    "aphi_2_1": {
                        "1.1": [0.08, 0.1, 0.12],
                        "1.2": [0.18, 0.2, 0.22],
                        "1.3": [0.28, 0.3, 0.32],
                        "1.4": [0.38, 0.4, 0.42],
                    },
                    "a_0_0": {
                        "1.1": [4.5, 5.0, 5.5],
                        "1.2": [7.5, 8.0, 8.5],
                        "1.3": [11.5, 12.0, 12.5],
                        "1.4": [14.5, 15.0, 15.5],
                    },
                },
                "h0_data": {
                    "mass": [1.1, 1.2, 1.3, 1.4],
                    "val": [100.0, 400.0, 900.0, 1600.0],
                    "err": [10.0, 20.0, 30.0, 40.0],
                },
            }
            with open(json_file, "w") as f:
                json.dump(dummy_data, f)

        plots = plot_amplitudes_across_tbins(
            cfg,
            tbins=["t010020", "t020032"],
            run_dir=str(self.test_dir),
            output_dir=out_dir,
            formats=("png",),
        )

        self.assertTrue(len(plots) > 0)
        # Check individual wave files
        self.assertTrue((out_dir / "individual" / "a_2_1_intensity_vs_t.png").is_file())
        self.assertTrue((out_dir / "individual" / "a_0_0_intensity_vs_t.png").is_file())
        # Reference amplitude phase is excluded from phase plots as required
        self.assertFalse((out_dir / "individual" / "aphi_2_1_phase_vs_t.png").exists())
        # Non-reference amplitude relative phase is plotted
        self.assertTrue((out_dir / "individual" / "aphi_0_0_phase_vs_t.png").is_file())
        # Check summed wave files
        self.assertTrue((out_dir / "total_intensity_vs_t.png").is_file())
        self.assertTrue((out_dir / "D_waves_total_intensity_vs_t.png").is_file())
        self.assertTrue((out_dir / "S_wave_total_intensity_vs_t.png").is_file())
        # Check overview files
        self.assertTrue((out_dir / "amplitudes_overview_all_tbins.png").is_file())
        # Verify copied to t-bin plot dirs
        for t in ["t010020", "t020032"]:
            self.assertTrue((self.test_dir / t / "plots" / "amplitudes_overview_all_tbins.png").is_file())

    def test_plot_amplitudes_across_tbins_sharey(self):
        import json
        out_dir = self.test_dir / "plots" / "amps_sharey"
        cfg = FitConfig()
        cfg.project.name = "TestProjectShareY"
        cfg.dataset.tbins = ["t010020", "t020032"]

        for t in ["t010020", "t020032"]:
            t_dir = self.test_dir / t
            cache_dir = t_dir / ".cache"
            cache_dir.mkdir(parents=True, exist_ok=True)
            p_dir = t_dir / "plots"
            p_dir.mkdir(parents=True, exist_ok=True)
            dummy_data = {
                "masses": [1.1, 1.2, 1.3],
                "nominal": {
                    "a_2_1": [10.0, 20.0, 30.0],
                },
            }
            with open(cache_dir / "amplitudes_data.json", "w") as f:
                json.dump(dummy_data, f)

        plots = plot_amplitudes_across_tbins(
            cfg,
            tbins=["t010020", "t020032"],
            run_dir=str(self.test_dir),
            output_dir=out_dir,
            formats=("png",),
            sharey=True,
        )
        self.assertTrue(len(plots) > 0)
        self.assertTrue((out_dir / "individual" / "a_2_1_intensity_vs_t.png").is_file())
        self.assertTrue((out_dir / "total_intensity_vs_t.png").is_file())

    def test_wrap_angle_pi(self):
        from moments2amplitudes.plotter import wrap_angle_pi
        import numpy as np

        self.assertAlmostEqual(wrap_angle_pi(0.0).item(), 0.0)
        self.assertAlmostEqual(wrap_angle_pi(np.pi / 2.0).item(), np.pi / 2.0)
        self.assertAlmostEqual(wrap_angle_pi(-np.pi / 2.0).item(), -np.pi / 2.0)
        self.assertAlmostEqual(wrap_angle_pi(3.0 * np.pi / 2.0).item(), -np.pi / 2.0)
        self.assertAlmostEqual(wrap_angle_pi(2.0 * np.pi).item(), 0.0)

    def test_split_phase_jumps(self):
        from moments2amplitudes.plotter import split_phase_jumps
        import numpy as np

        m = [1.0, 1.1, 1.2, 1.3]
        # Smooth phase -> no NaN
        y_smooth = [0.1, 0.2, 0.3, 0.4]
        m_out, y_out = split_phase_jumps(m, y_smooth)
        self.assertEqual(len(m_out), 4)
        self.assertFalse(np.any(np.isnan(y_out)))

        # Jump across branch cut (+3.1 to -3.1, jump ~ 6.2 > 3.0)
        y_jump = [2.9, 3.1, -3.1, -2.9]
        m_out, y_out = split_phase_jumps(m, y_jump)
        self.assertEqual(len(m_out), 5)
        self.assertTrue(np.isnan(y_out[2]))

    def test_format_relative_phase_labels(self):
        from moments2amplitudes.plotter import (
            format_relative_phase_symbol,
            format_relative_phase_label,
        )

        sym = format_relative_phase_symbol("aphi_1_1", "aphi_2_1")
        self.assertIn("P^{+}_{1}", sym)
        self.assertIn("D^{+}_{1}", sym)
        self.assertIn("\\Delta\\phi", sym)

        lbl = format_relative_phase_label("aphi_1_1", "aphi_2_1", plus_minus=True)
        self.assertTrue(lbl.startswith("$\\pm\\Delta\\phi"))
        self.assertIn("P^{+}_{1}", lbl)
        self.assertIn("D^{+}_{1}", lbl)

    def test_resolve_anchor_amplitudes(self):
        from moments2amplitudes.plotter import resolve_anchor_amplitudes
        from moments2amplitudes.schema import FitConfig, PlottingConfig, WaveConfig, ResonanceConfig

        # 1. Explicit config
        cfg_explicit = FitConfig(
            project=None,
            dataset=None,
            waves=[],
            plotting=PlottingConfig(reference_amplitude="a_2_1"),
        )
        nominal_amps = {
            "a_2_1": [1.0, 2.0],
            "aphi_2_1": [0.5, 0.5],
            "a_2_2": [10.0, 20.0],
            "aphi_2_2": [1.0, 1.0],
            "b_2_1": [1.0, 2.0],
            "bphi_2_1": [0.2, 0.2],
        }
        active_phases = ["aphi_2_1", "aphi_2_2", "bphi_2_1"]
        anchors = resolve_anchor_amplitudes(cfg_explicit, nominal_amps, active_phases)
        self.assertEqual(anchors["a"], "aphi_2_1")
        self.assertEqual(anchors["b"], "bphi_2_1")

        # 2. Auto-detection: D-wave has phase_fixed=True, a_2_2 is dominant
        cfg_auto = FitConfig(
            project=None,
            dataset=None,
            waves=[
                WaveConfig(
                    l=2,
                    name="D-wave",
                    resonances=[
                        ResonanceConfig(name="a2_1320", mass=1.32, width=0.11, phase_fixed=True)
                    ],
                )
            ],
            plotting=PlottingConfig(),
        )
        anchors_auto = resolve_anchor_amplitudes(cfg_auto, nominal_amps, active_phases)
        self.assertEqual(anchors_auto["a"], "aphi_2_2")

    def test_plot_resonance_parameters_vs_t(self):
        from moments2amplitudes.schema import ModelType, ResonanceConfig, MassBinsConfig, DatasetConfig

        cfg = FitConfig(
            project=None,
            dataset=DatasetConfig(moments_dir="", mass_bins=MassBinsConfig(min=1.0, max=1.6, step=0.2)),
            waves=[
                WaveConfig(l=0, name="S-wave", type=WaveType.PIECEWISE),
                WaveConfig(
                    l=1,
                    name="P-wave",
                    type=WaveType.MASS_DEPENDENT,
                    model=ModelType.BARRIER_BREIT_WIGNER,
                    poly_order=-1,
                    resonances=[
                        ResonanceConfig(
                            name="pi1_1400",
                            mass=1.354,
                            mass_fixed=False,
                            mass_limits=(1.1, 1.8),
                            width=0.33,
                            width_fixed=True,
                        )
                    ],
                ),
                WaveConfig(
                    l=2,
                    name="D-wave",
                    type=WaveType.MASS_DEPENDENT,
                    model=ModelType.BARRIER_TWO_BREIT_WIGNER,
                    poly_order=-1,
                    resonances=[
                        ResonanceConfig(
                            name="a2_1320",
                            mass=1.318,
                            mass_fixed=False,
                            mass_limits=(1.25, 1.40),
                            width=0.107,
                            width_fixed=False,
                            width_limits=(0.08, 0.15),
                            phase_fixed=True,
                        ),
                        ResonanceConfig(
                            name="a2_1700",
                            mass=1.706,
                            mass_fixed=True,
                            width=0.38,
                            width_fixed=False,
                            width_limits=(0.20, 0.40),
                            phase=0.5,
                            phase_fixed=False,
                        ),
                    ],
                ),
            ],
            plotting=PlottingConfig(formats=["png"]),
        )

        t_data = {
            "t010020": {
                "name": "t010020",
                "t_center": 0.15,
                "nominal_params": {
                    "MD_shared_pi1_1400_Mass": 1.38,
                    "MD_shared_pi1_1400_Width": 0.33,
                    "MD_shared_a2_1320_1700_Mass1": 1.321,
                    "MD_shared_a2_1320_1700_Width1": 0.109,
                    "MD_shared_a2_1320_1700_Mass2": 1.706,
                    "MD_shared_a2_1320_1700_Width2": 0.3999,
                    "MD_shared_a2_1700_GlobalPhase": 0.45,
                    "MD_a_2_1_a2_1320_1700_k": 120.0,
                    "MD_a_2_1_a2_1320_1700_a0": 0.0,
                },
                "nominal_param_errors": {
                    "MD_shared_pi1_1400_Mass": 0.03,
                    "MD_shared_a2_1320_1700_Mass1": 0.004,
                    "MD_shared_a2_1320_1700_Width1": 0.006,
                    "MD_shared_a2_1320_1700_Width2": 0.02,
                    "MD_shared_a2_1700_GlobalPhase": 0.12,
                    "MD_a_2_1_a2_1320_1700_k": 8.0,
                },
                "mcmc_params": {
                    "MD_shared_pi1_1400_Mass": [1.36, 1.38, 1.41],
                    "MD_shared_a2_1320_1700_Mass1": [1.318, 1.321, 1.325],
                },
            },
            "t020032": {
                "name": "t020032",
                "t_center": 0.26,
                "nominal_params": {
                    "MD_shared_pi1_1400_Mass": 1.42,
                    "MD_shared_pi1_1400_Width": 0.33,
                    "MD_shared_a2_1320_1700_Mass1": 1.317,
                    "MD_shared_a2_1320_1700_Width1": 0.105,
                    "MD_shared_a2_1320_1700_Mass2": 1.706,
                    "MD_shared_a2_1320_1700_Width2": 0.34,
                    "MD_shared_a2_1700_GlobalPhase": 0.62,
                    "MD_a_2_1_a2_1320_1700_k": 95.0,
                    "MD_a_2_1_a2_1320_1700_a0": 0.0,
                },
                "nominal_param_errors": {
                    "MD_shared_pi1_1400_Mass": 0.04,
                    "MD_shared_a2_1320_1700_Mass1": 0.005,
                    "MD_shared_a2_1320_1700_Width1": 0.007,
                    "MD_shared_a2_1320_1700_Width2": 0.03,
                    "MD_shared_a2_1700_GlobalPhase": 0.15,
                    "MD_a_2_1_a2_1320_1700_k": 6.0,
                },
                "mcmc_params": {},
            },
            "t032050": {
                "name": "t032050",
                "t_center": 0.41,
                "nominal_params": {
                    "MD_shared_pi1_1400_Mass": 1.39,
                    "MD_shared_pi1_1400_Width": 0.33,
                    "MD_shared_a2_1320_1700_Mass1": 1.320,
                    "MD_shared_a2_1320_1700_Width1": 0.111,
                    "MD_shared_a2_1320_1700_Mass2": 1.706,
                    "MD_shared_a2_1320_1700_Width2": 0.36,
                    "MD_shared_a2_1700_GlobalPhase": 0.51,
                    "MD_a_2_1_a2_1320_1700_k": 72.0,
                    "MD_a_2_1_a2_1320_1700_a0": 0.0,
                },
                "nominal_param_errors": {
                    "MD_shared_pi1_1400_Mass": 0.05,
                    "MD_shared_a2_1320_1700_Mass1": 0.006,
                    "MD_shared_a2_1320_1700_Width1": 0.008,
                    "MD_shared_a2_1320_1700_Width2": 0.04,
                    "MD_shared_a2_1700_GlobalPhase": 0.18,
                    "MD_a_2_1_a2_1320_1700_k": 5.0,
                },
                "mcmc_params": {},
            },
        }

        out_dir = self.test_dir / "res_params_vs_t"
        generated = plot_resonance_parameters_vs_t(
            cfg, t_data, ["t010020", "t020032", "t032050"], out_dir, formats=("png",)
        )
        self.assertGreater(len(generated), 0)
        self.assertTrue((out_dir / "resonance_parameters_vs_t.png").is_file())
        self.assertGreater((out_dir / "resonance_parameters_vs_t.png").stat().st_size, 1000)
        self.assertTrue((out_dir / "resonance_pi1_1400_params_vs_t.png").is_file())
        self.assertTrue((out_dir / "resonance_a2_1320_params_vs_t.png").is_file())
        self.assertTrue((out_dir / "resonance_a2_1700_params_vs_t.png").is_file())
        self.assertTrue((out_dir / "parameters_vs_t" / "MD_shared_a2_1320_1700_Mass1_vs_t.png").is_file())
        # Fixed shape parameters and unused poly parameters (poly_order=-1) must be excluded
        self.assertFalse((out_dir / "parameters_vs_t" / "MD_shared_pi1_1400_Width_vs_t.png").exists())
        self.assertFalse((out_dir / "parameters_vs_t" / "MD_shared_a2_1320_1700_Mass2_vs_t.png").exists())
        self.assertFalse((out_dir / "parameters_vs_t" / "MD_a_2_1_a2_1320_1700_a0_vs_t.png").exists())

        # Single t-bin should return empty cleanly
        self.assertEqual(
            plot_resonance_parameters_vs_t(cfg, t_data, ["t010020"], out_dir, formats=("png",)),
            [],
        )

    def test_plot_spd_total_intensity_vs_t(self):
        from moments2amplitudes.schema import MassBinsConfig, DatasetConfig

        masses = [1.0, 1.2, 1.4, 1.6]
        cfg = FitConfig(
            project=None,
            dataset=DatasetConfig(moments_dir="", mass_bins=MassBinsConfig(min=1.0, max=1.6, step=0.2)),
            waves=[
                WaveConfig(l=0, name="S-wave", type=WaveType.PIECEWISE),
                WaveConfig(l=1, name="P-wave", type=WaveType.PIECEWISE),
                WaveConfig(l=2, name="D-wave", type=WaveType.PIECEWISE),
            ],
            plotting=PlottingConfig(formats=["png"]),
        )
        mags = ["a_0_0", "a_1_1", "b_1_1", "a_2_1", "b_2_1", "a_2_2"]
        t_data = {
            "t010020": {
                "name": "t010020",
                "t_center": 0.15,
                "masses": masses,
                "nominal": {
                    "a_0_0": [20.0, 25.0, 22.0, 15.0],
                    "a_1_1": [5.0, 8.0, 7.0, 4.0],
                    "b_1_1": [3.0, 4.0, 3.5, 2.0],
                    "a_2_1": [10.0, 30.0, 45.0, 12.0],
                    "b_2_1": [6.0, 12.0, 15.0, 5.0],
                    "a_2_2": [4.0, 10.0, 18.0, 6.0],
                },
                "mcmc": {
                    "a_0_0": {1.0: [19.0, 21.0], 1.2: [24.0, 26.0], 1.4: [21.0, 23.0], 1.6: [14.0, 16.0]},
                    "a_1_1": {1.0: [4.5, 5.5], 1.2: [7.5, 8.5], 1.4: [6.5, 7.5], 1.6: [3.5, 4.5]},
                    "a_2_1": {1.0: [9.0, 11.0], 1.2: [28.0, 32.0], 1.4: [43.0, 47.0], 1.6: [11.0, 13.0]},
                },
                "h0_data": {"mass": masses, "val": [600.0, 1900.0, 3200.0, 480.0], "err": [30.0, 50.0, 60.0, 25.0]},
            },
            "t020032": {
                "name": "t020032",
                "t_center": 0.26,
                "masses": masses,
                "nominal": {
                    "a_0_0": [15.0, 18.0, 16.0, 10.0],
                    "a_1_1": [4.0, 6.0, 5.5, 3.0],
                    "b_1_1": [2.0, 3.0, 2.5, 1.5],
                    "a_2_1": [8.0, 22.0, 34.0, 9.0],
                    "b_2_1": [4.0, 9.0, 11.0, 3.0],
                    "a_2_2": [3.0, 8.0, 14.0, 4.0],
                },
                "mcmc": {},
                "h0_data": {"mass": masses, "val": [350.0, 1050.0, 1800.0, 290.0], "err": [20.0, 35.0, 45.0, 18.0]},
            },
        }

        out_dir = self.test_dir / "spd_vs_t"
        generated = plot_spd_total_intensity_vs_t(
            cfg, t_data, ["t010020", "t020032"], mags, out_dir, fit_mass_step=0.2, formats=("png",)
        )
        self.assertGreater(len(generated), 0)
        self.assertTrue((out_dir / "spd_waves_total_reflectivity_vs_t.png").is_file())
        self.assertGreater((out_dir / "spd_waves_total_reflectivity_vs_t.png").stat().st_size, 1000)
        self.assertTrue((out_dir / "spd_waves_vs_mass_all_tbins.png").is_file())
        self.assertGreater((out_dir / "spd_waves_vs_mass_all_tbins.png").stat().st_size, 1000)
        self.assertTrue((out_dir / "spd_waves_refl_separated_vs_mass_all_tbins.png").is_file())
        self.assertGreater((out_dir / "spd_waves_refl_separated_vs_mass_all_tbins.png").stat().st_size, 1000)
        self.assertTrue((out_dir / "spd_waves_refl_separated_vs_t.png").is_file())
        self.assertGreater((out_dir / "spd_waves_refl_separated_vs_t.png").stat().st_size, 1000)

        # Single t-bin returns empty
        self.assertEqual(
            plot_spd_total_intensity_vs_t(cfg, t_data, ["t010020"], mags, out_dir, formats=("png",)),
            [],
        )

    def test_compute_relative_phase_data(self):
        from moments2amplitudes.plotter import compute_relative_phase_data
        import numpy as np

        masses = [1.1, 1.2, 1.3]
        nom_phase = [1.0, 2.0, 3.0]
        anchor_phase = [0.5, 1.0, 1.5]
        mcmc_phase = {
            1.1: [1.1, 0.9],
            1.2: [2.1, 1.9],
            1.3: [3.1, 2.9],
        }
        anchor_mcmc = {
            1.1: [0.6, 0.4],
            1.2: [1.1, 0.9],
            1.3: [1.6, 1.4],
        }

        y_rel, mcmc_rel = compute_relative_phase_data(
            masses, nom_phase, anchor_phase, mcmc_phase, anchor_mcmc
        )
        self.assertEqual(len(y_rel), 3)
        self.assertAlmostEqual(y_rel[0], 0.5)
        self.assertAlmostEqual(y_rel[1], 1.0)
        self.assertAlmostEqual(y_rel[2], 1.5)

        self.assertIsNotNone(mcmc_rel)
        # Sample 0 diff: 1.1 - 0.6 = 0.5, Sample 1 diff: 0.9 - 0.4 = 0.5
        self.assertAlmostEqual(mcmc_rel[1.1][0], 0.5)
        self.assertAlmostEqual(mcmc_rel[1.1][1], 0.5)

    def test_resolve_reference_phase_info_two_breit_wigner(self):
        from moments2amplitudes.plotter import resolve_reference_phase_info
        from moments2amplitudes.schema import FitConfig, PlottingConfig, WaveConfig, ResonanceConfig, ModelType, WaveType

        masses = [1.1, 1.2, 1.3, 1.4, 1.5]
        cfg_two_bw = FitConfig(
            project=None,
            dataset=None,
            waves=[
                WaveConfig(
                    l=0,
                    name="S-wave",
                    type=WaveType.PIECEWISE,
                ),
                WaveConfig(
                    l=1,
                    name="P-wave",
                    type=WaveType.MASS_DEPENDENT,
                    model=ModelType.BARRIER_BREIT_WIGNER,
                    resonances=[ResonanceConfig(name="pi1_1400", mass=1.354, width=0.33)],
                ),
                WaveConfig(
                    l=2,
                    name="D-wave",
                    type=WaveType.MASS_DEPENDENT,
                    model=ModelType.BARRIER_TWO_BREIT_WIGNER,
                    resonances=[
                        ResonanceConfig(name="a2_1320", mass=1.318, width=0.107, phase_fixed=True),
                        ResonanceConfig(name="a2_1700", mass=1.706, width=0.38, phase_fixed=False),
                    ],
                ),
            ],
            plotting=PlottingConfig(),
        )

        nominal_amps = {
            "a_0_0": [10.0] * 5,
            "aphi_0_0": [-0.5] * 5,
            "a_1_1": [20.0] * 5,
            "aphi_1_1": [0.8] * 5,
            "a_2_1": [50.0] * 5,
            "aphi_2_1": [1.2] * 5,
            "a_2_2": [30.0] * 5,
            "aphi_2_2": [1.5] * 5,
        }
        active_phases = ["aphi_0_0", "aphi_1_1", "aphi_2_1", "aphi_2_2"]

        ref_info = resolve_reference_phase_info(
            config=cfg_two_bw,
            nominal_amps=nominal_amps,
            active_phases=active_phases,
            masses=masses,
        )

        self.assertTrue(ref_info["is_two_bw"])
        # Reference wave is the a2(1320) partial wave with the biggest magnitude (a_2_1)
        self.assertEqual(ref_info["ref_name"], "a_2_1")
        self.assertIn("aphi_2_1", ref_info["omit_phases"])
        self.assertEqual(len(ref_info["nominal_ref"]["a"]), len(masses))

    def test_resolve_reference_phase_info_single_breit_wigner(self):
        from moments2amplitudes.plotter import resolve_reference_phase_info
        from moments2amplitudes.schema import FitConfig, PlottingConfig, WaveConfig, ResonanceConfig, ModelType, WaveType

        masses = [1.1, 1.2, 1.3, 1.4, 1.5]
        cfg_single_bw = FitConfig(
            project=None,
            dataset=None,
            waves=[
                WaveConfig(
                    l=0,
                    name="S-wave",
                    type=WaveType.PIECEWISE,
                ),
                WaveConfig(
                    l=2,
                    name="D-wave",
                    type=WaveType.MASS_DEPENDENT,
                    model=ModelType.BARRIER_BREIT_WIGNER,
                    resonances=[
                        ResonanceConfig(name="a2_1320", mass=1.318, width=0.107, phase_fixed=True),
                    ],
                ),
            ],
            plotting=PlottingConfig(),
        )

        nominal_amps = {
            "a_0_0": [10.0] * 5,
            "aphi_0_0": [-0.5] * 5,
            "a_2_1": [50.0] * 5,
            "aphi_2_1": [1.2] * 5,
        }
        active_phases = ["aphi_0_0", "aphi_2_1"]

        ref_info = resolve_reference_phase_info(
            config=cfg_single_bw,
            nominal_amps=nominal_amps,
            active_phases=active_phases,
            masses=masses,
        )

        self.assertFalse(ref_info["is_two_bw"])
        self.assertTrue(ref_info["is_single_bw"])
        # USER RULE: For single Breit-Wigner, relative phase is identically 0, so D-waves are omitted
        self.assertIn("aphi_2_1", ref_info["omit_phases"])
        self.assertNotIn("aphi_0_0", ref_info["omit_phases"])

    def test_format_relative_phase_label_resonance_reference(self):
        from moments2amplitudes.plotter import (
            format_relative_phase_symbol,
            format_relative_phase_label,
        )

        sym_d = format_relative_phase_symbol("aphi_2_1", "a2_1320")
        self.assertEqual(sym_d, r"\Delta\phi(D^{+}_{1} - a_{2}(1320))")

        lbl_d = format_relative_phase_label("aphi_2_1", "a2_1320")
        self.assertEqual(lbl_d, r"$\Delta\phi(D^{+}_{1} - a_{2}(1320))$")

        lbl_p = format_relative_phase_label("aphi_1_1", "a2_1320")
        self.assertEqual(lbl_p, r"$\Delta\phi(P^{+}_{1} - a_{2}(1320))$")

        lbl_s = format_relative_phase_label("aphi_0_0", "a2_1320")
        self.assertEqual(lbl_s, r"$\Delta\phi(S^{+}_{0} - a_{2}(1320))$")

    def test_format_reference_partial_wave_label(self):
        from moments2amplitudes.plotter import format_reference_partial_wave_label

        ref_info = {
            "ref_name": "a_2_2",
            "ref_names": {"a": "a_2_2", "b": "b_2_2"},
        }

        # Specific reflectivity
        self.assertEqual(format_reference_partial_wave_label(ref_info, refl_char="a"), "$D^{+}_{2}$")
        self.assertEqual(format_reference_partial_wave_label(ref_info, refl_char="b"), "$D^{-}_{2}$")

        # Both reflectivities active
        phases_both = ["aphi_1_1", "bphi_1_1"]
        self.assertEqual(
            format_reference_partial_wave_label(ref_info, phases=phases_both),
            "$D^{+}_{2}$ (+), $D^{-}_{2}$ (-)",
        )

        # Only positive reflectivity active
        phases_pos = ["aphi_1_1", "aphi_2_1"]
        self.assertEqual(format_reference_partial_wave_label(ref_info, phases=phases_pos), "$D^{+}_{2}$")

        # Custom reference amplitude (e.g. D+_1)
        ref_info_d1 = {
            "ref_name": "a_2_1",
            "ref_names": {"a": "a_2_1", "b": "b_2_1"},
        }
        self.assertEqual(format_reference_partial_wave_label(ref_info_d1, refl_char="a"), "$D^{+}_{1}$")
        self.assertEqual(format_reference_partial_wave_label(ref_info_d1, refl_char="b"), "$D^{-}_{1}$")

    def test_spectroscopic_parsing_and_anchor_resolution(self):
        from moments2amplitudes.plotter import get_amplitude_symbol_parts, resolve_anchor_amplitudes
        from moments2amplitudes.schema import FitConfig, PlottingConfig

        # Spectroscopic parsing
        self.assertEqual(get_amplitude_symbol_parts("D+_2"), ("D", "+", "2", "amp"))
        self.assertEqual(get_amplitude_symbol_parts("D-_1"), ("D", "-", "1", "amp"))
        self.assertEqual(get_amplitude_symbol_parts("P+_1"), ("P", "+", "1", "amp"))
        self.assertEqual(get_amplitude_symbol_parts("$D^{+}_{2}$"), ("D", "+", "2", "amp"))

        # Anchor resolution with spectroscopic string in config
        cfg = FitConfig(
            project=None,
            dataset=None,
            waves=[],
            plotting=PlottingConfig(reference_amplitude="D+_2"),
        )
        nominal_amps = {
            "a_2_2": [10.0, 20.0],
            "aphi_2_2": [1.0, 1.0],
            "b_2_2": [10.0, 20.0],
            "bphi_2_2": [1.0, 1.0],
        }
        active_phases = ["aphi_2_2", "bphi_2_2"]
        anchors = resolve_anchor_amplitudes(cfg, nominal_amps, active_phases)
        self.assertEqual(anchors["a"], "aphi_2_2")
        self.assertEqual(anchors["b"], "bphi_2_2")

    def test_relative_phase_plots_explicit_reference_wave(self):
        from moments2amplitudes.plotter import plot_amplitudes_and_bands
        import json
        import matplotlib.pyplot as plt

        cfg = FitConfig()
        out_dir = self.test_dir / "plots_ref_check"
        out_dir.mkdir(parents=True, exist_ok=True)
        cache_dir = out_dir / ".cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        json_path = cache_dir / "amplitudes_data.json"

        data = {
            "masses": [1.1, 1.2, 1.3, 1.4],
            "nominal": {
                "a_2_1": [10.0, 20.0, 30.0, 40.0],
                "aphi_2_1": [0.1, 0.2, 0.3, 0.4],
                "a_0_0": [5.0, 8.0, 12.0, 15.0],
                "aphi_0_0": [0.5, 0.6, 0.7, 0.8],
            },
        }
        with open(json_path, "w") as f:
            json.dump(data, f)

        captured_titles = {}
        captured_ylabels = {}

        orig_savefig = plt.Figure.savefig
        def custom_savefig(self, fname, *args, **kwargs):
            fname_str = str(fname)
            for ax in self.axes:
                title = ax.get_title()
                ylabel = ax.get_ylabel()
                if title:
                    captured_titles.setdefault(fname_str, []).append(title)
                if ylabel:
                    captured_ylabels.setdefault(fname_str, []).append(ylabel)
            return orig_savefig(self, fname, *args, **kwargs)

        dummy_nom = self.test_dir / "dummy_nom_ref.root"
        dummy_nom.touch()

        with patch("subprocess.run") as mock_run, patch.object(plt.Figure, "savefig", custom_savefig):
            mock_run.return_value = MagicMock(returncode=0)
            plot_amplitudes_and_bands(dummy_nom, None, out_dir, cfg, formats=("png",))

        # Check amplitudes_phase_vs_mass.png
        phase_vs_mass_key = str(out_dir / "amplitudes_phase_vs_mass.png")
        self.assertIn(phase_vs_mass_key, captured_titles)
        self.assertTrue(any("Reference Partial-Wave: $D^{+}_{1}$" in t for t in captured_titles[phase_vs_mass_key]))
        self.assertTrue(any("Reference Partial-Wave: $D^{+}_{1}$" in y for y in captured_ylabels[phase_vs_mass_key]))

        # Check amplitudes_magnitude_vs_mass.png bottom panel
        mag_vs_mass_key = str(out_dir / "amplitudes_magnitude_vs_mass.png")
        self.assertIn(mag_vs_mass_key, captured_ylabels)
        self.assertTrue(any("Ref Partial-Wave: $D^{+}_{1}$" in y for y in captured_ylabels[mag_vs_mass_key]))

        # Check individual combo plot a_0_0_intensity_and_phase.png
        combo_key = str(out_dir / "amplitudes" / "individual" / "a_0_0_intensity_and_phase.png")
        self.assertIn(combo_key, captured_ylabels)
        self.assertTrue(any("Ref Partial-Wave: $D^{+}_{1}$" in y for y in captured_ylabels[combo_key]))

        # Check individual phase plot aphi_0_0_phase.png
        indiv_phi_key = str(out_dir / "amplitudes" / "individual" / "aphi_0_0_phase.png")
        self.assertIn(indiv_phi_key, captured_titles)
        self.assertTrue(any("Reference Partial-Wave: $D^{+}_{1}$" in t for t in captured_titles[indiv_phi_key]))
        self.assertTrue(any("Reference Partial-Wave: $D^{+}_{1}$" in y for y in captured_ylabels[indiv_phi_key]))


if __name__ == "__main__":
    unittest.main()


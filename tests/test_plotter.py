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
                "aphi_0_0": [0.0, 0.0, 0.0, 0.0],
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
        self.assertIn("a_2_1_intensity_and_phase.png", plot_names)
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

        # 2. Test continuous posterior density band (density_band=True)
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
        self.assertGreater(len(ax2.collections), 10)
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


if __name__ == "__main__":
    unittest.main()


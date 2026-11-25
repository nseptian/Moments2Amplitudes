"""
Tests for Moments2Amplitudes data models and schema definitions.
"""

import unittest
from moments2amplitudes.schema import (
    FitConfig,
    ProjectConfig,
    DatasetConfig,
    MassBinsConfig,
    Chi2Config,
    WaveConfig,
    ResonanceConfig,
    MultiStartConfig,
    parse_seed_range,
    MCMCPropagationConfig,
    PlottingConfig,
    EnvironmentConfig,
    WaveType,
    ModelType,
    InitStrategy,
)


class TestSchema(unittest.TestCase):
    def test_default_schema_instantiation(self):
        cfg = FitConfig()
        self.assertEqual(cfg.project.name, "PartialWaveFit")
        self.assertTrue(cfg.dataset.use_covariance)
        self.assertEqual(cfg.chi2.max_l, 4)
        self.assertEqual(cfg.multi_start.n_starts, 100)
        self.assertFalse(cfg.mcmc_propagation.enabled)

    def test_mass_bins_generation(self):
        bin_cfg = MassBinsConfig(min=1.0, max=1.2, step=0.05)
        bins = bin_cfg.get_bin_centers()
        self.assertEqual(len(bins), 5)
        self.assertAlmostEqual(bins[0], 1.00, places=4)
        self.assertAlmostEqual(bins[-1], 1.20, places=4)

        explicit_cfg = MassBinsConfig(bins=[1.0, 1.05, 1.15, 1.25])
        self.assertEqual(explicit_cfg.get_bin_centers(), [1.0, 1.05, 1.15, 1.25])

    def test_wave_configuration_types(self):
        w_pw = WaveConfig(l=0, name="S0", type=WaveType.PIECEWISE)
        self.assertEqual(w_pw.type, WaveType.PIECEWISE)
        self.assertEqual(w_pw.l, 0)

        w_md = WaveConfig(
            l=2,
            name="D1",
            type=WaveType.MASS_DEPENDENT,
            model=ModelType.BREIT_WIGNER,
            resonances=[ResonanceConfig(name="a2_1320", mass=1.318, width=0.107)],
        )
        self.assertEqual(w_md.type, WaveType.MASS_DEPENDENT)
        self.assertEqual(len(w_md.resonances), 1)
        self.assertEqual(w_md.resonances[0].name, "a2_1320")

    def test_barrier_wave_configuration(self):
        w_barrier = WaveConfig(
            l=2,
            name="D2",
            type=WaveType.MASS_DEPENDENT,
            model=ModelType.BARRIER_BREIT_WIGNER,
            barrier=True,
            barrier_l=2,
            resonances=[ResonanceConfig(name="a2_1320", mass=1.318, width=0.107)],
        )
        self.assertEqual(w_barrier.model, ModelType.BARRIER_BREIT_WIGNER)
        self.assertTrue(w_barrier.barrier)
        self.assertEqual(w_barrier.barrier_l, 2)

    def test_invalid_mass_binning(self):
        with self.assertRaises(ValueError):
            MassBinsConfig().get_bin_centers()

    def test_flatte_resonance_configuration(self):
        res = ResonanceConfig(
            name="a0_980",
            mass=0.980,
            g_etapi=0.5,
            g_etapi_fixed=False,
            g_etapi_limits=[0.1, 1.5],
            g_KK=-99.0,
            g_KK_fixed=True,
            g_etaprimepi=-99.0,
            g_etaprimepi_fixed=True,
        )
        self.assertEqual(res.mass, 0.980)
        self.assertEqual(res.g_etapi, 0.5)
        self.assertEqual(res.g_s, 0.5)
        self.assertFalse(res.g_etapi_fixed)
        self.assertEqual(res.g_etapi_limits, [0.1, 1.5])
        self.assertEqual(res.g_KK, -99.0)
        self.assertTrue(res.g_KK_fixed)
        self.assertEqual(res.g_etaprimepi, -99.0)
        self.assertTrue(res.g_etaprimepi_fixed)
        self.assertIsNone(res.width)

    def test_plotting_config(self):
        plot_cfg = PlottingConfig(n_lowest_chi2=15, moments=["H_0_0_0", "H_0_2_0"])
        self.assertTrue(plot_cfg.enabled)
        self.assertEqual(plot_cfg.n_lowest_chi2, 15)
        self.assertEqual(plot_cfg.moments, ["H_0_0_0", "H_0_2_0"])
        self.assertTrue(plot_cfg.draw_best_fits)
        self.assertEqual(plot_cfg.formats, ["png", "pdf"])
        self.assertTrue(plot_cfg.amplitudes_vs_t)
        self.assertFalse(plot_cfg.amplitudes_vs_t_sharey)


    def test_chi2_config_included_L(self):
        c1 = Chi2Config(included_L=[0, 1, 2, 3, 4])
        self.assertEqual(c1.included_L, [0, 1, 2, 3, 4])

        c2 = Chi2Config(included_L=[4])
        self.assertEqual(c2.included_L, [4])

        # Mutation of moments L list
        c2.included_L = [0, 2, 4]
        self.assertEqual(c2.included_L, [0, 2, 4])

    def test_parse_seed_range(self):
        self.assertEqual(parse_seed_range([5001, 10000]), (5001, 10000))
        self.assertEqual(parse_seed_range((5001, 10000)), (5001, 10000))
        self.assertEqual(parse_seed_range("5001:10000"), (5001, 10000))
        self.assertEqual(parse_seed_range("5001-10000"), (5001, 10000))
        self.assertEqual(parse_seed_range("5001..10000"), (5001, 10000))
        self.assertEqual(parse_seed_range("5001, 10000"), (5001, 10000))

        # Single seed range
        self.assertEqual(parse_seed_range("5001:5001"), (5001, 5001))

        # Errors
        with self.assertRaises(ValueError):
            parse_seed_range("10000:5001")  # inverted
        with self.assertRaises(ValueError):
            parse_seed_range("0:100")       # 0 not allowed (1-based)
        with self.assertRaises(ValueError):
            parse_seed_range("-5:100")      # negative
        with self.assertRaises(ValueError):
            parse_seed_range("invalid")     # non-numeric

    def test_multi_start_seed_range_generation(self):
        # Default 1..100
        ms_default = MultiStartConfig(n_starts=100)
        self.assertEqual(ms_default.get_seed_range(), (1, 100))
        self.assertEqual(len(ms_default.get_seeds()), 100)
        self.assertEqual(ms_default.get_seeds()[0], 1)
        self.assertEqual(ms_default.get_seeds()[-1], 100)

        # Explicit seed_start and seed_end
        ms_explicit = MultiStartConfig(seed_start=5001, seed_end=5010)
        self.assertEqual(ms_explicit.get_seed_range(), (5001, 5010))
        self.assertEqual(ms_explicit.get_seeds(), list(range(5001, 5011)))

        # Explicit seed_range string
        ms_str = MultiStartConfig(seed_range="5001:5005")
        self.assertEqual(ms_str.get_seed_range(), (5001, 5005))
        self.assertEqual(ms_str.get_seeds(), [5001, 5002, 5003, 5004, 5005])

        # Override n_starts takes precedence when no range specified
        ms_override = MultiStartConfig(n_starts=50)
        self.assertEqual(len(ms_override.get_seeds(n_starts_override=20)), 20)
        self.assertEqual(ms_override.get_seed_range(n_starts_override=20), (1, 20))

        # skip_existing
        ms_skip = MultiStartConfig(skip_existing=True)
        self.assertTrue(ms_skip.skip_existing)


if __name__ == "__main__":
    unittest.main()

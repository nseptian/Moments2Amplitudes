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
    MCMCPropagationConfig,
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

    def test_invalid_mass_binning(self):
        with self.assertRaises(ValueError):
            MassBinsConfig().get_bin_centers()


if __name__ == "__main__":
    unittest.main()

"""
Tests for reflectivity-separated resonance cross-section evaluation.
"""

import unittest
import math
from moments2amplitudes.schema import CrossSectionConfig
from moments2amplitudes.cross_section import (
    parse_t_bin_from_string,
    compute_bw_integral,
    WaveResonanceYield,
    CrossSectionResult,
    format_cross_section_table,
    propagate_mcmc_cross_section,
    find_mcmc_root_file,
)
import numpy as np


class TestCrossSection(unittest.TestCase):
    def test_parse_t_bin(self):
        self.assertEqual(parse_t_bin_from_string("t010020"), (0.10, 0.20))
        self.assertEqual(parse_t_bin_from_string("run/t020032/Mom2Amp"), (0.20, 0.325))
        self.assertEqual(parse_t_bin_from_string("MomentFit_m080200/t032050"), (0.325, 0.50))
        self.assertEqual(parse_t_bin_from_string("t050075"), (0.50, 0.75))
        self.assertEqual(parse_t_bin_from_string("t075100"), (0.75, 1.00))
        self.assertEqual(parse_t_bin_from_string("0.15:0.35"), (0.15, 0.35))
        self.assertIsNone(parse_t_bin_from_string("arbitrary_string"))

    def test_compute_bw_integral(self):
        # Integral of BW with k=10.0, M=1.318, Gamma=0.107 over [0.8, 2.0] / 0.04
        integral = compute_bw_integral(10.0, 1.318, 0.107, 0.8, 2.0, 0.04)
        self.assertGreater(integral, 1000.0)
        # Yield scales strictly as k^2
        integral_2x = compute_bw_integral(20.0, 1.318, 0.107, 0.8, 2.0, 0.04)
        self.assertAlmostEqual(integral_2x / integral, 4.0, places=4)

    def test_cross_section_result_math(self):
        res = CrossSectionResult(
            resonance="a2_1320",
            mass=1.318,
            mass_error=0.001,
            width=0.107,
            width_error=0.002,
            t_min=0.1,
            t_max=0.2,
            t_center=0.15,
            delta_t=0.1,
            luminosity=257.4,
            branching_ratio=0.0565,
            den_factor=1e6,
            pos_refl_yield=100000.0,
            pos_refl_yield_err=5000.0,
            neg_refl_yield=50000.0,
            neg_refl_yield_err=2500.0,
            total_yield=150000.0,
            total_yield_err=math.sqrt(5000**2 + 2500**2),
            cs_pos=0.06876,
            cs_pos_err=0.00344,
            cs_neg=0.03438,
            cs_neg_err=0.00172,
            cs_total=0.10314,
            cs_total_err=0.00385,
            pos_refl_fraction=100000.0 / 150000.0,
            wave_yields={
                "a_2_1": WaveResonanceYield("a_2_1", 1, 10.0, 0.5, 50000.0, 5000.0),
                "b_2_1": WaveResonanceYield("b_2_1", -1, 7.07, 0.35, 50000.0, 4950.0),
            }
        )

        self.assertAlmostEqual(res.pos_refl_fraction, 2.0 / 3.0, places=4)
        tbl = format_cross_section_table(res)
        self.assertIn("A2_1320", tbl)
        self.assertIn("Positive Reflectivity (+)", tbl)
        self.assertIn("Negative Reflectivity (-)", tbl)
        self.assertIn("Total Cross Section", tbl)

    def test_propagate_mcmc_cross_section(self):
        res = CrossSectionResult(
            resonance="a2_1320",
            mass=1.318,
            mass_error=0.001,
            width=0.107,
            width_error=0.002,
            t_min=0.1,
            t_max=0.2,
            t_center=0.15,
            delta_t=0.1,
            luminosity=257.4,
            branching_ratio=0.0565,
            den_factor=1e6,
            pos_refl_yield=100000.0,
            pos_refl_yield_err=5000.0,
            neg_refl_yield=50000.0,
            neg_refl_yield_err=2500.0,
            total_yield=150000.0,
            total_yield_err=5590.0,
            cs_pos=0.06876,
            cs_pos_err=0.00344,
            cs_neg=0.03438,
            cs_neg_err=0.00172,
            cs_total=0.10314,
            cs_total_err=0.00385,
            pos_refl_fraction=100000.0 / 150000.0,
            wave_yields={
                "a_2_1": WaveResonanceYield("a_2_1", 1, 10.0, 0.5, 50000.0, 5000.0),
                "b_2_1": WaveResonanceYield("b_2_1", -1, 7.07, 0.35, 50000.0, 4950.0),
            }
        )

        # Generate synthetic MCMC samples (N=100)
        np.random.seed(42)
        n_samples = 100
        mock_mcmc = {
            "resonance": "a2_1320",
            "entries_count": n_samples,
            "waves": ["a_2_1", "b_2_1"],
            "reflectivities": [1, -1],
            "is_valid": [1] * n_samples,
            "status": [0] * n_samples,
            "mass": [1.318] * n_samples,
            "width": [0.107] * n_samples,
            "couplings": {
                "a_2_1": list(np.random.normal(10.0, 0.4, n_samples)),
                "b_2_1": list(np.random.normal(7.07, 0.3, n_samples)),
            }
        }

        updated_res = propagate_mcmc_cross_section(res, mock_mcmc, mass_limits=(0.8, 2.0), bin_width=0.04)
        self.assertTrue(updated_res.has_mcmc)
        self.assertEqual(updated_res.mcmc_samples_count, 100)
        self.assertGreater(updated_res.cs_pos_mcmc_err, 0.0)
        self.assertGreater(updated_res.cs_neg_mcmc_err, 0.0)
        self.assertGreater(updated_res.cs_total_mcmc_err, 0.0)
        self.assertGreater(updated_res.pos_refl_fraction_mcmc_err, 0.0)
        self.assertGreater(updated_res.cs_total_mcmc_err_lo, 0.0)
        self.assertGreater(updated_res.cs_total_mcmc_err_hi, 0.0)

        # Check wave yield errors updated
        self.assertGreater(updated_res.wave_yields["a_2_1"].yield_mcmc_err, 0.0)
        self.assertGreater(updated_res.wave_yields["b_2_1"].yield_mcmc_err, 0.0)

        # Check ASCII table formatting includes MCMC lines
        tbl = format_cross_section_table(updated_res)
        self.assertIn("MCMC Propagation:     Enabled (100 samples", tbl)
        self.assertIn("(MCMC)", tbl)
        self.assertIn("MCMC 68% CI", tbl)

    def test_mcmc_invalid_filtering(self):
        res = CrossSectionResult(
            resonance="a2_1320", mass=1.318, mass_error=0.0, width=0.107, width_error=0.0,
            t_min=0.1, t_max=0.2, t_center=0.15, delta_t=0.1, luminosity=257.4,
            branching_ratio=0.0565, den_factor=1e6, pos_refl_yield=1000.0, pos_refl_yield_err=10.0,
            neg_refl_yield=500.0, neg_refl_yield_err=5.0, total_yield=1500.0, total_yield_err=15.0,
            cs_pos=0.001, cs_pos_err=0.0001, cs_neg=0.0005, cs_neg_err=0.00005,
            cs_total=0.0015, cs_total_err=0.00015, pos_refl_fraction=0.67,
            wave_yields={"a_2_1": WaveResonanceYield("a_2_1", 1, 10.0, 0.1, 1000.0, 10.0)}
        )

        mock_mcmc = {
            "resonance": "a2_1320",
            "entries_count": 5,
            "waves": ["a_2_1"],
            "reflectivities": [1],
            "is_valid": [1, 0, 1, 1, 0],
            "status": [0, 0, 0, 1, 0],
            "mass": [1.318] * 5,
            "width": [0.107] * 5,
            "couplings": {"a_2_1": [10.0, 999.0, 10.2, 888.0, 777.0]}
        }

        updated = propagate_mcmc_cross_section(res, mock_mcmc, filter_valid=True)
        # Only entry 0 and 2 have is_valid=1 and status=0 -> 2 valid samples
        self.assertEqual(updated.mcmc_samples_count, 2)


if __name__ == "__main__":
    unittest.main()


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
)


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


if __name__ == "__main__":
    unittest.main()

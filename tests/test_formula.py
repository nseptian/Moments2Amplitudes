"""
Tests for mathematical formula translation into C++ expressions.
"""

import unittest
from moments2amplitudes.driver import translate_formula_to_cpp


class TestFormula(unittest.TestCase):
    def test_translate_breit_wigner(self):
        formula = "k / (M^2 - s - i * M * width)"
        params = ["k", "M", "width"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertTrue(any("const double s = mass * mass;" in line for line in param_lines))
        self.assertTrue(any("const double k = (p.size() > 0)" in line for line in param_lines))
        self.assertEqual(expr, "k / std::complex<double>((M * M) - s, -M * width)")

    def test_translate_two_breit_wigner(self):
        formula = "k1 / (M1^2 - s - i*M1*w1) + (k2 / (M2^2 - s - i*M2*w2)) * exp(i*phase)"
        params = ["k1", "M1", "w1", "k2", "M2", "w2", "phase"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("std::complex<double>((M1 * M1) - s, -M1 * w1)", expr)
        self.assertIn("std::complex<double>((M2 * M2) - s, -M2 * w2)", expr)
        self.assertIn("std::polar(1.0, phase)", expr)

    def test_translate_flatte(self):
        formula = "k_fl / (M0^2 - s - i * (g_etapi^2 * rho_etapi + g_KK^2 * rho_KK + g_etaprimepi^2 * rho_etaprimepi))"
        params = ["k_fl", "M0", "g_etapi", "g_KK", "g_etaprimepi"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertTrue(any("PhaseSpaceFactor" in line for line in param_lines))
        self.assertTrue(any("gamma_total" in line for line in param_lines))
        self.assertIn("std::complex<double>((M0 * M0) - s, -gamma_total)", expr)

    def test_translate_breit_wigner_with_kinematic_factor(self):
        formula = "sqrt(rho_etapi) * k / (M^2 - s - i * M * width)"
        params = ["k", "M", "width"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("std::sqrt(m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.547853, 0.13957))", expr)
        self.assertIn("k / std::complex<double>((M * M) - s, -M * width)", expr)

    def test_translate_breit_wigner_with_kinematic_factor_parentheses(self):
        formula = "sqrt(rho_etapi) * (k / (M^2 - s - i * M * width))"
        params = ["k", "M", "width"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("std::sqrt(m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.547853, 0.13957))", expr)
        self.assertIn("k / std::complex<double>((M * M) - s, -M * width)", expr)

    def test_translate_two_breit_wigner_with_kinematic_factor(self):
        formula = "sqrt(rho_etapi) * (k1 / (M1^2 - s - i*M1*w1) + (k2 / (M2^2 - s - i*M2*w2)) * exp(i*phase))"
        params = ["k1", "M1", "w1", "k2", "M2", "w2", "phase"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("std::sqrt(m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.547853, 0.13957))", expr)
        self.assertIn("std::complex<double>((M1 * M1) - s, -M1 * w1)", expr)
        self.assertIn("std::complex<double>((M2 * M2) - s, -M2 * w2)", expr)
        self.assertIn("std::polar(1.0, phase)", expr)


if __name__ == "__main__":
    unittest.main()

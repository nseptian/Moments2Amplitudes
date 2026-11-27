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

    def test_translate_flatte_with_phase(self):
        formula = "k / (Mass^2 - s - i * (g_etapi^2 * rho_etapi + g_KK^2 * rho_KK + g_etaprimepi^2 * rho_etaprimepi)) * exp(i * phase)"
        params = ["k", "Mass", "g_etapi", "g_KK", "g_etaprimepi", "phase"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertTrue(any("gamma_total" in line for line in param_lines))
        self.assertIn("std::polar(1.0, phase)", expr)

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

    def test_translate_two_breit_wigner_two_phases(self):
        formula = "(k1 / (Mass1^2 - s - i * Mass1 * Width1)) * exp(i * phase1) + (k2 / (Mass2^2 - s - i * Mass2 * Width2)) * exp(i * phase2)"
        params = ["k1", "Mass1", "Width1", "phase1", "k2", "Mass2", "Width2", "phase2"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("std::complex<double>((Mass1 * Mass1) - s, -Mass1 * Width1)", expr)
        self.assertIn("std::polar(1.0, phase1)", expr)
        self.assertIn("std::complex<double>((Mass2 * Mass2) - s, -Mass2 * Width2)", expr)
        self.assertIn("std::polar(1.0, phase2)", expr)

    def test_translate_two_breit_wigner_two_phases_with_kinematic_factor(self):
        formula = "sqrt(rho_etapi) * ((k1 / (Mass1^2 - s - i * Mass1 * Width1)) * exp(i * phase1) + (k2 / (Mass2^2 - s - i * Mass2 * Width2)) * exp(i * phase2))"
        params = ["k1", "Mass1", "Width1", "phase1", "k2", "Mass2", "Width2", "phase2"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("std::sqrt(m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.547853, 0.13957))", expr)
        self.assertIn("std::polar(1.0, phase1)", expr)
        self.assertIn("std::polar(1.0, phase2)", expr)

    def test_translate_barrier_breit_wigner(self):
        formula = "barrier_bw_etapi(k, M, width, 2)"
        params = ["k", "M", "width"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("m2pw::MassDependentFunction::BarrierBreitWigner(mass, k, M, width, 2)", expr)

    def test_translate_barrier_breit_wigner_with_kinematic_factor(self):
        formula = "sqrt(rho_etapi) * barrier_bw_etapi(k, M, width, 2)"
        params = ["k", "M", "width"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("std::sqrt(m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.547853, 0.13957))", expr)
        self.assertIn("m2pw::MassDependentFunction::BarrierBreitWigner(mass, k, M, width, 2)", expr)

    def test_translate_two_barrier_breit_wigner_two_phases(self):
        formula = "barrier_bw_etapi(k1, Mass1, Width1, 2) * exp(i * phase1) + barrier_bw_etapi(k2, Mass2, Width2, 2) * exp(i * phase2)"
        params = ["k1", "Mass1", "Width1", "phase1", "k2", "Mass2", "Width2", "phase2"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("m2pw::MassDependentFunction::BarrierBreitWigner(mass, k1, Mass1, Width1, 2)", expr)
        self.assertIn("m2pw::MassDependentFunction::BarrierBreitWigner(mass, k2, Mass2, Width2, 2)", expr)
        self.assertIn("std::polar(1.0, phase1)", expr)
        self.assertIn("std::polar(1.0, phase2)", expr)

    def test_translate_two_barrier_breit_wigner_two_phases_with_kinematic(self):
        formula = "sqrt(rho_etapi) * (barrier_bw_etapi(k1, Mass1, Width1, 2) * exp(i * phase1) + barrier_bw_etapi(k2, Mass2, Width2, 2) * exp(i * phase2))"
        params = ["k1", "Mass1", "Width1", "phase1", "k2", "Mass2", "Width2", "phase2"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("std::sqrt(m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.547853, 0.13957))", expr)
        self.assertIn("m2pw::MassDependentFunction::BarrierBreitWigner(mass, k1, Mass1, Width1, 2)", expr)
        self.assertIn("m2pw::MassDependentFunction::BarrierBreitWigner(mass, k2, Mass2, Width2, 2)", expr)
        self.assertIn("std::polar(1.0, phase1)", expr)
        self.assertIn("std::polar(1.0, phase2)", expr)

    def test_translate_barrier_ratio_and_dynamic_width(self):
        formula = "barrier_ratio_etapi(M, 2) * k / (M^2 - s - i * M * dynamic_width_etapi(M, width, 2))"
        params = ["k", "M", "width"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        self.assertIn("m2pw::MassDependentFunction::BlattWeisskopfBarrierRatio(mass, M, 2)", expr)
        self.assertIn("m2pw::MassDependentFunction::DynamicBreitWignerWidth(mass, M, width, 2)", expr)

    def test_translate_peak_normalized_and_constant_width_barrier_bw(self):
        for formula, cpp_fn in [
            ("sqrt(rho_etapi) * peak_barrier_bw_etapi(k, M, width, 2)", "PeakNormalizedBarrierBreitWigner"),
            ("sqrt(rho_etapi) * const_barrier_bw_etapi(k, M, width, 2)", "ConstantWidthBarrierBreitWigner"),
            ("sqrt(rho_etapi) * peak_const_barrier_bw_etapi(k, M, width, 2)", "PeakNormalizedConstantWidthBarrierBreitWigner"),
        ]:
            param_lines, expr = translate_formula_to_cpp(formula, ["k", "M", "width"])
            self.assertIn(f"m2pw::MassDependentFunction::{cpp_fn}(mass, k, M, width, 2)", expr)
            self.assertIn("PhaseSpaceFactor(mass, 0.547853, 0.13957)", expr)

    def test_translate_mixed_two_barrier_breit_wigner_peak_const(self):
        formula = "sqrt(rho_etapi) * (barrier_bw_etapi(k1, Mass1, Width1, 2) * exp(i * phase1) + peak_const_barrier_bw_etapi(k2, Mass2, Width2, 2) * exp(i * phase2))"
        params = ["k1", "Mass1", "Width1", "phase1", "k2", "Mass2", "Width2", "phase2"]
        param_lines, expr = translate_formula_to_cpp(formula, params)
        joined = "\n".join(param_lines)
        self.assertIn("const double k2 = (p.size() > 3) ? p[3] : 0.0;", joined)
        self.assertIn("const double phase1 = (p.size() > 7) ? p[6] : 0.0;", joined)
        self.assertIn("m2pw::MassDependentFunction::BarrierBreitWigner(mass, k1, Mass1, Width1, 2)", expr)
        self.assertIn("m2pw::MassDependentFunction::PeakNormalizedConstantWidthBarrierBreitWigner(mass, k2, Mass2, Width2, 2)", expr)
        self.assertIn("std::polar(1.0, phase1)", expr)
        self.assertIn("std::polar(1.0, phase2)", expr)


if __name__ == "__main__":
    unittest.main()

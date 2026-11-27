"""
Tests for FitSummary formatting and label transformations.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from moments2amplitudes.summary import format_parameter_label, FitSummary
from moments2amplitudes.rundir import RunDirectory
from moments2amplitudes.schema import FitConfig


class TestSummary(unittest.TestCase):
    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_format_parameter_label(self):
        self.assertEqual(format_parameter_label("MD_a2_1320_Mass"), "Mass (a2_1320)")
        self.assertEqual(format_parameter_label("MD_a2_1320_Width"), "Width (a2_1320)")
        self.assertEqual(format_parameter_label("MD_a_2_1_a2_1320_k"), "k [a_2_1] (a2_1320)")
        self.assertEqual(format_parameter_label("MI_1.3200_a_0_0"), "a_0_0 @ 1.3200 GeV")
        self.assertEqual(format_parameter_label("simple_param"), "simple_param")

    def test_summary_ascii_and_markdown(self):
        cfg = FitConfig()
        cfg.project.name = "test_summary_proj"
        rundir = RunDirectory(self.temp_dir / "summary_run").ensure_directories()

        summary = FitSummary(cfg, rundir)
        summary.chi2 = 123.456
        summary.status = 0
        summary.is_valid = True
        summary.is_valid_error = True
        summary.parameters = [
            {
                "name": "MD_a2_1320_Mass",
                "label": "Mass (a2_1320)",
                "value": 1.3182,
                "error": 0.0012,
                "mcmc_median": 1.3181,
                "mcmc_err_lo": 0.0011,
                "mcmc_err_hi": 0.0013,
                "status": "Free",
            },
            {
                "name": "MD_a2_1320_Width",
                "label": "Width (a2_1320)",
                "value": 0.1070,
                "error": 0.0,
                "mcmc_median": None,
                "mcmc_err_lo": None,
                "mcmc_err_hi": None,
                "status": "Fixed",
            },
        ]

        ascii_out = summary.generate_ascii_table()
        self.assertIn("Moments2Amplitudes Fit Summary Report", ascii_out)
        self.assertIn("123.4560", ascii_out)
        self.assertIn("Mass (a2_1320)", ascii_out)
        self.assertIn("fixed", ascii_out)

        md_out = summary.generate_markdown_report()
        self.assertIn("# Fit Summary: test_summary_proj", md_out)
        self.assertIn("| **Mass (a2_1320)** | `MD_a2_1320_Mass` | `1.31820` |", md_out)

    def test_filter_poly_when_poly_order_minus_one(self):
        from moments2amplitudes.schema import WaveConfig, WaveType, ResonanceConfig
        cfg = FitConfig()
        cfg.waves = [
            WaveConfig(
                l=2,
                name="D-wave",
                type=WaveType.MASS_DEPENDENT,
                poly_order=-1,
                resonances=[ResonanceConfig(name="a2_1320")]
            )
        ]
        rundir = RunDirectory(self.temp_dir / "poly_test").ensure_directories()
        summary = FitSummary(cfg, rundir)

        # Poly parameters should be filtered out
        self.assertFalse(summary._should_display_parameter("MD_b_2_0_a2_1320_re_0"))
        self.assertFalse(summary._should_display_parameter("MD_b_2_0_a2_1320_im_0"))
        self.assertFalse(summary._should_display_parameter("MD_b_2_0_a2_1320_m_threshold"))
        self.assertFalse(summary._should_display_parameter("MD_b_2_0_a2_1320_m_expansion"))

        # Physical parameters should be kept
        self.assertTrue(summary._should_display_parameter("MD_b_2_0_a2_1320_k"))
        self.assertTrue(summary._should_display_parameter("MD_shared_a2_1320_Mass"))
        self.assertTrue(summary._should_display_parameter("MD_shared_a2_1320_Width"))

        txt_file, md_file = summary.save_reports()
        self.assertTrue(txt_file.is_file())
        self.assertTrue(md_file.is_file())

    def test_summary_seed_range_and_plateau(self):
        cfg = FitConfig()
        cfg.project.name = "seed_range_summary_proj"
        cfg.multi_start.seed_range = [5001, 10000]
        rundir = RunDirectory(self.temp_dir / "seed_summary_run").ensure_directories()

        summary = FitSummary(cfg, rundir)
        summary.chi2 = 85.2
        summary.total_starts = 10000
        summary.converged_starts = 9800
        summary.plateau_count = 450
        summary.status = 0
        summary.is_valid = True

        ascii_out = summary.generate_ascii_table()
        self.assertIn("10000", ascii_out)
        self.assertIn("450", ascii_out)

        md_out = summary.generate_markdown_report()
        self.assertIn("10000", md_out)
        self.assertIn("450", md_out)

    def test_summary_at_limit_marked_not_converged(self):
        cfg = FitConfig()
        cfg.project.name = "limit_summary_proj"
        rundir = RunDirectory(self.temp_dir / "limit_summary_run").ensure_directories()

        summary = FitSummary(cfg, rundir)
        summary.chi2 = 250.0
        summary.status = 0
        summary.is_valid = False
        summary.at_limit = True

        ascii_out = summary.generate_ascii_table()
        self.assertIn("NOT CONVERGED (AT LIMIT, Status 0)", ascii_out)

        md_out = summary.generate_markdown_report()
        self.assertIn("Failed (Hit Limit, Status 0)", md_out)

        # Verify load_from_log marks at_limit=1 as invalid even if is_valid=1 appeared
        log_file = self.temp_dir / "fit_seed_1.log"
        log_file.write_text(
            "[FIT][MIN] seed=1 is_valid=1 is_valid_error=0 status=0 at_limit=1\n"
            "[FIT][RESULT] param=MD_shared_a2_1320_Width value=0.25 error=0.01\n"
        )
        summary2 = FitSummary(cfg, rundir)
        self.assertTrue(summary2.load_from_log(log_file))
        self.assertTrue(summary2.at_limit)
        self.assertFalse(summary2.is_valid)


if __name__ == "__main__":
    unittest.main()

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

        txt_file, md_file = summary.save_reports()
        self.assertTrue(txt_file.is_file())
        self.assertTrue(md_file.is_file())


if __name__ == "__main__":
    unittest.main()

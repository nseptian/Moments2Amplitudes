"""
Tests for publication plotter engine and moment comparison logic.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from moments2amplitudes.schema import FitConfig, PlottingConfig
from moments2amplitudes.rundir import RunDirectory
from moments2amplitudes.plotter import discover_lowest_chi2_fits


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

        fits = discover_lowest_chi2_fits(self.rundir, n_lowest=3)
        self.assertEqual(len(fits), 3)
        self.assertAlmostEqual(fits[0][0], 75.1)
        self.assertAlmostEqual(fits[1][0], 85.2)
        self.assertAlmostEqual(fits[2][0], 95.0)


if __name__ == "__main__":
    unittest.main()

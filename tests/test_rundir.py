"""
Tests for RunDirectory management and filesystem layout.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from moments2amplitudes.rundir import RunDirectory, get_run_directory
from moments2amplitudes.schema import FitConfig


class TestRunDirectory(unittest.TestCase):
    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_rundir_hierarchy(self):
        target = self.temp_dir / "test_run"
        rundir = RunDirectory(target)
        self.assertEqual(rundir.root, target)
        self.assertEqual(rundir.fit_dir, target / "fit")
        self.assertEqual(rundir.fit_starts_dir, target / "fit" / "starts")
        self.assertEqual(rundir.mcmc_dir, target / "mcmc")
        self.assertEqual(rundir.mcmc_workers_dir, target / "mcmc" / "workers")
        self.assertEqual(rundir.plots_dir, target / "plots")
        self.assertEqual(rundir.reports_dir, target / "reports")

        rundir.ensure_directories()
        self.assertTrue(rundir.fit_starts_dir.is_dir())
        self.assertTrue(rundir.mcmc_workers_dir.is_dir())
        self.assertTrue(rundir.plots_dir.is_dir())
        self.assertTrue(rundir.reports_dir.is_dir())

    def test_rundir_file_resolution(self):
        rundir = RunDirectory(self.temp_dir / "run1")
        self.assertEqual(rundir.get_seed_root_file(42).name, "fit_seed_42.root")
        self.assertEqual(rundir.get_seed_log_file(42).name, "fit_seed_42.log")
        self.assertEqual(rundir.get_mcmc_worker_root_file(3).name, "mcmc_worker_3.root")
        self.assertEqual(rundir.get_mcmc_worker_log_file(3).name, "mcmc_worker_3.log")

    def test_get_run_directory_override(self):
        cfg = FitConfig()
        cfg.project.name = "my_custom_project"
        rundir = get_run_directory(cfg, str(self.temp_dir / "explicit_run"))
        self.assertEqual(rundir.root, (self.temp_dir / "explicit_run").resolve())

        rundir_default = get_run_directory(cfg)
        self.assertEqual(rundir_default.root, Path("runs/my_custom_project").resolve())


if __name__ == "__main__":
    unittest.main()

"""
Tests for RunDirectory management and filesystem layout.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from moments2amplitudes.rundir import RunDirectory, get_run_directory, resolve_init_file_for_tbin
from moments2amplitudes.driver import resolve_mcmc_tree_path
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

    def test_rundir_tbin_templating(self):
        cfg = FitConfig()
        cfg.project.name = "proj"
        cfg.dataset.t_bin = "t010020"
        cfg.multi_start.output_dir = "logs_dbw_{t_bin}"

        rundir = get_run_directory(cfg)
        self.assertEqual(rundir.root.name, "logs_dbw_t010020")

        # Test default runs/ naming when output_dir is default
        cfg.multi_start.output_dir = ""
        rundir_default = get_run_directory(cfg)
        self.assertEqual(rundir_default.root.name, "proj_t010020")

    def test_rundir_auto_nesting(self):
        cfg = FitConfig()
        cfg.project.name = "proj"
        cfg.dataset.t_bin = "t020032"
        cfg.multi_start.output_dir = "logs_dbw"  # No {tbin} template

        rundir = get_run_directory(cfg)
        self.assertEqual(rundir.root.name, "t020032")
        self.assertEqual(rundir.root.parent.name, "logs_dbw")

    def test_resolve_init_file_for_tbin_templated(self):
        res = resolve_init_file_for_tbin("logs/{tbin}/fit/nominal_best.root", "t020032", repo_root=self.temp_dir)
        self.assertEqual(res, self.temp_dir / "logs/t020032/fit/nominal_best.root")

    def test_resolve_init_file_for_tbin_directory(self):
        run_root = self.temp_dir / "logs_stage1"
        tbin_dir = run_root / "t020032" / "fit"
        tbin_dir.mkdir(parents=True)
        target_f = tbin_dir / "nominal_best.root"
        target_f.touch()

        res = resolve_init_file_for_tbin(run_root, "t020032", repo_root=self.temp_dir)
        self.assertEqual(res, target_f)

    def test_resolve_init_file_for_tbin_untemplated_path(self):
        run_root = self.temp_dir / "logs_stage1"
        tbin_dir = run_root / "t032050" / "fit"
        tbin_dir.mkdir(parents=True)
        target_f = tbin_dir / "nominal_best.root"
        target_f.touch()

        # User gave path pointing to a generic nominal_best.root
        user_path = run_root / "fit" / "nominal_best.root"
        res = resolve_init_file_for_tbin(user_path, "t032050", repo_root=self.temp_dir)
        self.assertEqual(res, target_f)

    def test_resolve_mcmc_tree_path(self):
        # 1. Templated path
        p1 = "/data/{tbin}/tree.root"
        self.assertEqual(resolve_mcmc_tree_path(p1, "t010020"), "/data/t010020/tree.root")

        # 2. Path with default t-bin substituted when real file exists
        p2_dir = self.temp_dir / "data" / "t020032"
        p2_dir.mkdir(parents=True)
        p2_file = p2_dir / "ResultsBruMcmc.root"
        p2_file.touch()

        input_path = str(self.temp_dir / "data" / "t010020" / "ResultsBruMcmc.root")
        res2 = resolve_mcmc_tree_path(input_path, "t020032")
        self.assertEqual(res2, str(p2_file))


if __name__ == "__main__":
    unittest.main()

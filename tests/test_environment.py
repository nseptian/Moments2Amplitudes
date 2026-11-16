import os
import unittest
from pathlib import Path
from moments2amplitudes.environment import (
    resolve_rootsys,
    resolve_brufit,
    get_root_binary,
    get_root_subprocess_env,
    is_pixi_or_conda_path,
    DEFAULT_ROOTSYS,
    DEFAULT_BRUFIT,
    get_default_max_workers,
)


class TestEnvironment(unittest.TestCase):
    def test_is_pixi_or_conda_path(self):
        self.assertTrue(is_pixi_or_conda_path("/home/user/repo/.pixi/envs/default"))
        self.assertTrue(is_pixi_or_conda_path("/opt/miniconda3/envs/myenv"))
        self.assertTrue(is_pixi_or_conda_path("/d/grid17/septian/Moments2Amplitudes/.pixi/envs/default/bin/root"))
        self.assertFalse(is_pixi_or_conda_path("/d/grid17/septian/root_v6-36-08_cuda_install"))
        self.assertFalse(is_pixi_or_conda_path("/usr/local/root"))

    def test_resolve_rootsys_defaults_or_env(self):
        rootsys = resolve_rootsys(None)
        self.assertIn("root", rootsys.lower())
        self.assertFalse(is_pixi_or_conda_path(rootsys))

    def test_resolve_brufit_defaults_or_env(self):
        brufit = resolve_brufit(None)
        self.assertIn("brufit", brufit.lower())

    def test_get_root_subprocess_env(self):
        env = get_root_subprocess_env(None)
        self.assertIn("ROOTSYS", env)
        self.assertIn("BRUFIT", env)
        self.assertIn("PATH", env)
        self.assertIn("LD_LIBRARY_PATH", env)
        self.assertIn("ROOT_INCLUDE_PATH", env)
        # Ensure rootsys bin is in PATH
        rootsys = env["ROOTSYS"]
        self.assertTrue(env["PATH"].startswith(f"{rootsys}/bin"))
        # Ensure rootsys lib is in LD_LIBRARY_PATH
        self.assertTrue(env["LD_LIBRARY_PATH"].startswith(f"{rootsys}/lib"))

    def test_get_default_max_workers(self):
        cpus = os.cpu_count() or 2
        expected = max(1, cpus // 2)
        self.assertEqual(get_default_max_workers(), expected)


if __name__ == "__main__":
    unittest.main()

"""
Tests for YAML configuration loading and validation across all configurations.
"""

import unittest
from pathlib import Path
from moments2amplitudes.config import load_config, validate_config
from moments2amplitudes.schema import FitConfig, ModelType, WaveType


class TestConfig(unittest.TestCase):
    def test_load_all_existing_configs(self):
        config_dir = Path("configs")
        config_files = list(config_dir.glob("*.yaml"))
        self.assertGreaterEqual(len(config_files), 4)

        for cfg_file in config_files:
            cfg = load_config(cfg_file)
            self.assertIsInstance(cfg, FitConfig)
            self.assertNotEqual(cfg.project.name, "")
            errors = validate_config(cfg)
            self.assertEqual(len(errors), 0, f"Config {cfg_file.name} had errors: {errors}")

    def test_breit_wigner_config(self):
        cfg = load_config(Path("configs/d_breit_wigner.yaml"))
        self.assertEqual(cfg.project.name, "EtaPi0_DBreitWigner_Stage1")
        d_wave = next(w for w in cfg.waves if w.l == 2)
        self.assertEqual(d_wave.model, ModelType.BREIT_WIGNER)
        self.assertEqual(len(d_wave.resonances), 1)
        self.assertEqual(d_wave.resonances[0].name, "a2_1320")

    def test_two_breit_wigner_config(self):
        cfg = load_config(Path("configs/d_two_breit_wigner.yaml"))
        self.assertEqual(cfg.project.name, "EtaPi0_DTwoBreitWigner")
        d_wave = next(w for w in cfg.waves if w.l == 2)
        self.assertEqual(d_wave.model, ModelType.TWO_BREIT_WIGNER)
        self.assertEqual(len(d_wave.resonances), 2)

    def test_flatte_config(self):
        cfg = load_config(Path("configs/s_flatte_d_bw.yaml"))
        self.assertEqual(cfg.project.name, "EtaPi0_SFlatte_DBreitWigner")
        s_wave = next(w for w in cfg.waves if w.l == 0)
        self.assertEqual(s_wave.model, ModelType.FLATTE)

    def test_piecewise_config(self):
        cfg = load_config(Path("configs/s_piecewise_d_bw.yaml"))
        self.assertEqual(cfg.project.name, "EtaPi0_SPiecewise_DBreitWigner_Stage2")
        s_wave = next(w for w in cfg.waves if w.l == 0)
        self.assertEqual(s_wave.type, WaveType.PIECEWISE)
        self.assertTrue(cfg.plotting.enabled)
        self.assertEqual(cfg.plotting.n_lowest_chi2, 10)
        self.assertEqual(cfg.plotting.moments, ["all"])

    def test_validation_errors(self):
        cfg = FitConfig()
        cfg.dataset.moments_dir = ""
        errors = validate_config(cfg)
        self.assertTrue(any("moments_dir" in e for e in errors))


if __name__ == "__main__":
    unittest.main()

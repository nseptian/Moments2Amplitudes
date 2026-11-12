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

    def test_t_bin_resolution_and_aliases(self):
        # Multiple t-bins from a2CS_D_BreitWigner.yaml
        cfg = load_config(Path("configs/a2CS_D_BreitWigner.yaml"))
        self.assertEqual(len(cfg.dataset.get_t_bins()), 5)
        self.assertEqual(cfg.dataset.get_t_bins()[0], "t010020")
        self.assertTrue(cfg.dataset.get_moments_dir().endswith("t010020"))
        self.assertTrue(cfg.dataset.get_moments_dir("t020032").endswith("t020032"))

        # Do not duplicate if moments_dir already ends with t_bin
        cfg.dataset.moments_dir = "/path/to/data/t010020"
        cfg.dataset.t_bin = "t010020"
        self.assertEqual(cfg.dataset.get_moments_dir(), "/path/to/data/t010020")

        # Absolute t_bin
        cfg.dataset.t_bin = "/custom/path/to/t020032"
        self.assertEqual(cfg.dataset.get_moments_dir(), "/custom/path/to/t020032")

    def test_parse_t_bins_forms(self):
        from moments2amplitudes.config import parse_t_bins_input, STANDARD_TBINS

        self.assertEqual(parse_t_bins_input("all"), STANDARD_TBINS)
        self.assertEqual(parse_t_bins_input(["all"]), STANDARD_TBINS)
        self.assertEqual(parse_t_bins_input([0, 1, 2]), ["t010020", "t020032", "t032050"])
        self.assertEqual(parse_t_bins_input("t010020,t020032"), ["t010020", "t020032"])
        self.assertEqual(parse_t_bins_input(["t010020", "t032050"]), ["t010020", "t032050"])
        self.assertEqual(parse_t_bins_input(3), ["t050075"])

    def test_optional_environment_block(self):
        cfg = load_config(Path("configs/a2CS_D_BreitWigner.yaml"))
        self.assertIsNotNone(cfg.environment)
        self.assertIsNotNone(cfg.environment.brufit_dir)
        self.assertGreater(len(cfg.environment.brufit_dir), 0)


if __name__ == "__main__":
    unittest.main()

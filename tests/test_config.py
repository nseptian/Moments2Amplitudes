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
            if cfg_file.name.startswith("momentfit_"):
                from moments2amplitudes.momentfit.config import load_momentfit_config, validate_momentfit_config
                p_cfg = load_momentfit_config(cfg_file)
                p_errors = validate_momentfit_config(p_cfg, check_files=False)
                self.assertEqual(len(p_errors), 0, f"MomentFit config {cfg_file.name} had errors: {p_errors}")
                continue
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

    def test_barrier_config_parsing(self):
        from moments2amplitudes.config import parse_wave_config
        # Single barrier BW
        wave_dict = {
            "l": 2,
            "name": "D-wave",
            "type": "mass_dependent",
            "model": "barrier_breit_wigner",
            "resonances": [{"name": "a2_1320", "mass": 1.318, "width": 0.107}],
        }
        w = parse_wave_config(wave_dict)
        self.assertEqual(w.model, ModelType.BARRIER_BREIT_WIGNER)
        self.assertTrue(w.barrier)
        self.assertEqual(w.barrier_l, 2)
        self.assertIn("barrier_bw_etapi", w.formula)

        # Two barrier BW
        wave2_dict = {
            "l": 2,
            "name": "D-wave",
            "type": "mass_dependent",
            "model": "barrier_two_breit_wigner",
            "resonances": [
                {"name": "a2_1320", "mass": 1.318, "width": 0.107},
                {"name": "a2_1700", "mass": 1.706, "width": 0.380},
            ],
        }
        w2 = parse_wave_config(wave2_dict)
        self.assertEqual(w2.model, ModelType.BARRIER_TWO_BREIT_WIGNER)
        self.assertTrue(w2.barrier)
        self.assertIn("barrier_bw_etapi", w2.formula)

    def test_n_threads_config_preserved(self):
        cfg = load_config(Path("configs/a2CS_SPD_Kinematic_Barrier.yaml"))
        self.assertEqual(cfg.multi_start.n_threads, 50)
        self.assertEqual(cfg.multi_start.n_tbin_workers, 5)

    def test_flatte_config(self):
        cfg = load_config(Path("configs/s_flatte_d_bw.yaml"))
        self.assertEqual(cfg.project.name, "EtaPi0_SFlatte_DBreitWigner")
        s_wave = next(w for w in cfg.waves if w.l == 0)
        self.assertEqual(s_wave.model, ModelType.FLATTE)
        self.assertEqual(s_wave.resonances[0].g_s, 0.5)
        self.assertIsNone(s_wave.resonances[0].width)

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

    def test_fixed_phase_and_shared_phase_options(self):
        # 1. Test from updated SPiecewise_DBreitWigner configs
        cfg = load_config(Path("configs/a2CS_SPiecewise_DBreitWigner.yaml"))
        d_wave = next(w for w in cfg.waves if w.l == 2)
        res = d_wave.resonances[0]
        self.assertEqual(res.name, "a2_1320")
        self.assertTrue(res.phase_fixed)
        self.assertEqual(res.fixed_phase, 0.0)
        self.assertEqual(res.phase, 0.0)
        self.assertEqual(d_wave.init.option, 2)
        self.assertTrue(bool(d_wave.init.shared_phase))

        # 2. Test shared_phase variations and string keywords
        from moments2amplitudes.config import parse_wave_config

        # reflectivity keyword
        w1 = parse_wave_config({
            "l": 2, "name": "D-wave", "type": "mass_dependent", "shared_phase": "reflectivity",
            "resonances": [{"name": "a2_1320"}]
        })
        self.assertEqual(w1.init.option, 2)

        # L keyword
        w2 = parse_wave_config({
            "l": 2, "name": "D-wave", "type": "mass_dependent", "shared_phase": "l",
            "resonances": [{"name": "a2_1320"}]
        })
        self.assertEqual(w2.init.option, 1)

        # Resonance with fixed phase automatically defaults to option 2 if unspecified
        w3 = parse_wave_config({
            "l": 2, "name": "D-wave", "type": "mass_dependent",
            "resonances": [{"name": "a2_1320", "phase": 0.0, "phase_fixed": True}]
        })
        self.assertEqual(w3.init.option, 2)
        self.assertTrue(w3.resonances[0].phase_fixed)
        self.assertEqual(w3.resonances[0].fixed_phase, 0.0)

        # Wave-level phase
        w4 = parse_wave_config({
            "l": 2, "name": "D-wave", "type": "mass_dependent", "phase": 0.0, "phase_fixed": True,
            "resonances": [{"name": "a2_1320"}, {"name": "a2_1700"}]
        })
        self.assertTrue(w4.resonances[0].phase_fixed)
        self.assertTrue(w4.resonances[1].phase_fixed)
        self.assertEqual(w4.init.option, 2)

        # General phase scalar 0.0 directly sets fixed
        w6 = parse_wave_config({
            "l": 2, "name": "D-wave", "type": "mass_dependent",
            "resonances": [{"name": "a2_1320", "phase": 0.0}]
        })
        self.assertTrue(w6.resonances[0].phase_fixed)

        # 3. Test resonance-level shared_phase options
        w7 = parse_wave_config({
            "l": 2, "name": "D-wave", "type": "mass_dependent",
            "resonances": [
                {"name": "a2_1320", "phase": 0.0, "phase_fixed": True, "shared_phase": True},
                {"name": "a2_1700", "phase": 0.0, "phase_fixed": False, "shared_phase": False}
            ]
        })
        self.assertTrue(w7.resonances[0].shared_phase)
        self.assertFalse(w7.resonances[1].shared_phase)
        self.assertEqual(w7.init.option, 2)

        w8 = parse_wave_config({
            "l": 2, "name": "D-wave", "type": "mass_dependent",
            "resonances": [
                {"name": "a2_1320", "shared_phase": "l"}
            ]
        })
        self.assertEqual(w8.resonances[0].shared_phase, "l")
        self.assertEqual(w8.init.option, 1)
        self.assertEqual(w6.resonances[0].fixed_phase, 0.0)

        # General floating phase
        w7 = parse_wave_config({
            "l": 2, "name": "D-wave", "type": "mass_dependent",
            "resonances": [{"name": "a2_1700", "phase": 1.57, "phase_fixed": False}]
        })
        self.assertFalse(w7.resonances[0].phase_fixed)
        self.assertIsNone(w7.resonances[0].fixed_phase)
        self.assertEqual(w7.resonances[0].phase, 1.57)

    def test_chi2_parsing_included_L(self):
        cfg = load_config(Path("configs/a2CS_SFlattePoly_PBreitWigner_DTwoBreitWigner.yaml"))
        self.assertEqual(cfg.chi2.included_L, [0, 1, 2, 3, 4])

    def test_multi_start_seed_range_config_parsing(self):
        import tempfile

        base_yaml = """
        project:
          name: "TestProject"
        dataset:
          moments_dir: "data/test"
          t_bin: "t010020"
        mass_bins:
          min: 1.0
          max: 1.2
          step: 0.1
        waves:
          - l: 0
            name: "S0"
            type: "piecewise"
        multi_start:
          n_starts: 100
          seed_range: [5001, 10000]
          skip_existing: true
        """
        with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
            f.write(base_yaml)
            f.flush()
            cfg = load_config(f.name)

        self.assertEqual(cfg.multi_start.seed_range, (5001, 10000))
        self.assertTrue(cfg.multi_start.skip_existing)
        self.assertEqual(cfg.multi_start.get_seed_range(), (5001, 10000))
        self.assertEqual(len(cfg.multi_start.get_seeds()), 5000)

        # String range
        yaml_str = base_yaml.replace("[5001, 10000]", '"5001:6000"')
        with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
            f.write(yaml_str)
            f.flush()
            cfg_str = load_config(f.name)
        self.assertEqual(cfg_str.multi_start.get_seed_range(), (5001, 6000))
        self.assertEqual(len(cfg_str.multi_start.get_seeds()), 1000)

        # seed_start & seed_end
        yaml_se = base_yaml.replace("seed_range: [5001, 10000]", "seed_start: 2001\n          seed_end: 2500")
        with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
            f.write(yaml_se)
            f.flush()
            cfg_se = load_config(f.name)
        self.assertEqual(cfg_se.multi_start.get_seed_range(), (2001, 2500))
        self.assertEqual(len(cfg_se.multi_start.get_seeds()), 500)

        # Validation failure on invalid range
        yaml_bad = base_yaml.replace("seed_range: [5001, 10000]", "seed_start: 3000\n          seed_end: 2000")
        with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
            f.write(yaml_bad)
            f.flush()
            cfg_bad = load_config(f.name)
        errors = validate_config(cfg_bad)
        self.assertTrue(any("seed" in e for e in errors))

    def test_cli_overrides_seed_range(self):
        import argparse
        from moments2amplitudes.cli import apply_cli_overrides

        cfg = FitConfig()
        args = argparse.Namespace(
            tbins=None,
            t_bin=None,
            style=None,
            distribution_bins=None,
            chi2_bins=None,
            stride=None,
            n_starts=None,
            trials=2500,
            seed_range="5001:7500",
            seed_start=None,
            seed_end=None,
            skip_existing=True,
        )
        apply_cli_overrides(cfg, args)
        self.assertEqual(cfg.multi_start.n_starts, 2500)
        self.assertEqual(cfg.multi_start.seed_range, "5001:7500")
        self.assertTrue(cfg.multi_start.skip_existing)
        self.assertEqual(cfg.multi_start.get_seed_range(), (5001, 7500))
        self.assertEqual(len(cfg.multi_start.get_seeds()), 2500)

    def test_plotting_reference_amplitude_config_and_cli(self):
        import tempfile
        import argparse
        from moments2amplitudes.cli import apply_cli_overrides

        yaml_content = """
project:
  name: "TestAnchor"
plotting:
  reference_amplitude: "a_2_1"
"""
        with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
            f.write(yaml_content)
            f.flush()
            cfg = load_config(f.name)

        self.assertEqual(cfg.plotting.reference_amplitude, "a_2_1")
        self.assertEqual(cfg.plotting.anchor_amplitude, "a_2_1")

        # Test CLI override
        args = argparse.Namespace(
            reference_amplitude="a_2_2",
            anchor_amplitude=None,
        )
        apply_cli_overrides(cfg, args)
        self.assertEqual(cfg.plotting.reference_amplitude, "a_2_2")
        self.assertEqual(cfg.plotting.anchor_amplitude, "a_2_2")

    def test_gaussian_constraints_parsing_and_validation(self):
        import tempfile

        valid_yaml = """
project:
  name: "TestGaussConstraints"
dataset:
  moments_dir: "data/test"
  t_bin: "t010020"
  mass_bins:
    min: 1.0
    max: 1.2
    step: 0.1
waves:
  - l: 2
    name: "D-wave"
    type: "mass_dependent"
    model: "barrier_two_breit_wigner"
    resonances:
      - name: "a2_1320"
        mass: 1.318
        width: 0.107
        gaussian_constraints:
          mass: [1.3182, 0.0006]
          width:
            mean: 0.107
            sigma: 0.005
      - name: "a2_1700"
        mass: 1.706
        width: 0.380
        gaussian_constraints:
  - l: 1
    name: "P-wave"
    type: "mass_dependent"
    model: "breit_wigner"
    resonances:
      - name: "pi1_1400"
        mass: 1.354
        width: 0.330
        gaussian_constraints:
          MD_shared_pi1_1400_Mass: [1.354, 0.02]
"""
        with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
            f.write(valid_yaml)
            f.flush()
            cfg = load_config(f.name)

        self.assertEqual(validate_config(cfg), [])
        d_wave = next(w for w in cfg.waves if w.l == 2)
        self.assertEqual(d_wave.resonances[0].gaussian_constraints["mass"], [1.3182, 0.0006])
        self.assertEqual(d_wave.resonances[0].gaussian_constraints["width"], [0.107, 0.005])
        self.assertEqual(d_wave.resonances[1].gaussian_constraints, {})

        # Invalid shorthand key, non-positive sigma, boolean, non-dict, mismatched model param
        bad_cases = [
            ("maas: [1.318, 0.005]", "recognized parameter"),
            ("mass: [1.318, 0.0]", "sigma must be > 0"),
            ("mass: [1.318, -0.005]", "sigma must be > 0"),
            ("mass: [1.318, true]", "numeric"),
            ("mass: [1.318]", "must be [mean, sigma]"),
            ("g_KK: [0.4, 0.05]", "Breit-Wigner"),
        ]
        base_bad_yaml = """
project:
  name: "BadGauss"
dataset:
  moments_dir: "data/test"
  mass_bins:
    min: 1.0
    max: 1.2
    step: 0.1
waves:
  - l: 2
    name: "D-wave"
    type: "mass_dependent"
    model: "breit_wigner"
    resonances:
      - name: "a2_1320"
        mass: 1.318
        width: 0.107
        gaussian_constraints:
          {gc_block}
"""
        for gc_block, expected_substr in bad_cases:
            with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
                f.write(base_bad_yaml.format(gc_block=gc_block))
                f.flush()
                bad_cfg = load_config(f.name)
            errs = validate_config(bad_cfg)
            self.assertTrue(any(expected_substr in e for e in errs), f"Expected '{expected_substr}' in {errs}")

        # Non-dict gaussian_constraints
        non_dict_yaml = base_bad_yaml.replace("gaussian_constraints:\n          {gc_block}", "gaussian_constraints: [1.318, 0.005]")
        with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
            f.write(non_dict_yaml)
            f.flush()
            nd_cfg = load_config(f.name)
        errs_nd = validate_config(nd_cfg)
        self.assertTrue(any("mapping" in e for e in errs_nd), f"Expected 'mapping' in {errs_nd}")

    def test_minimizer_reject_at_limit_options(self):
        import tempfile
        yaml_str = """
project:
  name: "TestMinLimit"
dataset:
  moments_dir: "data/test"
  mass_bins:
    min: 1.0
    max: 1.2
    step: 0.1
waves:
  - l: 0
    name: "S0"
    type: "piecewise"
minimizer:
  reject_at_limit: false
  limit_tolerance: 1.0e-4
"""
        with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
            f.write(yaml_str)
            f.flush()
            cfg = load_config(f.name)
        self.assertFalse(cfg.minimizer.reject_at_limit)
        self.assertAlmostEqual(cfg.minimizer.limit_tolerance, 1.0e-4)
        self.assertEqual(validate_config(cfg), [])

        cfg.minimizer.limit_tolerance = 0.0
        self.assertTrue(any("limit_tolerance" in e for e in validate_config(cfg)))

    def test_improved_init_defaults_and_unbounded_piecewise_phase_validation(self):
        import tempfile
        yaml_ok = """
project:
  name: "TestInitRules"
dataset:
  moments_dir: "data/test"
  t_bin: "t010020"
  mass_bins:
    min: 1.0
    max: 1.2
    step: 0.1
waves:
  - l: 0
    name: "S0"
    type: "piecewise"
    init:
      strategy: "random"
      phase:
        strategy: "random"
        min: -3.14159265
        max: 3.14159265
        fixed: false
  - l: 2
    name: "D-wave"
    type: "mass_dependent"
    model: "barrier_breit_wigner"
    resonances:
      - name: "a2_1320"
        mass: 1.318
        width: 0.107
    init:
      strategy: "from_file"
      file_path: "logs_a2CS_dbw/{tbin}/fit/nominal_best.root"
      spread: 0.25
      coupling:
        strategy: "random"
        min: -80.0
        max: 80.0
      poly:
        scale: 250.0
        growth: 3.8
"""
        with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
            f.write(yaml_ok)
            f.flush()
            cfg = load_config(f.name)
        self.assertEqual(validate_config(cfg), [])
        s_wave = next(w for w in cfg.waves if w.l == 0)
        d_wave = next(w for w in cfg.waves if w.l == 2)
        self.assertEqual(s_wave.init.magnitude.max_val, 350.0)
        self.assertEqual(d_wave.init.spread, 0.25)
        self.assertEqual(d_wave.init.coupling.min_val, -80.0)
        self.assertEqual(d_wave.init.coupling.max_val, 80.0)
        self.assertEqual(d_wave.init.poly_scale, 250.0)
        self.assertEqual(d_wave.init.poly_growth, 3.8)

        # Limiting piecewise phase creates hard walls in chi2 space and must fail validation
        yaml_bad_pw_phase_limit = yaml_ok.replace(
            "        fixed: false",
            "        fixed: false\n        limits: [-3.14159, 3.14159]",
        )
        with tempfile.NamedTemporaryFile("w", suffix=".yaml") as f:
            f.write(yaml_bad_pw_phase_limit)
            f.flush()
            cfg_bad = load_config(f.name)
        errs = validate_config(cfg_bad)
        self.assertTrue(any("hard wall" in e and "unbounded" in e for e in errs), f"Expected unbounded phase error, got {errs}")

        # Negative spread invalid
        cfg.waves[1].init.spread = -0.1
        self.assertTrue(any("spread" in e for e in validate_config(cfg)))


if __name__ == "__main__":
    unittest.main()

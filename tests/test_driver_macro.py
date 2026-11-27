"""
Tests for ROOT driver macro generation.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock
from moments2amplitudes.config import load_config
from moments2amplitudes.driver import generate_root_macro, run_mcmc_propagation_tbins
from moments2amplitudes.schema import FitConfig


class TestDriverMacro(unittest.TestCase):
    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_generate_macro_for_all_configs(self):
        configs = [
            "configs/d_breit_wigner.yaml",
            "configs/d_two_breit_wigner.yaml",
            "configs/s_flatte_d_bw.yaml",
            "configs/s_piecewise_d_bw.yaml",
            "configs/a2CS_SFlattePoly_PBreitWigner_DTwoBreitWigner.yaml",
            "configs/a2CS_SFlattePoly_PBreitWigner_DTwoBreitWigner_Barrier.yaml",
            "configs/a2CS_SPD_Kinematic_Barrier.yaml",
        ]

        for cpath in configs:
            cfg = load_config(Path(cpath))
            macro_path = self.temp_dir / f"macro_{cfg.project.name}.C"
            generate_root_macro(cfg, macro_path)
            self.assertTrue(macro_path.is_file())
            content = macro_path.read_text()
            self.assertIn("MassDependentMoments moms", content)
            self.assertIn("MassDependentFitter fitter", content)
            self.assertIn("void RunWorker", content)

    def test_macro_generation_fixed_phase_and_shared_phase(self):
        cfg = load_config(Path("configs/a2CS_SPiecewise_DBreitWigner.yaml"))
        macro_path = self.temp_dir / "macro_test_phase_ref.C"
        generate_root_macro(cfg, macro_path)
        self.assertTrue(macro_path.is_file())
        content = macro_path.read_text()

        # Check AddMassDependentParameter uses option 2
        self.assertIn("paramManager.AddMassDependentParameter({2}", content)
        self.assertIn(", 2);", content)

        # Check fixed phase for a2_1320
        self.assertIn("// Fixed phase for resonance a2_1320 (phases fixed to 0.0)", content)
        self.assertIn('"MD_shareda_a2_1320_globalphase"', content)
        self.assertIn('"MD_sharedb_a2_1320_globalphase"', content)
        self.assertIn("paramManager.FixParameter(pName);", content)

    def test_macro_generation_two_breit_wigner_published_model(self):
        cfg = load_config(Path("configs/a2CS_published_model.yaml"))
        macro_path = self.temp_dir / "macro_test_published_model.C"
        generate_root_macro(cfg, macro_path)
        self.assertTrue(macro_path.is_file())
        content = macro_path.read_text()

        # Check single combined wave registered in massDependentWaves[2]
        self.assertIn('massDepConfig.massDependentWaves[2].push_back({"a2_1320_1700"', content)
        self.assertEqual(content.count('massDepConfig.massDependentWaves[2].push_back('), 1)

        # Check explicit formula registered once for combined wave
        self.assertIn('fitter.SetWaveCustomFormula(2, "a2_1320_1700", eval_formula_2_a2_1320_1700', content)

        # Check AddMassDependentParameter uses option 2
        self.assertIn("paramManager.AddMassDependentParameter({2}, seed, 2);", content)

        # Check Mass1 and Width1 initialization and limits
        self.assertIn('paramManager.SetParameterInitialValue("MD_shared_a2_1320_1700_Mass1", 1.318);', content)
        self.assertIn('paramManager.SetParameterLimits("MD_shared_a2_1320_1700_Mass1", 1.2, 1.45);', content)
        self.assertIn('paramManager.SetParameterInitialValue("MD_shared_a2_1320_1700_Width1", 0.107);', content)
        self.assertIn('paramManager.SetParameterLimits("MD_shared_a2_1320_1700_Width1", 0.05, 0.25);', content)

        # Check Mass2 and Width2 initialization and fixed parameters
        self.assertIn('paramManager.SetParameterInitialValue("MD_shared_a2_1320_1700_Mass2", 1.706);', content)
        self.assertIn('paramManager.FixParameter("MD_shared_a2_1320_1700_Mass2");', content)
        self.assertIn('paramManager.SetParameterInitialValue("MD_shared_a2_1320_1700_Width2", 0.38);', content)
        self.assertIn('paramManager.FixParameter("MD_shared_a2_1320_1700_Width2");', content)

        # Check globalphase1 for resonance 1 (fixed to 0.0 as reference)
        self.assertIn('"MD_shareda_a2_1320_1700_globalphase1"', content)
        self.assertIn('"MD_sharedb_a2_1320_1700_globalphase1"', content)
        self.assertIn('Fixed phase for a2_1320_1700 globalphase1', content)

        # Check globalphase2 for resonance 2 (initial 0.0, floating relative phase)
        self.assertIn('"MD_shareda_a2_1320_1700_globalphase2"', content)
        self.assertIn('"MD_sharedb_a2_1320_1700_globalphase2"', content)
        self.assertIn('Initial phase for a2_1320_1700 globalphase2', content)

        # Check RunWorker has mcmcStride parameter and FitMCMCTree uses C++ adaptive stride
        self.assertIn("int mcmcStride = 1", content)
        self.assertIn("Long64_t totalAvailable = moms.GetMCMCEntries();", content)
        self.assertIn("if (stride <= 0) {", content)
        self.assertIn("fitter.FitMCMCTree(paramManager, moms, outFile, startEntry, nEntriesToFit, stride, outDir);", content)

    def test_macro_generation_barrier_breit_wigner(self):
        cfg = load_config(Path("configs/a2CS_D_BreitWigner_Kinematic.yaml"))
        d_wave = next(w for w in cfg.waves if w.l == 2)
        d_wave.formula = "sqrt(rho_etapi) * barrier_bw_etapi(k, M, width, 2)"
        macro_path = self.temp_dir / "macro_test_barrier_bw.C"
        generate_root_macro(cfg, macro_path)
        self.assertTrue(macro_path.is_file())
        content = macro_path.read_text()
        self.assertIn("m2pw::MassDependentFunction::BarrierBreitWigner(mass, k, M, width, 2)", content)
        self.assertIn("std::sqrt(m2pw::MassDependentFunction::PhaseSpaceFactor(mass, 0.547853, 0.13957))", content)

    def test_macro_generation_spd_kinematic_barrier_flatte_couplings(self):
        cfg = load_config(Path("configs/a2CS_SPD_Kinematic_Barrier.yaml"))
        macro_path = self.temp_dir / "macro_test_spd_couplings.C"
        generate_root_macro(cfg, macro_path)
        self.assertTrue(macro_path.is_file())
        content = macro_path.read_text()

        # Check moments L
        self.assertIn("momentsConfig.includedL = {0, 1, 2, 3, 4};", content)

        # Check Flatte couplings: g_etapi, g_KK, g_etaprimepi
        self.assertIn('paramManager.SetParameterInitialValue("MD_shared_a0_980_Mass", 0.987);', content)
        self.assertIn('paramManager.SetParameterLimits("MD_shared_a0_980_Mass", 0.9, 1.05);', content)
        self.assertIn('paramManager.SetParameterInitialValue("MD_shared_a0_980_g_etapi", 0.4);', content)
        self.assertIn('paramManager.SetParameterLimits("MD_shared_a0_980_g_etapi", 0.1, 1.5);', content)
        self.assertIn('paramManager.SetParameterInitialValue("MD_shared_a0_980_g_KK", -99.0);', content)
        self.assertIn('paramManager.FixParameter("MD_shared_a0_980_g_KK");', content)
        self.assertIn('paramManager.SetParameterInitialValue("MD_shared_a0_980_g_etaprimepi", -99.0);', content)
        self.assertIn('paramManager.FixParameter("MD_shared_a0_980_g_etaprimepi");', content)
        self.assertNotIn('MD_shared_a0_980_Width', content)

    def test_macro_generation_resonance_fixed_zero(self):
        cfg = load_config(Path("configs/a2CS_SFlattePoly_PBreitWigner_DTwoBreitWigner.yaml"))
        macro_path = self.temp_dir / "macro_test_fixed_zero.C"
        generate_root_macro(cfg, macro_path)
        self.assertTrue(macro_path.is_file())
        content = macro_path.read_text()
        self.assertIn("// Fixed zero coupling for mass-dependent wave L=1 (pi1_1600)", content)
        self.assertIn('TString kName = Form("MD_%s_%d_%d_pi1_1600_k", prefix, 1, m);', content)
        self.assertIn("paramManager.SetParameterInitialValue(kName, 0.0);", content)
        self.assertIn("paramManager.FixParameter(kName);", content)

    @patch("moments2amplitudes.driver.ensure_aclic_updated")
    @patch("moments2amplitudes.driver.execute_root_command")
    @patch("subprocess.run")
    def test_auto_stride_different_bin_totals(self, mock_subp, mock_exec, mock_aclic):
        from moments2amplitudes.driver import run_mcmc_propagation
        mock_subp.return_value.returncode = 0
        mock_exec.return_value = (0, "")
        cfg = load_config(Path("configs/a2CS_published_model.yaml"))
        cfg.mcmc_propagation.burn_in = 0

        # Scenario 1: By default auto-stride passes eff_stride=0 to trigger C++ runtime adaptive thinning
        from moments2amplitudes.rundir import get_run_directory
        with tempfile.TemporaryDirectory() as td:
            rundir = get_run_directory(cfg, td)
            rundir.ensure_directories()
            for w in range(4):
                (rundir.mcmc_workers_dir / f"mcmc_worker_{w}.root").touch()
            (rundir.root / "central_best_params.root").touch()
            with patch("moments2amplitudes.driver.run_all_diagnostics_and_plots"):
                run_mcmc_propagation(cfg, run_dir=td, n_workers=4, target_entries=200)
            self.assertEqual(mock_exec.call_count, 4)
            # Verify workers partitioned in draw space [0, 200) across 4 workers: [0, 50), [50, 100), [100, 150), [150, 200)
            calls = mock_exec.call_args_list
            expected_spans = [(0, 0, 50), (1, 50, 100), (2, 100, 150), (3, 150, 200)]
            for (call, (exp_w, exp_start, exp_end)) in zip(calls, expected_spans):
                args_tuple = call[0][1]
                self.assertEqual(args_tuple[0], "mcmc")
                self.assertEqual(args_tuple[1], exp_w)
                self.assertEqual(args_tuple[2], exp_start)
                self.assertEqual(args_tuple[3], exp_end)
                self.assertEqual(args_tuple[6], 0)  # eff_stride = 0 (auto-stride)

        # Scenario 2: Explicit manual stride override takes precedence
        mock_exec.reset_mock()
        with tempfile.TemporaryDirectory() as td:
            rundir = get_run_directory(cfg, td)
            rundir.ensure_directories()
            for w in range(4):
                (rundir.mcmc_workers_dir / f"mcmc_worker_{w}.root").touch()
            (rundir.root / "central_best_params.root").touch()
            with patch("moments2amplitudes.driver.run_all_diagnostics_and_plots"):
                run_mcmc_propagation(cfg, run_dir=td, n_workers=4, target_entries=200, stride=3)
            self.assertEqual(mock_exec.call_count, 4)
            for call in mock_exec.call_args_list:
                args_tuple = call[0][1]
                self.assertEqual(args_tuple[6], 3)  # eff_stride = 3

    @patch("moments2amplitudes.driver.ensure_aclic_updated")
    @patch("moments2amplitudes.driver.run_mcmc_propagation")
    def test_run_mcmc_propagation_tbins(self, mock_run, mock_aclic):
        mock_run.return_value = {"output_file": "results.root", "successful_fits": 10}
        cfg = FitConfig()
        cfg.dataset.t_bins = ["t010020", "t020032"]
        cfg.mcmc_propagation.mcmc_tree_path = "/path/{tbin}/tree.root"
        
        results = run_mcmc_propagation_tbins(cfg, n_workers=2, n_tbin_workers=2, target_entries=100)
        self.assertTrue(results["success"])
        self.assertEqual(len(results["tbins"]), 2)
        self.assertIn("t010020", results["tbins"])
        self.assertIn("t020032", results["tbins"])
        self.assertEqual(mock_run.call_count, 2)

    @patch("moments2amplitudes.cli.run_mcmc_propagation_tbins")
    def test_cli_propagate_command(self, mock_prop_tbins):
        import argparse
        from moments2amplitudes.cli import cmd_propagate
        mock_prop_tbins.return_value = {"success": True, "tbins": {}}
        args = argparse.Namespace(
            config="configs/d_breit_wigner.yaml",
            best_params=None,
            mcmc_tree=None,
            n_workers=4,
            tbin_workers=3,
            target_entries=50,
            output_file=None,
            run_dir=None,
            t_bin=None,
            tbins=["0", "1"],
        )
        res = cmd_propagate(args)
        self.assertEqual(res, 0)
        mock_prop_tbins.assert_called_once()
        _, kwargs = mock_prop_tbins.call_args
        self.assertEqual(kwargs["n_workers"], 4)
        self.assertEqual(kwargs["n_tbin_workers"], 3)
        self.assertEqual(kwargs["target_entries"], 50)
        self.assertEqual(kwargs["tbins"], ["t010020", "t020032"])

    @patch("moments2amplitudes.cli.plot_tbins_parallel")
    def test_cli_plot_only_propagated(self, mock_plot_parallel):
        import argparse
        from moments2amplitudes.cli import cmd_plot
        mock_plot_parallel.return_value = {"success": True, "total_plots": 15}
        args = argparse.Namespace(
            config="configs/d_breit_wigner.yaml",
            run_dir=None,
            formats="png",
            n_lowest_chi2=None,
            unnormalized_moments_file=None,
            moments=None,
            t_bin=None,
            tbins=["t010020", "t020032"],
            n_workers=2,
            only_propagated=True,
            only_distributions=False,
            only_amplitudes=False,
        )
        res = cmd_plot(args)
        self.assertEqual(res, 0)
        mock_plot_parallel.assert_called_once()
        _, kwargs = mock_plot_parallel.call_args
        self.assertTrue(kwargs["only_propagated"])
        self.assertFalse(kwargs["only_distributions"])

    @patch("moments2amplitudes.cli.plot_tbins_parallel")
    def test_cli_plot_only_distributions(self, mock_plot_parallel):
        import argparse
        from moments2amplitudes.cli import cmd_plot
        mock_plot_parallel.return_value = {"success": True, "total_plots": 8}
        args = argparse.Namespace(
            config="configs/d_breit_wigner.yaml",
            run_dir=None,
            formats="png",
            n_lowest_chi2=None,
            unnormalized_moments_file=None,
            moments=None,
            t_bin=None,
            tbins=["t010020", "t020032"],
            n_workers=2,
            only_propagated=False,
            only_distributions=True,
            only_amplitudes=False,
        )
        res = cmd_plot(args)
        self.assertEqual(res, 0)
        mock_plot_parallel.assert_called_once()
        _, kwargs = mock_plot_parallel.call_args
        self.assertTrue(kwargs["only_distributions"])
        self.assertFalse(kwargs["only_propagated"])

    def test_macro_generation_gaussian_constraints_single_and_two_bw(self):
        cfg = load_config(Path("configs/a2CS_SPiecewise_DBreitWigner_Kinematic.yaml"))
        d_wave = next(w for w in cfg.waves if w.l == 2)
        d_wave.resonances[0].gaussian_constraints = {
            "mass": [1.3182, 0.0006],
            "Width": [0.107, 0.005],
        }
        macro_path = self.temp_dir / "macro_test_gauss_single.C"
        generate_root_macro(cfg, macro_path)
        content = macro_path.read_text()
        self.assertIn("#include <stdexcept>", content)
        self.assertIn('if (!paramManager.SetGaussianConstraint("MD_shared_a2_1320_Mass", 1.3182, 0.0006))', content)
        self.assertIn('if (!paramManager.SetGaussianConstraint("MD_shared_a2_1320_Width", 0.107, 0.005))', content)

        # Two-Breit-Wigner: both shorthand keys and single-BW-style MD_shared_<res.name>_Mass should map to wave_name + index
        cfg2 = load_config(Path("configs/a2CS_published_model.yaml"))
        d_wave2 = next(w for w in cfg2.waves if w.l == 2)
        d_wave2.resonances[0].gaussian_constraints = {
            "mass": [1.3182, 0.0006],
            "width": [0.107, 0.005],
        }
        d_wave2.resonances[1].gaussian_constraints = {
            "MD_shared_a2_1700_Mass": [1.706, 0.02],
            "width": [0.38, 0.05],
        }
        macro_path2 = self.temp_dir / "macro_test_gauss_twobw.C"
        generate_root_macro(cfg2, macro_path2)
        content2 = macro_path2.read_text()
        self.assertIn('if (!paramManager.SetGaussianConstraint("MD_shared_a2_1320_1700_Mass1", 1.3182, 0.0006))', content2)
        self.assertIn('if (!paramManager.SetGaussianConstraint("MD_shared_a2_1320_1700_Width1", 0.107, 0.005))', content2)
        self.assertIn('if (!paramManager.SetGaussianConstraint("MD_shared_a2_1320_1700_Mass2", 1.706, 0.02))', content2)
        self.assertIn('if (!paramManager.SetGaussianConstraint("MD_shared_a2_1320_1700_Width2", 0.38, 0.05))', content2)

    def test_macro_generation_boundary_limit_checks_and_multistart_rejection(self):
        from moments2amplitudes.driver import run_multi_start
        cfg = load_config(Path("configs/d_breit_wigner.yaml"))
        macro_path = self.temp_dir / "macro_test_limit_check.C"
        generate_root_macro(cfg, macro_path)
        content = macro_path.read_text()
        self.assertIn("fitter.SetRejectBoundaryLimits(true);", content)
        self.assertIn("fitter.SetBoundaryLimitTolerance(1e-05);", content)
        self.assertIn("const bool atLimit = fitter.HasMinimizerHitLimit();", content)
        self.assertIn('<< " AT_LIMIT=" << (atLimit ? 1 : 0)', content)

        # Verify multi-start rejects seed that terminated at limit (AT_LIMIT=1) even if VALID=1 STATUS=0
        cfg.multi_start.n_starts = 2
        cfg.multi_start.n_threads = 1
        cfg.multi_start.output_dir = str(self.temp_dir / "ms_limit_run")

        def fake_exec(mpath, args, log_path, config, timeout=600):
            seed = args[1]
            out_dir = Path(args[4])
            out_dir.mkdir(parents=True, exist_ok=True)
            out_root = out_dir / f"MassDependentResults_{seed}.root"
            out_root.touch()
            if seed == 1:
                out = f"[RESULT] SEED=1 CHI2=100.0 VALID=0 AT_LIMIT=1 STATUS=0 OUT={out_root}\n"
            else:
                out = f"[RESULT] SEED=2 CHI2=120.0 VALID=0 AT_LIMIT=1 STATUS=0 OUT={out_root}\n"
            Path(log_path).write_text(out)
            return 0, out

        with patch("moments2amplitudes.driver.ensure_aclic_updated"), \
             patch("moments2amplitudes.driver.execute_root_command", side_effect=fake_exec), \
             patch("moments2amplitudes.driver.plot_multistart_landscape"):
            res = run_multi_start(cfg)
            self.assertFalse(res["converged"])

    def test_macro_generation_piecewise_unbounded_phase_and_perturbed_from_file(self):
        yaml_str = """
project:
  name: "TestImprovedInitMacro"
dataset:
  moments_dir: "data/test"
  t_bin: "t010020"
  mass_bins:
    min: 1.06
    max: 1.58
    step: 0.04
chi2:
  included_L: [0, 1, 2, 3, 4]
  max_l: 4
waves:
  - l: 0
    name: "S-wave"
    type: "piecewise"
    init:
      strategy: "from_file"
      file_path: "logs_best/{tbin}/fit/nominal_best.root"
      spread: 0.25
      magnitude:
        strategy: "random"
        min: 0.0
        max: 350.0
      phase:
        strategy: "random"
        min: -3.14159265
        max: 3.14159265
        fixed: false
  - l: 2
    name: "D-wave"
    type: "mass_dependent"
    model: "barrier_two_breit_wigner"
    resonances:
      - name: "a2_1320"
        mass: 1.318
        width: 0.107
        mass_limits: [1.20, 1.45]
        width_limits: [0.05, 0.25]
        phase: 0.0
        phase_fixed: true
        shared_phase: true
      - name: "a2_1700"
        mass: 1.706
        width: 0.380
        mass_limits: [1.50, 1.85]
        width_limits: [0.10, 0.40]
        phase: 0.0
        phase_fixed: false
        shared_phase: true
    init:
      strategy: "from_file"
      file_path: "logs_best/{tbin}/fit/nominal_best.root"
      spread: 0.25
      coupling:
        strategy: "random"
        min: -80.0
        max: 80.0
      phase:
        strategy: "random"
        fixed: false
"""
        cfg_file = self.temp_dir / "improved_init.yaml"
        cfg_file.write_text(yaml_str)
        cfg = load_config(cfg_file)
        macro_path = self.temp_dir / "macro_improved_init.C"
        generate_root_macro(cfg, macro_path)
        content = macro_path.read_text()

        # Piecewise wave uses RandomFullPhase [0, 350] with seed-aware file init and ZERO phase limits
        self.assertIn("WaveInitRule::Strategy::RandomFullPhase, 0.0, 350.0, 0.0, false, -3.14159265, 3.14159265, 0.25, 10.0,", content)
        self.assertIn('paramManager.AddMassIndependentParameter({0}, "', content)
        self.assertIn('logs_best/t010020/fit/nominal_best.root", seed);', content)
        self.assertNotIn('SetParameterLimits("MI_', content)
        self.assertNotIn("SetParameterLimits(aPhaseName", content)
        self.assertNotIn("SetParameterLimits(bPhaseName", content)

        # Mass-dependent wave emits MDInitRule + WasLoadedFromFile guards + keeps reference phase fixed
        self.assertIn("MassDependentFitter::MDInitRule mdRule_2;", content)
        self.assertIn("mdRule_2.couplingMin = -80.0;", content)
        self.assertIn("mdRule_2.couplingMax = 80.0;", content)
        self.assertIn("mdRule_2.polyScale = 250.0;", content)
        self.assertIn("mdRule_2.polyGrowth = 3.8;", content)
        self.assertIn("mdRule_2.fileRelSpread = 0.25;", content)
        self.assertIn("mdRule_2.randomizeLoadedCouplings = true;", content)
        self.assertIn("mdRule_2.randomizeFloatingPhases = true;", content)
        self.assertIn("paramManager.SetMDInitRule(2, mdRule_2);", content)
        self.assertIn('paramManager.AddMassDependentParameter({2}, "', content)
        self.assertIn('logs_best/t010020/fit/nominal_best.root", 2, seed);', content)
        self.assertIn('paramManager.WasLoadedFromFile("MD_shared_a2_1320_1700_Mass1")', content)
        self.assertIn("Fixed phase for a2_1320_1700 globalphase1", content)
        self.assertNotIn("Initial phase for a2_1320_1700 globalphase2", content)

    def test_macro_generation_fit_a21700_peak_const_lineshape_and_constraints(self):
        cfg = load_config(Path("configs/SPiecewise_DTwoBreitWigner_FitA21700.yaml"))
        macro_path = self.temp_dir / "macro_fit_a21700.C"
        generate_root_macro(cfg, macro_path)
        content = macro_path.read_text()

        self.assertIn("momentsConfig.includedL = {0, 2, 4};", content)
        self.assertIn("const double k2 = (p.size() > 3) ? p[3] : 0.0;", content)
        self.assertIn("WaveInitRule::Strategy::FixedZero", content)
        self.assertIn("WaveInitRule::Strategy::RandomFullPhase, 0.0, 350.0", content)
        self.assertIn("m2pw::MassDependentFunction::BarrierBreitWigner(mass, k1, Mass1, Width1, 2)", content)
        self.assertIn("m2pw::MassDependentFunction::PeakNormalizedConstantWidthBarrierBreitWigner(mass, k2, Mass2, Width2, 2)", content)
        self.assertIn("mdRule_2.couplingMin = -120.0;", content)
        self.assertIn("mdRule_2.couplingMax = 120.0;", content)
        self.assertIn('paramManager.SetGaussianConstraint("MD_shared_a2_1320_1700_Mass2", 1.706, 0.014)', content)
        self.assertIn('paramManager.SetGaussianConstraint("MD_shared_a2_1320_1700_Width2", 0.38, 0.055)', content)


if __name__ == "__main__":
    unittest.main()

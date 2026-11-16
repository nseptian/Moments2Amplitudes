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

        # Check RunWorker has mcmcStride parameter and FitMCMCTree uses it
        self.assertIn("int mcmcStride = 1", content)
        self.assertIn("fitter.FitMCMCTree(paramManager, moms, outFile, mcmcStart, nEntriesToFit, mcmcStride, outDir);", content)

    @patch("moments2amplitudes.driver.get_min_available_mcmc_entries")
    @patch("moments2amplitudes.driver.ensure_aclic_updated")
    @patch("moments2amplitudes.driver.execute_root_command")
    @patch("subprocess.run")
    def test_auto_stride_different_bin_totals(self, mock_subp, mock_exec, mock_aclic, mock_min_entries):
        from moments2amplitudes.driver import run_mcmc_propagation
        mock_subp.return_value.returncode = 0
        mock_exec.return_value = (0, "")
        cfg = load_config(Path("configs/a2CS_published_model.yaml"))
        cfg.mcmc_propagation.burn_in = 0

        # Scenario 1: Bin with 1582 entries, target 200 draws -> stride = 1582 // 200 = 7
        from moments2amplitudes.rundir import get_run_directory
        mock_min_entries.return_value = 1582
        with tempfile.TemporaryDirectory() as td:
            rundir = get_run_directory(cfg, td)
            rundir.ensure_directories()
            for w in range(4):
                (rundir.mcmc_workers_dir / f"mcmc_worker_{w}.root").touch()
            (rundir.root / "central_best_params.root").touch()
            with patch("moments2amplitudes.driver.run_all_diagnostics_and_plots"):
                run_mcmc_propagation(cfg, run_dir=td, n_workers=4, max_entries=200, auto_stride=True)
            for call in mock_exec.call_args_list:
                args_tuple = call[0][1]
                self.assertEqual(args_tuple[6], 7)  # eff_stride = 7

        # Scenario 2: Bin with 2500 entries, target 200 draws -> stride = 2500 // 200 = 12
        mock_exec.reset_mock()
        mock_min_entries.return_value = 2500
        with tempfile.TemporaryDirectory() as td:
            rundir = get_run_directory(cfg, td)
            rundir.ensure_directories()
            for w in range(4):
                (rundir.mcmc_workers_dir / f"mcmc_worker_{w}.root").touch()
            (rundir.root / "central_best_params.root").touch()
            with patch("moments2amplitudes.driver.run_all_diagnostics_and_plots"):
                run_mcmc_propagation(cfg, run_dir=td, n_workers=4, max_entries=200, auto_stride=True)
            for call in mock_exec.call_args_list:
                args_tuple = call[0][1]
                self.assertEqual(args_tuple[6], 12)  # eff_stride = 12

    @patch("moments2amplitudes.driver.ensure_aclic_updated")
    @patch("moments2amplitudes.driver.run_mcmc_propagation")
    def test_run_mcmc_propagation_tbins(self, mock_run, mock_aclic):
        mock_run.return_value = {"output_file": "results.root", "successful_fits": 10}
        cfg = FitConfig()
        cfg.dataset.t_bins = ["t010020", "t020032"]
        cfg.mcmc_propagation.mcmc_tree_path = "/path/{tbin}/tree.root"
        
        results = run_mcmc_propagation_tbins(cfg, n_workers=2, n_tbin_workers=2, max_entries=100)
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
            max_entries=50,
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
        self.assertEqual(kwargs["max_entries"], 50)
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


if __name__ == "__main__":
    unittest.main()

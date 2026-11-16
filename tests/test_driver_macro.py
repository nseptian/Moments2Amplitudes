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

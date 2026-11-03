"""
Tests for ROOT driver macro generation.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from moments2amplitudes.config import load_config
from moments2amplitudes.driver import generate_root_macro


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


if __name__ == "__main__":
    unittest.main()

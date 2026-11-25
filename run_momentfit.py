#!/usr/bin/env python3
"""
Top-level executable launcher for the MomentFit Analysis Workflow.

Stages executed in order:
  1. Linear algebra moments calculation (Moments.Main)
  2. Export moments value and covariance (SPD regularization for BruFit)
  3. MomentFit with BruFit MCMC using fixed covariance from linear algebra

Can be invoked directly:
  ./run_momentfit.py --config configs/momentfit_nominal.yaml
or via pixi:
  pixi run momentfit --config configs/momentfit_nominal.yaml
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def _is_inside_pixi_or_has_root() -> bool:
    """Checks if we are already inside a working ROOT / Moments.Main environment."""
    try:
        import ROOT
        # Check that TMatrixDSym is available
        _ = ROOT.TMatrixDSym
        return True
    except (ImportError, Exception):
        return False


def main():
    repo_root = Path(__file__).resolve().parent
    python_dir = repo_root / "python"
    if str(python_dir) not in sys.path:
        sys.path.insert(0, str(python_dir))

    # Detect if we need to dispatch via pixi
    if not _is_inside_pixi_or_has_root():
        pixi_bin = shutil.which("pixi") or os.path.expanduser("~/.pixi/bin/pixi")
        # Manifest candidate: project pixi.toml first, then fallbacks
        manifest_candidates = [
            str(repo_root / "pixi.toml"),
            "/d/grid17/septian/Moments-dev/Moments.Main/pixi.toml",
            "/d/grid17/septian/Moments-dev/Moments.Nizar/pixi.toml",
        ]
        manifest_path = next((m for m in manifest_candidates if os.path.isfile(m)), None)

        if os.path.isfile(pixi_bin) and manifest_path:
            print(f"[*] Dispatching execution into pixi environment ({manifest_path})...")
            cmd = [pixi_bin, "run", "--manifest-path", manifest_path, "python", str(Path(__file__).resolve())] + sys.argv[1:]
            ret = subprocess.run(cmd)
            sys.exit(ret.returncode)
        else:
            print("[!] Warning: Not running in pixi environment and ROOT could not be loaded directly.")

    from moments2amplitudes.momentfit.runner import main as momentfit_main
    momentfit_main()


if __name__ == "__main__":
    main()

"""
Run Directory Management for Moments2Amplitudes.

Provides standardized directory structure for fit workflows:
runs/<project_name>/
├── config.yaml          # Snapshot of configuration
├── worker.C             # Generated C++ ROOT worker macro
├── fit/                 # Multi-start fit results
│   ├── starts/          # Individual seed results & logs
│   ├── nominal_best.root# Best global minimum fit
│   └── all_starts.root  # Merged multi-start results
├── mcmc/                # MCMC uncertainty propagation
│   ├── workers/         # Parallel worker chunks & logs
│   └── mcmc_propagated.root # Merged posterior samples
├── plots/               # Publication figures (PDF & PNG)
└── reports/             # Summary tables (txt, md)
"""

from __future__ import annotations
import shutil
from pathlib import Path
from typing import Optional, List, Dict, Any

from .schema import FitConfig


class RunDirectory:
    """Manages the filesystem hierarchy for a specific fit project run."""

    def __init__(self, run_dir: Path | str):
        self.root = Path(run_dir).resolve()
        
        # Subdirectories
        self.fit_dir = self.root / "fit"
        self.fit_starts_dir = self.fit_dir / "starts"
        self.mcmc_dir = self.root / "mcmc"
        self.mcmc_workers_dir = self.mcmc_dir / "workers"
        self.plots_dir = self.root / "plots"
        self.moments_plots_dir = self.plots_dir / "moments"
        self.reports_dir = self.root / "reports"

        # Key file paths
        self.config_snapshot = self.root / "config.yaml"
        self.worker_macro = self.root / "worker.C"
        self.nominal_best_root = self.fit_dir / "nominal_best.root"
        self.all_starts_root = self.fit_dir / "all_starts.root"
        self.mcmc_propagated_root = self.mcmc_dir / "mcmc_propagated.root"
        self.summary_txt = self.reports_dir / "fit_summary.txt"
        self.summary_md = self.reports_dir / "fit_summary.md"

    def ensure_directories(self) -> RunDirectory:
        """Create all required directory paths if they do not exist."""
        self.root.mkdir(parents=True, exist_ok=True)
        self.fit_starts_dir.mkdir(parents=True, exist_ok=True)
        self.mcmc_workers_dir.mkdir(parents=True, exist_ok=True)
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.moments_plots_dir.mkdir(parents=True, exist_ok=True)
        self.reports_dir.mkdir(parents=True, exist_ok=True)
        return self

    def save_config_snapshot(self, config_source_path: Path | str) -> Path:
        """Save a snapshot of the input YAML configuration for reproducibility."""
        src = Path(config_source_path).resolve()
        if src.is_file():
            self.ensure_directories()
            shutil.copyfile(src, self.config_snapshot)
        return self.config_snapshot

    def save_nominal_best(self, best_seed_file: Path | str) -> Path:
        """Copy the best seed's ROOT output file to nominal_best.root."""
        src = Path(best_seed_file).resolve()
        if src.is_file():
            shutil.copyfile(src, self.nominal_best_root)
        return self.nominal_best_root

    def get_seed_root_file(self, seed: int) -> Path:
        return self.fit_starts_dir / f"fit_seed_{seed}.root"

    def get_seed_log_file(self, seed: int) -> Path:
        return self.fit_starts_dir / f"fit_seed_{seed}.log"

    def get_mcmc_worker_root_file(self, worker_id: int) -> Path:
        return self.mcmc_workers_dir / f"mcmc_worker_{worker_id}.root"

    def get_mcmc_worker_log_file(self, worker_id: int) -> Path:
        return self.mcmc_workers_dir / f"mcmc_worker_{worker_id}.log"


def get_run_directory(config: FitConfig, run_dir_override: Optional[str] = None) -> RunDirectory:
    """
    Resolve the RunDirectory for a configuration.
    If run_dir_override is provided, uses it.
    Otherwise defaults to runs/<project_name>.
    """
    if run_dir_override:
        target_dir = Path(run_dir_override)
    elif config.multi_start.output_dir and config.multi_start.output_dir not in ("logs", "run/central_fits", ""):
        target_dir = Path(config.multi_start.output_dir)
    else:
        proj_name = config.project.name or "fit_run"
        target_dir = Path("runs") / proj_name

    rundir = RunDirectory(target_dir)
    rundir.ensure_directories()
    return rundir

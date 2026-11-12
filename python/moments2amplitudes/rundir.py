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
    Supports {t_bin} and {tbin} templating in output directory paths,
    and automatically nests under the t-bin subdirectory when executing for a specific t-bin.
    """
    t_tag = config.dataset.t_bin or ""
    if run_dir_override:
        has_template = "{t_bin}" in run_dir_override or "{tbin}" in run_dir_override
        raw = run_dir_override.replace("{t_bin}", t_tag).replace("{tbin}", t_tag)
        target_dir = Path(raw)
        if t_tag and not has_template and target_dir.name != t_tag:
            target_dir = target_dir / t_tag
    elif config.multi_start.output_dir and config.multi_start.output_dir not in ("logs", "run/central_fits", ""):
        raw_out = config.multi_start.output_dir
        has_template = "{t_bin}" in raw_out or "{tbin}" in raw_out
        raw = raw_out.replace("{t_bin}", t_tag).replace("{tbin}", t_tag)
        target_dir = Path(raw)
        if t_tag and not has_template and target_dir.name != t_tag:
            target_dir = target_dir / t_tag
    else:
        proj_name = config.project.name or "fit_run"
        if t_tag:
            target_dir = Path("runs") / f"{proj_name}_{t_tag}"
        else:
            target_dir = Path("runs") / proj_name

    rundir = RunDirectory(target_dir)
    rundir.ensure_directories()
    return rundir


def resolve_init_file_for_tbin(
    file_path: str | Path,
    t_bin: Optional[str] = None,
    repo_root: Optional[Path] = None,
) -> Path:
    """
    Resolve initial parameters ROOT file path for a given t-bin.
    Supports:
    1. Direct templating: e.g. logs_a2CS_dbw/{tbin}/fit/nominal_best.root
    2. Multi-tbin directory: e.g. logs_a2CS_dbw, discovering logs_a2CS_dbw/<tbin>/fit/nominal_best.root
    3. Untemplated path: e.g. logs_a2CS_dbw/fit/nominal_best.root, checking if <dir>/<tbin>/fit/nominal_best.root exists
    4. Exact file fallback if single fit or global file.
    """
    if not file_path:
        raise ValueError("file_path cannot be empty.")

    root_dir = repo_root or Path(__file__).resolve().parent.parent.parent
    raw_str = str(file_path).strip()
    t_clean = str(t_bin).strip() if t_bin else ""

    # 1. Direct template replacement if {t_bin} or {tbin} or {t-bin} present
    if t_clean and ("{t_bin}" in raw_str or "{tbin}" in raw_str or "{t-bin}" in raw_str):
        formatted = raw_str.replace("{t_bin}", t_clean).replace("{tbin}", t_clean).replace("{t-bin}", t_clean)
        cand = Path(formatted)
        if not cand.is_absolute():
            cand = root_dir / cand
        return cand

    cand_path = Path(raw_str)
    if not cand_path.is_absolute():
        cand_path = root_dir / cand_path

    if not t_clean:
        return cand_path

    # 2. If cand_path is a directory (e.g. "logs_a2CS_dbw")
    if cand_path.is_dir():
        for sub_p in [
            cand_path / t_clean / "fit" / "nominal_best.root",
            cand_path / t_clean / "central_best_params.root",
            cand_path / t_clean / "fit" / "starts" / "MassDependentResults_1.root",
            cand_path / t_clean / "nominal_best.root",
            cand_path / "fit" / "nominal_best.root",
            cand_path / "central_best_params.root",
        ]:
            if sub_p.is_file():
                return sub_p

    # 3. If cand_path is a file path (e.g. "logs_a2CS_dbw/fit/nominal_best.root")
    # Check if a t-bin specific sibling / descendant exists
    parts = list(cand_path.parts)
    if t_clean not in parts:
        for idx in range(len(parts) - 1, 0, -1):
            cand_parts = parts[:idx] + [t_clean] + parts[idx:]
            t_cand = Path(*cand_parts)
            if t_cand.is_file():
                return t_cand

        for base_p in [cand_path.parent, cand_path.parent.parent]:
            for alt_name in [
                base_p / t_clean / "fit" / "nominal_best.root",
                base_p / t_clean / "central_best_params.root",
                base_p / t_clean / cand_path.name,
            ]:
                if alt_name.is_file():
                    return alt_name

    # 4. If exact file exists as given, return it
    if cand_path.is_file():
        return cand_path

    # Fallback candidate for clear error reporting
    if cand_path.is_dir():
        return cand_path / t_clean / "fit" / "nominal_best.root"
    return cand_path

"""
Fit Summary and Diagnostics Table Generator for Moments2Amplitudes.

Extracts fitted parameters, goodness-of-fit metrics, and MCMC credible
intervals, generating terminal ASCII tables and saved markdown reports.
"""

from __future__ import annotations
import os
import math
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple

from .schema import FitConfig, WaveType
from .rundir import RunDirectory
from .environment import is_pixi_or_conda_path


def format_parameter_label(raw_name: str) -> str:
    """Transform internal ROOT parameter names into human-readable physics notation."""
    if raw_name.startswith("MD_"):
        parts = raw_name[3:].split("_")
        if "Mass" in parts or "Width" in parts:
            param = parts[-1]
            res = "_".join(parts[:-1])
            return f"{param} ({res})"
        if len(parts) >= 3 and parts[0] in ("a", "b", "aphi", "bphi"):
            param = parts[-1]
            wave = f"{parts[0]}_{parts[1]}_{parts[2]}"
            res = "_".join(parts[3:-1]) if len(parts) > 4 else parts[3]
            return f"{param} [{wave}] ({res})"
        return raw_name[3:]
    elif raw_name.startswith("MI_"):
        parts = raw_name[3:].split("_")
        if len(parts) >= 4:
            mass = parts[0]
            wave = "_".join(parts[1:])
            return f"{wave} @ {mass} GeV"
        return raw_name[3:]
    return raw_name


class FitSummary:
    """Container for fit quality metrics and parameter values."""

    def __init__(self, config: FitConfig, rundir: RunDirectory):
        self.config = config
        self.rundir = rundir
        self.chi2: float = 0.0
        self.ndf: int = 1
        self.is_valid: bool = False
        self.is_valid_error: bool = False
        self.status: int = -1
        self.best_seed: int = -1
        self.total_starts: int = 0
        self.converged_starts: int = 0
        self.plateau_count: int = 0
        self.elapsed_seconds: float = 0.0
        
        # Parameters: list of dicts {name, label, value, error, mcmc_median, mcmc_lo, mcmc_hi, status}
        self.parameters: List[Dict[str, Any]] = []

    def _should_display_parameter(self, param_name: str) -> bool:
        """Filter out polynomial parameters when poly_order is set to -1."""
        is_poly = (
            param_name.endswith("m_threshold") or
            param_name.endswith("m_expansion") or
            "_re_" in param_name or
            "_im_" in param_name
        )
        if not is_poly:
            return True

        for w in self.config.waves:
            if w.type == WaveType.MASS_DEPENDENT and w.poly_order < 0:
                for res in w.resonances:
                    if res.name in param_name:
                        return False
                if len([mw for mw in self.config.waves if mw.type == WaveType.MASS_DEPENDENT]) == 1:
                    return False
        return True

    def load_from_root(
        self,
        root_file: Path | str,
        mcmc_file: Optional[Path | str] = None,
        cache_json: Optional[Path | str] = None,
    ) -> bool:
        """Load nominal best-fit parameters and optional MCMC posterior intervals."""
        root_path = Path(root_file)
        if not root_path.is_file():
            return False

        loaded = False
        # Attempt uproot
        try:
            import uproot
            loaded = self._load_with_uproot(root_path, mcmc_file)
        except (ImportError, Exception):
            pass

        # If MCMC intervals not populated yet, check json cache
        has_mcmc = any(p.get("mcmc_median") is not None for p in self.parameters)
        if not has_mcmc and (mcmc_file or cache_json or self.rundir):
            self.load_mcmc_from_cache(cache_json or mcmc_file or root_path)

        return loaded

    def load_mcmc_from_cache(self, cache_or_path: Optional[Path | str] = None) -> bool:
        """Populate MCMC posterior intervals from cache JSON if available."""
        import json
        import numpy as np

        cand_paths: List[Path] = []
        if cache_or_path:
            p = Path(cache_or_path)
            if p.is_file() and p.suffix == ".json":
                cand_paths.append(p)
            elif p.is_dir():
                cand_paths.append(p / ".cache" / "amplitudes_data.json")
                cand_paths.append(p / "plots" / ".cache" / "amplitudes_data.json")
            elif p.is_file():
                cand_paths.append(p.parent.parent / "plots" / ".cache" / "amplitudes_data.json")
                cand_paths.append(p.parent / ".cache" / "amplitudes_data.json")

        if self.rundir and getattr(self.rundir, "plots_dir", None):
            cand_paths.append(self.rundir.plots_dir / ".cache" / "amplitudes_data.json")

        for cp in cand_paths:
            if cp.is_file() and cp.stat().st_size > 10:
                try:
                    with open(cp, "r") as f:
                        jdata = json.load(f)
                    raw_params = jdata.get("mcmc_params", {})
                    if raw_params:
                        for p in self.parameters:
                            k = p["name"]
                            if k in raw_params:
                                arr = np.array(raw_params[k], dtype=float)
                                if len(arr) > 0 and np.std(arr) > 1e-12:
                                    med = float(np.median(arr))
                                    p16, p84 = float(np.percentile(arr, 16)), float(np.percentile(arr, 84))
                                    p["mcmc_median"] = med
                                    p["mcmc_err_lo"] = med - p16
                                    p["mcmc_err_hi"] = p84 - med
                        return True
                except Exception:
                    pass
        return False

    def load_from_log(
        self,
        log_path: Path | str,
        mcmc_file: Optional[Path | str] = None,
        cache_json: Optional[Path | str] = None,
    ) -> bool:
        """Parse fit parameters and errors directly from the driver log file."""
        log_file = Path(log_path)
        if not log_file.is_file():
            return False

        param_dict: Dict[str, Dict[str, float]] = {}
        with open(log_file, "r") as f:
            for line in f:
                if "[FIT][RESULT]" in line and "param=" in line:
                    parts = line.strip().split()
                    p_name = ""
                    val = None
                    err = None
                    for token in parts:
                        if token.startswith("param="):
                            p_name = token.split("=")[1]
                        elif token.startswith("value="):
                            try:
                                val = float(token.split("=")[1])
                            except ValueError:
                                pass
                        elif token.startswith("error="):
                            try:
                                err = float(token.split("=")[1])
                            except ValueError:
                                pass
                    if p_name:
                        if p_name not in param_dict:
                            param_dict[p_name] = {"value": 0.0, "error": 0.0}
                        if val is not None:
                            param_dict[p_name]["value"] = val
                        if err is not None:
                            param_dict[p_name]["error"] = err
                elif "[FIT][MIN]" in line:
                    for token in line.strip().split():
                        if token.startswith("is_valid="):
                            self.is_valid = (token.split("=")[1] == "1")
                        elif token.startswith("is_valid_error="):
                            self.is_valid_error = (token.split("=")[1] == "1")
                        elif token.startswith("status="):
                            try:
                                self.status = int(token.split("=")[1])
                            except ValueError:
                                pass

        if not param_dict:
            return False

        self.parameters = []
        for name, data in param_dict.items():
            if not self._should_display_parameter(name):
                continue
            val = data["value"]
            err = data["error"]
            self.parameters.append({
                "name": name,
                "label": format_parameter_label(name),
                "value": val,
                "error": err,
                "mcmc_median": None,
                "mcmc_err_lo": None,
                "mcmc_err_hi": None,
                "status": "Fixed" if err == 0.0 else "Free",
            })

        if mcmc_file or cache_json or self.rundir:
            self.load_mcmc_from_cache(cache_json or mcmc_file)

        return True

    def _load_with_uproot(self, root_path: Path, mcmc_file: Optional[Path | str]) -> bool:
        import uproot
        import numpy as np

        with uproot.open(root_path) as f:
            if "chi2" in f:
                c_tree = f["chi2"]
                self.chi2 = float(c_tree["chi2"].array()[0])
                self.is_valid = bool(c_tree["isValid"].array()[0])
                self.is_valid_error = bool(c_tree["isValidError"].array()[0])
                self.status = int(c_tree["Status"].array()[0])

            if "all_params" in f:
                p_tree = f["all_params"]
                keys = [k for k in p_tree.keys() if not k.endswith("_error")]
                for k in keys:
                    if not self._should_display_parameter(k):
                        continue
                    val = float(p_tree[k].array()[0])
                    err = float(p_tree[f"{k}_error"].array()[0]) if f"{k}_error" in p_tree else 0.0
                    self.parameters.append({
                        "name": k,
                        "label": format_parameter_label(k),
                        "value": val,
                        "error": err,
                        "mcmc_median": None,
                        "mcmc_err_lo": None,
                        "mcmc_err_hi": None,
                        "status": "Fixed" if err == 0.0 else "Free",
                    })

        # MCMC posteriors
        if mcmc_file and Path(mcmc_file).is_file():
            with uproot.open(mcmc_file) as mf:
                if "all_params" in mf or "mcmc_params" in mf:
                    tree_name = "all_params" if "all_params" in mf else "mcmc_params"
                    mtree = mf[tree_name]
                    mask = None
                    if "Status" in mtree and "isValid" in mtree:
                        mask = (np.array(mtree["Status"].array()) == 0) & (np.array(mtree["isValid"].array()) == 1)
                    elif "Status" in mtree:
                        mask = (np.array(mtree["Status"].array()) == 0)
                    elif "isValid" in mtree:
                        mask = (np.array(mtree["isValid"].array()) == 1)

                    for p in self.parameters:
                        k = p["name"]
                        if k in mtree:
                            arr = np.array(mtree[k].array())
                            if mask is not None and len(mask) == len(arr):
                                arr = arr[mask]
                            if len(arr) > 0:
                                med = float(np.median(arr))
                                p16, p84 = float(np.percentile(arr, 16)), float(np.percentile(arr, 84))
                                p["mcmc_median"] = med
                                p["mcmc_err_lo"] = med - p16
                                p["mcmc_err_hi"] = p84 - med

        return True

    def _load_with_pyroot(self, root_path: Path, mcmc_file: Optional[Path | str]) -> bool:
        import ROOT
        f = ROOT.TFile.Open(str(root_path), "READ")
        if not f or f.IsZombie():
            return False

        c_tree = f.Get("chi2")
        if c_tree and c_tree.GetEntries() > 0:
            c_tree.GetEntry(0)
            self.chi2 = float(getattr(c_tree, "chi2", 0.0))
            self.is_valid = bool(getattr(c_tree, "isValid", False))
            self.is_valid_error = bool(getattr(c_tree, "isValidError", False))
            self.status = int(getattr(c_tree, "Status", -1))

        p_tree = f.Get("all_params")
        if p_tree and p_tree.GetEntries() > 0:
            p_tree.GetEntry(0)
            branches = p_tree.GetListOfBranches()
            for b in branches:
                bname = b.GetName()
                if not bname.endswith("_error"):
                    if not self._should_display_parameter(bname):
                        continue
                    val = float(getattr(p_tree, bname, 0.0))
                    err = float(getattr(p_tree, f"{bname}_error", 0.0))
                    self.parameters.append({
                        "name": bname,
                        "label": format_parameter_label(bname),
                        "value": val,
                        "error": err,
                        "mcmc_median": None,
                        "mcmc_err_lo": None,
                        "mcmc_err_hi": None,
                        "status": "Fixed" if err == 0.0 else "Free",
                    })

        f.Close()
        return True

    def generate_ascii_table(self) -> str:
        """Produce an ASCII box-drawing summary table."""
        lines = []
        div = "=" * 80
        lines.append(div)
        lines.append(f" Moments2Amplitudes Fit Summary Report")
        lines.append(f" Project:      {self.config.project.name}")
        lines.append(f" Minimum Chi2: {self.chi2:.4f}")
        status_str = "SUCCESS (0)" if self.status == 0 and self.is_valid else f"STATUS {self.status}"
        lines.append(f" Status:       {status_str} | Valid Errors: {int(self.is_valid_error)}")
        lines.append(div)
        lines.append(f"{'Parameter':<32} {'Nominal':>12} {'Minuit Err':>12} {'MCMC (-/+)':>16}")
        lines.append("-" * 80)

        for p in self.parameters:
            lbl = p["label"][:31]
            val_s = f"{p['value']:12.5f}"
            err_s = f"± {p['error']:<9.5f}" if p["error"] > 0.0 else "fixed"
            if p["mcmc_median"] is not None:
                mcmc_s = f"-{p['mcmc_err_lo']:.4f}/+{p['mcmc_err_hi']:.4f}"
            else:
                mcmc_s = "-"
            lines.append(f"{lbl:<32} {val_s} {err_s:>12} {mcmc_s:>16}")

        lines.append(div)
        return "\n".join(lines)

    def generate_markdown_report(self) -> str:
        """Produce a GitHub-flavored Markdown summary report."""
        lines = [
            f"# Fit Summary: {self.config.project.name}",
            "",
            f"- **Minimum \\(\\chi^2\\)**: `{self.chi2:.4f}`",
            f"- **Minuit Status**: `{'Converged (0)' if self.status == 0 else f'Failed ({self.status})'}`",
            f"- **Valid Covariance / Errors**: `{self.is_valid_error}`",
            "",
            "## Fitted Parameters",
            "",
            "| Parameter | Raw Name | Nominal Value | Minuit Error | MCMC 68% Credible Interval | Status |",
            "| :--- | :--- | :---: | :---: | :---: | :---: |",
        ]

        for p in self.parameters:
            err_str = f"± {p['error']:.5f}" if p["error"] > 0 else "*fixed*"
            if p["mcmc_median"] is not None:
                mcmc_str = f"{p['mcmc_median']:.5f} (\\(-{p['mcmc_err_lo']:.4f}, +{p['mcmc_err_hi']:.4f}\\))"
            else:
                mcmc_str = "*N/A*"
            lines.append(f"| **{p['label']}** | `{p['name']}` | `{p['value']:.5f}` | {err_str} | {mcmc_str} | {p['status']} |")

        return "\n".join(lines)

    def save_reports(self) -> Tuple[Path, Path]:
        """Save text and markdown reports to the reports directory."""
        self.rundir.ensure_directories()
        txt = self.generate_ascii_table()
        md = self.generate_markdown_report()
        self.rundir.summary_txt.write_text(txt)
        self.rundir.summary_md.write_text(md)
        return self.rundir.summary_txt, self.rundir.summary_md


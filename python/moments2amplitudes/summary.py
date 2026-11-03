"""
Fit Summary and Diagnostics Table Generator for Moments2Amplitudes.

Extracts fitted parameters, goodness-of-fit metrics, and MCMC credible
intervals, generating terminal ASCII tables and saved markdown reports.
"""

from __future__ import annotations
import math
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple

from .schema import FitConfig
from .rundir import RunDirectory


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

    def load_from_root(self, root_file: Path | str, mcmc_file: Optional[Path | str] = None) -> bool:
        """Load nominal best-fit parameters and optional MCMC posterior intervals."""
        root_path = Path(root_file)
        if not root_path.is_file():
            return False

        # Attempt PyROOT or uproot
        try:
            import uproot
            return self._load_with_uproot(root_path, mcmc_file)
        except ImportError:
            pass

        try:
            import ROOT
            return self._load_with_pyroot(root_path, mcmc_file)
        except ImportError:
            pass

        return False

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
                    for p in self.parameters:
                        k = p["name"]
                        if k in mtree:
                            arr = np.array(mtree[k].array())
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

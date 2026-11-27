from __future__ import annotations
import os
import glob
import time
import re
try:
    import ROOT
except ImportError:
    ROOT = None
import numpy as np
import json
import gc
from datetime import datetime

try:
    import uproot
except ImportError:
    uproot = None


def _compute_autocorrelation_fft(data: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """
    Compute sample autocorrelation function acf(lag) = Cov(x_t, x_{t+lag}) / Var(x)
    matching BruFit's AutocorrPlot.cpp (parN / parD) via Wiener-Khinchin theorem FFT.
    """
    n = len(data)
    if n <= 1:
        return np.array([1.0], dtype=np.float64)
    xm = data - np.mean(data)
    par_d = np.sum(xm * xm)
    limit = min(n, (max_lag + 1) if max_lag is not None else (n // 2 + 1))
    if par_d <= 0.0:
        res = np.zeros(limit, dtype=np.float64)
        res[0] = 1.0
        return res
    n_fft = 2 ** int(np.ceil(np.log2(2 * n - 1)))
    fx = np.fft.fft(xm, n=n_fft)
    rx = np.fft.ifft(fx * np.conj(fx)).real[:limit]
    acf = np.asarray(rx / par_d, dtype=np.float64)
    acf[0] = 1.0
    return acf


class MomentsMCMC:

    """
    Helper class for loading, analyzing, and plotting BruFit MCMC moment results across mass bins.
    """

    def __init__(
        self,
        fit_result_dir: str,
        fit_bin_dir: list[str] | None = None,
        root_filename: str | None = None,
        burn_in: int = 2000
    ):
        """
        Initialize the MomentsMCMC container.

        Args:
            fit_result_dir: Directory containing mass bin fit outputs.
            fit_bin_dir: List of mass bin subdirectory names. If None, auto-discovered.
            root_filename: Name of the MCMC ROOT file (e.g. 'ResultsBruMcmcCovariance.root').
                           If None, auto-detected from candidates in the bin directories.
            burn_in: Number of initial MCMC steps to discard as burn-in (default: 2000).
        """
        self.fit_result_dir = fit_result_dir
        self.burn_in = burn_in
        self.__rdf = []
        self._moment_names = None
        self._bin_data_cache = {}
        self._diagnostics_cache = {}
        self.costheta_gj_canvas = []
        self.phi_gj_canvas = []
        self.acceptance_filename = "acceptance.root"
        self.acceptance = []
        self.is_acceptance_loaded = False

        # Auto-discover or set fit_bin_dir
        if fit_bin_dir is not None:
            self.fit_bin_dir = list(fit_bin_dir)
        else:
            self.fit_bin_dir = self._discover_mass_bin_dirs()

        # Parse mass bin center values from directory names (e.g. Mpi0eta1.020000 -> 1.02)
        self.fit_bin = []
        for d in self.fit_bin_dir:
            m = re.search(r"Mpi0eta([0-9.]+)", d)
            if m:
                self.fit_bin.append(float(m.group(1)))
            else:
                self.fit_bin.append(np.nan)

        # Auto-discover root filename if not explicitly provided
        if root_filename is not None:
            self.root_filename = root_filename
        else:
            self.root_filename = self._discover_root_filename()

        self.load_rdf()

    def set_burn_in(self, burn_in: int):
        """Set the burn-in period (number of initial samples to discard)."""
        self.burn_in = burn_in

    def _discover_mass_bin_dirs(self) -> list[str]:
        """Automatically find and sort mass bin subdirectories matching 'Mpi0eta*'."""
        if not os.path.isdir(self.fit_result_dir):
            return []
        subdirs = []
        for entry in os.listdir(self.fit_result_dir):
            full_p = os.path.join(self.fit_result_dir, entry)
            if os.path.isdir(full_p) and re.search(r"Mpi0eta[0-9.]+", entry):
                subdirs.append(entry)

        def _mass_key(d):
            m = re.search(r"Mpi0eta([0-9.]+)", d)
            return float(m.group(1)) if m else d

        subdirs.sort(key=_mass_key)
        return subdirs

    def _discover_root_filename(self) -> str:
        """Search candidate ROOT filenames in the first available bin directory or fit_result_dir."""
        candidates = [
            "ResultsBruMcmcFixedCovariance.root",
            "ResultsBruMcmcCovariance.root",
            "ResultsBruMcmc.root"
        ]
        check_dir = os.path.join(self.fit_result_dir, self.fit_bin_dir[0]) if self.fit_bin_dir else self.fit_result_dir
        for c in candidates:
            if os.path.exists(os.path.join(check_dir, c)):
                return c
        return "ResultsBruMcmcCovariance.root"

    # =========================================================================
    # Representation & Summary
    # =========================================================================

    def __repr__(self) -> str:
        n_bins = self.get_N_fit_bins()
        n_loaded = sum(1 for r in self.__rdf if r is not None)
        n_failed = n_bins - n_loaded
        return (f"<MomentsMCMC: {n_bins} mass bins ({n_loaded} loaded, {n_failed} failed/missing) "
                f"| burn_in={self.burn_in} | file='{self.root_filename}'>")

    def _repr_html_(self) -> str:
        """Rich HTML representation for Jupyter notebooks."""
        n_bins = self.get_N_fit_bins()
        n_loaded = sum(1 for r in self.__rdf if r is not None)
        n_failed = n_bins - n_loaded

        moments = self.get_moment_name_list()
        n_moments = len(moments)

        status_badge = (
            f"<span style='background-color: #28a745; color: white; padding: 2px 8px; border-radius: 4px; font-weight: bold;'>OK ({n_loaded}/{n_bins})</span>"
            if n_failed == 0 and n_bins > 0 else
            f"<span style='background-color: #e06c75; color: white; padding: 2px 8px; border-radius: 4px; font-weight: bold;'>{n_failed} Failed / {n_loaded} Loaded</span>"
        )

        return f"""
        <div style="border: 1px solid #d0d7de; border-radius: 8px; padding: 14px; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif; background: #ffffff; max-width: 680px; box-shadow: 0 1px 3px rgba(0,0,0,0.08);">
            <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 1px solid #eaecef; padding-bottom: 8px; margin-bottom: 10px;">
                <h4 style="margin: 0; color: #0969da; font-size: 16px;">MomentsMCMC Result Container</h4>
                {status_badge}
            </div>
            <table style="width: 100%; border-collapse: collapse; font-size: 13px; font-family: monospace;">
                <tr><td style="padding: 4px 0; font-weight: bold; width: 150px; color: #57606a;">Directory:</td><td style="padding: 4px 0; word-break: break-all;"><code>{self.fit_result_dir}</code></td></tr>
                <tr><td style="padding: 4px 0; font-weight: bold; color: #57606a;">ROOT File:</td><td style="padding: 4px 0;"><code>{self.root_filename}</code></td></tr>
                <tr><td style="padding: 4px 0; font-weight: bold; color: #57606a;">Mass Bins:</td><td style="padding: 4px 0;">{n_bins} bins ({n_loaded} successfully loaded)</td></tr>
                <tr><td style="padding: 4px 0; font-weight: bold; color: #57606a;">Burn-in:</td><td style="padding: 4px 0;"><b>{self.burn_in}</b> samples</td></tr>
                <tr><td style="padding: 4px 0; font-weight: bold; color: #57606a;">Moments:</td><td style="padding: 4px 0;">{n_moments} moments detected</td></tr>
            </table>
            <div style="margin-top: 10px; font-size: 11px; color: #6e7781; border-top: 1px solid #eaecef; padding-top: 6px;">
                Diagnostics: Call <code>mcmc.print_fit_status()</code> for fit status or <code>mcmc.print_mcmc_diagnostics()</code> for &tau;<sub>int</sub> & ESS.
            </div>
        </div>
        """

    def summary(self, print_table: bool = True) -> dict:
        """
        High-level summary of fit status and MCMC diagnostics.
        """
        status = self.get_fit_status(include_diagnostics=True)
        if print_table:
            self.print_fit_status()
        return status

    # =========================================================================
    # RDataFrame Loading & Management
    # =========================================================================

    def load_rdf(self):
        """
        Load RDataFrame for each mass bin directory. Failed or empty bins store None.
        """
        tree_name = "MCMCTree"
        if ROOT is None:
            print("ROOT is not available in this environment.")
            return

        def _try_load_rdf(fpath):
            if ROOT.gSystem.AccessPathName(fpath):
                print(f"File {fpath} not found.")
                return None
            try:
                tf = ROOT.TFile.Open(fpath)
                if not tf or tf.IsZombie():
                    print(f"File {fpath} is invalid or zombie.")
                    return None
                tree = tf.Get(tree_name)
                if not tree or tree.IsZombie():
                    tf.Close()
                    print(f"Tree '{tree_name}' not found in {fpath} (fit likely failed or incomplete).")
                    return None
                n_entries = tree.GetEntries()
                tf.Close()
                if n_entries <= 50:
                    print(f"Tree '{tree_name}' in {fpath} has only {n_entries} entries (incomplete chain).")
                    return None
                print(f"Loading file: {fpath}")
                return ROOT.RDataFrame(tree_name, fpath)
            except Exception as e:
                print(f"Failed to open/read {fpath}: {e}")
                return None

        self.__rdf = []
        if self.fit_bin_dir is None or len(self.fit_bin_dir) == 0:
            print("Loading single fit results from fit_result_dir.")
            filepath = f"{self.fit_result_dir}/{self.root_filename}"
            rdf = _try_load_rdf(filepath)
            if rdf is not None:
                self.__rdf.append(rdf)
            else:
                raise IndexError(f"{filepath} could not be loaded")
        else:
            for bdir in self.fit_bin_dir:
                filepath = f"{self.fit_result_dir}/{bdir}/{self.root_filename}"
                rdf = _try_load_rdf(filepath)
                self.__rdf.append(rdf)

    def print_rdf_fit_bin_dir(self):
        for bdir in self.fit_bin_dir:
            print(bdir)

    def get_N_fit_bins(self) -> int:
        return len(self.fit_bin_dir) if self.fit_bin_dir else len(self.__rdf)

    def get_N_rdf(self) -> int:
        return len(self.__rdf)

    def get_rdf(self, fit_bin_index: int):
        if 0 <= fit_bin_index < len(self.__rdf):
            return self.__rdf[fit_bin_index]
        raise IndexError(f"Mass bin index {fit_bin_index} out of range (0..{len(self.__rdf)-1}).")

    def get_moment_name_list(self) -> list[str]:
        """Get list of moment column names starting with 'H_'."""
        if self._moment_names is not None:
            return self._moment_names
        for rdf in self.__rdf:
            if rdf is not None:
                column_names = rdf.GetColumnNames()
                self._moment_names = [str(name) for name in column_names if str(name).startswith("H_")]
                return self._moment_names
        return []

    # =========================================================================
    # Normalized Moments Extraction & Plotting
    # =========================================================================

    def get_normalized_moments_distribution(self, fit_bin_index: int, var_name: str, burn_in: int = 0) -> ROOT.TH1D:
        """
        Get 1D distribution histogram of a normalized moment, discarding the first 'burn_in' steps,
        using 'weight' branch if present.
        """
        rdf = self.get_rdf(fit_bin_index)
        if rdf is not None:
            col_names = [str(c) for c in rdf.GetColumnNames()]
            has_wt = "weight" in col_names

            if burn_in > 0:
                if has_wt:
                    w_arr = np.asarray(rdf.AsNumpy(["weight"])["weight"], dtype=np.float64)
                    cum_w = np.cumsum(w_arr)
                    total_steps = cum_w[-1] if len(cum_w) > 0 else 0
                    eff_burn_steps = min(burn_in, total_steps // 2) if burn_in >= total_steps - 50 else burn_in
                    eff_burn_entry = int(np.searchsorted(cum_w, eff_burn_steps))
                    if 0 < eff_burn_entry < len(w_arr):
                        rdf = rdf.Range(eff_burn_entry, len(w_arr))
                else:
                    n_raw = rdf.Count().GetValue()
                    eff_burn = burn_in
                    if eff_burn >= n_raw - 50:
                        eff_burn = min(eff_burn, n_raw // 2)
                    if eff_burn > 0:
                        rdf = rdf.Range(eff_burn, n_raw)

            if has_wt:
                hist = rdf.Histo1D(var_name, "weight").GetValue()
            else:
                hist = rdf.Histo1D(var_name).GetValue()
            var_parts = var_name.split("_")
            if len(var_parts) == 4:
                H, alpha, L, M = var_parts
                hist.GetXaxis().SetTitle(f"H_{{{alpha}}}({L},{M})")
                hist.SetTitle("")
            else:
                hist.GetXaxis().SetTitle(var_name)
            hist.GetXaxis().SetTitleSize(0.05)
            hist.GetYaxis().SetTitleSize(0.05)
            return hist
        raise ValueError(f"RDF not loaded for mass bin index {fit_bin_index}.")

    def get_normalized_moments_stats(self, fit_bin_index: int, var_name: str, burn_in: int = None) -> tuple[float, float]:
        """
        Compute both mean and RMS error in a single histogram pass.
        
        Returns:
            (mean, error)
        """
        eff_burn_in = self.burn_in if burn_in is None else burn_in
        hist = self.get_normalized_moments_distribution(fit_bin_index, var_name, burn_in=eff_burn_in)
        if hist is not None:
            return float(hist.GetMean()), float(hist.GetRMS())
        raise ValueError(f"Failed to extract moment distribution for bin {fit_bin_index}")

    def get_normalized_moments_mean(self, fit_bin_index: int, var_name: str, burn_in: int = None) -> float:
        eff_burn_in = self.burn_in if burn_in is None else burn_in
        hist = self.get_normalized_moments_distribution(fit_bin_index, var_name, burn_in=eff_burn_in)
        return float(hist.GetMean())

    def get_normalized_moments_error(self, fit_bin_index: int, var_name: str, burn_in: int = None) -> float:
        eff_burn_in = self.burn_in if burn_in is None else burn_in
        hist = self.get_normalized_moments_distribution(fit_bin_index, var_name, burn_in=eff_burn_in)
        return float(hist.GetRMS())

    def get_normalized_moments_vs_mass(self, var_name: str, burn_in: int = None, skip_failed: bool = True):
        """
        Get normalized moment mean and error across all mass bins in a single pass per bin.

        Args:
            var_name: Moment name (e.g. 'H_0_1_0').
            burn_in: Burn-in steps to discard (default: self.burn_in).
            skip_failed: If True, omit bins where the fit failed or RDF is missing.

        Returns:
            (mass_bins, means, errors) as numpy arrays.
        """
        eff_burn_in = self.burn_in if burn_in is None else burn_in
        n_fit_bins = self.get_N_fit_bins()
        masses = []
        means = []
        errors = []

        for i in range(n_fit_bins):
            if i >= len(self.__rdf) or self.__rdf[i] is None:
                if skip_failed:
                    continue
                mass_val = self.fit_bin[i] if i < len(self.fit_bin) else i
                masses.append(mass_val)
                means.append(np.nan)
                errors.append(np.nan)
                continue

            try:
                mean, error = self.get_normalized_moments_stats(i, var_name, burn_in=eff_burn_in)
                mass_val = self.fit_bin[i] if i < len(self.fit_bin) else i
                masses.append(mass_val)
                means.append(mean)
                errors.append(error)
            except Exception as e:
                if not skip_failed:
                    raise e

        return np.array(masses, dtype=float), np.array(means, dtype=float), np.array(errors, dtype=float)

    def get_graph_normalized_moments_vs_mass(self, var_name: str, burn_in: int = None) -> ROOT.TGraphErrors:
        """Create a TGraphErrors of normalized moment vs mass, safely skipping failed bins."""
        fit_bins, means, error = self.get_normalized_moments_vs_mass(var_name, burn_in=burn_in, skip_failed=True)
        if len(fit_bins) == 0:
            raise ValueError(f"No valid bins available to build graph for {var_name}")

        gr = ROOT.TGraphErrors(len(fit_bins), fit_bins, means, np.zeros_like(fit_bins), error)
        var_parts = var_name.split("_")
        if len(var_parts) >= 4:
            alpha, L, M = var_parts[1:4]
            gr.SetName(f"gr_normalized_moments_{alpha}_{L}_{M}_vs_mass")
            gr.SetTitle(f"H_{{{alpha}}}({L},{M})")
        else:
            gr.SetName(f"gr_{var_name}_vs_mass")
            gr.SetTitle(var_name)
        gr.GetXaxis().SetTitle("M_{#eta#pi^{0}} (GeV/c^{2})")
        gr.GetYaxis().SetTitle("Normalized moment value")
        gr.SetMarkerStyle(20)
        return gr

    # =========================================================================
    # Acceptance & Unnormalized Moments
    # =========================================================================

    def load_acceptance(self):
        """Load acceptance values from acceptance.root, or default to 1.0."""
        self._bin_data_cache.clear()
        file_path = f"{self.fit_result_dir}/{self.acceptance_filename}"
        if not ROOT.gSystem.AccessPathName(file_path):
            print(f"Loading acceptance file: {file_path}")
            f = ROOT.TFile.Open(file_path, "READ")
            if f and not f.IsZombie():
                t = f.Get("acceptance")
                if t and not t.IsZombie():
                    n_entries = t.GetEntries()
                    self.acceptance = []
                    for i in range(n_entries):
                        t.GetEntry(i)
                        val = t.GetLeaf("acceptance").GetValue()
                        self.acceptance.append(float(val))
                    print(f"Loaded {len(self.acceptance)} acceptance values from {file_path}.")
                    self.is_acceptance_loaded = True
                    f.Close()
                else:
                    f.Close()
                    print(f"WARNING: Tree 'acceptance' not found in file {file_path}. Acceptance cannot be loaded.")
                    print("         The unnormalized moments are not yield acceptance corrected.")
                    self.acceptance = [1.0] * self.get_N_fit_bins()
                    self.is_acceptance_loaded = False
            else:
                print(f"WARNING: File {file_path} could not be opened. Acceptance cannot be loaded.")
                print("         The unnormalized moments are not yield acceptance corrected.")
                self.acceptance = [1.0] * self.get_N_fit_bins()
                self.is_acceptance_loaded = False
        else:
            print(f"WARNING: Acceptance file {file_path} not found. Acceptance cannot be loaded.")
            print("         The unnormalized moments are not yield acceptance corrected.")
            self.acceptance = [1.0] * self.get_N_fit_bins()
            self.is_acceptance_loaded = False

    def get_acceptance(self, fit_bin_index: int) -> float:
        if not self.acceptance:
            self.load_acceptance()
        if 0 <= fit_bin_index < len(self.acceptance):
            return self.acceptance[fit_bin_index]
        return 1.0

    def get_unnormalized_moments_distribution(self, fit_bin_index: int, var_name: str) -> ROOT.TH1D:
        """Get 1D distribution histogram of unnormalized moments (moment * Yld_Moments / 2), using weights if present."""
        rdf = self.get_rdf(fit_bin_index)
        if rdf is not None:
            col_names = [str(c) for c in rdf.GetColumnNames()]
            has_wt = "weight" in col_names
            product_col = f"Unnormalized_{var_name}"
            if var_name != "Yld_Moments" and var_name != "H_0_0_0":
                rdf_with_product = rdf.Define(product_col, f"{var_name}*Yld_Moments/2")
                if has_wt:
                    hist = rdf_with_product.Histo1D(product_col, "weight").GetValue()
                else:
                    hist = rdf_with_product.Histo1D(product_col).GetValue()
            else:
                if has_wt:
                    hist = rdf.Histo1D("Yld_Moments", "weight").GetValue()
                else:
                    hist = rdf.Histo1D("Yld_Moments").GetValue()
            var_parts = var_name.split("_")
            if len(var_parts) == 4:
                H, alpha, L, M = var_parts
                hist.GetXaxis().SetTitle(f"Unnormalized H_{{{alpha}}}({L},{M})")
                hist.SetTitle("")
            else:
                hist.GetXaxis().SetTitle(var_name)
            hist.GetXaxis().SetTitleSize(0.05)
            hist.GetYaxis().SetTitleSize(0.05)
            return hist
        raise ValueError(f"RDF not loaded for mass bin index {fit_bin_index}.")

    def _get_bin_numpy_data(self, fit_bin_index: int) -> dict | None:
        if fit_bin_index in self._bin_data_cache:
            return self._bin_data_cache[fit_bin_index]
        rdf = self.get_rdf(fit_bin_index)
        if rdf is None:
            return None
        col_names = self.get_moment_name_list()
        cols_to_load = list(col_names)
        if "Yld_Moments" not in cols_to_load:
            cols_to_load.append("Yld_Moments")
        col_all = [str(c) for c in rdf.GetColumnNames()]
        if "weight" in col_all and "weight" not in cols_to_load:
            cols_to_load.append("weight")
        data = rdf.AsNumpy(cols_to_load)
        self._bin_data_cache[fit_bin_index] = data
        return data

    def get_unnormalized_moments_stats(self, fit_bin_index: int, var_name: str) -> tuple[float, float]:
        """
        Compute both unnormalized moment mean and error in a single pass using fast cached numpy data.
        Uses sample weights if present.
        
        Returns:
            (mean, error)
        """
        acc = self.get_acceptance(fit_bin_index)
        acc_scale = acc if acc != 0 else 1.0

        data = self._get_bin_numpy_data(fit_bin_index)
        if data is not None and "Yld_Moments" in data:
            yld = data["Yld_Moments"]
            if var_name == "Yld_Moments" or var_name == "H_0_0_0":
                prod = yld
            elif var_name in data:
                prod = data[var_name] * yld / 2.0
            else:
                hist = self.get_unnormalized_moments_distribution(fit_bin_index, var_name)
                return float(hist.GetMean()) / acc_scale, float(hist.GetRMS()) / acc_scale

            weights = data.get("weight", None)
            if weights is not None:
                w = np.asarray(weights, dtype=np.float64)
                w_mean = np.average(prod, weights=w)
                w_var = np.average((prod - w_mean)**2, weights=w)
                return float(w_mean) / acc_scale, float(np.sqrt(w_var)) / acc_scale
            return float(np.mean(prod)) / acc_scale, float(np.std(prod)) / acc_scale

        hist = self.get_unnormalized_moments_distribution(fit_bin_index, var_name)
        return float(hist.GetMean()) / acc_scale, float(hist.GetRMS()) / acc_scale

    def get_unnormalized_moments_mean(self, fit_bin_index: int, var_name: str) -> float:
        mean, _ = self.get_unnormalized_moments_stats(fit_bin_index, var_name)
        return mean

    def get_unnormalized_moments_error(self, fit_bin_index: int, var_name: str) -> float:
        _, error = self.get_unnormalized_moments_stats(fit_bin_index, var_name)
        return error

    def get_unnormalized_moments_vs_mass(self, var_name: str, skip_failed: bool = True):
        """
        Get unnormalized moment mean and error across all mass bins in a single pass per bin.

        Args:
            var_name: Moment name (e.g. 'H_0_1_0' or 'Yld_Moments').
            skip_failed: If True, omit bins where the fit failed.

        Returns:
            (mass_bins, means, errors) as numpy arrays.
        """
        n_fit_bins = self.get_N_fit_bins()
        masses = []
        means = []
        errors = []

        for i in range(n_fit_bins):
            if i >= len(self.__rdf) or self.__rdf[i] is None:
                if skip_failed:
                    continue
                mass_val = self.fit_bin[i] if i < len(self.fit_bin) else i
                masses.append(mass_val)
                means.append(np.nan)
                errors.append(np.nan)
                continue

            try:
                mean, error = self.get_unnormalized_moments_stats(i, var_name)
                mass_val = self.fit_bin[i] if i < len(self.fit_bin) else i
                masses.append(mass_val)
                means.append(mean)
                errors.append(error)
            except Exception as e:
                if not skip_failed:
                    raise e

        return np.array(masses, dtype=float), np.array(means, dtype=float), np.array(errors, dtype=float)

    def get_graph_unnormalized_moments_vs_mass(self, var_name: str) -> ROOT.TGraphErrors:
        """Create a TGraphErrors of unnormalized moment vs mass, safely skipping failed bins."""
        query_var = "Yld_Moments" if var_name == "H_0_0_0" else var_name
        fit_bins, means, error = self.get_unnormalized_moments_vs_mass(query_var, skip_failed=True)
        if len(fit_bins) == 0:
            raise ValueError(f"No valid bins available to build graph for {var_name}")

        gr = ROOT.TGraphErrors(len(fit_bins), fit_bins, means, np.zeros_like(fit_bins), error)
        clean_name = "H_0_0_0" if var_name == "Yld_Moments" else var_name
        var_parts = clean_name.split("_")
        if len(var_parts) >= 4:
            alpha, L, M = var_parts[1:4]
            gr.SetName(f"gr_unnormalized_moments_{alpha}_{L}_{M}_vs_mass")
            gr.SetTitle(f"H_{{{alpha}}}({L},{M})")
        else:
            gr.SetName(f"gr_unnormalized_{clean_name}_vs_mass")
            gr.SetTitle(clean_name)
        gr.GetXaxis().SetTitle("M_{#eta#pi^{0}} (GeV/c^{2})")
        gr.GetYaxis().SetTitle("Unnormalized moment value")
        gr.SetMarkerStyle(20)
        return gr

    def save_unnormalized_moments_vs_mass(self, output_file: str):
        """
        Save unnormalized moments graphs vs mass to a ROOT file.
        Omit failed bins safely without raising unhandled exceptions.
        """
        f = ROOT.TFile.Open(output_file, "RECREATE")
        if not f or f.IsZombie():
            raise IOError(f"Cannot create file {output_file}")

        f.cd()
        moment_name_lists = self.get_moment_name_list()
        all_moments = list(moment_name_lists)
        if "H_0_0_0" not in all_moments:
            all_moments.append("H_0_0_0")

        graphs_to_save = []
        saved_count = 0
        for moment_name in all_moments:
            try:
                gr = self.get_graph_unnormalized_moments_vs_mass(moment_name)
                graphs_to_save.append(gr)
                f.WriteObject(gr, gr.GetName())
                saved_count += 1
            except Exception as e:
                print(f"Warning: Could not generate graph for {moment_name}: {e}")

        f.Write()
        f.Close()
        print(f"Successfully saved {saved_count} unnormalized moment graphs to {output_file}")

    def load_saved_unnormalized_moments_vs_mass(self, input_file: str, moment_name: str) -> ROOT.TGraphErrors | None:
        """
        Load a previously saved unnormalized moment graph from a ROOT file without leaking file descriptors.
        """
        f = ROOT.TFile.Open(input_file, "READ")
        if not f or f.IsZombie():
            raise FileNotFoundError(f"File {input_file} not found or could not be opened. Make sure the file was saved before.")

        target_name = moment_name.replace("H_", "gr_unnormalized_moments_") + "_vs_mass"
        gr = f.Get(target_name)
        if not gr or gr.IsZombie():
            gr = f.Get(f"gr_unnormalized_{moment_name}_vs_mass")

        if not gr or gr.IsZombie():
            f.Close()
            print(f"Graph '{target_name}' not found in file {input_file}.")
            return None

        gr_clone = gr.Clone()
        f.Close()
        return gr_clone

    def print_unnormalized_moments_vs_mass(self, output_file: str):
        """Print unnormalized moment vs mass plots to a multi-page PDF."""
        moment_name_lists = self.get_moment_name_list()
        if not moment_name_lists:
            print("No moments available to plot.")
            return

        plots_per_page = 3
        n_pages = (len(moment_name_lists) + plots_per_page - 1) // plots_per_page

        canvas = ROOT.TCanvas("c_unnorm", "Unnormalized Moments vs Mass", 1800, 600)
        canvas.Print(output_file + "[")

        for page in range(n_pages):
            canvas.Clear()
            canvas.Divide(plots_per_page, 1)
            start_idx = page * plots_per_page
            sublist = moment_name_lists[start_idx:start_idx + plots_per_page]
            for j, moment_name in enumerate(sublist):
                canvas.cd(j + 1)
                try:
                    gr = self.get_graph_unnormalized_moments_vs_mass(moment_name)
                    gr.SetMarkerColor(ROOT.kBlack)
                    gr.SetMarkerStyle(20)
                    gr.SetMarkerSize(1.5)
                    gr.SetLineColor(ROOT.kBlack)
                    gr.SetLineWidth(2)
                    gr.DrawClone("AP")
                except Exception as e:
                    print(f"Could not plot {moment_name}: {e}")
            canvas.Print(output_file)

        canvas.Print(output_file + "]")
        canvas.Close()
        print(f"Saved unnormalized moments plot to {output_file}")

    def print_saved_unnormalized_moments_vs_mass(self, input_file: str, output_file: str):
        """Print pre-saved unnormalized moment vs mass graphs from ROOT file to a multi-page PDF."""
        moment_name_lists = self.get_moment_name_list()
        if not moment_name_lists:
            print("No moments available to plot.")
            return

        plots_per_page = 3
        n_pages = (len(moment_name_lists) + plots_per_page - 1) // plots_per_page

        canvas = ROOT.TCanvas("c_saved_unnorm", "Saved Unnormalized Moments vs Mass", 1800, 600)
        canvas.Print(output_file + "[")

        for page in range(n_pages):
            canvas.Clear()
            canvas.Divide(plots_per_page, 1)
            start_idx = page * plots_per_page
            sublist = moment_name_lists[start_idx:start_idx + plots_per_page]
            for j, moment_name in enumerate(sublist):
                canvas.cd(j + 1)
                gr = self.load_saved_unnormalized_moments_vs_mass(input_file, moment_name)
                if gr is None:
                    continue
                gr.SetMarkerColor(ROOT.kBlack)
                gr.SetMarkerStyle(20)
                gr.SetMarkerSize(1.5)
                gr.SetLineColor(ROOT.kBlack)
                gr.SetLineWidth(2)
                gr.DrawClone("AP")
            canvas.Print(output_file)

        canvas.Print(output_file + "]")
        canvas.Close()
        gc.collect()
        print(f"Saved unnormalized moments plot to {output_file}")

    def get_graphs_unnormalized_moments_vs_mass_by_ess(
        self,
        var_name: str,
        ess_threshold: float = 100.0,
        burn_in: int = None,
        cached_diagnostics: dict = None
    ) -> tuple[ROOT.TGraphErrors, ROOT.TGraphErrors]:
        """
        Build dual TGraphErrors for unnormalized moments vs mass, categorizing points by ESS:
          - gr_good: Points with ESS >= ess_threshold (Black marker)
          - gr_bad:  Points with ESS < ess_threshold  (Red marker)

        Args:
            var_name: Moment name (e.g. 'H_0_1_0' or 'H_0_0_0').
            ess_threshold: Effective sample size threshold (default: 100.0).
            burn_in: Burn-in steps to discard.
            cached_diagnostics: Optional precomputed dict from self.get_all_mcmc_diagnostics().

        Returns:
            (gr_good, gr_bad) as ROOT.TGraphErrors
        """
        query_var = "Yld_Moments" if var_name == "H_0_0_0" else var_name
        clean_name = "H_0_0_0" if var_name == "Yld_Moments" else var_name
        n_fit_bins = self.get_N_fit_bins()

        if cached_diagnostics is not None:
            diags = cached_diagnostics
        else:
            diags = self.get_all_mcmc_diagnostics(burn_in=burn_in)

        good_masses, good_means, good_errors = [], [], []
        bad_masses, bad_means, bad_errors = [], [], []

        for i in range(n_fit_bins):
            if i >= len(self.__rdf) or self.__rdf[i] is None:
                continue

            try:
                mean, error = self.get_unnormalized_moments_stats(i, query_var)
                mass_val = self.fit_bin[i] if i < len(self.fit_bin) else float(i)

                # Look up ESS for this moment
                bin_diag = diags.get(i, {})
                per_m = bin_diag.get('per_moment', {})
                if clean_name in per_m:
                    ess = per_m[clean_name]['ess']
                else:
                    ess = self.get_mcmc_ess(i, clean_name, burn_in=burn_in)

                if ess >= ess_threshold:
                    good_masses.append(mass_val)
                    good_means.append(mean)
                    good_errors.append(error)
                else:
                    bad_masses.append(mass_val)
                    bad_means.append(mean)
                    bad_errors.append(error)
            except Exception:
                continue

        var_parts = clean_name.split("_")
        title_tex = f"H_{{{var_parts[1]}}}({var_parts[2]},{var_parts[3]})" if len(var_parts) >= 4 else clean_name

        n_good = len(good_masses)
        if n_good > 0:
            gr_good = ROOT.TGraphErrors(
                n_good,
                np.array(good_masses, dtype=float),
                np.array(good_means, dtype=float),
                np.zeros(n_good, dtype=float),
                np.array(good_errors, dtype=float)
            )
        else:
            gr_good = ROOT.TGraphErrors()
        gr_good.SetName(f"gr_unnorm_{clean_name}_good_ess")
        gr_good.SetTitle(f"{title_tex} (ESS #geq {int(ess_threshold)})")
        gr_good.SetMarkerStyle(20)
        gr_good.SetMarkerSize(1.2)
        gr_good.SetMarkerColor(ROOT.kBlack)
        gr_good.SetLineColor(ROOT.kBlack)
        gr_good.SetLineWidth(2)

        n_bad = len(bad_masses)
        if n_bad > 0:
            gr_bad = ROOT.TGraphErrors(
                n_bad,
                np.array(bad_masses, dtype=float),
                np.array(bad_means, dtype=float),
                np.zeros(n_bad, dtype=float),
                np.array(bad_errors, dtype=float)
            )
        else:
            gr_bad = ROOT.TGraphErrors()
        gr_bad.SetName(f"gr_unnorm_{clean_name}_bad_ess")
        gr_bad.SetTitle(f"{title_tex} (ESS < {int(ess_threshold)})")
        gr_bad.SetMarkerStyle(20)
        gr_bad.SetMarkerSize(1.2)
        gr_bad.SetMarkerColor(ROOT.kRed)
        gr_bad.SetLineColor(ROOT.kRed)
        gr_bad.SetLineWidth(2)

        return gr_good, gr_bad

    def print_unnormalized_moments_vs_mass_by_ess(
        self,
        output_file: str,
        ess_threshold: float = 100.0,
        burn_in: int = None,
        cached_diagnostics: dict = None
    ):
        """
        Print unnormalized moments vs mass to a multi-page PDF, distinguishing marker color by ESS:
          - Black: Success with good ESS (ESS >= threshold)
          - Red:   Success with bad ESS  (ESS < threshold)
        """
        moment_name_lists = self.get_moment_name_list()
        if not moment_name_lists:
            print("No moments available to plot.")
            return

        # Precompute diagnostics once for efficiency
        diags = cached_diagnostics if cached_diagnostics is not None else self.get_all_mcmc_diagnostics(burn_in=burn_in)

        plots_per_page = 3
        n_pages = (len(moment_name_lists) + plots_per_page - 1) // plots_per_page

        canvas = ROOT.TCanvas("c_unnorm_ess", "Unnormalized Moments vs Mass by ESS", 1800, 600)
        canvas.Print(output_file + "[")

        for page in range(n_pages):
            page_objects = []
            canvas.Clear()
            canvas.Divide(plots_per_page, 1)
            start_idx = page * plots_per_page
            sublist = moment_name_lists[start_idx:start_idx + plots_per_page]

            for j, moment_name in enumerate(sublist):
                canvas.cd(j + 1)
                ROOT.gPad.SetLeftMargin(0.15)
                ROOT.gPad.SetBottomMargin(0.15)
                ROOT.gPad.SetRightMargin(0.05)
                ROOT.gPad.SetTopMargin(0.10)
                ROOT.gPad.SetGridx(True)
                ROOT.gPad.SetGridy(True)

                gr_good, gr_bad = self.get_graphs_unnormalized_moments_vs_mass_by_ess(
                    moment_name,
                    ess_threshold=ess_threshold,
                    burn_in=burn_in,
                    cached_diagnostics=diags
                )

                has_good = gr_good.GetN() > 0
                has_bad = gr_bad.GetN() > 0

                if not has_good and not has_bad:
                    continue

                page_objects.extend([gr_good, gr_bad])

                var_parts = moment_name.split("_")
                title_tex = f"H_{{{var_parts[1]}}}({var_parts[2]},{var_parts[3]})" if len(var_parts) >= 4 else moment_name

                mg = ROOT.TMultiGraph()
                mg.SetName(f"mg_{moment_name}")
                mg.SetTitle(f"{title_tex};M_{{#eta#pi^{{0}}}} (GeV/c^{{2}});Unnormalized moment value")

                if has_good:
                    mg.Add(gr_good, "P")
                if has_bad:
                    mg.Add(gr_bad, "P")

                mg.Draw("A")
                page_objects.append(mg)

                mg.GetXaxis().SetTitleSize(0.05)
                mg.GetYaxis().SetTitleSize(0.05)
                mg.GetXaxis().SetTitleOffset(1.1)
                mg.GetYaxis().SetTitleOffset(1.2)

                # Reference zero line
                xmin = mg.GetXaxis().GetXmin()
                xmax = mg.GetXaxis().GetXmax()
                line_zero = ROOT.TLine(xmin, 0, xmax, 0)
                line_zero.SetLineColor(ROOT.kGray + 1)
                line_zero.SetLineStyle(2)
                line_zero.DrawClone("same")

                # Legend
                leg = ROOT.TLegend(0.60, 0.74, 0.92, 0.88)
                leg.SetBorderSize(1)
                leg.SetFillStyle(1001)
                leg.SetFillColor(ROOT.kWhite)
                leg.SetTextSize(0.038)
                if has_good:
                    leg.AddEntry(gr_good, f"ESS #geq {int(ess_threshold)} (Good)", "ep")
                if has_bad:
                    leg.AddEntry(gr_bad, f"ESS < {int(ess_threshold)} (Caution)", "ep")
                leg.DrawClone("same")

            canvas.Update()
            canvas.Print(output_file)

        canvas.Print(output_file + "]")
        canvas.Close()
        gc.collect()
        print(f"Saved ESS-colored unnormalized moments plot to {output_file}")

    def save_unnormalized_moments_vs_mass_by_ess(
        self,
        output_file: str,
        ess_threshold: float = 100.0,
        burn_in: int = None,
        cached_diagnostics: dict = None
    ):
        """
        Save unnormalized moments graphs vs mass categorized by ESS to a ROOT file.
        For each moment, saves:
          - gr_unnorm_<moment>_good_ess (Black points, ESS >= threshold)
          - gr_unnorm_<moment>_bad_ess  (Red points, ESS < threshold)
          - gr_unnormalized_moments_<alpha>_<L>_<M>_vs_mass (Combined graph for backward compatibility)
        """
        f = ROOT.TFile.Open(output_file, "RECREATE")
        if not f or f.IsZombie():
            raise IOError(f"Cannot create file {output_file}")

        f.cd()
        diags = cached_diagnostics if cached_diagnostics is not None else self.get_all_mcmc_diagnostics(burn_in=burn_in)
        moment_name_lists = self.get_moment_name_list()
        all_moments = list(moment_name_lists)
        if "H_0_0_0" not in all_moments:
            all_moments.append("H_0_0_0")

        saved_count = 0
        graphs_cache = []
        for moment_name in all_moments:
            try:
                # Combined graph
                gr_all = self.get_graph_unnormalized_moments_vs_mass(moment_name)
                graphs_cache.append(gr_all)
                f.WriteObject(gr_all, gr_all.GetName())

                # Dual graphs by ESS
                gr_good, gr_bad = self.get_graphs_unnormalized_moments_vs_mass_by_ess(
                    moment_name,
                    ess_threshold=ess_threshold,
                    burn_in=burn_in,
                    cached_diagnostics=diags
                )
                graphs_cache.extend([gr_good, gr_bad])
                f.WriteObject(gr_good, gr_good.GetName())
                f.WriteObject(gr_bad, gr_bad.GetName())
                saved_count += 1
            except Exception as e:
                print(f"Warning: Could not save graphs for {moment_name}: {e}")

        f.Write()
        f.Close()
        print(f"Successfully saved {saved_count} unnormalized moment graph sets to {output_file}")

    # =========================================================================
    # MCMC Trace & Autocorrelation Graphs
    # =========================================================================

    def get_mcmc_trace_graph(self, fit_bin_index: int, var_name: str) -> ROOT.TGraph:
        """Get MCMC trace graph for a specific moment, using weights to map steps if present."""
        rdf = self.get_rdf(fit_bin_index)
        if rdf is not None:
            col_names = [str(c) for c in rdf.GetColumnNames()]
            has_wt = "weight" in col_names
            if has_wt:
                data = rdf.AsNumpy([var_name, "weight"])
                vals = np.asarray(data[var_name], dtype=np.float64)
                wts = np.asarray(data["weight"], dtype=np.float64)

                # Check if weights can be represented as step durations
                if np.all(wts >= 1) and np.allclose(wts, np.round(wts)):
                    int_wts = np.maximum(1, np.round(wts)).astype(int)
                    total_steps = int(np.sum(int_wts))
                    if total_steps <= 2000000:
                        y_unrolled = np.repeat(vals, int_wts)
                        x_steps = np.arange(len(y_unrolled), dtype=np.float64)
                        graph = ROOT.TGraph(len(x_steps), x_steps, y_unrolled)
                    else:
                        cum_steps = np.concatenate(([0.0], np.cumsum(wts, dtype=np.float64)))
                        x_pts = np.empty(2 * len(vals), dtype=np.float64)
                        y_pts = np.empty(2 * len(vals), dtype=np.float64)
                        x_pts[0::2] = cum_steps[:-1]
                        x_pts[1::2] = cum_steps[1:]
                        y_pts[0::2] = vals
                        y_pts[1::2] = vals
                        graph = ROOT.TGraph(len(x_pts), x_pts, y_pts)
                else:
                    cum_steps = np.concatenate(([0.0], np.cumsum(wts, dtype=np.float64)))
                    x_pts = np.empty(2 * len(vals), dtype=np.float64)
                    y_pts = np.empty(2 * len(vals), dtype=np.float64)
                    x_pts[0::2] = cum_steps[:-1]
                    x_pts[1::2] = cum_steps[1:]
                    y_pts[0::2] = vals
                    y_pts[1::2] = vals
                    graph = ROOT.TGraph(len(x_pts), x_pts, y_pts)
            else:
                graph = rdf.Graph("entry", var_name)

            graph.SetTitle(f"MCMC Trace {var_name}")
            graph.GetXaxis().SetTitle("Steps")
            graph.GetXaxis().SetTitleSize(0.05)
            var_parts = var_name.split("_")
            if len(var_parts) == 4:
                H, alpha, L, M = var_parts
                graph.GetYaxis().SetTitle(f"H_{{{alpha}}}({L},{M})")
            else:
                graph.GetYaxis().SetTitle(var_name)
            graph.GetYaxis().SetTitleSize(0.05)
            return graph
        raise ValueError(f"RDF not loaded for mass bin index {fit_bin_index}.")

    def print_all_moment_trace_graph(self, output_file: str, true_value_file: str = None):
        """Print MCMC traces for all moments and bins, skipping failed bins safely."""
        moment_name_lists = self.get_moment_name_list()
        n_moments = len(moment_name_lists)
        n_fit_bins = self.get_N_fit_bins() if self.get_N_fit_bins() != 0 else 1

        print(f"Number of mass bins: {n_fit_bins}, Number of moments: {n_moments}")

        true_values = None
        if true_value_file:
            true_values = get_all_moment_true_values(true_value_file)

        plots_per_page = 12
        n_pages = (n_moments + plots_per_page - 1) // plots_per_page

        ROOT.TGaxis.SetMaxDigits(4)
        opened_pdf = False
        for i in range(n_fit_bins):
            if i >= len(self.__rdf) or self.__rdf[i] is None:
                print(f"Skipping MCMC trace plot for mass bin index {i} (no data)")
                continue

            mass_val = self.fit_bin[i] if i < len(self.fit_bin) else float(i)
            print(f"Plotting MCMC traces for mass bin index {i} (M={mass_val:.3f} GeV)")

            for page in range(n_pages):
                canvas = ROOT.TCanvas(f"c_page{page}_bin{i}", f"c_page{page}_bin{i}", 900, 1200)
                canvas.Divide(3, 4)

                if not opened_pdf:
                    canvas.Print(output_file + "[")
                    opened_pdf = True

                start_idx = page * plots_per_page
                end_idx = min(start_idx + plots_per_page, n_moments)

                for j in range(start_idx, end_idx):
                    pad_idx = j - start_idx + 1
                    canvas.cd(pad_idx)
                    ROOT.gPad.SetLeftMargin(0.18)
                    ROOT.gPad.SetBottomMargin(0.18)
                    ROOT.gPad.SetRightMargin(0.06)
                    ROOT.gPad.SetTopMargin(0.14)

                    m_name = moment_name_lists[j]
                    graph = self.get_mcmc_trace_graph(i, m_name)
                    var_parts = m_name.split("_")
                    if len(var_parts) >= 4:
                        alpha, L, M = var_parts[1:4]
                        moment_label = f"H_{{{alpha}}}({L},{M})"
                    else:
                        moment_label = m_name

                    graph.SetTitle("")
                    graph.GetXaxis().SetTitle("Steps")
                    graph.GetXaxis().SetTitleSize(0.052)
                    graph.GetXaxis().SetTitleOffset(1.25)
                    graph.GetXaxis().SetLabelSize(0.045)
                    graph.GetXaxis().SetLabelOffset(0.012)
                    graph.GetXaxis().SetNdivisions(504)

                    graph.GetYaxis().SetTitle(moment_label)
                    graph.GetYaxis().SetTitleSize(0.055)
                    graph.GetYaxis().SetTitleOffset(1.4)
                    graph.GetYaxis().SetLabelSize(0.045)
                    graph.GetYaxis().SetNdivisions(505)

                    graph.SetMarkerStyle(20)
                    graph.SetMarkerSize(0.3)
                    graph.SetLineColor(ROOT.kBlack)
                    graph.SetLineWidth(1)
                    graph.DrawClone("AL")

                    # Prominent moment, bin, and mass information in header
                    tex = ROOT.TLatex()
                    tex.SetNDC()
                    tex.SetTextFont(42)
                    tex.SetTextSize(0.046)
                    tex.SetTextColor(ROOT.kBlack)
                    tex.DrawLatex(0.18, 0.90, f"{moment_label}  -  Bin {i} (M = {mass_val:.3f} GeV)")

                    # Draw burn-in dashed vertical line if configured
                    eff_burn_in = getattr(self, 'burn_in', 2000)
                    if eff_burn_in and eff_burn_in > 0 and graph.GetN() > eff_burn_in:
                        ymin = graph.GetYaxis().GetXmin()
                        ymax = graph.GetYaxis().GetXmax()
                        bline = ROOT.TLine(eff_burn_in, ymin, eff_burn_in, ymax)
                        bline.SetLineColor(ROOT.kOrange + 7)
                        bline.SetLineStyle(2)
                        bline.SetLineWidth(1)
                        bline.DrawClone("same")

                    if true_values is not None:
                        true_value = true_values.get(m_name)
                        if true_value is not None:
                            line = ROOT.TLine(0, true_value * 2, graph.GetN(), true_value * 2)
                            line.SetLineColor(ROOT.kRed)
                            line.SetLineStyle(ROOT.kDashed)
                            line.SetLineWidth(2)
                            line.DrawClone("same")

                canvas.Update()
                canvas.Print(output_file)
                canvas.Close()

        if opened_pdf:
            temp_canvas = ROOT.TCanvas("c_close", "c_close", 100, 100)
            temp_canvas.Print(output_file + "]")
            temp_canvas.Close()
            print(f"Saved to {output_file}")

    def print_all_moment_posterior_distribution(
        self,
        fit_bin_index: int,
        output_file: str,
        burn_in: int = None
    ):
        """
        Print 1D posterior distributions of all moments for a single mass bin to a multi-page PDF.
        Includes Mean, RMS, tau_int, and ESS on each pad cleanly separated from the plot frame.
        """
        rdf = self.get_rdf(fit_bin_index)
        if rdf is None:
            print(f"Mass bin {fit_bin_index} has no valid RDF. Skipping posterior distribution plot.")
            return

        moment_name_lists = self.get_moment_name_list()
        n_moments = len(moment_name_lists)
        eff_burn_in = self.burn_in if burn_in is None else burn_in

        try:
            diag = self.get_mcmc_diagnostics(fit_bin_index, burn_in=eff_burn_in)
            per_moment = diag['per_moment']
        except Exception:
            per_moment = {}

        mass_val = self.fit_bin[fit_bin_index] if fit_bin_index < len(self.fit_bin) else float(fit_bin_index)
        plots_per_page = 12
        n_pages = (n_moments + plots_per_page - 1) // plots_per_page

        ROOT.TGaxis.SetMaxDigits(4)
        for page in range(n_pages):
            canvas = ROOT.TCanvas(f"c_post_p{page}_b{fit_bin_index}", f"Posterior Bin {fit_bin_index}", 900, 1200)
            canvas.Divide(3, 4)

            if page == 0:
                canvas.Print(output_file + "[")

            start_idx = page * plots_per_page
            end_idx = min(start_idx + plots_per_page, n_moments)

            for j in range(start_idx, end_idx):
                pad_idx = j - start_idx + 1
                canvas.cd(pad_idx)
                ROOT.gPad.SetLeftMargin(0.18)
                ROOT.gPad.SetBottomMargin(0.16)
                ROOT.gPad.SetRightMargin(0.06)
                ROOT.gPad.SetTopMargin(0.22)

                m_name = moment_name_lists[j]
                hist = self.get_normalized_moments_distribution(fit_bin_index, m_name, burn_in=eff_burn_in)
                hist.SetStats(0)
                hist.SetTitle("")

                var_parts = m_name.split("_")
                if len(var_parts) >= 4:
                    alpha, L, M = var_parts[1:4]
                    moment_label = f"H_{{{alpha}}}({L},{M})"
                else:
                    moment_label = m_name

                hist.GetXaxis().SetTitle(f"Normalized {moment_label}")
                hist.GetXaxis().SetTitleSize(0.052)
                hist.GetXaxis().SetTitleOffset(1.2)
                hist.GetXaxis().SetLabelSize(0.045)
                hist.GetXaxis().SetNdivisions(505)

                hist.GetYaxis().SetTitle("Posterior Samples")
                hist.GetYaxis().SetTitleSize(0.052)
                hist.GetYaxis().SetTitleOffset(1.4)
                hist.GetYaxis().SetLabelSize(0.045)
                hist.GetYaxis().SetNdivisions(505)

                hist.SetLineColor(ROOT.kAzure + 2)
                hist.SetFillColorAlpha(ROOT.kAzure - 9, 0.45)
                hist.SetLineWidth(2)
                hist.SetMaximum(hist.GetMaximum() * 1.15)
                hist.DrawClone("HIST")

                m_diag = per_moment.get(m_name, {})
                tau_val = m_diag.get('tau', 0.0)
                ess_val = m_diag.get('ess', 0.0)

                tex = ROOT.TLatex()
                tex.SetNDC()
                tex.SetTextFont(42)

                # Line 1: Moment, Bin, and Mass (well above plot frame at y=0.92)
                tex.SetTextSize(0.046)
                tex.SetTextColor(ROOT.kBlack)
                tex.DrawLatex(0.18, 0.92, f"{moment_label}  -  Bin {fit_bin_index} (M = {mass_val:.3f} GeV)")

                # Line 2: Mean +- RMS (left) and Tau | ESS (right) at y=0.84 (top frame border is at 0.78)
                tex.SetTextSize(0.040)
                tex.SetTextColor(ROOT.kBlack)
                tex.DrawLatex(0.18, 0.84, f"#mu = {hist.GetMean():.3f} #pm {hist.GetRMS():.3f}")

                ess_color = ROOT.kRed if ess_val < 100 else ROOT.kBlack
                tex.SetTextColor(ess_color)
                tex.DrawLatex(0.58, 0.84, f"#tau_{{int}} = {tau_val:.1f} | ESS = {ess_val:.1f}")

            canvas.Update()
            canvas.Print(output_file)
            canvas.Close()

        close_c = ROOT.TCanvas("c_close_post", "c_close_post", 100, 100)
        close_c.Print(output_file + "]")
        close_c.Close()
        print(f"Saved posterior distributions for bin {fit_bin_index} to {output_file}")

    def print_posterior_distribution_all_bins(self, output_file: str, burn_in: int = None):
        """
        Print 1D posterior distributions for all successful mass bins into a combined multi-page PDF.
        Includes Mean, RMS, tau_int, and ESS on each pad cleanly separated from the plot frame.
        """
        moment_name_lists = self.get_moment_name_list()
        n_moments = len(moment_name_lists)
        eff_burn_in = self.burn_in if burn_in is None else burn_in
        plots_per_page = 12
        n_pages_per_bin = (n_moments + plots_per_page - 1) // plots_per_page

        ROOT.TGaxis.SetMaxDigits(4)
        opened_pdf = False
        for i in range(self.get_N_fit_bins()):
            if i >= len(self.__rdf) or self.__rdf[i] is None:
                continue

            try:
                diag = self.get_mcmc_diagnostics(i, burn_in=eff_burn_in)
                per_moment = diag['per_moment']
            except Exception:
                per_moment = {}

            mass_val = self.fit_bin[i] if i < len(self.fit_bin) else float(i)
            print(f"Plotting posterior distributions for mass bin {i} (M={mass_val:.3f} GeV)")

            for page in range(n_pages_per_bin):
                canvas = ROOT.TCanvas(f"c_post_all_b{i}_p{page}", f"Bin {i} Page {page}", 900, 1200)
                canvas.Divide(3, 4)

                if not opened_pdf:
                    canvas.Print(output_file + "[")
                    opened_pdf = True

                start_idx = page * plots_per_page
                end_idx = min(start_idx + plots_per_page, n_moments)

                for j in range(start_idx, end_idx):
                    pad_idx = j - start_idx + 1
                    canvas.cd(pad_idx)
                    ROOT.gPad.SetLeftMargin(0.18)
                    ROOT.gPad.SetBottomMargin(0.16)
                    ROOT.gPad.SetRightMargin(0.06)
                    ROOT.gPad.SetTopMargin(0.22)

                    m_name = moment_name_lists[j]
                    hist = self.get_normalized_moments_distribution(i, m_name, burn_in=eff_burn_in)
                    hist.SetStats(0)
                    hist.SetTitle("")

                    var_parts = m_name.split("_")
                    if len(var_parts) >= 4:
                        alpha, L, M = var_parts[1:4]
                        moment_label = f"H_{{{alpha}}}({L},{M})"
                    else:
                        moment_label = m_name

                    hist.GetXaxis().SetTitle(f"Normalized {moment_label}")
                    hist.GetXaxis().SetTitleSize(0.052)
                    hist.GetXaxis().SetTitleOffset(1.2)
                    hist.GetXaxis().SetLabelSize(0.045)
                    hist.GetXaxis().SetNdivisions(505)

                    hist.GetYaxis().SetTitle("Posterior Samples")
                    hist.GetYaxis().SetTitleSize(0.052)
                    hist.GetYaxis().SetTitleOffset(1.4)
                    hist.GetYaxis().SetLabelSize(0.045)
                    hist.GetYaxis().SetNdivisions(505)

                    hist.SetLineColor(ROOT.kAzure + 2)
                    hist.SetFillColorAlpha(ROOT.kAzure - 9, 0.45)
                    hist.SetLineWidth(2)
                    hist.SetMaximum(hist.GetMaximum() * 1.15)
                    hist.DrawClone("HIST")

                    m_diag = per_moment.get(m_name, {})
                    tau_val = m_diag.get('tau', 0.0)
                    ess_val = m_diag.get('ess', 0.0)

                    tex = ROOT.TLatex()
                    tex.SetNDC()
                    tex.SetTextFont(42)

                    # Line 1: Moment, Bin, and Mass (well above plot frame at y=0.92)
                    tex.SetTextSize(0.046)
                    tex.SetTextColor(ROOT.kBlack)
                    tex.DrawLatex(0.18, 0.92, f"{moment_label}  -  Bin {i} (M = {mass_val:.3f} GeV)")

                    # Line 2: Mean +- RMS (left) and Tau | ESS (right) at y=0.84 (top frame border is at 0.78)
                    tex.SetTextSize(0.040)
                    tex.SetTextColor(ROOT.kBlack)
                    tex.DrawLatex(0.18, 0.84, f"#mu = {hist.GetMean():.3f} #pm {hist.GetRMS():.3f}")

                    ess_color = ROOT.kRed if ess_val < 100 else ROOT.kBlack
                    tex.SetTextColor(ess_color)
                    tex.DrawLatex(0.58, 0.84, f"#tau_{{int}} = {tau_val:.1f} | ESS = {ess_val:.1f}")

                canvas.Update()
                canvas.Print(output_file)
                canvas.Close()

        if opened_pdf:
            close_c = ROOT.TCanvas("c_close_all_post", "c_close_all_post", 100, 100)
            close_c.Print(output_file + "]")
            close_c.Close()
            print(f"Saved all posterior distributions to {output_file}")

    def _get_moment_chain_data(self, fit_bin_index: int, var_name: str, burn_in: int = None) -> tuple[np.ndarray, int, int]:
        """
        Extract moment MCMC chain data, unrolling by 'weight' if present and applying burn-in.

        Returns:
            (chain_data_after_burn, total_raw_steps, eff_burn_in_steps)
        """
        rdf = self.get_rdf(fit_bin_index)
        if rdf is None:
            raise ValueError(f"RDF not loaded for mass bin index {fit_bin_index}.")

        col_names = [str(c) for c in rdf.GetColumnNames()]
        has_wt = "weight" in col_names

        if has_wt:
            raw_dict = rdf.AsNumpy([var_name, "weight"])
            v = np.asarray(raw_dict[var_name], dtype=np.float64)
            w = np.asarray(raw_dict["weight"], dtype=np.float64)
            int_w = np.maximum(1, np.round(w)).astype(int)
            data = np.repeat(v, int_w)
            n_raw = len(data)
        else:
            data = np.asarray(rdf.AsNumpy([var_name])[var_name], dtype=np.float64)
            n_raw = len(data)

        eff_burn_in = self.burn_in if burn_in is None else burn_in
        if eff_burn_in >= n_raw - 50:
            eff_burn_in = min(eff_burn_in, n_raw // 2)
        if eff_burn_in > 0:
            data = data[eff_burn_in:]

        return data, n_raw, eff_burn_in

    def get_mcmc_autocorrelation_graph(self, fit_bin_index: int, var_name: str, max_lag: int = None, burn_in: int = None) -> ROOT.TGraph:
        """Get autocorrelation graph for a specific moment, using weights if present."""
        data, _, _ = self._get_moment_chain_data(fit_bin_index, var_name, burn_in=burn_in)
        n_samples = len(data)
        if max_lag is None:
            max_lag = min(5000, n_samples // 2)

        acf = _compute_autocorrelation_fft(data, max_lag=max_lag)
        n_pts = len(acf)
        graph = ROOT.TGraph(n_pts)
        graph.SetTitle(f"MCMC Autocorrelation {var_name}")
        graph.GetXaxis().SetTitle("Lag")
        graph.GetXaxis().SetTitleSize(0.05)
        graph.GetYaxis().SetTitle("Normalized Autocorrelation")
        graph.GetYaxis().SetTitleSize(0.05)
        for lag in range(n_pts):
            graph.SetPoint(lag, float(lag), float(acf[lag]))

        return graph

    def print_all_moment_mcmc_autocorrelation(self, fit_bin_index: int, output_file: str):
        """Print autocorrelation curves for all moments for a single bin to a multi-page PDF, using weights if present."""
        rdf = self.get_rdf(fit_bin_index)
        if rdf is None:
            print(f"Mass bin {fit_bin_index} has no valid RDF. Skipping autocorrelation plot.")
            return

        moment_name_lists = self.get_moment_name_list()
        n_moments = len(moment_name_lists)
        print(f"Plotting autocorrelation for {n_moments} moments in bin {fit_bin_index}")

        plots_per_page = 12
        n_pages = (n_moments + plots_per_page - 1) // plots_per_page

        for page in range(n_pages):
            canvas = ROOT.TCanvas(f"c_acf_page{page}", f"c_acf_page{page}", 900, 1200)
            canvas.Divide(3, 4)

            if page == 0:
                canvas.Print(output_file + "[")

            start_idx = page * plots_per_page
            end_idx = min(start_idx + plots_per_page, n_moments)

            for j in range(start_idx, end_idx):
                pad_idx = j - start_idx + 1
                canvas.cd(pad_idx)
                ROOT.gPad.SetLeftMargin(0.18)
                ROOT.gPad.SetBottomMargin(0.15)
                ROOT.gPad.SetRightMargin(0.05)
                ROOT.gPad.SetTopMargin(0.1)

                m_name = moment_name_lists[j]
                data, _, _ = self._get_moment_chain_data(fit_bin_index, m_name)
                n_samples = len(data)
                eff_max_lag = min(5000, n_samples // 2)
                acf = _compute_autocorrelation_fft(data, max_lag=eff_max_lag)
                n_pts = len(acf)

                graph = ROOT.TGraph(n_pts)
                var_parts = m_name.split("_")
                if len(var_parts) >= 4:
                    alpha, L, M = var_parts[1:4]
                    graph.SetTitle(f"H_{{{alpha}}}({L},{M})")
                graph.GetXaxis().SetTitle("Lag")
                graph.GetYaxis().SetTitle("Normalized Autocorrelation")
                graph.GetXaxis().SetTitleSize(0.05)
                graph.GetYaxis().SetTitleSize(0.05)

                for lag in range(n_pts):
                    graph.SetPoint(lag, float(lag), float(acf[lag]))

                graph.SetMarkerStyle(20)
                graph.SetMarkerSize(0.3)
                graph.SetLineColor(ROOT.kBlack)
                graph.SetLineWidth(1)
                graph.DrawClone("AL")

                line_zero = ROOT.TLine(0, 0, n_pts - 1, 0)
                line_zero.SetLineColor(ROOT.kGray)
                line_zero.SetLineStyle(2)
                line_zero.DrawClone("same")

            canvas.Update()
            canvas.Print(output_file)
            canvas.Close()

        close_canvas = ROOT.TCanvas("c_close_acf", "c_close_acf", 100, 100)
        close_canvas.Print(output_file + "]")
        close_canvas.Close()
        print(f"Saved to {output_file}")

    # =========================================================================
    # MCMC Chain Diagnostics (Integrated Autocorrelation & ESS)
    # =========================================================================

    def get_mcmc_integrated_autocorrelation_time(self, fit_bin_index: int, var_name: str, max_lag: int = None, burn_in: int = None) -> float:
        """
        Compute integrated autocorrelation time: tau_int = 0.5 + sum_{k=1} acf(k)
        Stops summation when acf(k) <= 0 (matching BruFit's Geyer sequence estimator).
        Uses sample weights if present.
        """
        data, _, _ = self._get_moment_chain_data(fit_bin_index, var_name, burn_in=burn_in)
        n_samples = len(data)
        if n_samples <= 1:
            return 0.5
        eff_max_lag = (n_samples // 2) if max_lag is None else min(max_lag, n_samples // 2)
        acf = _compute_autocorrelation_fft(data, max_lag=eff_max_lag)

        tau_int = 0.5
        for lag in range(1, len(acf)):
            if acf[lag] <= 0:
                break
            tau_int += acf[lag]

        return float(tau_int)

    def get_sum_mcmc_integrated_autocorrelation_time(self, fit_bin_index: int, max_lag: int = None, burn_in: int = None) -> float:
        """Compute sum of tau_int across all moments using fast vectorized diagnostics."""
        diag = self.get_mcmc_diagnostics(fit_bin_index, max_lag=max_lag, burn_in=burn_in)
        return float(sum(vals['tau'] for vals in diag['per_moment'].values()))

    def get_mcmc_ess(self, fit_bin_index: int, var_name: str, max_lag: int = None, burn_in: int = None) -> float:
        """
        Calculate Effective Sample Size: ESS = N_eff / (2 * tau_int)
        Uses sample weights if present.
        """
        data, _, _ = self._get_moment_chain_data(fit_bin_index, var_name, burn_in=burn_in)
        n_eff = len(data)
        if n_eff == 0:
            return 0.0

        tau_int = self.get_mcmc_integrated_autocorrelation_time(fit_bin_index, var_name, max_lag=max_lag, burn_in=burn_in)
        if tau_int > 0:
            ess = n_eff / (2.0 * tau_int)
            return min(float(n_eff), float(ess))
        return float(n_eff)

    def get_mcmc_diagnostics(self, fit_bin_index: int, max_lag: int = None, burn_in: int = None) -> dict:
        """
        Fast vectorized MCMC diagnostics across all moments in a mass bin (loads all columns once).
        Uses sample weights if present and caches computed results per (bin, burn_in, max_lag).

        Args:
            fit_bin_index: Mass bin index.
            max_lag: Maximum lag for autocorrelation calculation.
            burn_in: Initial samples to discard (default: self.burn_in).

        Returns:
            Dictionary containing:
                - 'bin': mass bin index
                - 'mass': mass value
                - 'n_raw': total MCMC samples (sum of weights if weighted)
                - 'burn_in': burn-in samples discarded
                - 'n_samples': effective samples evaluated
                - 'mean_tau', 'min_tau', 'max_tau'
                - 'mean_ess', 'min_ess', 'max_ess'
                - 'per_moment': dict mapping moment -> {'tau': tau, 'ess': ess}
        """
        eff_burn_in = self.burn_in if burn_in is None else burn_in
        cache_key = (fit_bin_index, eff_burn_in, max_lag)
        if cache_key in self._diagnostics_cache:
            return self._diagnostics_cache[cache_key]

        all_data = self._get_bin_numpy_data(fit_bin_index)
        if all_data is None:
            raise ValueError(f"RDF not loaded for mass bin index {fit_bin_index} (fit may have failed).")

        moment_name_lists = self.get_moment_name_list()
        if not moment_name_lists:
            raise ValueError("No moments found in RDF.")

        batch = compute_batch_chain_diagnostics(
            all_data,
            moment_name_lists,
            burn_in=eff_burn_in,
            max_lag=max_lag,
        )
        mass_val = self.fit_bin[fit_bin_index] if fit_bin_index < len(self.fit_bin) else None
        if batch is None:
            n_tree = len(all_data[moment_name_lists[0]])
            per_moment = {name: {'tau': 0.5, 'ess': float(n_tree)} for name in moment_name_lists}
            result = {
                'bin': fit_bin_index,
                'mass': mass_val,
                'n_raw': int(n_tree),
                'burn_in': int(min(eff_burn_in, n_tree // 2)),
                'n_samples': int(n_tree),
                'mean_tau': 0.5,
                'min_tau': 0.5,
                'max_tau': 0.5,
                'mean_ess': float(n_tree),
                'min_ess': float(n_tree),
                'max_ess': float(n_tree),
                'per_moment': per_moment,
            }
        else:
            result = {
                'bin': fit_bin_index,
                'mass': mass_val,
                'n_raw': batch['n_raw'],
                'burn_in': batch['burn_in'],
                'n_samples': batch['n_samples'],
                'mean_tau': batch['tau_mean'],
                'min_tau': batch['tau_min'],
                'max_tau': batch['tau_max'],
                'mean_ess': batch['ess_mean'],
                'min_ess': batch['ess_min'],
                'max_ess': batch['ess_max'],
                'per_moment': batch['per_moment'],
            }

        self._diagnostics_cache[cache_key] = result
        return result

    def get_all_mcmc_diagnostics(self, max_lag: int = None, burn_in: int = None) -> dict:
        """Compute MCMC diagnostics for all successfully loaded mass bins."""
        diagnostics = {}
        for i in range(self.get_N_fit_bins()):
            if i < len(self.__rdf) and self.__rdf[i] is not None:
                try:
                    diagnostics[i] = self.get_mcmc_diagnostics(i, max_lag=max_lag, burn_in=burn_in)
                except Exception as e:
                    print(f"Warning: Failed to compute diagnostics for bin {i}: {e}")
        return diagnostics

    def print_mcmc_diagnostics(self, fit_bin_index: int = None, max_lag: int = None, burn_in: int = None):
        """
        Print formatted integrated autocorrelation time and ESS diagnostics.
        If fit_bin_index is provided, prints detailed moment-by-moment table.
        If None, prints summary across all successful bins.
        """
        if fit_bin_index is not None:
            if fit_bin_index < 0 or fit_bin_index >= len(self.__rdf) or self.__rdf[fit_bin_index] is None:
                print(f"Cannot display diagnostics: mass bin {fit_bin_index} has no valid MCMC data (fit failed or missing).")
                return None
            diag = self.get_mcmc_diagnostics(fit_bin_index, max_lag=max_lag, burn_in=burn_in)
            mass_str = f"{diag['mass']:.4f} GeV" if diag['mass'] is not None else f"Bin {fit_bin_index}"
            print("=" * 80)
            print(f"MCMC Diagnostics for Mass Bin {fit_bin_index} ({mass_str})")
            print(f"Total Samples: {diag['n_raw']} | Burn-in: {diag['burn_in']} | Analyzed Samples: {diag['n_samples']}")
            print("=" * 80)
            print(f"{'Moment':<25} {'Integrated Tau (tau_int)':<30} {'Effective Sample Size (ESS)':<25}")
            print("-" * 80)
            for m_name, vals in diag['per_moment'].items():
                print(f"{m_name:<25} {vals['tau']:<30.2f} {vals['ess']:<25.1f}")
            print("-" * 80)
            print(f"Summary: Mean Tau = {diag['mean_tau']:.2f} (max: {diag['max_tau']:.2f}) | "
                  f"Min ESS = {diag['min_ess']:.1f} (mean: {diag['mean_ess']:.1f})")
            print("=" * 80)
            return diag

        all_diags = self.get_all_mcmc_diagnostics(max_lag=max_lag, burn_in=burn_in)
        eff_burn_in = self.burn_in if burn_in is None else burn_in
        print("=" * 105)
        print(f"MCMC Diagnostics Summary Across All Successful Mass Bins (configured burn_in={eff_burn_in})")
        print("=" * 105)
        print(f"{'Bin':<5} {'Mass (GeV)':<12} {'N Samples':<12} {'Tau (Mean)':<14} {'Tau (Max)':<14} {'ESS (Min)':<14} {'ESS (Mean)':<14}")
        print("-" * 105)
        for b_idx, d in sorted(all_diags.items()):
            m_str = f"{d['mass']:.4f}" if d['mass'] is not None else "N/A"
            print(f"{b_idx:<5} {m_str:<12} {d['n_samples']:<12} {d['mean_tau']:<14.2f} {d['max_tau']:<14.2f} {d['min_ess']:<14.1f} {d['mean_ess']:<14.1f}")
        print("=" * 105)
        return all_diags

    def get_fit_status(
        self,
        log_files_or_dir: str | list[str] = None,
        active_timeout_min: float = 30.0,
        as_dataframe: bool = False,
        include_diagnostics: bool = True
    ) -> dict:
        """Extract fit status for this MCMC dataset from log files."""
        if log_files_or_dir is None:
            log_files_or_dir = self.fit_result_dir
        return get_fit_status(
            log_files_or_dir=log_files_or_dir,
            active_timeout_min=active_timeout_min,
            as_dataframe=as_dataframe,
            mcmc=self if include_diagnostics else None,
            burn_in=self.burn_in,
            compute_diagnostics=include_diagnostics
        )

    def print_fit_status(
        self,
        log_files_or_dir: str | list[str] = None,
        filter_status: str = None,
        active_timeout_min: float = 30.0,
        show_diagnostics: bool = True
    ):
        """Print fit status summary table from log files."""
        if log_files_or_dir is None:
            log_files_or_dir = self.fit_result_dir
        return print_fit_status(
            log_files_or_dir=log_files_or_dir,
            filter_status=filter_status,
            active_timeout_min=active_timeout_min,
            mcmc=self if show_diagnostics else None,
            burn_in=self.burn_in,
            show_diagnostics=show_diagnostics
        )

    # =========================================================================
    # Angular Distributions
    # =========================================================================

    def load_angular_distribution_canvas(self):
        """Load angular distribution canvases from fit output ROOT files."""
        self.costheta_gj_canvas = []
        self.phi_gj_canvas = []
        for i, fit_bin_dir in enumerate(self.fit_bin_dir):
            filePath = f"{self.fit_result_dir}/{fit_bin_dir}/{self.root_filename}"
            if not os.path.exists(filePath):
                self.costheta_gj_canvas.append(None)
                self.phi_gj_canvas.append(None)
                continue

            f = ROOT.TFile.Open(filePath, "READ")
            if f and not f.IsZombie():
                cosTheta_gj_histname = None
                phi_gj_histname = None
                for key in f.GetListOfKeys():
                    name = key.GetName()
                    if "cosTheta_eta_gj" in name:
                        cosTheta_gj_histname = name
                    if "phi_eta_gj" in name:
                        phi_gj_histname = name

                c_costheta = f.Get(cosTheta_gj_histname).Clone() if cosTheta_gj_histname else None
                c_phi = f.Get(phi_gj_histname).Clone() if phi_gj_histname else None
                self.costheta_gj_canvas.append(c_costheta)
                self.phi_gj_canvas.append(c_phi)
                f.Close()
            else:
                self.costheta_gj_canvas.append(None)
                self.phi_gj_canvas.append(None)

    def get_phi_gj_data(self, fit_bin_index: int) -> ROOT.TH1 | None:
        """Extract phi_GJ data histogram without console output spam."""
        if fit_bin_index >= len(self.phi_gj_canvas) or self.phi_gj_canvas[fit_bin_index] is None:
            return None
        pad = self.phi_gj_canvas[fit_bin_index].GetPad(1)
        if not pad:
            return None
        primitives = pad.GetListOfPrimitives()
        for i in range(primitives.GetSize()):
            obj_name = primitives.At(i).GetName()
            if "h_Data" in obj_name:
                return pad.GetPrimitive(obj_name)
        return None

    def get_costheta_gj_data(self, fit_bin_index: int) -> ROOT.TH1 | None:
        """Extract cosTheta_GJ data histogram without console output spam."""
        if fit_bin_index >= len(self.costheta_gj_canvas) or self.costheta_gj_canvas[fit_bin_index] is None:
            return None
        pad = self.costheta_gj_canvas[fit_bin_index].GetPad(1)
        if not pad:
            return None
        primitives = pad.GetListOfPrimitives()
        for i in range(primitives.GetSize()):
            obj_name = primitives.At(i).GetName()
            if "h_Data" in obj_name:
                return pad.GetPrimitive(obj_name)
        return None

    def get_phi_gj_fit(self, fit_bin_index: int) -> ROOT.TH1 | None:
        """Extract phi_GJ fit histogram without console output spam."""
        if fit_bin_index >= len(self.phi_gj_canvas) or self.phi_gj_canvas[fit_bin_index] is None:
            return None
        pad = self.phi_gj_canvas[fit_bin_index].GetPad(1)
        if not pad:
            return None
        primitives = pad.GetListOfPrimitives()
        for i in range(primitives.GetSize()):
            obj_name = primitives.At(i).GetName()
            if "TotalPDF" in obj_name:
                return pad.GetPrimitive(obj_name)
        return None

    def get_costheta_gj_fit(self, fit_bin_index: int) -> ROOT.TH1 | None:
        """Extract cosTheta_GJ fit histogram without console output spam."""
        if fit_bin_index >= len(self.costheta_gj_canvas) or self.costheta_gj_canvas[fit_bin_index] is None:
            return None
        pad = self.costheta_gj_canvas[fit_bin_index].GetPad(1)
        if not pad:
            return None
        primitives = pad.GetListOfPrimitives()
        for i in range(primitives.GetSize()):
            obj_name = primitives.At(i).GetName()
            if "TotalPDF" in obj_name:
                return pad.GetPrimitive(obj_name)
        return None

    def print_phi_eta_gj(self, output_file: str):
        canvas = ROOT.TCanvas("c_phi_eta_gj", "phi_eta_gj", 1200, 800)
        canvas.Divide(6, 5)
        for i in range(self.get_N_fit_bins()):
            canvas.cd(i + 1)
            h_data = self.get_phi_gj_data(i)
            h_fit = self.get_phi_gj_fit(i)
            if h_data is None or h_fit is None:
                continue
            ROOT.gPad.SetLeftMargin(0.15)
            ROOT.gPad.SetRightMargin(0.02)
            ROOT.gPad.SetBottomMargin(0.15)
            mass_val = self.fit_bin[i] if i < len(self.fit_bin) else 0.0
            h_data.SetTitle(f"M_{{#eta#pi^{{0}}}} = {mass_val:0.2f} GeV")
            h_data.GetXaxis().SetTitle("#phi_{GJ} (rad)")
            h_data.GetYaxis().SetTitle("counts")
            h_data.GetXaxis().SetTitleSize(0.05)
            h_data.GetYaxis().SetTitleSize(0.05)
            h_data.GetXaxis().SetTitleOffset(1.0)
            h_data.GetXaxis().SetRangeUser(-3.2, 3.2)
            h_data.GetYaxis().SetTitleOffset(1.0)
            h_data.SetMarkerSize(0.2)
            h_data.DrawClone()
            h_fit.SetLineColor(ROOT.kRed)
            h_fit.SetLineWidth(2)
            h_fit.DrawClone("same")
            leg = ROOT.TLegend(0.5, 0.75, 0.65, 0.9)
            leg.AddEntry(h_data, "Data", "lep")
            leg.AddEntry(h_fit, "Fit", "l")
            leg.DrawClone()

        canvas.Print(output_file)
        canvas.Close()
        print(f"Saved to {output_file}")

    def print_costheta_eta_gj(self, output_file: str):
        canvas = ROOT.TCanvas("c_costheta_eta_gj", "costheta_eta_gj", 1200, 800)
        canvas.Divide(6, 5)
        for i in range(self.get_N_fit_bins()):
            canvas.cd(i + 1)
            h_data = self.get_costheta_gj_data(i)
            h_fit = self.get_costheta_gj_fit(i)
            if h_data is None or h_fit is None:
                continue
            ROOT.gPad.SetLeftMargin(0.15)
            ROOT.gPad.SetRightMargin(0.02)
            ROOT.gPad.SetBottomMargin(0.15)
            mass_val = self.fit_bin[i] if i < len(self.fit_bin) else 0.0
            h_data.SetTitle(f"M_{{#eta#pi^{{0}}}} = {mass_val:0.2f} GeV")
            h_data.GetXaxis().SetTitle("cos(#theta_{GJ})")
            h_data.GetYaxis().SetTitle("counts")
            h_data.GetXaxis().SetTitleSize(0.05)
            h_data.GetYaxis().SetTitleSize(0.05)
            h_data.GetXaxis().SetTitleOffset(1.0)
            h_data.GetXaxis().SetRangeUser(-1, 1)
            h_data.GetYaxis().SetTitleOffset(1.0)
            h_data.SetMarkerSize(0.2)
            h_data.DrawClone()
            h_fit.SetLineColor(ROOT.kRed)
            h_fit.SetLineWidth(2)
            h_fit.DrawClone("same")
            leg = ROOT.TLegend(0.5, 0.75, 0.65, 0.9)
            leg.AddEntry(h_data, "Data", "lep")
            leg.AddEntry(h_fit, "Fit", "l")
            leg.DrawClone()

        canvas.Print(output_file)
        canvas.Close()
        print(f"Saved to {output_file}")


# =============================================================================
# Standalone Functions & Helpers
# =============================================================================

def get_moment_true_value(json_file: str, var_name: str) -> float | None:
    """Get the true value of a moment from a JSON file."""
    with open(json_file, 'r') as f:
        true_values = json.load(f)
    split_var_name = var_name.split("_")
    if len(split_var_name) == 4:
        H, alpha, L, M = split_var_name
        for true_value_dict in true_values:
            true_alpha = true_value_dict["qn"]["momentIndex"]
            true_L = true_value_dict["qn"]["L"]
            true_M = true_value_dict["qn"]["M"]
            if (int(alpha) == int(true_alpha)) and (int(L) == int(true_L)) and (int(M) == int(true_M)):
                return true_value_dict["val_Im"] if int(alpha) == 2 else true_value_dict["val_Re"]
    return None


def get_all_moment_true_values(json_file: str) -> dict[str, float]:
    """Get all true moment values from a JSON file into a dict mapping 'H_alpha_L_M' -> value."""
    with open(json_file, 'r') as f:
        true_values = json.load(f)
    true_value_dict = {}
    for entry in true_values:
        alpha = entry["qn"]["momentIndex"]
        L = entry["qn"]["L"]
        M = entry["qn"]["M"]
        var_name = f"H_{alpha}_{L}_{M}"
        true_value_dict[var_name] = entry["val_Im"] if int(alpha) == 2 else entry["val_Re"]
    return true_value_dict


def find_last_mcmc_step(fpath: str) -> str:
    """Read backwards in chunks to find the last MCMC step progress indicator (e.g. '800/2000')."""
    pattern = re.compile(r'(\d+/\d+)')
    try:
        with open(fpath, 'rb') as f:
            f.seek(0, os.SEEK_END)
            pos = f.tell()
            remainder = b''
            while pos > 0:
                read_size = min(pos, 262144)
                pos -= read_size
                f.seek(pos)
                chunk = f.read(read_size) + remainder
                matches = pattern.findall(chunk.decode('ascii', errors='ignore'))
                if matches:
                    return matches[-1]
                remainder = chunk[:100]
    except Exception:
        pass
    return 'N/A'


def parse_single_log_status(fpath: str, active_timeout_min: float = 30.0) -> dict:
    """
    Parse a single logRooFit*.txt log file to determine its fit status.

    Returns dict with keys: bin, file, path, mass, target_bin, status, details, last_step, acceptance, mtime, size_bytes.
    """
    fname = os.path.basename(fpath)
    match_bin = re.search(r'logRooFit(\d+)\.txt', fname)
    bin_idx = int(match_bin.group(1)) if match_bin else -1

    if not os.path.exists(fpath):
        return {
            'bin': bin_idx, 'file': fname, 'path': fpath, 'mass': None, 'target_bin': None,
            'status': 'MISSING', 'details': 'File does not exist',
            'last_step': '-', 'acceptance': '-', 'mtime': '-', 'mtime_epoch': 0, 'size_bytes': 0
        }

    fsize = os.path.getsize(fpath)
    if fsize == 0:
        return {
            'bin': bin_idx, 'file': fname, 'path': fpath, 'mass': None, 'target_bin': None,
            'status': 'EMPTY', 'details': 'File is 0 bytes',
            'last_step': '-', 'acceptance': '-', 'mtime': '-', 'mtime_epoch': 0, 'size_bytes': 0
        }

    mtime = os.path.getmtime(fpath)
    time_diff_min = (time.time() - mtime) / 60.0
    mtime_str = datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')

    with open(fpath, 'r', errors='ignore') as f:
        if fsize > 300000:
            head = f.read(50000)
            f.seek(max(0, fsize - 250000))
            tail = f.read()
            sample_text = head + '\n' + tail
        else:
            sample_text = f.read()
            head = sample_text
            tail = sample_text

    mass_match = re.search(r'Mpi0eta([0-9.]+)', head)
    mass = float(mass_match.group(1)) if mass_match else None

    target_match = re.search(r'--> Target Bin\s*:\s*([A-Za-z0-9_.]+)', head)
    target_bin = target_match.group(1) if target_match else (f"Mpi0eta{mass:.6f}_" if mass is not None else None)

    is_failed = ('WARNING FitManager::Run: Fit failed' in tail or
                 'Skipping PlotDataModel() due to failed fit' in tail or
                 '1DStep Failed! Retrying (10/10)' in tail)

    if is_failed:
        reason = 'Fit failed (FitManager::Run returned kFALSE)'
        if '1DStep Failed! Retrying (10/10)' in sample_text:
            reason = '1DStep failed after 10 retries'
        return {
            'bin': bin_idx, 'file': fname, 'path': fpath, 'mass': mass, 'target_bin': target_bin,
            'status': 'FAILED', 'details': reason,
            'last_step': '-', 'acceptance': '-', 'mtime': mtime_str, 'mtime_epoch': mtime, 'size_bytes': fsize
        }

    is_success = ('PosteriorCovariance' in tail and 'Freeing MC Tree Memory' in tail) or \
                 ("--> Matrix successfully written to file as 'PosteriorCovariance'" in tail)

    if is_success:
        reason = 'Posterior covariance written successfully'
        return {
            'bin': bin_idx, 'file': fname, 'path': fpath, 'mass': mass, 'target_bin': target_bin,
            'status': 'SUCCESS', 'details': reason,
            'last_step': 'Done', 'acceptance': '-', 'mtime': mtime_str, 'mtime_epoch': mtime, 'size_bytes': fsize
        }

    last_step = find_last_mcmc_step(fpath)
    acc_match = re.findall(r'acceptance ([0-9.]+)', tail)
    last_acc = f'{float(acc_match[-1]):.3f}' if acc_match else 'N/A'

    if time_diff_min <= active_timeout_min:
        status = 'RUNNING'
        details = f'MCMC in progress (active {time_diff_min:.1f}m ago)'
    else:
        status = 'STALLED'
        details = f'Incomplete (modified {time_diff_min/60.0:.1f}h ago)'

    return {
        'bin': bin_idx, 'file': fname, 'path': fpath, 'mass': mass, 'target_bin': target_bin,
        'status': status, 'details': details,
        'last_step': last_step, 'acceptance': last_acc, 'mtime': mtime_str, 'mtime_epoch': mtime, 'size_bytes': fsize
    }


def check_projection_scaling(root_fpath: str, tolerance: float = 0.30) -> tuple[bool, float | None, str]:
    """
    Checks whether the RooPlot projection curve in ResultsBruMcmc*.root matches
    the data events normalization.
    
    When allow_negative_pdf=True and sampling wanders into negative-intensity
    phase space, BruComponentsPDF's quadratic penalty guard inflates the analytical
    normalization integral (denominator), while the 1D projection numerator is not
    inflated, causing the projection curve to be depressed by 8x-18x (Data/Curve >> 1.0).
    
    Returns:
        (is_valid: bool, ratio: float | None, message: str)
        where ratio = sum(data_points) / sum(curve_evals).
        A healthy physical fit has (1.0 - tolerance) <= ratio <= (1.0 + tolerance), typically ~1.00.
    """
    if ROOT is None or not os.path.exists(root_fpath):
        return True, None, "File missing or ROOT unavailable"

    try:
        f = ROOT.TFile.Open(root_fpath, "READ")
        if not f or f.IsZombie():
            if f:
                f.Close()
            return False, None, "File is zombie or cannot be opened"

        canvas = None
        for k in f.GetListOfKeys():
            name = k.GetName()
            if "__Phi" in name or "__Theta" in name or "phi_eta_gj" in name or "cosTheta_eta_gj" in name:
                obj = f.Get(name)
                if obj and (obj.InheritsFrom("TCanvas") or obj.InheritsFrom("RooPlot")):
                    canvas = obj
                    break

        if not canvas:
            f.Close()
            return True, None, "No projection canvas found in file"

        pad = canvas.GetPad(1) if hasattr(canvas, "GetPad") and canvas.GetPad(1) else canvas
        h_data = None
        curve = None

        if hasattr(pad, "GetListOfPrimitives"):
            for obj in pad.GetListOfPrimitives():
                if obj.InheritsFrom("RooHist") or "h_Data" in obj.GetName():
                    h_data = obj
                elif obj.InheritsFrom("RooCurve") and "paramBox" not in obj.GetName():
                    curve = obj
        elif hasattr(canvas, "getHist") and hasattr(canvas, "getCurve"):
            h_data = canvas.getHist()
            curve = canvas.getCurve()

        if not h_data or not curve:
            f.Close()
            return True, None, "Could not locate both data hist and curve in projection canvas"

        n_pts = h_data.GetN()
        if n_pts == 0:
            f.Close()
            return True, None, "Empty data histogram"

        s_data = sum(h_data.GetPointY(i) for i in range(n_pts))
        s_curve = sum(curve.Eval(h_data.GetPointX(i)) for i in range(n_pts))
        f.Close()

        if s_curve <= 0:
            return False, 999.0, "Projection curve integral is zero or negative"

        ratio = s_data / s_curve
        min_allowed = 1.0 - tolerance
        max_allowed = 1.0 + tolerance

        if ratio < min_allowed or ratio > max_allowed:
            msg = (
                f"Negative-PDF scaling distortion detected! Data/Curve ratio is {ratio:.2f} "
                f"(expected ~1.00, allowed [{min_allowed:.2f}, {max_allowed:.2f}]). "
                f"The negative PDF penalty guard inflated the normalization integral, "
                f"indicating walker trapped in unphysical negative-intensity phase space."
            )
            return False, ratio, msg

        return True, ratio, f"Normal projection curve scaling (Data/Curve = {ratio:.2f})"
    except Exception as e:
        return True, None, f"Error inspecting curve scaling: {e}"


def compute_batch_chain_diagnostics(
    all_data: dict,
    col_names: list[str],
    burn_in: int = 2000,
    max_lag: int | None = None
) -> dict | None:
    """
    Vectorized calculation of integrated autocorrelation time (tau_int) and ESS across all H_* moments
    using a single index expansion and batched 2D NumPy real FFTs (Wiener-Khinchin theorem).
    Matches BruFit's AutocorrPlot.cpp definition.
    """
    if not col_names:
        return None

    has_wt = "weight" in all_data
    n_tree = len(all_data[col_names[0]])
    if n_tree == 0:
        return None

    if has_wt:
        w_arr = np.asarray(all_data["weight"], dtype=np.float64)
        int_w = np.maximum(1, np.round(w_arr)).astype(np.int64)
        n_raw = int(np.sum(int_w))
    else:
        int_w = None
        n_raw = int(n_tree)

    if n_raw <= 50:
        return None

    eff_burn_in = burn_in
    if eff_burn_in >= n_raw - 50:
        eff_burn_in = min(eff_burn_in, n_raw // 2)

    n_eff = n_raw - eff_burn_in
    if n_eff <= 1:
        return None

    if int_w is not None:
        cum_w = np.cumsum(int_w)
        eff_burn_entry = int(np.searchsorted(cum_w, eff_burn_in))
        emp_acc = float(n_tree - eff_burn_entry) / float(max(1, n_eff))
    else:
        cum_w = None
        arr0 = np.asarray(all_data[col_names[0]], dtype=np.float64)
        emp_acc = float(np.mean(arr0[eff_burn_in + 1:] != arr0[eff_burn_in:-1])) if len(arr0) > eff_burn_in + 1 else None

    # Bound FFT array length and column chunk size to prevent OOM on ultra-low-acceptance chains
    max_fft_samples = 250_000
    stride = 1 if n_eff <= max_fft_samples else ((n_eff + max_fft_samples - 1) // max_fft_samples)
    if stride == 1:
        if int_w is not None:
            idx = np.repeat(np.arange(n_tree, dtype=np.int64), int_w)[eff_burn_in:]
        else:
            idx = slice(eff_burn_in, n_raw)
        n_grid = n_eff
    else:
        sample_pos = np.arange(eff_burn_in, n_raw, stride, dtype=np.int64)
        n_grid = len(sample_pos)
        if cum_w is not None:
            idx = np.minimum(np.searchsorted(cum_w, sample_pos, side="right"), n_tree - 1)
        else:
            idx = sample_pos

    n_cols = len(col_names)
    grid_max_lag = (n_grid // 2) if max_lag is None else min(max(1, max_lag // stride), n_grid // 2)
    limit = min(n_grid, grid_max_lag + 1)
    n_fft = 1 << int(np.ceil(np.log2(max(2, 2 * n_grid - 1))))

    taus = np.empty(n_cols, dtype=np.float64)
    esses = np.empty(n_cols, dtype=np.float64)
    per_moment = {}
    chunk_size = 16

    for start in range(0, n_cols, chunk_size):
        chunk_cols = col_names[start:start + chunk_size]
        c_len = len(chunk_cols)
        mat = np.empty((c_len, n_grid), dtype=np.float64)
        for j, name in enumerate(chunk_cols):
            mat[j] = np.asarray(all_data[name], dtype=np.float64)[idx]

        means = mat.mean(axis=1, keepdims=True)
        xm = mat - means
        par_d = np.einsum("ij,ij->i", xm, xm)

        fx = np.fft.rfft(xm, n=n_fft, axis=1)
        rx = np.fft.irfft(fx * np.conj(fx), n=n_fft, axis=1)[:, :limit]

        for j, name in enumerate(chunk_cols):
            pd = par_d[j]
            if pd <= 0.0:
                tau = 0.5
            else:
                acf_row = rx[j] / pd
                neg = np.flatnonzero(acf_row[1:] <= 0.0)
                stop = (1 + int(neg[0])) if len(neg) > 0 else limit
                tau_grid = 0.5 + float(np.sum(acf_row[1:stop]))
                tau = tau_grid if stride == 1 else max(0.5, float(stride) * tau_grid)
            ess = (n_eff / (2.0 * tau)) if tau > 0.0 else float(n_eff)
            ess = min(float(n_eff), ess)
            taus[start + j] = tau
            esses[start + j] = ess
            per_moment[name] = {'tau': float(tau), 'ess': float(ess)}

    return {
        'n_raw': int(n_raw),
        'burn_in': int(eff_burn_in),
        'n_eff': int(n_eff),
        'n_samples': int(n_eff),
        'tau_mean': float(np.mean(taus)),
        'tau_max': float(np.max(taus)),
        'tau_min': float(np.min(taus)),
        'ess_min': float(np.min(esses)),
        'ess_mean': float(np.mean(esses)),
        'ess_max': float(np.max(esses)),
        'acceptance': emp_acc,
        'per_moment': per_moment,
    }


def compute_chain_diagnostics_from_root(root_fpath: str, burn_in: int = 2000, tree_name: str = "MCMCTree") -> dict | None:
    """
    Load an MCMCTree from a ROOT file and calculate integrated autocorrelation time and ESS.
    Uses uproot when available for fast direct branch extraction, falling back to ROOT.RDataFrame.
    """
    if not root_fpath or not os.path.exists(root_fpath):
        return None

    if uproot is not None:
        try:
            with uproot.open(root_fpath) as uf:
                if tree_name in uf:
                    tree = uf[tree_name]
                    keys = [str(k) for k in tree.keys()]
                    all_h = [k for k in keys if k.startswith("H_") and not k.endswith("_err")]
                    strict_h = [k for k in all_h if len(k.split("_")) == 4]
                    col_names = strict_h if strict_h else all_h
                    if not col_names:
                        return None
                    cols_to_load = list(col_names)
                    if "weight" in keys:
                        cols_to_load.append("weight")
                    all_data = tree.arrays(cols_to_load, library="np")
                    return compute_batch_chain_diagnostics(all_data, col_names, burn_in=burn_in)
        except Exception:
            pass

    if ROOT is None:
        return None
    try:
        rdf = ROOT.RDataFrame(tree_name, root_fpath)
        all_cols = [str(c) for c in rdf.GetColumnNames()]
        all_h = [c for c in all_cols if c.startswith("H_") and not c.endswith("_err")]
        strict_h = [c for c in all_h if len(c.split("_")) == 4]
        col_names = strict_h if strict_h else all_h
        if not col_names:
            return None

        has_wt = "weight" in all_cols
        cols_to_load = list(col_names)
        if has_wt:
            cols_to_load.append("weight")

        all_data = rdf.AsNumpy(cols_to_load)
        return compute_batch_chain_diagnostics(all_data, col_names, burn_in=burn_in)
    except Exception:
        return None


def _resolve_bin_result_root(log_dir: str, bin_record: dict) -> str | None:
    """Find the MCMC result ROOT file path for a mass bin record."""
    target_bin = bin_record.get('target_bin')
    mass = bin_record.get('mass')
    candidate_dirs = []
    if target_bin:
        candidate_dirs.append(os.path.join(log_dir, target_bin))
    if mass is not None:
        candidate_dirs.append(os.path.join(log_dir, f"Mpi0eta{mass:.6f}_"))
        candidate_dirs.append(os.path.join(log_dir, f"Mpi0eta{mass:.6f}"))
        candidate_dirs.append(os.path.join(log_dir, f"Mpi0eta{mass:.4f}"))

    root_candidates = [
        "ResultsBruMcmcFixedCovariance.root",
        "ResultsBruMcmcCovariance.root",
        "ResultsBruMcmc.root",
    ]
    for cdir in candidate_dirs:
        for rname in root_candidates:
            p = os.path.join(cdir, rname)
            if os.path.isfile(p):
                return p
    return None


def _check_mcmc_root_valid(root_fpath: str, tree_name: str = "MCMCTree") -> tuple[bool, int]:
    """Return (has_valid_tree, n_entries) for an MCMC ROOT output file."""
    if not os.path.isfile(root_fpath) or os.path.getsize(root_fpath) == 0:
        return False, 0
    try:
        with open(root_fpath, "rb") as fb:
            if fb.read(4) != b"root":
                return False, 0
    except OSError:
        return False, 0

    if uproot is not None:
        try:
            with uproot.open(root_fpath) as uf:
                if len(uf.keys()) == 0 or tree_name not in uf:
                    return False, 0
                n_entries = int(uf[tree_name].num_entries)
                return (n_entries > 0), n_entries
        except Exception:
            return False, 0

    if ROOT is not None:
        try:
            f = ROOT.TFile.Open(root_fpath, "READ")
            if not f or f.IsZombie() or f.GetNkeys() == 0 or f.TestBit(ROOT.TFile.kRecovered):
                if f:
                    f.Close()
                return False, 0
            tree = f.Get(tree_name)
            if not tree:
                f.Close()
                return False, 0
            n_entries = int(tree.GetEntries())
            f.Close()
            return (n_entries > 0), n_entries
        except Exception:
            return False, 0

    return True, 1


def get_fit_status(
    log_files_or_dir: str | list[str],
    active_timeout_min: float = 30.0,
    as_dataframe: bool = False,
    mcmc: MomentsMCMC | None = None,
    burn_in: int = 2000,
    compute_diagnostics: bool = True
) -> dict:
    """
    Extract and list fit status across all log files matching a path or glob pattern.
    For successful fits, reports integrated autocorrelation time (tau_int) and ESS.
    """
    if isinstance(log_files_or_dir, (list, tuple)):
        file_list = list(log_files_or_dir)
    elif os.path.isdir(str(log_files_or_dir)):
        file_list = glob.glob(os.path.join(str(log_files_or_dir), 'logRooFit*.txt'))
    else:
        file_list = glob.glob(str(log_files_or_dir))

    def _sort_key(f):
        m = re.search(r'logRooFit(\d+)\.txt', os.path.basename(f))
        return int(m.group(1)) if m else f

    file_list.sort(key=_sort_key)
    bin_results = [parse_single_log_status(f, active_timeout_min=active_timeout_min) for f in file_list]

    for r in bin_results:
        r['tau_mean'] = None
        r['tau_max'] = None
        r['ess_min'] = None
        r['ess_mean'] = None
        r['n_eff'] = None

        found_root = _resolve_bin_result_root(os.path.dirname(r['path']), r)
        if found_root and r['status'] in ('RUNNING', 'STALLED', 'MISSING', 'EMPTY'):
            has_tree, n_entries = _check_mcmc_root_valid(found_root)
            if has_tree and n_entries > 0:
                r['status'] = 'SUCCESS'
                r['last_step'] = 'Done'
                r['details'] = f'MCMCTree verified ({n_entries} entries)'

        if compute_diagnostics and r['status'] in ('SUCCESS', 'RUNNING', 'STALLED'):
            if mcmc is not None:
                b_idx = r['bin']
                if b_idx is not None and 0 <= b_idx < mcmc.get_N_rdf() and mcmc.get_rdf(b_idx) is not None:
                    try:
                        diag = mcmc.get_mcmc_diagnostics(b_idx, burn_in=burn_in)
                        r['tau_mean'] = diag['mean_tau']
                        r['tau_max'] = diag['max_tau']
                        r['ess_min'] = diag['min_ess']
                        r['ess_mean'] = diag['mean_ess']
                        r['n_eff'] = diag['n_samples']
                        r['diagnostics'] = diag
                    except Exception:
                        pass
            elif found_root:
                if ROOT is not None:
                    scale_ok, scale_ratio, scale_msg = check_projection_scaling(found_root)
                    r['scale_ratio'] = scale_ratio
                    r['scale_ok'] = scale_ok
                    if not scale_ok:
                        r['status'] = 'SCALE_MISMATCH'
                        r['details'] = f"Negative-PDF scaling mismatch (ratio {scale_ratio:.2f} != 1.0)"
                if r['status'] == 'SUCCESS':
                    diag = compute_chain_diagnostics_from_root(found_root, burn_in=burn_in)
                    if diag:
                        r['tau_mean'] = diag['tau_mean']
                        r['tau_max'] = diag['tau_max']
                        r['ess_min'] = diag['ess_min']
                        r['ess_mean'] = diag['ess_mean']
                        r['n_eff'] = diag['n_eff']
                        r['diagnostics'] = diag

    summary = {}
    for r in bin_results:
        summary[r['status']] = summary.get(r['status'], 0) + 1

    failed_bins = [r['bin'] for r in bin_results if r['status'] == 'FAILED']
    scale_mismatch_bins = [r['bin'] for r in bin_results if r['status'] == 'SCALE_MISMATCH']
    success_bins = [r['bin'] for r in bin_results if r['status'] == 'SUCCESS']
    running_bins = [r['bin'] for r in bin_results if r['status'] == 'RUNNING']
    stalled_bins = [r['bin'] for r in bin_results if r['status'] == 'STALLED']
    missing_bins = [r['bin'] for r in bin_results if r['status'] in ('MISSING', 'EMPTY')]

    result_dict = {
        'summary': summary,
        'bins': bin_results,
        'failed_bins': failed_bins,
        'scale_mismatch_bins': scale_mismatch_bins,
        'success_bins': success_bins,
        'running_bins': running_bins,
        'stalled_bins': stalled_bins,
        'missing_bins': missing_bins,
    }

    if as_dataframe:
        try:
            import pandas as pd
            result_dict['df'] = pd.DataFrame(bin_results)
        except ImportError:
            pass

    return result_dict


def print_fit_status(
    log_files_or_dir: str | list[str],
    filter_status: str = None,
    active_timeout_min: float = 30.0,
    mcmc: MomentsMCMC | None = None,
    burn_in: int = 2000,
    show_diagnostics: bool = True
) -> dict:
    """Print a formatted summary table of fit statuses from log files."""
    status_info = get_fit_status(
        log_files_or_dir,
        active_timeout_min=active_timeout_min,
        mcmc=mcmc,
        burn_in=burn_in,
        compute_diagnostics=show_diagnostics
    )
    bin_results = status_info['bins']

    if filter_status is not None:
        filter_upper = filter_status.upper()
        bin_results = [r for r in bin_results if r['status'] == filter_upper]

    header_title = f"Fit Status Summary ({log_files_or_dir})"
    print("=" * 118)
    print(header_title)
    print("=" * 118)
    print(f"{'Bin':<5} {'Mass (GeV)':<12} {'Status':<10} {'Tau (mean/max)':<18} {'ESS (min/mean)':<18} {'Last Step':<12} {'Acceptance':<12} {'Details'}")
    print("-" * 118)

    for r in bin_results:
        m_str = f"{r['mass']:.4f}" if r['mass'] is not None else "N/A"
        tau_str = f"{r['tau_mean']:.1f} / {r['tau_max']:.1f}" if r.get('tau_mean') is not None else "-"
        ess_str = f"{r['ess_min']:.1f} / {r['ess_mean']:.1f}" if r.get('ess_min') is not None else "-"
        print(f"{r['bin']:<5} {m_str:<12} {r['status']:<10} {tau_str:<18} {ess_str:<18} {r['last_step']:<12} {r['acceptance']:<12} {r['details']}")

    print("=" * 118)
    summary_str = " | ".join(f"{k}: {v}" for k, v in sorted(status_info['summary'].items()))
    print(f"Total Bins: {len(status_info['bins'])} | {summary_str}")
    if status_info.get('scale_mismatch_bins'):
        print(f"⚠️  Negative-PDF Scaling Mismatch Bins ({len(status_info['scale_mismatch_bins'])}): {status_info['scale_mismatch_bins']}")
    if status_info['failed_bins']:
        print(f"Failed Bins ({len(status_info['failed_bins'])}): {status_info['failed_bins']}")
    if status_info['running_bins']:
        print(f"Running Bins ({len(status_info['running_bins'])}): {status_info['running_bins']}")
    if status_info['stalled_bins']:
        print(f"Stalled Bins ({len(status_info['stalled_bins'])}): {status_info['stalled_bins']}")
    if status_info['success_bins']:
        print(f"Successful Bins ({len(status_info['success_bins'])}): {status_info['success_bins']}")

    success_with_diag = [r for r in bin_results if r.get('tau_mean') is not None]
    if success_with_diag:
        mean_taus = [r['tau_mean'] for r in success_with_diag]
        max_taus = [r['tau_max'] for r in success_with_diag]
        min_esses = [r['ess_min'] for r in success_with_diag]
        mean_esses = [r['ess_mean'] for r in success_with_diag]
        print("-" * 118)
        print(f"MCMC Diagnostics for Successful Bins ({len(success_with_diag)} bins evaluated, configured burn_in={burn_in}):")
        print(f"  Integrated Autocorrelation Time (tau_int) -> Average mean: {np.mean(mean_taus):.2f} | Peak slowest: {np.max(max_taus):.2f}")
        print(f"  Effective Sample Size (ESS)              -> Average mean: {np.mean(mean_esses):.1f} | Worst-case min: {np.min(min_esses):.1f}")
    print("=" * 118)

    return status_info


# Alias for backward compatibility
list_fit_status = print_fit_status

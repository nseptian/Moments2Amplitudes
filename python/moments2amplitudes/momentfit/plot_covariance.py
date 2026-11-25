"""
Plotting utilities for MCMC Proposal Covariance Matrices.

Compares covariance and correlation matrices before alignment (raw linear algebra moments)
and after alignment (Symmetric Positive-Definite regularized covariance for BruFit).

Visualizations generated:
  1. Per-mass-bin 4-panel comparison:
     - Panel 1: Correlation Matrix Before Alignment (Raw Moments Covariance)
     - Panel 2: Correlation Matrix After Alignment (BruFit SPD Regularized)
     - Panel 3: Correlation Difference (Delta R = R_aligned - R_raw)
     - Panel 4: Eigenvalue Spectrum (lambda_k Before vs After Alignment, symlog scale)
  2. Multi-bin summary across mass:
     - Minimum eigenvalues before vs after alignment
     - Count of regularized/clamped modes vs mass
     - Condition number and max correlation shift vs mass
"""

from __future__ import annotations

import functools
import glob
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

print = functools.partial(print, flush=True)


def _load_matrix_uproot(file_path: str) -> Tuple[np.ndarray, List[str]]:
    """Loads MyCovMatrix and parameter names from a ROOT file using uproot."""
    import uproot

    with uproot.open(file_path) as root_f:
        if "MyCovMatrix;1" in root_f:
            cov_obj = root_f["MyCovMatrix;1"]
        elif "MyCovMatrix" in root_f:
            cov_obj = root_f["MyCovMatrix"]
        else:
            raise KeyError(f"'MyCovMatrix' not found in {file_path}. Available keys: {root_f.keys()}")

        elems = np.asarray(cov_obj.members["fElements"], dtype=np.float64)
        nrows = int(cov_obj.members["fNrows"])

        mat = np.zeros((nrows, nrows), dtype=np.float64)
        idx = 0
        for i in range(nrows):
            row_len = nrows - i
            mat[i, i:] = elems[idx : idx + row_len]
            mat[i:, i] = elems[idx : idx + row_len]
            idx += row_len

        # Extract parameter names
        names: List[str] = []
        if "MyCovMatrix_Names;1" in root_f or "MyCovMatrix_Names" in root_f:
            names_obj = root_f.get("MyCovMatrix_Names;1") or root_f.get("MyCovMatrix_Names")
            names = [str(x) for x in names_obj]
        elif "paramNames;1" in root_f or "paramNames" in root_f:
            names_obj = root_f.get("paramNames;1") or root_f.get("paramNames")
            names_str = str(names_obj.string() if hasattr(names_obj, "string") else names_obj)
            names = [x.strip() for x in names_str.split(",") if x.strip()]
        else:
            names = [f"Param_{i}" for i in range(nrows)]

    return mat, names


def _load_matrix_pyroot(file_path: str) -> Tuple[np.ndarray, List[str]]:
    """Fallback: loads MyCovMatrix and parameter names using PyROOT."""
    import ROOT

    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        raise IOError(f"Could not open ROOT file: {file_path}")

    m_obj = f.Get("MyCovMatrix")
    if not m_obj:
        f.Close()
        raise KeyError(f"'MyCovMatrix' not found in {file_path}")

    nrows = m_obj.GetNrows()
    mat = np.zeros((nrows, nrows), dtype=np.float64)
    for i in range(nrows):
        for j in range(nrows):
            mat[i, j] = m_obj[i][j]

    names: List[str] = []
    names_list = f.Get("MyCovMatrix_Names")
    if names_list:
        for i in range(names_list.GetSize()):
            names.append(str(names_list.At(i).GetString()))
    else:
        pnames_str = f.Get("paramNames")
        if pnames_str:
            names = [x.strip() for x in str(pnames_str.GetString()).split(",") if x.strip()]
        else:
            names = [f"Param_{i}" for i in range(nrows)]

    f.Close()
    return mat, names


def load_covariance_matrix(file_path: str) -> Tuple[np.ndarray, List[str]]:
    """Loads covariance matrix and parameter names using uproot (with PyROOT fallback)."""
    try:
        return _load_matrix_uproot(file_path)
    except Exception as e_uproot:
        try:
            return _load_matrix_pyroot(file_path)
        except Exception as e_pyroot:
            raise RuntimeError(
                f"Failed to load covariance matrix from {file_path}.\n"
                f"Uproot error: {e_uproot}\n"
                f"PyROOT error: {e_pyroot}"
            )


def cov_to_correlation(cov: np.ndarray, eps: float = 1e-12) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes standard correlation matrix R_ij = Cov_ij / (sigma_i * sigma_j).
    Returns (corr_matrix, sigmas).
    """
    diag = np.diag(cov)
    sigma = np.sqrt(np.maximum(np.abs(diag), eps))
    outer = np.outer(sigma, sigma)
    corr = cov / outer
    np.fill_diagonal(corr, 1.0)
    return np.clip(corr, -1.0, 1.0), sigma


def check_positive_definite(mat: np.ndarray) -> Tuple[bool, float, int]:
    """
    Checks if symmetric matrix is positive-definite via eigenvalues and Cholesky.
    Returns (is_spd, min_eigenval, num_non_positive_eigenvals).
    """
    evals = np.linalg.eigvalsh(mat)
    min_ev = float(np.min(evals))
    num_non_pos = int(np.sum(evals <= 1e-12))
    try:
        _ = np.linalg.cholesky(mat)
        is_spd = True
    except np.linalg.LinAlgError:
        is_spd = False
    return is_spd, min_ev, num_non_pos


def _find_moment_block_boundaries(param_names: List[str]) -> List[Tuple[str, int, int]]:
    """
    Identifies index ranges for H_0, H_1, H_2, and Yield parameters.
    Returns list of (label, start_idx, end_idx).
    """
    blocks = []
    n = len(param_names)
    h0_idx = [i for i, name in enumerate(param_names) if name.startswith("H_0")]
    h1_idx = [i for i, name in enumerate(param_names) if name.startswith("H_1")]
    h2_idx = [i for i, name in enumerate(param_names) if name.startswith("H_2")]
    yld_idx = [i for i, name in enumerate(param_names) if "Yld" in name or "Yield" in name]

    if h0_idx:
        blocks.append(("H_0", h0_idx[0], h0_idx[-1]))
    if h1_idx:
        blocks.append(("H_1", h1_idx[0], h1_idx[-1]))
    if h2_idx:
        blocks.append(("H_2", h2_idx[0], h2_idx[-1]))
    if yld_idx:
        blocks.append(("Yield", yld_idx[0], yld_idx[-1]))

    return blocks


def plot_covariance_comparison_single_bin(
    raw_path: str,
    aligned_path: str,
    output_dir: str,
    mass_bin_label: str = "",
    tbin: str = "",
    max_l: Union[int, str] = "",
    formats: Optional[List[str]] = None,
    min_eigen_floor: float = 1e-4,
) -> Optional[Dict[str, Any]]:
    """
    Generates a 4-panel comparison figure for a single mass bin:
      1. Raw correlation matrix (before alignment)
      2. Regularized SPD correlation matrix (after alignment)
      3. Correlation difference (Delta R)
      4. Eigenvalue spectrum (lambda_k before vs after on symlog scale)
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if formats is None:
        formats = ["png", "pdf"]

    cov_raw, names_raw = load_covariance_matrix(raw_path)
    cov_aligned, names_aligned = load_covariance_matrix(aligned_path)

    param_names = names_aligned if names_aligned else names_raw
    n = len(cov_raw)

    r_raw, sigma_raw = cov_to_correlation(cov_raw)
    r_aligned, sigma_aligned = cov_to_correlation(cov_aligned)
    r_diff = r_aligned - r_raw

    # Eigenvalues of correlation matrices
    evals_raw = np.linalg.eigvalsh(r_raw)
    evals_aligned = np.linalg.eigvalsh(r_aligned)

    is_spd_raw, min_ev_raw, n_neg_raw = check_positive_definite(r_raw)
    is_spd_aligned, min_ev_aligned, n_neg_aligned = check_positive_definite(r_aligned)

    # Also check full covariance eigenvalues
    cov_evals_raw = np.linalg.eigvalsh(cov_raw)
    cov_evals_aligned = np.linalg.eigvalsh(cov_aligned)
    cov_spd_raw, cov_min_ev_raw, _ = check_positive_definite(cov_raw)
    cov_spd_aligned, cov_min_ev_aligned, _ = check_positive_definite(cov_aligned)

    # Block boundaries for grid lines
    blocks = _find_moment_block_boundaries(param_names)
    dividers = [b[2] + 0.5 for b in blocks[:-1]] if len(blocks) > 1 else []

    fig, axes = plt.subplots(2, 2, figsize=(15, 13))
    tbin_str = f" [{tbin}]" if tbin else ""
    l_str = f" (maxL={max_l})" if max_l else ""
    mass_str = f" Mass Bin: {mass_bin_label}" if mass_bin_label else ""

    fig.suptitle(
        f"MCMC Proposal Covariance Matrix Alignment: Before vs After SPD Regularization\n"
        f"{tbin_str}{l_str}{mass_str}",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )

    # --------------------------------------------------------------------------
    # Panel 1: Correlation Matrix Before Alignment
    # --------------------------------------------------------------------------
    ax0 = axes[0, 0]
    im0 = ax0.imshow(r_raw, cmap="coolwarm", vmin=-1.0, vmax=1.0, aspect="equal")
    status_raw_str = "Cholesky: FAILED (Indefinite)" if not is_spd_raw else "Cholesky: PASSED"
    ax0.set_title(
        f"Before Alignment: Raw Linear Algebra Moments\n"
        f"Dim: {n}×{n} | Non-pos evals: {n_neg_raw} | {status_raw_str}",
        fontsize=11,
        fontweight="semibold",
        color="darkred" if not is_spd_raw else "black",
    )
    ax0.set_xlabel("Parameter Index", fontsize=10)
    ax0.set_ylabel("Parameter Index", fontsize=10)
    for div in dividers:
        ax0.axhline(div, color="black", linestyle="--", linewidth=0.7, alpha=0.7)
        ax0.axvline(div, color="black", linestyle="--", linewidth=0.7, alpha=0.7)
    for b_name, s, e in blocks:
        mid = (s + e) / 2.0
        ax0.text(mid, -2.5, b_name, ha="center", va="bottom", fontsize=8, fontweight="bold", color="navy")
    plt.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04, label="Correlation $R_{ij}$")

    # --------------------------------------------------------------------------
    # Panel 2: Correlation Matrix After Alignment
    # --------------------------------------------------------------------------
    ax1 = axes[0, 1]
    im1 = ax1.imshow(r_aligned, cmap="coolwarm", vmin=-1.0, vmax=1.0, aspect="equal")
    status_aligned_str = "Cholesky: PASSED (Strictly SPD)" if is_spd_aligned else "Cholesky: FAILED"
    ax1.set_title(
        f"After Alignment: BruFit SPD Regularized\n"
        f"Dim: {n}×{n} | Non-pos evals: {n_neg_aligned} | {status_aligned_str}",
        fontsize=11,
        fontweight="semibold",
        color="darkgreen" if is_spd_aligned else "darkred",
    )
    ax1.set_xlabel("Parameter Index", fontsize=10)
    ax1.set_ylabel("Parameter Index", fontsize=10)
    for div in dividers:
        ax1.axhline(div, color="black", linestyle="--", linewidth=0.7, alpha=0.7)
        ax1.axvline(div, color="black", linestyle="--", linewidth=0.7, alpha=0.7)
    for b_name, s, e in blocks:
        mid = (s + e) / 2.0
        ax1.text(mid, -2.5, b_name, ha="center", va="bottom", fontsize=8, fontweight="bold", color="navy")
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04, label="Correlation $R_{ij}$")

    # --------------------------------------------------------------------------
    # Panel 3: Correlation Difference (Delta R = After - Before)
    # --------------------------------------------------------------------------
    ax2 = axes[1, 0]
    max_d = max(0.01, float(np.max(np.abs(r_diff))))
    im2 = ax2.imshow(r_diff, cmap="seismic", vmin=-max_d, vmax=max_d, aspect="equal")
    mean_abs_shift = float(np.mean(np.abs(r_diff)))
    ax2.set_title(
        fr"Correlation Difference ($\Delta R = R_{{aligned}} - R_{{raw}}$)" + "\n"
        fr"Max |$\Delta R$|: {max_d:.3e} | Mean |$\Delta R$|: {mean_abs_shift:.3e}",
        fontsize=11,
        fontweight="semibold",
    )
    ax2.set_xlabel("Parameter Index", fontsize=10)
    ax2.set_ylabel("Parameter Index", fontsize=10)
    for div in dividers:
        ax2.axhline(div, color="black", linestyle="--", linewidth=0.7, alpha=0.7)
        ax2.axvline(div, color="black", linestyle="--", linewidth=0.7, alpha=0.7)
    for b_name, s, e in blocks:
        mid = (s + e) / 2.0
        ax2.text(mid, -2.5, b_name, ha="center", va="bottom", fontsize=8, fontweight="bold", color="navy")
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04, label=r"Shift $\Delta R_{ij}$")

    # --------------------------------------------------------------------------
    # Panel 4: Eigenvalue Spectrum (symlog scale)
    # --------------------------------------------------------------------------
    ax3 = axes[1, 1]
    k_idx = np.arange(1, n + 1)

    ax3.plot(
        k_idx,
        evals_raw,
        "o-",
        color="crimson",
        markersize=3.5,
        linewidth=1.2,
        label=fr"Raw $R$ ($\lambda_{{min}} = {min_ev_raw:.2e}$, $N_{{neg}} = {n_neg_raw}$)",
    )
    ax3.plot(
        k_idx,
        evals_aligned,
        "s-",
        color="royalblue",
        markersize=3.5,
        linewidth=1.2,
        label=fr"Aligned $R$ ($\lambda_{{min}} = {min_ev_aligned:.2e}$, $N_{{neg}} = 0$)",
    )

    ax3.axhline(min_eigen_floor, color="darkgreen", linestyle="--", linewidth=1.5, label=f"Floor ({min_eigen_floor:.0e})")
    ax3.axhline(0.0, color="gray", linestyle=":", linewidth=1.0)

    ax3.set_yscale("symlog", linthresh=min_eigen_floor)
    ax3.set_xlabel(r"Eigenvalue Index $k$ (sorted ascending)", fontsize=10)
    ax3.set_ylabel(r"Eigenvalue $\lambda_k$ (symlog scale)", fontsize=10)
    ax3.set_title(
        "Correlation Matrix Spectral Regularization\n"
        fr"Projected onto SPD cone with floor $\lambda \geq {min_eigen_floor:.0e}$",
        fontsize=11,
        fontweight="semibold",
    )
    ax3.grid(True, which="both", linestyle="--", alpha=0.4)
    ax3.legend(loc="lower right", fontsize=9, framealpha=0.9)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    os.makedirs(output_dir, exist_ok=True)
    clean_lbl = mass_bin_label.replace("#", "").replace("{", "").replace("}", "").replace("/", "_")
    base_name = f"covariance_alignment_{clean_lbl}" if clean_lbl else "covariance_alignment"

    saved_paths = []
    for fmt in formats:
        out_p = os.path.join(output_dir, f"{base_name}.{fmt}")
        fig.savefig(out_p, dpi=200 if fmt == "png" else None)
        saved_paths.append(out_p)

    plt.close(fig)

    return {
        "mass_label": mass_bin_label,
        "n_params": n,
        "raw_is_spd": is_spd_raw,
        "raw_min_eigenval": min_ev_raw,
        "raw_n_neg_eigenvals": n_neg_raw,
        "aligned_is_spd": is_spd_aligned,
        "aligned_min_eigenval": min_ev_aligned,
        "aligned_n_neg_eigenvals": n_neg_aligned,
        "max_corr_diff": max_d,
        "mean_abs_corr_diff": mean_abs_shift,
        "saved_plots": saved_paths,
    }


def plot_covariance_summary_across_mass(
    bin_results: List[Dict[str, Any]],
    output_dir: str,
    tbin: str = "",
    max_l: Union[int, str] = "",
    formats: Optional[List[str]] = None,
    min_eigen_floor: float = 1e-4,
) -> Optional[List[str]]:
    """
    Generates a multi-panel summary figure tracking regularization metrics across all mass bins.
    """
    if not bin_results:
        return None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if formats is None:
        formats = ["png", "pdf"]

    mass_centers = []
    for r in bin_results:
        lbl = r["mass_label"]
        # Match pattern e.g. Mpi0eta_1.300 or 1.300
        match = re.search(r"(\d+\.\d+)", lbl)
        if match:
            mass_centers.append(float(match.group(1)))
        else:
            mass_centers.append(len(mass_centers) + 1.0)

    order = np.argsort(mass_centers)
    masses = np.array(mass_centers)[order]
    results_ordered = [bin_results[i] for i in order]

    raw_min_ev = np.array([r["raw_min_eigenval"] for r in results_ordered])
    aligned_min_ev = np.array([r["aligned_min_eigenval"] for r in results_ordered])
    n_neg = np.array([r["raw_n_neg_eigenvals"] for r in results_ordered])
    max_d = np.array([r["max_corr_diff"] for r in results_ordered])
    mean_d = np.array([r["mean_abs_corr_diff"] for r in results_ordered])

    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    tbin_str = f" [{tbin}]" if tbin else ""
    l_str = f" (maxL={max_l})" if max_l else ""

    fig.suptitle(
        f"Covariance Matrix SPD Regularization Summary vs Mass\n{tbin_str}{l_str}",
        fontsize=13,
        fontweight="bold",
        y=0.98,
    )

    # Subplot 1: Minimum Eigenvalue vs Mass
    ax0 = axes[0]
    ax0.plot(masses, raw_min_ev, "o-", color="crimson", linewidth=1.5, markersize=5, label=r"Raw $R$ $\lambda_{min}$ (Before Alignment)")
    ax0.plot(masses, aligned_min_ev, "s-", color="royalblue", linewidth=1.5, markersize=5, label=r"Aligned $R$ $\lambda_{min}$ (After Alignment)")
    ax0.axhline(min_eigen_floor, color="darkgreen", linestyle="--", linewidth=1.5, label=f"Floor ({min_eigen_floor:.0e})")
    ax0.axhline(0.0, color="gray", linestyle=":", linewidth=1.0)
    ax0.set_yscale("symlog", linthresh=min_eigen_floor)
    ax0.set_ylabel(r"Min Eigenvalue $\lambda_{min}$", fontsize=10)
    ax0.set_title("Minimum Eigenvalue of Correlation Matrix", fontsize=11, fontweight="semibold")
    ax0.grid(True, which="both", linestyle="--", alpha=0.4)
    ax0.legend(loc="upper right", fontsize=9, framealpha=0.9)

    # Subplot 2: Number of Negative / Clamped Eigenvalues
    ax1 = axes[1]
    ax1.bar(masses, n_neg, width=(masses[1] - masses[0]) * 0.8 if len(masses) > 1 else 0.03, color="coral", edgecolor="darkred", alpha=0.85, label=r"Raw Non-Positive Eigenvalues ($N_{\leq 0}$)")
    ax1.set_ylabel("Count of Non-Positive Modes", fontsize=10)
    ax1.set_title("Number of Negative / Indefinite Modes Fixed by SPD Regularization", fontsize=11, fontweight="semibold")
    ax1.grid(True, linestyle="--", alpha=0.4)
    ax1.legend(loc="upper right", fontsize=9, framealpha=0.9)

    # Subplot 3: Correlation Matrix Shifts
    ax2 = axes[2]
    ax2.plot(masses, max_d, "d-", color="darkorange", linewidth=1.5, markersize=5, label=r"Max Absolute Shift $\max(|\Delta R_{ij}|)$")
    ax2.plot(masses, mean_d, "^-", color="purple", linewidth=1.5, markersize=5, label=r"Mean Absolute Shift $\mathrm{avg}(|\Delta R_{ij}|)$")
    ax2.set_xlabel(r"Invariant Mass $m(\pi^0\eta)$ [GeV/$c^2$]", fontsize=10)
    ax2.set_ylabel(r"Correlation Shift $|\Delta R|$", fontsize=10)
    ax2.set_title("Magnitude of Regularization Adjustments", fontsize=11, fontweight="semibold")
    ax2.grid(True, linestyle="--", alpha=0.4)
    ax2.legend(loc="upper right", fontsize=9, framealpha=0.9)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    os.makedirs(output_dir, exist_ok=True)
    saved_paths = []
    for fmt in formats:
        out_p = os.path.join(output_dir, f"covariance_alignment_summary_vs_mass.{fmt}")
        fig.savefig(out_p, dpi=200 if fmt == "png" else None)
        saved_paths.append(out_p)

    plt.close(fig)
    return saved_paths


def plot_covariance_for_directory(
    cov_dir: str,
    output_dir: Optional[str] = None,
    mass_bins_filter: Optional[List[str]] = None,
    tbin: str = "",
    max_l: Union[int, str] = "",
    formats: Optional[List[str]] = None,
) -> List[str]:
    """
    Finds and plots all pairs of raw and regularized covariance matrices in a directory:
      Raw:     covariance_matrix_<lbl>.root
      Aligned: BruCovMatrix_<lbl>.root
    """
    if not os.path.isdir(cov_dir):
        print(f"⚠️ Warning: Directory does not exist: {cov_dir}")
        return []

    if output_dir is None:
        output_dir = os.path.join(cov_dir, "plots_covariance")

    raw_files = sorted(glob.glob(os.path.join(cov_dir, "covariance_matrix_*.root")))
    if not raw_files:
        print(f"⚠️ No covariance_matrix_*.root files found in {cov_dir}")
        return []

    print(f"\n[Plotting Covariance Matrices] Found {len(raw_files)} covariance files in:\n  {cov_dir}")
    print(f"  Destination directory: {output_dir}")

    bin_results = []
    all_saved_plots = []

    for raw_p in raw_files:
        fname = os.path.basename(raw_p)
        # Match covariance_matrix_<lbl>.root
        lbl = fname.replace("covariance_matrix_", "").replace(".root", "")

        # Check mass bin filter if specified
        if mass_bins_filter:
            match = False
            for f in mass_bins_filter:
                if f in lbl or f in fname:
                    match = True
                    break
            if not match:
                continue

        aligned_p = os.path.join(cov_dir, f"BruCovMatrix_{lbl}.root")
        if not os.path.isfile(aligned_p):
            print(f"  ⚠️ Skipping {lbl}: matching BruCovMatrix_{lbl}.root not found.")
            continue

        res = plot_covariance_comparison_single_bin(
            raw_path=raw_p,
            aligned_path=aligned_p,
            output_dir=output_dir,
            mass_bin_label=lbl,
            tbin=tbin,
            max_l=max_l,
            formats=formats,
        )
        if res:
            bin_results.append(res)
            all_saved_plots.extend(res["saved_plots"])
            print(
                f"  ✅ Plotted {lbl}: Raw indef={res['raw_n_neg_eigenvals']} -> "
                f"Aligned SPD={res['aligned_is_spd']} (min eval {res['aligned_min_eigenval']:.2e})"
            )

    # Generate overall summary plot across all mass bins if multiple bins processed
    if len(bin_results) > 1:
        summary_plots = plot_covariance_summary_across_mass(
            bin_results=bin_results,
            output_dir=output_dir,
            tbin=tbin,
            max_l=max_l,
            formats=formats,
        )
        if summary_plots:
            all_saved_plots.extend(summary_plots)
            print(f"  📊 Generated multi-bin covariance alignment summary: {summary_plots[0]}")

    print(f"  Total plots generated: {len(all_saved_plots)}")
    return all_saved_plots


def plot_covariance_from_config(
    cfg: Any,
    tbins_filter: Optional[List[str]] = None,
    max_l_filter: Optional[List[int]] = None,
    mass_bins_filter: Optional[List[str]] = None,
    formats: Optional[List[str]] = None,
    output_dir: Optional[str] = None,
) -> Dict[str, List[str]]:
    """
    Plots covariance comparison across all target tbins and maxL configurations.
    """
    tbins = tbins_filter if tbins_filter is not None else cfg.analysis.tbins
    max_ls = max_l_filter if max_l_filter is not None else cfg.analysis.max_l
    beam_pol = cfg.analysis.beam_pol_label

    results: Dict[str, List[str]] = {}

    for tbin in tbins:
        for l in max_ls:
            cov_dir = cfg.linear_algebra.get_output_dir(tbin, l, beam_pol)
            key = f"{tbin}_maxL{l}"
            target_out = output_dir if output_dir else os.path.join(cov_dir, "plots_covariance")
            plots = plot_covariance_for_directory(
                cov_dir=cov_dir,
                output_dir=target_out,
                mass_bins_filter=mass_bins_filter,
                tbin=tbin,
                max_l=l,
                formats=formats,
            )
            results[key] = plots

    return results

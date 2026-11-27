#!/usr/bin/env python3
"""
MomentFitCollector.py

Automated diagnostic collector and next-fit generator for BruFit MomentFit pipeline.
Evaluates BruFit MCMC fit results across all t-bins and mass bins:
  1. Identifies unsuccessful bins (fit failure, stuck chains, or ESS < threshold).
  2. Diagnoses the specific failure mechanism:
       - Case 1: Zero acceptance rate in first chain (Phase 1) -> scale down initial value scale.
       - Case 2: Stuck MCMC in covariance chain (Phase 2)       -> scale down step size.
       - Case 3: Succeeded but mean ESS < threshold (low mixing) -> increase step size and iterations.
  3. Formulates parameter recommendations for the next fit.
  4. Generates a standalone, executable retry script (e.g. RunMomentFit_retry.sh) that:
       - Safely backs up prior failed bin results.
       - Re-fits only the specific unsuccessful mass bins via MomentFit.C in parallel.
       - Re-monitors the run directory on completion and triggers acceptance if all pass.
  5. Supports interactive confirmation and an --auto-retry flag for automated iteration loops.

Usage:
  ./MomentFitCollector.py --run-dir /path/to/run_dir [--tbin all|0|t010020] [--ess-threshold 15] [--ess-metric mean] [--auto-retry]
"""

from __future__ import annotations
import os
import sys
import glob
import re
import json
import shutil
import argparse
import subprocess
import concurrent.futures
import struct
import tempfile
from datetime import datetime

# Default maximum concurrent jobs: half of available CPUs
DEFAULT_MAX_PARALLEL = max(1, (os.cpu_count() or 2) // 2)
MAX_RETRY_ITERS = 200000


def _scale_retry_iters(current_iters: int, factor: float, floor_iters: int = 12000, max_iters: int = MAX_RETRY_ITERS) -> int:
    """Increase MCMC iteration count monotonically toward max_iters."""
    if current_iters >= max_iters:
        return current_iters
    scaled = int(round((current_iters * factor) / 1000.0) * 1000)
    return min(max_iters, max(floor_iters, current_iters + 1000, scaled))


def recommend_low_ess_params(
    bin_scale: float,
    bin_step: float,
    bin_iters: int,
    burn_in: int,
    eval_ess: float,
    ess_threshold: float,
    mean_tau: float | None,
    bin_acceptance: float | None,
    metric_name: str = "mean",
    max_iters: int = MAX_RETRY_ITERS,
    bin_gibbs: int = 0,
    low_acceptance_threshold: float = 0.03,
    gibbs_low_acceptance_threshold: float = 0.12,
    high_acceptance_threshold: float = 0.35,
) -> tuple[float, float, int, str]:
    """Compute retry (scale, step, iters, reason) for a converged bin with low ESS."""
    eff_post_samples = max(100, bin_iters - burn_in)
    tau_str = f"{mean_tau:.1f}" if mean_tau is not None else "N/A"
    eff_low_acc = gibbs_low_acceptance_threshold if bin_gibbs > 0 else low_acceptance_threshold

    if bin_acceptance is not None and bin_acceptance > high_acceptance_threshold:
        rec_scale = bin_scale
        rec_step = min(2.5, round(bin_step * 1.5, 3))
        rec_iters = _scale_retry_iters(bin_iters, factor=1.25, floor_iters=bin_iters, max_iters=max_iters)
        reason = (
            f"Low {metric_name} ESS ({eval_ess:.1f} < {ess_threshold:.1f}) with high acceptance ({bin_acceptance:.1%} > {high_acceptance_threshold:.1%}). "
            f"Proposals accepted too easily (tiny diffusive steps, tau={tau_str}). "
            f"Increasing step size ({bin_step:.3f} -> {rec_step:.3f}) and iterations ({bin_iters} -> {rec_iters}) to improve decorrelation."
        )
        return rec_scale, rec_step, rec_iters, reason

    if bin_acceptance is not None and bin_acceptance < eff_low_acc:
        rec_scale = bin_scale
        rec_step = max(0.05, round(bin_step * 0.6, 3))
        rec_iters = _scale_retry_iters(bin_iters, factor=1.25, floor_iters=bin_iters, max_iters=max_iters)
        mode_lbl = f"Gibbs-{bin_gibbs}" if bin_gibbs > 0 else "full-vector"
        reason = (
            f"Low {metric_name} ESS ({eval_ess:.1f} < {ess_threshold:.1f}) with low acceptance ({bin_acceptance:.1%} < {eff_low_acc:.1%}, {mode_lbl}). "
            f"Proposal jumps too large in covariance space (frequent rejections, tau={tau_str}). "
            f"Decreasing step size ({bin_step:.3f} -> {rec_step:.3f}) and increasing iterations ({bin_iters} -> {rec_iters}) to boost acceptance."
        )
        return rec_scale, rec_step, rec_iters, reason

    tau_ref = mean_tau if (mean_tau is not None and mean_tau > 1.0) else (eff_post_samples / max(1.0, eval_ess))
    needed_post = tau_ref * ess_threshold * 1.25
    needed_iters = int(round((burn_in + needed_post) / 1000.0) * 1000)
    if bin_iters >= max_iters:
        rec_iters = bin_iters
    else:
        rec_iters = min(max_iters, max(12000, needed_iters))
        if rec_iters <= bin_iters:
            rec_iters = _scale_retry_iters(bin_iters, factor=1.5, floor_iters=12000, max_iters=max_iters)
    rec_step = bin_step
    rec_scale = bin_scale
    h_acc_str = f"healthy acceptance ({bin_acceptance:.1%})" if bin_acceptance is not None else "measured autocorrelation"
    reason = (
        f"Low {metric_name} ESS ({eval_ess:.1f} < {ess_threshold:.1f}) with {h_acc_str}. "
        f"Autocorrelation time (tau={tau_ref:.1f}) requires scaling sampling iterations "
        f"({bin_iters} -> {rec_iters}) while maintaining calibrated step size ({rec_step:.3f})."
    )
    return rec_scale, rec_step, rec_iters, reason


def compute_retry_burnin(base_burn_in: int, rec_iters: int, burnin_fraction: float = 0.10) -> int:
    """Scale retry burn-in proportionally with iterations while keeping it strictly below half the chain."""
    if rec_iters <= 0:
        return max(0, base_burn_in)
    max_safe = max(0, rec_iters // 2)
    if burnin_fraction <= 0.0:
        return min(base_burn_in, max_safe) if max_safe > 0 else max(0, base_burn_in)
    prop = int(round((float(burnin_fraction) * rec_iters) / 500.0) * 500)
    scaled = max(base_burn_in, prop)
    return min(scaled, max_safe) if max_safe > 0 else scaled


def recommend_retry_gibbs_block_size(
    max_l: int,
    bin_gibbs: int,
    bin_acceptance: float | None,
    failure_type: str,
    gibbs_activation_acceptance: float = 0.05,
) -> int:
    """Enable or refine blocked Gibbs updates when global covariance proposals accept poorly or stall."""
    if bin_gibbs > 0:
        if failure_type == "STUCK_COVARIANCE" or (bin_acceptance is not None and bin_acceptance < 0.02):
            return max(4, min(bin_gibbs, max(4, bin_gibbs // 2)))
        return bin_gibbs

    needs_gibbs = (
        failure_type == "STUCK_COVARIANCE"
        or (max_l >= 8 and failure_type in ("LOW_ESS", "ZOMBIE_ROOT", "CORRUPT_ROOT", "FAILED_OTHER"))
        or (
            failure_type in ("LOW_ESS", "ZOMBIE_ROOT", "CORRUPT_ROOT", "FAILED_OTHER")
            and bin_acceptance is not None
            and ((bin_acceptance < gibbs_activation_acceptance and max_l >= 6) or bin_acceptance < 0.01)
        )
    )
    if not needs_gibbs:
        return 0
    if max_l >= 8:
        return 16
    if max_l >= 6:
        return 12
    return 8


def _compute_shrinkage_spd_covariance(
    samples,
    weights,
    param_names: list[str],
    prior_cov=None,
    shrinkage: float = 0.2,
    min_eigen_floor: float = 1e-4,
    yield_variance_floor: float = 10000.0,
):
    """Blend weighted posterior covariance with optional prior covariance and project onto SPD cone."""
    import numpy as np

    X = np.asarray(samples, dtype=np.float64)
    n_rows, n_cols = X.shape
    if n_rows < 2 or n_cols == 0:
        raise ValueError("Need at least 2 chain rows and 1 parameter to compute empirical covariance.")

    w = np.asarray(weights, dtype=np.float64) if weights is not None else np.ones(n_rows, dtype=np.float64)
    w = np.where(np.isfinite(w) & (w > 0.0), w, 0.0)
    w_sum = float(np.sum(w))
    if w_sum <= 0.0:
        w = np.ones(n_rows, dtype=np.float64)
        w_sum = float(n_rows)

    mu = np.sum(X * w[:, None], axis=0) / w_sum
    centered = X - mu
    denom = max(w_sum - 1.0, 1.0)
    emp_cov = (centered.T * w) @ centered / denom
    emp_cov = 0.5 * (emp_cov + emp_cov.T)

    if prior_cov is not None:
        prior_arr = np.asarray(prior_cov, dtype=np.float64)
        if prior_arr.shape != emp_cov.shape:
            raise ValueError(f"Prior shape {prior_arr.shape} != empirical shape {emp_cov.shape}")
        shrink = float(np.clip(shrinkage, 0.0, 1.0))
        blended = (1.0 - shrink) * emp_cov + shrink * (0.5 * (prior_arr + prior_arr.T))
    else:
        blended = emp_cov

    blended = 0.5 * (blended + blended.T)

    yield_idx = None
    for idx, pname in enumerate(param_names):
        if pname.startswith("Yld_") or "Yield" in pname:
            yield_idx = idx
            break

    if yield_idx is not None:
        y_var = abs(float(blended[yield_idx, yield_idx]))
        if y_var <= 1.0:
            y_var = float(yield_variance_floor)
        blended[yield_idx, :] = 0.0
        blended[:, yield_idx] = 0.0
        blended[yield_idx, yield_idx] = y_var

    diag_abs = np.abs(np.diag(blended))
    sigma = np.where(diag_abs > 1e-8, np.sqrt(diag_abs), 0.05)
    if yield_idx is not None:
        sigma[yield_idx] = float(np.sqrt(blended[yield_idx, yield_idx]))

    R = blended / np.outer(sigma, sigma)
    if yield_idx is not None:
        R[yield_idx, :] = 0.0
        R[:, yield_idx] = 0.0
    np.fill_diagonal(R, 1.0)
    R = 0.5 * (R + R.T)

    eigvals, eigvecs = np.linalg.eigh(R)
    eigvals = np.maximum(eigvals, float(min_eigen_floor))
    spd_R = (eigvecs * eigvals) @ eigvecs.T
    spd_R = 0.5 * (spd_R + spd_R.T)
    d = np.sqrt(np.maximum(np.diag(spd_R), 1e-12))
    spd_R = spd_R / np.outer(d, d)
    np.fill_diagonal(spd_R, 1.0)

    spd_cov = 0.5 * ((spd_R * np.outer(sigma, sigma)) + (spd_R * np.outer(sigma, sigma)).T)
    try:
        np.linalg.cholesky(spd_cov)
    except np.linalg.LinAlgError:
        spd_cov = spd_cov + np.diag((sigma ** 2) * 1e-4)
        np.linalg.cholesky(spd_cov)
    return spd_cov


def _load_bru_cov_matrix_root(cov_file: str, matrix_name: str = "MyCovMatrix"):
    """Load TMatrixDSym and ordered parameter names from a BruCovMatrix ROOT file."""
    import numpy as np

    root = _get_root()
    if root is None or not os.path.isfile(cov_file):
        return None, None
    tfile = root.TFile.Open(cov_file, "READ")
    if not tfile or tfile.IsZombie():
        if tfile:
            tfile.Close()
        return None, None
    try:
        mat = tfile.Get(matrix_name)
        name_list = tfile.Get(f"{matrix_name}_Names")
        if not mat or not name_list:
            return None, None
        n = int(mat.GetNrows())
        arr = np.empty((n, n), dtype=np.float64)
        for i in range(n):
            for j in range(n):
                arr[i, j] = float(mat[i][j])
        names = [str(obj.GetString()) for obj in name_list]
        return arr, names
    finally:
        tfile.Close()


def _save_bru_cov_matrix_root(out_path: str, cov_array, param_names: list[str], matrix_name: str = "MyCovMatrix") -> bool:
    """Save an SPD numpy covariance matrix and parameter name list in BruFit format."""
    root = _get_root()
    if root is None:
        return False
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    n = len(param_names)
    tmat = root.TMatrixDSym(n)
    for i in range(n):
        for j in range(n):
            tmat[i][j] = float(cov_array[i, j])
    name_list = root.TList()
    name_list.SetName(f"{matrix_name}_Names")
    for pname in param_names:
        name_list.Add(root.TObjString(pname))
    tfile = root.TFile.Open(out_path, "RECREATE")
    if not tfile or tfile.IsZombie():
        if tfile:
            tfile.Close()
        return False
    try:
        tmat.Write(matrix_name)
        name_list.Write(name_list.GetName(), root.TObject.kSingleKey)
    finally:
        tfile.Close()
    return True


def _resolve_source_cov_dir(base_cov_dir: str, tbin_name: str, max_l: int) -> str:
    """Resolve Stage 2 covariance directory matching MomentFit.C convention."""
    if not base_cov_dir:
        cand = os.path.join(repo_root, "Moments.out", "nominal", tbin_name, f"All.maxL_{max_l}")
        return cand if os.path.isdir(cand) else ""
    if tbin_name in base_cov_dir:
        return base_cov_dir.rstrip("/")
    return os.path.join(base_cov_dir, tbin_name, f"All.maxL_{max_l}")


def recondition_covariance_for_bin(
    bin_record: dict,
    run_dir: str,
    src_cov_dir: str,
    max_l: int,
    burn_in: int = 2000,
    shrinkage: float = 0.2,
    min_eigen_floor: float = 1e-4,
    min_post_samples: int = 200,
) -> str | None:
    """
    Build an opt-in shrinkage-regularized BruCovMatrix from post-burn-in MCMCTree samples
    under <run_dir>/<tbin>/reconditioned_cov without modifying Stage 2 files.
    """
    import numpy as np

    root_file = bin_record.get("root_file")
    tbin_name = bin_record.get("tbin_name", "")
    m_center = bin_record.get("mass_center")
    if not root_file or not os.path.isfile(root_file) or m_center is None or uproot is None:
        return None

    resolved_src = _resolve_source_cov_dir(src_cov_dir, tbin_name, max_l)
    src_bru = os.path.join(resolved_src, f"BruCovMatrix_Mpi0eta_{m_center:.3f}.root") if resolved_src else ""
    src_init = os.path.join(resolved_src, f"covariance_matrix_Mpi0eta_{m_center:.3f}.root") if resolved_src else ""
    prior_cov, prior_names = _load_bru_cov_matrix_root(src_bru) if src_bru else (None, None)

    try:
        with uproot.open(root_file) as f:
            if "MCMCTree" not in f:
                return None
            tree = f["MCMCTree"]
            all_keys = list(tree.keys())
            if prior_names:
                param_names = list(prior_names)
                read_cols = [n for n in param_names if n in all_keys]
            else:
                param_names = _select_h_columns(all_keys)
                if "Yld_Moments" in all_keys and "Yld_Moments" not in param_names:
                    param_names.append("Yld_Moments")
                read_cols = list(param_names)
            if not read_cols:
                return None

            weight_col = "weight" if "weight" in all_keys else None
            branches = read_cols + ([weight_col] if weight_col else [])
            arrays = tree.arrays(branches, library="np")
    except Exception as exc:
        print(f"  Warning: Could not read MCMCTree for covariance reconditioning ({tbin_name} bin {bin_record.get('bin_idx')}): {exc}")
        return None

    n_entries = len(arrays[read_cols[0]])
    if n_entries == 0:
        return None

    if weight_col:
        raw_w = np.asarray(arrays[weight_col], dtype=np.float64)
        int_w = np.maximum(1, np.round(np.where(np.isfinite(raw_w), raw_w, 1.0))).astype(np.int64)
    else:
        int_w = np.ones(n_entries, dtype=np.int64)

    cum_w = np.cumsum(int_w)
    if cum_w[-1] <= burn_in + min_post_samples:
        return None

    start_idx = int(np.searchsorted(cum_w, burn_in, side="right"))
    post_w = int_w[start_idx:].astype(np.float64)
    if cum_w[start_idx] > burn_in and start_idx < len(post_w) + start_idx:
        post_w[0] = float(cum_w[start_idx] - burn_in)
    if len(post_w) < 2 or float(np.sum(post_w)) < min_post_samples:
        return None

    if prior_names and prior_cov is not None and len(read_cols) < len(param_names):
        n_full = len(param_names)
        X_sub = np.column_stack([np.asarray(arrays[c][start_idx:], dtype=np.float64) for c in read_cols])
        idx_map = [param_names.index(c) for c in read_cols]
        sub_prior = prior_cov[np.ix_(idx_map, idx_map)]
        sub_spd = _compute_shrinkage_spd_covariance(
            X_sub, post_w, read_cols,
            prior_cov=sub_prior, shrinkage=shrinkage, min_eigen_floor=min_eigen_floor,
        )
        full_mat = np.array(prior_cov, dtype=np.float64, copy=True)
        for ii, orig_i in enumerate(idx_map):
            for jj, orig_j in enumerate(idx_map):
                full_mat[orig_i, orig_j] = sub_spd[ii, jj]
        spd_cov = full_mat
    else:
        X = np.column_stack([np.asarray(arrays[c][start_idx:], dtype=np.float64) for c in param_names])
        spd_cov = _compute_shrinkage_spd_covariance(
            X, post_w, param_names,
            prior_cov=prior_cov, shrinkage=shrinkage, min_eigen_floor=min_eigen_floor,
        )

    recond_dir = os.path.join(run_dir, tbin_name, "reconditioned_cov")
    out_bru = os.path.join(recond_dir, f"BruCovMatrix_Mpi0eta_{m_center:.3f}.root")
    if not _save_bru_cov_matrix_root(out_bru, spd_cov, param_names):
        return None

    out_init = os.path.join(recond_dir, f"covariance_matrix_Mpi0eta_{m_center:.3f}.root")
    if src_init and os.path.isfile(src_init) and not os.path.exists(out_init):
        try:
            shutil.copy2(src_init, out_init)
        except OSError:
            pass

    return recond_dir


def _resolve_pattern_for_tbin(pattern: str, tbin_name: str, data_dir: str) -> str:
    """Resolve {tbin}, {data_dir}, and {base_dir} placeholders for retry macro invocations."""
    if not pattern:
        return ""
    resolved = pattern.replace("{tbin}", tbin_name)
    if data_dir:
        resolved = resolved.replace("{data_dir}", data_dir).replace("{base_dir}", data_dir)
    return resolved

# Ensure local repository python libraries are found with highest priority
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
local_python = os.path.join(repo_root, "python")

search_paths = [
    local_python,
    "/d/home/septian/EtaPi0Analysis/python",
    "/d/grid17/septian/EtaPi0Analysis/python"
]
for p in reversed(search_paths):
    if os.path.isdir(p):
        while p in sys.path:
            sys.path.remove(p)
        sys.path.insert(0, p)

try:
    import uproot
except ImportError:
    uproot = None

# Lazy-load PyROOT only when needed (e.g., TCanvas projection scaling checks or uproot fallback).
# Deferring ROOT initialization avoids startup overhead on cache hits and keeps ProcessPoolExecutor forks clean.
ROOT = None
_ROOT_TRIED = False


def _get_root():
    global ROOT, _ROOT_TRIED
    if _ROOT_TRIED:
        return ROOT
    _ROOT_TRIED = True
    try:
        import ROOT as _root_mod
        _root_mod.gROOT.SetBatch(True)
        _root_mod.gErrorIgnoreLevel = _root_mod.kWarning
        ROOT = _root_mod
        return ROOT
    except Exception:
        root_lib_candidates = []
        if "ROOTSYS" in os.environ:
            root_lib_candidates.append(os.path.join(os.environ["ROOTSYS"], "lib"))
        root_lib_candidates.extend([
            "/d/grid17/septian/root_v6-36-08_cuda_install/lib",
            "/d/home/septian/root/lib",
            "/usr/local/root/lib"
        ])
        for rlib in root_lib_candidates:
            if os.path.isdir(rlib) and rlib not in sys.path:
                sys.path.append(rlib)
                try:
                    import ROOT as _root_mod
                    _root_mod.gROOT.SetBatch(True)
                    _root_mod.gErrorIgnoreLevel = _root_mod.kWarning
                    ROOT = _root_mod
                    return ROOT
                except Exception:
                    pass
    ROOT = None
    return None


FitResults = None

import numpy as np


TBIN_INDEX_MAP = {
    0: "t010020",
    1: "t020032",
    2: "t032050",
    3: "t050075",
    4: "t075100"
}
TBIN_NAME_TO_INDEX = {v: k for k, v in TBIN_INDEX_MAP.items()}


def _fast_binary_zombie_check(root_path: str) -> bool:
    """Fast pure-binary ROOT header check for missing/truncated/corrupt files (no ROOT/uproot open)."""
    if not root_path or not os.path.isfile(root_path):
        return False
    try:
        fsize = os.path.getsize(root_path)
    except OSError:
        return False
    if fsize < 10240:
        return True

    try:
        with open(root_path, "rb") as bf:
            magic = bf.read(4)
            if magic != b"root":
                return True
            version = struct.unpack(">i", bf.read(4))[0]
            bf.read(4)  # fBEGIN
            if version > 1000000:
                fEND = struct.unpack(">q", bf.read(8))[0]
                fSeekFree = struct.unpack(">q", bf.read(8))[0]
            else:
                fEND = struct.unpack(">i", bf.read(4))[0]
                fSeekFree = struct.unpack(">i", bf.read(4))[0]
            if fSeekFree == 0 or fEND < fsize - 4096:
                return True
    except Exception:
        return True
    return False


def is_root_file_zombie(root_path: str) -> bool:
    """
    Determine whether a ROOT file is a zombie, unclosed, truncated, or corrupted.
    Checks binary header first, then uproot keys / ROOT TFile.
    """
    if not root_path or not os.path.isfile(root_path):
        return False
    if _fast_binary_zombie_check(root_path):
        return True

    if uproot is not None:
        try:
            with uproot.open(root_path) as uf:
                return len(uf.keys()) == 0
        except Exception:
            return True

    root_mod = _get_root()
    if root_mod is not None:
        try:
            old_level = getattr(root_mod, "gErrorIgnoreLevel", 0)
            root_mod.gErrorIgnoreLevel = getattr(root_mod, "kFatal", 5000)
            f = root_mod.TFile.Open(root_path, "READ")
            root_mod.gErrorIgnoreLevel = old_level
            if not f:
                return True
            is_zomb = (
                f.IsZombie() or
                f.GetNkeys() == 0 or
                f.TestBit(getattr(root_mod.TFile, "kRecovered", 0x00000400))
            )
            try:
                f.Close()
            except Exception:
                pass
            if is_zomb:
                return True
        except Exception:
            return True

    return False


def check_mcmc_tree_entries(root_path: str, tree_name: str = "MCMCTree") -> tuple[bool, int]:
    """
    Check if root_path exists, is NOT a zombie/corrupt file, and has valid entries in tree_name (default MCMCTree).
    Returns (has_entries: bool, n_entries: int).
    """
    if not os.path.isfile(root_path):
        return False, 0

    if _fast_binary_zombie_check(root_path):
        return False, 0

    if uproot is not None:
        try:
            with uproot.open(root_path) as uf:
                if len(uf.keys()) == 0 or tree_name not in uf:
                    return False, 0
                n_entries = int(uf[tree_name].num_entries)
                return (n_entries > 0), n_entries
        except Exception:
            return False, 0

    root_mod = _get_root()
    if root_mod is None:
        return True, 1000

    try:
        f = root_mod.TFile.Open(root_path, "READ")
        if not f or f.IsZombie() or f.GetNkeys() == 0 or f.TestBit(getattr(root_mod.TFile, "kRecovered", 0x00000400)):
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


def check_projection_scaling(root_path: str, tolerance: float = 0.30, skip_zombie_check: bool = False) -> tuple[bool, float | None, str]:
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
    if not os.path.isfile(root_path) or (not skip_zombie_check and is_root_file_zombie(root_path)):
        return True, None, "File missing or zombie"

    root_mod = _get_root()
    if root_mod is None:
        return True, None, "ROOT not available to evaluate curve scaling"

    try:
        f = root_mod.TFile.Open(root_path, "READ")
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


def parse_run_dir_params(run_dir: str) -> dict:
    """Extract default parameters (maxL, scale, step, mass range) from run directory name."""
    params = {
        "max_l": 4,
        "init_scale": 1.0,
        "step": 1.0,
        "m_min": 0.8,
        "m_max": 2.0,
        "n_bins": 30,
        "mcmc_iterations": 10000,
        "mcmc_burnin": 2000
    }
    bname = os.path.basename(os.path.normpath(run_dir))

    # Match Lmax
    m_l = re.search(r"Lmax(\d+)", bname)
    if m_l:
        params["max_l"] = int(m_l.group(1))

    # Match InitValueScale
    m_scale = re.search(r"InitValueScale([0-9.]+)", bname)
    if m_scale:
        params["init_scale"] = float(m_scale.group(1))

    # Match step
    m_step = re.search(r"step([0-9.]+)", bname)
    if m_step:
        params["step"] = float(m_step.group(1))

    # Match 6-digit mass tag e.g. m080200 -> 0.80 to 2.00, m104172 -> 1.04 to 1.72
    m_mass6 = re.search(r"m(\d{3})(\d{3})", bname)
    if m_mass6:
        params["m_min"] = float(m_mass6.group(1)) / 100.0
        params["m_max"] = float(m_mass6.group(2)) / 100.0
        params["n_bins"] = int(round((params["m_max"] - params["m_min"]) / 0.04))
    else:
        # Match 4-digit mass tag e.g. m0820 -> 0.8 to 2.0
        m_mass4 = re.search(r"m(\d{2})(\d{2})", bname)
        if m_mass4:
            params["m_min"] = float(m_mass4.group(1)) / 10.0
            params["m_max"] = float(m_mass4.group(2)) / 10.0
            params["n_bins"] = int(round((params["m_max"] - params["m_min"]) / 0.04))

    return params


def discover_tbin_dirs(base_dir: str, tbin_filter: str = "all") -> list[tuple[int, str, str]]:
    """Discover matching t-bin folders: returns list of (tbin_idx, tbin_name, full_path)."""
    tbin_dirs = []
    if not os.path.isdir(base_dir):
        return []

    target_idx = None
    target_name = None
    if tbin_filter != "all":
        if tbin_filter.isdigit():
            target_idx = int(tbin_filter)
            target_name = TBIN_INDEX_MAP.get(target_idx)
        else:
            target_name = tbin_filter
            target_idx = TBIN_NAME_TO_INDEX.get(tbin_filter)

    # Check if base_dir is already a specific tbin folder
    bname = os.path.basename(os.path.normpath(base_dir))
    if re.match(r"^t\d{6}$", bname) and ".bak" not in bname:
        idx = TBIN_NAME_TO_INDEX.get(bname, 0)
        return [(idx, bname, os.path.abspath(base_dir))]

    for entry in sorted(os.listdir(base_dir)):
        p = os.path.join(base_dir, entry)
        if os.path.isdir(p) and re.match(r"^t\d{6}$", entry) and ".bak" not in entry:
            idx = TBIN_NAME_TO_INDEX.get(entry)
            if target_name is not None and entry != target_name and idx != target_idx:
                continue
            if idx is not None:
                tbin_dirs.append((idx, entry, os.path.abspath(p)))

    return tbin_dirs


def _safe_float(val, default=None):
    """Safely parse a string to float, returning default if None, N/A, -, nan, or invalid."""
    if val is None:
        return default
    if isinstance(val, (int, float)):
        return float(val)
    s = str(val).strip()
    if s in ("N/A", "-", "None", "nan", "NaN", "null", ""):
        return default
    try:
        return float(s)
    except (ValueError, TypeError):
        return default


def parse_diagnostics_file(diag_file: str) -> dict[int, dict]:
    """Parse existing MCMC_Diagnostics.txt or MomentFit_Diagnostics_<tbin>.txt."""
    if not os.path.isfile(diag_file):
        return {}

    parsed = {}
    try:
        with open(diag_file, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        # Parse Fit Status Summary Table
        in_status_table = False
        in_mcmc_summary_table = False
        in_caution_table = False

        for line in lines:
            line_str = line.strip()
            if "Fit Status Summary" in line_str:
                in_status_table = True
                in_mcmc_summary_table = False
                in_caution_table = False
                continue
            if in_status_table and line_str.startswith("Bin") and "Status" in line_str:
                continue
            if in_status_table and line_str.startswith("==="):
                continue
            if in_status_table and line_str.startswith("Total Bins:"):
                in_status_table = False
                continue

            if in_status_table and line_str and not line_str.startswith("---"):
                parts = line_str.split()
                if len(parts) >= 6 and parts[0].isdigit():
                    b_idx = int(parts[0])
                    mass_val = _safe_float(parts[1])
                    if mass_val is None:
                        mass_val = round(0.8 + 0.04 * b_idx + 0.02, 4)
                    status = parts[2]
                    # Parse Tau (mean/max)
                    tau_mean = None
                    tau_max = None
                    ess_min = None
                    ess_mean = None
                    # e.g. "1635.1 / 3430.2"
                    m_tau = re.search(r"([0-9.]+)\s*/\s*([0-9.]+)", line_str)
                    if m_tau:
                        tau_mean = _safe_float(m_tau.group(1))
                        tau_max = _safe_float(m_tau.group(2))
                    # Parse ESS (min/mean)
                    m_ess = re.findall(r"([0-9.]+)\s*/\s*([0-9.]+)", line_str)
                    if len(m_ess) >= 2:
                        ess_min = _safe_float(m_ess[1][0])
                        ess_mean = _safe_float(m_ess[1][1])

                    # Parse Acceptance if present (e.g. 0.234 or column 10)
                    acc_val = None
                    m_acc = re.search(r"\s+([0-9]\.[0-9]{3})\s+", line_str)
                    if m_acc:
                        val = _safe_float(m_acc.group(1))
                        if val is not None and 0.0 <= val <= 1.0:
                            acc_val = val
                    elif len(parts) >= 11:
                        val = _safe_float(parts[10])
                        if val is not None and 0.0 <= val <= 1.0:
                            acc_val = val

                    parsed[b_idx] = {
                        "bin_idx": b_idx,
                        "mass": mass_val,
                        "status": status,
                        "tau_mean": tau_mean,
                        "tau_max": tau_max,
                        "ess_min": ess_min,
                        "ess_mean": ess_mean,
                        "acceptance": acc_val,
                        "worst_moment": None,
                        "n_low_ess": 0,
                        "total_moments": 0
                    }

            # Parse MCMC Diagnostics Summary Across All Successful Mass Bins Table
            if "MCMC Diagnostics Summary Across All Successful Mass Bins" in line_str:
                in_mcmc_summary_table = True
                in_status_table = False
                in_caution_table = False
                continue
            if in_mcmc_summary_table and line_str.startswith("Bin") and "Tau" in line_str:
                continue
            if in_mcmc_summary_table and (line_str.startswith("===") or line_str.startswith("---")):
                continue
            if in_mcmc_summary_table and (line_str.startswith("[") or "Cautionary Bins" in line_str):
                in_mcmc_summary_table = False

            if in_mcmc_summary_table and line_str and not line_str.startswith("-") and not line_str.startswith("="):
                parts = line_str.split()
                if len(parts) >= 7 and parts[0].isdigit():
                    b_idx = int(parts[0])
                    m_val = _safe_float(parts[1])
                    if m_val is None:
                        m_val = round(0.8 + 0.04 * b_idx + 0.02, 4)
                    t_mean = _safe_float(parts[3])
                    t_max = _safe_float(parts[4])
                    e_min = _safe_float(parts[5])
                    e_mean = _safe_float(parts[6])
                    if b_idx not in parsed:
                        parsed[b_idx] = {
                            "bin_idx": b_idx,
                            "mass": m_val,
                            "status": "SUCCESS",
                            "tau_mean": t_mean,
                            "tau_max": t_max,
                            "ess_min": e_min,
                            "ess_mean": e_mean,
                            "worst_moment": None,
                            "n_low_ess": 0,
                            "total_moments": 0
                        }
                    else:
                        parsed[b_idx]["tau_mean"] = t_mean
                        parsed[b_idx]["tau_max"] = t_max
                        parsed[b_idx]["ess_min"] = e_min
                        parsed[b_idx]["ess_mean"] = e_mean

            # Parse Cautionary Bins table
            if "Cautionary Bins with ESS <" in line_str:
                in_caution_table = True
                continue
            if in_caution_table and line_str.startswith("Bin") and "Low ESS" in line_str:
                continue
            if in_caution_table and (line_str.startswith("Note:") or line_str.startswith("===") or line_str.startswith("---") and not line_str.startswith("----------------------------------------------------------------------------------------------------")):
                if line_str.startswith("Note:"):
                    in_caution_table = False
                continue

            if in_caution_table and line_str and not line_str.startswith("-") and not line_str.startswith("="):
                parts = line_str.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    b_idx = int(parts[0])
                    worst_m = parts[3]
                    worst_ess = _safe_float(parts[4])
                    # Low ESS Moments e.g. 39/39
                    m_frac = re.search(r"(\d+)/(\d+)", parts[2])
                    n_bad = int(m_frac.group(1)) if m_frac else 0
                    n_tot = int(m_frac.group(2)) if m_frac else 0

                    if b_idx in parsed:
                        parsed[b_idx]["worst_moment"] = worst_m
                        parsed[b_idx]["ess_min"] = worst_ess
                        parsed[b_idx]["n_low_ess"] = n_bad
                        parsed[b_idx]["total_moments"] = n_tot
    except Exception as e:
        print(f"Warning: Failed to parse diagnostic file {diag_file}: {e}")

    return parsed


def extract_hyperparameters_from_log(log_path: str) -> tuple[float | None, float | None, int | None]:
    """
    Extract (scale, step, iters) from a fit log file.
    Checks macro invocation line and MCMC configuration prints.
    """
    if not os.path.isfile(log_path):
        return None, None, None
    try:
        with open(log_path, "r", errors="ignore") as f:
            head = f.read(50000)
    except Exception:
        return None, None, None

    scale, step, iters = None, None, None

    # Macro call pattern: MomentFit.C(tbin, maxL, scale, step, bin, iters, ...)
    m_macro = re.search(r'MomentFit(?:_General)?\.C\s*\(\s*[-0-9]+\s*,\s*[-0-9]+\s*,\s*([0-9.]+)\s*,\s*([0-9.]+)\s*,\s*[-0-9]+\s*,\s*(\d+)', head)
    if m_macro:
        scale = float(m_macro.group(1))
        step = float(m_macro.group(2))
        iters = int(m_macro.group(3))
        return scale, step, iters

    # Legacy MomentFit_m080200.C(tbin, maxL, scale, step, ...)
    m_legacy = re.search(r'MomentFit_m080200\.C\s*\(\s*[-0-9]+\s*,\s*[-0-9]+\s*,\s*([0-9.]+)\s*,\s*([0-9.]+)', head)
    if m_legacy:
        scale = float(m_legacy.group(1))
        step = float(m_legacy.group(2))

    # Configuration prints
    m_samples = re.search(r'Number of samples:\s*(\d+)', head)
    if m_samples:
        iters = int(m_samples.group(1))

    m_step = re.search(r'Scale step size:\s*([0-9.]+)', head)
    if m_step and step is None:
        step = float(m_step.group(1))

    m_scale = re.search(r'(?:Initial Values Pattern:.*\(Scale:\s*([0-9.]+)\)|-->\s*Init Scale\s*:\s*([0-9.]+)|Initial value scale:\s*([0-9.]+))', head)
    if m_scale and scale is None:
        scale = float(m_scale.group(1) or m_scale.group(2) or m_scale.group(3))

    return scale, step, iters


def extract_gibbs_from_log(log_path: str) -> int | None:
    """Extract Gibbs block size from a MomentFit log file if recorded."""
    if not os.path.isfile(log_path):
        return None
    try:
        with open(log_path, "r", errors="ignore") as f:
            head = f.read(50000)
    except Exception:
        return None

    m_gibbs_print = re.search(r'Gibbs block size:\s*(\d+)', head)
    if m_gibbs_print:
        return int(m_gibbs_print.group(1))

    m_macro = re.search(
        r'MomentFit(?:_General)?\.C\s*\([^;)]*?,\s*(?:true|false)\s*,\s*(?:true|false)\s*,\s*(?:true|false)\s*,\s*(\d+)',
        head,
    )
    if m_macro:
        return int(m_macro.group(1))
    return None


def analyze_log_content(log_path: str) -> dict:
    """Analyze a single log file to detect failure stages, acceptance, and bailout."""
    result = {
        "exists": False,
        "status": "UNKNOWN",
        "phase1_zero_acc": False,
        "phase2_stuck": False,
        "last_acceptance": None,
        "bailout": False,
        "details": ""
    }
    if not os.path.isfile(log_path):
        result["status"] = "MISSING"
        result["details"] = "Log file not found"
        return result

    result["exists"] = True
    fsize = os.path.getsize(log_path)
    if fsize == 0:
        result["status"] = "EMPTY"
        result["details"] = "Log file is 0 bytes"
        return result

    try:
        with open(log_path, "r", errors="ignore") as f:
            if fsize > 400000:
                head = f.read(60000)
                f.seek(max(0, fsize - 300000))
                tail = f.read()
                content = head + "\n" + tail
            else:
                content = f.read()
                head = content
                tail = content
    except Exception as e:
        result["status"] = "READ_ERROR"
        result["details"] = str(e)
        return result

    # First check: If a completed ROOT file with MCMCTree exists, fit is strictly SUCCESS
    log_dir = os.path.dirname(os.path.abspath(log_path))
    m_idx = re.search(r'(?:logRooFit|bin)(\d+)', os.path.basename(log_path))
    b_idx = int(m_idx.group(1)) if m_idx else None
    m_val = round(0.8 + 0.04 * b_idx + 0.02, 4) if b_idx is not None else None
    m_head = re.search(r'Mpi0eta([0-9.]+)', head)
    if m_head:
        m_val = float(m_head.group(1))

    search_dirs = [log_dir]
    m_tbin = re.search(r'(t\d{6})', os.path.basename(log_path))
    if os.path.basename(log_dir) == "logs":
        parent = os.path.dirname(log_dir)
        if m_tbin and os.path.isdir(os.path.join(parent, m_tbin.group(1))):
            search_dirs = [os.path.join(parent, m_tbin.group(1))]
        elif os.path.isdir(parent):
            search_dirs = [parent]
    elif re.match(r"^t\d{6}$", os.path.basename(log_dir)):
        search_dirs = [log_dir]

    cand_dirs = []
    for sdir in search_dirs:
        if m_val is not None:
            for cname in [f"Mpi0eta{m_val:.6f}_", f"Mpi0eta{m_val:.6f}", f"Mpi0eta{m_val:.4f}_", f"Mpi0eta{m_val:.4f}"]:
                cp = os.path.join(sdir, cname)
                if os.path.isdir(cp) and ".bak" not in cp and cp not in cand_dirs:
                    cand_dirs.append(cp)
            if os.path.isdir(sdir):
                for entry in os.listdir(sdir):
                    cp = os.path.join(sdir, entry)
                    if os.path.isdir(cp) and ".bak" not in entry:
                        m_ent = re.search(r"Mpi0eta([0-9.]+)", entry)
                        if m_ent and abs(float(m_ent.group(1)) - m_val) < 0.005:
                            if cp not in cand_dirs:
                                cand_dirs.append(cp)

    stuck_in_cand = False
    zombie_in_cand = None
    for cdir in cand_dirs:
        for rname in ["ResultsBruMcmcFixedCovariance.root", "ResultsBruMcmcCovariance.root", "ResultsBruMcmc.root"]:
            rp = os.path.join(cdir, rname)
            if os.path.isfile(rp):
                if is_root_file_zombie(rp):
                    if not zombie_in_cand:
                        zombie_in_cand = rp
                    continue
                has_entries, n_ent = check_mcmc_tree_entries(rp)
                if has_entries:
                    result["status"] = "SUCCESS"
                    result["details"] = f"{rname} present (MCMCTree: {n_ent} entries)"
                    return result
                elif (0 < n_ent <= 50) or os.path.getsize(rp) < 120000:
                    stuck_in_cand = True

    # Check for fast bailout first (overrides any partial write indicators)
    has_bailout = ("Fast Bailout Triggered" in content or
                   "Bailout" in tail or "bailout" in tail.lower() or stuck_in_cand)

    zero_err_count = len(re.findall(r"H_\d_\d_\d\s+[-0-9.eE]+\s+\+-\s+0(?:\.0*)?(?!\d)", tail))
    has_zero_matrix = "-0          -0          -0          -0          -0" in tail
    if zero_err_count >= 10 or has_zero_matrix:
        has_bailout = True

    # Check for success indicators
    is_success = (not has_bailout) and (
        ("PosteriorCovariance" in tail and "Freeing MC Tree Memory" in tail) or
        ("--> Matrix successfully written to file as 'PosteriorCovariance'" in tail)
    )
    if is_success:
        result["status"] = "SUCCESS"
        result["details"] = "Posterior covariance written successfully"
        return result

    if zombie_in_cand and not is_success:
        result["status"] = "ZOMBIE_ROOT"
        result["details"] = f"Zombie/corrupted {os.path.basename(zombie_in_cand)} detected (fit interrupted before completion)"
        result["zombie_file"] = zombie_in_cand
        return result

    # Check for failure keywords
    is_failed = ("WARNING FitManager::Run: Fit failed" in tail or
                 "Skipping PlotDataModel() due to failed fit" in tail or
                 "1DStep Failed! Retrying (10/10)" in tail or
                 "Error in <" in tail or
                 "FATAL" in tail)

    # Check Phase 1 vs Phase 2
    in_phase1 = "*** Starting Phase 1" in head or "*** Starting Phase 1" in tail
    in_phase2 = "*** Launching Official Covariance Chain" in content or "*** Phase 2" in content

    # Acceptance parsing
    acc_matches = re.findall(r"acceptance\s+([0-9.]+)", tail)
    if acc_matches:
        try:
            result["last_acceptance"] = float(acc_matches[-1])
        except ValueError:
            pass

    # Extract attempted parameters from macro invocation or configuration
    scale_match = re.search(r"\(Scale:\s*([0-9.]+)\)", content) or re.search(r"Init Scale\s*:\s*([0-9.]+)", content)
    if scale_match:
        try:
            result["scale"] = float(scale_match.group(1))
        except ValueError:
            pass

    step_match = re.search(r"Scale step size:\s*([0-9.]+)", content)
    if step_match:
        try:
            result["step"] = float(step_match.group(1))
        except ValueError:
            pass

    iters_match = re.search(r"Number of samples:\s*([0-9]+)", content)
    if iters_match:
        try:
            result["iters"] = int(iters_match.group(1))
        except ValueError:
            pass

    # Zero acceptance in Phase 1
    ph1_zero = False
    if "1DStep Failed! Retrying (10/10)" in content:
        ph1_zero = True
    elif in_phase1 and not in_phase2:
        if "accepted 0 out of" in tail or (result["last_acceptance"] is not None and result["last_acceptance"] < 0.005):
            ph1_zero = True

    # Stuck in Phase 2
    ph2_stuck = False
    if in_phase2:
        if has_bailout:
            ph2_stuck = True
            result["bailout"] = True
        elif result["last_acceptance"] is not None and result["last_acceptance"] < 0.005:
            ph2_stuck = True
        elif is_failed or "BruSequentialProposal" in tail:
            ph2_stuck = True

    result["phase1_zero_acc"] = ph1_zero
    result["phase2_stuck"] = ph2_stuck

    if ph1_zero:
        result["status"] = "ZERO_ACC_PH1"
        result["details"] = "Initial burn-in chain stuck (acceptance near 0 / 1DStep failed)"
    elif ph2_stuck:
        result["status"] = "STUCK_COVARIANCE"
        result["details"] = "Covariance chain stuck or bailed out (acceptance near 0)"
    elif is_failed:
        result["status"] = "FAILED"
        result["details"] = "FitManager::Run returned failed"
    else:
        result["status"] = "INCOMPLETE"
        result["details"] = "MCMC interrupted or incomplete"

    return result


if FitResults is not None and hasattr(FitResults, "compute_brufit_chain_diagnostics"):
    compute_brufit_chain_diagnostics = FitResults.compute_brufit_chain_diagnostics
else:
    def _vehtari_norm_quantile(p):
        p_arr = np.atleast_1d(np.asarray(p, dtype=float))
        a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
        b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01]
        c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
        d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00]
        p_low, p_high = 0.02425, 1.0 - 0.02425
        q = np.zeros_like(p_arr)
        mask_mid = (p_arr >= p_low) & (p_arr <= p_high)
        if np.any(mask_mid):
            q_val = p_arr[mask_mid] - 0.5
            r = q_val * q_val
            q[mask_mid] = ((((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q_val) / ((((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0))
        mask_low = (p_arr < p_low) & (p_arr > 0.0)
        if np.any(mask_low):
            r = np.sqrt(-2.0 * np.log(p_arr[mask_low]))
            q[mask_low] = ((((((c[0]*r + c[1])*r + c[2])*r + c[3])*r + c[4])*r + c[5])) / (((((d[0]*r + d[1])*r + d[2])*r + d[3])*r + 1.0))
        mask_high = (p_arr > p_high) & (p_arr < 1.0)
        if np.any(mask_high):
            r = np.sqrt(-2.0 * np.log(1.0 - p_arr[mask_high]))
            q[mask_high] = -((((((c[0]*r + c[1])*r + c[2])*r + c[3])*r + c[4])*r + c[5])) / (((((d[0]*r + d[1])*r + d[2])*r + d[3])*r + 1.0))
        return float(q[0]) if np.isscalar(p) else q

    def _vehtari_rank_normalize(x):
        x_arr = np.asarray(x, dtype=float)
        n = len(x_arr)
        if n == 0:
            return np.array([], dtype=float)
        order = np.argsort(x_arr)
        ranks = np.empty(n, dtype=float)
        ranks[order] = np.arange(1, n + 1, dtype=float)
        p = (ranks - 0.375) / (n + 0.25)
        return _vehtari_norm_quantile(p)

    def _vehtari_diagnostics(z):
        n = len(z)
        half = n // 2
        if half < 2:
            return 0.0, 1.0, float(n)
        z1, z2 = z[:half], z[half:2 * half]
        m1, m2 = np.mean(z1), np.mean(z2)
        v1 = np.sum((z1 - m1)**2) / (half - 1.0)
        v2 = np.sum((z2 - m2)**2) / (half - 1.0)
        W = 0.5 * (v1 + v2)
        if W <= 0.0:
            return 1.0, 1.0, float(2 * half)
        grand_mean = (m1 + m2) / 2.0
        B = half * ((m1 - grand_mean)**2 + (m2 - grand_mean)**2)
        var_plus = ((half - 1.0) / half) * W + B / half
        rhat = np.sqrt(var_plus / W)

        if half >= 64:
            n_fft = 2 ** int(np.ceil(np.log2(2 * half - 1)))
            f1 = np.fft.fft(z1, n=n_fft)
            c1 = np.fft.ifft(f1 * np.conj(f1)).real[:half]
            f2 = np.fft.fft(z2, n=n_fft)
            c2 = np.fft.ifft(f2 * np.conj(f2)).real[:half]

            cs1_sq = np.cumsum(np.concatenate(([0.0], z1**2)))
            cs2_sq = np.cumsum(np.concatenate(([0.0], z2**2)))

            t = np.arange(1, half)
            s1 = (cs1_sq[half] - cs1_sq[t]) + (cs1_sq[half - t]) - 2.0 * c1[1:]
            s2 = (cs2_sq[half] - cs2_sq[t]) + (cs2_sq[half - t]) - 2.0 * c2[1:]
            vt = (s1 + s2) / (2.0 * (half - t))

            rho = 1.0 - (vt / (2.0 * var_plus))
            neg_idx = np.where(rho < 0.0)[0]
            if len(neg_idx) > 0:
                sum_rho = float(np.sum(rho[:neg_idx[0]]))
            else:
                sum_rho = float(np.sum(rho))
        else:
            sum_rho = 0.0
            for t in range(1, half):
                diff1, diff2 = z1[t:] - z1[:-t], z2[t:] - z2[:-t]
                vt = (np.sum(diff1**2) + np.sum(diff2**2)) / (2.0 * (half - t))
                rho_t = 1.0 - (vt / (2.0 * var_plus))
                if rho_t < 0.0:
                    break
                sum_rho += rho_t

        tau = 1.0 + 2.0 * sum_rho
        ess = min(float(2 * half), float(2 * half) / tau)
        return float(rhat), float(tau), float(ess)

    def compute_brufit_chain_diagnostics(data):
        arr = np.asarray(data, dtype=float)
        if len(arr) < 4:
            return 0.0, 1.0, float(len(arr))
        if np.all(arr == arr[0]):
            return 1.0, 1.0, float(len(arr))
        z = _vehtari_rank_normalize(arr)
        return _vehtari_diagnostics(z)


def _compute_batch_chain_diagnostics(data: dict, cols: list[str], burn_in: int = 2000) -> dict | None:
    """
    Vectorized tau_int and ESS across all H_* columns using a single weight-index expansion
    and batched 2D NumPy real FFTs (matching BruFit's AutocorrPlot.cpp definition).
    """
    if not cols:
        return None
    n_tree = len(data[cols[0]])
    if n_tree == 0:
        return None

    has_wt = "weight" in data
    if has_wt:
        w_arr = np.asarray(data["weight"], dtype=np.float64)
        int_w = np.maximum(1, np.round(w_arr)).astype(np.int64)
        n_raw = int(np.sum(int_w))
        eff_burn = burn_in if burn_in < n_raw - 50 else min(burn_in, n_raw // 2)
        cum_w = np.cumsum(int_w)
        eff_burn_entry = int(np.searchsorted(cum_w, eff_burn))
    else:
        int_w = None
        n_raw = int(n_tree)
        eff_burn = burn_in if burn_in < n_raw - 50 else min(burn_in, n_raw // 2)
        eff_burn_entry = eff_burn

    n_samples = n_raw - eff_burn
    if n_samples <= 50:
        return None

    if int_w is not None:
        emp_acc = float(n_tree - eff_burn_entry) / float(max(1, n_samples))
    else:
        cum_w = None
        arr0 = np.asarray(data[cols[0]], dtype=np.float64)
        emp_acc = float(np.mean(arr0[eff_burn + 1:] != arr0[eff_burn:-1])) if len(arr0) > eff_burn + 1 else None

    # Bound FFT array length and column chunk size to prevent OOM on ultra-low-acceptance chains
    max_fft_samples = 250_000
    stride = 1 if n_samples <= max_fft_samples else ((n_samples + max_fft_samples - 1) // max_fft_samples)
    if stride == 1:
        if int_w is not None:
            idx = np.repeat(np.arange(n_tree, dtype=np.int64), int_w)[eff_burn:]
        else:
            idx = slice(eff_burn, n_raw)
        n_grid = n_samples
    else:
        sample_pos = np.arange(eff_burn, n_raw, stride, dtype=np.int64)
        n_grid = len(sample_pos)
        if cum_w is not None:
            idx = np.minimum(np.searchsorted(cum_w, sample_pos, side="right"), n_tree - 1)
        else:
            idx = sample_pos

    n_cols = len(cols)
    calc_max_lag = n_grid // 2
    limit = min(n_grid, calc_max_lag + 1)
    n_fft = 1 << int(np.ceil(np.log2(max(2, 2 * n_grid - 1))))

    taus = np.empty(n_cols, dtype=np.float64)
    esses = np.empty(n_cols, dtype=np.float64)
    per_moment = {}
    chunk_size = 16

    for start in range(0, n_cols, chunk_size):
        chunk_cols = cols[start:start + chunk_size]
        c_len = len(chunk_cols)
        mat = np.empty((c_len, n_grid), dtype=np.float64)
        for j, c in enumerate(chunk_cols):
            mat[j] = np.asarray(data[c], dtype=np.float64)[idx]

        means = mat.mean(axis=1, keepdims=True)
        xm = mat - means
        par_d = np.einsum("ij,ij->i", xm, xm)

        fx = np.fft.rfft(xm, n=n_fft, axis=1)
        rx = np.fft.irfft(fx * np.conj(fx), n=n_fft, axis=1)[:, :limit]

        for j, c in enumerate(chunk_cols):
            pd = par_d[j]
            if pd <= 0.0:
                tau = 0.5
            else:
                acf_row = rx[j] / pd
                neg = np.flatnonzero(acf_row[1:] <= 0.0)
                stop = (1 + int(neg[0])) if len(neg) > 0 else limit
                tau_grid = 0.5 + float(np.sum(acf_row[1:stop]))
                tau = tau_grid if stride == 1 else max(0.5, float(stride) * tau_grid)
            ess = (n_samples / (2.0 * tau)) if tau > 0.0 else float(n_samples)
            ess = min(float(n_samples), ess)
            taus[start + j] = tau
            esses[start + j] = ess
            per_moment[c] = {"tau": float(tau), "ess": float(ess)}

    return {
        "n_raw": int(n_raw),
        "burn_in": int(eff_burn),
        "n_samples": int(n_samples),
        "tau_mean": float(np.mean(taus)),
        "tau_max": float(np.max(taus)),
        "ess_min": float(np.min(esses)),
        "ess_mean": float(np.mean(esses)),
        "acceptance": emp_acc,
        "per_moment": per_moment,
    }


def _select_h_columns(all_cols: list[str]) -> list[str]:
    all_h = [str(c) for c in all_cols if str(c).startswith("H_") and not str(c).endswith("_err")]
    strict_h = [c for c in all_h if len(c.split("_")) == 4]
    return strict_h if strict_h else all_h


def compute_bin_diagnostics_from_root(root_file: str, burn_in: int = 2000) -> dict | None:
    """Compute tau_int and ESS from ResultsBruMcmc*.root using fast uproot + vectorized 2D FFT."""
    if not root_file or not os.path.isfile(root_file):
        return None

    if uproot is not None:
        try:
            with uproot.open(root_file) as uf:
                if "MCMCTree" in uf:
                    tree = uf["MCMCTree"]
                    all_cols = [str(k) for k in tree.keys()]
                    cols = _select_h_columns(all_cols)
                    if not cols:
                        return None
                    cols_to_load = list(cols)
                    if "weight" in all_cols:
                        cols_to_load.append("weight")
                    data = tree.arrays(cols_to_load, library="np")
                    return _compute_batch_chain_diagnostics(data, cols, burn_in=burn_in)
        except Exception:
            pass

    root_mod = _get_root()
    if root_mod is None:
        return None

    try:
        rdf = root_mod.RDataFrame("MCMCTree", root_file)
        all_cols = [str(c) for c in rdf.GetColumnNames()]
        cols = _select_h_columns(all_cols)
        if not cols:
            return None

        cols_to_load = list(cols)
        if "weight" in all_cols:
            cols_to_load.append("weight")

        data = rdf.AsNumpy(cols_to_load)
        return _compute_batch_chain_diagnostics(data, cols, burn_in=burn_in)
    except Exception:
        return None


_ROOT_CANDIDATES = (
    "ResultsBruMcmcFixedCovariance.root",
    "ResultsBruMcmcCovariance.root",
    "ResultsBruMcmc.root",
)
_CACHE_FILENAME = ".collector_diag_cache.json"
_CACHE_VERSION = 1


def _stat_signature(path: str) -> tuple[int, int] | None:
    try:
        st = os.stat(path)
        return int(st.st_size), int(getattr(st, "st_mtime_ns", int(st.st_mtime * 1e9)))
    except OSError:
        return None


def _inspect_root_bin_task(task: tuple) -> tuple[int, int, dict]:
    """
    Worker-safe single-pass ROOT file health, projection scaling, and vectorized diagnostics check.
    Executes inside ProcessPoolExecutor workers with zero parent ROOT lock inheritance.
    """
    tbin_idx, b_idx, cand_paths, burn_in, need_diags, need_scaling = task

    root_file = None
    zombie_file = None
    found_dir = None
    diag = None
    sig = None

    for cdir, cand in cand_paths:
        if _fast_binary_zombie_check(cand):
            if not zombie_file:
                zombie_file = cand
                found_dir = cdir
            continue

        opened_ok = False
        has_tree = False
        n_entries = 0
        local_diag = None

        if uproot is not None:
            try:
                with uproot.open(cand) as uf:
                    keys = uf.keys()
                    if len(keys) == 0:
                        opened_ok = False
                    elif "MCMCTree" in uf:
                        tree = uf["MCMCTree"]
                        n_entries = int(tree.num_entries)
                        has_tree = n_entries > 0
                        opened_ok = True
                        if has_tree and need_diags and n_entries > 0:
                            all_cols = [str(k) for k in tree.keys()]
                            cols = _select_h_columns(all_cols)
                            if cols:
                                cols_to_load = list(cols)
                                if "weight" in all_cols:
                                    cols_to_load.append("weight")
                                data = tree.arrays(cols_to_load, library="np")
                                local_diag = _compute_batch_chain_diagnostics(data, cols, burn_in=burn_in)
                    else:
                        opened_ok = True
                        has_tree = False
            except Exception:
                opened_ok = False

        if not opened_ok:
            root_mod = _get_root()
            if root_mod is not None:
                try:
                    old_level = getattr(root_mod, "gErrorIgnoreLevel", 0)
                    root_mod.gErrorIgnoreLevel = getattr(root_mod, "kFatal", 5000)
                    f = root_mod.TFile.Open(cand, "READ")
                    root_mod.gErrorIgnoreLevel = old_level
                    is_z = (not f) or f.IsZombie() or f.GetNkeys() == 0 or f.TestBit(getattr(root_mod.TFile, "kRecovered", 0x00000400))
                    if is_z:
                        if f:
                            f.Close()
                        if not zombie_file:
                            zombie_file = cand
                            found_dir = cdir
                        continue
                    tree = f.Get("MCMCTree")
                    n_entries = int(tree.GetEntries()) if tree else 0
                    has_tree = n_entries > 0
                    f.Close()
                    opened_ok = True
                except Exception:
                    opened_ok = False
            elif uproot is None:
                opened_ok = True
                has_tree = True
                n_entries = 1000

        if not opened_ok:
            if not zombie_file:
                zombie_file = cand
                found_dir = cdir
            continue

        if has_tree:
            root_file = cand
            found_dir = cdir
            diag = local_diag
            sig = _stat_signature(cand)
            break

    scale_ok = True
    scale_ratio = None
    scale_msg = ""
    if root_file and need_scaling:
        scale_ok, scale_ratio, scale_msg = check_projection_scaling(root_file, skip_zombie_check=True)

    if root_file and scale_ok and need_diags and diag is None and uproot is None:
        diag = compute_bin_diagnostics_from_root(root_file, burn_in=burn_in)

    diag_compact = None
    if diag is not None:
        diag_compact = {
            "n_raw": diag.get("n_raw"),
            "n_samples": diag.get("n_samples"),
            "tau_mean": diag.get("tau_mean"),
            "tau_max": diag.get("tau_max"),
            "ess_min": diag.get("ess_min"),
            "ess_mean": diag.get("ess_mean"),
            "acceptance": diag.get("acceptance"),
            "per_moment_ess": {k: float(v["ess"]) for k, v in diag.get("per_moment", {}).items()},
        }

    return tbin_idx, b_idx, {
        "root_file": root_file,
        "zombie_file": zombie_file,
        "found_dir": found_dir,
        "sig": sig,
        "scale_ok": scale_ok,
        "scale_ratio": scale_ratio,
        "scale_msg": scale_msg,
        "diag": diag_compact,
    }


def _scan_tbin_directory(tbin_dir: str, run_dir: str, tbin_name: str, m_min: float, m_max: float, n_bins: int, burn_in: int, target_bin_idx: int | None = None):
    """Single-pass os.scandir + cache load for a t-bin directory."""
    bin_width = (m_max - m_min) / n_bins

    diag_alias = os.path.join(tbin_dir, "MCMC_Diagnostics.txt")
    diag_specific = os.path.join(tbin_dir, f"MomentFit_Diagnostics_{tbin_name}.txt")
    diag_file = None
    diag_mtime = None
    if os.path.isfile(diag_alias) and os.path.isfile(diag_specific):
        m_spec = os.path.getmtime(diag_specific)
        m_alias = os.path.getmtime(diag_alias)
        if m_spec >= m_alias:
            diag_file, diag_mtime = diag_specific, m_spec
        else:
            diag_file, diag_mtime = diag_alias, m_alias
    elif os.path.isfile(diag_specific):
        diag_file, diag_mtime = diag_specific, os.path.getmtime(diag_specific)
    elif os.path.isfile(diag_alias):
        diag_file, diag_mtime = diag_alias, os.path.getmtime(diag_alias)

    precomputed_diags = parse_diagnostics_file(diag_file) if diag_file else {}

    # Single directory scan for all mass-bin subdirs and logRooFit*.txt mtimes
    mass_subdirs = []
    dir_by_name = {}
    tbin_files_mtime = {}
    if os.path.isdir(tbin_dir):
        try:
            with os.scandir(tbin_dir) as it:
                for entry in it:
                    ename = entry.name
                    if entry.is_dir():
                        if ".bak" not in ename:
                            dir_by_name[ename] = entry.path
                            m = re.search(r"Mpi0eta([0-9.]+)", ename)
                            if m:
                                try:
                                    mass_subdirs.append((float(m.group(1)), entry.path))
                                except ValueError:
                                    pass
                    elif ename.startswith("logRooFit") and ename.endswith(".txt"):
                        try:
                            tbin_files_mtime[ename] = entry.stat().st_mtime
                        except OSError:
                            pass
        except OSError:
            pass

    retry_dir = os.path.join(run_dir, "logs")
    retry_files_mtime = {}
    retry_prefix = f"MomentFit_retry_{tbin_name}_bin"
    if os.path.isdir(retry_dir):
        try:
            with os.scandir(retry_dir) as it:
                for entry in it:
                    if entry.name.startswith(retry_prefix) and entry.name.endswith(".log") and entry.is_file():
                        try:
                            retry_files_mtime[entry.name] = entry.stat().st_mtime
                        except OSError:
                            pass
        except OSError:
            pass

    cache_path = os.path.join(tbin_dir, _CACHE_FILENAME)
    cache_data = {"version": _CACHE_VERSION, "bins": {}}
    if os.path.isfile(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as cf:
                loaded = json.load(cf)
            if isinstance(loaded, dict) and loaded.get("version") == _CACHE_VERSION and isinstance(loaded.get("bins"), dict):
                cache_data = loaded
        except Exception:
            pass

    cached_bins = cache_data["bins"]
    bin_indices = [target_bin_idx] if (target_bin_idx is not None and 0 <= target_bin_idx < n_bins) else list(range(n_bins))

    prepared_bins = []
    pending_tasks = []
    for b_idx in bin_indices:
        b_low = m_min + b_idx * bin_width
        b_high = m_min + (b_idx + 1) * bin_width
        b_center = (b_low + b_high) / 2.0
        default_name = f"Mpi0eta{b_center:.6f}_"
        bin_dir_path = os.path.join(tbin_dir, default_name)
        bin_dir_name = default_name

        cand_dirs = [p for mv, p in mass_subdirs if abs(mv - b_center) < 0.005]
        for cname in (
            f"Mpi0eta{b_center:.6f}_", f"Mpi0eta{b_center:.6f}",
            f"Mpi0eta{b_center:.4f}_", f"Mpi0eta{b_center:.4f}",
            f"Mpi0eta{b_center:.2f}_", f"Mpi0eta{b_center:.2f}",
        ):
            cp = dir_by_name.get(cname)
            if cp and cp not in cand_dirs:
                cand_dirs.append(cp)

        if cand_dirs:
            bin_dir_path = cand_dirs[0]
            bin_dir_name = os.path.basename(bin_dir_path)

        cand_files = []
        for cdir in cand_dirs:
            for rname in _ROOT_CANDIDATES:
                cand = os.path.join(cdir, rname)
                if os.path.isfile(cand):
                    cand_files.append((cdir, cand))

        log_name = f"logRooFit{b_idx}.txt"
        log_path = os.path.join(tbin_dir, log_name)
        retry_name = f"{retry_prefix}{b_idx}.log"
        retry_log = os.path.join(retry_dir, retry_name)
        retry_mtime = retry_files_mtime.get(retry_name)
        log_mtime = tbin_files_mtime.get(log_name)
        if retry_mtime is not None and (log_mtime is None or retry_mtime >= log_mtime):
            active_log_path = retry_log
        else:
            active_log_path = log_path

        # Check cache hit for primary candidate file
        cache_entry = cached_bins.get(str(b_idx))
        cache_hit = None
        if cand_files and isinstance(cache_entry, dict) and cache_entry.get("burn_in") == burn_in:
            cached_rname = cache_entry.get("root_name")
            for cdir, cand in cand_files:
                if os.path.basename(cand) == cached_rname and cache_entry.get("bin_dir_name") == os.path.basename(cdir):
                    sig = _stat_signature(cand)
                    if sig is not None and sig[0] == cache_entry.get("size") and sig[1] == cache_entry.get("mtime_ns"):
                        cache_hit = {
                            "root_file": cand,
                            "zombie_file": None,
                            "found_dir": cdir,
                            "sig": sig,
                            "scale_ok": cache_entry.get("scale_ok", True),
                            "scale_ratio": cache_entry.get("scale_ratio"),
                            "scale_msg": cache_entry.get("scale_msg", ""),
                            "diag": cache_entry.get("diag"),
                        }
                    break

        is_diag_stale = False
        if diag_mtime is not None:
            if retry_mtime is not None and retry_mtime > diag_mtime:
                is_diag_stale = True
            elif cand_files:
                sig0 = _stat_signature(cand_files[0][1])
                if sig0 is not None and (sig0[1] / 1e9) > diag_mtime:
                    is_diag_stale = True

        pd = precomputed_diags.get(b_idx) if (b_idx in precomputed_diags and not is_diag_stale) else None
        need_diags = not (
            pd is not None
            and pd.get("status") == "SUCCESS"
            and pd.get("ess_min") is not None
            and pd.get("ess_mean") is not None
            and pd.get("acceptance") is not None
        )

        if cache_hit is None and cand_files:
            pending_tasks.append((b_idx, cand_files, burn_in, need_diags, True))

        prepared_bins.append({
            "b_idx": b_idx,
            "b_low": b_low,
            "b_high": b_high,
            "b_center": b_center,
            "bin_dir_path": bin_dir_path,
            "bin_dir_name": bin_dir_name,
            "cand_files": cand_files,
            "cache_hit": cache_hit,
            "active_log_path": active_log_path,
            "retry_log": retry_log if retry_mtime is not None else None,
            "diag_mtime": diag_mtime,
            "pd": pd,
        })

    return {
        "cache_path": cache_path,
        "cache_data": cache_data,
        "prepared_bins": prepared_bins,
        "pending_tasks": pending_tasks,
    }


def _save_tbin_cache(cache_path: str, cache_data: dict) -> None:
    try:
        cache_dir = os.path.dirname(cache_path)
        fd, tmp_path = tempfile.mkstemp(prefix=".tmp_collector_cache_", suffix=".json", dir=cache_dir)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as tf:
                json.dump(cache_data, tf)
            os.replace(tmp_path, cache_path)
        except Exception:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
            raise
    except OSError:
        pass


def _run_root_inspection_tasks(tasks: list[tuple]) -> list[tuple[int, int, dict]]:
    """Run bin ROOT health + diagnostics tasks in parallel worker processes (with thread fallback)."""
    if not tasks:
        return []
    if len(tasks) <= 2:
        return [_inspect_root_bin_task(t) for t in tasks]

    n_workers = min(len(tasks), max(4, min(os.cpu_count() or 8, 16)))
    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
            return list(executor.map(_inspect_root_bin_task, tasks))
    except Exception:
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
            return list(executor.map(_inspect_root_bin_task, tasks))


def evaluate_tbin(
    tbin_idx: int,
    tbin_name: str,
    tbin_dir: str,
    run_dir: str,
    m_min: float = 0.8,
    m_max: float = 2.0,
    n_bins: int = 30,
    ess_threshold: float = 20.0,
    ess_metric: str = "mean",
    min_ess_threshold: float | None = 10.0,
    burn_in: int = 2000,
    current_scale: float = 1.0,
    current_step: float = 1.0,
    current_iters: int = 10000,
    current_gibbs: int = 0,
    max_l: int = 4,
    target_bin_idx: int | None = None,
    prepared_scan: dict | None = None,
    precomputed_root_results: dict[int, dict] | None = None,
    max_retry_iters: int = MAX_RETRY_ITERS,
    burnin_fraction: float = 0.10,
    low_acceptance_threshold: float = 0.03,
    gibbs_low_acceptance_threshold: float = 0.12,
    high_acceptance_threshold: float = 0.35,
    gibbs_activation_acceptance: float = 0.05,
) -> list[dict]:
    """
    Evaluate mass bins for a given t-bin.
    Uses single-pass os.scandir, atomic per-tbin JSON cache (.collector_diag_cache.json),
    uproot + vectorized 2D FFT chain diagnostics, and multi-process inspection for cache misses.
    """
    if prepared_scan is None:
        prepared_scan = _scan_tbin_directory(
            tbin_dir=tbin_dir,
            run_dir=run_dir,
            tbin_name=tbin_name,
            m_min=m_min,
            m_max=m_max,
            n_bins=n_bins,
            burn_in=burn_in,
            target_bin_idx=target_bin_idx,
        )

    cache_path = prepared_scan["cache_path"]
    cache_data = prepared_scan["cache_data"]
    cached_bins = cache_data["bins"]
    prepared_bins = prepared_scan["prepared_bins"]
    pending_tasks = prepared_scan["pending_tasks"]

    root_results = dict(precomputed_root_results) if precomputed_root_results is not None else {}
    if pending_tasks and precomputed_root_results is None:
        missing_tasks = [(tbin_idx, *task) for task in pending_tasks if task[0] not in root_results]
        for _, res_b, res_dict in _run_root_inspection_tasks(missing_tasks):
            root_results[res_b] = res_dict

    cache_dirty = False
    bin_evals = []
    for pbin in prepared_bins:
        b_idx = pbin["b_idx"]
        b_low = pbin["b_low"]
        b_high = pbin["b_high"]
        b_center = pbin["b_center"]
        bin_dir_path = pbin["bin_dir_path"]
        bin_dir_name = pbin["bin_dir_name"]
        log_path = pbin["active_log_path"]
        retry_log = pbin["retry_log"]

        inspect_res = pbin["cache_hit"]
        if inspect_res is None:
            inspect_res = root_results.get(b_idx)
            if inspect_res is not None and inspect_res.get("root_file") and inspect_res.get("sig"):
                sig = inspect_res["sig"]
                cached_bins[str(b_idx)] = {
                    "root_name": os.path.basename(inspect_res["root_file"]),
                    "bin_dir_name": os.path.basename(inspect_res.get("found_dir") or bin_dir_path),
                    "size": sig[0],
                    "mtime_ns": sig[1],
                    "burn_in": burn_in,
                    "scale_ok": inspect_res.get("scale_ok", True),
                    "scale_ratio": inspect_res.get("scale_ratio"),
                    "scale_msg": inspect_res.get("scale_msg", ""),
                    "diag": inspect_res.get("diag"),
                }
                cache_dirty = True

        root_file = inspect_res.get("root_file") if inspect_res else None
        zombie_file = inspect_res.get("zombie_file") if inspect_res else None
        if inspect_res and inspect_res.get("found_dir"):
            bin_dir_path = inspect_res["found_dir"]
            bin_dir_name = os.path.basename(bin_dir_path)

        parsed_scale, parsed_step, parsed_iters = extract_hyperparameters_from_log(log_path)
        parsed_gibbs = extract_gibbs_from_log(log_path)
        if (parsed_scale is None or parsed_step is None or parsed_iters is None) and retry_log and retry_log != log_path:
            r_scale, r_step, r_iters = extract_hyperparameters_from_log(retry_log)
            if parsed_scale is None: parsed_scale = r_scale
            if parsed_step is None: parsed_step = r_step
            if parsed_iters is None: parsed_iters = r_iters
        if parsed_gibbs is None and retry_log and retry_log != log_path:
            parsed_gibbs = extract_gibbs_from_log(retry_log)

        bin_scale = parsed_scale if parsed_scale is not None else current_scale
        bin_step = parsed_step if parsed_step is not None else current_step
        bin_iters = parsed_iters if parsed_iters is not None else current_iters
        bin_gibbs = parsed_gibbs if parsed_gibbs is not None else current_gibbs

        scale_ok = True
        scale_ratio = None
        scale_msg = ""
        if root_file:
            status = "SUCCESS"
            log_info = {
                "status": "SUCCESS", "phase1_zero_acc": False, "phase2_stuck": False,
                "scale": bin_scale, "step": bin_step, "iters": bin_iters,
                "details": f"{os.path.basename(root_file)} present (MCMCTree detected)"
            }
            scale_ok = bool(inspect_res.get("scale_ok", True))
            scale_ratio = inspect_res.get("scale_ratio")
            scale_msg = inspect_res.get("scale_msg", "")
        elif zombie_file:
            status = "ZOMBIE_ROOT"
            log_aux = analyze_log_content(log_path)
            log_info = {
                "status": "ZOMBIE_ROOT",
                "phase1_zero_acc": log_aux.get("phase1_zero_acc", False),
                "phase2_stuck": log_aux.get("phase2_stuck", False),
                "last_acceptance": log_aux.get("last_acceptance"),
                "scale": bin_scale, "step": bin_step, "iters": bin_iters,
                "details": f"Zombie/corrupted {os.path.basename(zombie_file)} detected (fit interrupted before completion)"
            }
        else:
            log_info = analyze_log_content(log_path)
            status = log_info["status"]
            if parsed_scale is None: bin_scale = log_info.get("scale", current_scale)
            if parsed_step is None: bin_step = log_info.get("step", current_step)
            if parsed_iters is None: bin_iters = log_info.get("iters", current_iters)

        min_ess = None
        mean_ess = None
        mean_tau = None
        max_tau = None
        bin_acceptance = None
        worst_moment = None
        n_low_ess = 0
        total_moments = 0

        pd = pbin["pd"]
        if zombie_file:
            status = "ZOMBIE_ROOT"
            pd = None

        if pd is not None:
            if not (pd.get("status") == "SUCCESS" and not root_file):
                min_ess = pd.get("ess_min")
                mean_ess = pd.get("ess_mean")
                mean_tau = pd.get("tau_mean")
                max_tau = pd.get("tau_max")
                bin_acceptance = pd.get("acceptance")
                worst_moment = pd.get("worst_moment")
                n_low_ess = pd.get("n_low_ess", 0)
                total_moments = pd.get("total_moments", 0)
                if not root_file and pd.get("status"):
                    status = pd["status"]

        if root_file and scale_ok and (min_ess is None or mean_ess is None or bin_acceptance is None):
            diag = inspect_res.get("diag") if inspect_res else None
            if diag is None and uproot is None:
                diag = compute_bin_diagnostics_from_root(root_file, burn_in=burn_in)
            if diag:
                if min_ess is None or mean_ess is None:
                    min_ess = diag.get("ess_min")
                    mean_ess = diag.get("ess_mean")
                    mean_tau = diag.get("tau_mean")
                    max_tau = diag.get("tau_max")
                    per_m_ess = diag.get("per_moment_ess")
                    if per_m_ess is None and isinstance(diag.get("per_moment"), dict):
                        per_m_ess = {k: v["ess"] for k, v in diag["per_moment"].items()}
                    if per_m_ess:
                        total_moments = len(per_m_ess)
                        bad_moments = [(k, val) for k, val in per_m_ess.items() if val < ess_threshold]
                        n_low_ess = len(bad_moments)
                        if bad_moments:
                            bad_moments.sort(key=lambda x: x[1])
                            worst_moment = bad_moments[0][0]
                if bin_acceptance is None:
                    bin_acceptance = diag.get("acceptance")
        if bin_acceptance is None and log_info.get("last_acceptance") is not None:
            bin_acceptance = log_info.get("last_acceptance")

        # Determine effective ESS value according to chosen metric (mean vs min)
        eval_ess = mean_ess if ess_metric == "mean" else min_ess
        if eval_ess is None:
            eval_ess = min_ess if ess_metric == "mean" else mean_ess

        # Classification & Diagnosis
        is_successful = False
        failure_type = None
        diagnosis_reason = ""
        rec_scale = bin_scale
        rec_step = bin_step
        rec_iters = bin_iters

        metric_name = "mean" if ess_metric == "mean" else "min"

        if status in ("ZOMBIE_ROOT", "CORRUPT_ROOT") or (zombie_file and not root_file):
            is_successful = False
            failure_type = "ZOMBIE_ROOT"
            zname = os.path.basename(zombie_file if zombie_file else (root_file or "ROOT file"))
            diagnosis_reason = f"Zombie/corrupted {zname} detected (process interrupted before completion; must be refit)"
            rec_scale = bin_scale
            rec_step = bin_step
            rec_iters = max(bin_iters, current_iters)
        elif root_file and not scale_ok:
            is_successful = False
            failure_type = "NORM_SCALE_MISMATCH"
            diagnosis_reason = (
                f"Negative-PDF normalisation scaling mismatch (Data/Curve ratio = {scale_ratio:.2f} != 1.0). "
                f"The negative intensity penalty guard inflated the normalization integral, indicating walker trapped in unphysical phase space. "
                f"Requires tighter step size or regularized physical initialization."
            )
            rec_scale = max(0.1, round(bin_scale * 0.5, 2))
            rec_step = max(0.05, round(bin_step * 0.5, 3))
            rec_iters = bin_iters
        elif status == "SUCCESS" and eval_ess is not None:
            acc_str = f", acc={bin_acceptance:.1%}" if bin_acceptance is not None else ""
            tau_str = f", tau={mean_tau:.1f}" if mean_tau is not None else ""
            mean_ok = (eval_ess >= ess_threshold)
            min_ok = (
                min_ess_threshold is None
                or min_ess_threshold <= 0.0
                or min_ess is None
                or min_ess >= min_ess_threshold
            )
            if mean_ok and min_ok:
                is_successful = True
                failure_type = "NONE"
                min_str = f", min={min_ess:.1f}" if min_ess is not None else ""
                diagnosis_reason = f"Converged with good {metric_name} ESS ({eval_ess:.1f} >= {ess_threshold:.1f}{min_str}{tau_str}{acc_str})"
            elif not mean_ok:
                is_successful = False
                failure_type = "LOW_ESS"
                rec_scale, rec_step, rec_iters, diagnosis_reason = recommend_low_ess_params(
                    bin_scale=bin_scale,
                    bin_step=bin_step,
                    bin_iters=bin_iters,
                    burn_in=burn_in,
                    eval_ess=eval_ess,
                    ess_threshold=ess_threshold,
                    mean_tau=mean_tau,
                    bin_acceptance=bin_acceptance,
                    metric_name=metric_name,
                    max_iters=max_retry_iters,
                    bin_gibbs=bin_gibbs,
                    low_acceptance_threshold=low_acceptance_threshold,
                    gibbs_low_acceptance_threshold=gibbs_low_acceptance_threshold,
                    high_acceptance_threshold=high_acceptance_threshold,
                )
            else:
                is_successful = False
                failure_type = "LOW_ESS"
                rec_scale, rec_step, rec_iters, diagnosis_reason = recommend_low_ess_params(
                    bin_scale=bin_scale,
                    bin_step=bin_step,
                    bin_iters=bin_iters,
                    burn_in=burn_in,
                    eval_ess=min_ess,
                    ess_threshold=min_ess_threshold,
                    mean_tau=max_tau if (max_tau is not None and max_tau > 1.0) else mean_tau,
                    bin_acceptance=bin_acceptance,
                    metric_name="min",
                    max_iters=max_retry_iters,
                    bin_gibbs=bin_gibbs,
                    low_acceptance_threshold=low_acceptance_threshold,
                    gibbs_low_acceptance_threshold=gibbs_low_acceptance_threshold,
                    high_acceptance_threshold=high_acceptance_threshold,
                )
                worst_note = f" (worst moment: {worst_moment})" if worst_moment else ""
                diagnosis_reason = (
                    f"{metric_name.capitalize()} ESS passed ({eval_ess:.1f} >= {ess_threshold:.1f}), "
                    f"but floor failed{worst_note}: {diagnosis_reason}"
                )
        elif status == "SUCCESS" and eval_ess is None:
            if root_file and (uproot is not None or _get_root() is not None):
                is_successful = False
                failure_type = "CORRUPT_ROOT"
                diagnosis_reason = f"Corrupted or unreadable ROOT file '{os.path.basename(root_file)}' (no valid MCMCTree found; needs refit)"
                rec_scale = bin_scale
                rec_step = bin_step
                rec_iters = bin_iters
            else:
                is_successful = True
                failure_type = "NONE"
                diagnosis_reason = "Posterior covariance written successfully (ESS not evaluated)"
        elif log_info["phase1_zero_acc"] or status == "ZERO_ACC_PH1":
            is_successful = False
            failure_type = "ZERO_ACC_PH1"
            diagnosis_reason = "Zero acceptance in Phase 1 burn-in. Initial moment scale too large or flat likelihood."
            # Rule 1: Scale down initial value scale and step
            rec_scale = max(0.1, round(bin_scale * 0.5, 2))
            rec_step = max(0.1, round(bin_step * 0.5, 2))
            rec_iters = bin_iters
        elif log_info["phase2_stuck"] or status == "STUCK_COVARIANCE":
            is_successful = False
            failure_type = "STUCK_COVARIANCE"
            diagnosis_reason = "MCMC stuck or bailed out in Phase 2 covariance chain. Step size too large for proposal acceptance."
            # Rule 2: Scale down step
            rec_scale = bin_scale
            rec_step = max(0.05, round(bin_step * 0.5, 3))
            rec_iters = bin_iters
        else:
            is_successful = False
            failure_type = "FAILED_OTHER"
            diagnosis_reason = f"Fit incomplete or aborted ({status}: {log_info['details']})"
            rec_scale = max(0.1, round(bin_scale * 0.5, 2))
            rec_step = max(0.05, round(bin_step * 0.5, 3))
            rec_iters = max(10000, bin_iters)

        rec_gibbs = (
            bin_gibbs
            if is_successful
            else recommend_retry_gibbs_block_size(
                max_l=max_l,
                bin_gibbs=bin_gibbs,
                bin_acceptance=bin_acceptance,
                failure_type=failure_type,
                gibbs_activation_acceptance=gibbs_activation_acceptance,
            )
        )
        rec_burnin = burn_in if is_successful else compute_retry_burnin(burn_in, rec_iters, burnin_fraction)

        rec = {
            "tbin_idx": tbin_idx,
            "tbin_name": tbin_name,
            "bin_idx": b_idx,
            "mass_low": b_low,
            "mass_high": b_high,
            "mass_center": b_center,
            "bin_dir_name": bin_dir_name,
            "bin_dir_path": bin_dir_path,
            "root_file": root_file,
            "zombie_file": zombie_file,
            "log_path": log_path,
            "status": status,
            "is_successful": is_successful,
            "failure_type": failure_type,
            "diagnosis_reason": diagnosis_reason,
            "min_ess": min_ess,
            "mean_ess": mean_ess,
            "eval_ess": eval_ess,
            "ess_metric": ess_metric,
            "mean_tau": mean_tau,
            "acceptance": bin_acceptance,
            "worst_moment": worst_moment,
            "n_low_ess": n_low_ess,
            "total_moments": total_moments,
            "current_scale": bin_scale,
            "current_step": bin_step,
            "current_iters": bin_iters,
            "current_gibbs": bin_gibbs,
            "rec_scale": rec_scale,
            "rec_step": rec_step,
            "rec_iters": rec_iters,
            "rec_burnin": rec_burnin,
            "rec_gibbs": rec_gibbs,
            "scale_ratio": scale_ratio
        }
        bin_evals.append(rec)

    if cache_dirty:
        _save_tbin_cache(cache_path, cache_data)

    return bin_evals


def print_evaluation_summary(
    evaluations: list[dict],
    ess_threshold: float,
    ess_metric: str = "mean",
    min_ess_threshold: float | None = None,
):
    """Print a clean diagnostic and recommendation table to the console."""
    total_bins = len(evaluations)
    successful_bins = [b for b in evaluations if b["is_successful"]]
    unsuccessful_bins = [b for b in evaluations if not b["is_successful"]]
    metric_str = "Mean ESS" if ess_metric == "mean" else "Min ESS"
    floor_str = (
        f", Min ESS >= {min_ess_threshold:.1f}"
        if (min_ess_threshold is not None and min_ess_threshold > 0.0 and ess_metric == "mean")
        else ""
    )

    print("\n" + "=" * 124)
    print(f" MOMENTFIT DIAGNOSTICS & SUCCESS EVALUATION (Threshold: {metric_str} >= {ess_threshold:.1f}{floor_str})")
    print(f" Total Bins Evaluated: {total_bins} | Succeeded: {len(successful_bins)} | Unsuccessful: {len(unsuccessful_bins)}")
    print("=" * 124)

    if not unsuccessful_bins:
        print(f"  ✅ All bins across all t-bins successfully converged with {metric_str} >= {ess_threshold:.1f}{floor_str}! No retries needed.")
        print("=" * 124 + "\n")
        return

    # Group by tbin
    tbins_present = sorted(list(set(b["tbin_name"] for b in evaluations)))
    for t_name in tbins_present:
        t_unsucc = [b for b in unsuccessful_bins if b["tbin_name"] == t_name]
        t_all = [b for b in evaluations if b["tbin_name"] == t_name]
        print(f"\n>>> t-bin: {t_name} ({len(t_unsucc)} / {len(t_all)} bins require retry)")
        if not t_unsucc:
            print(f"    ✅ All mass bins converged with good {metric_str}!")
            continue

        header = f"{'Bin':<5} {'Mass (GeV)':<12} {'Status / Failure':<18} {'Mean ESS':<10} {'Min ESS':<10} {'Tau':<8} {'Acc':<8} {'Curr (Scl/Stp)':<16} {'Recommended (Scl/Stp/Iter)':<26}"
        print("-" * 124)
        print(header)
        print("-" * 124)

        for b in t_unsucc:
            m_str = f"{b['mass_center']:.4f}"
            mean_ess_str = f"{b['mean_ess']:.1f}" if b.get("mean_ess") is not None else "-"
            min_ess_str = f"{b['min_ess']:.1f}" if b.get("min_ess") is not None else "-"
            tau_str = f"{b['mean_tau']:.1f}" if b.get("mean_tau") is not None else "-"
            acc_str = f"{b['acceptance']:.1%}" if b.get("acceptance") is not None else "-"
            curr_str = f"{b['current_scale']:.2f} / {b['current_step']:.3f}"
            gibbs_tag = f" / G{b.get('rec_gibbs', 0)}" if b.get("rec_gibbs", 0) > 0 else ""
            rec_str = f"{b['rec_scale']:.2f} / {b['rec_step']:.3f} / {b['rec_iters']}{gibbs_tag}"
            fail_str = b["failure_type"]
            print(f"{b['bin_idx']:<5} {m_str:<12} {fail_str:<18} {mean_ess_str:<10} {min_ess_str:<10} {tau_str:<8} {acc_str:<8} {curr_str:<16} {rec_str:<26}")

        print("-" * 124)
        print("Detailed Diagnoses & Action Plan:")
        for b in t_unsucc:
            print(f"  * Bin {b['bin_idx']} (m={b['mass_center']:.3f} GeV): {b['diagnosis_reason']}")
            cur_g = b.get("current_gibbs", 0)
            rec_g = b.get("rec_gibbs", 0)
            gibbs_note = f", enable blocked Gibbs (block={rec_g})" if (rec_g != cur_g and rec_g > 0) else (
                f", keep Gibbs block={rec_g}" if rec_g > 0 else ""
            )
            if b['failure_type'] == 'LOW_ESS':
                if b['rec_step'] > b['current_step']:
                    print(f"    -> Action: High acceptance -> INCREASE step size ({b['current_step']:.3f} -> {b['rec_step']:.3f}) and scale iterations ({b['current_iters']} -> {b['rec_iters']}){gibbs_note}.")
                elif b['rec_step'] < b['current_step']:
                    print(f"    -> Action: Low acceptance -> DECREASE step size ({b['current_step']:.3f} -> {b['rec_step']:.3f}) and scale iterations ({b['current_iters']} -> {b['rec_iters']}){gibbs_note}.")
                else:
                    print(f"    -> Action: Healthy acceptance -> KEEP step size ({b['current_step']:.3f}), SCALE iterations ({b['current_iters']} -> {b['rec_iters']}){gibbs_note} to reach target ESS.")
            elif b['failure_type'] == 'ZERO_ACC_PH1':
                print(f"    -> Action: SCALE DOWN initial value scale ({b['current_scale']:.2f} -> {b['rec_scale']:.2f}) to find valid likelihood basin in burn-in.")
            elif b['failure_type'] == 'STUCK_COVARIANCE':
                print(f"    -> Action: SCALE DOWN step size ({b['current_step']:.3f} -> {b['rec_step']:.3f}){gibbs_note} to improve proposal acceptance in covariance space.")
            elif b['failure_type'] in ('ZOMBIE_ROOT', 'CORRUPT_ROOT'):
                print(f"    -> Action: REFIT corrupted/zombie mass bin (prior interrupted directory backed up; running full sampling).")
            elif b['failure_type'] == 'NORM_SCALE_MISMATCH':
                print(f"    -> Action: NORM SCALE MISMATCH (Data/Curve = {b.get('scale_ratio', 0):.2f}): trapped in negative-intensity phase space. Decrease step/scale or regularize linear algebra moments (scale={b['rec_scale']:.2f}, step={b['rec_step']:.3f}).")
            else:
                print(f"    -> Action: Safe parameter scaling (scale={b['rec_scale']:.2f}, step={b['rec_step']:.3f}, iters={b['rec_iters']}).")

    print("\n" + "=" * 124 + "\n")


def evaluate_single_job(
    run_dir: str,
    tbin_str: str,
    bin_idx: int,
    m_min: float = 0.8,
    m_max: float = 2.0,
    n_bins: int = 30,
    ess_threshold: float = 20.0,
    ess_metric: str = "mean",
    min_ess_threshold: float | None = 10.0,
    burn_in: int = 2000,
    current_scale: float = 1.0,
    current_step: float = 1.0,
    current_iters: int = 10000,
    current_gibbs: int = 0,
    max_l: int = 4,
    max_retry_iters: int = MAX_RETRY_ITERS,
    burnin_fraction: float = 0.10,
    low_acceptance_threshold: float = 0.03,
    gibbs_low_acceptance_threshold: float = 0.12,
    high_acceptance_threshold: float = 0.35,
    gibbs_activation_acceptance: float = 0.05,
) -> int:
    """
    Evaluate outcome for a single completed refit job and print clean status banner.
    Returns:
      0: SUCCESS with Good ESS (>= ess_threshold)
      1: SUCCESS with Low ESS (< ess_threshold)
      2: STUCK in 1DStep (Phase 1)
      3: STUCK in Covariance (Phase 2)
      4: FAILED / OTHER
    """
    discovered_tbins = discover_tbin_dirs(run_dir, tbin_filter=tbin_str)
    tbin_dir = None
    tbin_name = None
    tbin_idx = 0
    if discovered_tbins:
        tbin_idx, tbin_name, tbin_dir = discovered_tbins[0]
    else:
        cand_dir = os.path.join(run_dir, tbin_str)
        if os.path.isdir(cand_dir):
            tbin_dir = cand_dir
            tbin_name = os.path.basename(cand_dir)
            tbin_idx = 0

    if not tbin_dir or not os.path.isdir(tbin_dir):
        print(f"  >>> [Outcome] ❌ ERROR: t-bin directory '{tbin_str}' not found in {run_dir}")
        return 4

    bin_evals = evaluate_tbin(
        tbin_idx=tbin_idx,
        tbin_name=tbin_name,
        tbin_dir=tbin_dir,
        run_dir=run_dir,
        m_min=m_min,
        m_max=m_max,
        n_bins=n_bins,
        ess_threshold=ess_threshold,
        ess_metric=ess_metric,
        min_ess_threshold=min_ess_threshold,
        burn_in=burn_in,
        current_scale=current_scale,
        current_step=current_step,
        current_iters=current_iters,
        current_gibbs=current_gibbs,
        max_l=max_l,
        target_bin_idx=bin_idx,
        max_retry_iters=max_retry_iters,
        burnin_fraction=burnin_fraction,
        low_acceptance_threshold=low_acceptance_threshold,
        gibbs_low_acceptance_threshold=gibbs_low_acceptance_threshold,
        high_acceptance_threshold=high_acceptance_threshold,
        gibbs_activation_acceptance=gibbs_activation_acceptance,
    )

    matching = [b for b in bin_evals if b["bin_idx"] == bin_idx]
    if not matching:
        print(f"  >>> [Outcome] ❌ ERROR: Bin {bin_idx} evaluation failed (index out of range).")
        return 4

    b = matching[0]
    m_center = b["mass_center"]
    mean_ess = b.get("mean_ess")
    min_ess = b.get("min_ess")
    mean_tau = b.get("mean_tau")
    acc_val = b.get("acceptance")
    fail_type = b["failure_type"]
    is_success = b["is_successful"]
    diagnosis = b.get("diagnosis_reason", "")

    if is_success or fail_type == "NONE":
        print(f"  >>> [Outcome] ✅ SUCCESS with Good ESS")
        ess_parts = []
        if mean_ess is not None:
            ess_parts.append(f"Mean ESS = {mean_ess:.1f} >= {ess_threshold:.1f}")
        if min_ess is not None:
            ess_parts.append(f"Min ESS = {min_ess:.1f}")
        if mean_tau is not None:
            ess_parts.append(f"Tau = {mean_tau:.1f}")
        if acc_val is not None:
            ess_parts.append(f"Acc = {acc_val:.1%}")
        if ess_parts:
            print(f"      Metrics: {', '.join(ess_parts)}")
        print(f"      Status:  Chain converged and met publication threshold ({ess_metric} ESS >= {ess_threshold:.1f}).")
        return 0

    elif fail_type == "LOW_ESS":
        print(f"  >>> [Outcome] ⚠️  SUCCESS with Low ESS")
        ess_parts = []
        if mean_ess is not None:
            ess_parts.append(f"Mean ESS = {mean_ess:.1f} < {ess_threshold:.1f}")
        if min_ess is not None:
            ess_parts.append(f"Min ESS = {min_ess:.1f}")
        if mean_tau is not None:
            ess_parts.append(f"Tau = {mean_tau:.1f}")
        if acc_val is not None:
            ess_parts.append(f"Acc = {acc_val:.1%}")
        if ess_parts:
            print(f"      Metrics: {', '.join(ess_parts)}")
        print(f"      Status:  Fit converged, but high autocorrelation / slow exploration.")
        if b['rec_step'] > b['current_step']:
            print(f"      Action:  High acceptance ({acc_val:.1%}) -> Recommend increasing step ({b['current_step']:.3f} -> {b['rec_step']:.3f}) and iters ({b['current_iters']} -> {b['rec_iters']}).")
        elif b['rec_step'] < b['current_step']:
            print(f"      Action:  Low acceptance ({acc_val:.1%}) -> Recommend decreasing step ({b['current_step']:.3f} -> {b['rec_step']:.3f}) and iters ({b['current_iters']} -> {b['rec_iters']}).")
        else:
            print(f"      Action:  Healthy acceptance ({acc_val:.1%}) -> Recommend scaling iterations ({b['current_iters']} -> {b['rec_iters']}) at step {b['rec_step']:.3f}.")
        return 1

    elif fail_type == "ZERO_ACC_PH1":
        print(f"  >>> [Outcome] ❌ STUCK in 1DStep")
        print(f"      Details: Phase 1 initial burn-in failed with near-zero acceptance (initial moments out of likelihood basin).")
        print(f"      Action:  Recommend retrying with smaller initial scale ({b['current_scale']:.2f} -> {b['rec_scale']:.2f}) and step ({b['current_step']:.3f} -> {b['rec_step']:.3f}).")
        return 2

    elif fail_type == "STUCK_COVARIANCE":
        print(f"  >>> [Outcome] ❌ STUCK in Covariance")
        print(f"      Details: Phase 2 proposal rejected in covariance space / chain bailed out early (acceptance near 0).")
        print(f"      Action:  Recommend retrying with smaller step size ({b['current_step']:.3f} -> {b['rec_step']:.4f}).")
        return 3

    elif fail_type in ("ZOMBIE_ROOT", "CORRUPT_ROOT"):
        print(f"  >>> [Outcome] ❌ ZOMBIE / CORRUPTED ROOT FILE")
        print(f"      Details: {diagnosis}")
        print(f"      Action:  Refit required (backup prior interrupted folder and refit with iters={b['rec_iters']}).")
        return 4

    else:
        print(f"  >>> [Outcome] ❌ FAILED / BAILED OUT")
        print(f"      Details: {diagnosis}")
        print(f"      Action:  Recommend safe scaling (scale={b['rec_scale']:.2f}, step={b['rec_step']:.4f}, iters={b['rec_iters']}).")
        return 4


def generate_retry_script(
    run_dir: str,
    unsuccessful_bins: list[dict],
    macro_path: str,
    max_l: int = 4,
    m_min: float = 0.8,
    m_max: float = 2.0,
    n_bins: int = 30,
    burn_in: int = 2000,
    ess_threshold: float = 20.0,
    min_ess_threshold: float | None = 10.0,
    max_parallel: int | None = None,
    output_script_path: str | None = None,
    cov_dir: str = "",
    data_dir: str = "",
    apply_mc_variance: bool = False,
    allow_negative_pdf: bool = True,
    use_tuning: bool = True,
    cold_warmup_acc: float = 0.15,
    data_pattern: str = "",
    mcrecon_pattern: str = "",
    mcgen_pattern: str = "",
    weight_pattern: str = "",
    ess_metric: str = "mean",
    min_acceptance: float = 0.05,
    max_acceptance: float = 0.15,
    target_acceptance: float | None = None,
    tune_mode: str = "MappedRhat",
    rhat_window: int = 1000,
    rhat_target: float = 1.05,
    ess_fraction: float = 0.02,
    rhat_max_retries: int = 5,
    bailout_threshold: int = 10000,
    bailout_min_acceptance: float = 0.0001,
    mc_integration_bins: int = 50,
    fixed_covariance: bool = True,
    use_covariance: bool = True,
    set_init_value: bool = True,
    use_gj_frame: bool = True,
    num_cpu: int = 1,
    seed: int = 111,
    max_retry_iters: int = MAX_RETRY_ITERS,
    burnin_fraction: float = 0.10,
    low_acceptance_threshold: float = 0.03,
    gibbs_low_acceptance_threshold: float = 0.12,
    high_acceptance_threshold: float = 0.35,
    gibbs_activation_acceptance: float = 0.05,
    recondition_covariance: bool = False,
    recondition_shrinkage: float = 0.2,
    recondition_min_eigen_floor: float = 1e-4,
) -> str:
    """Generate a self-contained shell script that refits only the unsuccessful mass bins."""
    if output_script_path is None:
        output_script_path = os.path.join(run_dir, "RunMomentFit_retry.sh")

    effective_max_parallel = max_parallel if max_parallel is not None else DEFAULT_MAX_PARALLEL
    logs_dir = os.path.join(run_dir, "logs")

    script_lines = [
        "#!/bin/bash",
        "# ==============================================================================",
        f"# Automatically generated MomentFit Retry Script",
        f"# Target Run Directory: {run_dir}",
        f"# Unsuccessful Bins to Refit: {len(unsuccessful_bins)}",
        f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "# ==============================================================================",
        "set -e",
        "",
        "# 1. Setup BruFit & Analysis Environment (skip if already set up)",
        'if [ -z "$BRUFIT" ]; then',
        '    CURRENT_DIR="$(pwd)"',
        '    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        '    if [ -f "$SCRIPT_DIR/../../set_env.sh" ]; then',
        '        source "$SCRIPT_DIR/../../set_env.sh"',
        '    elif [ -f "$SCRIPT_DIR/../set_env.sh" ]; then',
        '        source "$SCRIPT_DIR/../set_env.sh"',
        '    elif [ -f "$SCRIPT_DIR/../../set_ana_env.sh" ]; then',
        '        source "$SCRIPT_DIR/../../set_ana_env.sh"',
        '    elif [ -f "$SCRIPT_DIR/../set_ana_env.sh" ]; then',
        '        source "$SCRIPT_DIR/../set_ana_env.sh"',
        '    elif [ -f "/d/home/septian/EtaPi0Analysis/set_ana_env.sh" ]; then',
        '        source "/d/home/septian/EtaPi0Analysis/set_ana_env.sh"',
        '    elif [ -f "/d/grid17/septian/EtaPi0Analysis/set_ana_env.sh" ]; then',
        '        source "/d/grid17/septian/EtaPi0Analysis/set_ana_env.sh"',
        '    fi',
        '    cd "$CURRENT_DIR" 2>/dev/null || true',
        'fi',
        "",
        'if ! command -v root &>/dev/null; then',
        '    if [ -f "/d/grid17/septian/root_v6-36-08_cuda_install/bin/thisroot.sh" ]; then',
        '        export PATH="/d/grid17/septian/root_v6-36-08_cuda_install/bin:$PATH"',
        '        export LD_LIBRARY_PATH="/d/grid17/septian/root_v6-36-08_cuda_install/lib:$LD_LIBRARY_PATH"',
        '        source "/d/grid17/septian/root_v6-36-08_cuda_install/bin/thisroot.sh"',
        '    fi',
        'fi',
        "",
        f'PYTHON_BIN="{os.path.join(repo_root, ".pixi/envs/default/bin/python")}"',
        '[ ! -f "$PYTHON_BIN" ] && PYTHON_BIN="/d/grid17/septian/EtaPi0Analysis/.pixi/envs/default/bin/python"',
        '[ ! -f "$PYTHON_BIN" ] && PYTHON_BIN="/d/home/septian/.pixi/envs/default/bin/python"',
        '[ ! -f "$PYTHON_BIN" ] && PYTHON_BIN="python3"',
        "",
        f'RUN_DIR="{run_dir}"',
        f'LOGS_DIR="{logs_dir}"',
        f'MACRO_PATH="{macro_path}"',
        f'MONITOR_SCRIPT="{os.path.join(repo_root, "macros/MonitorMomentFit.py")}"',
        '[ ! -f "$MONITOR_SCRIPT" ] && MONITOR_SCRIPT="/d/home/septian/EtaPi0Analysis/macros/MonitorMomentFit.py"',
        f'COLLECTOR_SCRIPT="{os.path.join(repo_root, "macros/MomentFitCollector.py")}"',
        '[ ! -f "$COLLECTOR_SCRIPT" ] && COLLECTOR_SCRIPT="/d/home/septian/EtaPi0Analysis/macros/MomentFitCollector.py"',
        "",
        'mkdir -p "$LOGS_DIR"',
        "",
        '# Concurrency limit (supports optional numeric $1 or environment variable MAX_PARALLEL; preserves remaining flags in "$@")',
        'NUM_CPUS=$(nproc 2>/dev/null || grep -c ^processor /proc/cpuinfo 2>/dev/null || echo 2)',
        'DEFAULT_PARALLEL=$(( NUM_CPUS / 2 ))',
        '[ "$DEFAULT_PARALLEL" -lt 1 ] && DEFAULT_PARALLEL=1',
        f'MAX_PARALLEL="${{MAX_PARALLEL:-{"$DEFAULT_PARALLEL" if max_parallel is None else max_parallel}}}"',
        'if [[ "${1:-}" =~ ^[0-9]+$ ]]; then',
        '    MAX_PARALLEL="$1"',
        '    shift',
        'fi',
        'TIMESTAMP=$(date +"%Y%m%d_%H%M%S")',
        "",
        'echo "=========================================================="',
        'echo " Starting MomentFit Targeted Refit for Unsuccessful Bins"',
        f'echo " Total Bins to Refit: {len(unsuccessful_bins)}"',
        'echo " Concurrency Limit:   $MAX_PARALLEL jobs"',
        'echo "=========================================================="',
        "",
        "# 2. Backup prior bin folders to preserve history before refitting",
        'echo ""',
        'echo "[1/4] Creating backups of prior failed bin folders..."'
    ]

    # Add backup commands
    backed_up = set()
    for b in unsuccessful_bins:
        b_path = b["bin_dir_path"]
        if b_path not in backed_up:
            backed_up.add(b_path)
            script_lines.append(f'if [ -d "{b_path}" ]; then')
            script_lines.append(f'    cp -r "{b_path}" "{b_path}.bak_${{TIMESTAMP}}" 2>/dev/null || true')
            script_lines.append('fi')

    script_lines.extend([
        "",
        "# 3. Dynamic rolling job execution pool for targeted refits",
        'echo ""',
        'echo "[2/4] Launching targeted single-bin refits..."',
        "RUNNING_JOBS=0",
        "FAILED_JOBS=0",
        "declare -A PID_TO_JOB",
        "",
        "reap_one_job() {",
        "    local finished_pid",
        "    if wait -n -p finished_pid 2>/dev/null; then",
        '        local desc="${PID_TO_JOB[$finished_pid]:-PID $finished_pid}"',
        '        echo "  [Slot Freed] ✅ Finished: $desc"',
        "    else",
        "        FAILED_JOBS=$((FAILED_JOBS + 1))",
        '        local desc="${PID_TO_JOB[$finished_pid]:-PID $finished_pid}"',
        '        echo "  [Slot Freed] ⚠️  Finished with bailout/error: $desc"',
        "    fi",
        '    [ -n "$finished_pid" ] && unset "PID_TO_JOB[$finished_pid]"',
        "    RUNNING_JOBS=$((RUNNING_JOBS - 1))",
        "}",
        "",
        "wait_for_slot() {",
        '    while [ "$RUNNING_JOBS" -ge "$MAX_PARALLEL" ]; do',
        "        reap_one_job",
        "    done",
        "}",
        "",
        "cleanup() {",
        '    echo ""',
        '    echo "Caught interrupt signal! Terminating running refit jobs..."',
        '    for pid in "${!PID_TO_JOB[@]}"; do',
        '        kill -TERM "$pid" 2>/dev/null || true',
        '    done',
        "    exit 130",
        "}",
        "trap cleanup SIGINT SIGTERM",
        ""
    ])

    # Add each retry command
    mc_var_str = "true" if apply_mc_variance else "false"
    neg_pdf_str = "true" if allow_negative_pdf else "false"
    tuning_str = "true" if use_tuning else "false"
    eff_target_acc = (
        target_acceptance
        if target_acceptance is not None
        else 0.5 * (float(min_acceptance) + float(max_acceptance))
    )

    script_lines.extend([
        "# Export MomentFit.C configuration environment overrides",
        f'export MOM2AMP_MIN_ACCEPTANCE="{min_acceptance}"',
        f'export MOM2AMP_MAX_ACCEPTANCE="{max_acceptance}"',
        f'export MOM2AMP_TARGET_ACCEPTANCE="{eff_target_acc}"',
        f'export MOM2AMP_TUNE_MODE="{tune_mode}"',
        f'export MOM2AMP_RHAT_WINDOW="{rhat_window}"',
        f'export MOM2AMP_RHAT_TARGET="{rhat_target}"',
        f'export MOM2AMP_ESS_FRACTION="{ess_fraction}"',
        f'export MOM2AMP_RHAT_MAX_RETRIES="{rhat_max_retries}"',
        f'export MOM2AMP_BAILOUT_THRESHOLD="{bailout_threshold}"',
        f'export MOM2AMP_BAILOUT_MIN_ACC="{bailout_min_acceptance}"',
        f'export MOM2AMP_MC_INT_BINS="{mc_integration_bins}"',
        f'export MOM2AMP_FIXED_COV="{"true" if fixed_covariance else "false"}"',
        f'export MOM2AMP_USE_COV="{"true" if use_covariance else "false"}"',
        f'export MOM2AMP_SET_INIT_VALUE="{"true" if set_init_value else "false"}"',
        f'export MOM2AMP_USE_GJ_FRAME="{"true" if use_gj_frame else "false"}"',
        f'export MOM2AMP_NUM_CPU="{num_cpu}"',
        f'export MOM2AMP_SEED="{seed}"',
        "",
    ])

    for idx, b in enumerate(unsuccessful_bins, 1):
        t_idx = b["tbin_idx"]
        t_name = b["tbin_name"]
        bin_idx = b["bin_idx"]
        scale = b["rec_scale"]
        step = b["rec_step"]
        iters = b["rec_iters"]
        bin_burnin = b.get("rec_burnin", burn_in)
        gibbs = b.get("rec_gibbs", 0)
        bin_cov_dir = b.get("rec_cov_dir") or cov_dir
        m_center = b["mass_center"]
        fail_type = b["failure_type"]
        log_f = f"${{LOGS_DIR}}/MomentFit_retry_{t_name}_bin{bin_idx}.log"

        data_f = _resolve_pattern_for_tbin(data_pattern, t_name, data_dir)
        mcrecon_f = _resolve_pattern_for_tbin(mcrecon_pattern, t_name, data_dir)
        mcgen_f = _resolve_pattern_for_tbin(mcgen_pattern, t_name, data_dir)
        weight_f = _resolve_pattern_for_tbin(weight_pattern, t_name, data_dir)

        cmd_root = (
            f'root -l -b -q "${{BRUFIT}}/macros/LoadBru.C" '
            f'"{macro_path}({t_idx}, {max_l}, {scale}, {step}, {bin_idx}, {iters}, {bin_burnin}, '
            f'{m_min}, {m_max}, {n_bins}, \\"${{RUN_DIR}}\\", \\"{bin_cov_dir}\\", \\"{data_dir}\\", '
            f'{mc_var_str}, {neg_pdf_str}, {tuning_str}, {gibbs}, {cold_warmup_acc}, '
            f'\\"{data_f}\\", \\"{mcrecon_f}\\", \\"{mcgen_f}\\", \\"{weight_f}\\")" > "{log_f}" 2>&1'
        )
        gibbs_msg = f", gibbs={gibbs}" if gibbs > 0 else ""
        recond_msg = ", reconditioned_cov" if b.get("rec_cov_dir") else ""

        script_lines.extend([
            "wait_for_slot",
            f'echo "  [{idx}/{len(unsuccessful_bins)}] Refitting {t_name} Bin {bin_idx} (m={m_center:.3f} GeV, {fail_type}) -> scale={scale}, step={step}, iters={iters}, burnin={bin_burnin}{gibbs_msg}{recond_msg}"',
            f'( ',
            f'    set +e',
            f'    {cmd_root}',
            f'    JOB_EXIT=$?',
            f'    cp -f "{log_f}" "${{RUN_DIR}}/{t_name}/logRooFit{bin_idx}.txt" 2>/dev/null || true',
            f'    echo ""',
            f'    echo "  ----------------------------------------------------------------------"',
            f'    echo "  [Refit Completed] {t_name} Bin {bin_idx} (m={m_center:.3f} GeV)"',
            f'    "$PYTHON_BIN" "$COLLECTOR_SCRIPT" \\',
            f'        --run-dir "$RUN_DIR" \\',
            f'        --eval-job "{t_name}" "{bin_idx}" \\',
            f'        --max-l "{max_l}" \\',
            f'        --init-scale "{scale}" \\',
            f'        --step "{step}" \\',
            f'        --iters "{iters}" \\',
            f'        --gibbs-block-size "{gibbs}" \\',
            f'        --ess-threshold "{ess_threshold}" \\',
            f'        --ess-metric "{ess_metric}" \\',
            f'        --min-ess-threshold "{min_ess_threshold if min_ess_threshold is not None else 0}" \\',
            f'        --burn-in "{bin_burnin}" \\',
            f'        --max-retry-iters "{max_retry_iters}" \\',
            f'        --burnin-fraction "{burnin_fraction}" \\',
            f'        --low-acceptance-threshold "{low_acceptance_threshold}" \\',
            f'        --gibbs-low-acceptance-threshold "{gibbs_low_acceptance_threshold}" \\',
            f'        --high-acceptance-threshold "{high_acceptance_threshold}" \\',
            f'        --gibbs-activation-acceptance "{gibbs_activation_acceptance}" \\',
            f'        --macro-path "$MACRO_PATH" 2>/dev/null || true',
            f'    echo "  ----------------------------------------------------------------------"',
            f'    exit $JOB_EXIT',
            f') &',
            'JOB_PID=$!',
            f'PID_TO_JOB[$JOB_PID]="{t_name} Bin {bin_idx} (m={m_center:.3f} GeV)"',
            "RUNNING_JOBS=$((RUNNING_JOBS + 1))",
            ""
        ])

    script_lines.extend([
        "# Wait for remaining running jobs in the pool",
        'echo ""',
        'echo "Waiting for remaining $RUNNING_JOBS refit job(s) in pool to complete..."',
        'while [ "$RUNNING_JOBS" -gt 0 ]; do',
        '    reap_one_job',
        'done',
        "",
        'if [ "$FAILED_JOBS" -gt 0 ]; then',
        '    echo "⚠️  $FAILED_JOBS refit job(s) completed with non-zero exit code or fast bailout."',
        'else',
        '    echo "✅ All single-bin refits completed successfully!"',
        'fi',
        "",
        "# 4. Re-run post-fit monitoring and acceptance verification",
        'echo ""',
        'echo "[3/4] Re-running monitoring and diagnostic pipeline..."',
        '"$PYTHON_BIN" "$MONITOR_SCRIPT" \\',
        '    --run-dir "$RUN_DIR" \\',
        '    --tbin "all" \\',
        f'    --ess-threshold "{ess_threshold}" \\',
        f'    --max-l "{max_l}"',
        "",
        "# 5. Re-evaluate fit quality across all mass bins and check convergence",
        'echo ""',
        'echo "[4/4] Evaluating overall fit quality and convergence..."',
        '"$PYTHON_BIN" "$COLLECTOR_SCRIPT" \\',
        '    --run-dir "$RUN_DIR" \\',
        f'    --max-l "{max_l}" \\',
        f'    --ess-threshold "{ess_threshold}" \\',
        f'    --ess-metric "{ess_metric}" \\',
        f'    --min-ess-threshold "{min_ess_threshold if min_ess_threshold is not None else 0}" \\',
        f'    --burn-in "{burn_in}" \\',
        f'    --cov-dir "{cov_dir}" \\',
        f'    --data-dir "{data_dir}" \\',
        f'    --cold-warmup-acc "{cold_warmup_acc}" \\',
        f'    --min-acceptance "{min_acceptance}" \\',
        f'    --max-acceptance "{max_acceptance}" \\',
        f'    --target-acceptance "{eff_target_acc}" \\',
        f'    --tune-mode "{tune_mode}" \\',
        f'    --rhat-window "{rhat_window}" \\',
        f'    --rhat-target "{rhat_target}" \\',
        f'    --ess-fraction "{ess_fraction}" \\',
        f'    --rhat-max-retries "{rhat_max_retries}" \\',
        f'    --bailout-threshold "{bailout_threshold}" \\',
        f'    --bailout-min-acceptance "{bailout_min_acceptance}" \\',
        f'    --mc-integration-bins "{mc_integration_bins}" \\',
        f'    --num-cpu "{num_cpu}" \\',
        f'    --seed "{seed}" \\',
        f'    --max-retry-iters "{max_retry_iters}" \\',
        f'    --burnin-fraction "{burnin_fraction}" \\',
        f'    --low-acceptance-threshold "{low_acceptance_threshold}" \\',
        f'    --gibbs-low-acceptance-threshold "{gibbs_low_acceptance_threshold}" \\',
        f'    --high-acceptance-threshold "{high_acceptance_threshold}" \\',
        f'    --gibbs-activation-acceptance "{gibbs_activation_acceptance}" \\',
        f'    {"--mc-variance" if apply_mc_variance else "--no-mc-variance"} \\',
        f'    {"--negative-pdf" if allow_negative_pdf else "--no-negative-pdf"} \\',
        f'    {"--tuning" if use_tuning else "--no-tuning"} \\',
        f'    {"--fixed-covariance" if fixed_covariance else "--no-fixed-covariance"} \\',
        f'    {"--use-covariance" if use_covariance else "--no-use-covariance"} \\',
        f'    {"--set-init-value" if set_init_value else "--no-set-init-value"} \\',
        f'    {"--gj-frame" if use_gj_frame else "--hel-frame"} \\',
        f'    {"--recondition-covariance" if recondition_covariance else "--no-recondition-covariance"} \\',
        f'    --recondition-shrinkage "{recondition_shrinkage}" \\',
        f'    --recondition-min-eigen-floor "{recondition_min_eigen_floor}" \\',
        f'    --data-pattern "{data_pattern}" \\',
        f'    --mcrecon-pattern "{mcrecon_pattern}" \\',
        f'    --mcgen-pattern "{mcgen_pattern}" \\',
        f'    --weight-pattern "{weight_pattern}" \\',
        '    --macro-path "$MACRO_PATH" "$@"',
        "",
        'echo "=========================================================="',
        'echo " Refit Pipeline Complete for: $RUN_DIR"',
        'echo "=========================================================="',
        ""
    ])

    try:
        with open(output_script_path, "w", encoding="utf-8") as f:
            f.write("\n".join(script_lines))
        os.chmod(output_script_path, 0o755)
    except OSError as e:
        alt_path = os.path.join("/tmp", os.path.basename(output_script_path))
        with open(alt_path, "w", encoding="utf-8") as f:
            f.write("\n".join(script_lines))
        os.chmod(alt_path, 0o755)
        print(f"Notice: Target script path was read-only ({output_script_path}: {e}); retry script written to fallback: {alt_path}")
        output_script_path = alt_path

    return output_script_path


def main():
    parser = argparse.ArgumentParser(
        description="MomentFitCollector: diagnose BruFit MomentFit runs and generate targeted retries for unsuccessful bins."
    )
    parser.add_argument(
        "--run-dir", "-d", type=str, default=None,
        help="Path to MomentFit run directory (e.g. /d/grid17/septian/EtaPi0Analysis/run/MomentFit_m080200_...)."
    )
    parser.add_argument(
        "--tbin", "-t", type=str, default="all",
        help="t-bin to process: 'all', an index (0..4), or a name like 't010020' (default: 'all')."
    )
    parser.add_argument(
        "--ess-threshold", type=float, default=20.0,
        help="Effective Sample Size threshold for a successful fit (default: 20.0)."
    )
    parser.add_argument(
        "--ess-metric", type=str, default="mean", choices=["mean", "min"],
        help="ESS metric to evaluate against threshold: 'mean' or 'min' (default: 'mean')."
    )
    parser.add_argument(
        "--min-ess-threshold", type=float, default=None,
        help="Safety floor on worst-moment ESS even when --ess-metric=mean (<= 0 disables; default: min(10.0, ess-threshold))."
    )
    parser.add_argument(
        "--burn-in", "-b", type=int, default=2000,
        help="Number of initial MCMC samples to discard as burn-in (default: 2000)."
    )
    parser.add_argument(
        "--max-l", type=int, default=None,
        help="Max L multipole (default: auto-detected or 4)."
    )
    parser.add_argument(
        "--init-scale", type=float, default=None,
        help="Current initial value scale (default: auto-detected or 1.0)."
    )
    parser.add_argument(
        "--step", type=float, default=None,
        help="Current MCMC step size (default: auto-detected or 1.0)."
    )
    parser.add_argument(
        "--iters", "--mcmc-iterations", type=int, default=None,
        help="Current MCMC iterations (default: auto-detected or 10000)."
    )
    parser.add_argument(
        "--m-min", type=float, default=None,
        help="Lower mass boundary in GeV (default: auto-detected or 0.8)."
    )
    parser.add_argument(
        "--m-max", type=float, default=None,
        help="Upper mass boundary in GeV (default: auto-detected or 2.0)."
    )
    parser.add_argument(
        "--n-bins", type=int, default=None,
        help="Total mass bins (default: auto-detected or 30)."
    )
    parser.add_argument(
        "--macro-path", type=str, default=None,
        help="Path to MomentFit.C macro (default: auto-detected)."
    )
    parser.add_argument(
        "--output-script", "-o", type=str, default=None,
        help="Output path for the generated retry shell script (default: <run-dir>/RunMomentFit_retry.sh)."
    )
    parser.add_argument(
        "--summary-json", type=str, default=None,
        help="Output path for JSON summary (default: <run-dir>/unsuccessful_bins_summary.json)."
    )
    parser.add_argument(
        "--max-parallel", "-j", type=int, default=None,
        help=f"Maximum concurrent single-bin refit jobs in generated script (default: half of available CPUs, currently {DEFAULT_MAX_PARALLEL})."
    )
    parser.add_argument(
        "--non-interactive", "-y", action="store_true",
        help="Accept recommended parameters without prompting."
    )
    parser.add_argument(
        "--auto-retry", action="store_true",
        help="Automatically execute the generated retry script after creation."
    )
    parser.add_argument(
        "--no-exec", action="store_true",
        help="Do not execute the retry script from within Python (script generation only)."
    )
    parser.add_argument(
        "--eval-job", nargs=2, metavar=("TBIN", "BIN_IDX"), default=None,
        help="Evaluate outcome for a single completed refit job (t-bin and mass bin index) and report status/diagnosis."
    )
    parser.add_argument(
        "--cov-dir", type=str, default="",
        help="Linear algebra covariance directory forwarded to retry MomentFit.C calls."
    )
    parser.add_argument(
        "--data-dir", type=str, default="",
        help="Input data directory forwarded to retry MomentFit.C calls."
    )
    parser.add_argument(
        "--mc-variance", dest="apply_mc_variance", action="store_true", default=False,
        help="Enable MC variance weighting in retry MomentFit.C calls."
    )
    parser.add_argument(
        "--no-mc-variance", dest="apply_mc_variance", action="store_false",
        help="Disable MC variance weighting in retry MomentFit.C calls."
    )
    parser.add_argument(
        "--negative-pdf", dest="allow_negative_pdf", action="store_true", default=True,
        help="Allow negative PDF evaluations in retry MomentFit.C calls (default: True)."
    )
    parser.add_argument(
        "--no-negative-pdf", dest="allow_negative_pdf", action="store_false",
        help="Disallow negative PDF evaluations in retry MomentFit.C calls."
    )
    parser.add_argument(
        "--tuning", dest="use_tuning", action="store_true", default=True,
        help="Enable BruFit MCMC covariance tuning in retry MomentFit.C calls (default: True)."
    )
    parser.add_argument(
        "--no-tuning", dest="use_tuning", action="store_false",
        help="Disable BruFit MCMC covariance tuning in retry MomentFit.C calls."
    )
    parser.add_argument(
        "--gibbs-block-size", type=int, default=0,
        help="Default Gibbs block size for MomentFit.C retries (0 = full-vector proposals unless enabled on low acceptance)."
    )
    parser.add_argument(
        "--cold-warmup-acc", type=float, default=0.15,
        help="Cold warm-up acceptance rate forwarded to retry MomentFit.C calls (default: 0.15)."
    )
    parser.add_argument("--min-acceptance", type=float, default=0.05, help="MCMC min acceptance bound.")
    parser.add_argument("--max-acceptance", type=float, default=0.15, help="MCMC max acceptance bound.")
    parser.add_argument("--target-acceptance", type=float, default=None, help="MCMC target acceptance rate.")
    parser.add_argument("--tune-mode", type=str, default="MappedRhat", help="MCMC tuning mode (MappedRhat or Acceptance).")
    parser.add_argument("--rhat-window", type=int, default=1000, help="MappedRhat diagnostic window.")
    parser.add_argument("--rhat-target", type=float, default=1.05, help="Target R-hat for adaptive tuning.")
    parser.add_argument("--ess-fraction", type=float, default=0.02, help="Minimum ESS fraction during tuning.")
    parser.add_argument("--rhat-max-retries", type=int, default=5, help="Maximum adaptive R-hat retries.")
    parser.add_argument("--bailout-threshold", type=int, default=10000, help="Bailout threshold proposals.")
    parser.add_argument("--bailout-min-acceptance", type=float, default=0.0001, help="Bailout minimum acceptance.")
    parser.add_argument("--mc-integration-bins", type=int, default=50, help="Angular MC integration grid bins.")
    parser.add_argument("--fixed-covariance", dest="fixed_covariance", action="store_true", default=True)
    parser.add_argument("--no-fixed-covariance", dest="fixed_covariance", action="store_false")
    parser.add_argument("--use-covariance", dest="use_covariance", action="store_true", default=True)
    parser.add_argument("--no-use-covariance", dest="use_covariance", action="store_false")
    parser.add_argument("--set-init-value", dest="set_init_value", action="store_true", default=True)
    parser.add_argument("--no-set-init-value", dest="set_init_value", action="store_false")
    parser.add_argument("--gj-frame", dest="use_gj_frame", action="store_true", default=True)
    parser.add_argument("--hel-frame", dest="use_gj_frame", action="store_false")
    parser.add_argument("--num-cpu", type=int, default=1, help="RooFit NumCPU per job.")
    parser.add_argument("--seed", type=int, default=111, help="gRandom RNG seed.")
    parser.add_argument("--max-retry-iters", type=int, default=MAX_RETRY_ITERS, help="Maximum retry iteration cap.")
    parser.add_argument("--burnin-fraction", type=float, default=0.10, help="Retry burn-in fraction of iterations.")
    parser.add_argument("--low-acceptance-threshold", type=float, default=0.03, help="Low-acceptance threshold for full-vector retries.")
    parser.add_argument("--gibbs-low-acceptance-threshold", type=float, default=0.12, help="Low-acceptance threshold for blocked Gibbs retries.")
    parser.add_argument("--high-acceptance-threshold", type=float, default=0.35, help="High-acceptance threshold for step increase.")
    parser.add_argument("--gibbs-activation-acceptance", type=float, default=0.05, help="Acceptance below which maxL>=6 enables Gibbs.")
    parser.add_argument("--recondition-covariance", dest="recondition_covariance", action="store_true", default=False,
                        help="Opt-in: empirical posterior covariance reconditioning for LOW_ESS retries (default: False).")
    parser.add_argument("--no-recondition-covariance", dest="recondition_covariance", action="store_false")
    parser.add_argument("--recondition-shrinkage", type=float, default=0.2, help="Shrinkage blend toward Stage 2 covariance.")
    parser.add_argument("--recondition-min-eigen-floor", type=float, default=1e-4, help="Correlation eigenvalue floor for reconditioning.")
    parser.add_argument("--data-pattern", type=str, default="", help="Optional custom data file pattern.")
    parser.add_argument("--mcrecon-pattern", type=str, default="", help="Optional custom recon MC file pattern.")
    parser.add_argument("--mcgen-pattern", type=str, default="", help="Optional custom gen MC file pattern.")
    parser.add_argument("--weight-pattern", type=str, default="", help="Optional custom weight file pattern.")
    args = parser.parse_args()

    # Resolve run directory
    run_dir = args.run_dir
    if run_dir is None:
        candidates = glob.glob(os.path.join(repo_root, "run/MomentFit_*")) + \
                     glob.glob("/d/grid17/septian/EtaPi0Analysis/run/MomentFit_*") + \
                     glob.glob("/d/home/septian/EtaPi0Analysis/run/MomentFit_*")
        candidates = [c for c in candidates if os.path.isdir(c)]
        if candidates:
            candidates.sort(key=lambda p: os.path.getmtime(p), reverse=True)
            run_dir = candidates[0]
            print(f"Auto-selected newest run directory: {run_dir}")
        else:
            print("ERROR: No run directory specified and no MomentFit_* run folders found.")
            sys.exit(1)

    run_dir = os.path.abspath(run_dir)
    dir_params = parse_run_dir_params(run_dir)

    max_l = args.max_l if args.max_l is not None else dir_params["max_l"]
    init_scale = args.init_scale if args.init_scale is not None else dir_params["init_scale"]
    step = args.step if args.step is not None else dir_params["step"]
    m_min = args.m_min if args.m_min is not None else dir_params["m_min"]
    m_max = args.m_max if args.m_max is not None else dir_params["m_max"]
    n_bins = args.n_bins if args.n_bins is not None else dir_params["n_bins"]
    current_iters = args.iters if args.iters is not None else dir_params["mcmc_iterations"]

    # Resolve macro path
    macro_path = args.macro_path
    if macro_path is None:
        for cand in [
            os.path.join(repo_root, "macros/MomentFit.C"),
            "/d/home/septian/Moments2Amplitudes/macros/MomentFit.C",
            "/d/home/septian/EtaPi0Analysis/macros/MomentFit.C",
            "/d/grid17/septian/EtaPi0Analysis/macros/MomentFit.C",
            "/d/home/septian/EtaPi0Analysis/macros/MomentFit_General.C",
            "/d/grid17/septian/EtaPi0Analysis/macros/MomentFit_General.C"
        ]:
            if os.path.isfile(cand):
                macro_path = cand
                break
    if not macro_path:
        macro_path = "/d/home/septian/EtaPi0Analysis/macros/MomentFit.C"

    if args.min_ess_threshold is None:
        min_ess_threshold = min(10.0, args.ess_threshold) if args.ess_metric == "mean" else None
    elif args.min_ess_threshold <= 0.0:
        min_ess_threshold = None
    else:
        min_ess_threshold = args.min_ess_threshold

    # Fast evaluation of a single completed refit job
    if args.eval_job:
        tbin_arg, bin_idx_arg = args.eval_job
        try:
            b_idx = int(bin_idx_arg)
        except ValueError:
            print(f"Error: Invalid bin index '{bin_idx_arg}'. Must be an integer.")
            sys.exit(4)
        code = evaluate_single_job(
            run_dir=run_dir,
            tbin_str=tbin_arg,
            bin_idx=b_idx,
            m_min=m_min,
            m_max=m_max,
            n_bins=n_bins,
            ess_threshold=args.ess_threshold,
            ess_metric=args.ess_metric,
            min_ess_threshold=min_ess_threshold,
            burn_in=args.burn_in,
            current_scale=init_scale,
            current_step=step,
            current_iters=current_iters,
            current_gibbs=args.gibbs_block_size,
            max_l=max_l,
            max_retry_iters=args.max_retry_iters,
            burnin_fraction=args.burnin_fraction,
            low_acceptance_threshold=args.low_acceptance_threshold,
            gibbs_low_acceptance_threshold=args.gibbs_low_acceptance_threshold,
            high_acceptance_threshold=args.high_acceptance_threshold,
            gibbs_activation_acceptance=args.gibbs_activation_acceptance,
        )
        sys.exit(code)

    # Discover t-bins
    tbin_dirs = discover_tbin_dirs(run_dir, args.tbin)
    if not tbin_dirs:
        print(f"ERROR: No matching t-bin directories found in '{run_dir}'.")
        sys.exit(1)

    print(f"\nEvaluating {len(tbin_dirs)} t-bin directory/directories:")
    for idx, name, p in tbin_dirs:
        print(f"  - [{idx}] {name}: {p}")

    # Fast Phase A: scan each t-bin directory once and reuse cached diagnostics
    scans_by_tbin = {}
    all_pending_tasks = []
    for t_idx, t_name, t_path in tbin_dirs:
        scan = _scan_tbin_directory(
            tbin_dir=t_path,
            run_dir=run_dir,
            tbin_name=t_name,
            m_min=m_min,
            m_max=m_max,
            n_bins=n_bins,
            burn_in=args.burn_in,
        )
        scans_by_tbin[t_idx] = (t_name, t_path, scan)
        for task in scan["pending_tasks"]:
            all_pending_tasks.append((t_idx, *task))

    # Fast Phase B: process uncached ROOT inspections across all (t-bin, mass-bin) tasks in parallel processes
    root_results_by_tbin: dict[int, dict[int, dict]] = {t_idx: {} for t_idx, _, _ in tbin_dirs}
    for res_t, res_b, res_dict in _run_root_inspection_tasks(all_pending_tasks):
        root_results_by_tbin[res_t][res_b] = res_dict

    all_evaluations = []
    for t_idx, t_name, t_path in tbin_dirs:
        _, _, scan = scans_by_tbin[t_idx]
        evals = evaluate_tbin(
            tbin_idx=t_idx,
            tbin_name=t_name,
            tbin_dir=t_path,
            run_dir=run_dir,
            m_min=m_min,
            m_max=m_max,
            n_bins=n_bins,
            ess_threshold=args.ess_threshold,
            ess_metric=args.ess_metric,
            min_ess_threshold=min_ess_threshold,
            burn_in=args.burn_in,
            current_scale=init_scale,
            current_step=step,
            current_iters=current_iters,
            current_gibbs=args.gibbs_block_size,
            max_l=max_l,
            prepared_scan=scan,
            precomputed_root_results=root_results_by_tbin.get(t_idx),
            max_retry_iters=args.max_retry_iters,
            burnin_fraction=args.burnin_fraction,
            low_acceptance_threshold=args.low_acceptance_threshold,
            gibbs_low_acceptance_threshold=args.gibbs_low_acceptance_threshold,
            high_acceptance_threshold=args.high_acceptance_threshold,
            gibbs_activation_acceptance=args.gibbs_activation_acceptance,
        )
        all_evaluations.extend(evals)

    # Print summary table
    print_evaluation_summary(
        all_evaluations,
        ess_threshold=args.ess_threshold,
        ess_metric=args.ess_metric,
        min_ess_threshold=min_ess_threshold,
    )

    unsuccessful_bins = [b for b in all_evaluations if not b["is_successful"]]

    # Save JSON summary
    summary_path = args.summary_json if args.summary_json else os.path.join(run_dir, "unsuccessful_bins_summary.json")
    json_data = {
        "run_dir": run_dir,
        "timestamp": datetime.now().isoformat(),
        "ess_threshold": args.ess_threshold,
        "ess_metric": args.ess_metric,
        "min_ess_threshold": min_ess_threshold,
        "burn_in": args.burn_in,
        "total_bins": len(all_evaluations),
        "successful_count": len(all_evaluations) - len(unsuccessful_bins),
        "unsuccessful_count": len(unsuccessful_bins),
        "parameters": {
            "max_l": max_l,
            "init_scale": init_scale,
            "step": step,
            "m_min": m_min,
            "m_max": m_max,
            "n_bins": n_bins,
            "gibbs_block_size": args.gibbs_block_size,
            "max_retry_iters": args.max_retry_iters,
            "burnin_fraction": args.burnin_fraction,
            "low_acceptance_threshold": args.low_acceptance_threshold,
            "gibbs_low_acceptance_threshold": args.gibbs_low_acceptance_threshold,
            "high_acceptance_threshold": args.high_acceptance_threshold,
            "gibbs_activation_acceptance": args.gibbs_activation_acceptance,
            "recondition_covariance": args.recondition_covariance,
        },
        "unsuccessful_bins": [
            {
                "tbin_idx": b["tbin_idx"],
                "tbin_name": b["tbin_name"],
                "bin_idx": b["bin_idx"],
                "mass_center": b["mass_center"],
                "failure_type": b["failure_type"],
                "status": b["status"],
                "mean_ess": b.get("mean_ess"),
                "min_ess": b.get("min_ess"),
                "eval_ess": b.get("eval_ess"),
                "worst_moment": b["worst_moment"],
                "diagnosis_reason": b["diagnosis_reason"],
                "current_scale": b["current_scale"],
                "current_step": b["current_step"],
                "current_gibbs": b.get("current_gibbs", 0),
                "rec_scale": b["rec_scale"],
                "rec_step": b["rec_step"],
                "rec_iters": b["rec_iters"],
                "rec_burnin": b.get("rec_burnin", args.burn_in),
                "rec_gibbs": b.get("rec_gibbs", 0),
                "rec_cov_dir": b.get("rec_cov_dir"),
            }
            for b in unsuccessful_bins
        ]
    }
    try:
        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(json_data, f, indent=2)
        print(f"Diagnostics summary JSON saved to: {summary_path}")
    except OSError as e:
        alt_summary = os.path.join("/tmp", os.path.basename(summary_path))
        try:
            with open(alt_summary, "w", encoding="utf-8") as f:
                json.dump(json_data, f, indent=2)
            print(f"Notice: Target path was read-only ({summary_path}: {e}); summary JSON written to fallback: {alt_summary}")
        except OSError:
            print(f"Notice: Could not write diagnostics summary JSON to {summary_path} ({e}).")

    if not unsuccessful_bins:
        print("\nAll mass bins are successful! No retry script generation needed.")
        return

    n_ph1 = sum(1 for b in unsuccessful_bins if b["failure_type"] == "ZERO_ACC_PH1")
    n_cov = sum(1 for b in unsuccessful_bins if b["failure_type"] == "STUCK_COVARIANCE")
    n_low = sum(1 for b in unsuccessful_bins if b["failure_type"] == "LOW_ESS")
    n_other = sum(1 for b in unsuccessful_bins if b["failure_type"] not in ("ZERO_ACC_PH1", "STUCK_COVARIANCE", "LOW_ESS"))

    print("\n" + "=" * 80)
    print(f" Targeted Retry Summary Breakdown ({len(unsuccessful_bins)} / {len(all_evaluations)} mass bins require refit):")
    if n_ph1 > 0:
        print(f"   ❌ STUCK in 1DStep (Phase 1):      {n_ph1:>2d} bin(s) -> Suggestion: reduce initial scale to fit likelihood basin")
    if n_cov > 0:
        print(f"   ❌ STUCK in Covariance (Phase 2):  {n_cov:>2d} bin(s) -> Suggestion: reduce step size to improve proposal acceptance")
    if n_low > 0:
        print(f"   ⚠️  SUCCESS with Low ESS:           {n_low:>2d} bin(s) -> Suggestion: increase MCMC iterations or adjust step size")
    if n_other > 0:
        print(f"   ❌ FAILED / BAILED OUT:            {n_other:>2d} bin(s) -> Suggestion: safe parameters / check ROOT log")
    print("=" * 80)

    # Check for user confirmation / parameter adjustment if interactive
    if not args.non_interactive and not args.auto_retry and not args.no_exec and sys.stdin.isatty():
        try:
            print("")
            user_input = input("Accept recommended parameters and generate retry script? [Y/n/custom]: ").strip()
            if user_input.lower() == "n":
                print("Aborting retry script generation.")
                return
            elif user_input.lower() == "custom":
                custom_scale_in = input(f"Enter global retry initial scale [default {init_scale}]: ").strip()
                custom_step_in = input(f"Enter global retry step size [default {step}]: ").strip()
                custom_iters_in = input(f"Enter retry MCMC iterations [default {current_iters}]: ").strip()

                for b in unsuccessful_bins:
                    if custom_scale_in:
                        b["rec_scale"] = float(custom_scale_in)
                    if custom_step_in:
                        b["rec_step"] = float(custom_step_in)
                    if custom_iters_in:
                        b["rec_iters"] = int(custom_iters_in)
                print("Updated retry parameters with custom user inputs.")
        except (EOFError, KeyboardInterrupt):
            pass

    # Optional empirical posterior covariance reconditioning for converged LOW_ESS bins
    if args.recondition_covariance and args.fixed_covariance:
        n_recond = 0
        for b in unsuccessful_bins:
            if b["failure_type"] == "LOW_ESS":
                recond_dir = recondition_covariance_for_bin(
                    bin_record=b,
                    run_dir=run_dir,
                    src_cov_dir=args.cov_dir,
                    max_l=max_l,
                    burn_in=args.burn_in,
                    shrinkage=args.recondition_shrinkage,
                    min_eigen_floor=args.recondition_min_eigen_floor,
                )
                if recond_dir:
                    b["rec_cov_dir"] = recond_dir
                    n_recond += 1
        if n_recond > 0:
            print(f"Reconditioned empirical SPD covariance matrices for {n_recond} LOW_ESS bin(s).")

    # Generate retry shell script
    retry_script = generate_retry_script(
        run_dir=run_dir,
        unsuccessful_bins=unsuccessful_bins,
        macro_path=macro_path,
        max_l=max_l,
        m_min=m_min,
        m_max=m_max,
        n_bins=n_bins,
        burn_in=args.burn_in,
        ess_threshold=args.ess_threshold,
        min_ess_threshold=min_ess_threshold,
        max_parallel=args.max_parallel,
        output_script_path=args.output_script,
        cov_dir=args.cov_dir,
        data_dir=args.data_dir,
        apply_mc_variance=args.apply_mc_variance,
        allow_negative_pdf=args.allow_negative_pdf,
        use_tuning=args.use_tuning,
        cold_warmup_acc=args.cold_warmup_acc,
        data_pattern=args.data_pattern,
        mcrecon_pattern=args.mcrecon_pattern,
        mcgen_pattern=args.mcgen_pattern,
        weight_pattern=args.weight_pattern,
        ess_metric=args.ess_metric,
        min_acceptance=args.min_acceptance,
        max_acceptance=args.max_acceptance,
        target_acceptance=args.target_acceptance,
        tune_mode=args.tune_mode,
        rhat_window=args.rhat_window,
        rhat_target=args.rhat_target,
        ess_fraction=args.ess_fraction,
        rhat_max_retries=args.rhat_max_retries,
        bailout_threshold=args.bailout_threshold,
        bailout_min_acceptance=args.bailout_min_acceptance,
        mc_integration_bins=args.mc_integration_bins,
        fixed_covariance=args.fixed_covariance,
        use_covariance=args.use_covariance,
        set_init_value=args.set_init_value,
        use_gj_frame=args.use_gj_frame,
        num_cpu=args.num_cpu,
        seed=args.seed,
        max_retry_iters=args.max_retry_iters,
        burnin_fraction=args.burnin_fraction,
        low_acceptance_threshold=args.low_acceptance_threshold,
        gibbs_low_acceptance_threshold=args.gibbs_low_acceptance_threshold,
        high_acceptance_threshold=args.high_acceptance_threshold,
        gibbs_activation_acceptance=args.gibbs_activation_acceptance,
        recondition_covariance=args.recondition_covariance,
        recondition_shrinkage=args.recondition_shrinkage,
        recondition_min_eigen_floor=args.recondition_min_eigen_floor,
    )
    print(f"\nGenerated ready-to-run retry script: {retry_script}")
    print(f"To manually inspect and launch the retries:")
    print(f"  bash {retry_script}\n")

    # If --auto-retry was specified, execute the retry script immediately
    if args.auto_retry and not args.no_exec:
        print("==========================================================")
        print(" --auto-retry specified: Automatically executing retry script...")
        print("==========================================================")
        subprocess.run(["bash", retry_script], check=True)
    elif not args.non_interactive and not args.no_exec and sys.stdin.isatty():
        try:
            print("-" * 80)
            prompt_exec = input(f"Would you like to execute the retry script now for these {len(unsuccessful_bins)} unsuccessful bin(s)? [Y/n]: ").strip()
            if prompt_exec.lower() not in ("n", "no"):
                print("==========================================================")
                print(f" Executing targeted retry script for {len(unsuccessful_bins)} bin(s)...")
                print("==========================================================")
                subprocess.run(["bash", retry_script], check=True)
            else:
                print("Skipping retry execution. You can manually launch it later:")
                print(f"  bash {retry_script}\n")
        except (EOFError, KeyboardInterrupt):
            pass


if __name__ == "__main__":
    main()

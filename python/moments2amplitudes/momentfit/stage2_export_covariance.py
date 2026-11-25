"""
Stage 2: Export moments values, format covariance matrices, and regularize to SPD for BruFit.
"""

from __future__ import annotations

import functools
import glob
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .config import PipelineConfig

print = functools.partial(print, flush=True)


def _get_py39_env(mom2amp_dir: str) -> dict:
    """Constructs a clean isolated environment for Python 3.9 + ROOT 6.28."""
    repo_root = Path(mom2amp_dir)
    py39_env_dir = repo_root / ".pixi" / "envs" / "py39"
    if not (py39_env_dir / "bin" / "python").is_file():
        py39_env_dir = Path("/d/grid17/septian/Moments-dev/Moments.Nizar/.pixi/envs/default")

    clean_pypath = [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p and not any(k in p for k in ("root_v6", "cuda_install", "cuda_build", "brufit"))]
    clean_ldpath = [p for p in os.environ.get("LD_LIBRARY_PATH", "").split(os.pathsep) if p and not any(k in p for k in ("root_v6", "cuda_install", "cuda_build", "brufit"))]

    env = os.environ.copy()
    for k in list(env.keys()):
        if k.startswith("ROOT"):
            env.pop(k, None)
    env["PATH"] = f"{py39_env_dir}/bin:/usr/bin:/bin"
    env["LD_LIBRARY_PATH"] = f"{py39_env_dir}/lib"
    env["PYTHONPATH"] = f"/d/grid17/septian/Moments-dev/Moments.Main:{repo_root}/python"
    env["PYTHONNOUSERSITE"] = "1"
    env["MPLCONFIGDIR"] = "/tmp"
    return env


def _find_py39_runner(mom2amp_dir: str) -> Tuple[List[str], dict]:
    """
    Locates the Python 3.9 / ROOT 6.28 runner command prefix for linear algebra moment unpickling.
    """
    repo_root = Path(mom2amp_dir)
    env = _get_py39_env(mom2amp_dir)

    # 1. Direct py39 environment binary if present
    py39_bin = repo_root / ".pixi" / "envs" / "py39" / "bin" / "python"
    if py39_bin.is_file():
        return [str(py39_bin)], env

    # 2. Moments.Nizar fallback
    nizar_bin = Path("/d/grid17/septian/Moments-dev/Moments.Nizar/.pixi/envs/default/bin/python")
    if nizar_bin.is_file():
        return [str(nizar_bin)], env

    # 3. Pixi run
    pixi_bin = shutil.which("pixi") or os.path.expanduser("~/.pixi/bin/pixi")
    manifest_path = repo_root / "pixi.toml"
    if os.path.isfile(pixi_bin) and manifest_path.is_file():
        return [pixi_bin, "run", "--manifest-path", str(manifest_path), "-e", "py39", "python"], env

    raise RuntimeError(
        "Could not locate a valid Python 3.9 + ROOT 6.28 environment for Moments.Main pickle extraction.\n"
        "Please ensure '.pixi/envs/py39' or '/d/grid17/septian/Moments-dev/Moments.Nizar' is available."
    )


def _init_export_environment(cfg: PipelineConfig) -> None:
    """Initialize environment for Python 3.9 worker."""
    moments_main_dir = cfg.environment.moments_main_dir
    if moments_main_dir not in sys.path:
        sys.path.insert(0, moments_main_dir)

    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetMacroPath(f"{moments_main_dir}:{ROOT.gROOT.GetMacroPath()}")

    # Ensure pickle compatibility
    import moments.MomentCalculator
    sys.modules["MomentCalculator"] = moments.MomentCalculator


def extract_covariance_matrices_py39(
    cfg: PipelineConfig,
    tbin: str,
    max_l: int,
) -> bool:
    """
    Extracts raw covariance matrices and moment values from pickle using Moments.Main (Python 3.9).
    """
    _init_export_environment(cfg)

    import ROOT
    from moments.MomentCalculator import (
        MomentCalculator,
        MomentResultsKinematicBinning,
        MomentValue,
        QnMomentIndex,
        binLabel,
    )

    out_dir = cfg.linear_algebra.get_output_dir(tbin, max_l, cfg.analysis.beam_pol_label)
    prefix = "norm" if cfg.linear_algebra.normalize_moments else "unnorm"
    pkl_path = os.path.join(out_dir, f"{prefix}_moments_phys.pkl")

    if not os.path.isfile(pkl_path):
        print(f"❌ Error: Moments pickle file not found for tbin={tbin}, maxL={max_l}: {pkl_path}")
        return False

    print(f"\n[Extracting Moments] tbin={tbin}, maxL={max_l}")
    print(f"  Input: {pkl_path}")

    moment_results: MomentResultsKinematicBinning = MomentResultsKinematicBinning.loadPickle(pkl_path)
    n_bins = len(moment_results)
    print(f"  Loaded {n_bins} moment result bins.")

    if n_bins == 0:
        print(f"❌ Error: No moment results found in {pkl_path}")
        return False

    first_result = moment_results[0]
    param_names: List[str] = []

    counter = 0
    for moment_value in first_result:
        qn: QnMomentIndex = moment_value.qn
        idx = qn.momentIndex
        L = qn.L
        M = qn.M
        if idx == 0 and L == 0 and M == 0:
            continue
        param_names.append(f"H_{idx}_{L}_{M}")
        counter += 1
        if counter == (len(first_result) - 1):
            break

    param_names.append("Yld_Moments")
    n_moments = len(param_names)

    bin_width = cfg.analysis.mass_binning.bin_width

    for bin_idx, moment_result in enumerate(moment_results):
        h000_qn = QnMomentIndex(momentIndex=0, L=0, M=0)
        h000_flat = moment_result.indices[h000_qn]
        h000_val: MomentValue = moment_result[h000_flat]
        h000_real, h000_uncert = h000_val.real

        mass_center = None
        for var, center in moment_result.binCenters.items():
            if var.name == cfg.analysis.mass_binning.variable:
                mass_center = center
                break
        if mass_center is None:
            mass_center = cfg.analysis.mass_binning.get_bin_center(bin_idx)

        bin_low = mass_center - (bin_width / 2.0)
        bin_high = mass_center + (bin_width / 2.0)

        moment_result.normalize()

        cov_matrix = ROOT.TMatrixDSym(n_moments)
        flat_indices: List[int] = []
        use_real_parts: List[bool] = []
        h1_norm_moments: Dict[str, ROOT.TH1D] = {}

        for label in param_names:
            if label.startswith("H_"):
                parts = label.split("_")
                m_idx = int(parts[1])
                qn = QnMomentIndex(
                    momentIndex=m_idx,
                    L=int(parts[2]),
                    M=int(parts[3]),
                )
                flat_idx = moment_result.indices[qn]
                flat_indices.append(flat_idx)
                use_real_part = m_idx != 2
                use_real_parts.append(use_real_part)

                val = moment_result[flat_idx]
                val_num = val.real if use_real_part else val.imag
                if hasattr(val_num, "nominal_value"):
                    mean_val = float(val_num.nominal_value)
                    err_val = float(val_num.std_dev)
                elif isinstance(val_num, (tuple, list)):
                    mean_val = float(val_num[0])
                    err_val = float(val_num[1])
                else:
                    mean_val = float(val_num)
                    err_val = 0.0

                h1 = ROOT.TH1D(label, label, 1, bin_low, bin_high)
                h1.SetBinContent(1, mean_val * 2.0)
                h1.SetBinError(1, err_val * 2.0)
                h1_norm_moments[label] = h1
            else:
                flat_indices.append(h000_flat)
                use_real_parts.append(True)

        scales = [float(h000_real) if name == "Yld_Moments" else 2.0 for name in param_names]

        for i in range(n_moments):
            for j in range(i, n_moments):
                cov_est = moment_result.covariance((flat_indices[i], flat_indices[j]), (use_real_parts[i], use_real_parts[j]))
                orig = float(cov_est[0, 0] if i == j else cov_est[0, 1])
                scaled = orig * scales[i] * scales[j]
                cov_matrix[i][j] = scaled
                cov_matrix[j][i] = scaled

        bin_lbl = binLabel(moment_result)
        raw_cov_path = os.path.join(out_dir, f"covariance_matrix_{bin_lbl}.root")
        root_file = ROOT.TFile(raw_cov_path, "RECREATE")
        cov_matrix.Write("MyCovMatrix")

        for label in param_names:
            if label in h1_norm_moments:
                h1_norm_moments[label].Write()

        param_names_str = ROOT.TObjString(",".join(param_names))
        param_names_str.Write("paramNames")
        root_file.Close()

    print(f"  Extracted {n_bins} raw covariance matrices to {out_dir}")
    return True


def export_single_task(
    cfg: PipelineConfig,
    tbin: str,
    max_l: int,
    plot_cov: bool = False,
    plot_cov_bins: Optional[List[str]] = None,
) -> bool:
    """
    Exports covariance and moment values for a single (tbin, maxL) task:
      1. Extract raw covariance matrix via Python 3.9 (Moments.Main)
      2. Regularize to SPD and save in BruFit format via ROOT 6.36
      3. Optionally plot proposal covariance comparisons before and after alignment
    """
    out_dir = cfg.linear_algebra.get_output_dir(tbin, max_l, cfg.analysis.beam_pol_label)

    # 1. Extraction step
    is_py312 = sys.version_info >= (3, 10)
    if is_py312:
        cmd_prefix, py39_env = _find_py39_runner(cfg.environment.mom2amp_dir)
        cmd = cmd_prefix + [
            "-m", "moments2amplitudes.momentfit.stage2_export_covariance",
            "--extract-only",
            "--config", cfg.config_path,
            "--tbin", tbin,
            "--max-l", str(max_l),
        ]
        ret = subprocess.run(cmd, env=py39_env)
        if ret.returncode != 0:
            print(f"❌ Error: Python 3.9 extraction failed for tbin={tbin}, maxL={max_l} (code {ret.returncode})")
            return False
    else:
        success = extract_covariance_matrices_py39(cfg, tbin, max_l)
        if not success:
            return False

    # 2. SPD regularization in BruFit format using ROOT 6.36
    sh_script = os.path.join(cfg.environment.mom2amp_dir, "macros", "ExportCovariance.sh")
    if not os.path.isfile(sh_script):
        sh_script = os.path.join(cfg.environment.mom2amp_dir, "macros", "legacy", "ExportCovariance.sh")

    if os.path.isfile(sh_script):
        print(f"  Exporting BruFit-format SPD covariance matrices using ROOT 6.36 via {sh_script}...")
        ret = subprocess.run(["bash", sh_script, out_dir])
        if ret.returncode != 0:
            print(f"❌ Error: ExportCovariance.sh failed with exit code {ret.returncode}")
            return False
    else:
        import ROOT
        macro_path = cfg.export_covariance.macro_path
        ROOT.gROOT.LoadMacro(f"{macro_path}+")
        ROOT.ExportCovarianceC(out_dir)

    # Verify that BruCovMatrix files exist
    bru_files = glob.glob(os.path.join(out_dir, "BruCovMatrix_*.root"))
    if not bru_files:
        print(f"❌ Error: No BruCovMatrix files generated in {out_dir}")
        return False

    print(f"  ✅ Successfully exported {len(bru_files)} BruCovMatrix bins for tbin={tbin}, maxL={max_l}.")

    # 3. Optional proposal covariance comparison plotting
    do_plot = plot_cov or getattr(cfg.export_covariance, "plot_covariance", False)
    if do_plot:
        try:
            from .plot_covariance import plot_covariance_for_directory
            plot_covariance_for_directory(
                cov_dir=out_dir,
                mass_bins_filter=plot_cov_bins,
                tbin=tbin,
                max_l=max_l,
            )
        except Exception as e:
            print(f"  ⚠️ Warning: Covariance plotting failed for tbin={tbin}, maxL={max_l}: {e}")

    return True


def run_stage2(
    cfg: PipelineConfig,
    tbins_filter: Optional[List[str]] = None,
    max_l_filter: Optional[List[int]] = None,
    dry_run: bool = False,
    plot_cov: bool = False,
    plot_cov_bins: Optional[List[str]] = None,
) -> bool:
    """
    Executes Stage 2 across all configured t-bins and max_l multipoles.
    """
    tbins = tbins_filter if tbins_filter is not None else cfg.analysis.tbins
    max_ls = max_l_filter if max_l_filter is not None else cfg.analysis.max_l

    tasks: List[Tuple[str, int]] = [(tbin, l) for tbin in tbins for l in max_ls]

    print("\n" + "=" * 70)
    print(" STAGE 2: Export Moments Value and Covariance (SPD Regularization)")
    print("=" * 70)
    print(f" Target t-bins:      {', '.join(tbins)}")
    print(f" Multipoles Max L:   {max_ls}")
    print(f" Total Tasks:        {len(tasks)}")
    print(f" Macro:              {cfg.export_covariance.macro_path}")
    if plot_cov or getattr(cfg.export_covariance, "plot_covariance", False):
        print(f" Plot Covariance:    Enabled")
    print("=" * 70 + "\n")

    if dry_run:
        print("[DRY-RUN] Tasks to export:")
        for tbin, l in tasks:
            out = cfg.linear_algebra.get_output_dir(tbin, l, cfg.analysis.beam_pol_label)
            print(f"  * tbin={tbin}, maxL={l} -> {out}/BruCovMatrix_*.root")
        return True

    failed_tasks: List[Tuple[str, int]] = []
    for tbin, l in tasks:
        success = export_single_task(cfg, tbin, l, plot_cov=plot_cov, plot_cov_bins=plot_cov_bins)
        if not success:
            failed_tasks.append((tbin, l))

    print("\n" + "-" * 70)
    if failed_tasks:
        print(f"❌ Stage 2 encountered errors in {len(failed_tasks)} task(s):")
        for tbin, l in failed_tasks:
            print(f"  - tbin={tbin}, maxL={l}")
        return False
    else:
        print(f"✅ Stage 2 completed successfully! All {len(tasks)} tasks exported.")
        return True


def main() -> None:
    import argparse
    from .config import load_pipeline_config

    parser = argparse.ArgumentParser(description="Stage 2 Covariance Extraction Worker (Python 3.9 / ROOT 6.28)")
    parser.add_argument("--extract-only", action="store_true", help="Extract raw covariance matrices only")
    parser.add_argument("--config", "-c", required=True, help="Path to YAML configuration")
    parser.add_argument("--tbin", "-t", required=True, help="Target t-bin")
    parser.add_argument("--max-l", "-l", type=int, required=True, help="Target max L")
    args = parser.parse_args()

    cfg = load_pipeline_config(args.config)
    if args.extract_only:
        success = extract_covariance_matrices_py39(cfg, args.tbin, args.max_l)
    else:
        success = export_single_task(cfg, args.tbin, args.max_l)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()

"""Unit tests for BruFitResults status checks, diagnostics cache, and monitor helpers."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "python"))
sys.path.insert(0, str(REPO_ROOT / "macros"))

import BruFitResults as bfr
from moments2amplitudes.momentfit.config import load_momentfit_config
from moments2amplitudes.momentfit.runner import run_monitor


class TestMonitorAndBinnerFixes(unittest.TestCase):
    """Tests for MCMCTree status fallback, diagnostics caching, and monitor CLI plumbing."""

    def test_get_fit_status_recovers_stalled_log_with_valid_mcmctree(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bin_dir = Path(tmpdir) / "Mpi0eta1.060000_"
            bin_dir.mkdir(parents=True)
            log_path = Path(tmpdir) / "logRooFit0.txt"
            log_path.write_text(
                "--> Target Bin : Mpi0eta1.060000_\n"
                "Step 15800/15811 acceptance 0.077\n",
                encoding="utf-8",
            )
            old_time = 1_700_000_000.0
            os.utime(log_path, (old_time, old_time))

            root_path = bin_dir / "ResultsBruMcmcFixedCovariance.root"
            root_path.write_bytes(b"root_fake_valid_header")

            with mock.patch.object(bfr, "_check_mcmc_root_valid", return_value=(True, 15812)):
                res = bfr.get_fit_status(tmpdir, active_timeout_min=30.0, compute_diagnostics=False)

            self.assertEqual(res["summary"], {"SUCCESS": 1})
            self.assertEqual(res["success_bins"], [0])
            self.assertEqual(res["stalled_bins"], [])
            self.assertIn("MCMCTree verified (15812 entries)", res["bins"][0]["details"])

    def test_get_fit_status_does_not_override_explicit_failed_log(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bin_dir = Path(tmpdir) / "Mpi0eta1.060000_"
            bin_dir.mkdir(parents=True)
            log_path = Path(tmpdir) / "logRooFit0.txt"
            log_path.write_text(
                "--> Target Bin : Mpi0eta1.060000_\n"
                "WARNING FitManager::Run: Fit failed\n",
                encoding="utf-8",
            )
            root_path = bin_dir / "ResultsBruMcmcFixedCovariance.root"
            root_path.write_bytes(b"root_fake_valid_header")

            with mock.patch.object(bfr, "_check_mcmc_root_valid", return_value=(True, 100)):
                res = bfr.get_fit_status(tmpdir, active_timeout_min=30.0, compute_diagnostics=False)

            self.assertEqual(res["summary"], {"FAILED": 1})
            self.assertEqual(res["failed_bins"], [0])
            self.assertEqual(res["success_bins"], [])

    def test_moments_mcmc_get_mcmc_diagnostics_uses_batch_fft_and_cache(self):
        rng = np.random.default_rng(7)
        n_rows = 2500
        all_data = {
            "H_0_0_0": rng.normal(1.0, 0.1, size=n_rows),
            "H_0_2_0": rng.normal(0.0, 0.2, size=n_rows),
            "Yld_Moments": rng.normal(1000.0, 20.0, size=n_rows),
            "weight": np.ones(n_rows, dtype=np.float64),
        }

        mcmc = object.__new__(bfr.MomentsMCMC)
        mcmc.burn_in = 500
        mcmc.fit_bin = [1.06]
        mcmc._moment_names = ["H_0_0_0", "H_0_2_0"]
        mcmc._bin_data_cache = {}
        mcmc._diagnostics_cache = {}

        with mock.patch.object(mcmc, "_get_bin_numpy_data", return_value=all_data) as mock_get:
            d1 = mcmc.get_mcmc_diagnostics(0)
            d2 = mcmc.get_mcmc_diagnostics(0)

        mock_get.assert_called_once_with(0)
        self.assertEqual(d1, d2)
        self.assertEqual(d1["n_raw"], n_rows)
        self.assertEqual(d1["burn_in"], 500)
        self.assertEqual(d1["n_samples"], 2000)
        self.assertIn("H_0_2_0", d1["per_moment"])
        self.assertGreater(d1["mean_ess"], 0.0)

    def test_infer_mass_binning_and_runner_monitor_dry_run(self):
        import MonitorMomentFit as mmf

        with tempfile.TemporaryDirectory() as tmpdir:
            tbin_dir = Path(tmpdir) / "MomentFit_m104160_FixedCov_Lmax4" / "t010020"
            for i in range(14):
                center = 1.04 + (i + 0.5) * (1.60 - 1.04) / 14.0
                (tbin_dir / f"Mpi0eta{center:.6f}_").mkdir(parents=True)

            m_min, m_max, n_bins = mmf._infer_mass_binning(str(tbin_dir))
            self.assertAlmostEqual(m_min, 1.04)
            self.assertAlmostEqual(m_max, 1.60)
            self.assertEqual(n_bins, 14)

            cfg = load_momentfit_config(str(REPO_ROOT / "configs" / "momentfit_test.yaml"))
            with mock.patch("builtins.print") as mock_print:
                ok = run_monitor(cfg=cfg, tbins_filter=["t010020"], max_l_filter=[4], dry_run=True)

            self.assertTrue(ok)
            printed = "\n".join(" ".join(str(a) for a in call.args) for call in mock_print.call_args_list)
            self.assertIn("--max-l 4", printed)
            self.assertIn("--m-min 1.04", printed)
            self.assertIn("--m-max 1.6", printed)
            self.assertIn("--n-bins 14", printed)

    def test_low_ess_retry_iterations_increase_up_to_200000(self):
        import MomentFitCollector as mfc

        self.assertEqual(mfc.MAX_RETRY_ITERS, 200000)

        # Healthy acceptance with large tau scales above old 40000 cap up to 200000
        _, step_h, iters_h, _ = mfc.recommend_low_ess_params(
            bin_scale=1.0, bin_step=1.0, bin_iters=40000, burn_in=2000,
            eval_ess=8.0, ess_threshold=20.0, mean_tau=5000.0, bin_acceptance=0.20,
        )
        self.assertAlmostEqual(step_h, 1.0)
        self.assertEqual(iters_h, 127000)

        # Further retry at 150000 hits 200000 cap cleanly
        _, _, iters_cap, _ = mfc.recommend_low_ess_params(
            bin_scale=1.0, bin_step=1.0, bin_iters=150000, burn_in=2000,
            eval_ess=15.0, ess_threshold=20.0, mean_tau=200.0, bin_acceptance=0.20,
        )
        self.assertEqual(iters_cap, 200000)

        # Low acceptance (< 0.03 full-vector, or < 0.12 blocked Gibbs) decreases step and scales iters up to 200000
        _, step_low, iters_low, _ = mfc.recommend_low_ess_params(
            bin_scale=1.0, bin_step=1.0, bin_iters=40000, burn_in=2000,
            eval_ess=10.0, ess_threshold=20.0, mean_tau=300.0, bin_acceptance=0.02,
            bin_gibbs=0,
        )
        self.assertAlmostEqual(step_low, 0.6)
        self.assertEqual(iters_low, 50000)

        _, step_low_g, iters_low_g, _ = mfc.recommend_low_ess_params(
            bin_scale=1.0, bin_step=1.0, bin_iters=40000, burn_in=2000,
            eval_ess=10.0, ess_threshold=20.0, mean_tau=300.0, bin_acceptance=0.05,
            bin_gibbs=12,
        )
        self.assertAlmostEqual(step_low_g, 0.6)
        self.assertEqual(iters_low_g, 50000)

        # High acceptance increases step and scales iters above 40000 up to 200000
        _, step_high, iters_high, _ = mfc.recommend_low_ess_params(
            bin_scale=1.0, bin_step=1.0, bin_iters=80000, burn_in=2000,
            eval_ess=10.0, ess_threshold=20.0, mean_tau=300.0, bin_acceptance=0.50,
        )
        self.assertAlmostEqual(step_high, 1.5)
        self.assertEqual(iters_high, 100000)

        # At 200000, iters stays at 200000 without overflowing
        _, _, iters_at_max, _ = mfc.recommend_low_ess_params(
            bin_scale=1.0, bin_step=1.0, bin_iters=200000, burn_in=2000,
            eval_ess=10.0, ess_threshold=20.0, mean_tau=5000.0, bin_acceptance=0.20,
        )
        self.assertEqual(iters_at_max, 200000)

        # Acceptance == 0.0 triggers low-acceptance step reduction, not healthy branch
        _, step_zero, iters_zero, reason_zero = mfc.recommend_low_ess_params(
            bin_scale=1.0, bin_step=1.0, bin_iters=10000, burn_in=2000,
            eval_ess=13.0, ess_threshold=20.0, mean_tau=13000.0, bin_acceptance=0.0,
        )
        self.assertAlmostEqual(step_zero, 0.6)
        self.assertEqual(iters_zero, 12000)
        self.assertIn("low acceptance", reason_zero)

    def test_gibbs_recommendation_min_ess_floor_and_forwarding(self):
        import MomentFitCollector as mfc
        from moments2amplitudes.momentfit.runner import run_collector

        self.assertEqual(mfc.recommend_retry_gibbs_block_size(8, 0, 0.0, "LOW_ESS"), 16)
        self.assertEqual(mfc.recommend_retry_gibbs_block_size(8, 0, None, "ZOMBIE_ROOT"), 16)
        self.assertEqual(mfc.recommend_retry_gibbs_block_size(6, 0, 0.02, "LOW_ESS"), 12)
        self.assertEqual(mfc.recommend_retry_gibbs_block_size(4, 0, 0.0, "STUCK_COVARIANCE"), 8)
        self.assertEqual(mfc.recommend_retry_gibbs_block_size(4, 0, 0.07, "LOW_ESS"), 0)

        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            tbin_dir = run_dir / "t010020"
            tbin_dir.mkdir()
            (tbin_dir / "Mpi0eta1.060000_").mkdir()
            root_path = tbin_dir / "Mpi0eta1.060000_" / "ResultsBruMcmcFixedCovariance.root"
            root_path.write_bytes(b"root")

            prepared_scan = {
                "cache_path": str(tbin_dir / mfc._CACHE_FILENAME),
                "cache_data": {"version": mfc._CACHE_VERSION, "bins": {}},
                "pending_tasks": [],
                "prepared_bins": [
                    {
                        "b_idx": 0,
                        "b_low": 1.04,
                        "b_high": 1.08,
                        "b_center": 1.06,
                        "bin_dir_path": str(tbin_dir / "Mpi0eta1.060000_"),
                        "bin_dir_name": "Mpi0eta1.060000_",
                        "cand_files": [(str(tbin_dir / "Mpi0eta1.060000_"), str(root_path))],
                        "cache_hit": {
                            "root_file": str(root_path),
                            "zombie_file": None,
                            "found_dir": str(tbin_dir / "Mpi0eta1.060000_"),
                            "sig": (4, 1),
                            "scale_ok": True,
                            "scale_ratio": 1.0,
                            "scale_msg": "ok",
                            "diag": {
                                "ess_mean": 72.4,
                                "ess_min": 4.4,
                                "tau_mean": 23952.4,
                                "tau_max": 150352.4,
                                "acceptance": 0.009,
                                "per_moment_ess": {"H_0_0_0": 150.0, "H_0_4_0": 4.4},
                            },
                        },
                        "active_log_path": str(tbin_dir / "logRooFit0.txt"),
                        "retry_log": None,
                        "diag_mtime": None,
                        "pd": None,
                    }
                ],
            }

            # Passes mean >= 20, fails min_ess_threshold=10 -> flagged LOW_ESS with Gibbs=12 at max_l=6
            evals = mfc.evaluate_tbin(
                tbin_idx=0, tbin_name="t010020", tbin_dir=str(tbin_dir), run_dir=str(run_dir),
                m_min=1.04, m_max=1.08, n_bins=1, ess_threshold=20.0, ess_metric="mean",
                min_ess_threshold=10.0, burn_in=2000, current_scale=1.0, current_step=1.0,
                current_iters=10000, current_gibbs=0, max_l=6, prepared_scan=prepared_scan,
                precomputed_root_results={},
            )
            self.assertEqual(len(evals), 1)
            self.assertFalse(evals[0]["is_successful"])
            self.assertEqual(evals[0]["failure_type"], "LOW_ESS")
            self.assertEqual(evals[0]["rec_gibbs"], 12)
            self.assertIn("floor failed", evals[0]["diagnosis_reason"])

            script_path = mfc.generate_retry_script(
                run_dir=str(run_dir),
                unsuccessful_bins=evals,
                macro_path=str(REPO_ROOT / "macros" / "MomentFit.C"),
                max_l=6, m_min=1.04, m_max=1.08, n_bins=1, burn_in=2000,
                ess_threshold=20.0, min_ess_threshold=10.0,
                cov_dir="/cov/dir", data_dir="/data/dir",
                apply_mc_variance=True, allow_negative_pdf=True, use_tuning=True,
                cold_warmup_acc=0.2,
                data_pattern="{data_dir}/data_merged_{tbin}.root",
            )
            script_text = Path(script_path).read_text(encoding="utf-8")
            self.assertIn('/cov/dir', script_text)
            self.assertIn('/data/dir/data_merged_t010020.root', script_text)
            self.assertIn(', true, true, true, 12, 0.2,', script_text)
            self.assertIn('--min-ess-threshold "10.0"', script_text)

            cfg = load_momentfit_config(str(REPO_ROOT / "configs" / "momentfit_test.yaml"))
            with mock.patch("builtins.print") as mock_print:
                ok = run_collector(cfg=cfg, tbins_filter=["t010020"], max_l_filter=[6], no_exec=True, dry_run=True)
            self.assertTrue(ok)
            printed = "\n".join(" ".join(str(a) for a in call.args) for call in mock_print.call_args_list)
            self.assertIn("--min-ess-threshold 10.0", printed)
            self.assertIn("--cov-dir", printed)
            self.assertIn("--data-dir", printed)
            self.assertIn("--gibbs-block-size 0", printed)

    def test_huge_weight_chain_diagnostics_bounded_memory(self):
        import MomentFitCollector as mfc

        rng = np.random.default_rng(99)
        n_tree = 2000
        weights = np.full(n_tree, 5000.0, dtype=np.float64)  # n_raw = 10,000,000
        h0 = np.cumsum(rng.normal(0.0, 0.1, size=n_tree))
        data = {f"H_0_{k}_0": h0 + 0.01 * k for k in range(25)}
        data["weight"] = weights
        cols = [f"H_0_{k}_0" for k in range(25)]

        d_bfr = bfr.compute_batch_chain_diagnostics(data, cols, burn_in=2000)
        d_mfc = mfc._compute_batch_chain_diagnostics(data, cols, burn_in=2000)

        self.assertIsNotNone(d_bfr)
        self.assertIsNotNone(d_mfc)
        self.assertEqual(d_bfr["n_raw"], 10_000_000)
        self.assertEqual(d_mfc["n_raw"], 10_000_000)
        self.assertAlmostEqual(d_bfr["tau_mean"], d_mfc["tau_mean"], places=6)
        self.assertAlmostEqual(d_bfr["ess_mean"], d_mfc["ess_mean"], places=6)
        self.assertGreater(d_mfc["tau_mean"], 1000.0)


if __name__ == "__main__":
    unittest.main()

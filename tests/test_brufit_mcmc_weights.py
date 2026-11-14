import os
import tempfile
import unittest
import numpy as np

try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    HAS_ROOT = True
except ImportError:
    HAS_ROOT = False

import BruFitResults as BFR


@unittest.skipUnless(HAS_ROOT, "ROOT is required for BruFitResults tests")
class TestBruFitMCMCWeights(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root_file = os.path.join(self.temp_dir.name, "ResultsBruMcmcFixedCovariance.root")

        # Generate synthetic MCMC data with non-trivial weights
        np.random.seed(42)
        self.n_entries = 200
        self.weights = np.random.randint(1, 6, size=self.n_entries).astype(np.float64)
        self.h010 = np.random.normal(0.5, 0.1, size=self.n_entries).astype(np.float64)
        self.h020 = np.random.normal(-0.3, 0.08, size=self.n_entries).astype(np.float64)
        self.yld = np.random.normal(1000.0, 50.0, size=self.n_entries).astype(np.float64)

        # Create ROOT file with MCMCTree
        f = ROOT.TFile.Open(self.root_file, "RECREATE")
        tree = ROOT.TTree("MCMCTree", "MCMCTree")

        entry_arr = np.zeros(1, dtype=np.int64)
        w_arr = np.zeros(1, dtype=np.float64)
        yld_arr = np.zeros(1, dtype=np.float64)
        h010_arr = np.zeros(1, dtype=np.float64)
        h020_arr = np.zeros(1, dtype=np.float64)

        tree.Branch("entry", entry_arr, "entry/L")
        tree.Branch("weight", w_arr, "weight/D")
        tree.Branch("Yld_Moments", yld_arr, "Yld_Moments/D")
        tree.Branch("H_0_1_0", h010_arr, "H_0_1_0/D")
        tree.Branch("H_0_2_0", h020_arr, "H_0_2_0/D")

        for i in range(self.n_entries):
            entry_arr[0] = i
            w_arr[0] = self.weights[i]
            yld_arr[0] = self.yld[i]
            h010_arr[0] = self.h010[i]
            h020_arr[0] = self.h020[i]
            tree.Fill()

        tree.Write()
        f.Close()

        self.mcmc = BFR.MomentsMCMC(
            fit_result_dir=self.temp_dir.name,
            root_filename="ResultsBruMcmcFixedCovariance.root",
            burn_in=10
        )

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_weighted_mean_and_standard_deviation(self):
        """Test that get_normalized_moments_stats returns weighted mean and std."""
        mean, rms = self.mcmc.get_normalized_moments_stats(0, "H_0_1_0", burn_in=0)

        # Expected weighted mean and RMS
        expected_mean = float(np.average(self.h010, weights=self.weights))
        expected_var = float(np.average((self.h010 - expected_mean)**2, weights=self.weights))
        expected_rms = float(np.sqrt(expected_var))

        self.assertAlmostEqual(mean, expected_mean, places=3)
        self.assertAlmostEqual(rms, expected_rms, places=3)

    def test_weighted_posterior_distribution(self):
        """Test that get_normalized_moments_distribution uses weights."""
        hist = self.mcmc.get_normalized_moments_distribution(0, "H_0_1_0", burn_in=0)
        self.assertIsNotNone(hist)
        # Integral should equal total sum of weights
        total_w = float(np.sum(self.weights))
        self.assertAlmostEqual(hist.Integral(), total_w, delta=1e-5)

        expected_mean = float(np.average(self.h010, weights=self.weights))
        self.assertAlmostEqual(hist.GetMean(), expected_mean, places=3)

    def test_weighted_trace_graph(self):
        """Test that get_mcmc_trace_graph builds steps based on weights."""
        graph = self.mcmc.get_mcmc_trace_graph(0, "H_0_1_0")
        self.assertIsNotNone(graph)
        total_steps = int(np.sum(self.weights))
        self.assertEqual(graph.GetN(), total_steps)

        # First entry has weight self.weights[0]
        w0 = int(self.weights[0])
        for step in range(w0):
            x = graph.GetPointX(step)
            y = graph.GetPointY(step)
            self.assertEqual(x, float(step))
            self.assertAlmostEqual(y, self.h010[0], places=5)

    def test_weighted_ess_and_autocorrelation(self):
        """Test that ESS and integrated autocorrelation time respect weights."""
        tau = self.mcmc.get_mcmc_integrated_autocorrelation_time(0, "H_0_1_0", burn_in=0)
        ess = self.mcmc.get_mcmc_ess(0, "H_0_1_0", burn_in=0)
        total_steps = float(np.sum(self.weights))

        self.assertGreater(tau, 0.0)
        self.assertGreater(ess, 0.0)
        self.assertLessEqual(ess, total_steps)

    def test_weighted_unnormalized_moments_stats(self):
        """Test that get_unnormalized_moments_stats computes weighted mean and error."""
        mean, err = self.mcmc.get_unnormalized_moments_stats(0, "H_0_1_0")
        prod = self.h010 * self.yld / 2.0
        expected_mean = float(np.average(prod, weights=self.weights))
        expected_var = float(np.average((prod - expected_mean)**2, weights=self.weights))
        expected_err = float(np.sqrt(expected_var))

        self.assertAlmostEqual(mean, expected_mean, places=2)
        self.assertAlmostEqual(err, expected_err, places=2)

    def test_compute_chain_diagnostics_from_root(self):
        """Test compute_chain_diagnostics_from_root with weighted tree."""
        diags = BFR.compute_chain_diagnostics_from_root(self.root_file, burn_in=0)
        self.assertIsNotNone(diags)
        total_steps = int(np.sum(self.weights))
        self.assertEqual(diags["n_raw"], total_steps)
        self.assertIn("H_0_1_0", diags["per_moment"])
        self.assertIn("H_0_2_0", diags["per_moment"])


if __name__ == "__main__":
    unittest.main()

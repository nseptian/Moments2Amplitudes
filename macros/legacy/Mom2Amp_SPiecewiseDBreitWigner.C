#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include "TStopwatch.h"
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>

const Double_t FIRST_MASS_BIN_CENTER = 0.82; // in GeV
const Double_t MASS_BIN_WIDTH = 0.04;
const Int_t N_MASS_BINS = 20;

TString FIT_RESULTS_DIR = "/d/grid17/septian/EtaPi0Analysis/run/MomentFit_m080200/t010020/";

/**
 * @brief Helper to locate the best D-wave result ROOT file
 */
TString ResolveDWaveFile(TString userPath) {
    if (!userPath.IsNull() && gSystem->AccessPathName(userPath) == 0) {
        return userPath;
    }
    std::vector<TString> candidates = {
        userPath,
        "result_tree_DBreitWigner_best.root",
        "macros/result_tree_DBreitWigner_best.root",
        "logs/central_best_params.root",
        "macros/logs/central_best_params.root",
        "result_tree_test_0.root",
        "macros/result_tree_test_0.root"
    };
    for (const auto& path : candidates) {
        if (!path.IsNull() && gSystem->AccessPathName(path) == 0) {
            return path;
        }
    }
    return "";
}

/**
 * @brief S-wave Piecewise + D-wave Breit-Wigner fit with randomized starts & MCMC uncertainty propagation
 * @param nRandomStarts Number of randomized starts for central S-wave (default: 100). If <= 1, runs single seed.
 * @param burnIn Number of initial MCMC samples to skip (default: 0)
 * @param maxEntries Maximum number of MCMC samples to evaluate (-1 for all, 0 to skip MCMC)
 * @param stride Thinning stride (1 = every sample, 5 = every 5th sample, etc.)
 * @param dWaveFile ROOT file containing best D-wave fit parameters from Mom2Amp_DBreitWigner.C
 * @param outputFile Output ROOT file containing central trees and MCMC propagation trees
 * @param logDir Directory to store logs and diagnostic summaries (default: "logs_spiecewise")
 * @param nParallelJobs Number of parallel worker jobs (default: 40)
 */
void Mom2Amp_SPiecewiseDBreitWigner(int nRandomStarts = 100,
                                   Long64_t burnIn = 0,
                                   Long64_t maxEntries = -1,
                                   int stride = 1,
                                   TString dWaveFile = "result_tree_DBreitWigner_best.root",
                                   TString outputFile = "mcmc_propagation_spiecewise_results.root",
                                   TString logDir = "logs_spiecewise",
                                   int nParallelJobs = 40) {

    auto& setup = ConfigureAmpsNoValues(2, 2, 2); // (Lmax, MMax, Nref)
    
    // 1. Load moments data with error handling
    MassDependentMoments massDepMoments;
    try {
        massDepMoments = MassDependentMoments(N_MASS_BINS, 
                                             MASS_BIN_WIDTH, 
                                             FIRST_MASS_BIN_CENTER, 
                                             FIT_RESULTS_DIR,
                                             "ResultsBruMcmcCovariance.root", 
                                             "Mpi0eta");
    } catch (const std::exception& e) {
        std::cerr << "\n=== ERROR ===" << std::endl;
        std::cerr << "Failed to load moments data: " << e.what() << std::endl;
        std::cerr << "Please check directory: " << FIT_RESULTS_DIR << std::endl;
        return;
    }

    std::vector<double> massBins = massDepMoments.GetMassBins();

    using namespace m2pw;
    using MDModel = MassDependentFitter::MassDependenceConfig::ModelType;
    using WaveModelConfig = MassDependentFitter::MassDependenceConfig::WaveModelConfig;
    
    auto a2_1320_config = MassDependentFitter::CreateCustomConfig(
        {{2, {WaveModelConfig{"a2_1320", MDModel::BreitWignerPlusPolynomial, -1}}}}, 
        {0, 1}
    );
    auto H024_config = MassDependentFitter::CreateL0L2L4OnlyConfig();

    MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, a2_1320_config, H024_config};
    fitter.SetEquationValues(massDepMoments);
    fitter.PrintIncludedMoments();

    // 2. Resolve and verify D-wave initial parameter file
    TString resolvedDWaveFile = ResolveDWaveFile(dWaveFile);
    if (!resolvedDWaveFile.IsNull()) {
        std::cout << "\n========================================================" << std::endl;
        std::cout << " Found D-wave best fit file: " << resolvedDWaveFile << std::endl;
        std::cout << " Initializing D-wave parameters from this file." << std::endl;
        std::cout << " S-wave amplitudes randomized; S-wave phases set to zero." << std::endl;
        std::cout << "========================================================\n" << std::endl;
    } else {
        std::cerr << "\n[Warning] D-wave result file not found: " << dWaveFile << std::endl;
        std::cerr << "Will use standard Breit-Wigner default parameters for D-wave.\n" << std::endl;
    }

    gSystem->mkdir(logDir.Data(), true);
    TString logDirCentral = Form("%s/central_fits", logDir.Data());
    gSystem->mkdir(logDirCentral.Data(), true);

    // ========================================================================
    // Step 1: Finding Global Minimum on Central Moments
    // ========================================================================
    std::cout << "\n==========================================" << std::endl;
    std::cout << " Step 1: Finding Global Minimum on Central Moments (S-Piecewise + D-BW)" << std::endl;
    std::cout << "==========================================" << std::endl;

    double bestChi2 = 1e12;
    int bestSeed = -1;
    std::map<TString, double> bestParams;

    struct FitRecord {
        int seed = 0;
        double chi2 = 0.0;
        int status = -1;
        bool valid = false;
        std::map<TString, double> params;
    };
    std::vector<FitRecord> fitRecords;

    fitter.SetMinimizerTolerance(0.01);

    if (nRandomStarts > 1) {
        if (nParallelJobs > 1) {
            std::cout << Form("Running %d S-piecewise randomized fits in parallel across %d workers...",
                              nRandomStarts, nParallelJobs) << std::endl;
            std::cout << "Fit logs are being saved to: " << logDirCentral << "/\n" << std::endl;

            TString runScriptDir = gSystem->DirName(gSystem->Which(gROOT->GetMacroPath(), "Mom2Amp_SPiecewiseDBreitWigner.C"));
            if (runScriptDir.IsNull() || gSystem->AccessPathName(runScriptDir + "/RunParallelSPiecewiseFits.sh") != 0) {
                if (gSystem->AccessPathName("macros/RunParallelSPiecewiseFits.sh") == 0) {
                    runScriptDir = "macros";
                } else {
                    runScriptDir = ".";
                }
            }
            TString cmd = Form("bash %s/RunParallelSPiecewiseFits.sh %d %d %s %s",
                               runScriptDir.Data(), nRandomStarts, nParallelJobs,
                               (resolvedDWaveFile.IsNull() ? "\"\"" : resolvedDWaveFile.Data()),
                               logDirCentral.Data());

            int ret = gSystem->Exec(cmd.Data());
            if (ret != 0) {
                std::cerr << "Warning: Parallel fit script returned non-zero exit code: " << ret << std::endl;
            }

            for (int seed = 1; seed <= nRandomStarts; ++seed) {
                TString resFile = Form("%s/result_seed_%d.root", logDirCentral.Data(), seed);
                TString logFile = Form("%s/log_seed_%d.txt", logDirCentral.Data(), seed);
                TFile* fRes = TFile::Open(resFile);
                if (fRes && !fRes->IsZombie()) {
                    TTree* tChi2 = fRes->Get<TTree>("chi2");
                    TTree* tAll = fRes->Get<TTree>("all_params");
                    if (!tAll) tAll = fRes->Get<TTree>("mass_dependent_params");

                    if (tChi2) {
                        double chi2 = 0.0;
                        int status = -1;
                        bool valid = false;
                        tChi2->SetBranchAddress("chi2", &chi2);
                        tChi2->SetBranchAddress("Status", &status);
                        tChi2->SetBranchAddress("isValid", &valid);
                        tChi2->GetEntry(0);

                        FitRecord rec;
                        rec.seed = seed;
                        rec.chi2 = chi2;
                        rec.status = status;
                        rec.valid = valid;

                        if (tAll) {
                            TObjArray* branches = tAll->GetListOfBranches();
                            for (int b = 0; b < branches->GetEntries(); ++b) {
                                TString bName = branches->At(b)->GetName();
                                if (!bName.EndsWith("_error")) {
                                    double val = 0.0;
                                    tAll->SetBranchAddress(bName, &val);
                                    tAll->GetEntry(0);
                                    rec.params[bName] = val;
                                }
                            }
                        }
                        fitRecords.push_back(rec);

                        std::cout << Form("  [Seed %3d] Chi2 = %10.2f | Status = %d | Valid = %d | Log: %s%s",
                                          seed, chi2, status, (int)valid, logFile.Data(),
                                          (valid && chi2 < bestChi2 ? "  <-- new lowest chi2!" : ""))
                                  << std::endl;

                        if (valid && chi2 < bestChi2) {
                            bestChi2 = chi2;
                            bestSeed = seed;
                            bestParams = rec.params;
                        }
                    }
                    fRes->Close();
                    delete fRes;
                } else {
                    if (fRes) delete fRes;
                    std::cout << Form("  [Seed %3d] Failed to read result file %s (check log %s)",
                                      seed, resFile.Data(), logFile.Data()) << std::endl;
                }
            }
        } else {
            // Sequential execution
            std::cout << Form("Running %d S-piecewise randomized fits sequentially...", nRandomStarts) << std::endl;
            std::cout << "Fit logs are being saved to: " << logDirCentral << "/\n" << std::endl;

            for (int seed = 1; seed <= nRandomStarts; ++seed) {
                TString logFile = Form("%s/log_seed_%d.txt", logDirCentral.Data(), seed);
                gSystem->RedirectOutput(logFile.Data(), "w");

                MassDependentFitter::ParameterManager trialManager(fitter);
                trialManager.AddMassIndependentParameter(std::vector<int>{0, 1}, seed);

                for (const double massBin : massBins) {
                    for (const TString& parName : MassDependentFitter::GetParNames()) {
                        if (parName.Contains("_1_")) {
                            TString fixedName = Form("MI_%1.6f_%s", massBin, parName.Data());
                            trialManager.SetInitialValue(fixedName, 0.0, true);
                        }
                    }
                }

                if (!resolvedDWaveFile.IsNull()) {
                    trialManager.AddMassDependentParameter(std::vector<int>{2}, resolvedDWaveFile, 0);
                    double dMass = 1.318;
                    double dWidth = 0.107;
                    if (trialManager.parsList.find("MD_shared_a2_1320_Mass") != trialManager.parsList.end()) {
                        dMass = trialManager.parsList["MD_shared_a2_1320_Mass"];
                    }
                    if (trialManager.parsList.find("MD_shared_a2_1320_Width") != trialManager.parsList.end()) {
                        dWidth = trialManager.parsList["MD_shared_a2_1320_Width"];
                    }
                    trialManager.SetInitialValue("MD_shared_a2_1320_Mass", dMass, true);
                    trialManager.SetInitialValue("MD_shared_a2_1320_Width", dWidth, true);
                } else {
                    trialManager.AddMassDependentParameter(std::vector<int>{2}, 0, 0);
                    trialManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
                    trialManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);
                }

                TRandom3 rng(seed > 0 ? seed : 12345);
                for (const double massBin : massBins) {
                    for (const TString& parName : MassDependentFitter::GetParNames()) {
                        if (parName.Contains("_0_0")) {
                            TString sName = Form("MI_%1.6f_%s", massBin, parName.Data());
                            if (trialManager.parsList.find(sName) != trialManager.parsList.end()) {
                                if (sName.Contains("phi")) {
                                    trialManager.SetInitialValue(sName, 0.0, false);
                                } else {
                                    double rMag = rng.Uniform(0.0, 200.0);
                                    trialManager.SetInitialValue(sName, rMag, false);
                                }
                            }
                        }
                    }
                }

                fitter.SetMinimizerPrintLevel(1);
                fitter.MinimizeChi2(trialManager);

                double chi2 = fitter.GetLastChi2();
                int status = fitter.GetMinimizerStatus();
                bool valid = fitter.IsMinimizerValid();

                gSystem->RedirectOutput(0);

                FitRecord rec;
                rec.seed = seed;
                rec.chi2 = chi2;
                rec.status = status;
                rec.valid = valid;
                rec.params = trialManager.parsList;
                fitRecords.push_back(rec);

                std::cout << Form("  [Seed %3d] Chi2 = %10.2f | Status = %d | Valid = %d | Log: %s%s",
                                  seed, chi2, status, (int)valid, logFile.Data(),
                                  (valid && chi2 < bestChi2 ? "  <-- new lowest chi2!" : ""))
                          << std::endl;

                if (valid && chi2 < bestChi2) {
                    bestChi2 = chi2;
                    bestSeed = seed;
                    bestParams = trialManager.parsList;
                }
            }
        }

        // Summary analysis of convergence and global minimum cluster
        std::vector<FitRecord> validFits;
        for (const auto& r : fitRecords) {
            if (r.valid) validFits.push_back(r);
        }
        std::sort(validFits.begin(), validFits.end(), [](const FitRecord& a, const FitRecord& b) {
            return a.chi2 < b.chi2;
        });

        const double clusterThreshold = 0.05;
        std::vector<int> globalMinimaSeeds;
        for (const auto& r : validFits) {
            if (std::abs(r.chi2 - bestChi2) <= clusterThreshold) {
                globalMinimaSeeds.push_back(r.seed);
            }
        }

        std::cout << "\n========================================================" << std::endl;
        std::cout << "         S-PIECEWISE RANDOMIZED FITS SUMMARY            " << std::endl;
        std::cout << "========================================================" << std::endl;
        std::cout << Form(" Total starts:     %d", nRandomStarts) << std::endl;
        std::cout << Form(" Valid converged:  %zu (%.1f%%)", validFits.size(), (validFits.size() * 100.0 / nRandomStarts)) << std::endl;
        std::cout << Form(" Lowest Chi2:      %.2f (Seed %d)", bestChi2, bestSeed) << std::endl;
        std::cout << Form(" Global Minima Cluster: %zu fits within Delta Chi2 <= %.2f (%.1f%% of valid)",
                          globalMinimaSeeds.size(), clusterThreshold,
                          (validFits.empty() ? 0.0 : globalMinimaSeeds.size() * 100.0 / validFits.size())) << std::endl;
        std::cout << " Seeds in cluster: ";
        for (size_t i = 0; i < globalMinimaSeeds.size(); ++i) {
            std::cout << globalMinimaSeeds[i] << (i + 1 < globalMinimaSeeds.size() ? ", " : "\n");
        }
        std::cout << "--------------------------------------------------------" << std::endl;
        std::cout << " Top Lowest Chi2 Fits:" << std::endl;
        std::cout << "  Rank |  Seed |       Chi2 | Delta Chi2 | Status | Interpretation" << std::endl;
        std::cout << " --------------------------------------------------------" << std::endl;
        const size_t nToShow = std::min((size_t)10, validFits.size());
        for (size_t i = 0; i < nToShow; ++i) {
            double dChi2 = validFits[i].chi2 - bestChi2;
            TString interp = (dChi2 <= clusterThreshold) ? "Global Minimum" : "Local Minimum";
            std::cout << Form("   %3zu |   %4d | %10.2f |    %+7.2f |    %3d | %s",
                              i + 1, validFits[i].seed, validFits[i].chi2, dChi2,
                              validFits[i].status, interp.Data())
                      << std::endl;
        }
        std::cout << "========================================================\n" << std::endl;

        // Save detailed text summary
        TString summaryPath = Form("%s/spiecewise_fits_summary.txt", logDir.Data());
        std::ofstream sumFile(summaryPath.Data());
        if (sumFile.is_open()) {
            sumFile << "========================================================" << std::endl;
            sumFile << "         S-PIECEWISE RANDOMIZED FITS SUMMARY            " << std::endl;
            sumFile << "========================================================" << std::endl;
            sumFile << " Total starts:     " << nRandomStarts << std::endl;
            sumFile << " Valid converged:  " << validFits.size() << " (" << std::fixed << std::setprecision(1) << (validFits.size() * 100.0 / nRandomStarts) << "%)" << std::endl;
            sumFile << " Lowest Chi2:      " << std::fixed << std::setprecision(2) << bestChi2 << " (Seed " << bestSeed << ")" << std::endl;
            sumFile << " Global Minima Cluster: " << globalMinimaSeeds.size() << " fits within Delta Chi2 <= " << clusterThreshold << std::endl;
            sumFile << " Seeds in cluster: ";
            for (size_t i = 0; i < globalMinimaSeeds.size(); ++i) {
                sumFile << globalMinimaSeeds[i] << (i + 1 < globalMinimaSeeds.size() ? ", " : "\n");
            }
            sumFile << "--------------------------------------------------------" << std::endl;
            sumFile << " Ranked Fits:" << std::endl;
            for (size_t i = 0; i < validFits.size(); ++i) {
                double dChi2 = validFits[i].chi2 - bestChi2;
                std::string interp = (dChi2 <= clusterThreshold) ? "Global Minimum" : "Local Minimum";
                sumFile << std::setw(6) << (i + 1) << " | "
                        << std::setw(5) << validFits[i].seed << " | "
                        << std::setw(10) << std::fixed << std::setprecision(2) << validFits[i].chi2 << " | "
                        << std::setw(10) << std::showpos << dChi2 << std::noshowpos << " | "
                        << std::setw(6) << validFits[i].status << " | "
                        << interp << std::endl;
            }
            sumFile << "========================================================" << std::endl;
            sumFile.close();
            std::cout << "[Diagnostic] Full summary log saved to " << summaryPath << std::endl;
        }

        // Draw Chi2 Distribution and Plateau Plot
        if (!validFits.empty()) {
            gStyle->SetOptStat(1111);
            TCanvas* cChi2 = new TCanvas("c_spiecewise_chi2_dist", "S-Piecewise Starts Chi2 Distribution", 1200, 500);
            cChi2->Divide(2, 1);

            // Pad 1: Sorted Chi2 Plateau Plot
            cChi2->cd(1);
            gPad->SetGrid();
            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.06);

            const int nValid = validFits.size();
            std::vector<double> ranks(nValid);
            std::vector<double> chi2Vals(nValid);
            for (int i = 0; i < nValid; ++i) {
                ranks[i] = i + 1;
                chi2Vals[i] = validFits[i].chi2;
            }

            TGraph* grPlateau = new TGraph(nValid, ranks.data(), chi2Vals.data());
            grPlateau->SetName("gr_spiecewise_chi2_plateau");
            grPlateau->SetTitle("Sorted #chi^{2} vs Fit Rank (S-Piecewise Plateau);Fit Rank (ordered by #chi^{2});#chi^{2}");
            grPlateau->SetMarkerStyle(20);
            grPlateau->SetMarkerSize(0.9);
            grPlateau->SetMarkerColor(kAzure + 1);
            grPlateau->SetLineColor(kAzure + 1);
            grPlateau->SetLineWidth(2);
            grPlateau->GetYaxis()->SetTitleOffset(1.3);
            grPlateau->Draw("ALP");

            TLine* lineBest = new TLine(1, bestChi2, nValid, bestChi2);
            lineBest->SetLineColor(kRed + 1);
            lineBest->SetLineStyle(2);
            lineBest->SetLineWidth(2);
            lineBest->Draw("same");

            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.035);
            latex.SetTextColor(kRed + 2);
            latex.DrawLatex(0.18, 0.82, Form("#chi^{2}_{min} = %.2f (Cluster: %zu / %d fits)", 
                                           bestChi2, globalMinimaSeeds.size(), nValid));

            // Pad 2: Histogram of Chi2
            cChi2->cd(2);
            gPad->SetGrid();
            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.06);

            double histMax = validFits.back().chi2;
            if (nValid >= 5 && histMax > validFits[std::min((size_t)nValid - 1, (size_t)25)].chi2 * 1.5) {
                histMax = validFits[std::min((size_t)nValid - 1, (size_t)25)].chi2 * 1.5;
            }
            if (histMax <= bestChi2) histMax = bestChi2 + 50.0;

            TH1D* hChi2 = new TH1D("h_spiecewise_chi2_dist", "Distribution of #chi^{2} Across S-Piecewise Starts;#chi^{2};Number of Fits", 
                                   35, bestChi2 * 0.98, histMax);
            hChi2->SetFillColorAlpha(kTeal + 3, 0.45);
            hChi2->SetLineColor(kTeal + 3);
            hChi2->SetLineWidth(2);
            hChi2->GetYaxis()->SetTitleOffset(1.3);

            for (const auto& r : validFits) {
                hChi2->Fill(r.chi2);
            }
            hChi2->Draw("HIST");

            cChi2->Update();
            cChi2->SaveAs(Form("%s/spiecewise_chi2_distribution.png", logDir.Data()));
            cChi2->SaveAs(Form("%s/spiecewise_chi2_distribution.pdf", logDir.Data()));
            cChi2->SaveAs("spiecewise_chi2_distribution.png");
            cChi2->SaveAs("spiecewise_chi2_distribution.pdf");
            std::cout << "[Diagnostic] Diagnostic plot saved: " << logDir << "/spiecewise_chi2_distribution.png and .pdf\n" << std::endl;
        }
    }

    // Run best fit on central moments and save central_best_params.root
    std::cout << "Setting up nominal central fit with best parameters (Seed " << (bestSeed > 0 ? bestSeed : 0) << ")..." << std::endl;
    MassDependentFitter::ParameterManager bestManager(fitter);
    bestManager.AddMassIndependentParameter(std::vector<int>{0, 1}, (bestSeed > 0 ? bestSeed : 0));

    for (const double massBin : massBins) {
        for (const TString& parName : MassDependentFitter::GetParNames()) {
            if (parName.Contains("_1_")) {
                TString fixedName = Form("MI_%1.6f_%s", massBin, parName.Data());
                bestManager.SetInitialValue(fixedName, 0.0, true);
            }
        }
    }

    if (!resolvedDWaveFile.IsNull()) {
        bestManager.AddMassDependentParameter(std::vector<int>{2}, resolvedDWaveFile, 0);
        double dMass = 1.318;
        double dWidth = 0.107;
        if (bestManager.parsList.find("MD_shared_a2_1320_Mass") != bestManager.parsList.end()) {
            dMass = bestManager.parsList["MD_shared_a2_1320_Mass"];
        }
        if (bestManager.parsList.find("MD_shared_a2_1320_Width") != bestManager.parsList.end()) {
            dWidth = bestManager.parsList["MD_shared_a2_1320_Width"];
        }
        bestManager.SetInitialValue("MD_shared_a2_1320_Mass", dMass, true);
        bestManager.SetInitialValue("MD_shared_a2_1320_Width", dWidth, true);
    } else {
        bestManager.AddMassDependentParameter(std::vector<int>{2}, 0, 0);
        bestManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
        bestManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);
    }

    if (!bestParams.empty()) {
        bestManager.SetParameterInitialValues(bestParams);
    } else {
        TRandom3 rng(12345);
        for (const double massBin : massBins) {
            for (const TString& parName : MassDependentFitter::GetParNames()) {
                if (parName.Contains("_0_0")) {
                    TString sName = Form("MI_%1.6f_%s", massBin, parName.Data());
                    if (bestManager.parsList.find(sName) != bestManager.parsList.end()) {
                        if (sName.Contains("phi")) {
                            bestManager.SetInitialValue(sName, 0.0, false);
                        } else {
                            double rMag = rng.Uniform(0.0, 200.0);
                            bestManager.SetInitialValue(sName, rMag, false);
                        }
                    }
                }
            }
        }
    }

    TString bestCentralLog = Form("%s/log_best_central_fit.txt", logDir.Data());
    gSystem->RedirectOutput(bestCentralLog.Data(), "w");
    fitter.SetMinimizerPrintLevel(2);
    fitter.MinimizeChi2(bestManager);
    gSystem->RedirectOutput(0);

    TString centralBestFile = Form("%s/central_best_params.root", logDir.Data());
    fitter.MakeResultTree(bestManager, centralBestFile);

    std::cout << "Nominal central fit finished with Chi2 = " << fitter.GetLastChi2() 
              << " | Log: " << bestCentralLog 
              << " | Central params file: " << centralBestFile << std::endl;

    // ========================================================================
    // Step 2: Propagating Uncertainty over MCMC Tree
    // ========================================================================
    if (maxEntries == 0) {
        std::cout << "\n[Info] maxEntries=0: Skipping MCMC uncertainty propagation." << std::endl;
        return;
    }

    std::cout << "\n==========================================" << std::endl;
    std::cout << " Step 2: Propagating Uncertainty over MCMC Tree (S-Piecewise + D-BW)" << std::endl;
    std::cout << " (Using best central fit as initial values for each sample)" << std::endl;
    std::cout << " burnIn=" << burnIn << ", maxEntries=" << maxEntries 
              << ", stride=" << stride << ", out=" << outputFile << std::endl;
    std::cout << "==========================================" << std::endl;

    TString logDirMCMC = Form("%s/mcmc_fits", logDir.Data());
    gSystem->mkdir(logDirMCMC.Data(), true);

    const Long64_t totalAvailable = massDepMoments.GetMCMCEntries();
    Long64_t actualEntriesToFit = (maxEntries > 0 && (burnIn + maxEntries <= totalAvailable))
                                  ? maxEntries
                                  : (totalAvailable - burnIn);
    if (actualEntriesToFit <= 0) actualEntriesToFit = 0;

    if (nParallelJobs > 1 && actualEntriesToFit > 0) {
        std::cout << Form("Running S-Piecewise MCMC propagation in parallel across %d workers (%lld entries)...",
                          nParallelJobs, actualEntriesToFit) << std::endl;
        std::cout << "Worker logs: " << logDirMCMC << "/log_worker_*.txt" << std::endl;
        std::cout << "Fit logs:    " << logDirMCMC << "/log_entry_*.txt\n" << std::endl;

        TString runScriptDir = gSystem->DirName(gSystem->Which(gROOT->GetMacroPath(), "Mom2Amp_SPiecewiseDBreitWigner.C"));
        if (runScriptDir.IsNull() || gSystem->AccessPathName(runScriptDir + "/RunParallelMCMCFitsSPiecewise.sh") != 0) {
            if (gSystem->AccessPathName("macros/RunParallelMCMCFitsSPiecewise.sh") == 0) {
                runScriptDir = "macros";
            } else {
                runScriptDir = ".";
            }
        }
        TString cmd = Form("bash %s/RunParallelMCMCFitsSPiecewise.sh %d %lld %lld %d %s %s",
                           runScriptDir.Data(), nParallelJobs, burnIn, actualEntriesToFit, stride,
                           centralBestFile.Data(), logDirMCMC.Data());

        int ret = gSystem->Exec(cmd.Data());
        if (ret != 0) {
            std::cerr << "Warning: Parallel S-piecewise MCMC script returned exit code: " << ret << std::endl;
        }

        // Merge chunk files into outputFile
        std::cout << "\nMerging parallel worker chunks into " << outputFile << "..." << std::endl;
        TFile* fOut = TFile::Open(outputFile, "RECREATE");
        if (!fOut || fOut->IsZombie()) {
            std::cerr << "Error: Cannot open output file " << outputFile << " for merging!" << std::endl;
            if (fOut) delete fOut;
            return;
        }

        // Copy central trees
        TFile* fCent = TFile::Open(centralBestFile, "READ");
        if (fCent && !fCent->IsZombie()) {
            TTree* tResult = fCent->Get<TTree>("result");
            if (tResult) {
                fOut->cd();
                TTree* cResult = tResult->CloneTree(-1, "fast");
                cResult->SetName("central_result");
                cResult->SetTitle("Nominal Central Amplitudes and Moments per Mass Bin");
                cResult->Write();
            }
            TTree* tAll = fCent->Get<TTree>("all_params");
            if (!tAll) tAll = fCent->Get<TTree>("mass_dependent_params");
            if (tAll) {
                fOut->cd();
                TTree* cParams = tAll->CloneTree(-1, "fast");
                cParams->SetName("central_params");
                cParams->SetTitle("Nominal Central Fitted Parameters");
                cParams->Write();
            }
            TTree* tChi2 = fCent->Get<TTree>("chi2");
            if (tChi2) {
                fOut->cd();
                TTree* cChi2 = tChi2->CloneTree(-1, "fast");
                cChi2->SetName("central_chi2");
                cChi2->SetTitle("Nominal Central Chi2 and Status");
                cChi2->Write();
            }
            fCent->Close();
            delete fCent;
        }

        // Chain and merge MCMC trees
        TChain chainParams("mcmc_params");
        TChain chainChi2("mcmc_chi2");
        TChain chainResult("mcmc_result");
        int foundChunks = 0;

        for (int w = 0; w < nParallelJobs; ++w) {
            TString chunkPath = Form("%s/chunk_worker_%d.root", logDirMCMC.Data(), w);
            if (gSystem->AccessPathName(chunkPath) == 0) {
                chainParams.Add(chunkPath);
                chainChi2.Add(chunkPath);
                chainResult.Add(chunkPath);
                foundChunks++;
            }
        }

        std::cout << Form("Found %d chunk files with %lld total MCMC fitted entries.",
                          foundChunks, chainParams.GetEntries()) << std::endl;

        fOut->cd();
        if (chainParams.GetEntries() > 0) {
            TTree* mParams = chainParams.CloneTree(0);
            mParams->CopyEntries(&chainParams);
            mParams->Write("mcmc_params");
        }
        if (chainChi2.GetEntries() > 0) {
            TTree* mChi2 = chainChi2.CloneTree(0);
            mChi2->CopyEntries(&chainChi2);
            mChi2->Write("mcmc_chi2");
        }
        if (chainResult.GetEntries() > 0) {
            TTree* mResult = chainResult.CloneTree(0);
            mResult->CopyEntries(&chainResult);
            mResult->Write("mcmc_result");
        }

        // Summary analysis on MCMC fits
        Long64_t nEntriesFitted = chainChi2.GetEntries();
        if (nEntriesFitted > 0) {
            double c_val = 0.0;
            bool v_val = false;
            int s_val = -1;
            Long64_t e_val = 0;
            chainChi2.SetBranchAddress("chi2", &c_val);
            chainChi2.SetBranchAddress("isValid", &v_val);
            chainChi2.SetBranchAddress("Status", &s_val);
            chainChi2.SetBranchAddress("mcmc_entry", &e_val);

            Long64_t nValid = 0;
            double sumChi2 = 0.0;
            double minChi2 = 1e12;
            double maxChi2 = -1e12;

            TString mcmcSumPath = Form("%s/mcmc_fits_summary.txt", logDir.Data());
            std::ofstream mSumFile(mcmcSumPath.Data());
            if (mSumFile.is_open()) {
                mSumFile << "========================================================" << std::endl;
                mSumFile << "     S-PIECEWISE MCMC PROPAGATION FITS SUMMARY          " << std::endl;
                mSumFile << "========================================================" << std::endl;
                mSumFile << " Total entries fitted: " << nEntriesFitted << std::endl;
                mSumFile << "  Entry |       Chi2 | Status | Valid | Log File" << std::endl;
                mSumFile << " --------------------------------------------------------" << std::endl;
            }

            for (Long64_t i = 0; i < nEntriesFitted; ++i) {
                chainChi2.GetEntry(i);
                if (v_val) {
                    nValid++;
                    sumChi2 += c_val;
                    if (c_val < minChi2) minChi2 = c_val;
                    if (c_val > maxChi2) maxChi2 = c_val;
                }
                if (mSumFile.is_open()) {
                    mSumFile << std::setw(7) << e_val << " | "
                             << std::setw(10) << std::fixed << std::setprecision(2) << c_val << " | "
                             << std::setw(6) << s_val << " | "
                             << std::setw(5) << (int)v_val << " | "
                             << Form("%s/log_entry_%lld.txt", logDirMCMC.Data(), e_val)
                             << std::endl;
                }
            }
            if (mSumFile.is_open()) {
                mSumFile << "========================================================" << std::endl;
                mSumFile.close();
            }

            std::cout << "\n========================================================" << std::endl;
            std::cout << "     S-PIECEWISE MCMC PROPAGATION FITS SUMMARY          " << std::endl;
            std::cout << "========================================================" << std::endl;
            std::cout << Form(" Total fitted:    %lld", nEntriesFitted) << std::endl;
            std::cout << Form(" Valid converged: %lld (%.1f%%)", nValid, (nValid * 100.0 / nEntriesFitted)) << std::endl;
            if (nValid > 0) {
                std::cout << Form(" Chi2 range:      [%.2f, %.2f] (mean: %.2f)", minChi2, maxChi2, sumChi2 / nValid) << std::endl;
            }
            std::cout << Form(" Detailed summary saved to: %s", mcmcSumPath.Data()) << std::endl;
            std::cout << "========================================================\n" << std::endl;
        }

        // Save diagnostic canvas into the output ROOT file as well
        TCanvas* cChi2 = dynamic_cast<TCanvas*>(gROOT->FindObject("c_spiecewise_chi2_dist"));
        if (cChi2) {
            fOut->cd();
            cChi2->Write("c_chi2_dist", TObject::kOverwrite);
            std::cout << "[Diagnostic] Diagnostic canvas written into " << outputFile << std::endl;
        }
        fOut->Close();
        delete fOut;

    } else {
        // Sequential fallback
        fitter.FitMCMCTree(bestManager, massDepMoments, outputFile, burnIn, maxEntries, stride, logDirMCMC);

        TFile* fOut = TFile::Open(outputFile, "UPDATE");
        if (fOut && !fOut->IsZombie()) {
            TCanvas* cChi2 = dynamic_cast<TCanvas*>(gROOT->FindObject("c_spiecewise_chi2_dist"));
            if (cChi2) {
                cChi2->Write("c_chi2_dist", TObject::kOverwrite);
            }
            fOut->Close();
            delete fOut;
        }
    }
}

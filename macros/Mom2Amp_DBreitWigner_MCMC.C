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
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>

const Double_t FIRST_MASS_BIN_CENTER = 1.06; // in GeV
const Double_t MASS_BIN_WIDTH = 0.04;
const Int_t N_MASS_BINS = 14;

TString FIT_RESULTS_DIR = "/d/grid17/septian/EtaPi0Analysis/run/MomentFit_m080200/t010020/";

/**
 * @brief Run nominal central fit and propagate MCMC tree uncertainties
 * @param nRandomStarts Number of random seeds to search for global minimum (default: 100). If 1, runs single seed 0.
 * @param burnIn Number of initial MCMC samples to skip (default: 0)
 * @param maxEntries Maximum number of MCMC samples to evaluate (-1 for all)
 * @param stride Thinning stride (1 = every sample, 5 = every 5th sample, etc.)
 * @param outputFile Output ROOT file containing central and MCMC trees
 * @param logDir Directory to store individual fit logs and diagnostic summaries (default: "logs")
 * @param nParallelJobs Number of parallel worker jobs for the randomized starts (default: 40)
 */
void Mom2Amp_DBreitWigner_MCMC(int nRandomStarts = 100,
                              Long64_t burnIn = 0, 
                              Long64_t maxEntries = -1, 
                              int stride = 1, 
                              TString outputFile = "mcmc_propagation_results.root",
                              TString logDir = "logs",
                              int nParallelJobs = 40) {

    auto& setup = ConfigureAmpsNoValues(2, 2, 2); // (Lmax, MMax, Nref)
    
    // Load moments data with error handling (loads central moments and MCMC samples)
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

    std::cout << "\nLoaded " << massDepMoments.GetNMassBins() << " mass bins." << std::endl;
    std::cout << "Available MCMC entries: " << massDepMoments.GetMCMCEntries() << std::endl;

    if (!massDepMoments.HasMCMCSamples()) {
        std::cerr << "Error: No MCMC samples found in loaded moment files!" << std::endl;
        return;
    }

    std::vector<double> massBins = massDepMoments.GetMassBins();

    using namespace m2pw;
    using MDModel = MassDependentFitter::MassDependenceConfig::ModelType;
    using WaveModelConfig = MassDependentFitter::MassDependenceConfig::WaveModelConfig;
    
    auto a2_1320_config = MassDependentFitter::CreateCustomConfig(
        {{2, {WaveModelConfig{"a2_1320", MDModel::BreitWignerPlusPolynomial, -1}}}}, 
        {0}
    );
    auto H4_config = MassDependentFitter::CreateL4OnlyConfig();

    MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, a2_1320_config, H4_config};
    
    // Set nominal central equation values
    fitter.SetEquationValues(massDepMoments);
    fitter.PrintIncludedMoments();

    // 1. Search for global minimum across random seeds (or run single seed)
    std::cout << "\n==========================================" << std::endl;
    std::cout << " Step 1: Finding Global Minimum on Central Moments" << std::endl;
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

    // Create log directories
    gSystem->mkdir(logDir.Data(), true);
    TString logDirCentral = Form("%s/central_fits", logDir.Data());
    gSystem->mkdir(logDirCentral.Data(), true);

    if (nRandomStarts > 1) {
        if (nParallelJobs > 1) {
            std::cout << Form("Running %d randomized initial starts in parallel across %d workers...",
                              nRandomStarts, nParallelJobs) << std::endl;
            std::cout << "Fit logs are being saved to: " << logDirCentral << "/\n" << std::endl;

            TString runScript = gSystem->DirName(gSystem->Which(gROOT->GetMacroPath(), "Mom2Amp_DBreitWigner_MCMC.C"));
            if (runScript.IsNull() || gSystem->AccessPathName(runScript + "/RunParallelRandomizedFits.sh") != 0) {
                if (gSystem->AccessPathName("macros/RunParallelRandomizedFits.sh") == 0) {
                    runScript = "macros";
                } else {
                    runScript = ".";
                }
            }
            TString cmd = Form("bash %s/RunParallelRandomizedFits.sh %d %d %s",
                               runScript.Data(), nRandomStarts, nParallelJobs, logDirCentral.Data());
            
            int ret = gSystem->Exec(cmd.Data());
            if (ret != 0) {
                std::cerr << "Warning: Parallel fit script returned non-zero exit code: " << ret << std::endl;
            }

            // Collect results from output root files
            for (int seed = 1; seed <= nRandomStarts; ++seed) {
                TString resFile = Form("%s/result_seed_%d.root", logDirCentral.Data(), seed);
                TString logFile = Form("%s/log_seed_%d.txt", logDirCentral.Data(), seed);
                TFile* fRes = TFile::Open(resFile);
                if (fRes && !fRes->IsZombie()) {
                    TTree* tChi2 = fRes->Get<TTree>("chi2");
                    TTree* tMD = fRes->Get<TTree>("mass_dependent_params");
                    if (tChi2 && tMD) {
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

                        TObjArray* branches = tMD->GetListOfBranches();
                        for (int b = 0; b < branches->GetEntries(); ++b) {
                            TString bName = branches->At(b)->GetName();
                            if (!bName.EndsWith("_error")) {
                                double val = 0.0;
                                tMD->SetBranchAddress(bName, &val);
                                tMD->GetEntry(0);
                                rec.params[bName] = val;
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
            std::cout << "Running " << nRandomStarts << " randomized initial starts sequentially..." << std::endl;
            std::cout << "Fit logs are being saved to: " << logDirCentral << "/\n" << std::endl;

            for (int seed = 1; seed <= nRandomStarts; ++seed) {
                TString logFile = Form("%s/log_seed_%d.txt", logDirCentral.Data(), seed);

                // Redirect stdout and stderr to per-fit log file
                gSystem->RedirectOutput(logFile.Data(), "w");

                std::cout << "==================================================" << std::endl;
                std::cout << " Central Fit Run Log: Seed " << seed << std::endl;
                std::cout << "==================================================" << std::endl;

                MassDependentFitter::ParameterManager trialManager(fitter);
                trialManager.AddMassDependentParameter(std::vector<int>{2}, seed, 0);
                trialManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
                trialManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);

                trialManager.PrintParameters();

                fitter.SetMinimizerPrintLevel(1); // Full Minuit details in log file
                fitter.MinimizeChi2(trialManager);

                double chi2 = fitter.GetLastChi2();
                bool valid = fitter.IsMinimizerValid();
                int status = fitter.GetMinimizerStatus();

                std::cout << "\nFit Finished: Chi2 = " << chi2 
                          << " | Status = " << status 
                          << " | Valid = " << (valid ? "true" : "false") << std::endl;

                // Restore terminal output
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

        // Analysis of convergence and global minimum cluster
        std::vector<FitRecord> validFits;
        for (const auto& r : fitRecords) {
            if (r.valid) validFits.push_back(r);
        }
        std::sort(validFits.begin(), validFits.end(), [](const FitRecord& a, const FitRecord& b) {
            return a.chi2 < b.chi2;
        });

        // Threshold for identical minima (within 0.05 units of chi2)
        const double clusterThreshold = 0.05;
        std::vector<int> globalMinimaSeeds;
        for (const auto& r : validFits) {
            if (std::abs(r.chi2 - bestChi2) <= clusterThreshold) {
                globalMinimaSeeds.push_back(r.seed);
            }
        }

        std::cout << "\n========================================================" << std::endl;
        std::cout << "           RANDOMIZED FITS SUMMARY                      " << std::endl;
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

        // Save detailed text summary to file
        TString summaryPath = Form("%s/central_fits_summary.txt", logDir.Data());
        std::ofstream sumFile(summaryPath.Data());
        if (sumFile.is_open()) {
            sumFile << "========================================================" << std::endl;
            sumFile << "           RANDOMIZED FITS SUMMARY                      " << std::endl;
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
            sumFile << "  Rank |  Seed |       Chi2 | Delta Chi2 | Status | Interpretation" << std::endl;
            sumFile << " --------------------------------------------------------" << std::endl;
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
            TCanvas* cChi2 = new TCanvas("c_chi2_dist", "Randomized Starts Chi2 Distribution", 1200, 500);
            cChi2->Divide(2, 1);

            // Pad 1: Sorted Chi2 Plateau Plot (shows flat plateau for global minimum)
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
            grPlateau->SetName("gr_chi2_plateau");
            grPlateau->SetTitle("Sorted #chi^{2} vs Fit Rank (Global Minimum Plateau);Fit Rank (ordered by #chi^{2});#chi^{2}");
            grPlateau->SetMarkerStyle(20);
            grPlateau->SetMarkerSize(0.9);
            grPlateau->SetMarkerColor(kAzure + 1);
            grPlateau->SetLineColor(kAzure + 1);
            grPlateau->SetLineWidth(2);
            grPlateau->GetYaxis()->SetTitleOffset(1.3);
            grPlateau->Draw("ALP");

            // Dashed line at lowest chi2
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

            TH1D* hChi2 = new TH1D("h_chi2_dist", "Distribution of #chi^{2} Across Random Starts;#chi^{2};Number of Fits", 
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
            cChi2->SaveAs(Form("%s/central_chi2_distribution.png", logDir.Data()));
            cChi2->SaveAs(Form("%s/central_chi2_distribution.pdf", logDir.Data()));
            cChi2->SaveAs("central_chi2_distribution.png");
            cChi2->SaveAs("central_chi2_distribution.pdf");
            std::cout << "[Diagnostic] Diagnostic plot saved: " << logDir << "/central_chi2_distribution.png and .pdf\n" << std::endl;
        }
    }

    // Set up paramManager with the best solution
    MassDependentFitter::ParameterManager paramManager(fitter);
    paramManager.AddMassDependentParameter(std::vector<int>{2}, (bestSeed > 0 ? bestSeed : 0), 0);
    paramManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
    paramManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);

    if (!bestParams.empty()) {
        paramManager.SetParameterInitialValues(bestParams);
    }

    // Save final central fit log
    TString bestCentralLog = Form("%s/log_best_central_fit.txt", logDir.Data());
    gSystem->RedirectOutput(bestCentralLog.Data(), "w");
    fitter.SetMinimizerPrintLevel(2);
    fitter.MinimizeChi2(paramManager);
    gSystem->RedirectOutput(0);

    // Save best central parameters ROOT file for parallel MCMC workers
    TString centralBestFile = Form("%s/central_best_params.root", logDir.Data());
    fitter.MakeResultTree(paramManager, centralBestFile);

    std::cout << "Nominal central fit finished with Chi2 = " << fitter.GetLastChi2() 
              << " | Log: " << bestCentralLog 
              << " | Best parameters file: " << centralBestFile << std::endl;

    // 2. Propagate uncertainty across MCMC tree starting from the global minimum
    std::cout << "\n==========================================" << std::endl;
    std::cout << " Step 2: Propagating Uncertainty over MCMC Tree" << std::endl;
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
        std::cout << Form("Running MCMC propagation in parallel across %d workers (%lld entries)...",
                          nParallelJobs, actualEntriesToFit) << std::endl;
        std::cout << "Worker logs: " << logDirMCMC << "/log_worker_*.txt" << std::endl;
        std::cout << "Fit logs:    " << logDirMCMC << "/log_entry_*.txt\n" << std::endl;

        TString runScriptDir = gSystem->DirName(gSystem->Which(gROOT->GetMacroPath(), "Mom2Amp_DBreitWigner_MCMC.C"));
        if (runScriptDir.IsNull() || gSystem->AccessPathName(runScriptDir + "/RunParallelMCMCFits.sh") != 0) {
            if (gSystem->AccessPathName("macros/RunParallelMCMCFits.sh") == 0) {
                runScriptDir = "macros";
            } else {
                runScriptDir = ".";
            }
        }
        TString cmd = Form("bash %s/RunParallelMCMCFits.sh %d %lld %lld %d %s %s",
                           runScriptDir.Data(), nParallelJobs, burnIn, actualEntriesToFit, stride,
                           centralBestFile.Data(), logDirMCMC.Data());

        int ret = gSystem->Exec(cmd.Data());
        if (ret != 0) {
            std::cerr << "Warning: Parallel MCMC script returned exit code: " << ret << std::endl;
        }

        // Merge chunk ROOT files into outputFile
        std::cout << "\nMerging parallel worker chunks into " << outputFile << "..." << std::endl;
        TFile* fOut = TFile::Open(outputFile, "RECREATE");
        if (!fOut || fOut->IsZombie()) {
            std::cerr << "Error: Cannot open output file " << outputFile << " for merging!" << std::endl;
            if (fOut) delete fOut;
            return;
        }

        // Copy central trees from centralBestFile
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
            TTree* tMD = fCent->Get<TTree>("mass_dependent_params");
            if (tMD) {
                fOut->cd();
                TTree* cParams = tMD->CloneTree(-1, "fast");
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

        // Chain and merge MCMC trees from all worker chunks
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
            std::vector<double> allChi2;

            TString mcmcSumPath = Form("%s/mcmc_fits_summary.txt", logDir.Data());
            std::ofstream mSumFile(mcmcSumPath.Data());
            if (mSumFile.is_open()) {
                mSumFile << "========================================================" << std::endl;
                mSumFile << "         MCMC PROPAGATION FITS SUMMARY                  " << std::endl;
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
                    allChi2.push_back(c_val);
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
            std::cout << "         MCMC PROPAGATION FITS SUMMARY                  " << std::endl;
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
        TCanvas* cChi2 = dynamic_cast<TCanvas*>(gROOT->FindObject("c_chi2_dist"));
        if (cChi2) {
            fOut->cd();
            cChi2->Write("c_chi2_dist", TObject::kOverwrite);
            std::cout << "[Diagnostic] Diagnostic canvas written into " << outputFile << std::endl;
        }
        fOut->Close();
        delete fOut;
    } else {
        // Sequential fallback
        fitter.FitMCMCTree(paramManager, massDepMoments, outputFile, burnIn, maxEntries, stride, logDirMCMC);

        TFile* fOut = TFile::Open(outputFile, "UPDATE");
        if (fOut && !fOut->IsZombie()) {
            TCanvas* cChi2 = dynamic_cast<TCanvas*>(gROOT->FindObject("c_chi2_dist"));
            if (cChi2) {
                cChi2->Write("c_chi2_dist", TObject::kOverwrite);
                std::cout << "[Diagnostic] Diagnostic canvas written into " << outputFile << std::endl;
            }
            fOut->Close();
            delete fOut;
        }
    }
}


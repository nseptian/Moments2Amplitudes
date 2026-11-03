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
#include <iostream>

const Double_t FIRST_MASS_BIN_CENTER = 1.06; // in GeV
const Double_t MASS_BIN_WIDTH = 0.04;
const Int_t N_MASS_BINS = 14;

TString FIT_RESULTS_DIR = "/d/grid17/septian/EtaPi0Analysis/run/MomentFit_m080200/t010020/";

/**
 * @brief D-wave Breit-Wigner fit to H4 moments with randomized starts
 * @param nRandomStarts Number of random seeds to search for global minimum (default: 100). If <= 1, runs single seed.
 * @param outputFile Best result ROOT file (default: "result_tree_DBreitWigner_best.root")
 * @param logDir Directory to store logs and diagnostic plots (default: "logs_dbreitwigner")
 * @param nParallelJobs Number of parallel worker jobs (default: 40)
 */
void Mom2Amp_DBreitWigner(int nRandomStarts = 100,
                          TString outputFile = "result_tree_DBreitWigner_best.root",
                          TString logDir = "logs_dbreitwigner",
                          int nParallelJobs = 40) {

    auto& setup = ConfigureAmpsNoValues(2, 2, 2); // (Lmax, MMax, Nref)
    
    // Load moments data with error handling
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
        std::cerr << "Please check that directory exists: " << FIT_RESULTS_DIR << std::endl;
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
    fitter.SetEquationValues(massDepMoments);
    fitter.PrintIncludedMoments();

    gSystem->mkdir(logDir.Data(), true);
    fitter.SetMinimizerTolerance(0.01);

    if (nRandomStarts > 1) {
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

        if (nParallelJobs > 1) {
            std::cout << Form("Running %d randomized D-wave fits in parallel across %d workers...",
                              nRandomStarts, nParallelJobs) << std::endl;
            std::cout << "Fit logs are being saved to: " << logDir << "/\n" << std::endl;

            TString runScriptDir = gSystem->DirName(gSystem->Which(gROOT->GetMacroPath(), "Mom2Amp_DBreitWigner.C"));
            if (runScriptDir.IsNull() || gSystem->AccessPathName(runScriptDir + "/RunParallelRandomizedFits.sh") != 0) {
                if (gSystem->AccessPathName("macros/RunParallelRandomizedFits.sh") == 0) {
                    runScriptDir = "macros";
                } else {
                    runScriptDir = ".";
                }
            }
            TString cmd = Form("bash %s/RunParallelRandomizedFits.sh %d %d %s",
                               runScriptDir.Data(), nRandomStarts, nParallelJobs, logDir.Data());

            int ret = gSystem->Exec(cmd.Data());
            if (ret != 0) {
                std::cerr << "Warning: Parallel fit script returned non-zero exit code: " << ret << std::endl;
            }

            for (int seed = 1; seed <= nRandomStarts; ++seed) {
                TString resFile = Form("%s/result_seed_%d.root", logDir.Data(), seed);
                TString logFile = Form("%s/log_seed_%d.txt", logDir.Data(), seed);
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
            // Sequential
            for (int seed = 1; seed <= nRandomStarts; ++seed) {
                TString logFile = Form("%s/log_seed_%d.txt", logDir.Data(), seed);
                gSystem->RedirectOutput(logFile.Data(), "w");

                MassDependentFitter::ParameterManager trialManager(fitter);
                trialManager.AddMassDependentParameter(std::vector<int>{2}, seed, 0);
                trialManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
                trialManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);

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

        // Summary analysis
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
        std::cout << "           D-BREIT-WIGNER FITS SUMMARY                  " << std::endl;
        std::cout << "========================================================" << std::endl;
        std::cout << Form(" Total starts:     %d", nRandomStarts) << std::endl;
        std::cout << Form(" Valid converged:  %zu (%.1f%%)", validFits.size(), (validFits.size() * 100.0 / nRandomStarts)) << std::endl;
        std::cout << Form(" Lowest Chi2:      %.2f (Seed %d)", bestChi2, bestSeed) << std::endl;
        std::cout << Form(" Global Minima Cluster: %zu fits within Delta Chi2 <= %.2f",
                          globalMinimaSeeds.size(), clusterThreshold) << std::endl;
        std::cout << "========================================================\n" << std::endl;

        // Save detailed text summary
        TString summaryPath = Form("%s/dbreitwigner_fits_summary.txt", logDir.Data());
        std::ofstream sumFile(summaryPath.Data());
        if (sumFile.is_open()) {
            sumFile << "========================================================" << std::endl;
            sumFile << "           D-BREIT-WIGNER FITS SUMMARY                  " << std::endl;
            sumFile << "========================================================" << std::endl;
            sumFile << " Total starts:     " << nRandomStarts << std::endl;
            sumFile << " Valid converged:  " << validFits.size() << std::endl;
            sumFile << " Lowest Chi2:      " << std::fixed << std::setprecision(2) << bestChi2 << " (Seed " << bestSeed << ")" << std::endl;
            sumFile << " Global Minima Cluster: " << globalMinimaSeeds.size() << " fits" << std::endl;
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
        }

        // Draw Chi2 Distribution and Plateau Plot
        if (!validFits.empty()) {
            gStyle->SetOptStat(1111);
            TCanvas* cChi2 = new TCanvas("c_dbreitwigner_chi2_dist", "D-BreitWigner Starts Chi2 Distribution", 1200, 500);
            cChi2->Divide(2, 1);

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
            grPlateau->SetName("gr_dbreitwigner_chi2_plateau");
            grPlateau->SetTitle("Sorted #chi^{2} vs Fit Rank (D-BreitWigner Plateau);Fit Rank (ordered by #chi^{2});#chi^{2}");
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

            cChi2->cd(2);
            gPad->SetGrid();
            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.06);

            double histMax = validFits.back().chi2;
            if (nValid >= 5 && histMax > validFits[std::min((size_t)nValid - 1, (size_t)25)].chi2 * 1.5) {
                histMax = validFits[std::min((size_t)nValid - 1, (size_t)25)].chi2 * 1.5;
            }
            if (histMax <= bestChi2) histMax = bestChi2 + 50.0;

            TH1D* hChi2 = new TH1D("h_dbreitwigner_chi2_dist", "Distribution of #chi^{2} Across D-BreitWigner Starts;#chi^{2};Number of Fits", 
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
            cChi2->SaveAs(Form("%s/dbreitwigner_chi2_distribution.png", logDir.Data()));
            cChi2->SaveAs(Form("%s/dbreitwigner_chi2_distribution.pdf", logDir.Data()));
            cChi2->SaveAs("dbreitwigner_chi2_distribution.png");
            cChi2->SaveAs("dbreitwigner_chi2_distribution.pdf");
        }

        // Final best fit
        MassDependentFitter::ParameterManager bestManager(fitter);
        bestManager.AddMassDependentParameter(std::vector<int>{2}, (bestSeed > 0 ? bestSeed : 0), 0);
        bestManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
        bestManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);
        if (!bestParams.empty()) {
            bestManager.SetParameterInitialValues(bestParams);
        }

        TString bestLog = Form("%s/log_best_dbreitwigner_fit.txt", logDir.Data());
        gSystem->RedirectOutput(bestLog.Data(), "w");
        fitter.SetMinimizerPrintLevel(2);
        fitter.MinimizeChi2(bestManager);
        gSystem->RedirectOutput(0);

        fitter.MakeResultTree(bestManager, outputFile);
        fitter.MakeResultTree(bestManager, "result_tree_test_0.root");
        std::cout << "Best D-wave fit saved to " << outputFile 
                  << " (Chi2 = " << fitter.GetLastChi2() << ", Log: " << bestLog << ")" << std::endl;

        TFile* fOut = TFile::Open(outputFile, "UPDATE");
        if (fOut && !fOut->IsZombie()) {
            TCanvas* cChi2 = dynamic_cast<TCanvas*>(gROOT->FindObject("c_dbreitwigner_chi2_dist"));
            if (cChi2) {
                cChi2->Write("c_chi2_dist", TObject::kOverwrite);
            }
            fOut->Close();
            delete fOut;
        }

    } else {
        // Single seed run
        int seed = (nRandomStarts <= 0 ? 0 : nRandomStarts);
        MassDependentFitter::ParameterManager paramManager(fitter);
        paramManager.AddMassDependentParameter(std::vector<int>{2}, seed, 0);
        paramManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
        paramManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);

        paramManager.PrintParameters();
        fitter.SetMinimizerPrintLevel(2);
        fitter.MinimizeChi2(paramManager);

        fitter.MakeResultTree(paramManager, outputFile);
        fitter.MakeResultTree(paramManager, "result_tree_test_0.root");
        std::cout << "Single fit saved to " << outputFile << " (Chi2 = " << fitter.GetLastChi2() << ")" << std::endl;
    }
}

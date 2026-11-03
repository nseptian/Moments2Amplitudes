#ifndef MCMC_FITTER_CPP
#define MCMC_FITTER_CPP

#include "MassDependentFitter.h"
#include "FitInternalUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include <iostream>
#include <vector>
#include <map>
#include <set>

namespace m2pw {

    void MassDependentFitter::SetEquationValuesFromMCMCEntry(const MassDependentMoments& moms, Long64_t entryIndex) {
        for (auto& [mass_bin, equations] : eqns_) {
            const MomentHelper moments = moms.GetMoments(mass_bin);
            for (auto& eqn : equations) {
                double val = moments.GetMCMCUnnormalizedVal(eqn.GetName(), entryIndex);
                eqn.SetTargetValue(val);
            }
        }
        cacheValid_ = false;
    }

    void MassDependentFitter::FitMCMCTree(ParameterManager& paramManager, 
                                          MassDependentMoments& moms, 
                                          const TString& outputFileName,
                                          Long64_t burnIn, 
                                          Long64_t maxEntries, 
                                          int stride,
                                          const TString& logDir) {
        if (eqns_.empty() || pars_.empty()) {
            LogFitWarn("MCMC_FIT", "reason=no_equations_or_parameters_configured");
            return;
        }

        if (!moms.HasMCMCSamples()) {
            std::cerr << "Error in FitMCMCTree: No MCMC samples found in MassDependentMoments." << std::endl;
            return;
        }

        const Long64_t totalAvailable = moms.GetMCMCEntries();
        if (burnIn < 0) burnIn = 0;
        if (stride < 1) stride = 1;

        const Long64_t startEntry = burnIn;
        Long64_t endEntry = totalAvailable;
        if (maxEntries > 0 && (startEntry + maxEntries) < totalAvailable) {
            endEntry = startEntry + maxEntries;
        }

        if (startEntry >= endEntry) {
            std::cerr << "Error in FitMCMCTree: startEntry (" << startEntry 
                      << ") >= endEntry (" << endEntry << "). No entries to fit." << std::endl;
            return;
        }

        const Long64_t totalFits = (endEntry - startEntry + stride - 1) / stride;

        std::cout << "\n========================================================" << std::endl;
        std::cout << " Starting MCMC Tree Uncertainty Propagation" << std::endl;
        std::cout << " Total available MCMC entries: " << totalAvailable << std::endl;
        std::cout << " Burn-in (start entry):        " << startEntry << std::endl;
        std::cout << " End entry:                    " << endEntry << std::endl;
        std::cout << " Stride (thinning):            " << stride << std::endl;
        std::cout << " Total fits to perform:        " << totalFits << std::endl;
        std::cout << " Output ROOT file:             " << outputFileName << std::endl;
        std::cout << "========================================================\n" << std::endl;

        // 1. Save central parameter baseline
        const std::map<TString, double> centralParamValues = paramManager.parsList;
        const std::vector<double> centralParErrors = paramManager.parErrors;
        const double centralChi2 = lastChi2_;
        const int centralStatus = minimizerStatus_;
        const bool centralIsValid = minimizerIsValid_;
        const bool centralIsValidError = minimizerHasValidErrors_;
        const int centralSeed = paramManager.randomSeed;

        const int nPars = pars_.begin()->second.Nvars();
        std::vector<TString> ampNames(nPars);
        for (int i = 0; i < nPars; ++i) {
            ampNames[i] = pars_.begin()->second.GetParName(i);
        }

        std::map<double, std::map<TString, double>> centralAmplitudes;
        for (const auto& [mass_bin, pars] : pars_) {
            for (const auto& ampName : ampNames) {
                centralAmplitudes[mass_bin][ampName] = pars.GetCurrentVal(ampName);
            }
        }

        std::set<TString> allMomentNames;
        for (const auto& [mass_bin, equations] : eqns_) {
            for (const auto& eqn : equations) {
                const TString& eqnName = eqn.GetName();
                if (eqnName.BeginsWith("H_")) {
                    allMomentNames.insert(eqnName);
                }
            }
        }
        std::vector<TString> momentNamesList(allMomentNames.begin(), allMomentNames.end());

        std::map<double, std::map<TString, double>> centralMoments;
        for (const auto& [mass_bin, equations] : eqns_) {
            for (const auto& eqn : equations) {
                if (eqn.GetName().BeginsWith("H_")) {
                    centralMoments[mass_bin][eqn.GetName()] = eqn.GetOrigFormulaValue();
                }
            }
        }

        // 2. Open output file
        TFile* outFile = TFile::Open(outputFileName, "RECREATE");
        if (!outFile || outFile->IsZombie()) {
            std::cerr << "Error in FitMCMCTree: Failed to create output file " << outputFileName << std::endl;
            if (outFile) delete outFile;
            return;
        }

        const std::vector<TString>& registeredParNames = paramManager.parIndexNames;

        // Central Trees
        TTree* centralParamsTree = new TTree("central_params", "Nominal Best-fit Parameters");
        std::vector<double> centralAllParVals(registeredParNames.size(), 0.0);
        std::vector<double> centralAllParErrs(registeredParNames.size(), 0.0);
        for (size_t i = 0; i < registeredParNames.size(); ++i) {
            const TString& pName = registeredParNames[i];
            auto valIt = centralParamValues.find(pName);
            centralAllParVals[i] = (valIt != centralParamValues.end()) ? valIt->second : 0.0;
            double err = 0.0;
            auto it = paramManager.nameToIndex.find(pName);
            if (it != paramManager.nameToIndex.end() && it->second < static_cast<int>(centralParErrors.size())) {
                err = centralParErrors[it->second];
            }
            centralAllParErrs[i] = err;
            centralParamsTree->Branch(pName, &centralAllParVals[i]);
            centralParamsTree->Branch(pName + "_error", &centralAllParErrs[i]);
        }
        centralParamsTree->Fill();

        TTree* centralChi2Tree = new TTree("central_chi2", "Nominal Fit Chi2 and Status");
        double c_chi2 = centralChi2;
        int c_status = centralStatus;
        bool c_isValid = centralIsValid;
        bool c_isValidError = centralIsValidError;
        int c_seed = centralSeed;
        centralChi2Tree->Branch("chi2", &c_chi2);
        centralChi2Tree->Branch("Status", &c_status);
        centralChi2Tree->Branch("isValid", &c_isValid);
        centralChi2Tree->Branch("isValidError", &c_isValidError);
        centralChi2Tree->Branch("seed", &c_seed);
        centralChi2Tree->Fill();

        TTree* centralResultTree = new TTree("central_result", "Nominal Amplitudes and Moments per Mass Bin");
        double c_mass_bin = 0.0;
        std::vector<double> c_ampVals(ampNames.size(), 0.0);
        std::vector<double> c_momVals(momentNamesList.size(), 0.0);
        centralResultTree->Branch("mass_bin", &c_mass_bin);
        for (size_t a = 0; a < ampNames.size(); ++a) {
            centralResultTree->Branch(ampNames[a], &c_ampVals[a]);
        }
        for (size_t m = 0; m < momentNamesList.size(); ++m) {
            centralResultTree->Branch(momentNamesList[m], &c_momVals[m]);
        }
        for (const auto& mass_bin : massBins_) {
            c_mass_bin = mass_bin;
            for (size_t a = 0; a < ampNames.size(); ++a) {
                c_ampVals[a] = centralAmplitudes[mass_bin][ampNames[a]];
            }
            for (size_t m = 0; m < momentNamesList.size(); ++m) {
                c_momVals[m] = centralMoments[mass_bin][momentNamesList[m]];
            }
            centralResultTree->Fill();
        }

        // MCMC Trees
        TTree* mcmcParamsTree = new TTree("mcmc_params", "MCMC Posterior Sample Fitted Parameters");
        Long64_t mcmc_entry = 0;
        double m_chi2 = 0.0;
        int m_status = -1;
        bool m_isValid = false;
        bool m_isValidError = false;
        mcmcParamsTree->Branch("mcmc_entry", &mcmc_entry);
        std::vector<double> mcmcParVals(registeredParNames.size(), 0.0);
        for (size_t i = 0; i < registeredParNames.size(); ++i) {
            mcmcParamsTree->Branch(registeredParNames[i], &mcmcParVals[i]);
        }
        mcmcParamsTree->Branch("chi2", &m_chi2);
        mcmcParamsTree->Branch("Status", &m_status);
        mcmcParamsTree->Branch("isValid", &m_isValid);
        mcmcParamsTree->Branch("isValidError", &m_isValidError);

        TTree* mcmcChi2Tree = new TTree("mcmc_chi2", "MCMC Posterior Sample Chi2 and Status");
        mcmcChi2Tree->Branch("mcmc_entry", &mcmc_entry);
        mcmcChi2Tree->Branch("chi2", &m_chi2);
        mcmcChi2Tree->Branch("Status", &m_status);
        mcmcChi2Tree->Branch("isValid", &m_isValid);
        mcmcChi2Tree->Branch("isValidError", &m_isValidError);

        TTree* mcmcResultTree = new TTree("mcmc_result", "MCMC Amplitudes and Moments per Mass Bin");
        double m_mass_bin = 0.0;
        std::vector<double> m_ampVals(ampNames.size(), 0.0);
        std::vector<double> m_momVals(momentNamesList.size(), 0.0);
        mcmcResultTree->Branch("mcmc_entry", &mcmc_entry);
        mcmcResultTree->Branch("mass_bin", &m_mass_bin);
        for (size_t a = 0; a < ampNames.size(); ++a) {
            mcmcResultTree->Branch(ampNames[a], &m_ampVals[a]);
        }
        for (size_t m = 0; m < momentNamesList.size(); ++m) {
            mcmcResultTree->Branch(momentNamesList[m], &m_momVals[m]);
        }

        // 3. MCMC Loop
        int originalPrintLevel = minimizerPrintLevel_;
        SetMinimizerPrintLevel(0);

        TStopwatch sw;
        sw.Start();

        Long64_t nProcessed = 0;
        Long64_t nConverged = 0;

        for (Long64_t entry = startEntry; entry < endEntry; entry += stride) {
            if (!logDir.IsNull()) {
                gSystem->mkdir(logDir.Data(), true);
                TString entryLog = Form("%s/log_entry_%lld.txt", logDir.Data(), entry);
                gSystem->RedirectOutput(entryLog.Data(), "w");
                SetMinimizerPrintLevel(1);
                std::cout << "==================================================" << std::endl;
                std::cout << " MCMC Entry Fit Run Log: Entry " << entry << std::endl;
                std::cout << "==================================================" << std::endl;
            }

            SetEquationValuesFromMCMCEntry(moms, entry);
            paramManager.SetParameterInitialValues(centralParamValues);
            MinimizeChi2(paramManager);

            if (!logDir.IsNull()) {
                std::cout << "\nFit Finished: Chi2 = " << lastChi2_
                          << " | Status = " << minimizerStatus_
                          << " | Valid = " << (minimizerIsValid_ ? "true" : "false") << std::endl;
                gSystem->RedirectOutput(0);
                SetMinimizerPrintLevel(0);
            }

            Chi2Function finalEvalFunction(*this, paramManager.parIndexNames, paramManager);
            std::vector<double> finalParVals = paramManager.GetValues();
            finalEvalFunction(finalParVals.data());

            mcmc_entry = entry;
            m_chi2 = lastChi2_;
            m_status = minimizerStatus_;
            m_isValid = minimizerIsValid_;
            m_isValidError = minimizerHasValidErrors_;

            if (m_isValid) nConverged++;

            for (size_t i = 0; i < registeredParNames.size(); ++i) {
                mcmcParVals[i] = paramManager.parsList[registeredParNames[i]];
            }
            mcmcParamsTree->Fill();
            mcmcChi2Tree->Fill();

            for (const auto& mass_bin : massBins_) {
                m_mass_bin = mass_bin;
                const auto parsIt = pars_.find(mass_bin);
                if (parsIt != pars_.end()) {
                    for (size_t a = 0; a < ampNames.size(); ++a) {
                        m_ampVals[a] = parsIt->second.GetCurrentVal(ampNames[a]);
                    }
                }
                const auto eqnIt = eqns_.find(mass_bin);
                for (size_t m = 0; m < momentNamesList.size(); ++m) {
                    m_momVals[m] = 0.0;
                    if (eqnIt != eqns_.end()) {
                        for (const auto& eqn : eqnIt->second) {
                            if (eqn.GetName() == momentNamesList[m]) {
                                m_momVals[m] = eqn.GetOrigFormulaValue();
                                break;
                            }
                        }
                    }
                }
                mcmcResultTree->Fill();
            }

            nProcessed++;

            if (nProcessed % 50 == 0 || nProcessed == totalFits) {
                double elapsed = sw.RealTime();
                sw.Continue();
                double rate = nProcessed / (elapsed > 0.0 ? elapsed : 1.0);
                double eta = (totalFits - nProcessed) / (rate > 0.0 ? rate : 1.0);
                std::cout << Form("[FitMCMCTree] Sample %lld/%lld (%lld%%) | Chi2=%.2f Status=%d Valid=%d | Rate: %.1f fits/s | ETA: %.1fs",
                                  nProcessed, totalFits, (nProcessed * 100 / totalFits),
                                  m_chi2, m_status, (int)m_isValid, rate, eta)
                          << std::endl;
            }
        }

        // 4. Restore state
        SetMinimizerPrintLevel(originalPrintLevel);
        SetEquationValues(moms);
        paramManager.SetParameterInitialValues(centralParamValues);

        // 5. Write and close output file
        outFile->Write();
        outFile->Close();
        delete outFile;

        std::cout << "\n========================================================" << std::endl;
        std::cout << Form(" FitMCMCTree Finished: %lld/%lld (%lld%%) converged.", nConverged, totalFits, (totalFits > 0 ? nConverged * 100 / totalFits : 0)) << std::endl;
        std::cout << " Results saved to: " << outputFileName << std::endl;
        std::cout << "========================================================\n" << std::endl;
    }

} // namespace m2pw

#endif // MCMC_FITTER_CPP


#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include <iostream>
#include <map>

/**
 * @brief Worker macro to process a chunk of MCMC entries in parallel
 * @param workerId Unique integer ID for this worker
 * @param startEntry Global MCMC start entry for this chunk
 * @param nEntries Number of entries this worker should process
 * @param stride Step size between entries
 * @param centralBestFile ROOT file containing nominal best-fit parameters (with tree mass_dependent_params)
 * @param outDir Directory to write chunk ROOT files and per-entry log files
 */
void RunSingleMCMCWorker(int workerId,
                         Long64_t startEntry,
                         Long64_t nEntries,
                         int stride,
                         TString centralBestFile,
                         TString outDir = "logs/mcmc_fits") {

    // 1. Ensure output directory exists
    gSystem->mkdir(outDir.Data(), true);

    std::cout << Form("========================================================\n");
    std::cout << Form(" [Worker %2d] Starting chunk: startEntry=%lld, nEntries=%lld, stride=%d\n",
                      workerId, startEntry, nEntries, stride);
    std::cout << Form(" Central best parameters: %s\n", centralBestFile.Data());
    std::cout << Form(" Output directory:        %s\n", outDir.Data());
    std::cout << Form("========================================================\n");

    // 2. Configure amplitudes and moments
    auto& setup = ConfigureAmpsNoValues(2, 2, 2);
    MassDependentMoments massDepMoments(14, 0.04, 1.06,
                                       "/d/grid17/septian/EtaPi0Analysis/run/MomentFit_m080200/t010020/",
                                       "ResultsBruMcmcCovariance.root", "Mpi0eta");

    auto a2_1320_config = m2pw::MassDependentFitter::CreateCustomConfig(
        {{2, {m2pw::MassDependentFitter::MassDependenceConfig::WaveModelConfig{"a2_1320", 
             m2pw::MassDependentFitter::MassDependenceConfig::ModelType::BreitWignerPlusPolynomial, -1}}}}, 
        {0}
    );
    auto H4_config = m2pw::MassDependentFitter::CreateL4OnlyConfig();
    
    m2pw::MassDependentFitter fitter{setup, massDepMoments.GetMassBins(), 2, 0.0, {"H_3"}, a2_1320_config, H4_config};
    fitter.SetEquationValues(massDepMoments);

    // 3. Set up ParameterManager
    m2pw::MassDependentFitter::ParameterManager paramManager(fitter);
    paramManager.AddMassDependentParameter(std::vector<int>{2}, 0, 0);
    paramManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
    paramManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);

    // 4. Load central best-fit parameters from centralBestFile
    std::map<TString, double> centralParamValues;
    TFile* fCentral = TFile::Open(centralBestFile);
    if (!fCentral || fCentral->IsZombie()) {
        std::cerr << "Error in RunSingleMCMCWorker: Cannot open central file " << centralBestFile << std::endl;
        if (fCentral) delete fCentral;
        return;
    }
    TTree* tMD = fCentral->Get<TTree>("mass_dependent_params");
    if (tMD) {
        TObjArray* branches = tMD->GetListOfBranches();
        for (int b = 0; b < branches->GetEntries(); ++b) {
            TString bName = branches->At(b)->GetName();
            if (!bName.EndsWith("_error")) {
                double val = 0.0;
                tMD->SetBranchAddress(bName, &val);
                tMD->GetEntry(0);
                centralParamValues[bName] = val;
            }
        }
    } else {
        std::cerr << "Error in RunSingleMCMCWorker: Tree mass_dependent_params not found in " << centralBestFile << std::endl;
        fCentral->Close();
        delete fCentral;
        return;
    }
    fCentral->Close();
    delete fCentral;

    paramManager.SetParameterInitialValues(centralParamValues);
    fitter.SetMinimizerTolerance(0.01);

    // 5. Run MCMC propagation for this chunk using FitMCMCTree
    TString chunkFileName = Form("%s/chunk_worker_%d.root", outDir.Data(), workerId);
    fitter.FitMCMCTree(paramManager, massDepMoments, chunkFileName, startEntry, nEntries, stride, outDir);

    std::cout << Form("[Worker %2d] Finished chunk. Output saved to: %s\n", workerId, chunkFileName.Data());
}

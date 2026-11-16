#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include <iostream>
#include <map>

/**
 * @brief Worker macro to process a chunk of MCMC entries for S-Piecewise + D-BreitWigner
 * @param workerId Unique integer ID for this worker
 * @param startEntry Global MCMC start entry for this chunk
 * @param nEntries Number of entries this worker should process
 * @param stride Step size between entries
 * @param centralBestFile ROOT file containing nominal best-fit parameters
 * @param outDir Directory to write chunk ROOT files and per-entry log files
 */
void RunSingleMCMCWorkerSPiecewise(int workerId,
                                   Long64_t startEntry,
                                   Long64_t nEntries,
                                   int stride,
                                   TString centralBestFile,
                                   TString outDir = "logs_spiecewise/mcmc_fits") {

    // 1. Ensure output directory exists
    gSystem->mkdir(outDir.Data(), true);

    std::cout << Form("========================================================\n");
    std::cout << Form(" [Worker %2d] Starting S-Piecewise MCMC chunk: startEntry=%lld, nEntries=%lld, stride=%d\n",
                      workerId, startEntry, nEntries, stride);
    std::cout << Form(" Central best parameters: %s\n", centralBestFile.Data());
    std::cout << Form(" Output directory:        %s\n", outDir.Data());
    std::cout << Form("========================================================\n");

    const Double_t FIRST_MASS_BIN_CENTER = 0.82; // in GeV
    const Double_t MASS_BIN_WIDTH = 0.04;
    const Int_t N_MASS_BINS = 20;
    TString FIT_RESULTS_DIR = "/d/grid17/septian/EtaPi0Analysis/run/MomentFit_m080200/t010020/";

    // 2. Configure amplitudes and moments
    auto& setup = ConfigureAmpsNoValues(2, 2, 2);
    MassDependentMoments massDepMoments(N_MASS_BINS, MASS_BIN_WIDTH, FIRST_MASS_BIN_CENTER,
                                       FIT_RESULTS_DIR, "ResultsBruMcmcCovariance.root", "Mpi0eta");

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

    // 3. Set up ParameterManager
    MassDependentFitter::ParameterManager paramManager(fitter);

    paramManager.AddMassIndependentParameter(std::vector<int>{0, 1}, 0);

    for (const double massBin : massBins) {
        for (const TString& parName : MassDependentFitter::GetParNames()) {
            if (parName.Contains("_1_")) {
                TString fixedName = Form("MI_%1.6f_%s", massBin, parName.Data());
                paramManager.SetInitialValue(fixedName, 0.0, true);
            }
        }
    }

    paramManager.AddMassDependentParameter(std::vector<int>{2}, 0, 0);

    // 4. Load central best-fit parameters from centralBestFile
    std::map<TString, double> centralParamValues;
    TFile* fCentral = TFile::Open(centralBestFile);
    if (!fCentral || fCentral->IsZombie()) {
        std::cerr << "Error in RunSingleMCMCWorkerSPiecewise: Cannot open central file " << centralBestFile << std::endl;
        if (fCentral) delete fCentral;
        return;
    }

    TTree* tAll = fCentral->Get<TTree>("all_params");
    if (!tAll) tAll = fCentral->Get<TTree>("central_params");
    if (tAll) {
        TObjArray* branches = tAll->GetListOfBranches();
        for (int b = 0; b < branches->GetEntries(); ++b) {
            TString bName = branches->At(b)->GetName();
            if (!bName.EndsWith("_error")) {
                double val = 0.0;
                tAll->SetBranchAddress(bName, &val);
                tAll->GetEntry(0);
                centralParamValues[bName] = val;
            }
        }
    } else {
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
        }
    }
    fCentral->Close();
    delete fCentral;

    paramManager.SetParameterInitialValues(centralParamValues);

    // Fix D-wave mass and width to central best-fit values
    if (centralParamValues.find("MD_shared_a2_1320_Mass") != centralParamValues.end()) {
        paramManager.SetInitialValue("MD_shared_a2_1320_Mass", centralParamValues["MD_shared_a2_1320_Mass"], true);
    } else {
        paramManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
    }
    if (centralParamValues.find("MD_shared_a2_1320_Width") != centralParamValues.end()) {
        paramManager.SetInitialValue("MD_shared_a2_1320_Width", centralParamValues["MD_shared_a2_1320_Width"], true);
    } else {
        paramManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);
    }

    fitter.SetMinimizerTolerance(0.01);

    // 5. Run MCMC propagation for this chunk using FitMCMCTree
    TString chunkFileName = Form("%s/chunk_worker_%d.root", outDir.Data(), workerId);
    fitter.FitMCMCTree(paramManager, massDepMoments, chunkFileName, startEntry, nEntries, stride, outDir);

    std::cout << Form("[Worker %2d] Finished chunk. Output saved to: %s\n", workerId, chunkFileName.Data());
}


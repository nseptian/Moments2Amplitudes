#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include <TRandom3.h>
#include <iostream>

TString ResolveDWaveFileInWorker(TString userPath) {
    if (!userPath.IsNull() && gSystem->AccessPathName(userPath) == 0) {
        return userPath;
    }
    std::vector<TString> candidates = {
        userPath,
        "result_tree_DBreitWigner_best.root",
        "macros/result_tree_DBreitWigner_best.root",
        "result_tree_test_0.root",
        "macros/result_tree_test_0.root",
        "logs/central_best_params.root",
        "macros/logs/central_best_params.root"
    };
    for (const auto& path : candidates) {
        if (!path.IsNull() && gSystem->AccessPathName(path) == 0) {
            return path;
        }
    }
    return "";
}

void RunSingleSeedSPiecewise(int seed, 
                             TString dWaveFile = "result_tree_DBreitWigner_best.root", 
                             TString outDir = "logs_spiecewise") {
    // 1. Ensure output directory exists
    gSystem->mkdir(outDir.Data(), true);

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

    // Add S-wave (ell=0) and P-wave (ell=1)
    paramManager.AddMassIndependentParameter(std::vector<int>{0, 1}, seed);

    // Fix P-wave parameters to 0
    for (const double massBin : massBins) {
        for (const TString& parName : MassDependentFitter::GetParNames()) {
            if (parName.Contains("_1_")) {
                TString fixedName = Form("MI_%1.6f_%s", massBin, parName.Data());
                paramManager.SetInitialValue(fixedName, 0.0, true);
            }
        }
    }

    // Initialize ONLY the D-waves from dWaveFile
    TString resolvedFile = ResolveDWaveFileInWorker(dWaveFile);
    if (!resolvedFile.IsNull()) {
        paramManager.AddMassDependentParameter(std::vector<int>{2}, resolvedFile, 0);
        double dMass = 1.318;
        double dWidth = 0.107;
        if (paramManager.parsList.find("MD_shared_a2_1320_Mass") != paramManager.parsList.end()) {
            dMass = paramManager.parsList["MD_shared_a2_1320_Mass"];
        }
        if (paramManager.parsList.find("MD_shared_a2_1320_Width") != paramManager.parsList.end()) {
            dWidth = paramManager.parsList["MD_shared_a2_1320_Width"];
        }
        paramManager.SetInitialValue("MD_shared_a2_1320_Mass", dMass, true);
        paramManager.SetInitialValue("MD_shared_a2_1320_Width", dWidth, true);
    } else {
        paramManager.AddMassDependentParameter(std::vector<int>{2}, 0, 0);
        paramManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
        paramManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);
    }

    // For the piecewise S-waves: initialize magnitude randomly (positive) and phases to zero
    TRandom3 rng(seed > 0 ? seed : 12345);
    for (const double massBin : massBins) {
        for (const TString& parName : MassDependentFitter::GetParNames()) {
            if (parName.Contains("_0_0")) {
                TString sName = Form("MI_%1.6f_%s", massBin, parName.Data());
                if (paramManager.parsList.find(sName) != paramManager.parsList.end()) {
                    if (sName.Contains("phi")) {
                        paramManager.SetInitialValue(sName, 0.0, false);
                    } else {
                        double rMag = rng.Uniform(0.0, 200.0);
                        paramManager.SetInitialValue(sName, rMag, false);
                    }
                }
            }
        }
    }

    std::cout << "==================================================" << std::endl;
    std::cout << " S-Piecewise D-BreitWigner Fit Run Log: Seed " << seed << std::endl;
    std::cout << " D-wave input file: " << resolvedFile << std::endl;
    std::cout << "==================================================" << std::endl;
    paramManager.PrintParameters();

    // 4. Run minimization
    fitter.SetMinimizerTolerance(0.01);
    fitter.SetMinimizerPrintLevel(1);
    fitter.MinimizeChi2(paramManager);

    double chi2 = fitter.GetLastChi2();
    int status = fitter.GetMinimizerStatus();
    bool valid = fitter.IsMinimizerValid();

    // 5. Save result tree for this seed
    TString outRoot = Form("%s/result_seed_%d.root", outDir.Data(), seed);
    fitter.MakeResultTree(paramManager, outRoot);

    std::cout << Form("\n[RunSingleSeedSPiecewise] Seed %d Finished: Chi2 = %.2f | Status = %d | Valid = %d | Output = %s\n",
                      seed, chi2, status, (int)valid, outRoot.Data()) << std::endl;
}

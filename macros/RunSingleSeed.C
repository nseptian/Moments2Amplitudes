#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include <iostream>

void RunSingleSeed(int seed, TString outDir = "logs/central_fits") {
    // 1. Ensure output directory exists
    gSystem->mkdir(outDir.Data(), true);

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

    // 3. Initialize parameters with randomized seed
    m2pw::MassDependentFitter::ParameterManager pm(fitter);
    pm.AddMassDependentParameter(std::vector<int>{2}, seed, 0);
    pm.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
    pm.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);

    std::cout << "==================================================" << std::endl;
    std::cout << " Central Fit Run Log: Seed " << seed << std::endl;
    std::cout << "==================================================" << std::endl;
    pm.PrintParameters();

    // 4. Run minimization
    fitter.SetMinimizerTolerance(0.01);
    fitter.SetMinimizerPrintLevel(1);
    fitter.MinimizeChi2(pm);

    double chi2 = fitter.GetLastChi2();
    int status = fitter.GetMinimizerStatus();
    bool valid = fitter.IsMinimizerValid();

    // 5. Save result tree for this seed
    TString outRoot = Form("%s/result_seed_%d.root", outDir.Data(), seed);
    fitter.MakeResultTree(pm, outRoot);

    std::cout << Form("\n[RunSingleSeed] Seed %d Finished: Chi2 = %.2f | Status = %d | Valid = %d | Output = %s\n",
                      seed, chi2, status, (int)valid, outRoot.Data()) << std::endl;
}


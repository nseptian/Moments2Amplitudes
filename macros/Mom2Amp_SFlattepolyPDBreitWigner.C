#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include "TStopwatch.h"
#include <memory>
#include <vector>
#include <map>

const Double_t FIRST_MASS_BIN_CENTER = 0.82; // in GeV
const Double_t MASS_BIN_WIDTH = 0.04;
const Int_t N_MASS_BINS = 30;

TString FIT_RESULTS_DIR = "/d/grid17/septian/EtaPi0Analysis/run/MomentFit_m080200/t010020/";

void Mom2Amp_SFlattepolyPDBreitWigner(int seed=0){
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
        std::cerr << "Please check that the following directory exists and contains the required files:" << std::endl;
        std::cerr << "  " << FIT_RESULTS_DIR << std::endl;
        return;
    }

    // Generate mass bins
    std::vector<double> massBins = massDepMoments.GetMassBins();

    using namespace m2pw;
    using MDModel = MassDependentFitter::MassDependenceConfig::ModelType;
    using WaveModelConfig = MassDependentFitter::MassDependenceConfig::WaveModelConfig;
    
    // auto a0_980_pi1_1400_a2_1320_config = MassDependentFitter::CreateCustomConfig(
    //     {
    //         {0, {WaveModelConfig{"a0_980", MDModel::FlattePlusPolynomial, 2}}},
    //         {1, {WaveModelConfig{"pi1_1400", MDModel::BreitWigner}}},
    //         {2, {WaveModelConfig{"a2_1320", MDModel::BreitWigner}}}
    //     },
    //     {}
    // );
    auto a0_980_pi1_1400_a2_1320_1700_config = MassDependentFitter::CreateCustomConfig(
        {
            {0, {WaveModelConfig{"a0_980", MDModel::FlattePlusPolynomial, 2}}},
            {1, {WaveModelConfig{"pi1_1400", MDModel::BreitWigner}}},
            {2, {WaveModelConfig{"a2_1320_1700", MDModel::TwoBreitWigner}}}
        },
        {}
    );

    auto allMoments_config = MassDependentFitter::CreateIncludeAllConfig();

    MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, a0_980_pi1_1400_a2_1320_1700_config, allMoments_config};
    fitter.SetEquationValues(massDepMoments);
    
    fitter.PrintIncludedMoments();

    // const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_L4_only_0_sPlot.root";
    // const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_L4_only_0_nominal_t010030.root";
    // const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_BWa2_freeSWaves_ZeroPWaves_322_nominal_t010030.root";

    m2pw::MassDependentFitter::ParameterManager paramManager(fitter);
    
    // paramManager.AddMassDependentParameter(std::vector<int>{2}, 0, 0);
    paramManager.AddMassDependentParameter(std::vector<int>{0,1,2}, 0, 2);
    // paramManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
    // paramManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);
    paramManager.SetInitialValue("MD_shared_a2_1320_1700_Mass1", 1.318, true);
    paramManager.SetInitialValue("MD_shared_a2_1320_1700_Width1", 0.107, true);
    paramManager.SetInitialValue("MD_shared_a2_1320_1700_Mass2", 1.706, true);
    paramManager.SetInitialValue("MD_shared_a2_1320_1700_Width2", 0.380, true);
    paramManager.SetInitialValue("MD_shared_a0_980_Mass", 0.980, true);
    paramManager.SetInitialValue("MD_shared_a0_980_g_etapi", 0.353, true);
    paramManager.SetInitialValue("MD_shared_a0_980_g_KK", 0.0, true);
    paramManager.SetInitialValue("MD_a_0_0_a0_980_m_threshold", 0.547853 + 0.13957, true);
    paramManager.SetInitialValue("MD_b_0_0_a0_980_m_threshold", 0.547853 + 0.13957, true);
    paramManager.SetParameterLimits("MD_a_0_0_a0_980_m_expansion", 0.95, 1.2);
    paramManager.SetParameterLimits("MD_b_0_0_a0_980_m_expansion", 0.95, 1.2);
    paramManager.SetInitialValue("MD_shared_pi1_1400_Mass", 1.354, true);
    paramManager.SetInitialValue("MD_shared_pi1_1400_Width", 0.330, true);



    // TString a_1_m1_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_allMoments_BWa2_initializedSWaves_Pi1BWPWaves_139_nominal_t010030.root";

    // paramManager.AddMassIndependentParametersForPhase(a_1_m1_file,
    //                                             std::vector<TString>{"aphi_1_-1"});

    // paramManager.AddMassDependentParameter(std::vector<int>{2}, D_waves_model_file, 0);
    // paramManager.AddMassDependentParameter(std::vector<int>{1}, seed, 0);

    // paramManager.AddMassIndependentParametersForPhase(fitter.GetMassBins(),
                                                // fitter.GetParNames(),
                                                // std::vector<TString>{"bphi_1_-1"});

    // paramManager.AddMassDependentParametersForL(fitter.GetParNames(), 
    //                                       std::vector<TString>{"b_1_-1"}, 
    //                                       seed, 
    //                                       fitter.GetMassDependenceConfig(), 
    //                                       fitter.GetMomentsConfig(), 
    //                                       fitter,
    //                                       false,
    //                                       false,
    //                                       true);

    // paramManager.AddMassDependentParameters(fitter.GetParNames(), 
    //                                         D_waves_model_file, 
    //                                         fitter.GetMassDependenceConfig(), 
    //                                         fitter.GetMomentsConfig(), 
    //                                         fitter,0);

    // std::vector<string> fixedL = {"1-"};
    // std::map<string, double> fixedValues;
    // fixedValues["1-"] = 0.0; // Fix L=1 to 0

    // paramManager.AddFixedMassIndependentParametersForL(fitter.GetMassBins(),
    //                                                   fitter.GetParNames(),
    //                                                   fixedL,
    //                                                   fixedValues,
    //                                                   fitter.GetMomentsConfig(),
    //                                                   fitter,
    //                                                   false);
    
    // std::vector<int> initializedL = {0};
    // const TString S_waves_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_BWa2_freeSWaves_ZeroPWaves_322_nominal_t010030.root";
    // paramManager.AddMassIndependentParameter(initializedL, S_waves_file);

    // paramManager.AddMassIndependentParameters(fitter.GetMassBins(), 
    //                                           fitter.GetParNames(), 
    //                                           seed, 
    //                                           fitter.GetMassDependenceConfig(), 
    //                                           fitter.GetMomentsConfig(), 
    //                                           fitter);                                                      


    // fitter.PrintParNameIndices();
    // fitter.PrintIncludedMoments();

    paramManager.PrintParameters();
    fitter.SetMinimizerPrintLevel(2);
    fitter.SetEnableParallelization(true);
    fitter.SetNumThreads(2);

    fitter.MinimizeChi2(paramManager);

    fitter.MakeResultTree(paramManager, "result_SFlattepolyO2_PBW_D2BW_0.root");

}
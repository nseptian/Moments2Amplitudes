#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include "TStopwatch.h"
#include <memory>
#include <vector>
#include <map>

const Double_t FIRST_MASS_BIN_CENTER = 0.82; // in GeV
const Double_t MASS_BIN_WIDTH = 0.04;
const Int_t N_MASS_BINS = 20;

TString FIT_RESULTS_DIR = "/d/grid17/septian/EtaPi0Analysis/run/MomentFit_m080200/t010020/";

void Mom2Amp_SPiecewiseDBreitWigner(int seed=0){
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
    
    auto a2_1320_config = MassDependentFitter::CreateCustomConfig({{2, {WaveModelConfig{"a2_1320", MDModel::BreitWigner}}}}, {0, 1});
    // auto pi1_1400_a2_1320_config = MassDependentFitter::CreateCustomConfig(
    //     {
    //         {1, {WaveModelConfig{"pi1_1400", MDModel::BreitWigner}}},
    //         {2, {WaveModelConfig{"a2_1320", MDModel::BreitWigner}}}
    //     },
    //     {0});
    // auto H4_config = MassDependentFitter::CreateL4OnlyConfig();
    // auto H24_config = MassDependentFitter::CreateL2L4OnlyConfig();
    auto H024_config = MassDependentFitter::CreateL0L2L4OnlyConfig();
    // auto allMoments_config = MassDependentFitter::CreateIncludeAllConfig();

    MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, a2_1320_config, H024_config};
    fitter.SetEquationValues(massDepMoments);
    
    fitter.PrintIncludedMoments();

    // const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_L4_only_0_sPlot.root";
    // const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_L4_only_0_nominal_t010030.root";
    // const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_BWa2_freeSWaves_ZeroPWaves_322_nominal_t010030.root";

    m2pw::MassDependentFitter::ParameterManager paramManager(fitter);

    paramManager.AddMassIndependentParameter(std::vector<int>{0, 1}, seed);

    // Keep the ell=1 mass-independent parameters in the list, but fix them to zero.
    for (const double massBin : massBins) {
        for (const TString& parName : m2pw::MassDependentFitter::GetParNames()) {
            if (!parName.BeginsWith("a_1_") && !parName.BeginsWith("b_1_") &&
                !parName.BeginsWith("aphi_1_") && !parName.BeginsWith("bphi_1_")) {
                continue;
            }

            TString fixedName = Form("MI_%1.6f_%s", massBin, parName.Data());
            if (paramManager.SetInitialValue(fixedName, 0.0, true)) {
                continue;
            }
        }
    }

    paramManager.AddMassDependentParameter(std::vector<int>{2}, 0, 0);
    paramManager.SetInitialValue("MD_shared_a2_1320_Mass", 1.318, true);
    paramManager.SetInitialValue("MD_shared_a2_1320_Width", 0.107, true);

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
    fitter.MinimizeChi2(paramManager);

    // fitter.MakeResultTree(paramManager, "result_trsee_test_0.root");

}
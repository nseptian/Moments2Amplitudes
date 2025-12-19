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
// const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m080200_MCMCN6000BI1000S08WCOV_R6.34/";
// const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_Phase2/fitMoment_2019_11_polALL_sPlot_t010100_m080200_MCMC_R6.34/";
// const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase2_MomentsFit_t010020_Mpi0eta080200_800000_mandelstam_t010020_100200_MCMC_R6.34/";
const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase1_Phase2_MomentsFit_t010030_Mpi0eta080200_800000_mandelstam_t010030_100300_MCMC_R6.34/";

void H_odd_model(int seed=0){
    auto& setup = ConfigureAmpsNoValues(2, 2, 2); // (Lmax, MMax, Nref)
    
    // Load moments data with error handling
    MassDependentMoments massDepMoments;
    try {
        massDepMoments = MassDependentMoments(N_MASS_BINS, 
                                            MASS_BIN_WIDTH, 
                                            FIRST_MASS_BIN_CENTER, 
                                            FIT_RESULTS_DIR,
                                            "ResultsBruMcmcCovariance.root", 
                                            "Mpi0eta",
                                            "mandelstam_t0.200000");
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
    
    // auto a2_1320_config = MassDependentFitter::CreateDefaultConfig();
    auto pi1_1400_a2_1320_config = MassDependentFitter::CreateCustomConfig({{1, {"pi1_1400"}}, {2, {"a2_1320"}}}, {0});
    // auto H24_config = MassDependentFitter::CreateL2L4OnlyConfig();
    // auto H024_config = MassDependentFitter::CreateL0L2L4OnlyConfig();
    auto allMoments_config = MassDependentFitter::CreateIncludeAllConfig();

    MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, pi1_1400_a2_1320_config, allMoments_config};
    fitter.SetEquationValues(massDepMoments);
    
    fitter.PrintIncludedMoments();

    // const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_L4_only_0_sPlot.root";
    // const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_L4_only_0_nominal_t010030.root";
    const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_BWa2_freeSWaves_ZeroPWaves_322_nominal_t010030.root";

    m2pw::MassDependentFitter::ParameterManager paramManager(fitter);

    // TString a_1_m1_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_allMoments_BWa2_initializedSWaves_Pi1BWPWaves_139_nominal_t010030.root";

    // paramManager.AddMassIndependentParametersForPhase(a_1_m1_file,
    //                                             std::vector<TString>{"aphi_1_-1"});
    
    paramManager.AddMassDependentParametersForL(std::vector<int>{2}, 
                                          D_waves_model_file, 
                                          false,
                                          false,
                                          true);

    paramManager.AddMassDependentParametersForL(std::vector<int>{1}, 
                                          seed, 
                                          false,
                                          false,
                                          true);

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
    
    std::vector<int> initializedL = {0};
    const TString S_waves_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_BWa2_freeSWaves_ZeroPWaves_322_nominal_t010030.root";
    paramManager.AddMassIndependentParametersForL(fitter.GetMassBins(),
                                                fitter.GetParNames(),
                                                initializedL,
                                                S_waves_file);

    // paramManager.AddMassIndependentParameters(fitter.GetMassBins(), 
    //                                           fitter.GetParNames(), 
    //                                           seed, 
    //                                           fitter.GetMassDependenceConfig(), 
    //                                           fitter.GetMomentsConfig(), 
    //                                           fitter);                                                      


    fitter.PrintParNameIndices();
    fitter.PrintIncludedMoments();

    fitter.MinimizeChi2(paramManager);

    fitter.MakeResultTree(paramManager, "result_tree_allMoments_BWa2_initializedSWaves_Fit_AllBWPi1_"+std::to_string(seed)+"_nominal_t010030.root");

}
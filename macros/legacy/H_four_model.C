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
// const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase1_Phase2_MomentsFit_t010030_Mpi0eta080200_800000_mandelstam_t010030_100300_MCMC_R6.34/";
const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase1_Phase2_MomentsFit_t010030_Mpi0eta080200_800000_mandelstam_t010030_100300_MCMC_R6.34/";


void H_four_model(int seed = 0) {
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
        std::cerr << "Failed to load moments data!" << e.what() << std::endl;
        return;
    }

    // massDepMoments.PrintMoments();

    std::vector<double> massBins = massDepMoments.GetMassBins();

    using namespace m2pw;
    
    // Setup massdependentconfig default = a2_1320_config
    auto a2_1320_config = MassDependentFitter::CreateDefaultConfig();

    // Include only H4 moments in chi2
    auto H4_config = MassDependentFitter::CreateL4OnlyConfig();

    MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, a2_1320_config, H4_config};


    fitter.SetEquationValues(massDepMoments);

    fitter.PrintIncludedMoments();

    m2pw::MassDependentFitter::ParameterManager paramManager(fitter);
    paramManager.AddMassDependentParameter(std::vector<int>{2}, seed, 3);

    fitter.MinimizeChi2(paramManager);

    fitter.MakeResultTree(paramManager, "result_tree_L4_only_" + std::to_string(seed) + "_nominal_t010030.root");

}
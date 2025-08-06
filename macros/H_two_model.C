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
const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m080200_MCMCN6000BI1000S08WCOV_R6.34/";

void H_two_model(int seed=0){
    auto& setup = ConfigureAmpsNoValues(2, 2, 2); // (Lmax, MMax, Nref)
    
    // Load moments data
    MassDependentMoments massDepMoments(N_MASS_BINS, 
                                        MASS_BIN_WIDTH, 
                                        FIRST_MASS_BIN_CENTER, 
                                        FIT_RESULTS_DIR, 
                                        "Mpi0eta");
    
    // Generate mass bins
    std::vector<double> massBins = GenerateMassBins();

    using namespace m2pw;
    
    auto a2_1320_config = MassDependentFitter::CreateDefaultConfig();
    auto H24_config = MassDependentFitter::CreateL2L4OnlyConfig();

    MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, a2_1320_config, H24_config};
    fitter.SetEquationValues(massDepMoments);
    
    fitter.PrintIncludedMoments();

    const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/brufit/result_tree_L4_only_0.root";

    m2pw::MassDependentFitter::ParameterManager paramManager;
    paramManager.AddMassDependentParameters(fitter.GetParNames(), 
                                            D_waves_model_file, 
                                            fitter.GetMassDependenceConfig(), 
                                            fitter.GetHMomentsConfig(), 
                                            fitter,0);
    
    paramManager.AddMassIndependentParameters(fitter.GetMassBins(), 
                                              fitter.GetParNames(), 
                                              seed, 
                                              fitter.GetMassDependenceConfig(), 
                                              fitter.GetHMomentsConfig(), 
                                              fitter);

    // fitter.PrintParNameIndices();

    fitter.MinimizeChi2(paramManager, seed);

    fitter.MakeResultTree(paramManager, "result_tree_BWa2_freeSWaves_"+std::to_string(seed)+".root", seed);

}
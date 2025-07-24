#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include "TStopwatch.h"
#include <memory>
#include <vector>
#include <map>

namespace Config {
    const Double_t FIRST_MASS_BIN_CENTER = 0.82; // in GeV
    const Double_t MASS_BIN_WIDTH = 0.04;
    const Int_t N_MASS_BINS = 20;
    const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m080200_MCMCN6000BI1000S08WCOV_R6.34/";
    // const TString FIT_RESULTS_FILENAME = "ResultsBruMcmcCovariance.root";
    // const TString BW_NAMES[3] = {"k", "M", "width"};
} // namespace Config

MassDependentMoments LoadMomentsData() {
    MassDependentMoments massDepMoments;
    
    for (Int_t i = 0; i < Config::N_MASS_BINS; ++i) {
        Double_t massValue = Config::FIRST_MASS_BIN_CENTER + i * Config::MASS_BIN_WIDTH;
        const TString massBinDirName = Form("Mpi0eta%1.6f_/", massValue);
        const TString filePath = Config::FIT_RESULTS_DIR + massBinDirName;

        std::cout << "Loading moments from file: " << filePath << std::endl;

        MomentHelper moments;
        moments.Set(filePath, 1, 1);
        massDepMoments.SetMoments(massValue, moments);
    }
    
    return massDepMoments;
}

std::vector<double> GenerateMassBins() {
    std::vector<double> massBins;
    massBins.reserve(Config::N_MASS_BINS);
    
    for (Int_t i = 0; i < Config::N_MASS_BINS; ++i) {
        Double_t massValue = Config::FIRST_MASS_BIN_CENTER + i * Config::MASS_BIN_WIDTH;
        massBins.push_back(massValue);
    }
    
    return massBins;
}

void H_two_model(int seed=0){
    auto& setup = ConfigureAmpsNoValues(2, 2, 2); // (Lmax, MMax, Nref)
    
    // Load moments data
    MassDependentMoments massDepMoments = LoadMomentsData();
    
    // Generate mass bins
    std::vector<double> massBins = GenerateMassBins();

    using namespace m2pw;
    
    // Setup fitter
    // m2pw::MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}};
    // auto l4OnlyConfig = m2pw::MassDependentFitter::CreateL4OnlyConfig();
    auto a2_1320_config = MassDependentFitter::CreateDefaultConfig();

    std::cout << "Mass dependence configuration:" << std::endl;
    std::cout << "   - L=0 mass dependent: " << (a2_1320_config.IsMassDependent(0) ? "YES" : "NO") << std::endl;
    std::cout << "   - L=1 mass dependent: " << (a2_1320_config.IsMassDependent(1) ? "YES" : "NO") << std::endl;
    std::cout << "   - L=2 mass dependent: " << (a2_1320_config.IsMassDependent(2) ? "YES" : "NO") << std::endl;
    
    std::cout << "   - L=0 mass independent: " << (a2_1320_config.IsMassIndependent(0) ? "YES" : "NO") << std::endl;
    std::cout << "   - L=1 mass independent: " << (a2_1320_config.IsMassIndependent(1) ? "YES" : "NO") << std::endl;
    std::cout << "   - L=2 mass independent: " << (a2_1320_config.IsMassIndependent(2) ? "YES" : "NO") << std::endl;

    auto H4_config = MassDependentFitter::CreateL4OnlyConfig();
    auto H24_config = MassDependentFitter::CreateL2L4OnlyConfig();
    auto HAll_config = MassDependentFitter::CreateIncludeAllConfig();

    // MassDependentFitter::MassDependenceConfig a2_a0_config;
    // a2_a0_config.massDependentWaves[0] = {"a0_980"};
    // a2_a0_config.massDependentWaves[2] = {"a2_1320"};
    // a2_a0_config.massIndependentL = {1};

    MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, a2_1320_config, H24_config};
    fitter.SetEquationValues(massDepMoments);
    
    std::cout << "\n=== Include L=2,4 only configuration ===" << std::endl;
    // fitter.SetHMomentsConfig(H24_config);
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
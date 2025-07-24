#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include "TStopwatch.h"
#include <memory>
#include <vector>
#include <map>

// Configuration constants
namespace Config {
    const Double_t FIRST_MASS_BIN_CENTER = 0.82; // in GeV
    const Double_t MASS_BIN_WIDTH = 0.04;
    const Int_t N_MASS_BINS = 3;
    const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m080200_MCMCN6000BI1000S08WCOV_R6.34/";
    const TString FIT_RESULTS_FILENAME = "ResultsBruMcmcCovariance.root";
}

MassDependentMoments LoadMomentsData() {
    MassDependentMoments massDepMoments;
    
    for (Int_t i = 0; i < Config::N_MASS_BINS; ++i) {
        Double_t massValue = Config::FIRST_MASS_BIN_CENTER + i * Config::MASS_BIN_WIDTH;
        const TString massBinDirName = Form("Mpi0eta%1.6f_/", massValue);
        const TString filePath = Config::FIT_RESULTS_DIR + massBinDirName;

        std::cout << "Loading moments from file: " << filePath << std::endl;

        // std::unique_ptr<TFile> file(TFile::Open(filePath));
        // if (!file || file->IsZombie()) {
        //     std::cerr << "Error opening file: " << filePath << std::endl;
        //     continue;
        // }

        // Use stack allocation instead of new
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

void TestMassDependentMoments_Clean(int seed = 0) {

    auto& setup = ConfigureAmpsNoValues(2, 2, 2); // (Lmax, MMax, Nref)
    
    // Load moments data
    MassDependentMoments massDepMoments = LoadMomentsData();
    
    // Generate mass bins
    std::vector<double> massBins = GenerateMassBins();

    // Setup solver
    // m2pw::MassDependentEquationSolver solver{setup, massBins, 2, 0.0, {"H_3"}};

    auto l4OnlyConfig = m2pw::MassDependentEquationSolver::CreateL4OnlyConfig();
    auto L2L4OnlyConfig = m2pw::MassDependentEquationSolver::CreateL2L4OnlyConfig();
    auto massDepConfig = m2pw::MassDependentEquationSolver::CreateDefaultConfig();

    m2pw::MassDependentEquationSolver solver{setup, massBins, 2, 0.0, {"H_3"}, massDepConfig, L2L4OnlyConfig};
    solver.SetEquationValues(massDepMoments);

    const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/brufit/result_tree_L4_only.root";

    // --- New interface demonstration ---
    // 1. Load initial values from a result tree (if available)
    m2pw::MassDependentEquationSolver::ParameterManager paramManager;
    paramManager.AddMassDependentParameters(solver.GetParNames(), 
                                            D_waves_model_file, 
                                            solver.GetMassDependenceConfig(), 
                                            solver.GetHMomentsConfig(), 
                                            solver,1);
    
    paramManager.AddMassIndependentParameters(solver.GetMassBins(), 
                                              solver.GetParNames(), 
                                              seed, 
                                              solver.GetMassDependenceConfig(), 
                                              solver.GetHMomentsConfig(), 
                                              solver);

    solver.MinimizeChi2(paramManager, seed);

    // solver.RandomInit(paramManager, seed);
    // Uncomment and set the correct path if you want to test loading initial values
    // solver.LoadInit("/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m080200_MCMCN6000BI1000S08WCOV_R6.34/ResultsBruMcmcCovariance.root");

    // 2. (Parameter fixing interface not yet implemented)

    // 3. Run minimization with the new interface
    // solver.MinimizeChi2(paramManager);
}

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

void H_four_model(int seed = 0) {
    auto& setup = ConfigureAmpsNoValues(2, 2, 2); // (Lmax, MMax, Nref)
    
    // Load moments data
    MassDependentMoments massDepMoments = LoadMomentsData();
    
    // Generate mass bins
    std::vector<double> massBins = GenerateMassBins();

    using namespace m2pw;
    
    // Setup solver
    // m2pw::MassDependentEquationSolver solver{setup, massBins, 2, 0.0, {"H_3"}};
    // auto l4OnlyConfig = m2pw::MassDependentEquationSolver::CreateL4OnlyConfig();
    auto a2_1320_config = MassDependentEquationSolver::CreateDefaultConfig();
    auto H4_config = MassDependentEquationSolver::CreateL4OnlyConfig();

    // std::cout << "Mass dependence configuration:" << std::endl;
    // std::cout << "   - L=0 mass dependent: " << (a2_1320_config.IsMassDependent(0) ? "YES" : "NO") << std::endl;
    // std::cout << "   - L=1 mass dependent: " << (a2_1320_config.IsMassDependent(1) ? "YES" : "NO") << std::endl;
    // std::cout << "   - L=2 mass dependent: " << (a2_1320_config.IsMassDependent(2) ? "YES" : "NO") << std::endl;
    
    // std::cout << "   - L=0 mass independent: " << (a2_1320_config.IsMassIndependent(0) ? "YES" : "NO") << std::endl;
    // std::cout << "   - L=1 mass independent: " << (a2_1320_config.IsMassIndependent(1) ? "YES" : "NO") << std::endl;
    // std::cout << "   - L=2 mass independent: " << (a2_1320_config.IsMassIndependent(2) ? "YES" : "NO") << std::endl;

    MassDependentEquationSolver solver{setup, massBins, 2, 0.0, {"H_3"}, a2_1320_config, H4_config};
    solver.SetEquationValues(massDepMoments);

    // auto H4_config = MassDependentEquationSolver::CreateL4OnlyConfig();
    // auto H24_config = MassDependentEquationSolver::CreateL2L4OnlyConfig();
    // auto HAll_config = MassDependentEquationSolver::CreateIncludeAllConfig();
    
    // Test with L=4 only
    std::cout << "\n=== Testing with L=4 only configuration ===" << std::endl;
    solver.SetHMomentsConfig(H4_config);
    // solver.PrintIncludedMoments();

    m2pw::MassDependentEquationSolver::ParameterManager paramManager;
    paramManager.AddMassDependentParameters(solver.GetParNames(), 
                                            seed, 
                                            solver.GetMassDependenceConfig(), 
                                            solver.GetHMomentsConfig(), 
                                            solver);

    solver.MinimizeChi2(paramManager, seed);

    solver.MakeResultTree(paramManager, "result_tree_L4_only_" + std::to_string(seed) + ".root", seed);

    // TStopwatch timer;
    // timer.Start();
    // std::cout << "\n=== Evaluating chi2 with L=4 only ===" << std::endl;
    // Double_t chi2_L4_only = solver.EvalChi2();
    // Double_t eval_time_L4 = timer.RealTime();
    // std::cout << "Chi2 (L=4 only): " << chi2_L4_only << " (evaluation time: " << eval_time_L4 << " seconds)" << std::endl;
    
    // // Test with L=2,4 configuration
    // std::cout << "\n=== Testing with L=2,4 configuration ===" << std::endl;
    // solver.SetHMomentsConfig(H24_config);
    // solver.PrintIncludedMoments();
    
    // timer.Start();
    // std::cout << "\n=== Evaluating chi2 with L=2,4 ===" << std::endl;
    // Double_t chi2_L24 = solver.EvalChi2();
    // Double_t eval_time_L24 = timer.RealTime();
    // std::cout << "Chi2 (L=2,4): " << chi2_L24 << " (evaluation time: " << eval_time_L24 << " seconds)" << std::endl;
    
    // Test with all L values
    // std::cout << "\n=== Testing with ALL L values ===" << std::endl;
    // solver.SetHMomentsConfig(HAll_config);
    // solver.PrintIncludedMoments();
    
    // timer.Start();
    // std::cout << "\n=== Evaluating chi2 with ALL L values ===" << std::endl;
    // Double_t chi2_all = solver.EvalChi2();
    // Double_t eval_time_all = timer.RealTime();
    // std::cout << "Chi2 (ALL): " << chi2_all << " (evaluation time: " << eval_time_all << " seconds)" << std::endl;
    
    // Summary comparison
    // std::cout << "\n=== Chi2 Comparison ===" << std::endl;
    // std::cout << "L=4 only:  " << chi2_L4_only << std::endl;
    // std::cout << "L=2,4:     " << chi2_L24 << std::endl;
    // std::cout << "ALL L:     " << chi2_all << std::endl;

    // // Reset timer for minimization
    // timer.Start();
    
    // std::cout << "\n=== Running minimization ===" << std::endl;
    // // Run minimization
    // solver.MinimizeChi2(seed);
    
    // Double_t chi2_final = solver.GetLastChi2();
    // Double_t minimization_time = timer.RealTime();
    
    // std::cout << "\n=== Results ===" << std::endl;
    // std::cout << "Chi2 before minimization: " << chi2_initial << std::endl;
    // std::cout << "Chi2 after minimization:  " << chi2_final << std::endl;
    // std::cout << "Chi2 improvement:         " << (chi2_initial - chi2_final) << std::endl;
    // std::cout << "Minimization time:        " << minimization_time << " seconds" << std::endl;


}
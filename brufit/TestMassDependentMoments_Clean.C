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
    const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m010200_MCMCN6000BI1000S08WCOV_R6.34/";
    const TString FIT_RESULTS_FILENAME = "ResultsBruMcmcCovariance.root";
}

MassDependentMoments LoadMomentsData() {
    MassDependentMoments massDepMoments;
    
    for (Int_t i = 0; i < Config::N_MASS_BINS; ++i) {
        Double_t massValue = Config::FIRST_MASS_BIN_CENTER + i * Config::MASS_BIN_WIDTH;
        const TString massBinDirName = Form("Mpi0eta%1.6f_/", massValue);
        const TString filePath = Config::FIT_RESULTS_DIR + massBinDirName + Config::FIT_RESULTS_FILENAME;

        std::cout << "Loading moments from file: " << filePath << std::endl;

        std::unique_ptr<TFile> file(TFile::Open(filePath));
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filePath << std::endl;
            continue;
        }

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
    m2pw::MassDependentEquationSolver solver{setup, massBins, 2, 0.0, {"H_3"}};
    solver.SetEquationValues(massDepMoments);
    
    // Run minimization with new interface - all functionality is now in the solver
    solver.MinimizeChi2(seed);
}

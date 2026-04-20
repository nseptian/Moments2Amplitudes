#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include "TStopwatch.h"
#include <memory>
#include <vector>
#include <map>
#include <fstream>
#include <ROOT/TProcessExecutor.hxx>

// Struct to hold seed task parameters
struct SeedTask {
    int seed;
    Double_t lastMassBinCenter;
    Int_t nMassBins;
};

// const Double_t FIRST_MASS_BIN_CENTER = 0.82; // in GeV
const Double_t MASS_BIN_WIDTH = 0.04;
// const Int_t N_MASS_BINS = 20;
// const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m080200_MCMCN6000BI1000S08WCOV_R6.34/";
// const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_Phase2/fitMoment_2019_11_polALL_sPlot_t010100_m080200_MCMC_R6.34/";
// const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase2_MomentsFit_t010020_Mpi0eta080200_800000_mandelstam_t010020_100200_MCMC_R6.34/";
const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/FitMoments/fitMoment_Phase1_Phase2_MomentsFit_t010030_Mpi0eta080200_800000_mandelstam_t010030_100300_MCMC_R6.34/";
const TString OUTPUT_DIR = "/d/home/septian/Moments2Amplitudes/run/Mom2Amps_H_two_model_seq_init_forward/";
const TString LOG_DIR = "/d/home/septian/Moments2Amplitudes/run/Mom2Amps_H_two_model_seq_init_forward/logs/";

void H_two_model(const Double_t FIRST_MASS_BIN_CENTER, const Int_t N_MASS_BINS, const TString init_model_file="", const TString result_file="") {
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
    std::vector<double> initializedMassBins;
    std::vector<double> loadedMassBins;

    for (int iMassBin = 0; iMassBin < massBins.size(); ++iMassBin) {
        if (iMassBin == massBins.size()-1) {
            initializedMassBins.push_back(massBins[iMassBin]);
        } else {
            cout << "Added Loaded Mass bin " << iMassBin << " center: " << massBins[iMassBin] << endl;
            loadedMassBins.push_back(massBins[iMassBin]);
        }
    }

    using namespace m2pw;
    
    auto a2_1320_config = MassDependentFitter::CreateDefaultConfig();
    // auto a0_flatte_a2_1320_config = MassDependentFitter::CreateCustomConfig({{0, {"a0_980"}}, {2, {"a2_1320"}}}, {});
    // auto S_waves_conformalPol_a2_1320_config = MassDependentFitter::CreateCustomConfig({{0, {"conformal_poly_order3"}}, {2, {"a2_1320"}}}, {});
    // auto H24_config = MassDependentFitter::CreateL2L4OnlyConfig();
    auto H024_config = MassDependentFitter::CreateL0L2L4OnlyConfig();
    // auto allMoments_config = MassDependentFitter::CreateIncludeAllConfig();

    MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, a2_1320_config, H024_config};
    // MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, a0_flatte_a2_1320_config, H024_config};
    // MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}, S_waves_conformalPol_a2_1320_config, H024_config};

    fitter.SetEquationValues(massDepMoments);
    
    fitter.PrintIncludedMoments();

    // const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_L4_only_0_sPlot.root";
    // const TString D_waves_model_file = "/d/home/septian/Moments2Amplitudes/macros/result_tree_L4_only_0_nominal_t010030.root";

    // const TString init_model_file = "/d/home/septian/Moments2Amplitudes/run/H2_model_BWa2_freeSWaves_ZeroPWaves/result_tree_BWa2_freeSWaves_ZeroPWaves_322_nominal_t010030.root";

    m2pw::MassDependentFitter::ParameterManager paramManager(fitter);

    paramManager.AddMassDependentParameter(std::vector<int>{2}, init_model_file, 0);

    // paramManager.AddMassDependentParametersForL(std::vector<int>{0}, 
    //                               seed, 
    //                               false,
    //                               false,
    //                               true);

    // std::vector<int> fixedL = {1};
    // std::map<int, double> fixedValues;
    // fixedValues[1] = 0.0;        // Fix L=1 to 0

    // paramManager.AddFixedMassIndependentParametersForL(fitter.GetMassBins(),
    //                                                   fitter.GetParNames(),
    //                                                   fixedL,
    //                                                   fixedValues,
    //                                                   fitter.GetMomentsConfig(),
    //                                                   fitter,
    //                                                   false);

    std::vector<int> initializedMassIndependentL = {0};


    paramManager.AddMassIndependentParameter(initializedMassIndependentL, init_model_file);

    // std::cout << "\n=== Running H_two_model ===" << std::endl;
    // std::cout << "FIRST_MASS_BIN_CENTER = " << FIRST_MASS_BIN_CENTER << std::endl;
    // std::cout << "N_MASS_BINS = " << N_MASS_BINS << std::endl;

    fitter.PrintParNameIndices();
    fitter.PrintIncludedMoments();

    fitter.MinimizeChi2(paramManager);

    fitter.MakeResultTree(paramManager, result_file);

    std::cout << "=== Completed ===" << std::endl;

    // Close all open ROOT files
    TSeqCollection* files = gROOT->GetListOfFiles();
    if (files) {
        TIter next(files);
        TFile* file;
        while ((file = (TFile*)next())) {
            if (file && file->IsOpen()) {
                file->Close();
            }
        }
    }
    
}

void H_two_model_seq_init_forward() {

    std::vector<Double_t> lastMassBinCenters;

    cout << "Generating last mass bin centers from 1.3 to 1.58 with step 0.04" << endl;

    for (Double_t center = 1.34; center <= 1.6; center += 0.04) {
        lastMassBinCenters.push_back(center);
        cout << "  " << center << endl;
    }
    // Create log directory
    // if (gSystem->AccessPathName(LOG_DIR)) {
    //     gSystem->Exec(Form("mkdir -p %s", LOG_DIR.Data()));
    // }

    TString init_model_file = "/d/home/septian/Moments2Amplitudes/run/Mom2Amps_H_two_model_seq_init/result_tree_H_two_model_seq_init_nominal_t010030_firstbin_0.860000.root";

    if (gSystem->AccessPathName(OUTPUT_DIR)) {
        gSystem->Exec(Form("mkdir -p %s", OUTPUT_DIR.Data()));
    }

    TString output_file = "";

    for (Int_t iMassBin = 0; iMassBin < lastMassBinCenters.size(); ++iMassBin) {
        
        if (iMassBin > 0) {
            init_model_file = output_file;
        }

        Double_t lastMassBinCenter = lastMassBinCenters[iMassBin];
        Int_t N_MASS_BINS = TMath::FloorNint((lastMassBinCenter - 0.86)/MASS_BIN_WIDTH) + 1;

        std::cout << "\n========== Starting H_two_model_seq_init ==========" << std::endl;
        std::cout << "Last Mass Bin Center: " << lastMassBinCenter << std::endl;
        std::cout << "N_MASS_BINS: " << N_MASS_BINS << std::endl;
        std::cout << "================================================" << std::endl;

        output_file = OUTPUT_DIR + TString("result_tree_H_two_model_seq_init_forward_nominal_t010030_lastbin_" + std::to_string(lastMassBinCenters[iMassBin]) + ".root");
        
        H_two_model(0.86, N_MASS_BINS, init_model_file, output_file);
        
    }
}
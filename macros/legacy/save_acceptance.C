

namespace Config {
    const Double_t FIRST_MASS_BIN_CENTER = 0.82; // in GeV
    const Double_t MASS_BIN_WIDTH = 0.04;
    const Int_t N_MASS_BINS = 30;
    // const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m080200_MCMCN6000BI1000S08WCOV_R6.34/";
    // const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_Phase2/fitMoment_2019_11_polALL_sPlot_t010100_m080200_MCMC_R6.34/";
    // const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase2_MomentsFit_t010020_Mpi0eta080200_800000_mandelstam_t010020_100200_MCMC_R6.34/";
    // const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase1_Phase2_MomentsFit_t010030_Mpi0eta080200_800000_mandelstam_t010030_100300_MCMC_R6.34/";
    const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase1_Phase2_MomentsFit_nominal_Mpi0eta080200_800000_mandelstam_t010020_100200_MCMC_R6.34/";
} // namespace Config

void save_acceptance() {

    vector<double> massBinCenters;
    vector<double> acceptance;

    for (Int_t i = 0; i < Config::N_MASS_BINS; ++i) {
        Double_t massValue = Config::FIRST_MASS_BIN_CENTER + i * Config::MASS_BIN_WIDTH;
        const TString massBinDirName = Form("Mpi0eta%1.6f_mandelstam_t0.150000_/", massValue);
        const TString filePath = Config::FIT_RESULTS_DIR + massBinDirName + "ResultsCrossSection.root";

        std::cout << "Loading acceptance from file: " << filePath << std::endl;

        TFile *f = TFile::Open(filePath);
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening file: " << filePath << std::endl;
            continue;

        }
        HS::FIT::CrossSection *cross_section = (HS::FIT::CrossSection *)f->Get("cs");
        if (!cross_section) {
            std::cerr << "Error: CrossSection not found in file " << filePath << std::endl;
            continue;
        }

        double acc = cross_section->GetAcceptance();
        if (acc <= 0) {
            std::cerr << "Error: Invalid acceptance value in file " << filePath << std::endl;
            continue;
        }
        massBinCenters.push_back(massValue);
        acceptance.push_back(acc);
    }

    // Save the acceptance data to a file
    TString outputFilePath = Config::FIT_RESULTS_DIR + "acceptance.root";
    TFile *outputFile = TFile::Open(outputFilePath, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outputFilePath << std::endl;
        return;
    }

    TTree *tree = new TTree("acceptance", "Acceptance Data");
    double mass_bin, acceptance_value;

    tree->Branch("mass", &mass_bin, "mass/D");
    tree->Branch("acceptance", &acceptance_value, "acceptance/D");


    for (size_t i = 0; i < massBinCenters.size(); ++i) {
        mass_bin = massBinCenters[i];
        acceptance_value = acceptance[i];
        tree->Fill();
    }

    outputFile->Write();
    outputFile->Close();

    std::cout << "Acceptance data saved to: " << outputFilePath << std::endl;

}
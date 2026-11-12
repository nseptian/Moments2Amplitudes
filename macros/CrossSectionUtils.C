// CrossSectionEvaluator is already provided by BruFit session
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <iostream>
#include <vector>
#include <memory>

struct CrossSectionConfig {
    TString result_file_name = "ResultsBruMcmcCovariance.root";
    TString analysis_tag;
    TString input_data_dir;
    TString fit_result_dir;
    TString luminosity_file;
    Double_t total_branching_fraction = 0.3936*0.98823;  // BR(eta->gamma gamma) * BR(pi0->gamma gamma)
    std::vector<Double_t> beam_energy_limit = {8.0, 8.8};

    TString weights_file;
    TString weight_species = "Acc";
    TString weight_name = "weight";

    std::vector<TString> bin_variables = {"Mpi0eta"};
    std::vector<Int_t> n_bins = {30};
    std::vector<Double_t> bin_min = {0.8};
    std::vector<Double_t> bin_max = {2.0};

    int lmax = 2;
    int mmax = 2;
    int nrefl = 2;

    // Polarization configuration (do not change)
    double pol_value = 0.35;
    double pol_min = 0.34;
    double pol_max = 0.38;
};

void CalculateCrossSection(const CrossSectionConfig& config) {
    using namespace HS::FIT;

    CrossSectionEvaluator RF;
    RF.SetUp().SetOutDir(config.fit_result_dir);
    RF.SetResultDir(config.fit_result_dir);
    RF.SetResultFileName(config.result_file_name);

    PhotoTwoSpin0Amps amp_config("Moments");
    amp_config.SetManager(&RF);
    amp_config.SetDecayAngleCosTh("cosTheta_eta_gj[0.21,-1,1]");
    amp_config.SetDecayAnglePhi("phi_eta_gj[0.0,-3.2,3.2]");
    amp_config.SetPolPhi("Phi[0.2,-3.14159,3.14159]");
    amp_config.SetPolarisation(Form("Pol[%.2f,%.2f,%.2f]", config.pol_value, config.pol_min, config.pol_max));

    RF.SetUp().SetIDBranchName("UID");
    // Setup binning
    for (size_t i = 0; i < config.bin_variables.size(); ++i) {
        RF.Bins().LoadBinVar(config.bin_variables[i], config.n_bins[i], config.bin_min[i], config.bin_max[i]);
    }

    RF.Data().LoadWeights(config.weight_species, config.input_data_dir + config.weights_file, config.weight_name);
    RF.ReloadData(config.fit_result_dir);
    RF.ReloadSimulated(config.fit_result_dir, amp_config.GetName());
    RF.ReloadGenerated(config.fit_result_dir, amp_config.GetName());

    amp_config.SetLmax(config.lmax);
    amp_config.SetMmax(config.mmax);
    amp_config.SetNrefl(config.nrefl);
    amp_config.ConfigureMoments();
    amp_config.LoadModelPDF();

    // Integrated flux calculation from tagged_lumi histogram
    Double_t integrated_flux = 1.0;
    if (!config.luminosity_file.IsNull()) {
        std::unique_ptr<TFile> fluxFile(TFile::Open(config.luminosity_file));
        if (fluxFile && !fluxFile->IsZombie()) {
            auto* hFlux = dynamic_cast<TH1*>(fluxFile->Get("tagged_lumi"));
            if (hFlux) {
                Int_t b1 = hFlux->FindBin(config.beam_energy_limit[0]);
                Int_t b2 = hFlux->FindBin(config.beam_energy_limit[1]) - 1;
                integrated_flux = hFlux->Integral(b1, b2);
                std::cout << "Integrated flux [" << config.beam_energy_limit[0] << " - " 
                          << config.beam_energy_limit[1] << " GeV]: " << integrated_flux << std::endl;
            } else {
                std::cerr << "Warning: 'tagged_lumi' histogram not found in " << config.luminosity_file << std::endl;
            }
        } else {
            std::cerr << "Warning: Could not open luminosity file " << config.luminosity_file << std::endl;
        }
    }

    Double_t br = config.total_branching_fraction;
    Double_t target_thickness = 1.0; // inverse barn

    // Define cross-section formula: sigma = Yield / (Acceptance * Flux * TargetThickness * BR * BinWidth)
    RF.SetCrossSectionFormula([=](const CSData& d) -> std::pair<Double_t, Double_t> {
        if (d.acceptance <= 0.0 || integrated_flux <= 0.0) return {0.0, 0.0};
        Double_t bin_width = 1.0;
        if (d.kinematics.count("Mpi0eta")) {
            bin_width = d.kinematics.at("Mpi0eta").width;
        }
        Double_t denom = d.acceptance * integrated_flux * target_thickness * br * bin_width;
        Double_t cs = d.yield / denom;
        Double_t cs_err = cs * d.GetStatRelError(); // Combines relative errors of yield and acceptance
        return {cs, cs_err};
    });

    RF.SetNThreads(4);

    Here::Go(&RF);
}

void SaveAcceptance(const CrossSectionConfig& config) {
    std::vector<double> massBinCenters;
    std::vector<double> acceptance;
    
    // Calculate first_mass_bin_center from bin_min
    double bin_min = config.bin_min[0];
    double bin_max = config.bin_max[0];
    int n_bins = config.n_bins[0];
    double bin_width = (bin_max - bin_min) / n_bins;
    for (int i = 0; i < n_bins; ++i) {
        double massValue = bin_min + i * bin_width + bin_width / 2.0;
        const TString massBinDirName = Form("Mpi0eta%1.6f_/", massValue);
        const TString filePath = config.fit_result_dir + massBinDirName + "ResultsCrossSection.root";
        std::cout << "Loading acceptance from file: " << filePath << std::endl;
        TFile *f = TFile::Open(filePath);
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening file: " << filePath << std::endl;
            if (f) delete f;
            continue;
        }

        double acc = 0.0;
        // Try reading from modern CrossSectionEvaluator TTree first
        if (auto* tree = f->Get<TTree>("CrossSection")) {
            tree->SetBranchAddress("Acceptance", &acc);
            if (tree->GetEntries() > 0) tree->GetEntry(0);
        } else if (auto* cross_section = (HS::FIT::CrossSection*)f->Get("cs")) {
            // Fallback to legacy CrossSection object
            acc = cross_section->GetAcceptance();
        } else {
            std::cerr << "Error: Neither 'CrossSection' tree nor 'cs' object found in " << filePath << std::endl;
            f->Close();
            delete f;
            continue;
        }
        f->Close();
        delete f;

        if (acc <= 0) {
            std::cerr << "Error: Invalid acceptance value (" << acc << ") in file " << filePath << std::endl;
            continue;
        }
        massBinCenters.push_back(massValue);
        acceptance.push_back(acc);
    }

    TString outputFilePath = config.fit_result_dir + "acceptance.root";
    TFile *outputFile = TFile::Open(outputFilePath, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outputFilePath << std::endl;
        if (outputFile) delete outputFile;
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
    delete outputFile;
    std::cout << "Acceptance data saved to: " << outputFilePath << std::endl;
}

// Configuration structure for Moment fit
struct MomentFitConfig {
    // Input/Output paths
    TString analysis_tag;
    TString input_data_dir;
    TString output_dir;
    TString data_file;
    TString simulated_file;
    TString generated_file;
    TString weights_file;
    TString weight_species = "Acc"; 
    TString weight_name = "weight";  
    TString tree_name = "kin";
    
    // Processing options
    bool use_generated_mc = true;  // Optional acceptance correction for cross-section
    
    // Fit configuration
    vector<TString> bin_variables = {"Mpi0eta"};  // Branch names for binning
    vector<int> n_bins = {30};
    vector<double> bin_min = {0.8};    
    vector<double> bin_max = {2.0};

    int lmax = 2;
    int mmax = 2;
    int nrefl = 2;
    bool use_gj_frame = false;  // if false will use helicity frame, if true will use Gottfried-Jackson frame
    
    // Polarization configuration (do not change)
    double pol_value = 0.35;
    double pol_min = 0.34;
    double pol_max = 0.38;
    
    // MCMC parameters
    vector<int> mcmc_iterations = {10000, 10000, 3000};
    int mcmc_burnin = 2000;
    double mcmc_step = 0.8;
    double mcmc_target_acceptance = 0.14;
    double mcmc_min_acceptance = 0.1;
    double mcmc_max_acceptance = 0.18;
    bool mcmc_use_covariance = true;  // Whether to use covariance matrix in MCMC proposal
    bool mcmc_use_covariance_tuning = true;
    TString mcmc_tune_mode = "MappedRhat";
    int mcmc_rhat_window = 1000;
    double mcmc_rhat_target = 1.05;
    double mcmc_ess_fraction = 0.02;
    int mcmc_rhat_max_retries = 5;
    int mcmc_gibbs_block_size = 0;
    double mcmc_cold_warmup_acceptance = 0.15;
    
    // Bailout parameters (fast-fail if acceptance stays below threshold)
    int mcmc_bailout_threshold = 10000;
    double mcmc_bailout_min_acceptance = 0.0001; // Fast bailout if zero acceptance after 10,000 sampling proposals
    
    // New likelihood & PDF features
    bool apply_mc_variance = true;
    bool allow_negative_pdf = true;
    int mc_integration_bins = 50;

    // Use fixed covariance matrix from linear algebra methods
    bool fixed_covariance = false;
    TString covariance_file = "covariance.root";
    TString covariance_matrix_name = "mymatrix";

    bool set_init_value = true;
    TString init_value_file = "file.root";
    double init_value_scale = 0.5;
    
    // Processing options
    int num_cpu = 1;
    int Multi_workers = 80;
    bool run_fit = true;
    bool use_multi = true;
    bool verbose = false;
};

inline TString get_full_path(const TString& dir, const TString& file) {
    if (file.BeginsWith("/") || file.BeginsWith("~")) {
        return file;
    }
    TString d = dir;
    if (!d.EndsWith("/")) d += "/";
    return d + file;
}

void print_config(const MomentFitConfig& config) {
    cout << "\n=== BRUFIT MOMENT FIT CONFIGURATION ===" << endl;
    cout << "Analysis Tag: " << config.analysis_tag << endl;
    cout << "Output Directory: " << config.output_dir << endl;
    cout << "Input Files Directory: " << config.input_data_dir << endl;
    cout << "\nInput Files:" << endl;
    if (config.use_generated_mc) {
        cout << "  Generated MC (mcgen):    " << get_full_path(config.input_data_dir, config.generated_file) << endl;
    } else {
        cout << "  Generated MC (mcgen):    Not used (acceptance correction disabled)" << endl;
    }
    cout << "  Simulated MC (mcrecon):  " << get_full_path(config.input_data_dir, config.simulated_file) << endl;
    cout << "  Data (data):             " << get_full_path(config.input_data_dir, config.data_file) << endl;
    cout << "  Weights (weight):        " << get_full_path(config.input_data_dir, config.weights_file) << endl;
    cout << "  Tree Name:               " << (config.tree_name.IsNull() ? "kin" : config.tree_name) << endl;
    cout << "  Weight Species / Name:   " << config.weight_species << " / " << config.weight_name << endl;
    
    cout << "\nBinning variables: " << endl;
    for (size_t i = 0; i < config.bin_variables.size(); ++i) {
        cout << config.bin_variables[i] << " (" << config.n_bins[i] 
             << " bins, [" << config.bin_min[i] << ", " << config.bin_max[i] << "])";
        if (i < config.bin_variables.size() - 1) cout << ", ";
    }
    cout << endl;

    cout << "Fit Configuration:" << endl;
    cout << "  L_max: " << config.lmax << ", M_max: " << config.mmax 
         << ", N_refl: " << config.nrefl << endl;
    cout << endl;

    cout << "Markov Chain Monte Carlo Configuration: " << endl;
    cout << "  Number of samples: ";
    for (size_t i = 0; i < config.mcmc_iterations.size(); ++i) {
        cout << config.mcmc_iterations[i];
        if (i < config.mcmc_iterations.size() - 1) cout << ", ";
    }
    cout << endl;
    cout << "  Burn-in: " << config.mcmc_burnin << endl;
    cout << "  Scale step size: " << config.mcmc_step << endl;
    cout << "  Acceptance: target=" << config.mcmc_target_acceptance 
         << ", range=[" << config.mcmc_min_acceptance << ", " << config.mcmc_max_acceptance << "]" << endl;
    cout << "  Bailout: threshold=" << config.mcmc_bailout_threshold 
         << ", min_acceptance=" << config.mcmc_bailout_min_acceptance << endl;
    cout << "  Tuning: " << (config.mcmc_use_covariance_tuning ? "Enabled" : "Disabled")
         << " (Mode: " << config.mcmc_tune_mode << ")" << endl;
    cout << "  MC Variance: " << (config.apply_mc_variance ? "Enabled" : "Disabled") << endl;
    cout << "  Allow Negative PDF: " << (config.allow_negative_pdf ? "Enabled" : "Disabled") << endl;

    cout << endl << "Processing options:" << endl;
    cout << "  Use Multi: " << (config.use_multi ? "Yes" : "No") << endl;
    // if (config.use_multi) {
    //     cout << "  Multi workers: " << config.Multi_workers << endl;
    // }

    cout << "=======================================" << endl;
}

bool check_input_files(const MomentFitConfig& config) {
    cout << "\n=== CHECKING INPUT FILES ===" << endl;
    
    vector<pair<TString, TString>> files = {
        {get_full_path(config.input_data_dir, config.simulated_file), "Simulated MC (mcrecon)"}, 
        {get_full_path(config.input_data_dir, config.data_file), "Data (data)"},
        {get_full_path(config.input_data_dir, config.weights_file), "Data Weights (weight)"}
    };
    
    // Add generated MC only if enabled
    if (config.use_generated_mc) {
        files.insert(files.begin(), {get_full_path(config.input_data_dir, config.generated_file), "Generated MC (mcgen)"});
    }
    
    bool all_exist = true;
    for (const auto& file_pair : files) {
        TString filepath = file_pair.first;
        TString description = file_pair.second;
        
        if (gSystem->AccessPathName(filepath)) {
            cout << "  MISSING: " << description << endl;
            cout << "  Path: " << filepath << endl;
            all_exist = false;
        } else {
            cout << "  Found: " << description << endl;
            
            // Check if ROOT file can be opened
            TFile* test_file = TFile::Open(filepath, "READ");
            if (test_file && !test_file->IsZombie()) {
                cout << "    Size: " << (test_file->GetSize() / (1024*1024)) << " MB" << endl;
                test_file->Close();
            } else {
                cout << "    ERROR: Cannot open ROOT file" << endl;
                all_exist = false;
                if (test_file) test_file->Close();
            }
        }
    }
    
    cout << "=============================" << endl;
    return all_exist;
}

void configure_variables(PhotoTwoSpin0Amps& config, const MomentFitConfig& fit_config) {
    cout << "\n=== CONFIGURING AMPLITUDE MODEL ===" << endl;
    
    // Set decay angles and polarization
    if (fit_config.use_gj_frame) {
        config.SetDecayAngleCosTh("cosTheta_eta_gj[0.21,-1,1]");
        config.SetDecayAnglePhi("phi_eta_gj[0.0,-3.2,3.2]");  
        cout << "Using Gottfried-Jackson frame for decay angles" << endl;
        cout << "  Decay angles: cosTheta_eta_gj, phi_eta_gj" << endl;
    } else {
        config.SetDecayAngleCosTh("cosTheta_eta_hel[0.21,-1,1]");
        config.SetDecayAnglePhi("phi_eta_hel[0.0,-3.2,3.2]");  
        cout << "Using helicity frame for decay angles" << endl;
        cout << "  Decay angles: cosTheta_eta_hel, phi_eta_hel" << endl;
    }
    config.SetPolPhi("Phi[0.2,-3.2,3.2]");
    
    // Set polarization parameter
    TString pol_param = Form("Pol[%.3f,%.3f,%.3f]", 
                            fit_config.pol_value,
                            fit_config.pol_min, 
                            fit_config.pol_max);
    
    config.SetPolarisation(pol_param);
    // config.SetPolarisation("Pol[0.35,0.34,0.38]");
    cout << "  Polarization: " << pol_param << endl;
}

void configure_amplitude_model(PhotoTwoSpin0Amps& config, const MomentFitConfig& fit_config) {
    cout << "\n=== CONFIGURING AMPLITUDE MODEL ===" << endl;
    
    config.SetLmax(fit_config.lmax);
    config.SetMmax(fit_config.mmax); 
    config.SetNrefl(fit_config.nrefl);
    
    cout << "  Multipoles: L_max=" << fit_config.lmax 
         << ", M_max=" << fit_config.mmax << ", N_refl=" << fit_config.nrefl << endl;
    
    // Configure the model
    config.ConfigureMoments();
    config.LoadModelPDF();
    
    cout << "  Model configured and PDF loaded" << endl;
    cout << "==================================" << endl;
}

void load_data_files(FitManager& fitter, PhotoTwoSpin0Amps& amplitude_config, const MomentFitConfig& config) {
    cout << "\n=== LOADING DATA FILES ===" << endl;
    
    // Set UID branch for BruFit identification
    fitter.SetUp().SetIDBranchName("UID");
    
    TString tree = config.tree_name.IsNull() ? "kin" : config.tree_name;

    // Load generated MC (acceptance correction) - optional
    if (config.use_generated_mc) {
        cout << "  Loading generated MC (acceptance correction)..." << endl;
        fitter.LoadGenerated(tree, get_full_path(config.input_data_dir, config.generated_file), amplitude_config.GetName());
    } else {
        cout << "  Skipping generated MC (no acceptance correction)" << endl;
    }
    
    // Load simulated MC (detector response)
    cout << "  Loading simulated MC..." << endl;
    fitter.LoadSimulated(tree, get_full_path(config.input_data_dir, config.simulated_file), amplitude_config.GetName());

    // Load real data
    cout << "  Loading data..." << endl;
    fitter.LoadData(tree, get_full_path(config.input_data_dir, config.data_file));

    // Load weights
    cout << "  Loading weights..." << endl;
    fitter.Data().LoadWeights(config.weight_species, get_full_path(config.input_data_dir, config.weights_file), config.weight_name);

    cout << "  All required data files loaded successfully" << endl;
    cout << "============================" << endl;

    Int_t n_mc_bins = config.mc_integration_bins > 0 ? config.mc_integration_bins : 50;
    if (config.use_gj_frame) {
        fitter.SetUp().GetVar("cosTheta_eta_gj")->setBins(n_mc_bins);
        fitter.SetUp().GetVar("phi_eta_gj")->setBins(n_mc_bins);
    } else {
        fitter.SetUp().GetVar("cosTheta_eta_hel")->setBins(n_mc_bins);
        fitter.SetUp().GetVar("phi_eta_hel")->setBins(n_mc_bins);
    }

    fitter.SetUp().GetVar("Phi")->setBins(n_mc_bins);
}

void setup_binning(FitManager& fitter, const MomentFitConfig& config) {
    cout << "\n=== SETTING UP BINNING ===" << endl;
    
    // Validate that all binning vectors have the same size
    size_t n_variables = config.bin_variables.size();
    if (config.n_bins.size() != n_variables || 
        config.bin_min.size() != n_variables || 
        config.bin_max.size() != n_variables) {
        cout << "ERROR: Binning configuration vectors must have the same size!" << endl;
        cout << "  Variables: " << config.bin_variables.size() << endl;
        cout << "  n_bins: " << config.n_bins.size() << endl;
        cout << "  bin_min: " << config.bin_min.size() << endl;
        cout << "  bin_max: " << config.bin_max.size() << endl;
        return;
    }
    
    // Load binning for each variable
    for (size_t i = 0; i < n_variables; ++i) {
        
        fitter.Bins().LoadBinVar(config.bin_variables[i], 
                                config.n_bins[i], 
                                config.bin_min[i], 
                                config.bin_max[i]);
        
        cout << "  Variable " << (i+1) << ": " << config.bin_variables[i] << endl;
        cout << "    Bins: " << config.n_bins[i] << endl;
        cout << "    Range: [" << config.bin_min[i] << ", " << config.bin_max[i] << "]" << endl;
        cout << "    Bin width: " << (config.bin_max[i] - config.bin_min[i]) / config.n_bins[i] << endl;
    }
    
    cout << "  Total dimensions: " << n_variables << endl;
    
    // Calculate total number of bins
    int total_bins = 1;
    for (size_t i = 0; i < n_variables; ++i) {
        total_bins *= config.n_bins[i];
    }
    cout << "  Total bins: " << total_bins << endl;
    cout << "=========================" << endl;
}

void configure_minimizer(FitManager& fitter, const MomentFitConfig& config) {
    cout << "\n=== CONFIGURING MINIMIZER ===" << endl;
    
    // Setup MCMC with covariance
    // auto mcmc = new BruMcmcCovariance(config.mcmc_iterations, 
    //                                  config.mcmc_burnin,
    //                                  config.mcmc_step,
    //                                  config.mcmc_target_acceptance,
    //                                  config.mcmc_min_acceptance, 
    //                                  config.mcmc_max_acceptance);

    // Create null MCMC pointer
    BruMcmcCovariance* mcmc = nullptr;
    if (!config.fixed_covariance){
        mcmc = new BruMcmcCovariance(config.mcmc_iterations, 
                                    config.mcmc_burnin,
                                    config.mcmc_step,
                                    config.mcmc_target_acceptance,
                                    config.mcmc_min_acceptance, 
                                    config.mcmc_max_acceptance);
        
    } else {
        // Use fixed covariance MCMC
        mcmc = new BruMcmcFixedCovariance(config.covariance_file,
                                          config.covariance_matrix_name,
                                          std::vector<Int_t>{config.mcmc_iterations[0], config.mcmc_iterations[2]},
                                          config.mcmc_burnin,
                                          config.mcmc_step,
                                          config.mcmc_target_acceptance,
                                          config.mcmc_min_acceptance,
                                          config.mcmc_max_acceptance);
        
        cout << "  Using fixed covariance matrix from file: " << config.covariance_file << endl;
        cout << "  Covariance matrix name: " << config.covariance_matrix_name << endl;
    }
                                
    if (!config.mcmc_use_covariance) {
        mcmc->TurnOffCovariance();
        cout << "  MCMC covariance matrix: OFF (using simple N-dimensional Gaussian proposal)" << endl;
    }

    if (config.mcmc_use_covariance_tuning) {
        if (config.mcmc_tune_mode == "Acceptance" || config.mcmc_tune_mode == "kAcceptance") {
            mcmc->SetTuningMode(HS::FIT::McmcTuneMode::kAcceptance);
        } else {
            mcmc->SetTuningMode(HS::FIT::McmcTuneMode::kMappedRhat);
        }
        mcmc->TuneCovarianceStep();
        mcmc->SetDiagnosticTuningParameters(config.mcmc_rhat_window, config.mcmc_rhat_target, config.mcmc_ess_fraction, config.mcmc_rhat_max_retries);
        cout << "  Activated MCMC covariance tuning (Mode: " 
             << config.mcmc_tune_mode 
             << ", Window: " << config.mcmc_rhat_window << ", Target R-hat: " << config.mcmc_rhat_target << ")" << endl;
    }

    if (config.mcmc_gibbs_block_size > 0) {
        mcmc->SetGibbsBlockSize(config.mcmc_gibbs_block_size);
        cout << "  Gibbs block size: " << config.mcmc_gibbs_block_size << endl;
    }

    if (config.mcmc_cold_warmup_acceptance > 0) {
        mcmc->SetColdWarmUpAcceptance(config.mcmc_cold_warmup_acceptance);
        cout << "  Cold warm-up acceptance: " << config.mcmc_cold_warmup_acceptance << endl;
    }
    
    if (auto* fixedMcmc = dynamic_cast<HS::FIT::BruMcmcFixedCovariance*>(mcmc)) {
        fixedMcmc->SetBailoutParameters(config.mcmc_bailout_threshold, config.mcmc_bailout_min_acceptance);
        fixedMcmc->BruMcmc::SetBailoutParameters(config.mcmc_bailout_threshold, config.mcmc_bailout_min_acceptance);
        cout << "  Bailout parameters: threshold=" << config.mcmc_bailout_threshold 
             << ", min_acceptance=" << config.mcmc_bailout_min_acceptance << endl;
    }
    
    fitter.SetMinimiser(mcmc);
    
    // Add fit options
    fitter.SetUp().AddFitOption(RooFit::NumCPU(config.num_cpu));
    if (config.set_init_value) {
        auto* fixedMcmc = dynamic_cast<BruMcmcFixedCovariance*>(mcmc);
        if (config.fixed_covariance && fixedMcmc) {
            cout << "  Initializing per-bin moments via BruMcmcFixedCovariance from: " << config.init_value_file << endl;
            fixedMcmc->SetInitValues(config.init_value_file, config.init_value_scale);
        } else {
            cout << "  Initializing the moments value from Linear Algebra method calculation..." << endl;
            cout << "  Opening file: " << config.init_value_file << endl;
            TFile* file = TFile::Open(config.init_value_file);
            if (!file || file->IsZombie()) {
                cerr << "Error opening file: " << config.init_value_file << endl;
                return;
            }

            TObjString* paramNamesObj = dynamic_cast<TObjString*>(file->Get("paramNames"));
            if (!paramNamesObj) {
                cerr << "  Error: paramNames not found in file." << endl;
                file->Close();
                return;
            }

            cout << "  Parameter names from the file:" << endl;
            vector<string> paramNames;
            stringstream namesStream(paramNamesObj->GetString().Data());
            string paramName;
            while (getline(namesStream, paramName, ',')) {
                if (!paramName.empty()) {
                    paramNames.push_back(paramName);
                    cout << paramNames.back() << endl;
                }
            }

            for (auto name: paramNames){
                if ((name != "Yld_Moments")) {
                    TString hist_name = "h1_norm_" + TString(name);
                    TH1* h1 = dynamic_cast<TH1*>(file->Get(hist_name.Data()));
                    if (!h1) {
                        cerr << "  Error: histogram not found or invalid type: " << hist_name << endl;
                        continue;
                    }
                    double moment_value = h1->GetBinContent(1);
                    moment_value *= config.init_value_scale;
                    cout << "  Setting " << name << " = " << moment_value << endl;
                    fitter.SetUp().SetParVal(name, moment_value, kFALSE);
                }
            }
        }
    }
    
    cout << endl;
    cout << "  Minimizer: BruMcmcCovariance" << endl;
    cout << "  Iterations: ";
    for (size_t i = 0; i < config.mcmc_iterations.size(); ++i) {
        cout << config.mcmc_iterations[i];
        if (i < config.mcmc_iterations.size() - 1) cout << ", ";
    }
    cout << endl;
    cout << "  Burnin: " << config.mcmc_burnin << endl;
    cout << "  Step: " << config.mcmc_step << endl;
    cout << "  Target acceptance: " << config.mcmc_target_acceptance << endl;
    cout << "  Acceptance range: [" << config.mcmc_min_acceptance 
         << ", " << config.mcmc_max_acceptance << "]" << endl;
    cout << "  CPUs: " << config.num_cpu << endl;
    cout << "==============================" << endl;
}

void run_fit(FitManager& fitter, const MomentFitConfig& config) {
    if (!config.run_fit) {
        cout << "\n=== FIT EXECUTION SKIPPED ===" << endl;
        cout << "  run_fit is set to false. No fit will be executed." << endl;
        cout << "  You can set run_fit to true in the configuration to execute the fit." << endl;
        cout << "=============================" << endl;
        return;
    }

    cout << "\n=== STARTING Moment FIT ===" << endl;
    
    if (config.use_multi) {
        cout << "  Using multiprocessing " << endl;
        Multi::Go(&fitter);
    } else {
        cout << "  Using local processing" << endl;
        Here::Go(&fitter);
    }
    
    cout << "  Fit completed!" << endl;
    cout << "===========================" << endl;
}

// Main function with configurable parameters
void MomentFitUtils(TString analysis_tag = "",
                    TString input_data_dir = "",
                    vector<TString> bin_variables = {"Mpi0eta"},
                    vector<int> n_bins = {30},
                    vector<double> bin_min = {0.8},
                    vector<double> bin_max = {2.0},
                    bool use_generated_mc = false,
                    bool use_Multi = true,
                    int Multi_workers = 30) {
    
    cout << "\n" << string(60, '=') << endl;
    cout << "      Moment Fit - EtaPi0->4Gamma" << endl;  
    cout << string(60, '=') << endl;
    
    // Setup configuration
    MomentFitConfig config;
    // config.analysis_tag = analysis_tag;
    // config.run_period = run_period;
    // config.bin_variables = bin_variables;
    // config.n_bins = n_bins;
    // config.bin_min = bin_min;
    // config.bin_max = bin_max;
    // config.use_generated_mc = use_generated_mc;
    // config.use_Multi = use_Multi;
    // config.Multi_workers = Multi_workers;
    
    // Auto-detect input directories if not provided
    if (input_data_dir != "") {
        config.input_data_dir = input_data_dir;
    }
    
    // if (input_data_dir == "") {
    //     config.input_data_dir = "/d/home/septian/EtaPi0Analysis/test_splot_output/eta_splot/";
    // } else {
    //     config.input_data_dir = input_data_dir;
    // }
    
    // Set input file paths
    // TString file_prefix = Form("polALL_t010100_m080200_%s", analysis_tag.Data());
    
    // Generated MC file (optional - only if use_generated_mc is true)
    // if (use_generated_mc) {
    //     config.generated_file = "/d/home/septian/EtaPi0Analysis/run_merged/rootFiles/t010100_m080200_nominal_MomentFit_2019_11/polALL_t010100_m080200_nominal_MomentFit_2019_11_FTOT_generated_gen_amp_V2_t010100_data_flat.root";
    // }
    
    // Use processed files from our workflow
    // config.simulated_file = config.input_data_dir + file_prefix + "_FTOT_selected_WideMass_acc_flat.root";
    // config.data_file = config.input_data_dir + file_prefix + "_DTOT_selected_WideMass_acc_flat.root";
    
    // weights file
    // config.weights_file = config.input_data_dir + "Mpi0eta1.400000_/ResultsHSMinuit2.root";
    
    // Setup output directory
    // setup_output_directory(config);
    
    // Print configuration
    print_config(config);
    
    // Check input files
    if (!check_input_files(config)) {
        cout << "\nERROR: Missing input files. Please check the paths above." << endl;
        return;
    }
    
    // Initialize BruFit components
    cout << "\n=== INITIALIZING BRUFIT ===" << endl;
    FitManager fitter;
    fitter.SetUp().SetOutDir(config.output_dir);
    
    // Setup amplitude configuration
    PhotoTwoSpin0Amps amplitude_config("Moments");
    amplitude_config.SetManager(&fitter);
    
    // Configure the model
    configure_amplitude_model(amplitude_config, config);

    // Setup binning
    setup_binning(fitter, config);
    
    // Load data files
    load_data_files(fitter, amplitude_config, config);
    
    // Configure minimizer
    configure_minimizer(fitter, config);
    
    // Run the fit
    if (config.run_fit) run_fit(fitter, config);

}
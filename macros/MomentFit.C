#if defined(__CLING__)
R__ADD_INCLUDE_PATH(/d/home/septian/Moments2Amplitudes/macros)
R__ADD_INCLUDE_PATH(/d/home/septian/EtaPi0Analysis/macros)
R__ADD_INCLUDE_PATH(/d/grid17/septian/EtaPi0Analysis/macros)
#endif

#include "MomentFitUtils.C"
#include <TSystem.h>
#include <cstdlib>
#include <map>

// ==============================================================================
// MomentFit.C
//
// Generalized BruFit MomentFit Macro supporting:
//   1. Full mass range fits (all bins in parallel via Multi::Go)
//   2. Single mass bin fits / retries (synchronous via Here::Go)
//   3. Flexible mass ranges [m_min, m_max] and arbitrary bin count n_bins
//   4. Configurable MCMC iterations, burn-in, step size, and initial scale
//   5. Explicit mapping to linear algebra covariance & initial value files
// ==============================================================================

void MomentFit(
    int tbin_index = 0,
    int maxL = 4,
    double init_value_scale = 1.0,
    double mcmc_step = 1.0,
    int target_bin_index = -1,      // -1: fit all bins; >=0: fit only this specific mass bin
    int mcmc_iterations = 10000,
    int mcmc_burnin = 2000,
    double m_min = 0.8,
    double m_max = 2.0,
    int n_bins = 30,
    TString custom_run_dir = "",
    TString custom_cov_dir = "",
    TString custom_data_dir = "",
    bool apply_mc_variance = true,
    bool allow_negative_pdf = true,
    bool use_tuning = true,
    int gibbs_block_size = 0,
    double cold_warmup_acc = 0.15,
    TString custom_data_file = "",
    TString custom_mcrecon_file = "",
    TString custom_mcgen_file = "",
    TString custom_weight_file = "",
    bool run_fit_flag = true
) {
    if (tbin_index < 0 || tbin_index > 4) {
        cout << "Invalid tbin_index. Please provide a value between 0 and 4:" << endl;
        cout << "  0: t010020" << endl;
        cout << "  1: t020032" << endl;
        cout << "  2: t032050" << endl;
        cout << "  3: t050075" << endl;
        cout << "  4: t075100" << endl;
        return;
    }

    map<int, TString> tbins = {
        {0, "t010020"},
        {1, "t020032"},
        {2, "t032050"},
        {3, "t050075"},
        {4, "t075100"}
    };
    TString tbin_name = tbins[tbin_index];

    MomentFitConfig config;

    int tag_min = (int)TMath::Nint(m_min * 100);
    int tag_max = (int)TMath::Nint(m_max * 100);
    config.analysis_tag = Form("m%03d%03d_maxL%d", tag_min, tag_max, maxL);

    // 1. Resolve input data directory
    TString input_dir = custom_data_dir;
    if (input_dir.IsNull()) {
        input_dir = "/d/grid17/septian/EtaPi0Analysis/data/nominal_tbins/";
        if (gSystem->AccessPathName(input_dir + Form("data_merged_%s.root", tbin_name.Data()))) {
            input_dir = "/d/grid17/septian/EtaPi0Analysis/data/nominal/";
        }
    }
    if (!input_dir.EndsWith("/")) input_dir += "/";
    config.input_data_dir = input_dir;

    // 2. Resolve covariance matrix & initial value directory
    TString moments_base_dir = custom_cov_dir;
    if (moments_base_dir.IsNull()) {
        moments_base_dir = Form("/d/home/septian/Moments2Amplitudes/Moments.out/nominal/%s/All.maxL_%d",
                                tbin_name.Data(), maxL);
    } else if (!moments_base_dir.Contains(tbin_name)) {
        if (!moments_base_dir.EndsWith("/")) moments_base_dir += "/";
        moments_base_dir = Form("%s%s/All.maxL_%d", moments_base_dir.Data(), tbin_name.Data(), maxL);
    }
    if (moments_base_dir.EndsWith("/")) moments_base_dir.Remove(moments_base_dir.Length() - 1);

    // 3. Resolve output directory
    TString base_output_dir = custom_run_dir;
    if (base_output_dir.IsNull()) {
        base_output_dir = Form("/d/grid17/septian/EtaPi0Analysis/run/MomentFit_m%03d%03d_FixedCov_Lmax%d_InitValueScale%.1f_step%0.1f/",
                               tag_min, tag_max, maxL, init_value_scale, mcmc_step);
    }
    if (!base_output_dir.EndsWith("/")) base_output_dir += "/";

    if (base_output_dir.EndsWith(tbin_name + "/")) {
        config.output_dir = base_output_dir;
    } else {
        config.output_dir = base_output_dir + tbin_name + "/";
    }

    // 4. Configure mass binning
    config.bin_variables = {"Mpi0eta"};
    double full_bin_width = (m_max - m_min) / n_bins;

    if (target_bin_index >= 0) {
        if (target_bin_index >= n_bins) {
            cerr << "ERROR: target_bin_index (" << target_bin_index << ") out of range [0, " << (n_bins - 1) << "]" << endl;
            return;
        }
        double bin_low = m_min + target_bin_index * full_bin_width;
        double bin_high = m_min + (target_bin_index + 1) * full_bin_width;
        double bin_center = (bin_low + bin_high) / 2.0;

        config.n_bins = {1};
        config.bin_min = {bin_low};
        config.bin_max = {bin_high};
        config.use_multi = false; // Run single bin synchronously via Here::Go

        cout << "==========================================================" << endl;
        cout << Form(" Mode: SINGLE MASS BIN REFIT (Bin index: %d)", target_bin_index) << endl;
        cout << Form(" Mass Range: [%.4f, %.4f] GeV (Center: %.4f GeV)", bin_low, bin_high, bin_center) << endl;
        cout << "==========================================================" << endl;
    } else {
        config.n_bins = {n_bins};
        config.bin_min = {m_min};
        config.bin_max = {m_max};
        config.use_multi = true;  // Run all bins in parallel via Multi::Go

        cout << "==========================================================" << endl;
        cout << Form(" Mode: ALL MASS BINS FIT (%d bins)", n_bins) << endl;
        cout << Form(" Mass Range: [%.4f, %.4f] GeV (Bin width: %.4f GeV)", m_min, m_max, full_bin_width) << endl;
        cout << "==========================================================" << endl;
    }

    auto env_double = [](const char* name, double fallback) -> double {
        const char* val = gSystem->Getenv(name);
        return (val && val[0] != '\0') ? atof(val) : fallback;
    };
    auto env_int = [](const char* name, int fallback) -> int {
        const char* val = gSystem->Getenv(name);
        return (val && val[0] != '\0') ? atoi(val) : fallback;
    };
    auto env_bool = [](const char* name, bool fallback) -> bool {
        const char* val = gSystem->Getenv(name);
        if (!val || val[0] == '\0') return fallback;
        TString s(val);
        s.ToLower();
        return (s == "1" || s == "true" || s == "yes" || s == "on");
    };
    auto env_str = [](const char* name, const char* fallback) -> TString {
        const char* val = gSystem->Getenv(name);
        return (val && val[0] != '\0') ? TString(val) : TString(fallback);
    };

    // 5. Multipoles and physics configuration
    config.lmax = TMath::Floor(maxL / 2);
    config.mmax = TMath::Floor(maxL / 2);
    config.nrefl = 2;
    config.num_cpu = env_int("MOM2AMP_NUM_CPU", 1);
    config.run_fit = run_fit_flag;
    config.use_gj_frame = env_bool("MOM2AMP_USE_GJ_FRAME", true);

    // 6. MCMC Configuration
    config.mcmc_iterations = {mcmc_iterations, mcmc_iterations, mcmc_iterations};
    config.mcmc_burnin = mcmc_burnin;
    config.mcmc_step = mcmc_step;
    config.mcmc_min_acceptance = env_double("MOM2AMP_MIN_ACCEPTANCE", 0.05);
    config.mcmc_max_acceptance = env_double("MOM2AMP_MAX_ACCEPTANCE", 0.15);
    double default_target_acc = 0.5 * (config.mcmc_min_acceptance + config.mcmc_max_acceptance);
    config.mcmc_target_acceptance = env_double("MOM2AMP_TARGET_ACCEPTANCE", default_target_acc);

    // Likelihood & PDF options
    config.apply_mc_variance = apply_mc_variance;
    config.allow_negative_pdf = allow_negative_pdf;
    config.mc_integration_bins = env_int("MOM2AMP_MC_INT_BINS", 50);

    // Tuning options
    config.mcmc_use_covariance_tuning = use_tuning;
    config.mcmc_tune_mode = env_str("MOM2AMP_TUNE_MODE", "MappedRhat");
    config.mcmc_rhat_window = env_int("MOM2AMP_RHAT_WINDOW", 1000);
    config.mcmc_rhat_target = env_double("MOM2AMP_RHAT_TARGET", 1.05);
    config.mcmc_ess_fraction = env_double("MOM2AMP_ESS_FRACTION", 0.02);
    config.mcmc_rhat_max_retries = env_int("MOM2AMP_RHAT_MAX_RETRIES", 5);
    config.mcmc_gibbs_block_size = gibbs_block_size;
    config.mcmc_cold_warmup_acceptance = cold_warmup_acc;

    // Fast-fail bailout: cancel if after bailout_threshold sampling proposals acceptance falls short
    config.mcmc_bailout_threshold = env_int("MOM2AMP_BAILOUT_THRESHOLD", 10000);
    config.mcmc_bailout_min_acceptance = env_double("MOM2AMP_BAILOUT_MIN_ACC", 0.0001);

    // 7. Fixed covariance & Initial values
    config.fixed_covariance = env_bool("MOM2AMP_FIXED_COV", true);
    config.mcmc_use_covariance = env_bool("MOM2AMP_USE_COV", true);
    config.set_init_value = env_bool("MOM2AMP_SET_INIT_VALUE", true);
    config.init_value_scale = init_value_scale;
    config.covariance_matrix_name = "MyCovMatrix";
    int rng_seed = env_int("MOM2AMP_SEED", 111);

    if (target_bin_index >= 0) {
        double bin_low = m_min + target_bin_index * full_bin_width;
        double bin_high = m_min + (target_bin_index + 1) * full_bin_width;
        double center = (bin_low + bin_high) / 2.0;
        config.covariance_file = Form("%s/BruCovMatrix_Mpi0eta_%.3f.root", moments_base_dir.Data(), center);
        config.init_value_file = Form("%s/covariance_matrix_Mpi0eta_%.3f.root", moments_base_dir.Data(), center);
    } else {
        config.covariance_file = moments_base_dir + "/BruCovMatrix_Mpi0eta_%.3f.root";
        config.init_value_file = moments_base_dir + "/covariance_matrix_Mpi0eta_%.3f.root";
    }

    // 8. Data files
    auto resolve_data_file = [&](TString pattern, const char* default_pattern) -> TString {
        if (pattern.IsNull()) {
            return Form(default_pattern, tbin_name.Data());
        }
        pattern.ReplaceAll("{tbin}", tbin_name);
        pattern.ReplaceAll("{data_dir}", input_dir);
        pattern.ReplaceAll("{base_dir}", input_dir);
        if (pattern.Contains("%s")) {
            pattern = Form(pattern.Data(), tbin_name.Data());
        }
        return pattern;
    };

    config.data_file = resolve_data_file(custom_data_file, "data_merged_%s.root");
    config.simulated_file = resolve_data_file(custom_mcrecon_file, "reconmc_merged_%s.root");
    config.generated_file = resolve_data_file(custom_mcgen_file, "genmc_merged_%s.root");
    config.weights_file = resolve_data_file(custom_weight_file, "weights_data_merged_%s.root");

    print_config(config);
    cout << "  Fixed Covariance Directory: " << moments_base_dir << endl;
    cout << "  Covariance File Pattern:    " << config.covariance_file << endl;
    cout << "  Initial Values Pattern:     " << config.init_value_file << " (Scale: " << config.init_value_scale << ")" << endl;

    if (!check_input_files(config)) {
        cerr << "ERROR: Missing input data files. Please verify paths above." << endl;
        return;
    }

    // 9. BruFit setup
    FitManager fitter;
    gRandom->SetSeed(rng_seed);
    fitter.SetUp().SetOutDir(config.output_dir);
    fitter.SetUp().AddFitOption(RooFit::EvalBackend("cpu"));
    fitter.SetUp().ApplyMCVariance(config.apply_mc_variance);
    fitter.SetUp().SetAllowNegativePDF(config.allow_negative_pdf);
    fitter.SetUp().SetIDBranchName("UID");
    PhotoTwoSpin0Amps amplitude_config("Moments");
    amplitude_config.SetManager(&fitter);
    configure_variables(amplitude_config, config);
    setup_binning(fitter, config);
    load_data_files(fitter, amplitude_config, config);
    if (!config.run_fit) {
        cout << "\n=== BINNING ONLY MODE (run_fit=false) ===" << endl;
        cout << "Restored BinsConfig files in " << config.output_dir << "; skipping MCMC fit execution." << endl;
        cout << "=========================================" << endl;
        return;
    }
    configure_amplitude_model(amplitude_config, config);

    cout << "\n=== CONFIGURING MCMC MINIMIZER ===" << endl;
    cout << "  Fixed covariance:    " << (config.fixed_covariance ? "Enabled" : "Disabled") << endl;
    cout << "  Init values enabled: " << (config.set_init_value ? "Enabled" : "Disabled") << endl;
    cout << "  MC integration bins: " << config.mc_integration_bins << endl;
    cout << "  Use GJ frame:        " << (config.use_gj_frame ? "Yes" : "No (Helicity)") << endl;
    cout << "  RNG seed / NumCPU:   " << rng_seed << " / " << config.num_cpu << endl;

    HS::FIT::BruMcmcCovariance* mcmc = nullptr;
    if (config.fixed_covariance) {
        TString ctor_init_file = config.set_init_value ? config.init_value_file : TString("");
        auto* fixedMcmc = new HS::FIT::BruMcmcFixedCovariance(
            config.covariance_file,
            config.covariance_matrix_name,
            ctor_init_file,
            config.init_value_scale,
            std::vector<Int_t>{config.mcmc_iterations[0], config.mcmc_iterations[2]},
            config.mcmc_burnin,
            config.mcmc_step,
            config.mcmc_target_acceptance,
            config.mcmc_min_acceptance,
            config.mcmc_max_acceptance
        );

        if (target_bin_index >= 0) {
            double bin_low = m_min + target_bin_index * full_bin_width;
            double bin_high = m_min + (target_bin_index + 1) * full_bin_width;
            double center = (bin_low + bin_high) / 2.0;
            TString cov_file = Form("%s/BruCovMatrix_Mpi0eta_%.3f.root", moments_base_dir.Data(), center);
            TString init_file = Form("%s/covariance_matrix_Mpi0eta_%.3f.root", moments_base_dir.Data(), center);
            TString bName = fitter.Bins().BinName(0);

            fixedMcmc->SetBinFile(bName, cov_file);
            fixedMcmc->SetBinFile(0, cov_file);
            if (config.set_init_value) {
                fixedMcmc->SetBinInitFile(bName, init_file);
                fixedMcmc->SetBinInitFile(0, init_file);
            }
            cout << Form("  Mapped single bin: %s -> %s (init: %s)",
                         bName.Data(), cov_file.Data(),
                         config.set_init_value ? init_file.Data() : "disabled") << endl;
        } else {
            int fit_bins = config.n_bins[0];
            double b_width = (config.bin_max[0] - config.bin_min[0]) / fit_bins;
            for (int i = 0; i < fit_bins; ++i) {
                double center = config.bin_min[0] + (i + 0.5) * b_width;
                TString cov_file = Form("%s/BruCovMatrix_Mpi0eta_%.3f.root", moments_base_dir.Data(), center);
                TString init_file = Form("%s/covariance_matrix_Mpi0eta_%.3f.root", moments_base_dir.Data(), center);
                TString bName = fitter.Bins().BinName(i);

                fixedMcmc->SetBinFile(bName, cov_file);
                fixedMcmc->SetBinFile(i, cov_file);
                if (config.set_init_value) {
                    fixedMcmc->SetBinInitFile(bName, init_file);
                    fixedMcmc->SetBinInitFile(i, init_file);
                }
            }
            cout << "  Explicitly mapped " << fit_bins << " bins into BruMcmcFixedCovariance." << endl;
        }

        fixedMcmc->SetBailoutParameters(config.mcmc_bailout_threshold, config.mcmc_bailout_min_acceptance);
        fixedMcmc->BruMcmc::SetBailoutParameters(config.mcmc_bailout_threshold, config.mcmc_bailout_min_acceptance);
        mcmc = fixedMcmc;
    } else {
        mcmc = new HS::FIT::BruMcmcCovariance(
            config.mcmc_iterations,
            config.mcmc_burnin,
            config.mcmc_step,
            config.mcmc_target_acceptance,
            config.mcmc_min_acceptance,
            config.mcmc_max_acceptance
        );
        mcmc->SetBailoutParameters(config.mcmc_bailout_threshold, config.mcmc_bailout_min_acceptance);

        if (config.set_init_value) {
            TString init_file = config.init_value_file;
            if (init_file.Contains("%.3f")) {
                int mid_idx = config.n_bins[0] / 2;
                double mid_center = config.bin_min[0] + (mid_idx + 0.5) * ((config.bin_max[0] - config.bin_min[0]) / config.n_bins[0]);
                init_file = Form("%s/covariance_matrix_Mpi0eta_%.3f.root", moments_base_dir.Data(), mid_center);
            }
            TFile* ifile = TFile::Open(init_file);
            if (ifile && !ifile->IsZombie()) {
                TObjString* paramNamesObj = dynamic_cast<TObjString*>(ifile->Get("paramNames"));
                if (paramNamesObj) {
                    stringstream namesStream(paramNamesObj->GetString().Data());
                    string paramName;
                    while (getline(namesStream, paramName, ',')) {
                        if (paramName.empty() || paramName == "Yld_Moments") continue;
                        TString hist_name = "h1_norm_" + TString(paramName);
                        TH1* h1 = dynamic_cast<TH1*>(ifile->Get(hist_name.Data()));
                        if (h1) {
                            double val = h1->GetBinContent(1) * config.init_value_scale;
                            fitter.SetUp().SetParVal(paramName, val, kFALSE);
                        }
                    }
                    cout << "  Initialized moments from " << init_file << " (Scale: " << config.init_value_scale << ")" << endl;
                }
                ifile->Close();
            } else {
                cout << "  Notice: Init value file not opened (" << init_file << "); using default initial values." << endl;
                if (ifile) ifile->Close();
            }
        }
    }

    mcmc->SetDesiredAcceptance(config.mcmc_min_acceptance, config.mcmc_max_acceptance, config.mcmc_target_acceptance);

    if (config.mcmc_use_covariance_tuning) {
        if (config.mcmc_tune_mode == "Acceptance" || config.mcmc_tune_mode == "kAcceptance") {
            mcmc->SetTuningMode(HS::FIT::McmcTuneMode::kAcceptance);
        } else {
            mcmc->SetTuningMode(HS::FIT::McmcTuneMode::kMappedRhat);
        }
        mcmc->TuneCovarianceStep();
        mcmc->SetDiagnosticTuningParameters(config.mcmc_rhat_window, config.mcmc_rhat_target, config.mcmc_ess_fraction, config.mcmc_rhat_max_retries);
        cout << "  Activated MCMC covariance tuning (Mode: " << config.mcmc_tune_mode
             << ", Window: " << config.mcmc_rhat_window << ", Target R-hat: " << config.mcmc_rhat_target
             << ", ESS fraction: " << config.mcmc_ess_fraction << ", Max retries: " << config.mcmc_rhat_max_retries << ")" << endl;
    } else {
        mcmc->TuneCovarianceStep(kFALSE);
    }

    if (config.mcmc_gibbs_block_size > 0) {
        mcmc->SetGibbsBlockSize(config.mcmc_gibbs_block_size);
        cout << "  Gibbs block size: " << config.mcmc_gibbs_block_size << endl;
    }

    if (config.mcmc_cold_warmup_acceptance > 0) {
        mcmc->SetColdWarmUpAcceptance(config.mcmc_cold_warmup_acceptance);
        cout << "  Cold warm-up acceptance: " << config.mcmc_cold_warmup_acceptance << endl;
    }

    if (!config.mcmc_use_covariance) {
        mcmc->TurnOffCovariance();
        cout << "  MCMC covariance matrix: OFF" << endl;
    }

    cout << "  Bailout parameters: threshold=" << config.mcmc_bailout_threshold
         << ", min_acceptance=" << config.mcmc_bailout_min_acceptance << endl;
    fitter.SetMinimiser(mcmc);
    fitter.SetUp().AddFitOption(RooFit::NumCPU(config.num_cpu));
    cout << "===================================================" << endl;

    run_fit(fitter, config);
}

// Backward-compatibility wrapper for MomentFit_General
void MomentFit_General(
    int tbin_index = 0,
    int maxL = 4,
    double init_value_scale = 0.2,
    double mcmc_step = 0.1,
    int target_bin_index = -1,
    int mcmc_iterations = 6000,
    int mcmc_burnin = 1000,
    double m_min = 0.8,
    double m_max = 2.0,
    int n_bins = 30,
    TString custom_run_dir = "",
    TString custom_cov_dir = "",
    TString custom_data_dir = "",
    bool apply_mc_variance = true,
    bool allow_negative_pdf = true,
    bool use_tuning = true,
    int gibbs_block_size = 0,
    double cold_warmup_acc = 0.15,
    TString custom_data_file = "",
    TString custom_mcrecon_file = "",
    TString custom_mcgen_file = "",
    TString custom_weight_file = "",
    bool run_fit_flag = true
) {
    MomentFit(tbin_index, maxL, init_value_scale, mcmc_step, target_bin_index,
              mcmc_iterations, mcmc_burnin, m_min, m_max, n_bins,
              custom_run_dir, custom_cov_dir, custom_data_dir,
              apply_mc_variance, allow_negative_pdf, use_tuning,
              gibbs_block_size, cold_warmup_acc,
              custom_data_file, custom_mcrecon_file, custom_mcgen_file, custom_weight_file,
              run_fit_flag);
}

// Backward-compatibility wrapper for MomentFit_m080200
void MomentFit_m080200(
    int tbin_index = -1,
    int maxL = 4,
    double init_value_scale = 0.5,
    double mcmc_step = 0.8,
    bool apply_mc_variance = true,
    bool allow_negative_pdf = true,
    bool use_tuning = true,
    int gibbs_block_size = 0,
    double cold_warmup_acc = 0.15
) {
    MomentFit(tbin_index, maxL, init_value_scale, mcmc_step, -1, 6000, 1000, 0.8, 2.0, 30,
              "", "", "", apply_mc_variance, allow_negative_pdf, use_tuning,
              gibbs_block_size, cold_warmup_acc);
}

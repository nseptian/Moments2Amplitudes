
#if defined(__CLING__)
R__ADD_INCLUDE_PATH(/d/home/septian/Moments2Amplitudes/macros)
R__ADD_INCLUDE_PATH(/d/home/septian/EtaPi0Analysis/macros)
R__ADD_INCLUDE_PATH(/d/grid17/septian/EtaPi0Analysis/macros)
#endif
#include "CrossSectionUtils.C"
#include <TSystem.h>

void CalculateAcceptance(
    int tbin_index = 0,
    TString custom_run_dir = "",
    TString result_file_name = "",
    int maxL = 4,
    double m_min = 0.8,
    double m_max = 2.0,
    int n_bins = 30,
    TString custom_data_dir = "",
    bool force_rebin = false
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

    CrossSectionConfig config;

    // Resolve fit result directory
    TString fit_dir = custom_run_dir;
    if (fit_dir.IsNull()) {
        fit_dir = Form("/d/home/septian/Moments2Amplitudes/run/MomentFit_m%03d%03d_FixedCov_Lmax%d_InitValueScale1.0_step1.0/%s/",
                       (int)TMath::Nint(m_min * 100), (int)TMath::Nint(m_max * 100), maxL, tbin_name.Data());
    } else {
        if (!fit_dir.EndsWith("/")) fit_dir += "/";
        if (!fit_dir.Contains(tbin_name)) {
            fit_dir = fit_dir + tbin_name + "/";
        }
    }
    config.fit_result_dir = fit_dir;

    // Auto-detect result ROOT filename if not provided
    if (!result_file_name.IsNull()) {
        config.result_file_name = result_file_name;
    } else {
        double first_bin_center = m_min + 0.5 * (m_max - m_min) / n_bins;
        TString first_fixed_cov = config.fit_result_dir + Form("Mpi0eta%.6f_/ResultsBruMcmcFixedCovariance.root", first_bin_center);
        if (!gSystem->AccessPathName(first_fixed_cov) ||
            config.fit_result_dir.Contains("FixedCov")) {
            config.result_file_name = "ResultsBruMcmcFixedCovariance.root";
        } else {
            config.result_file_name = "ResultsBruMcmcCovariance.root";
        }
    }

    config.analysis_tag = Form("m%03d%03d_maxL%d", (int)TMath::Nint(m_min * 100), (int)TMath::Nint(m_max * 100), maxL);
    config.lmax = TMath::Floor(maxL / 2);
    config.mmax = TMath::Floor(maxL / 2);
    config.nrefl = 2;

    // Data directory with merged weights and binner inputs (prioritize nominal_tbins)
    TString weights_fname = Form("weights_data_merged_%s.root", tbin_name.Data());
    TString input_dir;
    if (!custom_data_dir.IsNull()) {
        input_dir = EnsureTrailingSlash(custom_data_dir);
    } else {
        std::vector<TString> candidate_dirs = {
            "/d/home/septian/Moments2Amplitudes/data/data/nominal_tbins/",
            "/d/grid17/septian/EtaPi0Analysis/data/nominal_tbins/"
        };
        for (const auto& cand : candidate_dirs) {
            if (!gSystem->AccessPathName(cand + weights_fname)) {
                input_dir = cand;
                break;
            }
        }
        if (input_dir.IsNull()) input_dir = candidate_dirs.front();
    }
    config.input_data_dir = input_dir;
    config.data_file = Form("data_merged_%s.root", tbin_name.Data());
    config.simulated_file = Form("reconmc_merged_%s.root", tbin_name.Data());
    config.generated_file = Form("genmc_merged_%s.root", tbin_name.Data());
    config.weights_file = weights_fname;
    config.force_rebin = force_rebin;
    config.weight_species = "Acc";
    config.weight_name = "weight";

    config.bin_variables = {"Mpi0eta"};
    config.n_bins = {n_bins};
    config.bin_min = {m_min};
    config.bin_max = {m_max};
    config.luminosity_file = "/d/home/septian/Moments2Amplitudes/data/flux_GlueXI_GlueX2020.root";

    cout << "==========================================================" << endl;
    cout << " Calculating Acceptance for " << tbin_name << " (maxL=" << maxL << ")" << endl;
    cout << " Fit Result Directory: " << config.fit_result_dir << endl;
    cout << " Result File:          " << config.result_file_name << endl;
    cout << " Input Data Directory: " << config.input_data_dir << endl;
    cout << " Weights File:         " << config.weights_file << endl;
    cout << " Luminosity File:      " << config.luminosity_file << endl;
    cout << "==========================================================" << endl;

    CalculateCrossSection(config);
    SaveAcceptance(config);
}

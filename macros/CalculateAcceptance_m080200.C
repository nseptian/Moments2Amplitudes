
#if defined(__CLING__)
R__ADD_INCLUDE_PATH(/d/home/septian/EtaPi0Analysis/macros)
R__ADD_INCLUDE_PATH(/d/grid17/septian/EtaPi0Analysis/macros)
R__ADD_INCLUDE_PATH(/d/home/septian/Moments2Amplitudes/macros)
#endif
#include "CrossSectionUtils.C"
#include <TSystem.h>

void CalculateAcceptance_m080200(
    int tbin_index = 0,
    TString custom_run_dir = "",
    TString result_file_name = "",
    int maxL = 4
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
        fit_dir = "/d/grid17/septian/EtaPi0Analysis/run/MomentFit_m080200/" + tbin_name + "/";
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
        if (!gSystem->AccessPathName(config.fit_result_dir + "Mpi0eta0.820000_/ResultsBruMcmcFixedCovariance.root") ||
            !gSystem->AccessPathName(config.fit_result_dir + "Mpi0eta0.940000_/ResultsBruMcmcFixedCovariance.root")) {
            config.result_file_name = "ResultsBruMcmcFixedCovariance.root";
        } else {
            config.result_file_name = "ResultsBruMcmcCovariance.root";
        }
    }

    config.analysis_tag = Form("m080200_maxL%d", maxL);
    config.lmax = TMath::Floor(maxL / 2);
    config.mmax = TMath::Floor(maxL / 2);
    config.nrefl = 2;

    // Data directory with merged weights
    TString input_dir = "/d/grid17/septian/EtaPi0Analysis/data/m080200/";
    if (gSystem->AccessPathName(input_dir + Form("weights_data_merged_%s.root", tbin_name.Data()))) {
        input_dir = "/d/grid17/septian/EtaPi0Analysis/data/nominal_tbins/";
    }
    config.input_data_dir = input_dir;
    config.weights_file = Form("weights_data_merged_%s.root", tbin_name.Data());
    config.weight_species = "Acc";
    config.weight_name = "weight";

    config.bin_variables = {"Mpi0eta"};
    config.n_bins = {30};
    config.bin_min = {0.8};
    config.bin_max = {2.0};
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
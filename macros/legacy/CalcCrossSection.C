

void CalcCrossSection(){

    CrossSection RF;

    // TString fitdir = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m080200_MCMCN6000BI1000S08WCOV_R6.34/";
    // TString fitdir = "/d/home/septian/EtaPi0Analysis/run_Phase2/fitMoment_2019_11_polALL_sPlot_t010100_m080200_MCMC_R6.34/";
    // TString weightsDir = "/d/home/septian/EtaPi0Analysis/run_Phase2/rootFiles/";

    // TString fitdir = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase2_MomentsFit_t010020_Mpi0eta080200_800000_mandelstam_t010020_100200_MCMC_R6.34/";
    // TString weightsDir = "/d/home/septian/EtaPi0Analysis/test_splot_output/eta_splot/";
    // TString weightFileName = "WeightsMpi0eta1.400000_.root";

    // TString fitdir = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase1_Phase2_MomentsFit_t010030_Mpi0eta080200_800000_mandelstam_t010030_100300_MCMC_R6.34/";
    TString fitdir = "/d/home/septian/EtaPi0Analysis/fitMoment_Phase1_Phase2_MomentsFit_nominal_Mpi0eta080200_800000_mandelstam_t010020_100200_MCMC_R6.34/";
    TString weightsDir = "/d/home/septian/EtaPi0Analysis/data/nominal/processed/";
    TString weightFileName = "weights_data_phase1_phase2_merged.root";

    RF.SetUp().SetOutDir(fitdir);
    RF.SetResultDir(fitdir);
    RF.SetResultFileName("ResultsBruMcmcCovariance.root");

    PhotoTwoSpin0Amps config("Moments");

    config.SetManager(&RF);
    
    config.SetDecayAngleCosTh("cosTheta_eta_gj[0.21,-1,1]");
    config.SetDecayAnglePhi("phi_eta_gj[0.0,-3.14159,3.14159]");
    config.SetPolPhi("Phi[0.2,-3.14159,3.14159]");

    config.SetPolarisation("Pol[0.35,0.34,0.38]");

    RF.SetUp().SetIDBranchName("UID");
    RF.Bins().LoadBinVar("Mpi0eta",30,0.8,2.0);
    RF.Bins().LoadBinVar("mandelstam_t",1,0.1,0.2);

    // RF.Data().LoadWeights("Acc",weightsDir + "weights_polALL_t010100_m080200_nominal_MomentsFit_PhaseI_2019_11_DTOT_selected_nominal_MomentsFit_acc_flat.root", "weight");
    // RF.Data().LoadWeights("EtaSignal",weightsDir + weightFileName,"HSsWeights");
    RF.Data().LoadWeights("Acc",weightsDir + weightFileName,"weight");
    RF.ReloadData(fitdir);
    RF.ReloadSimulated(fitdir,"Moments");
    RF.ReloadGenerated(fitdir,"Moments");

    config.SetLmax(2);
    config.SetMmax(2);
    config.SetNrefl(2);

    config.ConfigureMoments();
    config.LoadModelPDF();

    vector<Double_t > energybin = {8.0,8.8};
    RF.SetBeamEnergyBinLimits(energybin);

    // eta -> 2gamma BF = 0.3936
    // pi0 -> 2gamma BF = 0.98823
    double BR = 0.3936 * 0.98823;
    RF.SetBranchingRatio(BR);
    RF.SetTargetThickness(1);
    RF.SetFlux("/d/home/septian/Moments2Amplitudes/data/flux_GlueXI_GlueX2020.root","tagged_lumi");
    // RF.SetFlux("/d/home/septian/Moments2Amplitudes/data/flux_71350_73266.root","tagged_lumi");

    Here::Go(&RF);

    // RF.DrawResults(fitdir + "/CrossSection.root");
}
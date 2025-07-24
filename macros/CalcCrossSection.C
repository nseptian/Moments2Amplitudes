

void CalcCrossSection(){

    CrossSection RF;

    TString fitdir = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m080200_MCMCN6000BI1000S08WCOV_R6.34/";
    TString weightsDir = "/d/home/septian/EtaPi0Analysis/run_merged/rootFiles/";

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

    RF.Data().LoadWeights("Acc",weightsDir + "weights_polALL_t010100_m080200_nominal_MomentsFit_PhaseI_2019_11_DTOT_selected_nominal_MomentsFit_acc_flat.root", "weight");
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
    RF.SetFlux("/d/home/septian/Moments2Amplitudes/brufit/flux_GlueXI_GlueX2020.root","tagged_lumi");

    Here::Go(&RF);

    // RF.DrawResults(fitdir + "/CrossSection.root");
}
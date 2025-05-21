
void DoMomentFit(TString tag, TString inputData, TString inputSimulated, TString inputWeight, TString outDir)
{
    FitManager Fitter;// manage the fitting
    //set the output directory for the fit results files Results*.root
    Fitter.SetUp().SetOutDir(outDir);

    //Use amlitude configue class to define model
    PhotoTwoSpin0Amps config(tag.Data());
    config.SetManager(&Fitter);
    //set data variables which must be in the input tree
    // config.SetDecayAngleCosTh("cosTheta[0.21,-1,1]");
    // config.SetDecayAnglePhi("phi[0.1,-3.14159,3.14159]");
    // config.SetPolPhi("Phi[0.2,-3.14159,3.14159]");
    // config.SetConstPolarisation("Pol[1.0]");

    config.SetDecayAngleCosTh("cosTheta_eta_hel[0.21,-1,1]");
    config.SetDecayAnglePhi("phi_eta_hel[0.1,-3.14159,3.14159]");
    config.SetPolPhi("Phi[0.2,-3.14159,3.14159]");
    config.SetConstPolarisation("Pol[0.3]");

    // Fitter.SetUp().SetIDBranchName("UID");
    // Fitter.Bins().LoadBinVar("Mpi0eta",17,1.04,1.72);
    // Fitter.LoadSimulated("data","/d/home/septian/EtaPi0Analysis/run/fitBruMoments/plotsTestNizar/efficiencyReco.photoProd.root",config.GetName());
    // Fitter.LoadData("data","/d/home/septian/EtaPi0Analysis/run/fitBruMoments/plotsTestNizar/intensity.photoProd.root");
    // Fitter.Data().LoadWeights("Acc","/d/grid17/septian/EtaPi0Analysis/run/rootFiles/t010020_m104172_Phase1_brufit/weights_polALL_t010020_m104172_Phase1_brufit_DTOT_selected_nominal_wPhotonSyst_acc_flat.root","weight");

    // test fit with GlueX MC with GlueX acceptance
    // Fitter.LoadSimulated("kin","/d/grid17/septian/EtaPi0Analysis/MomentMCStudy/accepted.root",config.GetName());
    // Fitter.LoadData("kin","/d/grid17/septian/EtaPi0Analysis/MomentMCStudy/data.root");

    // test fit with GlueX MC with signal and background
    // Fitter.LoadData("kin","/d/grid17/septian/EtaPi0Analysis/MomentMCStudy/data_sig_bkg.root");
    // Fitter.LoadData("kin","/d/grid17/septian/EtaPi0Analysis/MomentMCStudy/data_sig_bkg_pol0.5.root");
    // Fitter.LoadSimulated("kin","/d/grid17/septian/EtaPi0Analysis/MomentMCStudy/accepted.root",config.GetName());
    // Fitter.Data().LoadWeights("Acc","/d/grid17/septian/EtaPi0Analysis/MomentMCStudy/weights_data_sig_bkg_pol0.5.root","weight");

    // test fit with GlueX MC with signal and background and seed variation
    Fitter.LoadData("kin",inputData);
    Fitter.LoadSimulated("kin",inputSimulated,config.GetName());
    // Fitter.Data().LoadWeights("Acc",inputWeight,"weight");

    config.SetLmax(2);
    config.SetMmax(2);
    config.SetNrefl(2);

    config.ConfigureMoments();
    config.LoadModelPDF();

    Fitter.SetUp().AddFitOption(RooFit::PrintEvalErrors(-1));
    // Fitter.SetUp().AddFitOption(RooFit::EvalBackend("legacy"));
    Fitter.SetUp().AddFitOption(RooFit::NumCPU(4));

    // Fitter.SetUp().ErrorsWrong();
    // Fitter.SetUp().ErrorsSumW2();
    // Fitter.SetUp().ErrorsAsymp();

    // auto mcmc=new BruMcmcSeqHelper(5000,1000,0.1,0.23,0.16,0.3);
    // auto mcmc=new BruMcmcCovariance(10000,2000,0.2,0.14,0.1,0.18);
    // auto mcmc = new BruMcmcCovariance(4000,1000,0.8,0.14,0.1,0.18);
    // auto mcmc = new BruMcmcCovariance(4000,1000,0.7,0.14,0.1,0.18);
    // auto mcmc = new BruMcmcCovariance(6000,1000,0.8,0.14,0.1,0.18);

    // Fitter.SetPlotOptions("MCMC:AUTOCORR");
    // auto mcmc = new BruMcmcCovariance({6000,6000,1000},1000,0.8,0.14,0.1,0.18);
    

    // _0
    // auto mcmc=new BruMcmcCovariance(5000,1000,0.1,0.07,0.04,0.12); // stuck

    // auto mcmc=new AmpMcmc(&config,15000,1000,0.01,20);

    // mcmc->TurnOffCovariance();
    Fitter.SetMinimiser(mcmc);

    Here::Go(&Fitter);


    // Fitter.InitPrevResult(Fitter.SetUp().GetOutDir(),Fitter.GetMinimiserType());
    auto minuit = new Minuit2();
    Fitter.SetMinimiser(minuit);
    Fitter.SetUp().ErrorsWrong();
    Here::Go(&Fitter);

    // delete mcmc;
    // delete minuit;

}

void MomentFit()
{

    TString inputDir="/d/grid17/septian/Moments2Amplitudes/samples/";
    TString inputData = Form("%sdata_sig_pol0.3_SD_waves.root",inputDir.Data());
    TString inputSimulated = Form("%s/accepted.root",inputDir.Data());
    TString inputWeight = Form("%s/weights_data_sig_bkg_pol0.5.root",inputDir.Data()); //not used
    TString outDir = "fitMoments_GlueXAcceptanceSignal_SDwaves_R6.34";
    TString tag = "Moment_GlueXAcceptanceSignal_SDwaves_R6.34";
    DoMomentFit(tag,inputData,inputSimulated,inputWeight,outDir);

    // Moment fit with seed from 1000 to 1100
    // for (Int_t seed=1005; seed<=1006; seed++){
    //     TString inputDir="/d/grid17/septian/EtaPi0Analysis/MomentMCStudy/data_seed_variation/";
    //     TString inputSimulated="/d/grid17/septian/EtaPi0Analysis/MomentMCStudy/accepted.root";
    //     TString inputData = Form("%s/data_sig_bkg_pol0.5_seed%d.root",inputDir.Data(),seed);
    //     TString inputWeight = Form("%s/weights_data_sig_bkg_pol0.5_seed%d.root",inputDir.Data(),seed);
    //     TString outDir = Form("fitMoment_seed%d/",seed);
    //     TString tag = Form("Moment_seed%d",seed);
    //     DoMomentFit(tag,inputData,inputSimulated,inputWeight,outDir);

    //     // reset root session
    //     gROOT->Reset();
    // }       
}
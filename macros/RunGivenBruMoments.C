#include "ConfigureAmpsNoValues.C"

void RunGivenBruMoments(){
  
  gRandom->SetSeed(0);
  
  //Define amplitudes (in additional file)
  auto& setup = ConfigureAmpsNoValues(2,2,2);//(Lmax,MMax,Nref)

  setup.SetParVal("aphi_0_0",0,kTRUE); //fix S real
  setup.SetParVal("bphi_0_0",0,kTRUE);

  // setup.SetParVal("b_0_0",1.0,kFALSE); //fix S real
  // setup.SetParVal("b_0_0",0.25,kFALSE);
  // setup.SetParVal("a_0_0",0.125,kFALSE);
  // setup.SetParVal("a_1_0",-0.499,kFALSE); //fix S real
  // setup.SetParVal("b_1_0",-0.999,kFALSE);

  // Set P-waves to 0
  // setup.SetParVal("a_1_0",0,kTRUE);
  // setup.SetParVal("b_1_0",0,kTRUE);
  // setup.SetParVal("a_1_1",0,kTRUE);
  // setup.SetParVal("b_1_1",0,kTRUE);
  // setup.SetParVal("a_1_-1",0,kTRUE);
  // setup.SetParVal("b_1_-1",0,kTRUE);
  // setup.SetParVal("aphi_1_0",0,kTRUE);
  // setup.SetParVal("bphi_1_0",0,kTRUE);
  // setup.SetParVal("aphi_1_1",0,kTRUE);
  // setup.SetParVal("bphi_1_1",0,kTRUE);
  // setup.SetParVal("aphi_1_-1",0,kTRUE);
  // setup.SetParVal("bphi_1_-1",0,kTRUE);

  //Fix with moments from ConfigureAmpsSPJune.C
  MomentHelper moments;
  // moments.Set("/d/home/septian/EtaPi0Analysis/MomentMCStudy/fitMoments_GlueXAcceptanceSigBkg_MCMCN4000BI1000S08/ResultsBruMcmcCovariance.root");

  // moments.Set("/d/grid17/septian/Moments2Amplitudes/brufit/fitMoments_GlueXAcceptanceSignal_SDwaves_R6.34/ResultsHSMinuit2.root",1,0);
  moments.Set("/d/home/septian/EtaPi0Analysis/run_Phase1/fitMoment_GlueXI_t010100_m080200_MCMCN6000BI1000S08WCOV/Mpi0eta1.420000_/ResultsBruMcmcCovariance.root",1,1);
  moments.PrintVals();
  //setup the fitter, arguments :
  //  setup = BruFit setup
  //  resolution = smear moments by resolution
  //  ignore_observables = do not include the following polarised moments
  //                       i.e. alpha = 0,1,2, or 3 => H_0,H_1,H_2,H_3
  m2pw::EquationSolver fitter{setup,0.0,{"H_3"}}; //ignore H_3
  //  m2pw::EquationSolver fitter{setup,0.05,{"H_3","H_2","H_1"}};
  fitter.SetEquationValues(moments);
  fitter.Print("v");

  fitter.GetPars().Randomise();
  fitter.Solve();

  //create output tree
  // fitter.MakeResultTree("resultsGivenBruMomentsNoH3_data_SDWaves_10000Sample_NoPWaves.root");
  // fitter.MakeResultTree("resultsGivenBruMomentsNoH3_data_MCMC_test_10000Sample.root");
  
  // gBenchmark->Start("fitter");

  // //loop and perform 10,000 minimisations with random starting amplitudes
  // for(int i = 0; i<10000;i++){
  //   if(i%100==0) cout<<i<<" "<<endl;
  //   fitter.GetPars().Randomise();
  //   fitter.Solve();
  //   fitter.FillTree();
  // }

  // gBenchmark->Stop("fitter");
  // gBenchmark->Print("fitter");

  // //Save results tree
  // fitter.GetPars().CloseTree();

  // fitter.PrintResult();
  
}
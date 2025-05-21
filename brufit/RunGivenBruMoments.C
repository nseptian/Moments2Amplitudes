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

  moments.Set("/d/grid17/septian/Moments2Amplitudes/brufit/fitMoments_GlueXAcceptanceSignal_SDwaves_R6.34/ResultsHSMinuit2.root");
  
  //setup the solver, arguments :
  //  setup = BruFit setup
  //  resolution = smear moments by resolution
  //  ignore_observables = do not include the following polarised moments
  //                       i.e. alpha = 0,1,2, or 3 => H_0,H_1,H_2,H_3
  m2pw::EquationSolver solver{setup,0.0,{"H_3"}}; //ignore H_3
  //  m2pw::EquationSolver solver{setup,0.05,{"H_3","H_2","H_1"}};
  solver.SetEquationValues(moments);
  solver.Print("v");

  //create output tree
  solver.MakeResultTree("resultsGivenBruMomentsNoH3_data_SDWaves_10000Sample_NoPWaves.root");
  
  gBenchmark->Start("solver");

  //loop and perform 10,000 minimisations with random starting amplitudes
  for(int i = 0; i<10000;i++){
    if(i%100==0) cout<<i<<" "<<endl;
    solver.GetPars().Randomise();
    solver.Solve();
    solver.FillTree();
  }

  gBenchmark->Stop("solver");
  gBenchmark->Print("solver");

  //Save results tree
  solver.GetPars().CloseTree();

  solver.PrintResult();
  
}
//Run with
//root --hsfit FitHSAsymmetryBins.C
void FitSDME(TString dataset = "2017-01", Bool_t Ebinned = kFALSE){
	
	TString datafile, sigfile, outdir, genfile, weightsDir;
	if(dataset=="2017-01"){ 
		datafile = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/trees/data_cuts_lineshape/2017-01_tree_baseline.root";
		sigfile = "/w/work5/home/pauli/Lambda1520/MC/lineshape_pKmKp__t_30_07__s17_ver40_4038/2017_01_tree_baseline.root";
		genfile  = "/w/work5/home/pauli/Lambda1520/MC/lineshape_pKmKp__t_30_07__s17_ver40_4038/hs_tree_gen_kpkm.root";
		if(Ebinned){
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/lineshape_fine_Ebins/weights/";
			outdir     = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/lineshape_fine_Ebins/fit/";
		}
		else{
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/lineshape_fine/weights/";
			outdir     = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/lineshape_fine/fit/";
		}

	}
	if(dataset=="2018-01"){
		datafile = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/trees/data_cuts_lineshape/2018-01_tree_baseline.root";
		sigfile = "/w/work5/home/pauli/Lambda1520/MC/lineshape_pKmKp__t_30_07__s18_ver32_4039/2018_01_tree_baseline.root";
		genfile  = "/w/work5/home/pauli/Lambda1520/MC/lineshape_pKmKp__t_30_07__s18_ver32_4039/hs_tree_gen_kpkm.root";
		if(Ebinned){
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/lineshape_fine_Ebins/weights/";
			outdir     = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/lineshape_fine_Ebins/fit/";
		}
		else{
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/lineshape_fine/weights/";
			outdir     = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/lineshape_fine/fit/";
		}
	}
	if(dataset=="2018-08"){
		datafile = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/trees/data_cuts_lineshape/2018-08_tree_baseline.root";
		sigfile = "/w/work5/home/pauli/Lambda1520/MC/lineshape_pKmKp__t_30_07__f18_ver31_4040/2018_08_tree_baseline.root";
		genfile  = "/w/work5/home/pauli/Lambda1520/MC/lineshape_pKmKp__t_30_07__f18_ver31_4040/hs_tree_gen_kpkm.root";
		if(Ebinned){
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine_Ebins/weights/";
			outdir     = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine_Ebins/fit/";
		}
		else{
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine/weights/";
			outdir     = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine/fit/";
		}
	}

	weightsDir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine_t1.0..1.5/weights/";
	outdir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine_t1.0..1.5/fit/";


	//JLAB
	// TString datafile = "/volatile/halld/home/ppauli/lineshape/kpkm/hs_tree_kpkm_2018-08.root";
	// TString sigfile  = "/work/halld/home/ppauli/data/MC/lineshape_pKmKp__t_30_07__f18_ver25/hs_tree_kpkm.root";
	// TString genfile  = "/work/halld/home/ppauli/data/MC/lineshape_pKmKp__t_30_07__f18_ver25/hs_tree_gen_kpkm.root";
	// TString weightsDir = "/volatile/halld/home/ppauli/lineshape/kpkm/brufit/2018-08/weights/";
	// TString outdir     = "/volatile/halld/home/ppauli/lineshape/kpkm/brufit/2018-08/fit/";


	/****************************************/
	/************Create FitManager***********/
	/****************************************/
	FitManager RF;
	RF.SetUp().SetOutDir(outdir);

	/****************************************/
	/*************Load Variables*************/
	/****************************************/
	RF.SetUp().LoadVariable("Theta_GJ_KMinus[0,3.2]");
	RF.SetUp().LoadVariable("phiAngle_GJ_KMinus[-3.2,3.2]");
	RF.SetUp().LoadVariable("phiAngle_KPlus[-3.2,3.2]");
	RF.SetUp().LoadVariable("polarisationAngle[-1.5,3.2]");
	RF.SetUp().LoadVariable("polarisationDegree[-0.2,1]");

	RF.SetUp().LoadAuxVar("ConfidenceLevel[0,1]");
	RF.SetUp().AddCut("ConfidenceLevel>1e-6");
	// RF.SetUp().LoadAuxVar("momentumTransferMinustmin[0,1.5]");
	RF.SetUp().LoadAuxVar("momentumTransferMinustmin[1.0,1.5]");

	RF.SetUp().LoadAuxVar("runNumber[0,1E6]");

	RF.SetUp().LoadFormula("A=3./(8.*3.14159265359)*TMath::Sin(@Theta_GJ_KMinus[])*TMath::Sin(@Theta_GJ_KMinus[])");
	RF.SetUp().LoadFormula("B=1./(2.*3.14159265359)*(3.*TMath::Cos(@Theta_GJ_KMinus[])*TMath::Cos(@Theta_GJ_KMinus[])-1.)");
	RF.SetUp().LoadFormula("C=(-1.)*TMath::Sqrt(3.)/(2.*3.14159265359)*TMath::Cos(@phiAngle_GJ_KMinus[])*TMath::Sin(2.*@Theta_GJ_KMinus[])");
	RF.SetUp().LoadFormula("D=(-1.)*TMath::Sqrt(3.)/(2.*3.14159265359)*TMath::Cos(2.*@phiAngle_GJ_KMinus[])*TMath::Sin(@Theta_GJ_KMinus[])*TMath::Sin(@Theta_GJ_KMinus[])");
	RF.SetUp().LoadFormula("E=(-1.)*3./(4.*3.14159265359)*@polarisationDegree[]*TMath::Cos(2*(@polarisationAngle[]-@phiAngle_KPlus[]))*TMath::Sin(@Theta_GJ_KMinus[])*TMath::Sin(@Theta_GJ_KMinus[])");
	RF.SetUp().LoadFormula("F=(-1.)*3./(4.*3.14159265359)*@polarisationDegree[]*TMath::Cos(2*(@polarisationAngle[]-@phiAngle_KPlus[]))*(1./3.+TMath::Cos(@Theta_GJ_KMinus[])*TMath::Cos(@Theta_GJ_KMinus[]))");
	RF.SetUp().LoadFormula("G=TMath::Sqrt(3.)/(2.*3.14159265359)*@polarisationDegree[]*TMath::Cos(2*(@polarisationAngle[]-@phiAngle_KPlus[]))*TMath::Cos(@phiAngle_GJ_KMinus[])*TMath::Sin(2.*@Theta_GJ_KMinus[])");
	RF.SetUp().LoadFormula("H=TMath::Sqrt(3.)/(2.*3.14159265359)*@polarisationDegree[]*TMath::Cos(2*(@polarisationAngle[]-@phiAngle_KPlus[]))*TMath::Cos(2.*@phiAngle_GJ_KMinus[])*TMath::Sin(@Theta_GJ_KMinus[])*TMath::Sin(@Theta_GJ_KMinus[])");
	RF.SetUp().LoadFormula("I=(-1.)*TMath::Sqrt(3.)/(2.*3.14159265359)*@polarisationDegree[]*TMath::Sin(2*(@polarisationAngle[]-@phiAngle_KPlus[]))*TMath::Sin(@phiAngle_GJ_KMinus[])*TMath::Sin(2.*@Theta_GJ_KMinus[])");
	RF.SetUp().LoadFormula("J=(-1.)*TMath::Sqrt(3.)/(2.*3.14159265359)*@polarisationDegree[]*TMath::Sin(2*(@polarisationAngle[]-@phiAngle_KPlus[]))*TMath::Sin(2.*@phiAngle_GJ_KMinus[])*TMath::Sin(@Theta_GJ_KMinus[])*TMath::Sin(@Theta_GJ_KMinus[])");

	RF.SetUp().LoadParameter("Rho011[0,-1,1]");
	RF.SetUp().LoadParameter("Rho031[0,-1,1]");
	RF.SetUp().LoadParameter("Rho03m1[0,-1,1]");
	RF.SetUp().LoadParameter("Rho111[0,-1,1]");
	RF.SetUp().LoadParameter("Rho133[0,-1,1]");
	RF.SetUp().LoadParameter("Rho131[0,-1,1]");
	RF.SetUp().LoadParameter("Rho13m1[0,-1,1]");
	RF.SetUp().LoadParameter("Rho231[0,-1,1]");
	RF.SetUp().LoadParameter("Rho23m1[0,-1,1]");

	RF.SetUp().SetIDBranchName("uniqueComboID");


	/**************************************************/
	/********* Need to create sideband weights ********/
	/**************************************************/

	TDirectory* saveDir=gDirectory;
	if(true){//make new sideband weights
		TFile* treefileSig=new TFile(sigfile);
		TTree* treeSig=(TTree*) treefileSig->Get("haspect");
		saveDir->cd();
		Weights* wgtsSBSig=new Weights("RndTime");
		wgtsSBSig->SetFile(outdir+"/SideBandSig.root");
		wgtsSBSig->SetSpecies("SideBand");
		wgtsSBSig->SetIDName("uniqueComboID");
		// wgtsSBSig->WeightBySelection(treeSig,"(1)","timingWeight");
		wgtsSBSig->WeightBySelection(treeSig,"(scaledtimingWeight==comboWeight)","scaledtimingWeight");
	// 	wgtsSBSig->WeightBySelection(treeSig,"(deltaRFTiming>-6.012&&deltaRFTiming<6.012)",1);
	// 	wgtsSBSig->WeightBySelection(treeSig,"(deltaRFTiming>-18.036&&deltaRFTiming<-6.012)||(deltaRFTiming>6.012&&deltaRFTiming<18.036)",-0.5); //Note : sideband weight is proportional to the time windows range ratio
		wgtsSBSig->SortWeights();
		wgtsSBSig->Save();//Save to disc
		delete wgtsSBSig; wgtsSBSig=nullptr;
	}
	
	/****************************************/
	/*************Make model PDF*************/
	/****************************************/
	//1 + A*COS2 + B*SIN2 (first argument = the 1 )
	//Note : seperates different products
	//     ; seperates different terms in the produt
	RF.SetUp().FactoryPDF("RooComponentsPDF::SDMES(0,{Theta_GJ_KMinus,phiAngle_GJ_KMinus,phiAngle_KPlus,polarisationAngle,polarisationDegree},=A:B;Rho011:C;Rho031:D;Rho03m1:E;Rho133:F;Rho111:G;Rho131:H;Rho13m1:I;Rho231:J;Rho23m1)WEIGHTS@SideBand,"+RF.SetUp().GetOutDir()+"/SideBandSig.root,RndTime");
// 	RF.SetUp().FactoryPDF("RooComponentsPDF::SDMES(0,{Theta_GJ_KMinus,phiAngle_GJ_KMinus,phiAngle_KPlus,polarisationAngle,polarisationDegree},=A:B;Rho011:C;Rho031:D;Rho03m1:E;Rho133:F;Rho111:G;Rho131:H;Rho13m1:I;Rho231:J;Rho23m1)WEIGHTS@LikeData,"+weightsDir+"/impWeights.root,impWeights");

	RF.SetUp().LoadSpeciesPDF("SDMES",1);

	/**************************************************/
	/********************Make bins*********************/
	/**************************************************/

	// Double_t binLimitst[] = {0.0,0.35,0.6,1.5}; //final binning
	// RF.Bins().LoadBinVar("momentumTransferMinustmin",3,binLimitst);

	RF.Bins().LoadBinVar("invMass_PKMinus",54,1.43,1.70);

	Double_t binLimitsEbinned[] = {6.5,7.0,7.6,8.2,8.8,9.6,10.4,11.2}; //final binning
	Double_t binLimitsE[] = {6.5,11.6}; //final binning
	if(Ebinned)
		RF.Bins().LoadBinVar("photonBeamEnergy",7,binLimitsEbinned);
	else
		RF.Bins().LoadBinVar("photonBeamEnergy",1,binLimitsE);


// 	RF.Bins().AddCut("invMass_PKMinus<1.6");


	/**************************************************/
	/****************Load data and MC******************/
	/**************************************************/
	RF.LoadData("haspect",datafile);
	RF.LoadSimulated("haspect",sigfile,"SDMES");
	RF.LoadGenerated("haspect",genfile,"SDMES",kTRUE);

	RF.Data().LoadWeights("Sig",weightsDir+"/Tweights.root");


	/**************************************************/
	/***********Choose minimiser and run***************/
	/**************************************************/
	RF.SetUp().AddFitOption(RooFit::Optimize(1));
	RF.SetUp().AddFitOption(RooFit::PrintEvalErrors(-1));//suppress error messaages
	// RF.SetUp().AddFitOption(RooFit::NumCPU(6)); //number of CPUs to split likelihood calc.

// 	RF.SetMinimiser(new Minuit2());

  //Or try an mcmc minimser 1000-># of points, 200->burnin 200 ~ 1/step size
	// RF.SetMinimiser(new RooMcmcSeq(50000,100,6)); //publication

	// RF.SetMinimiser(new RooMcmcSeq(50,0,6));
// 	RF.SetMinimiser(new RooMcmcUniform2Seq(1000,200,10));
// 	RF.SetMinimiser(new RooMcmcGaus(10,0,1000));

	///////////////////////////////(Nsamples,burnin,step size, desired acceptance, min acceptance, max acceptance)
	auto mcmc=new BruMcmcCovariance(5000,   1000,  0.1,          0.23,              0.16,           0.3);
	////mcmc->TurnOffCovariance();//BruMcmcCovariance only, do not proceed with covariance based sampling, just perform basic stepping
	RF.SetMinimiser(mcmc);
	cout << "Start Fit" << endl;
	// Here::Go(&RF);
	if (Ebinned)
		Proof::Go(&RF,378);
	else
		Proof::Go(&RF,54);

}

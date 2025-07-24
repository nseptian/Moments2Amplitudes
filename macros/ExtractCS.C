//Run with
//root --hsfit FitHSAsymmetryBins.C
void ExtractCS(TString dataset = "2017-01", Bool_t Ebinned = kFALSE){
	
	TString fitdir, weightsDir;
	if(dataset=="2017-01"){
		if(Ebinned){
			fitdir     = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/lineshape_fine_Ebins/fit/";
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/lineshape_fine_Ebins/weights/";
		}
		else{
			fitdir     = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/lineshape_fine/fit/";
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/lineshape_fine/weights/";
		}
	}
	if(dataset=="2018-01"){
		if(Ebinned){
			fitdir     = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/lineshape_fine_Ebins/fit/";
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/lineshape_fine_Ebins/weights/";
		}
		else{
			fitdir     = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/lineshape_fine/fit/";
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/lineshape_fine/weights/";
		}
	}
	if(dataset=="2018-08"){
		if(Ebinned){
			fitdir     = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine_Ebins/fit/";
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine_Ebins/weights/";
		}
		else{
			fitdir     = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine/fit/";
			weightsDir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine/weights/";
		}
	}

	weightsDir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine_t1.0..1.5/weights/";
	fitdir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine_t1.0..1.5/fit/";

	//JLAB
	// TString fitdir     = "/volatile/halld/home/ppauli/lineshape/kpkm/brufit/2018-08/fit/";
	// TString weightsDir = "/volatile/halld/home/ppauli/lineshape/kpkm/brufit/2018-08/weights/";

	/****************************************/
	/************Create CrossSection*********/
	/****************************************/
	CrossSection RF;
	RF.SetUp().SetIDBranchName("uniqueComboID");

	RF.SetUp().SetOutDir(fitdir);
	RF.SetResultDir(fitdir);
	// RF.SetResultFileName("ResultsHSRooMcmcSeq.root");
	RF.SetResultFileName("ResultsBruMcmcCovariance.root");


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

	/***************************************/
	/*************Get model PDF*************/
	/***************************************/
	//1 + A*COS2 + B*SIN2 (first argument = the 1 )
	//Note : seperates different products
	//     ; seperates different terms in the produt
	RF.SetUp().FactoryPDF("RooComponentsPDF::SDMES(0,{Theta_GJ_KMinus,phiAngle_GJ_KMinus,phiAngle_KPlus,polarisationAngle,polarisationDegree},=A:B;Rho011:C;Rho031:D;Rho03m1:E;Rho133:F;Rho111:G;Rho131:H;Rho13m1:I;Rho231:J;Rho23m1)WEIGHTS@SideBand,"+RF.SetUp().GetOutDir()+"/SideBandSig.root,RndTime");

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

	/**************************************************/
	/****************Load data and MC******************/
	/**************************************************/
	RF.ReloadData(fitdir);
	RF.ReloadSimulated(fitdir,"SDMES");
	RF.ReloadGenerated(fitdir,"SDMES");

	RF.Data().LoadWeights("Sig",weightsDir+"/Tweights.root");

	/**************************************************/
	/***********Calculate cross section****************/
	/**************************************************/
	cout << "Start cross section calculation" << endl;
// 	vector<Double_t> energybin = {8.2,8.8};
// 	RF.SetBeamEnergyBinLimits(energybin); //set directly when no binning in E was done
	RF.SetBeamEnergyBinLimits("photonBeamEnergy");
	if(dataset=="2017-01") RF.SetFlux("/home/pauli/Documents/GlueX/Lambda1520/cs/flux/flux_30274_31057.root","tagged_lumi"); // spring 17
	if(dataset=="2018-01") RF.SetFlux("/home/pauli/Documents/GlueX/Lambda1520/cs/flux/flux_40856_42559.root","tagged_lumi"); //spring 18
	if(dataset=="2018-08") RF.SetFlux("/home/pauli/Documents/GlueX/Lambda1520/cs/flux/flux_50677_51768.root","tagged_lumi"); //fall 18
 	// RF.SetFlux("/work/halld/home/ppauli/data/flux/flux_50677_51768.root","tagged_lumi"); // fall 18 at JLab
	RF.SetTargetThickness(1e6); //convert lumi in 1/pb to 1/µb
	RF.SetBranchingRatio(1);
	// RF.SampleAcceptance(); //sample acceptance from MCMC instead of doing the integral with the fit results

	// Here::Go(&RF);
	if (Ebinned)
		Proof::Go(&RF,378);
	else
		Proof::Go(&RF,54);


	/*********************************************/
	/***********Draw cross section****************/
	/*********************************************/
	RF.DrawResults(fitdir+"/CrossSectionAllE.root");

}

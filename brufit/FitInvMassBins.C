//brufit  FitInvMassBins.C
void FitInvMassBins(TString dataset = "2017-01", Bool_t Ebinned = kFALSE){

	// Loader::Compile("$HSCODE/hsfit/RelBreitWigner.cxx+");

	TString datafile, sigfile, outdir;
	if(dataset=="2017-01"){ // Spring '17
		datafile = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/trees/data_cuts_lineshape/2017-01_tree_baseline.root";
		sigfile = "/w/work5/home/pauli/Lambda1520/MC/lineshape_pKmKp__t_30_07__s17_ver40_4038/2017_01_tree_baseline.root";
		if(Ebinned)
			outdir = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/lineshape_fine_Ebins/weights/";
		else
			outdir = "/w/work5/home/pauli/Lambda1520/data/spring17_ver46/lineshape_fine/weights/";
	}
	if(dataset=="2018-01"){ // Spring '18
		datafile = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/trees/data_cuts_lineshape/2018-01_tree_baseline.root";
		sigfile = "/w/work5/home/pauli/Lambda1520/MC/lineshape_pKmKp__t_30_07__s18_ver32_4039/2018_01_tree_baseline.root";
		if(Ebinned)
			outdir = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/lineshape_fine_Ebins/weights/";
		else
			outdir = "/w/work5/home/pauli/Lambda1520/data/spring18_ver15/lineshape_fine/weights/";
	}
	if(dataset=="2018-08"){ // // Fall '18
		datafile = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/trees/data_cuts_lineshape/2018-08_tree_baseline.root";
		sigfile = "/w/work5/home/pauli/Lambda1520/MC/lineshape_pKmKp__t_30_07__f18_ver31_4040/2018_08_tree_baseline.root";
		if(Ebinned)
			outdir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine_Ebins/weights/";
		else
			outdir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine/weights/";
	}

	outdir = "/w/work5/home/pauli/Lambda1520/data/fall18_ver15/lineshape_fine_t1.0..1.5/weights/";

	//JLAB
	// TString datafile = "/volatile/halld/home/ppauli/lineshape/kpkm/hs_tree_kpkm_2018-08.root";
	// TString sigfile  = "/work/halld/home/ppauli/data/MC/lineshape_pKmKp__t_30_07__f18_ver25/hs_tree_kpkm.root";
	// TString outdir   = "/volatile/halld/home/ppauli/lineshape/kpkm/brufit/2018-08/weights/";

	sPlot RF;
	RF.SetUp().SetOutDir(outdir);

	/****************************************/
	/*************Load Variables*************/
	/****************************************/
	// RF.SetUp().LoadVariable("invMass_PKMinus[1.43,1.70]");//should be same name as variable in tree
	RF.SetUp().LoadVariable("momentumTransferMinustmin[0,1.5]");
	// RF.SetUp().LoadVariable("momentumTransferMinustmin[1.0,1.5]");
// 	RF.SetUp().LoadAuxVar("photonBeamEnergy[8.2,8.8]");
	// RF.SetUp().LoadAuxVar("deltaRFTiming[-20,20]");

	RF.SetUp().LoadAuxVar("ConfidenceLevel[0,1]");
	RF.SetUp().AddCut("ConfidenceLevel>1e-6");

	RF.SetUp().LoadAuxVar("runNumber[0,1E6]");

	RF.SetUp().SetIDBranchName("uniqueComboID");


	/**************************************************/
	/********* Need to create sideband weights ********/
	/**************************************************/
	TDirectory* saveDir=gDirectory;
	TFile* treefileData=new TFile(datafile);
	TTree* treeData=(TTree*) treefileData->Get("haspect");
	saveDir->cd();
	Weights* wgtsSBData=new Weights("HSsWeights");
	wgtsSBData->SetFile(RF.SetUp().GetOutDir()+"/SideBandData.root");
	wgtsSBData->SetSpecies("SideBand");
	wgtsSBData->SetIDName("uniqueComboID");
// 	wgtsSBData->WeightBySelection(treeData,"(1)","timingWeight");
	wgtsSBData->WeightBySelection(treeData,"(scaledtimingWeight==comboWeight)","scaledtimingWeight");
	wgtsSBData->SortWeights();
	wgtsSBData->Save();//Save to disc
	delete wgtsSBData; wgtsSBData=nullptr;

	saveDir=gDirectory;
	TFile* treefileSig=new TFile(sigfile);
	TTree* treeSig=(TTree*) treefileSig->Get("haspect");
	saveDir->cd();
	Weights* wgtsSBSig=new Weights("RndTime");
	wgtsSBSig->SetFile(RF.SetUp().GetOutDir()+"SideBandSig.root");
	wgtsSBSig->SetSpecies("SideBand");
	wgtsSBSig->SetIDName("uniqueComboID");
// 	wgtsSBSig->WeightBySelection(treeSig,"(1)","timingWeight");
	wgtsSBSig->WeightBySelection(treeSig,"(scaledtimingWeight==comboWeight)","scaledtimingWeight");
	wgtsSBSig->SortWeights();
	wgtsSBSig->Save();//Save to disc
	delete wgtsSBSig; wgtsSBSig=nullptr;

	
	/****************************************/
	/*************Make signal PDF************/
	/****************************************/
	RF.SetUp().LoadParameter("A[0.7,0,2]");
	RF.SetUp().LoadParameter("B[-3,-10,0]");
	RF.SetUp().FactoryPDF("EXPR::Sig('momentumTransferMinustmin**A * exp(B*momentumTransferMinustmin)', momentumTransferMinustmin, A, B) ");
	RF.SetUp().LoadSpeciesPDF("Sig",1);
	
	// RF.SetUp().FactoryPDF("RooHSEventsHistPDF::Sig(momentumTransferMinustmin,alpha[0.0025,0.0,1.0],off[0,-1,1],scale[1,0.5,1.5])WEIGHTS@SideBand,"+RF.SetUp().GetOutDir()+"/SideBandSig.root,RndTime");
	// RooHSEventsHistPDF* pdfSig=dynamic_cast<RooHSEventsHistPDF*>(RF.SetUp().WS().pdf("Sig"));
	// RF.SetUp().LoadSpeciesPDF("Sig",1);

	// RF.SetUp().FactoryPDF("RelBreitWigner::Sig(invMass_PKMinus,mK[0.493677,0.493677,0.493677],mp[0.938272046,0.938272046,0.938272046],L[2,2,2],mean[1.519,1.5,1.54],width[0.02,0,0.1])");
	// RF.SetUp().LoadSpeciesPDF("Sig",1);


	/**************************************************/
	/*************Make background PDFs*****************/
	/**************************************************/
	// RF.SetUp().FactoryPDF("Chebychev::Cheb(invMass_PKMinus,{ch0[0,-1,1],ch1[0,-1,1]})");
	// RF.SetUp().LoadSpeciesPDF("Cheb",1);


	/****************************************************************/
	/*************Add constraints to PDF fudge parameters************/
	/****************************************************************/
	//These are derived from the alpha off and scale parameter initial val and limits
	//i.e. intital val = gaussian mean; range = 10*sigma
	// RF.SetUp().AddGausConstraint(pdfSig->AlphaConstraint());
	// RF.SetUp().AddGausConstraint(pdfSig->OffConstraint());
	// RF.SetUp().AddGausConstraint(pdfSig->ScaleConstraint());

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
	// RF.LoadData("haspect",datafile, "Data");
	RF.LoadData("haspect",datafile); //new brufit doesn't allow data label
	RF.LoadSimulated("haspect",sigfile, "Sig");
// 	RF.ReloadData(datafile, "Data");
// 	RF.ReloadSimulated(sigfile, "Sig");

	RF.Data().LoadWeights("SideBand",RF.SetUp().GetOutDir()+"/SideBandData.root");

	gBenchmark->Start("Binned");

	// Here::Go(&RF);
	if (Ebinned)
		Proof::Go(&RF,378);
	else
		Proof::Go(&RF,54);
	gBenchmark->Stop("Binned");
	gBenchmark->Print("Binned");

	new TCanvas;
	RF.DrawWeighted("missingMassSq>>(80,-0.04,0.04)","Sig");
	// RF.DrawWeighted("missingMassSq>>(80,-0.04,0.04)","Cheb","1","same");

  //make sure weighted tree is written properly
  RF.DeleteWeightedTree();

}


int MakeMomentFitData(TString inFile="",TString outFile=""){
    TFile *f = TFile::Open(inFile);
    if (!f->IsOpen()) {
        std::cerr << "Could not open the file!" << std::endl;
        return 1;
    }
    Float_t in_cosTheta_eta_hel, in_phi_eta_hel, in_Phi;
    TTree *t = (TTree*)f->Get("kin");
    t->SetBranchAddress("cosTheta_eta_hel",&in_cosTheta_eta_hel);
    t->SetBranchAddress("phi_eta_hel",&in_phi_eta_hel);
    t->SetBranchAddress("Phi",&in_Phi);

    Double_t out_cosTheta_eta_hel, out_phi_eta_hel, out_Phi;
    TFile *fout = new TFile(outFile,"RECREATE");
    TTree *tOut = new TTree("kin","kin");
    tOut->Branch("cosTheta_eta_hel",&out_cosTheta_eta_hel);
    tOut->Branch("phi_eta_hel",&out_phi_eta_hel);
    tOut->Branch("Phi",&out_Phi);

    for (Int_t i=0; i<t->GetEntries(); i++){
        t->GetEntry(i);
        out_cosTheta_eta_hel = in_cosTheta_eta_hel;
        out_phi_eta_hel = TMath::DegToRad()*in_phi_eta_hel;
        out_Phi = TMath::DegToRad()*in_Phi;
        tOut->Fill();
    }

    fout->Write();
    fout->Close();
    f->Close();
    return 0;
}
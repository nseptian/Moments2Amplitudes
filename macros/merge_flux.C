
const TString fInput[4] = {"flux_30274_31057.root","flux_40856_42577.root","flux_50677_51768.root","flux_71350_73266.root"};

void merge_flux() {
    const TString& outputFile = "flux_GlueXI_GlueX2020.root";
    const TString hname = "tagged_lumi";

    TFile *output = TFile::Open(outputFile, "RECREATE");
    if (!output || output->IsZombie()) {
        std::cerr << "Error opening output file: " << outputFile << std::endl;
        return;
    }
    output->cd();
    TH1D *mergedFlux = new TH1D("tagged_lumi", "Tagged Luminosity", 1, 8.0, 8.8);

    Double_t tagged_lumi = 0.0;
    for (const auto& inputFile : fInput) {
        TFile *input = TFile::Open(inputFile);
        if (!input || input->IsZombie()) {
            std::cerr << "Error opening input file: " << inputFile << std::endl;
            continue;
        }

        TH1D *fluxHist = dynamic_cast<TH1D*>(input->Get(hname));
        if (!fluxHist) {
            std::cerr << "Histogram " << hname << " not found in file: " << inputFile << std::endl;
            input->Close();
            continue;
        }

        tagged_lumi += fluxHist->Integral();
        std::cout << "Adding " << inputFile << ": " << fluxHist->Integral() << " to total luminosity." << std::endl;
    }

    std::cout << "Total tagged luminosity: " << tagged_lumi << std::endl;
    mergedFlux->SetBinContent(1, tagged_lumi);

    output->cd();

    mergedFlux->Write();
    output->Close();

}
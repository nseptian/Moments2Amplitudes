#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"

void TestMassDependentMoments(){
    auto& setup = ConfigureAmpsNoValues(2,2,2); // (Lmax,MMax,Nref)
    setup.SetParVal("aphi_0_0",0,kTRUE); //fix S real
    setup.SetParVal("bphi_0_0",0,kTRUE);

    MassDependentMoments massDepMoments;

    const Double_t firstMassBinCenter = 0.82; // in GeV
    const Double_t massBinWidth = 0.04;
    const Int_t nMassBins = 5;
    const TString fitResultsDir = "/d/home/septian/EtaPi0Analysis/run_Phase1/fitMoment_GlueXI_t010100_m080200_MCMCN6000BI1000S08WCOV/";
    const TString fitResultsFilename = "ResultsBruMcmcCovariance.root";

    TFile *f[nMassBins];

    for (Int_t i = 0; i < nMassBins; ++i) {
        // cout << "Mass bin " << i << ": ";
        Double_t massValue = firstMassBinCenter + i * massBinWidth;
        // cout << "Mass value = " << massValue << endl;

        const TString massBinDirName = Form("Mpi0eta%1.6f_/", massValue);
        // cout << "Mass bin directory name: " << massBinDirName << endl;
        
        const TString filePath = fitResultsDir + massBinDirName + fitResultsFilename;

        cout << "Loading moments from file: " << filePath << endl;

        f[i] = TFile::Open(filePath);
        if (!f[i] || f[i]->IsZombie()) {
            cerr << "Error opening file: " << filePath << endl;
            continue;
        }

        MomentHelper* moments = new MomentHelper();
        
        moments->Set(filePath, 1, 1);

        // cout << "Moments loaded for mass value: " << massValue << endl;
        // moments->PrintVals();
        // cout << endl;

        
        massDepMoments.SetMoments(massValue, *moments);
        // cout << "Moments set for mass value: " << massValue << endl;

    }
    // Print all moments for each mass value
    // cout << "Mass-dependent moments:" << endl;
    massDepMoments.PrintMoments();
}
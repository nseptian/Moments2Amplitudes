// Helper macro to export resonance parameters and couplings to JSON for cross-section evaluation
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>

void ExportResonanceCouplings(
    const char* rootFilePath,
    const char* resonanceName,
    const char* outJsonPath
) {
    TFile* f = TFile::Open(rootFilePath, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open ROOT file " << rootFilePath << std::endl;
        if (f) delete f;
        return;
    }

    TTree* t = (TTree*)f->Get("mass_dependent_params");
    if (!t) {
        std::cerr << "Error: tree 'mass_dependent_params' not found in " << rootFilePath << std::endl;
        f->Close();
        delete f;
        return;
    }

    std::ofstream out(outJsonPath);
    if (!out.is_open()) {
        std::cerr << "Error: cannot open output file " << outJsonPath << std::endl;
        f->Close();
        delete f;
        return;
    }

    t->GetEntry(0);
    TString res(resonanceName);

    // 1. Mass and Width
    double massVal = 0.0, massErr = 0.0;
    double widthVal = 0.0, widthErr = 0.0;

    TString massBranch = Form("MD_shared_%s_Mass", res.Data());
    TString massErrBranch = Form("MD_shared_%s_Mass_error", res.Data());
    TString widthBranch = Form("MD_shared_%s_Width", res.Data());
    TString widthErrBranch = Form("MD_shared_%s_Width_error", res.Data());

    if (t->GetLeaf(massBranch)) massVal = t->GetLeaf(massBranch)->GetValue();
    if (t->GetLeaf(massErrBranch)) massErr = t->GetLeaf(massErrBranch)->GetValue();
    if (t->GetLeaf(widthBranch)) widthVal = t->GetLeaf(widthBranch)->GetValue();
    if (t->GetLeaf(widthErrBranch)) widthErr = t->GetLeaf(widthErrBranch)->GetValue();

    out << "{\n";
    out << "  \"resonance\": \"" << res.Data() << "\",\n";
    out << "  \"mass\": " << std::setprecision(8) << massVal << ",\n";
    out << "  \"mass_error\": " << std::setprecision(8) << massErr << ",\n";
    out << "  \"width\": " << std::setprecision(8) << widthVal << ",\n";
    out << "  \"width_error\": " << std::setprecision(8) << widthErr << ",\n";
    out << "  \"couplings\": {\n";

    // 2. Iterate branches to find MD_{wave}_{resonance}_k
    TObjArray* branches = t->GetListOfBranches();
    bool first = true;

    for (int i = 0; i < branches->GetEntries(); ++i) {
        TString bName = branches->At(i)->GetName();
        // Pattern: MD_{amp}_{res}_k, but NOT ending in _error
        if (bName.BeginsWith("MD_") && bName.EndsWith("_k") && !bName.EndsWith("_error")) {
            TString sub = bName;
            sub.ReplaceAll("MD_", "");
            sub.ReplaceAll("_k", "");
            // Now sub is {amp}_{res}
            if (sub.EndsWith(Form("_%s", res.Data()))) {
                TString waveName = sub;
                waveName.ReplaceAll(Form("_%s", res.Data()), "");
                
                double kVal = t->GetLeaf(bName)->GetValue();
                double kErr = 0.0;
                TString errBranch = Form("%s_error", bName.Data());
                if (t->GetLeaf(errBranch)) {
                    kErr = t->GetLeaf(errBranch)->GetValue();
                }

                int refl = 0;
                if (waveName.BeginsWith("a_")) refl = 1;
                else if (waveName.BeginsWith("b_")) refl = -1;

                if (!first) out << ",\n";
                first = false;

                out << "    \"" << waveName.Data() << "\": {\n";
                out << "      \"k\": " << std::setprecision(8) << kVal << ",\n";
                out << "      \"k_error\": " << std::setprecision(8) << kErr << ",\n";
                out << "      \"reflectivity\": " << refl << "\n";
                out << "    }";
            }
        }
    }

    out << "\n  }\n}\n";
    out.close();
    f->Close();
    delete f;

    std::cout << "[OK] Exported resonance couplings to: " << outJsonPath << std::endl;
}

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
    if (!t) t = (TTree*)f->Get("central_params");
    if (!t) t = (TTree*)f->Get("all_params");
    if (!t) {
        std::cerr << "Error: tree 'mass_dependent_params' (or 'central_params'/'all_params') not found in " << rootFilePath << std::endl;
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

    // Fallbacks for two-breit-wigner or alternate naming schemes
    if (massVal <= 0.0) {
        std::vector<TString> altMass;
        if (res.Contains("1700") || res.BeginsWith("res2")) {
            altMass = {
                "MD_shared_a2_1320_1700_Mass2",
                Form("MD_shared_%s_Mass2", res.Data()),
                Form("MD_shared_%s_Mass", res.Data()),
                Form("MD_shared_%s_Mass1", res.Data())
            };
        } else {
            altMass = {
                "MD_shared_a2_1320_1700_Mass1",
                Form("MD_shared_%s_Mass1", res.Data()),
                Form("MD_shared_%s_Mass", res.Data()),
                "MD_shared_a2_1320_Mass1"
            };
        }
        for (const auto& alt : altMass) {
            if (t->GetLeaf(alt)) {
                massVal = t->GetLeaf(alt)->GetValue();
                TString altErr = alt + "_error";
                if (t->GetLeaf(altErr)) massErr = t->GetLeaf(altErr)->GetValue();
                break;
            }
        }
    }
    if (widthVal <= 0.0) {
        std::vector<TString> altWidth;
        if (res.Contains("1700") || res.BeginsWith("res2")) {
            altWidth = {
                "MD_shared_a2_1320_1700_Width2",
                Form("MD_shared_%s_Width2", res.Data()),
                Form("MD_shared_%s_Width", res.Data()),
                Form("MD_shared_%s_Width1", res.Data())
            };
        } else {
            altWidth = {
                "MD_shared_a2_1320_1700_Width1",
                Form("MD_shared_%s_Width1", res.Data()),
                Form("MD_shared_%s_Width", res.Data()),
                "MD_shared_a2_1320_Width1"
            };
        }
        for (const auto& alt : altWidth) {
            if (t->GetLeaf(alt)) {
                widthVal = t->GetLeaf(alt)->GetValue();
                TString altErr = alt + "_error";
                if (t->GetLeaf(altErr)) widthErr = t->GetLeaf(altErr)->GetValue();
                break;
            }
        }
    }

    out << "{\n";
    out << "  \"resonance\": \"" << res.Data() << "\",\n";
    out << "  \"mass\": " << std::setprecision(8) << massVal << ",\n";
    out << "  \"mass_error\": " << std::setprecision(8) << massErr << ",\n";
    out << "  \"width\": " << std::setprecision(8) << widthVal << ",\n";
    out << "  \"width_error\": " << std::setprecision(8) << widthErr << ",\n";
    out << "  \"couplings\": {\n";

    // 2. Iterate branches to find MD_{wave}_{resonance}_k or _k1/_k2
    TObjArray* branches = t->GetListOfBranches();
    bool first = true;

    for (int i = 0; i < branches->GetEntries(); ++i) {
        TString bName = branches->At(i)->GetName();
        if (bName.EndsWith("_error")) continue;
        if (!bName.BeginsWith("MD_")) continue;

        TString kSuffix = "";
        if (bName.EndsWith("_k")) {
            kSuffix = "_k";
        } else if (bName.EndsWith("_k1") && (res == "a2_1320" || res.EndsWith("1320") || res.BeginsWith("res1"))) {
            kSuffix = "_k1";
        } else if (bName.EndsWith("_k2") && (res == "a2_1700" || res.EndsWith("1700") || res.BeginsWith("res2"))) {
            kSuffix = "_k2";
        } else if (bName.EndsWith("_k1")) {
            kSuffix = "_k1";
        }

        if (kSuffix.Length() == 0) continue;

        TString sub = bName;
        sub.ReplaceAll("MD_", "");
        sub.ReplaceAll(kSuffix, "");
        // Now sub is {amp}_{res}
        TString waveName = "";
        if (sub.EndsWith(Form("_%s", res.Data()))) {
            waveName = sub;
            waveName.ReplaceAll(Form("_%s", res.Data()), "");
        } else if (sub.EndsWith("_a2_1320_1700")) {
            waveName = sub;
            waveName.ReplaceAll("_a2_1320_1700", "");
        } else if (sub.BeginsWith("a_") || sub.BeginsWith("b_")) {
            waveName = sub;
        }

        if (waveName.Length() > 0) {
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

    out << "\n  }\n}\n";
    out.close();
    f->Close();
    delete f;

    std::cout << "[OK] Exported resonance couplings to: " << outJsonPath << std::endl;
}


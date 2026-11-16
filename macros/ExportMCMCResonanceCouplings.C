// Helper macro to export MCMC sample resonance couplings and parameters to JSON for cross-section evaluation
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

void ExportMCMCResonanceCouplings(
    const char* rootFilePath,
    const char* resonanceName,
    const char* outJsonPath,
    const char* preferredTreeName = "mcmc_params"
) {
    TFile* f = TFile::Open(rootFilePath, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open ROOT file " << rootFilePath << std::endl;
        if (f) delete f;
        return;
    }

    TTree* t = nullptr;
    if (preferredTreeName && strlen(preferredTreeName) > 0) {
        t = (TTree*)f->Get(preferredTreeName);
    }
    if (!t) t = (TTree*)f->Get("mcmc_params");
    if (!t) t = (TTree*)f->Get("mass_dependent_params");
    if (!t) t = (TTree*)f->Get("mcmc_result");
    if (!t) t = (TTree*)f->Get("central_params");
    if (!t) t = (TTree*)f->Get("all_params");

    if (!t) {
        std::cerr << "Error: no valid parameter tree found in " << rootFilePath << std::endl;
        f->Close();
        delete f;
        return;
    }

    Long64_t nEntries = t->GetEntries();
    if (nEntries <= 0) {
        std::cerr << "Error: tree '" << t->GetName() << "' has no entries in " << rootFilePath << std::endl;
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

    TString res(resonanceName);

    // Identify branches
    TString massBranch = Form("MD_shared_%s_Mass", res.Data());
    TString widthBranch = Form("MD_shared_%s_Width", res.Data());

    bool hasMass = (t->GetLeaf(massBranch) != nullptr);
    bool hasWidth = (t->GetLeaf(widthBranch) != nullptr);

    // Fallbacks for two-breit-wigner or alternate naming schemes
    if (!hasMass) {
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
                massBranch = alt;
                hasMass = true;
                break;
            }
        }
    }
    if (!hasWidth) {
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
                widthBranch = alt;
                hasWidth = true;
                break;
            }
        }
    }

    bool hasValid = (t->GetLeaf("isValid") != nullptr);
    bool hasStatus = (t->GetLeaf("Status") != nullptr);
    bool hasChi2 = (t->GetLeaf("chi2") != nullptr);
    bool hasEntry = (t->GetLeaf("mcmc_entry") != nullptr);
    bool hasWeight = (t->GetLeaf("weight") != nullptr);

    // Identify coupling branches: MD_{wave}_{res}_k or _k1/_k2
    TObjArray* branches = t->GetListOfBranches();
    std::vector<std::string> waveNames;
    std::vector<TString> couplingBranches;
    std::vector<int> reflectivities;

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
            int refl = 0;
            if (waveName.BeginsWith("a_")) refl = 1;
            else if (waveName.BeginsWith("b_")) refl = -1;

            waveNames.push_back(waveName.Data());
            couplingBranches.push_back(bName);
            reflectivities.push_back(refl);
        }
    }

    if (waveNames.empty()) {
        std::cerr << "Error: no coupling branches found for resonance '" << res.Data() << "' in tree '" << t->GetName() << "'" << std::endl;
        out.close();
        f->Close();
        delete f;
        return;
    }

    // Write JSON header
    out << "{\n";
    out << "  \"resonance\": \"" << res.Data() << "\",\n";
    out << "  \"tree_name\": \"" << t->GetName() << "\",\n";
    out << "  \"entries_count\": " << nEntries << ",\n";
    
    // Wave names
    out << "  \"waves\": [";
    for (size_t i = 0; i < waveNames.size(); ++i) {
        if (i > 0) out << ", ";
        out << "\"" << waveNames[i] << "\"";
    }
    out << "],\n";

    // Reflectivities
    out << "  \"reflectivities\": [";
    for (size_t i = 0; i < reflectivities.size(); ++i) {
        if (i > 0) out << ", ";
        out << reflectivities[i];
    }
    out << "],\n";

    // Columnar arrays for entry, is_valid, status, chi2, mass, width
    std::vector<Long64_t> vecEntries;
    std::vector<int> vecValid;
    std::vector<int> vecStatus;
    std::vector<double> vecChi2;
    std::vector<double> vecMass;
    std::vector<double> vecWidth;
    std::vector<double> vecWeights;
    std::map<std::string, std::vector<double>> mapCouplings;
    for (const auto& w : waveNames) {
        mapCouplings[w] = std::vector<double>();
    }

    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);

        Long64_t entryIdx = hasEntry ? (Long64_t)t->GetLeaf("mcmc_entry")->GetValue() : i;
        int isValidVal = hasValid ? (int)t->GetLeaf("isValid")->GetValue() : 1;
        int statusVal = hasStatus ? (int)t->GetLeaf("Status")->GetValue() : 0;
        double chi2Val = hasChi2 ? t->GetLeaf("chi2")->GetValue() : 0.0;
        double massVal = hasMass ? t->GetLeaf(massBranch)->GetValue() : 0.0;
        double widthVal = hasWidth ? t->GetLeaf(widthBranch)->GetValue() : 0.0;
        double weightVal = hasWeight ? t->GetLeaf("weight")->GetValue() : 1.0;

        vecEntries.push_back(entryIdx);
        vecValid.push_back(isValidVal);
        vecStatus.push_back(statusVal);
        vecChi2.push_back(chi2Val);
        vecMass.push_back(massVal);
        vecWidth.push_back(widthVal);
        vecWeights.push_back(weightVal);

        for (size_t w = 0; w < waveNames.size(); ++w) {
            double kVal = t->GetLeaf(couplingBranches[w])->GetValue();
            mapCouplings[waveNames[w]].push_back(kVal);
        }
    }

    // Output is_valid
    out << "  \"is_valid\": [";
    for (size_t i = 0; i < vecValid.size(); ++i) {
        if (i > 0) out << ", ";
        out << vecValid[i];
    }
    out << "],\n";

    // Output status
    out << "  \"status\": [";
    for (size_t i = 0; i < vecStatus.size(); ++i) {
        if (i > 0) out << ", ";
        out << vecStatus[i];
    }
    out << "],\n";

    // Output weights
    out << "  \"weights\": [";
    for (size_t i = 0; i < vecWeights.size(); ++i) {
        if (i > 0) out << ", ";
        out << std::setprecision(6) << vecWeights[i];
    }
    out << "],\n";

    // Output chi2
    out << "  \"chi2\": [";
    for (size_t i = 0; i < vecChi2.size(); ++i) {
        if (i > 0) out << ", ";
        out << std::setprecision(6) << vecChi2[i];
    }
    out << "],\n";

    // Output mass
    out << "  \"mass\": [";
    for (size_t i = 0; i < vecMass.size(); ++i) {
        if (i > 0) out << ", ";
        out << std::setprecision(6) << vecMass[i];
    }
    out << "],\n";

    // Output width
    out << "  \"width\": [";
    for (size_t i = 0; i < vecWidth.size(); ++i) {
        if (i > 0) out << ", ";
        out << std::setprecision(6) << vecWidth[i];
    }
    out << "],\n";

    // Output couplings
    out << "  \"couplings\": {\n";
    for (size_t w = 0; w < waveNames.size(); ++w) {
        const std::string& wName = waveNames[w];
        out << "    \"" << wName << "\": [";
        const auto& vals = mapCouplings[wName];
        for (size_t i = 0; i < vals.size(); ++i) {
            if (i > 0) out << ", ";
            out << std::setprecision(8) << vals[i];
        }
        out << "]";
        if (w + 1 < waveNames.size()) out << ",\n";
        else out << "\n";
    }
    out << "  }\n";

    out << "}\n";
    out.close();
    f->Close();
    delete f;

    std::cout << "[OK] Exported MCMC resonance couplings (" << nEntries << " entries) to: " << outJsonPath << std::endl;
}


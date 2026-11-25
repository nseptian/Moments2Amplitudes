// Helper macro to export unnormalized MCMC moments and fit results moments to JSON for plotting
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>

void ExportMomentComparison(
    const char* unnormMomentsFile,
    const char* bestFitFilesList,  // Comma-separated list of ROOT files for best fits (or nominal best)
    const char* outJsonFile,
    const char* filterMoments = "" // Comma-separated moment names, or empty for all
) {
    std::ofstream out(outJsonFile);
    if (!out.is_open()) {
        std::cerr << "Error: cannot open output file " << outJsonFile << std::endl;
        return;
    }

    out << "{\n";

    // 1. Export Data Moments from UnnormalizedMoments_vs_Mass.root
    out << "  \"data\": {\n";
    TFile* f_data = TFile::Open(unnormMomentsFile);
    bool firstMoment = true;

    if (f_data && !f_data->IsZombie()) {
        TList* keys = f_data->GetListOfKeys();
        for (int i = 0; i < keys->GetSize(); ++i) {
            TKey* k = (TKey*)keys->At(i);
            TString kname = k->GetName();
            if (kname.BeginsWith("gr_unnormalized_moments_")) {
                TString sub = kname;
                sub.ReplaceAll("gr_unnormalized_moments_", "");
                sub.ReplaceAll("_vs_mass", "");
                TString momName = Form("H_%s", sub.Data());

                TGraphErrors* gr = dynamic_cast<TGraphErrors*>(f_data->Get(kname));
                if (!gr) continue;

                if (!firstMoment) out << ",\n";
                firstMoment = false;

                out << "    \"" << momName << "\": {\n";
                out << "      \"mass\": [";
                for (int p = 0; p < gr->GetN(); ++p) {
                    if (p > 0) out << ", ";
                    out << std::setprecision(6) << gr->GetPointX(p);
                }
                out << "],\n";

                out << "      \"val\": [";
                for (int p = 0; p < gr->GetN(); ++p) {
                    if (p > 0) out << ", ";
                    out << std::setprecision(6) << gr->GetPointY(p);
                }
                out << "],\n";

                out << "      \"err\": [";
                for (int p = 0; p < gr->GetN(); ++p) {
                    if (p > 0) out << ", ";
                    out << std::setprecision(6) << gr->GetErrorY(p);
                }
                out << "]\n    }";
            }
        }
        f_data->Close();
        delete f_data;
    }
    out << "\n  },\n";

    // 2. Export Model Fits Moments
    out << "  \"fits\": [\n";

    TString filesStr(bestFitFilesList);
    TObjArray* fileTokens = filesStr.Tokenize(",");
    bool firstFit = true;

    for (int fIdx = 0; fIdx < fileTokens->GetEntries(); ++fIdx) {
        TString fPath = ((TObjString*)fileTokens->At(fIdx))->GetString().Strip(TString::kBoth);
        if (fPath.IsNull()) continue;

        TFile* f_fit = TFile::Open(fPath);
        if (!f_fit || f_fit->IsZombie()) {
            if (f_fit) delete f_fit;
            continue;
        }

        TTree* chi2Tree = (TTree*)f_fit->Get("chi2");
        double chi2Val = 0.0;
        int statusVal = 0;
        bool isValidVal = true;
        bool atLimitVal = false;
        if (chi2Tree && chi2Tree->GetEntries() > 0) {
            if (chi2Tree->GetBranch("chi2")) {
                chi2Tree->SetBranchAddress("chi2", &chi2Val);
            } else if (chi2Tree->GetBranch("Chi2")) {
                chi2Tree->SetBranchAddress("Chi2", &chi2Val);
            }
            if (chi2Tree->GetBranch("Status")) {
                chi2Tree->SetBranchAddress("Status", &statusVal);
            } else if (chi2Tree->GetBranch("status")) {
                chi2Tree->SetBranchAddress("status", &statusVal);
            }
            if (chi2Tree->GetBranch("isValid")) {
                chi2Tree->SetBranchAddress("isValid", &isValidVal);
            } else if (chi2Tree->GetBranch("valid")) {
                chi2Tree->SetBranchAddress("valid", &isValidVal);
            }
            if (chi2Tree->GetBranch("atLimit")) {
                chi2Tree->SetBranchAddress("atLimit", &atLimitVal);
            }
            chi2Tree->GetEntry(0);
        }

        // Only include converged fits (Status == 0, isValid == true, and not at parameter limits)
        if (statusVal != 0 || !isValidVal || atLimitVal) {
            f_fit->Close();
            delete f_fit;
            continue;
        }

        TTree* resTree = (TTree*)f_fit->Get("result");
        if (!resTree) {
            f_fit->Close();
            delete f_fit;
            continue;
        }

        if (!firstFit) out << ",\n";
        firstFit = false;

        out << "    {\n";
        out << "      \"file\": \"" << fPath << "\",\n";
        out << "      \"chi2\": " << chi2Val << ",\n";
        out << "      \"status\": " << statusVal << ",\n";
        out << "      \"moments\": {\n";

        // Find all H_ branches
        TObjArray* branches = resTree->GetListOfBranches();
        std::vector<TString> hBranches;
        for (int b = 0; b < branches->GetEntries(); ++b) {
            TString bName = branches->At(b)->GetName();
            if (bName.BeginsWith("H_")) {
                hBranches.push_back(bName);
            }
        }

        Double_t massBin = 0.0;
        resTree->SetBranchAddress("mass_bin", &massBin);
        std::map<TString, Double_t> branchVals;
        for (const auto& bName : hBranches) {
            branchVals[bName] = 0.0;
            resTree->SetBranchAddress(bName, &branchVals[bName]);
        }

        // Read all entries
        std::vector<Double_t> masses;
        std::map<TString, std::vector<Double_t>> fitMoments;
        for (Long64_t entry = 0; entry < resTree->GetEntries(); ++entry) {
            resTree->GetEntry(entry);
            masses.push_back(massBin);
            for (const auto& bName : hBranches) {
                // Model moments are divided by 2 to match BruFit unnormalized moment convention
                fitMoments[bName].push_back(branchVals[bName] / 2.0);
            }
        }

        bool firstBranch = true;
        for (const auto& bName : hBranches) {
            if (!firstBranch) out << ",\n";
            firstBranch = false;

            out << "        \"" << bName << "\": {\n";
            out << "          \"mass\": [";
            for (size_t p = 0; p < masses.size(); ++p) {
                if (p > 0) out << ", ";
                out << std::setprecision(6) << masses[p];
            }
            out << "],\n";

            out << "          \"val\": [";
            for (size_t p = 0; p < fitMoments[bName].size(); ++p) {
                if (p > 0) out << ", ";
                out << std::setprecision(6) << fitMoments[bName][p];
            }
            out << "]\n        }";
        }

        out << "\n      }\n    }";
        f_fit->Close();
        delete f_fit;
    }
    delete fileTokens;

    out << "\n  ]\n}\n";
    out.close();
    std::cout << "[OK] Exported moment comparison data to: " << outJsonFile << std::endl;
}


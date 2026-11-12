// Helper macro to export nominal and MCMC amplitude branches to JSON for plotting
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>

void ExportAmplitudes(
    const char* nominalRootFile,
    const char* mcmcRootFile,
    const char* outJsonFile
) {
    std::ofstream out(outJsonFile);
    if (!out.is_open()) {
        std::cerr << "Error: cannot open output JSON file " << outJsonFile << std::endl;
        return;
    }

    TFile* f_nom = TFile::Open(nominalRootFile);
    if (!f_nom || f_nom->IsZombie()) {
        std::cerr << "Error: cannot open nominal ROOT file " << nominalRootFile << std::endl;
        out << "{}\n";
        return;
    }

    TTree* t_nom = dynamic_cast<TTree*>(f_nom->Get("result"));
    if (!t_nom) {
        t_nom = dynamic_cast<TTree*>(f_nom->Get("central_result"));
    }
    if (!t_nom) {
        std::cerr << "Error: TTree 'result' or 'central_result' not found in " << nominalRootFile << std::endl;
        f_nom->Close();
        out << "{}\n";
        return;
    }

    const Long64_t nEntries = t_nom->GetEntries();
    double mass_val = 0.0;
    t_nom->SetBranchAddress("mass_bin", &mass_val);

    std::vector<std::string> ampNames;
    std::map<std::string, double> branchBuffers;

    TObjArray* branches = t_nom->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* br = dynamic_cast<TBranch*>(branches->At(i));
        if (!br) continue;
        TString bName = br->GetName();
        if (bName == "mass_bin") continue;
        if (bName.BeginsWith("H_")) continue; // Skip moment projections

        std::string nameStr = bName.Data();
        ampNames.push_back(nameStr);
        branchBuffers[nameStr] = 0.0;
        t_nom->SetBranchAddress(bName, &branchBuffers[nameStr]);
    }

    std::vector<double> massBins;
    std::map<std::string, std::vector<double>> nominalValues;
    for (const auto& name : ampNames) {
        nominalValues[name] = std::vector<double>();
    }

    for (Long64_t e = 0; e < nEntries; ++e) {
        t_nom->GetEntry(e);
        massBins.push_back(mass_val);
        for (const auto& name : ampNames) {
            nominalValues[name].push_back(branchBuffers[name]);
        }
    }
    f_nom->Close();

    // Start JSON output
    out << "{\n";
    out << "  \"masses\": [";
    for (size_t i = 0; i < massBins.size(); ++i) {
        if (i > 0) out << ", ";
        out << std::setprecision(6) << massBins[i];
    }
    out << "],\n";

    out << "  \"nominal\": {\n";
    for (size_t i = 0; i < ampNames.size(); ++i) {
        const std::string& name = ampNames[i];
        out << "    \"" << name << "\": [";
        const auto& vals = nominalValues[name];
        for (size_t j = 0; j < vals.size(); ++j) {
            if (j > 0) out << ", ";
            out << std::setprecision(6) << vals[j];
        }
        out << "]" << (i + 1 < ampNames.size() ? ",\n" : "\n");
    }
    out << "  },\n";

    // MCMC samples if available
    out << "  \"mcmc\": {\n";
    if (mcmcRootFile && std::string(mcmcRootFile) != "" && std::string(mcmcRootFile) != "none") {
        TFile* f_mcmc = TFile::Open(mcmcRootFile);
        if (f_mcmc && !f_mcmc->IsZombie()) {
            TTree* t_mcmc = dynamic_cast<TTree*>(f_mcmc->Get("mcmc_result"));
            if (!t_mcmc) {
                t_mcmc = dynamic_cast<TTree*>(f_mcmc->Get("result"));
            }
            if (t_mcmc) {
                double mcmc_mass = 0.0;
                t_mcmc->SetBranchAddress("mass_bin", &mcmc_mass);
                std::map<std::string, double> mcmcBuffers;
                std::vector<std::string> validMcmcAmps;
                for (const auto& name : ampNames) {
                    if (t_mcmc->GetBranch(name.c_str())) {
                        validMcmcAmps.push_back(name);
                        mcmcBuffers[name] = 0.0;
                        t_mcmc->SetBranchAddress(name.c_str(), &mcmcBuffers[name]);
                    }
                }

                // Map: ampName -> mass_bin -> list of samples
                std::map<std::string, std::map<double, std::vector<double>>> samps;
                const Long64_t nMcmc = t_mcmc->GetEntries();
                for (Long64_t e = 0; e < nMcmc; ++e) {
                    t_mcmc->GetEntry(e);
                    for (const auto& name : validMcmcAmps) {
                        samps[name][mcmc_mass].push_back(mcmcBuffers[name]);
                    }
                }

                for (size_t i = 0; i < validMcmcAmps.size(); ++i) {
                    const std::string& name = validMcmcAmps[i];
                    out << "    \"" << name << "\": {\n";
                    size_t mCount = 0;
                    for (const auto& [m, sList] : samps[name]) {
                        if (mCount > 0) out << ",\n";
                        out << "      \"" << std::fixed << std::setprecision(4) << m << "\": [";
                        for (size_t sIdx = 0; sIdx < sList.size(); ++sIdx) {
                            if (sIdx > 0) out << ", ";
                            out << std::setprecision(6) << sList[sIdx];
                        }
                        out << "]";
                        mCount++;
                    }
                    out << "\n    }" << (i + 1 < validMcmcAmps.size() ? ",\n" : "\n");
                }
            }

            out << "  }";

            // Also export mcmc_params tree if present
            TTree* t_params = dynamic_cast<TTree*>(f_mcmc->Get("mcmc_params"));
            if (t_params) {
                out << ",\n  \"mcmc_params\": {\n";
                TObjArray* pBranches = t_params->GetListOfBranches();
                std::vector<std::string> paramNames;
                std::map<std::string, double> pBuffers;
                for (int i = 0; i < pBranches->GetEntries(); ++i) {
                    TBranch* br = dynamic_cast<TBranch*>(pBranches->At(i));
                    if (!br) continue;
                    TString bName = br->GetName();
                    if (bName == "mcmc_entry" || bName == "Status" || bName == "isValid" || bName == "isValidError") continue;
                    std::string pStr = bName.Data();
                    paramNames.push_back(pStr);
                    pBuffers[pStr] = 0.0;
                    t_params->SetBranchAddress(bName, &pBuffers[pStr]);
                }
                std::map<std::string, std::vector<double>> paramVals;
                const Long64_t nParEntries = t_params->GetEntries();
                for (Long64_t e = 0; e < nParEntries; ++e) {
                    t_params->GetEntry(e);
                    for (const auto& p : paramNames) {
                        paramVals[p].push_back(pBuffers[p]);
                    }
                }
                for (size_t i = 0; i < paramNames.size(); ++i) {
                    const std::string& p = paramNames[i];
                    out << "    \"" << p << "\": [";
                    const auto& vec = paramVals[p];
                    for (size_t j = 0; j < vec.size(); ++j) {
                        if (j > 0) out << ", ";
                        out << std::setprecision(6) << vec[j];
                    }
                    out << "]" << (i + 1 < paramNames.size() ? ",\n" : "\n");
                }
                out << "  }\n";
            } else {
                out << "\n";
            }

            f_mcmc->Close();
        } else {
            out << "  }\n";
        }
    } else {
        out << "  }\n";
    }
    out << "}\n";
}

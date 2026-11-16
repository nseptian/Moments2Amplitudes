// Helper macro to export nominal and MCMC amplitude branches to JSON for plotting
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <iomanip>

void ExportAmplitudes(
    const char* nominalRootFile,
    const char* mcmcRootFile,
    const char* outJsonFile,
    const char* unnormMomentsFile = ""
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

    bool hasH0Model = (t_nom->GetBranch("H_0_0_0") != nullptr);
    double h0ModelBuf = 0.0;
    if (hasH0Model) {
        t_nom->SetBranchAddress("H_0_0_0", &h0ModelBuf);
    }

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
    std::vector<double> h0ModelValues;
    std::map<std::string, std::vector<double>> nominalValues;
    for (const auto& name : ampNames) {
        nominalValues[name] = std::vector<double>();
    }

    for (Long64_t e = 0; e < nEntries; ++e) {
        t_nom->GetEntry(e);
        massBins.push_back(mass_val);
        if (hasH0Model) {
            h0ModelValues.push_back(h0ModelBuf);
        }
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
            // Identify converged MCMC entries (Status == 0 and isValid == true)
            std::set<Long64_t> convergedEntries;
            bool hasConvCheck = false;
            TTree* t_conv_check = dynamic_cast<TTree*>(f_mcmc->Get("mcmc_params"));
            if (!t_conv_check) t_conv_check = dynamic_cast<TTree*>(f_mcmc->Get("mcmc_chi2"));
            if (t_conv_check) {
                hasConvCheck = true;
                Long64_t m_entry = 0;
                int m_status = 0;
                bool m_valid = true;
                bool has_entry_b = (t_conv_check->GetBranch("mcmc_entry") != nullptr);
                bool has_status_b = (t_conv_check->GetBranch("Status") != nullptr);
                bool has_valid_b = (t_conv_check->GetBranch("isValid") != nullptr);

                if (has_entry_b) t_conv_check->SetBranchAddress("mcmc_entry", &m_entry);
                if (has_status_b) t_conv_check->SetBranchAddress("Status", &m_status);
                if (has_valid_b) t_conv_check->SetBranchAddress("isValid", &m_valid);

                for (Long64_t e = 0; e < t_conv_check->GetEntries(); ++e) {
                    t_conv_check->GetEntry(e);
                    bool ok = true;
                    if (has_status_b && m_status != 0) ok = false;
                    if (has_valid_b && !m_valid) ok = false;
                    if (ok) {
                        convergedEntries.insert(has_entry_b ? m_entry : e);
                    }
                }
            }

            TTree* t_mcmc = dynamic_cast<TTree*>(f_mcmc->Get("mcmc_result"));
            if (!t_mcmc) {
                t_mcmc = dynamic_cast<TTree*>(f_mcmc->Get("result"));
            }
            if (t_mcmc) {
                double mcmc_mass = 0.0;
                t_mcmc->SetBranchAddress("mass_bin", &mcmc_mass);
                Long64_t mcmc_entry = 0;
                bool has_mcmc_entry = (t_mcmc->GetBranch("mcmc_entry") != nullptr);
                if (has_mcmc_entry) t_mcmc->SetBranchAddress("mcmc_entry", &mcmc_entry);
                double mcmc_weight = 1.0;
                bool has_mcmc_weight = (t_mcmc->GetBranch("weight") != nullptr);
                if (has_mcmc_weight) t_mcmc->SetBranchAddress("weight", &mcmc_weight);

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
                std::map<double, std::vector<double>> weightList;
                const Long64_t nMcmc = t_mcmc->GetEntries();
                const size_t nMassBins = (massBins.empty() ? 1 : massBins.size());
                for (Long64_t e = 0; e < nMcmc; ++e) {
                    t_mcmc->GetEntry(e);
                    if (hasConvCheck) {
                        Long64_t cur_id = has_mcmc_entry ? mcmc_entry : (e / nMassBins);
                        if (convergedEntries.find(cur_id) == convergedEntries.end()) {
                            continue; // Only include converged MCMC fits!
                        }
                    }
                    if (has_mcmc_weight) {
                        weightList[mcmc_mass].push_back(mcmc_weight);
                    }
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

                out << "  }";

                if (has_mcmc_weight && !weightList.empty()) {
                    out << ",\n  \"mcmc_weights\": {\n";
                    for (size_t i = 0; i < validMcmcAmps.size(); ++i) {
                        const std::string& name = validMcmcAmps[i];
                        out << "    \"" << name << "\": {\n";
                        size_t mCount = 0;
                        for (const auto& [m, wList] : weightList) {
                            if (mCount > 0) out << ",\n";
                            out << "      \"" << std::fixed << std::setprecision(4) << m << "\": [";
                            for (size_t sIdx = 0; sIdx < wList.size(); ++sIdx) {
                                if (sIdx > 0) out << ", ";
                                out << std::setprecision(6) << wList[sIdx];
                            }
                            out << "]";
                            mCount++;
                        }
                        out << "\n    }" << (i + 1 < validMcmcAmps.size() ? ",\n" : "\n");
                    }
                    out << "  }";
                }
            } else {
                out << "  }";
            }

            // Also export mcmc_params tree if present (converged fits only)
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

                Long64_t par_entry = 0;
                int par_status = 0;
                bool par_valid = true;
                bool has_par_entry = (t_params->GetBranch("mcmc_entry") != nullptr);
                bool has_par_status = (t_params->GetBranch("Status") != nullptr);
                bool has_par_valid = (t_params->GetBranch("isValid") != nullptr);

                if (has_par_entry) t_params->SetBranchAddress("mcmc_entry", &par_entry);
                if (has_par_status) t_params->SetBranchAddress("Status", &par_status);
                if (has_par_valid) t_params->SetBranchAddress("isValid", &par_valid);

                std::map<std::string, std::vector<double>> paramVals;
                const Long64_t nParEntries = t_params->GetEntries();
                for (Long64_t e = 0; e < nParEntries; ++e) {
                    t_params->GetEntry(e);
                    bool ok = true;
                    if (has_par_status && par_status != 0) ok = false;
                    if (has_par_valid && !par_valid) ok = false;
                    if (!ok) continue; // Only include converged fits!

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

    if (hasH0Model) {
        out << ",\n  \"h0_model\": {\n";
        out << "    \"mass\": [";
        for (size_t i = 0; i < massBins.size(); ++i) {
            if (i > 0) out << ", ";
            out << std::setprecision(6) << massBins[i];
        }
        out << "],\n    \"val\": [";
        for (size_t i = 0; i < h0ModelValues.size(); ++i) {
            if (i > 0) out << ", ";
            out << std::setprecision(6) << h0ModelValues[i];
        }
        out << "]\n  }";
    }

    if (unnormMomentsFile && strlen(unnormMomentsFile) > 0) {
        TFile* f_unnorm = TFile::Open(unnormMomentsFile);
        if (f_unnorm && !f_unnorm->IsZombie()) {
            TGraphErrors* gr = dynamic_cast<TGraphErrors*>(f_unnorm->Get("gr_unnormalized_moments_0_0_0_vs_mass"));
            if (!gr) {
                TList* ukeys = f_unnorm->GetListOfKeys();
                for (int i = 0; i < ukeys->GetSize(); ++i) {
                    TKey* uk = (TKey*)ukeys->At(i);
                    TString ukname = uk->GetName();
                    if (ukname.Contains("unnormalized_moments_0_0_0") || ukname.Contains("H_0_0_0")) {
                        gr = dynamic_cast<TGraphErrors*>(f_unnorm->Get(ukname));
                        if (gr) break;
                    }
                }
            }
            if (gr) {
                out << ",\n  \"h0_data\": {\n";
                out << "    \"mass\": [";
                for (int p = 0; p < gr->GetN(); ++p) {
                    if (p > 0) out << ", ";
                    out << std::setprecision(6) << gr->GetPointX(p);
                }
                out << "],\n    \"val\": [";
                for (int p = 0; p < gr->GetN(); ++p) {
                    if (p > 0) out << ", ";
                    out << std::setprecision(6) << gr->GetPointY(p);
                }
                out << "],\n    \"err\": [";
                for (int p = 0; p < gr->GetN(); ++p) {
                    if (p > 0) out << ", ";
                    out << std::setprecision(6) << gr->GetErrorY(p);
                }
                out << "]\n  }";
            }
            f_unnorm->Close();
        }
    }

    out << "\n}\n";
}

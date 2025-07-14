#include "MassDependentEquationSolver.h"
#include "TObjString.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace m2pw {

    MassDependentEquationSolver::MassDependentEquationSolver(const Setup& setup, 
                                                           std::vector<double> mass_bins, 
                                                           double l_max, double noise, 
                                                           std::vector<TString> noUse) 
        : massBins_(std::move(mass_bins)), massDependenceConfig_(MassDependenceConfig()), hMomentsConfig_(HMomentsConfig()) {
        
        InitializeMassBins(massBins_);
        InitializeMDFunctions();
        
        // Initialize equations for each mass bin
        const auto validEqn = [&noUse](const TString& name) {
            // if (name.BeginsWith("a_")) return false; // Skip amplitude equations
            if (name == "normalise") return false;
            return std::none_of(noUse.begin(), noUse.end(), 
                               [&name](const TString& match) { return name.Contains(match); });
        };

        for (const double& mass_bin : massBins_) {
            pars_[mass_bin] = ParameterHelper(setup);
            
            for (const auto& form : setup.ParameterFormulas()) {
                if (!validEqn(form->GetName())) continue;
                
                if (auto* var = dynamic_cast<RooFormulaVar*>(form)) {
                    eqns_[mass_bin].emplace_back(var, &pars_[mass_bin], noise);
                }
            }

            // Initialize equations
            for (auto& eqn : eqns_[mass_bin]) {
                eqn.FindL();
                eqn.DoEval(pars_[mass_bin].CurrentVals());
            }
        }

        // Initialize parameter name to index mapping
        int index = 0;
        for (int l = 0; l <= static_cast<int>(l_max); ++l) {
            for (int m = -l; m <= l; ++m) {
                parNameToIndex_[Form("a_%d_%d", l, m)] = index++;
                parNameToIndex_[Form("b_%d_%d", l, m)] = index++;
            }
        }
    }

    MassDependentEquationSolver::MassDependentEquationSolver(const Setup& setup, 
                                                           std::vector<double> mass_bins, 
                                                           double l_max, double noise, 
                                                           std::vector<TString> noUse,
                                                           const MassDependenceConfig& config) 
        : massBins_(std::move(mass_bins)), massDependenceConfig_(config), hMomentsConfig_(HMomentsConfig()) {
        
        InitializeMassBins(massBins_);
        InitializeMDFunctions();
        
        // Initialize equations for each mass bin
        const auto validEqn = [&noUse](const TString& name) {
            // if (name.BeginsWith("a_")) return false; // Skip amplitude equations
            if (name == "normalise") return false;
            return std::none_of(noUse.begin(), noUse.end(), 
                               [&name](const TString& match) { return name.Contains(match); });
        };

        for (const double& mass_bin : massBins_) {
            pars_[mass_bin] = ParameterHelper(setup);
            
            for (const auto& form : setup.ParameterFormulas()) {
                if (!validEqn(form->GetName())) continue;
                
                if (auto* var = dynamic_cast<RooFormulaVar*>(form)) {
                    eqns_[mass_bin].emplace_back(var, &pars_[mass_bin], noise);
                }
            }

            // Initialize equations
            for (auto& eqn : eqns_[mass_bin]) {
                eqn.FindL();
                eqn.DoEval(pars_[mass_bin].CurrentVals());
            }
        }

        // Initialize parameter name to index mapping
        int index = 0;
        for (int l = 0; l <= static_cast<int>(l_max); ++l) {
            for (int m = -l; m <= l; ++m) {
                parNameToIndex_[Form("a_%d_%d", l, m)] = index++;
                parNameToIndex_[Form("b_%d_%d", l, m)] = index++;
            }
        }
    }

    MassDependentEquationSolver::MassDependentEquationSolver(const Setup& setup, 
                                                           std::vector<double> mass_bins, 
                                                           double l_max, double noise, 
                                                           std::vector<TString> noUse,
                                                           const MassDependenceConfig& massDepConfig,
                                                           const HMomentsConfig& hMomentsConfig) 
        : massBins_(std::move(mass_bins)), massDependenceConfig_(massDepConfig), hMomentsConfig_(hMomentsConfig) {
        
        InitializeMassBins(massBins_);
        InitializeMDFunctions();
        
        // Initialize equations for each mass bin
        const auto validEqn = [&noUse](const TString& name) {
            // if (name.BeginsWith("a_")) return false; // Skip amplitude equations
            if (name == "normalise") return false;
            return std::none_of(noUse.begin(), noUse.end(), 
                               [&name](const TString& match) { return name.Contains(match); });
        };

        for (const double& mass_bin : massBins_) {
            pars_[mass_bin] = ParameterHelper(setup);
            
            for (const auto& form : setup.ParameterFormulas()) {
                if (!validEqn(form->GetName())) continue;
                
                if (auto* var = dynamic_cast<RooFormulaVar*>(form)) {
                    eqns_[mass_bin].emplace_back(var, &pars_[mass_bin], noise);
                }
            }

            // Initialize equations
            for (auto& eqn : eqns_[mass_bin]) {
                eqn.FindL();
                eqn.DoEval(pars_[mass_bin].CurrentVals());
            }
        }

        // Initialize parameter name to index mapping - FILTER BASED ON H MOMENTS CONFIG
        int index = 0;
        for (int l = 0; l <= static_cast<int>(l_max); ++l) {
            // Only include parameters for L values that are needed for the H moments
            if (!ParameterNeededForHMoments(l, hMomentsConfig)) {
                std::cout << "Skipping parameters for l=" << l << " (not needed for selected H moments)" << std::endl;
                continue;
            }
            
            for (int m = -l; m <= l; ++m) {
                parNameToIndex_[Form("a_%d_%d", l, m)] = index++;
                parNameToIndex_[Form("b_%d_%d", l, m)] = index++;
                std::cout << "Including parameters: a_" << l << "_" << m << " and b_" << l << "_" << m << std::endl;
            }
        }
        
        std::cout << "Total parameters to fit: " << parNameToIndex_.size() << std::endl;
    }

    void MassDependentEquationSolver::InitializeMassBins(const std::vector<double>& mass_bins) {
        if (mass_bins.empty()) {
            throw std::invalid_argument("Mass bins cannot be empty");
        }
        // Note: std::map doesn't have reserve(), it automatically manages memory
    }

    void MassDependentEquationSolver::InitializeMDFunctions() {
        if (massBins_.size() < 2) {
            throw std::invalid_argument("At least two mass bins required");
        }
        
        const double bin_width = massBins_[1] - massBins_[0];
        
        // Initialize mass-dependent functions based on waves in configuration
        massDepFuncs_.clear();
        
        // Collect all unique wave names from the configuration
        std::set<TString> allWaves;
        for (const auto& [l_value, waves] : massDependenceConfig_.massDependentWaves) {
            for (const TString& wave : waves) {
                allWaves.insert(wave);
            }
        }
        
        // Initialize functions for each unique wave
        int funcIndex = 0;
        waveToFunctionIndex_.clear();
        
        for (const TString& wave : allWaves) {
            if (wave == "a2_1320") {
                massDepFuncs_.emplace_back(static_cast<int>(massBins_.size()), 
                                          massBins_.front(), bin_width, 0);  // Breit-Wigner for a2(1320)
                waveToFunctionIndex_[wave] = funcIndex;
                std::cout << "  - Function " << funcIndex << ": " << wave << " Breit-Wigner" << std::endl;
                funcIndex++;
            } else if (wave == "a2_1700") {
                massDepFuncs_.emplace_back(static_cast<int>(massBins_.size()), 
                                          massBins_.front(), bin_width, 1);  // Breit-Wigner for a2(1700)
                waveToFunctionIndex_[wave] = funcIndex;
                std::cout << "  - Function " << funcIndex << ": " << wave << " Breit-Wigner" << std::endl;
                funcIndex++;
            } else if (wave == "a0_980") {
                massDepFuncs_.emplace_back(static_cast<int>(massBins_.size()), 
                                          massBins_.front(), bin_width, 2);  // Flatté for a0(980)
                waveToFunctionIndex_[wave] = funcIndex;
                std::cout << "  - Function " << funcIndex << ": " << wave << " Flatté" << std::endl;
                funcIndex++;
            } else {
                std::cerr << "Warning: Unknown wave " << wave << " in configuration. Skipping." << std::endl;
            }
        }
        
        std::cout << "Initialized " << massDepFuncs_.size() << " mass-dependent functions based on configuration." << std::endl;
    }

    void MassDependentEquationSolver::SetMassDependenceConfig(const MassDependenceConfig& config) {
        massDependenceConfig_ = config;
    }

    void MassDependentEquationSolver::SetHMomentsConfig(const HMomentsConfig& config) {
        hMomentsConfig_ = config;
    }

    // Core evaluation method - simplified and optimized
    double MassDependentEquationSolver::DoEval(const double* mass_dep_pars) {
        // Check cache first
        if (IsParameterCached(mass_dep_pars)) {
            return lastChi2_;
        }
        
        double total_chi2 = 0.0;
        
        for (size_t i = 0; i < massBins_.size(); ++i) {
            const double mass_bin = massBins_[i];
            auto& pars = pars_.at(mass_bin);
            
            // Set parameter values for this mass bin
            for (int j = 0; j < pars.Nvars(); ++j) {
                const TString par_name = pars.GetParName(j);
                std::unique_ptr<TObjArray> parts(par_name.Tokenize("_"));
                
                if (parts->GetEntries() < 3) continue;
                
                const TString refl_str = ((TObjString*)parts->At(0))->GetString();
                const TString l_str = ((TObjString*)parts->At(1))->GetString();
                const TString m_str = ((TObjString*)parts->At(2))->GetString();
                
                double value = 0.0;
                
                if (refl_str == "a" || refl_str == "b") {
                    const auto it = parNameToIndex_.find(par_name);
                    if (it != parNameToIndex_.end()) {
                        const double mass_dep_par = mass_dep_pars[it->second];
                        
                        if (l_str == "0" || l_str == "1") {
                            value = mass_dep_par;
                        } else if (l_str == "2") {
                            value = massDepFuncs_[0].GetPWMagnitude(mass_bin, mass_dep_par);
                        }
                    }
                } else if (refl_str == "aphi" || refl_str == "bphi") {
                    TString base_name = refl_str;
                    base_name.ReplaceAll("phi", "");
                    const TString base_par_name = Form("%s_%s_%s", base_name.Data(), l_str.Data(), m_str.Data());
                    
                    const auto it = parNameToIndex_.find(base_par_name);
                    if (it != parNameToIndex_.end()) {
                        const double mass_dep_par = mass_dep_pars[it->second];
                        if (l_str == "2") {
                            value = massDepFuncs_[0].GetPWPhase(mass_bin, mass_dep_par);
                        }
                    }
                }
                
                pars.SetCurrentVal(par_name, value);
            }
            
            // Evaluate chi2 for this mass bin
            total_chi2 += EvaluateChi2ForMassBin(mass_bin, pars);
        }
        
        // Update cache
        lastChi2_ = total_chi2;
        if (cachedParams_.size() != NDim()) {
            cachedParams_.resize(NDim());
        }
        std::copy(mass_dep_pars, mass_dep_pars + NDim(), cachedParams_.begin());
        cacheValid_ = true;
        
        return total_chi2;
    }

    unsigned int MassDependentEquationSolver::NDim() const {
        if (pars_.empty()) return 0;
        
        unsigned int n = 0;
        for (const auto& [mass_bin, pars] : pars_) {
            n += pars.Nvars();
        }
        return n;
    }

    double MassDependentEquationSolver::EvaluateChi2ForMassBin(double massBin, ParameterHelper& pars) const {
        double chi2 = 0.0;
        const auto current_vals = pars.CurrentVals();
        
        int includedEquations = 0;
        int excludedEquations = 0;
        
        for (auto& eqn : eqns_.at(massBin)) {
            // Extract L value from equation name to check if it should be included
            int L_value = ExtractLFromEquationName(eqn.GetName());
            
            // Apply H moments filtering
            bool shouldInclude = true;
            if (L_value >= 0) { // This is an H moment equation
                shouldInclude = hMomentsConfig_.ShouldIncludeL(L_value);
            }
            
            if (shouldInclude) {
                eqn.SetNeedsRecalc();
                double eqn_chi2 = eqn.DoEvalSq(current_vals);
                chi2 += eqn_chi2;
                includedEquations++;
                
                // Optional: detailed output for debugging
                // std::cout << "  Including " << eqn.GetName() << " (L=" << L_value << "): chi2 = " << eqn_chi2 << std::endl;
            } else {
                excludedEquations++;
                // std::cout << "  Excluding " << eqn.GetName() << " (L=" << L_value << ")" << std::endl;
            }
        }
        
        // Print summary for this mass bin
        // std::cout << "  Mass bin " << massBin << ": included " << includedEquations 
        //           << " equations, excluded " << excludedEquations << " equations" << std::endl;
        
        return chi2;
    }

    std::unique_ptr<MassDependentEquationSolver::ParameterInfo> 
    MassDependentEquationSolver::ParseParameterName(const TString& parName) const {
        auto info = std::make_unique<ParameterInfo>(parName);
        
        std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
        if (parts->GetEntries() < 3) {
            return nullptr;
        }
        
        return info;
    }

    // Optimized evaluation with parameter maps
    double MassDependentEquationSolver::DoEval(
        const std::map<TString, std::vector<double>>& massDepPars,
        const std::map<TString, std::map<TString, std::vector<double>>>& massIndepPars) {
        
        double total_chi2 = 0.0;
        
        for (size_t i = 0; i < massBins_.size(); ++i) {
            const double mass_bin = massBins_[i];
            const TString mass_bin_str = TString::Format("%1.2f", mass_bin);
            auto& pars = pars_.at(mass_bin);
            
            // Set parameter values
            for (int j = 0; j < pars.Nvars(); ++j) {
                const TString par_name = pars.GetParName(j);
                const auto param_info = ParseParameterName(par_name);
                
                if (!param_info) continue;
                
                // Extract l value from parameter name
                std::unique_ptr<TObjArray> parts(par_name.Tokenize("_"));
                if (parts->GetEntries() < 3) continue;
                
                TString l_str = ((TObjString*)parts->At(1))->GetString();
                int l_value = l_str.Atoi();
                
                double value = 0.0;
                if (!param_info->isPhase) {
                    // Handle magnitude parameters
                    if (massDependenceConfig_.IsMassDependent(l_value)) {
                        // Check if we need coherent sum or single wave
                        std::vector<TString> waves = massDependenceConfig_.GetWavesForL(l_value);
                        if (waves.size() > 1) {
                            // Multiple waves: use coherent sum
                            value = GetCoherentMagnitudeForL(mass_bin, par_name, massDepPars, l_value);
                        } else if (waves.size() == 1) {
                            // Single wave: direct evaluation
                            value = GetSingleWaveMagnitude(mass_bin, par_name, massDepPars, waves[0]);
                        }
                    } else {
                        // Use mass-independent parameters
                        const auto bin_it = massIndepPars.find(mass_bin_str);
                        if (bin_it != massIndepPars.end()) {
                            const auto par_it = bin_it->second.find(par_name);
                            if (par_it != bin_it->second.end() && !par_it->second.empty()) {
                                value = par_it->second[0];
                            }
                        }
                    }
                } else {
                    // Handle phase parameters
                    TString base_name = param_info->name;
                    base_name.ReplaceAll("phi", "");
                    
                    if (massDependenceConfig_.IsMassDependent(l_value)) {
                        // Check if we need coherent sum or single wave
                        std::vector<TString> waves = massDependenceConfig_.GetWavesForL(l_value);
                        if (waves.size() > 1) {
                            // Multiple waves: use coherent sum
                            value = GetCoherentPhaseForL(mass_bin, base_name, massDepPars, l_value);
                        } else if (waves.size() == 1) {
                            // Single wave: direct evaluation
                            value = GetSingleWavePhase(mass_bin, base_name, massDepPars, waves[0]);
                        }
                    } else {
                        // Use mass-independent parameters
                        const auto bin_it = massIndepPars.find(mass_bin_str);
                        if (bin_it != massIndepPars.end()) {
                            const auto par_it = bin_it->second.find(par_name);
                            if (par_it != bin_it->second.end() && !par_it->second.empty()) {
                                value = par_it->second[0];
                            }
                        }
                    }
                }
                
                pars.SetCurrentVal(par_name, value);
            }
            
            total_chi2 += EvaluateChi2ForMassBin(mass_bin, pars);
        }
        
        lastChi2_ = total_chi2;
        return total_chi2;
    }
    void MassDependentEquationSolver::PrintEquations(const TString opt, const double mass_bin_center) const {
        // std::cout << "MassDependentEquationSolver::PrintEquations() : mass_bin_center = " << mass_bin_center << std::endl;
        
        const auto it = eqns_.find(mass_bin_center);
        if (it != eqns_.end()) {
            for (const auto& eqn : it->second) {
                eqn.Print(opt);
            }
        } else {
            std::cout << "No equations found for mass bin center: " << mass_bin_center << std::endl;
        }
    }

    void MassDependentEquationSolver::SetEquationValues(MassDependentMoments& moms) {
        for (auto& [mass_bin, equations] : eqns_) {
            for (auto& eqn : equations) {
                eqn.FindL();
                // Get a copy instead of a reference to avoid binding to temporary
                MomentHelper moments = moms.GetMoments(mass_bin);
                eqn.SetEquationValue(moments.GetUnnormalizedVal(eqn.GetName()),
                                   moments.GetUnnormalizedError(eqn.GetName()));
                eqn.DoEval(pars_[mass_bin].CurrentVals());
            }
        }
    }

    void MassDependentEquationSolver::PrintParNameIndices() const {
        for (const auto& [name, index] : parNameToIndex_) {
            std::cout << "Parameter: " << name << ", Index: " << index << std::endl;
        }
    }

    void MassDependentEquationSolver::PrintParCurrentVals(double mass_bin_center) const {
        const auto it = pars_.find(mass_bin_center);
        if (it != pars_.end()) {
            std::cout << "Current parameter values for mass bin center " << mass_bin_center << ":" << std::endl;
            it->second.Print("v");
        } else {
            std::cout << "No parameters found for mass bin center: " << mass_bin_center << std::endl;
        }
    }

    void MassDependentEquationSolver::MakeResultTree(const TString& fileName) const {
        if (pars_.empty()) {
            std::cerr << "No parameters to save" << std::endl;
            return;
        }

        std::unique_ptr<TFile> file(TFile::Open(fileName, "RECREATE"));
        if (!file || file->IsZombie()) {
            std::cerr << "Error creating file: " << fileName << std::endl;
            return;
        }

        auto tree = std::make_unique<TTree>("result", "result");
        const int nMassBins = static_cast<int>(massBins_.size());
        const int nPars = pars_.begin()->second.Nvars();
        
        // Use vectors for better memory management
        std::vector<double> par_vals(nPars);
        double mass_bin_center = 0.0;
        
        // Collect all unique H moment names from all mass bins
        std::set<TString> allHMomentNames;
        for (const auto& [mass_bin, equations] : eqns_) {
            for (const auto& eqn : equations) {
                TString eqnName = eqn.GetName();
                if (eqnName.BeginsWith("H_")) {
                    allHMomentNames.insert(eqnName);
                }
            }
        }
        
        tree->Branch("mass_bin", &mass_bin_center);
        for (int i = 0; i < nPars; ++i) {
            const TString par_name = pars_.begin()->second.GetParName(i);
            tree->Branch(par_name, &par_vals[i]);
        }
        
        // Add branches for H moment values only (no errors)
        std::map<TString, double> hMomentBranchValues;
        for (const TString& momentName : allHMomentNames) {
            hMomentBranchValues[momentName] = 0.0;
            tree->Branch(momentName, &hMomentBranchValues[momentName]);
        }

        std::cout << "Filling the result tree with parameter values and H moment values..." << std::endl;

        for (const auto& [mass_bin, pars] : pars_) {
            mass_bin_center = mass_bin;
            
            // Fill parameter values
            for (int i = 0; i < nPars; ++i) {
                const TString par_name = pars.GetParName(i);
                par_vals[i] = pars.GetCurrentVal(par_name);
            }
            
            // Fill H moment values for this mass bin
            for (const TString& momentName : allHMomentNames) {
                hMomentBranchValues[momentName] = 0.0;
                
                // Find the equation for this moment in this mass bin
                const auto eqnIt = eqns_.find(mass_bin);
                if (eqnIt != eqns_.end()) {
                    for (const auto& eqn : eqnIt->second) {
                        if (eqn.GetName() == momentName) {
                            hMomentBranchValues[momentName] = eqn.GetOrigFormulaValue();
                            break;
                        }
                    }
                }
            }
            
            tree->Fill();
        }

        // Save chi2 information (without seed for this version)
        auto chi2_tree = std::make_unique<TTree>("chi2", "chi2");
        double chi2 = lastChi2_;
        chi2_tree->Branch("chi2", &chi2);
        chi2_tree->Fill();

        file->Write();
        std::cout << "Result tree saved to " << fileName << " with H moment values as individual branches" << std::endl;
    }

    void MassDependentEquationSolver::MakeResultTree(const TString& fileName, int seed) const {
        if (pars_.empty()) {
            std::cerr << "No parameters to save" << std::endl;
            return;
        }

        std::unique_ptr<TFile> file(TFile::Open(fileName, "RECREATE"));
        if (!file || file->IsZombie()) {
            std::cerr << "Error creating file: " << fileName << std::endl;
            return;
        }

        auto tree = std::make_unique<TTree>("result", "result");
        const int nMassBins = static_cast<int>(massBins_.size());
        const int nPars = pars_.begin()->second.Nvars();
        
        // Use vectors for better memory management
        std::vector<double> par_vals(nPars);
        double mass_bin_center = 0.0;
        // Note: seed is NOT stored in the main result tree, only in chi2 tree
        
        // Collect all unique H moment names from all mass bins
        std::set<TString> allHMomentNames;
        for (const auto& [mass_bin, equations] : eqns_) {
            for (const auto& eqn : equations) {
                TString eqnName = eqn.GetName();
                if (eqnName.BeginsWith("H_")) {
                    allHMomentNames.insert(eqnName);
                }
            }
        }
        
        tree->Branch("mass_bin", &mass_bin_center);
        for (int i = 0; i < nPars; ++i) {
            const TString par_name = pars_.begin()->second.GetParName(i);
            tree->Branch(par_name, &par_vals[i]);
        }
        
        // Add branches for H moment values only (no errors)
        std::map<TString, double> hMomentBranchValues;
        for (const TString& momentName : allHMomentNames) {
            hMomentBranchValues[momentName] = 0.0;
            tree->Branch(momentName, &hMomentBranchValues[momentName]);
        }

        std::cout << "Filling the result tree with parameter values and H moment values (seed: " << seed << ")..." << std::endl;

        for (const auto& [mass_bin, pars] : pars_) {
            mass_bin_center = mass_bin;
            
            // Fill parameter values
            for (int i = 0; i < nPars; ++i) {
                const TString par_name = pars.GetParName(i);
                par_vals[i] = pars.GetCurrentVal(par_name);
            }
            
            // Fill H moment values for this mass bin
            for (const TString& momentName : allHMomentNames) {
                hMomentBranchValues[momentName] = 0.0;
                
                // Find the equation for this moment in this mass bin
                const auto eqnIt = eqns_.find(mass_bin);
                if (eqnIt != eqns_.end()) {
                    for (const auto& eqn : eqnIt->second) {
                        if (eqn.GetName() == momentName) {
                            hMomentBranchValues[momentName] = eqn.EqnValue();
                            break;
                        }
                    }
                }
            }
            
            tree->Fill();
        }

        // Save chi2 and seed information (seed only stored in chi2 tree)
        auto chi2_tree = std::make_unique<TTree>("chi2", "chi2");
        double chi2 = lastChi2_;
        int seed_val = seed;
        chi2_tree->Branch("chi2", &chi2);
        chi2_tree->Branch("seed", &seed_val);
        chi2_tree->Fill();

        file->Write();
        std::cout << "Result tree saved to " << fileName << " with seed " << seed << " and H moment values as individual branches" << std::endl;
    }

    // ParameterManager implementation
    void MassDependentEquationSolver::ParameterManager::AddMassIndependentParameters(
        const std::vector<double>& massBins, 
        const std::vector<TString>& parNames, 
        int seed) {
        
        TRandom3 rng(seed);
        for (const double massBin : massBins) {
            for (const TString& parName : parNames) {
                if (parName.Contains("2")) continue;

                TString name = Form("MI_%1.6f_%s", massBin, parName.Data());
                double initialValue = name.Contains("_1_") ? 0.0 : rng.Uniform(-100, 100);
                
                parsList[name] = initialValue;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                totalNpars++;
            }
        }
    }

    void MassDependentEquationSolver::ParameterManager::AddMassDependentParameters(
        const std::vector<TString>& parNames, 
        int seed,
        const std::vector<int>& massDependentL) {  // NEW: Add l values to treat as mass dependent
        
        for (const TString& parName : parNames) {
            // Step 1: Parse the parameter name to extract l value
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;  // Need at least prefix_l_m format
            
            // Extract l value from parameter name (assuming format: prefix_l_m)
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int l_value = l_str.Atoi();
            
            // Step 2: Check if this l value should be mass dependent
            bool isMassDependent = std::find(massDependentL.begin(), massDependentL.end(), l_value) 
                                  != massDependentL.end();
            
            if (!isMassDependent) continue;  // Skip if not in mass dependent list
            
            // Step 3: Skip phase parameters for now (can be extended later)
            if (parName.Contains("phi_0") || parName.Contains("phi_1") || parName.Contains("phi_2")) {
                continue;
            }

            // Step 4: Add parameter with generalized naming
            TString name = Form("MD_%s_%s", parName.Data(), Config::BW_NAMES[0].Data());
            parsList[name] = TRandom3(seed).Uniform(5, 20);  // Random initial value
            nameToIndex[name] = totalNpars;
            parIndexNames.push_back(name);
            totalNpars++;
        }
    }

    // Keep the old interface for backward compatibility
    void MassDependentEquationSolver::ParameterManager::AddMassDependentParameters(
        const std::vector<TString>& parNames, 
        int seed) {
        
        // Default behavior: only l=2 is mass dependent
        std::vector<int> defaultMassDependentL = {2};
        AddMassDependentParameters(parNames, seed, defaultMassDependentL);
    }

    // void MassDependentEquationSolver::ParameterManager::AddMassDependentParameters(
    //     const std::vector<TString>& parNames, 
    //     int seed) {
        
    //     for (const TString& parName : parNames) {
    //         if (!(parName.Contains("2")) || parName.Contains("phi_0") || 
    //             parName.Contains("phi_1") || parName.Contains("phi_2")) continue;

    //         TString name = Form("MD_%s_%s", parName.Data(), Config::BW_NAMES[0].Data());
    //         parsList[name] = 10.0;
    //         nameToIndex[name] = totalNpars;
    //         parIndexNames.push_back(name);
    //         totalNpars++;
    //     }
    // }

    std::vector<double> MassDependentEquationSolver::ParameterManager::GetInitialValues() const {
        std::vector<double> values(totalNpars);
        for (int i = 0; i < totalNpars; ++i) {
            values[i] = parsList.at(parIndexNames[i]);
        }
        return values;
    }

    // Chi2Function implementation
    double MassDependentEquationSolver::Chi2Function::operator()(const double* mass_dep_pars) const {
        MapParameters(mass_dep_pars);
        return solver_.DoEval(massDepPars_, massIndepPars_);
    }

    void MassDependentEquationSolver::Chi2Function::MapParameters(const double* mass_dep_pars) const {
        massDepPars_.clear();
        massIndepPars_.clear();
        
        for (size_t i = 0; i < parIndexNames_.size(); ++i) {
            const TString& parName = parIndexNames_[i];
            
            if (parName.BeginsWith("MI_")) {
                MapMassIndependentParameter(parName, mass_dep_pars[i]);
            } else if (parName.BeginsWith("MD_")) {
                MapMassDependentParameter(parName, mass_dep_pars[i]);
            }
        }
    }

    void MassDependentEquationSolver::Chi2Function::MapMassIndependentParameter(
        const TString& parName, double value) const {
        
        int firstUnderscore = parName.Index("_", 3);
        int secondUnderscore = parName.Index("_", firstUnderscore + 1);
        int thirdUnderscore = parName.Index("_", secondUnderscore + 1);
        
        TString parKey = parName(firstUnderscore + 1, thirdUnderscore + 2 - firstUnderscore);
        TString massBinStr = parName(3, firstUnderscore - 7);
        
        massIndepPars_[massBinStr][parKey].push_back(value);
    }

    void MassDependentEquationSolver::Chi2Function::MapMassDependentParameter(
        const TString& parName, double value) const {
        
        int lastUnderscore = parName.Last('_');
        TString parKey = parName(3, lastUnderscore - 3);
        massDepPars_[parKey].push_back(value);
    }

    // Static utility method
    std::vector<TString> MassDependentEquationSolver::GetParNames() {
        std::vector<TString> parNames;
        for (int l = 0; l <= 2; ++l) {
            for (int m = -l; m <= l; ++m) {
                parNames.push_back(Form("a_%d_%d", l, m));
                parNames.push_back(Form("b_%d_%d", l, m));
                parNames.push_back(Form("aphi_%d_%d", l, m));
                parNames.push_back(Form("bphi_%d_%d", l, m));
            }
        }
        return parNames;
    }

    // Minimization methods
    void MassDependentEquationSolver::MinimizeChi2(int seed) {
        // Setup parameter manager
        ParameterManager paramManager;
        std::vector<TString> parNames = GetParNames();
        
        paramManager.AddMassIndependentParameters(massBins_, parNames, seed, massDependenceConfig_, hMomentsConfig_, *this);
        paramManager.AddMassDependentParameters(parNames, seed, massDependenceConfig_, hMomentsConfig_, *this);
        
        MinimizeChi2(paramManager, seed);
    }

    void MassDependentEquationSolver::MinimizeChi2(const ParameterManager& paramManager, int seed) {
        std::cout << "Total number of parameters: " << paramManager.totalNpars << std::endl;

        // Get initial parameter values
        std::vector<double> initialValues = paramManager.GetInitialValues();

        // Setup minimizer with Chi2Function wrapper
        auto minimizer = std::unique_ptr<ROOT::Math::Minimizer>(
            ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
        
        if (!minimizer) {
            std::cerr << "Error: Cannot create Minuit2 minimizer" << std::endl;
            return;
        }

        Chi2Function chi2Function(*this, paramManager.parIndexNames);
        ROOT::Math::Functor functor(chi2Function, paramManager.totalNpars);
        
        TStopwatch timer;
        timer.Start();

        std::cout << "Chi2 before minimization: " << chi2Function(initialValues.data()) << std::endl;

        minimizer->SetFunction(functor);
        
        // Set parameters
        for (int i = 0; i < paramManager.totalNpars; ++i) {
            const TString& parName = paramManager.parIndexNames[i];
            std::cout << "Setting parameter " << i << ": " << parName 
                 << " with initial value " << initialValues[i] << std::endl;
            minimizer->SetVariable(i, parName.Data(), initialValues[i], 0.1);
        }

        minimizer->SetPrintLevel(2);
        minimizer->SetMaxFunctionCalls(10000);
        minimizer->SetTolerance(Config::DEFAULT_CHI2_TOLERANCE);

        bool isValid = minimizer->Minimize();

        double minChi2 = minimizer->MinValue();
        std::cout << "Chi2 after minimization: " << minChi2 << std::endl;

        timer.Stop();
        std::cout << "Minimization time: " << timer.RealTime() << " seconds" << std::endl;

        minimizer->PrintResults();

        if (isValid) {
            MakeResultTree(Form("MassDependentResults_%d.root", seed), seed);
        }
        
        // Print results for first mass bin
        if (!massBins_.empty()) {
            PrintEquations("v", massBins_[0]);
            PrintParCurrentVals(massBins_[0]);
        }
    }

    // Helper method for single wave magnitude (no coherent sum needed)
    double MassDependentEquationSolver::GetSingleWaveMagnitude(
        double mass_bin, 
        const TString& par_name,
        const std::map<TString, std::vector<double>>& massDepPars,
        const TString& waveName) const {
        
        TString key = Form("%s_%s", par_name.Data(), waveName.Data());
        const auto it = massDepPars.find(key);
        
        if (it != massDepPars.end() && !it->second.empty()) {
            int funcIndex = GetFunctionIndexForWave(waveName);
            
            if (funcIndex >= 0 && funcIndex < static_cast<int>(massDepFuncs_.size())) {
                double coupling = it->second[0];
                double magnitude = massDepFuncs_[funcIndex].GetPWMagnitude(mass_bin, coupling);
                
                // std::cout << "Single wave " << waveName << " magnitude for " << par_name 
                //          << ": " << magnitude << " (coupling=" << coupling << ")" << std::endl;
                
                return magnitude;
            }
        }
        
        return 0.0;
    }

    // Helper method for single wave phase (no coherent sum needed)
    double MassDependentEquationSolver::GetSingleWavePhase(
        double mass_bin, 
        const TString& base_name,
        const std::map<TString, std::vector<double>>& massDepPars,
        const TString& waveName) const {
        
        TString key = Form("%s_%s", base_name.Data(), waveName.Data());
        const auto it = massDepPars.find(key);
        
        if (it != massDepPars.end() && !it->second.empty()) {
            int funcIndex = GetFunctionIndexForWave(waveName);
            
            if (funcIndex >= 0 && funcIndex < static_cast<int>(massDepFuncs_.size())) {
                double coupling = it->second[0];
                double phase = massDepFuncs_[funcIndex].GetPWPhase(mass_bin, coupling);
                
                // std::cout << "Single wave " << waveName << " phase for " << base_name 
                //          << ": " << phase << " (coupling=" << coupling << ")" << std::endl;
                
                return phase;
            }
        }
        
        return 0.0;
    }

    // Helper method to get coherent magnitude for specific l value
    double MassDependentEquationSolver::GetCoherentMagnitudeForL(
        double mass_bin, 
        const TString& par_name,
        const std::map<TString, std::vector<double>>& massDepPars,
        int l_value) const {
        
        std::complex<double> total_amplitude(0.0, 0.0);
        
        // Get wave names for this l value from configuration
        std::vector<TString> waveNames = massDependenceConfig_.GetWavesForL(l_value);
        
        std::cout << "Computing coherent magnitude for L=" << l_value << " with " 
                 << waveNames.size() << " waves" << std::endl;
        
        for (const TString& waveName : waveNames) {
            TString key = Form("%s_%s", par_name.Data(), waveName.Data());
            const auto it = massDepPars.find(key);
            
            if (it != massDepPars.end() && !it->second.empty()) {
                int funcIndex = GetFunctionIndexForWave(waveName);
                
                if (funcIndex >= 0 && funcIndex < static_cast<int>(massDepFuncs_.size())) {
                    double coupling = it->second[0];
                    double magnitude = massDepFuncs_[funcIndex].GetPWMagnitude(mass_bin, coupling);
                    double phase = massDepFuncs_[funcIndex].GetPWPhase(mass_bin, coupling);
                    
                    // Add complex amplitude
                    std::complex<double> wave_amplitude = magnitude * std::complex<double>(std::cos(phase), std::sin(phase));
                    total_amplitude += wave_amplitude;
                    
                    std::cout << "  Wave " << waveName << " contribution: |A|=" << magnitude 
                             << ", phase=" << phase << ", coupling=" << coupling << std::endl;
                }
            }
        }
        
        double total_magnitude = std::abs(total_amplitude);
        std::cout << "Total coherent magnitude for L=" << l_value << " " << par_name << ": " << total_magnitude << std::endl;
        
        return total_magnitude;
    }

    // Helper method to get coherent phase for specific l value
    double MassDependentEquationSolver::GetCoherentPhaseForL(
        double mass_bin, 
        const TString& base_name,
        const std::map<TString, std::vector<double>>& massDepPars,
        int l_value) const {
        
        std::complex<double> total_amplitude(0.0, 0.0);
        
        // Get wave names for this l value from configuration
        std::vector<TString> waveNames = massDependenceConfig_.GetWavesForL(l_value);
        
        std::cout << "Computing coherent phase for L=" << l_value << " with " 
                 << waveNames.size() << " waves" << std::endl;
        
        for (const TString& waveName : waveNames) {
            TString key = Form("%s_%s", base_name.Data(), waveName.Data());
            const auto it = massDepPars.find(key);
            
            if (it != massDepPars.end() && !it->second.empty()) {
                int funcIndex = GetFunctionIndexForWave(waveName);
                
                if (funcIndex >= 0 && funcIndex < static_cast<int>(massDepFuncs_.size())) {
                    double coupling = it->second[0];
                    double magnitude = massDepFuncs_[funcIndex].GetPWMagnitude(mass_bin, coupling);
                    double phase = massDepFuncs_[funcIndex].GetPWPhase(mass_bin, coupling);
                    
                    // Add complex amplitude
                    std::complex<double> wave_amplitude = magnitude * std::complex<double>(std::cos(phase), std::sin(phase));
                    total_amplitude += wave_amplitude;
                }
            }
        }
        
        double total_phase = std::arg(total_amplitude);
        std::cout << "Total coherent phase for L=" << l_value << " " << base_name << ": " << total_phase << std::endl;
        
        return total_phase;
    }

    // Helper method to map wave names to function indices
    int MassDependentEquationSolver::GetFunctionIndexForWave(const TString& waveName) const {
        const auto it = waveToFunctionIndex_.find(waveName);
        if (it != waveToFunctionIndex_.end()) {
            return it->second;
        }
        
        std::cerr << "Warning: Unknown wave name " << waveName << std::endl;
        return -1;
    }

    // Update ParameterManager methods
    void MassDependentEquationSolver::ParameterManager::AddMassIndependentParameters(
        const std::vector<double>& massBins, 
        const std::vector<TString>& parNames, 
        int seed,
        const MassDependenceConfig& config,
        const HMomentsConfig& hConfig,
        const MassDependentEquationSolver& solver) {
        
        TRandom3 rng(seed);
        for (const double massBin : massBins) {
            for (const TString& parName : parNames) {
                // Parse parameter name to extract l value
                std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
                if (parts->GetEntries() < 3) continue;
                
                TString l_str = ((TObjString*)parts->At(1))->GetString();
                int l_value = l_str.Atoi();
                
                // Check if this parameter is needed for the H moments configuration  
                if (!solver.ParameterNeededForHMoments(l_value, hConfig)) {
                    std::cout << "Skipping parameter " << parName << " (L=" << l_value << " not needed for H moments)" << std::endl;
                    continue;
                }
                
                // Only add if this l value should be mass independent
                if (!config.IsMassIndependent(l_value)) continue;

                TString name = Form("MI_%1.6f_%s", massBin, parName.Data());
                double initialValue = name.Contains("_1_") ? 0.0 : rng.Uniform(-100, 100);
                
                parsList[name] = initialValue;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                totalNpars++;
                
                std::cout << "Added mass-independent parameter: " << name << " (L=" << l_value << ")" << std::endl;
            }
        }
    }

    void MassDependentEquationSolver::ParameterManager::AddMassDependentParameters(
        const std::vector<TString>& parNames, 
        int seed,
        const MassDependenceConfig& config,
        const HMomentsConfig& hConfig,
        const MassDependentEquationSolver& solver) {
        
        for (const TString& parName : parNames) {
            // Parse parameter name to extract l value
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int l_value = l_str.Atoi();
            
            // Check if this parameter is needed for the H moments configuration
            if (!solver.ParameterNeededForHMoments(l_value, hConfig)) {
                std::cout << "Skipping parameter " << parName << " (L=" << l_value << " not needed for H moments)" << std::endl;
                continue;
            }
            
            // Check if this l value should be mass dependent
            if (!config.IsMassDependent(l_value)) continue;
            
            // Skip phase parameters for now
            if (parName.Contains("phi_0") || parName.Contains("phi_1") || parName.Contains("phi_2")) {
                continue;
            }

            // Add parameters for each wave that contributes to this l value
            std::vector<TString> waveNames = config.GetWavesForL(l_value);
            for (const TString& waveName : waveNames) {
                TString name = Form("MD_%s_%s_%s", parName.Data(), waveName.Data(), Config::BW_NAMES[0].Data());
                parsList[name] = 10.0;  // Default initial value
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                totalNpars++;
                
                std::cout << "Added mass-dependent parameter for L=" << l_value 
                         << ": " << name << " (index: " << totalNpars-1 << ")" << std::endl;
            }
        }
    }

    double MassDependentEquationSolver::EvalChi2() const {
        double total_chi2 = 0.0;
        
        for (const auto& [mass_bin, pars] : pars_) {
            // Use current parameter values to evaluate chi2
            double chi2_for_bin = EvaluateChi2ForMassBin(mass_bin, const_cast<ParameterHelper&>(pars));
            total_chi2 += chi2_for_bin;
            
            std::cout << "Mass bin " << mass_bin << " chi2: " << chi2_for_bin << std::endl;
        }
        
        // Update the cached chi2 value
        const_cast<MassDependentEquationSolver*>(this)->lastChi2_ = total_chi2;
        
        std::cout << "Total chi2: " << total_chi2 << std::endl;
        return total_chi2;
    }

    void MassDependentEquationSolver::PrintIncludedMoments() const {
        std::cout << "\n=== H Moments Configuration ===" << std::endl;
        
        if (hMomentsConfig_.includeAll) {
            std::cout << "Including ALL H moments in chi2 calculation" << std::endl;
            return;
        }
        
        if (!hMomentsConfig_.includedL.empty()) {
            std::cout << "Including H moments for L values: ";
            for (size_t i = 0; i < hMomentsConfig_.includedL.size(); ++i) {
                std::cout << hMomentsConfig_.includedL[i];
                if (i < hMomentsConfig_.includedL.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
        
        if (!hMomentsConfig_.excludedL.empty()) {
            std::cout << "Excluding H moments for L values: ";
            for (size_t i = 0; i < hMomentsConfig_.excludedL.size(); ++i) {
                std::cout << hMomentsConfig_.excludedL[i];
                if (i < hMomentsConfig_.excludedL.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
        
        // Show which specific moments will be included
        std::cout << "Specific H moments included in chi2:" << std::endl;
        for (int L = 0; L <= 4; ++L) {
            bool shouldInclude = hMomentsConfig_.ShouldIncludeL(L);
            std::cout << "  L=" << L << ": " << (shouldInclude ? "YES" : "NO") << std::endl;
        }
    }

    int MassDependentEquationSolver::ExtractLFromEquationName(const TString& eqnName) const {
        // Equation names are typically like "H_alpha_L_M" (e.g., "H_0_0_0", "H_1_2_1", "H_2_4_2")
        // std::cout << "Extracting L from equation name: " << eqnName << std::endl;
        if (!eqnName.BeginsWith("H_")) {
            return -1; // Not an H moment equation
        }
        
        std::unique_ptr<TObjArray> parts(eqnName.Tokenize("_"));
        if (parts->GetEntries() >= 3) {
            // Format is H_alpha_L_M, so L is at index 2
            TString L_str = ((TObjString*)parts->At(2))->GetString();
            int L_value = L_str.Atoi();
            // std::cout << "  Parsed L=" << L_value << " from equation " << eqnName << std::endl;
            return L_value;
        }
        
        // std::cout << "  Could not parse L value from equation " << eqnName << std::endl;
        return -1; // Could not parse L value
    }

    bool MassDependentEquationSolver::ParameterNeededForHMoments(int l, const HMomentsConfig& hConfig) const {
        // Check what H moment L values are included
        std::set<int> includedHMomentL;
        
        // Determine which H moment L values are actually included
        for (int H_L = 0; H_L <= 6; ++H_L) {
            if (hConfig.ShouldIncludeL(H_L)) {
                includedHMomentL.insert(H_L);
            }
        }
        
        // Rule 1: If only L=4 H moments are included -> need only l=2 amplitude parameters
        if (includedHMomentL.size() == 1 && includedHMomentL.count(4) == 1) {
            return (l == 2);
        }
        
        // Rule 2: If L=2 and L=4 H moments are included -> need l=0 and l=2 amplitude parameters
        if (includedHMomentL.size() == 2 && 
            includedHMomentL.count(2) == 1 && includedHMomentL.count(4) == 1) {
            return (l == 0 || l == 2);
        }
        
        // Rule 3: If all L moments are included -> need all parameters l=0,1,2 (assuming lmax=2)
        if (hConfig.includeAll) {
            return (l <= 2);  // Include l=0,1,2 (assuming lmax=2)
        }
        
        // Default: if none of the specific cases match, include all parameters
        return (l <= 2);
    }

}
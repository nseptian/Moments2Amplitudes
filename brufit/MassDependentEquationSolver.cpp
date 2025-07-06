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
        : massBins_(std::move(mass_bins)) {
        
        InitializeMassBins(massBins_);
        InitializeFunctions(l_max);
        
        // Initialize equations for each mass bin
        const auto validEqn = [&noUse](const TString& name) {
            if (name.BeginsWith("a_")) return false; // Skip amplitude equations
            if (name == "normalise") return false; // Skip a_0_0
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

    void MassDependentEquationSolver::InitializeMassBins(const std::vector<double>& mass_bins) {
        if (mass_bins.empty()) {
            throw std::invalid_argument("Mass bins cannot be empty");
        }
        // Note: std::map doesn't have reserve(), it automatically manages memory
    }

    void MassDependentEquationSolver::InitializeFunctions(double l_max) {
        if (massBins_.size() < 2) {
            throw std::invalid_argument("At least two mass bins required");
        }
        
        const double bin_width = massBins_[1] - massBins_[0];
        massDepFuncs_.emplace_back(static_cast<int>(massBins_.size()), 
                                  massBins_.front(), bin_width, 0);
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
        
        for (auto& eqn : eqns_.at(massBin)) {
            eqn.SetNeedsRecalc();
            chi2 += eqn.DoEvalSq(current_vals);
        }
        
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
                
                double value = 0.0;
                if (!param_info->isPhase) {
                    // Handle magnitude parameters
                    if (param_info->l_str == "2") {
                        const auto it = massDepPars.find(par_name);
                        if (it != massDepPars.end() && !it->second.empty()) {
                            value = massDepFuncs_[0].GetPWMagnitude(mass_bin, it->second[0]);
                        }
                    } else {
                        // Mass-independent parameters
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
                    
                    if (param_info->l_str == "2") {
                        const auto it = massDepPars.find(base_name);
                        if (it != massDepPars.end() && !it->second.empty()) {
                            value = massDepFuncs_[0].GetPWPhase(mass_bin, it->second[0]);
                        }
                    } else {
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
        std::cout << "MassDependentEquationSolver::PrintEquations() : mass_bin_center = " << mass_bin_center << std::endl;
        
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
        
        tree->Branch("mass_bin", &mass_bin_center);
        for (int i = 0; i < nPars; ++i) {
            const TString par_name = pars_.begin()->second.GetParName(i);
            tree->Branch(par_name, &par_vals[i]);
        }

        std::cout << "Filling the result tree with parameter values..." << std::endl;

        for (const auto& [mass_bin, pars] : pars_) {
            mass_bin_center = mass_bin;
            for (int i = 0; i < nPars; ++i) {
                const TString par_name = pars.GetParName(i);
                par_vals[i] = pars.GetCurrentVal(par_name);
            }
            tree->Fill();
        }

        // Save chi2 information
        auto chi2_tree = std::make_unique<TTree>("chi2", "chi2");
        double chi2 = lastChi2_;
        chi2_tree->Branch("chi2", &chi2);
        chi2_tree->Fill();

        file->Write();
        std::cout << "Result tree saved to " << fileName << std::endl;
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
        int seed_val = seed;
        
        tree->Branch("mass_bin", &mass_bin_center);
        tree->Branch("seed", &seed_val);
        for (int i = 0; i < nPars; ++i) {
            const TString par_name = pars_.begin()->second.GetParName(i);
            tree->Branch(par_name, &par_vals[i]);
        }

        std::cout << "Filling the result tree with parameter values (seed: " << seed << ")..." << std::endl;

        for (const auto& [mass_bin, pars] : pars_) {
            mass_bin_center = mass_bin;
            for (int i = 0; i < nPars; ++i) {
                const TString par_name = pars.GetParName(i);
                par_vals[i] = pars.GetCurrentVal(par_name);
            }
            tree->Fill();
        }

        // Save chi2 and seed information
        auto chi2_tree = std::make_unique<TTree>("chi2", "chi2");
        double chi2 = lastChi2_;
        chi2_tree->Branch("chi2", &chi2);
        chi2_tree->Branch("seed", &seed_val);
        chi2_tree->Fill();

        file->Write();
        std::cout << "Result tree saved to " << fileName << " with seed " << seed << std::endl;
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
        int seed) {
        
        for (const TString& parName : parNames) {
            if (!(parName.Contains("2")) || parName.Contains("phi_0") || 
                parName.Contains("phi_1") || parName.Contains("phi_2")) continue;

            TString name = Form("MD_%s_%s", parName.Data(), Config::BW_NAMES[0].Data());
            parsList[name] = 10.0;
            nameToIndex[name] = totalNpars;
            parIndexNames.push_back(name);
            totalNpars++;
        }
    }

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
        
        paramManager.AddMassIndependentParameters(massBins_, parNames, seed);
        paramManager.AddMassDependentParameters(parNames, seed);
        
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
            minimizer->SetVariable(i, parName.Data(), initialValues[i], 1.0);
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

}
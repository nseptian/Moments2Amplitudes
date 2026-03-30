#include "MassDependentFitter.h"
#include "TObjString.h"
#include "TRegexp.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace m2pw {

    namespace {
        void LogMDInfo(const TString& stage, const TString& message) {
            std::cout << "[MD][" << stage << "] " << message << std::endl;
        }

        void LogMDWarn(const TString& stage, const TString& message) {
            std::cerr << "[MD][" << stage << "] " << message << std::endl;
        }

        void LogMIInfo(const TString& stage, const TString& message) {
            std::cout << "[MI][" << stage << "] " << message << std::endl;
        }

        void LogMIWarn(const TString& stage, const TString& message) {
            std::cerr << "[MI][" << stage << "] " << message << std::endl;
        }

        void LogFitInfo(const TString& stage, const TString& message) {
            std::cout << "[FIT][" << stage << "] " << message << std::endl;
        }

        void LogFitWarn(const TString& stage, const TString& message) {
            std::cerr << "[FIT][" << stage << "] " << message << std::endl;
        }

        TString GlobalPhaseKeyForWave(const TString& waveName) {
            return waveName;
        }

        TString SharedResonanceKeyForWave(const TString& waveName) {
            return Form("shared_%s", waveName.Data());
        }
    }

    MassDependentFitter::MassDependentFitter(const Setup& setup, 
                                                           std::vector<double> mass_bins, 
                                                           double l_max, double noise, 
                                                           std::vector<TString> noUse,
                                                           const MassDependenceConfig& massDepConfig,
                                                           const MomentsConfig& hMomentsConfig) 
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
    }

    void MassDependentFitter::InitializeMassBins(const std::vector<double>& mass_bins) {
        if (mass_bins.empty()) {
            throw std::invalid_argument("Mass bins cannot be empty");
        }
        // Note: std::map doesn't have reserve(), it automatically manages memory
    }

    void MassDependentFitter::InitializeMDFunctions() {
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
                allWaves.insert(MassDependenceConfig::StripPhaseTag(wave));
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
                LogFitInfo("MDFUNC", Form("index=%d wave=%s model=BreitWigner", funcIndex, wave.Data()));
                funcIndex++;
            } else if (wave == "a2_1700") {
                massDepFuncs_.emplace_back(static_cast<int>(massBins_.size()), 
                                          massBins_.front(), bin_width, 1);  // Breit-Wigner for a2(1700)
                waveToFunctionIndex_[wave] = funcIndex;
                LogFitInfo("MDFUNC", Form("index=%d wave=%s model=BreitWigner", funcIndex, wave.Data()));
                funcIndex++;
            } else if (wave == "a0_980") {
                massDepFuncs_.emplace_back(static_cast<int>(massBins_.size()), 
                                          massBins_.front(), bin_width, 2);  // Flatté for a0(980)
                waveToFunctionIndex_[wave] = funcIndex;
                LogFitInfo("MDFUNC", Form("index=%d wave=%s model=Flatte", funcIndex, wave.Data()));
                funcIndex++;
            } else if (wave == "pi1_1400") {
                massDepFuncs_.emplace_back(static_cast<int>(massBins_.size()), 
                                          massBins_.front(), bin_width, 3);  // Breit-Wigner for pi1(1400)
                waveToFunctionIndex_[wave] = funcIndex;
                LogFitInfo("MDFUNC", Form("index=%d wave=%s model=BreitWigner", funcIndex, wave.Data()));
                funcIndex++;
            } else if (wave.Contains("conformal") || wave.Contains("poly")) {
                // Extract polynomial order from wave name if specified (e.g., "conformal_S0_order3")
                int polyOrder = 2;  // Default order
                if (wave.Contains("order")) {
                    TString orderStr = wave;
                    orderStr.Remove(0, orderStr.Index("order") + 5);
                    orderStr = orderStr(TRegexp("[0-9]+"));
                    if (orderStr.IsDigit()) {
                        polyOrder = orderStr.Atoi();
                    }
                }
                massDepFuncs_.emplace_back(static_cast<int>(massBins_.size()), 
                                          massBins_.front(), bin_width, 4, polyOrder);  // Conformal polynomial
                waveToFunctionIndex_[wave] = funcIndex;
                LogFitInfo("MDFUNC", Form("index=%d wave=%s model=ConformalPolynomial order=%d", funcIndex, wave.Data(), polyOrder));
                funcIndex++;
            } else {
                LogFitWarn("MDFUNC", Form("wave=%s reason=unknown_wave skipped=1", wave.Data()));
            }
        }
        
        LogFitInfo("MDFUNC", Form("initialized_count=%zu", massDepFuncs_.size()));
    }

    void MassDependentFitter::SetMassDependenceConfig(const MassDependenceConfig& config) {
        massDependenceConfig_ = config;
    }

    void MassDependentFitter::SetMomentsConfig(const MomentsConfig& config) {
        hMomentsConfig_ = config;
    }

    unsigned int MassDependentFitter::NDim() const {
        if (pars_.empty()) return 0;
        
        unsigned int n = 0;
        for (const auto& [mass_bin, pars] : pars_) {
            n += pars.Nvars();
        }
        return n;
    }

    double MassDependentFitter::EvaluateChi2ForMassBin(double massBin, ParameterHelper& pars) const {
        double chi2 = 0.0;
        const auto current_vals = pars.CurrentVals();
        
        int includedEquations = 0;
        int excludedEquations = 0;
        
        for (auto& eqn : eqns_.at(massBin)) {
            // Extract L value from equation name to check if it should be included
            int L_value = ExtractLFromEquationName(eqn.GetName());
            
            // Apply H(L,M)s filtering
            bool shouldInclude = true;
            if (L_value >= 0) { // This is an H(L,M) equation
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

    std::unique_ptr<MassDependentFitter::ParameterInfo> 
    MassDependentFitter::ParseParameterName(const TString& parName) const {
        auto info = std::make_unique<ParameterInfo>(parName);
        
        std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
        if (parts->GetEntries() < 3) {
            return nullptr;
        }
        
        return info;
    }

    double MassDependentFitter::DoEval(
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
                            value = GetSingleWaveMagnitude(mass_bin, par_name, massDepPars, l_value, waves[0]);
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
                    // Handle phase parameters - check mass-independent first
                    TString base_name = param_info->name;
                    base_name.ReplaceAll("phi", "");
                    
                    // First check if we have mass-independent phase parameters for this parameter
                    bool hasPhaseInMassIndepPars = false;
                    const auto bin_it = massIndepPars.find(mass_bin_str);
                    if (bin_it != massIndepPars.end()) {
                        const auto par_it = bin_it->second.find(par_name);
                        if (par_it != bin_it->second.end() && !par_it->second.empty()) {
                            hasPhaseInMassIndepPars = true;
                            value = par_it->second[0];
                        }
                    }
                    
                    // If no mass-independent phase found, check for mass-dependent
                    if (!hasPhaseInMassIndepPars && massDependenceConfig_.IsMassDependent(l_value)) {
                        // Check if we need coherent sum or single wave
                        std::vector<TString> waves = massDependenceConfig_.GetWavesForL(l_value);
                        if (waves.size() > 1) {
                            // Multiple waves: use coherent sum
                            value = GetCoherentPhaseForL(mass_bin, base_name, massDepPars, l_value);
                        } else if (waves.size() == 1) {
                            // Single wave: direct evaluation
                            value = GetSingleWavePhase(mass_bin, base_name, massDepPars, l_value, waves[0]);
                        }
                    }
                    // If hasPhaseInMassIndepPars is true, value is already set above
                }
                
                pars.SetCurrentVal(par_name, value);
            }
            
            total_chi2 += EvaluateChi2ForMassBin(mass_bin, pars);
        }
        
        lastChi2_ = total_chi2;
        return total_chi2;
    }
    void MassDependentFitter::PrintEquations(const TString opt, const double mass_bin_center) const {
        // std::cout << "MassDependentFitter::PrintEquations() : mass_bin_center = " << mass_bin_center << std::endl;
        
        const auto it = eqns_.find(mass_bin_center);
        if (it != eqns_.end()) {
            for (const auto& eqn : it->second) {
                eqn.Print(opt);
            }
        } else {
            LogFitWarn("EQUATION", Form("mass_bin=%g reason=no_equations_found", mass_bin_center));
        }
    }

    void MassDependentFitter::SetEquationValues(MassDependentMoments& moms) {
        for (auto& [mass_bin, equations] : eqns_) {
            for (auto& eqn : equations) {
                eqn.FindL();
                // Get a copy instead of a reference to avoid binding to temporary
                MomentHelper moments = moms.GetMoments(mass_bin);
                // cout << "Setting equation value for " << eqn.GetName() 
                //      << " at mass bin " << mass_bin 
                //      << ": value = " << moments.GetUnnormalizedVal(eqn.GetName())
                //      << ", error = " << moments.GetUnnormalizedError(eqn.GetName()) << endl;
                eqn.SetEquationValue(moments.GetUnnormalizedVal(eqn.GetName()),
                                   moments.GetUnnormalizedError(eqn.GetName()));
                eqn.DoEval(pars_[mass_bin].CurrentVals());
            }
        }
    }

    void MassDependentFitter::PrintParNameIndices() const {
        for (const auto& [name, index] : parNameToIndex_) {
            LogFitInfo("PARAM", Form("name=%s index=%d", name.Data(), index));
        }
    }

    void MassDependentFitter::PrintParCurrentVals(double mass_bin_center) const {
        const auto it = pars_.find(mass_bin_center);
        if (it != pars_.end()) {
            LogFitInfo("PARAM", Form("dump_current_values mass_bin=%g", mass_bin_center));
            it->second.Print("v");
        } else {
            LogFitWarn("PARAM", Form("mass_bin=%g reason=no_parameters_found", mass_bin_center));
        }
    }

    void MassDependentFitter::MakeResultTree(const ParameterManager& paramManager, const TString& fileName) const {
        if (pars_.empty()) {
            LogFitWarn("RESULT", "reason=no_parameters_to_save");
            return;
        }

        std::unique_ptr<TFile> file(TFile::Open(fileName, "RECREATE"));
        if (!file || file->IsZombie()) {
            LogFitWarn("RESULT", Form("failed_to_create_file path=%s", fileName.Data()));
            return;
        }

        auto tree = std::make_unique<TTree>("result", "result");
        const int nMassBins = static_cast<int>(massBins_.size());
        const int nPars = pars_.begin()->second.Nvars();

        std::vector<double> par_vals(nPars);
        std::vector<double> par_errors(nPars);
        double mass_bin_center = 0.0;

        // Collect all unique H(L,M) names from all mass bins
        std::set<TString> allMomentNames;
        for (const auto& [mass_bin, equations] : eqns_) {
            for (const auto& eqn : equations) {
                TString eqnName = eqn.GetName();
                if (eqnName.BeginsWith("H_")) {
                    allMomentNames.insert(eqnName);
                }
            }
        }

        tree->Branch("mass_bin", &mass_bin_center);
        for (int i = 0; i < nPars; ++i) {
            const TString par_name = pars_.begin()->second.GetParName(i);
            tree->Branch(par_name, &par_vals[i]);
            tree->Branch(par_name + "_error", &par_errors[i]);
        }

        // Add branches for H(L,M) values only (no errors)
        std::map<TString, double> hMomentBranchValues;
        for (const TString& momentName : allMomentNames) {
            hMomentBranchValues[momentName] = 0.0;
            tree->Branch(momentName, &hMomentBranchValues[momentName]);
        }
        
        LogFitInfo("RESULT", "filling_tree=result content=amplitude_and_moments");

        for (const auto& [mass_bin, pars] : pars_) {
            mass_bin_center = mass_bin;

            // Fill parameter values
            for (int i = 0; i < nPars; ++i) {
                const TString par_name = pars.GetParName(i);
                par_vals[i] = pars.GetCurrentVal(par_name);
            }

            // Fill H(L,M) values for this mass bin
            for (const TString& momentName : allMomentNames) {
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

        // vector<double> ParameterValues = paramManager.GetValues();

        LogFitInfo("RESULT", "filling_tree=mass_dependent_params");

        // create new tree to store mass dependent parameters
        auto mdTree = std::make_unique<TTree>("mass_dependent_params", "Mass Dependent Parameters");
        
        const int nMDPars = paramManager.totalNpars;
        std::vector<double> md_par_vals(nMDPars);
        int idx = 0;
        for (const auto& [parName, value] : paramManager.parsList) {
            // double currentValue = value;
            if (parName.BeginsWith("MD_")) {  // Only mass-dependent parameters
                LogFitInfo("RESULT", Form("tree=mass_dependent_params param=%s value=%g", parName.Data(), value));
                md_par_vals[idx] = value;
                mdTree->Branch(parName, &md_par_vals[idx]);
                idx++;
            }
        }

        mdTree->Fill();

        
        auto chi2_tree = std::make_unique<TTree>("chi2", "chi2");
        double chi2 = lastChi2_;
        chi2_tree->Branch("chi2", &chi2);
        int seed_val = paramManager.randomSeed;
        bool isValid = minimizerIsValid_;
        int status = minimizerStatus_;

        LogFitInfo("MIN", Form("seed=%d is_valid=%d status=%d", seed_val, isValid, status));

        chi2_tree->Branch("isValid", &isValid);
        chi2_tree->Branch("Status", &status);
        chi2_tree->Fill();

        file->Write();
        LogFitInfo("RESULT", Form("saved path=%s", fileName.Data()));
    }

    std::vector<double> MassDependentFitter::ParameterManager::GetValues() const {
        std::vector<double> values(totalNpars);
        for (int i = 0; i < totalNpars; ++i) {
            values[i] = parsList.at(parIndexNames[i]);
        }
        return values;
    }

    bool MassDependentFitter::ParameterManager::SetParameterInitialValue(const TString& parName, double value) {
        auto it = parsList.find(parName);
        if (it == parsList.end()) {
            LogMDWarn("CONFIG", Form("set_initial_failed param=%s reason=not_found", parName.Data()));
            return false;
        }

        const double oldValue = it->second;
        it->second = value;
        LogMDInfo("CONFIG", Form("set_initial_ok param=%s old=%g new=%g", parName.Data(), oldValue, value));
        return true;
    }

    bool MassDependentFitter::ParameterManager::SetParameterLimits(const TString& parName, double lower, double upper) {
        if (!(lower < upper)) {
            LogMDWarn("CONFIG", Form("set_limits_failed param=%s reason=invalid_range lower=%g upper=%g", parName.Data(), lower, upper));
            return false;
        }

        auto it = parsList.find(parName);
        if (it == parsList.end()) {
            LogMDWarn("CONFIG", Form("set_limits_failed param=%s reason=not_found", parName.Data()));
            return false;
        }

        parameterLimits[parName] = {lower, upper};
        LogMDInfo("CONFIG", Form("set_limits_ok param=%s lower=%g upper=%g", parName.Data(), lower, upper));
        return true;
    }

    void MassDependentFitter::ParameterManager::SetParameterInitialValues(const std::map<TString, double>& values) {
        int successCount = 0;
        for (const auto& [parName, value] : values) {
            if (SetParameterInitialValue(parName, value)) {
                successCount++;
            }
        }
        LogMDInfo("CONFIG", Form("set_initial_batch total=%zu success=%d failed=%zu",
            values.size(), successCount, values.size() - successCount));
    }

    void MassDependentFitter::ParameterManager::SetParameterLimitsBatch(const std::map<TString, std::pair<double, double>>& limits) {
        int successCount = 0;
        for (const auto& [parName, range] : limits) {
            if (SetParameterLimits(parName, range.first, range.second)) {
                successCount++;
            }
        }
        LogMDInfo("CONFIG", Form("set_limits_batch total=%zu success=%d failed=%zu",
            limits.size(), successCount, limits.size() - successCount));
    }

    void MassDependentFitter::ParameterManager::StoreResults(const ROOT::Math::Minimizer* minimizer) {
        if (!minimizer) return;
        
        const double* errors = minimizer->Errors();
        parErrors.resize(totalNpars);
        
        // Store all parameter errors
        for (int i = 0; i < totalNpars; ++i) {
            parErrors[i] = errors[i];
        }
        
        // Propagate errors to the pars_ ParameterHelper objects for mass-independent parameters
        for (int i = 0; i < totalNpars; ++i) {
            const TString& parName = parIndexNames[i];
            
            if (parName.BeginsWith("MI_")) {
                // Extract mass bin and parameter info from mass-independent parameter name
                // Format: MI_1.420000_a_0_0
                int firstUnderscore = parName.Index("_", 3);  // Find first underscore after "MI_"
                if (firstUnderscore > 0) {
                    TString massBinStr = parName(3, firstUnderscore - 3);
                    double massBin = massBinStr.Atof();
                    
                    // Extract the actual parameter name (after the mass bin part)
                    TString actualParName = parName(firstUnderscore + 1, parName.Length() - firstUnderscore - 1);
                    
                    // Find the corresponding ParameterHelper and set the error
                    auto parsIt = fitter.pars_.find(massBin);
                    if (parsIt != fitter.pars_.end()) {
                        parsIt->second.SetCurrentError(actualParName, errors[i]);
                    }
                }
            }
            // Note: Mass-dependent parameters don't directly correspond to individual 
            // ParameterHelper entries, so we don't propagate their errors there
        }
    }

    // Chi2Function implementation
    double MassDependentFitter::Chi2Function::operator()(const double* mass_dep_pars) const {
        MapParameters(mass_dep_pars);
        return fitter_.DoEval(massDepPars_, massIndepPars_);
    }

    void MassDependentFitter::Chi2Function::MapParameters(const double* mass_dep_pars) const {
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

    void MassDependentFitter::Chi2Function::MapMassIndependentParameter(
        const TString& parName, double value) const {
        
        int firstUnderscore = parName.Index("_", 3);
        int secondUnderscore = parName.Index("_", firstUnderscore + 1);
        int thirdUnderscore = parName.Index("_", secondUnderscore + 1);
        
        TString parKey = parName(firstUnderscore + 1, thirdUnderscore + 2 - firstUnderscore);
        TString massBinStr = parName(3, firstUnderscore - 7);
        
        massIndepPars_[massBinStr][parKey].push_back(value);
    }

    void MassDependentFitter::Chi2Function::MapMassDependentParameter(
        const TString& parName, double value) const {
        
        int lastUnderscore = parName.Last('_');
        TString parKey = parName(3, lastUnderscore - 3);
        massDepPars_[parKey].push_back(value);
    }

    // Static utility method
    std::vector<TString> MassDependentFitter::GetParNames() {
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

    void MassDependentFitter::MinimizeChi2(ParameterManager& paramManager) {
        LogFitInfo("MIN", Form("total_parameters=%d", paramManager.totalNpars));

        // Get initial parameter values
        std::vector<double> initialValues = paramManager.GetValues();

        auto minimizer = std::unique_ptr<ROOT::Math::Minimizer>(
            ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));

        if (!minimizer) {
            LogFitWarn("MIN", "reason=minuit2_creation_failed");
            return;
        }

        Chi2Function chi2Function(*this, paramManager.parIndexNames);
        ROOT::Math::Functor functor(chi2Function, paramManager.totalNpars);

        TStopwatch timer;
        timer.Start();

        LogFitInfo("MIN", Form("chi2_before=%g", chi2Function(initialValues.data())));

        minimizer->SetFunction(functor);

        auto isFixed = [&paramManager](const TString& parName) {
            return std::find(paramManager.fixedParNames.begin(), paramManager.fixedParNames.end(), parName) != paramManager.fixedParNames.end();
        };

        for (int i = 0; i < paramManager.totalNpars; ++i) {
            const TString& parName = paramManager.parIndexNames[i];
            
            if (isFixed(parName)) {
                LogFitInfo("MIN", Form("set_param index=%d name=%s state=FIXED value=%g", i, parName.Data(), initialValues[i]));
                minimizer->SetFixedVariable(i, parName.Data(), initialValues[i]);
            }
            else {
                LogFitInfo("MIN", Form("set_param index=%d name=%s state=FREE init=%g", i, parName.Data(), initialValues[i]));
                
                double step = 0.1;
                if (parName.Contains("phi") || parName.Contains("globalphase")) {
                    step = 0.01;
                }
                minimizer->SetVariable(i, parName.Data(), initialValues[i], step);
                
                // Add phi parameter constraints
                if (parName.Contains("phi") || parName.Contains("globalphase")) {
                    minimizer->SetVariableLimits(i, -TMath::Pi(), TMath::Pi());
                }

                const auto limitsIt = paramManager.parameterLimits.find(parName);
                if (limitsIt != paramManager.parameterLimits.end() && limitsIt->second.first < limitsIt->second.second) {
                    minimizer->SetVariableLimits(i, limitsIt->second.first, limitsIt->second.second);
                }
            }
        }

        minimizer->SetPrintLevel(2);
        minimizer->SetMaxFunctionCalls(1000000);
        // minimizer->SetTolerance(1e-6);

        minimizerIsValid_ = minimizer->Minimize();
        minimizerStatus_ = minimizer->Status();
        lastChi2_ = minimizer->MinValue();

        LogFitInfo("MIN", Form("chi2_after=%g", lastChi2_));

        timer.Stop();
        LogFitInfo("MIN", Form("elapsed_seconds=%g", timer.RealTime()));

        minimizer->PrintResults();

        const double* parValues = minimizer->X();
        
        // update parameter manager values
        for (int i = 0; i < paramManager.totalNpars; ++i) {
            const TString& parName = paramManager.parIndexNames[i];
            paramManager.parsList[parName] = parValues[i];
        }
        
        // Store the errors in the parameter manager AND propagate to pars_
        paramManager.StoreResults(minimizer.get());
    }

    // Helper method for single wave magnitude (no coherent sum needed)
    double MassDependentFitter::GetSingleWaveMagnitude(
        double mass_bin, 
        const TString& par_name,
        const std::map<TString, std::vector<double>>& massDepPars,
        int l_value,
        const TString& waveName) const {
        
        TString key = Form("%s_%s", par_name.Data(), waveName.Data());
        const auto it = massDepPars.find(key);
        
        if (it != massDepPars.end() && !it->second.empty()) {
            int funcIndex = GetFunctionIndexForWave(waveName);
            
            if (funcIndex >= 0 && funcIndex < static_cast<int>(massDepFuncs_.size())) {
                std::vector<double> params = it->second;
                if (params.size() == 1) {
                    auto sharedIt = massDepPars.find(SharedResonanceKeyForWave(waveName));
                    if (sharedIt != massDepPars.end()) {
                        if (sharedIt->second.size() >= 1) params.push_back(sharedIt->second[0]);
                        if (sharedIt->second.size() >= 2) params.push_back(sharedIt->second[1]);
                    }
                }
                const bool includeGlobalPhase = massDependenceConfig_.UseGlobalPhaseForWave(l_value, waveName);
                if (includeGlobalPhase) {
                    const auto phaseIt = massDepPars.find(GlobalPhaseKeyForWave(waveName));
                    const double globalPhase = (phaseIt != massDepPars.end() && !phaseIt->second.empty()) ? phaseIt->second[0] : 0.0;
                    params.push_back(globalPhase);
                }

                double magnitude = massDepFuncs_[funcIndex].GetPWMagnitude(mass_bin, params, includeGlobalPhase);
                
                // std::cout << "Single wave " << waveName << " magnitude for " << par_name 
                //          << ": " << magnitude << std::endl;
                
                return magnitude;
            }
        }
        
        return 0.0;
    }

    // Helper method for single wave phase (no coherent sum needed)
    double MassDependentFitter::GetSingleWavePhase(
        double mass_bin, 
        const TString& base_name,
        const std::map<TString, std::vector<double>>& massDepPars,
        int l_value,
        const TString& waveName) const {
        
        TString key = Form("%s_%s", base_name.Data(), waveName.Data());
        const auto it = massDepPars.find(key);
        
        if (it != massDepPars.end() && !it->second.empty()) {
            int funcIndex = GetFunctionIndexForWave(waveName);
            
            if (funcIndex >= 0 && funcIndex < static_cast<int>(massDepFuncs_.size())) {
                std::vector<double> params = it->second;
                if (params.size() == 1) {
                    auto sharedIt = massDepPars.find(SharedResonanceKeyForWave(waveName));
                    if (sharedIt != massDepPars.end()) {
                        if (sharedIt->second.size() >= 1) params.push_back(sharedIt->second[0]);
                        if (sharedIt->second.size() >= 2) params.push_back(sharedIt->second[1]);
                    }
                }
                const bool includeGlobalPhase = massDependenceConfig_.UseGlobalPhaseForWave(l_value, waveName);
                if (includeGlobalPhase) {
                    const auto phaseIt = massDepPars.find(GlobalPhaseKeyForWave(waveName));
                    const double globalPhase = (phaseIt != massDepPars.end() && !phaseIt->second.empty()) ? phaseIt->second[0] : 0.0;
                    params.push_back(globalPhase);
                }

                double phase = massDepFuncs_[funcIndex].GetPWPhase(mass_bin, params, includeGlobalPhase);
                
                // std::cout << "Single wave " << waveName << " phase for " << base_name 
                //          << ": " << phase << std::endl;
                
                return phase;
            }
        }
        
        return 0.0;
    }

    // Helper method to get coherent magnitude for specific l value
    double MassDependentFitter::GetCoherentMagnitudeForL(
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
                    std::vector<double> params = it->second;
                    if (params.size() == 1) {
                        auto sharedIt = massDepPars.find(SharedResonanceKeyForWave(waveName));
                        if (sharedIt != massDepPars.end()) {
                            if (sharedIt->second.size() >= 1) params.push_back(sharedIt->second[0]);
                            if (sharedIt->second.size() >= 2) params.push_back(sharedIt->second[1]);
                        }
                    }
                    const bool includeGlobalPhase = massDependenceConfig_.UseGlobalPhaseForWave(l_value, waveName);
                    if (includeGlobalPhase) {
                        const auto phaseIt = massDepPars.find(GlobalPhaseKeyForWave(waveName));
                        const double globalPhase = (phaseIt != massDepPars.end() && !phaseIt->second.empty()) ? phaseIt->second[0] : 0.0;
                        params.push_back(globalPhase);
                    }

                    double magnitude = massDepFuncs_[funcIndex].GetPWMagnitude(mass_bin, params, includeGlobalPhase);
                    double phase = massDepFuncs_[funcIndex].GetPWPhase(mass_bin, params, includeGlobalPhase);
                    
                    // Add complex amplitude
                    std::complex<double> wave_amplitude = magnitude * std::complex<double>(std::cos(phase), std::sin(phase));
                    total_amplitude += wave_amplitude;
                    
                    std::cout << "  Wave " << waveName << " contribution: |A|=" << magnitude 
                             << ", phase=" << phase << std::endl;
                }
            }
        }
        
        double total_magnitude = std::abs(total_amplitude);
        std::cout << "Total coherent magnitude for L=" << l_value << " " << par_name << ": " << total_magnitude << std::endl;
        
        return total_magnitude;
    }

    // Helper method to get coherent phase for specific l value
    double MassDependentFitter::GetCoherentPhaseForL(
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
                    std::vector<double> params = it->second;
                    if (params.size() == 1) {
                        auto sharedIt = massDepPars.find(SharedResonanceKeyForWave(waveName));
                        if (sharedIt != massDepPars.end()) {
                            if (sharedIt->second.size() >= 1) params.push_back(sharedIt->second[0]);
                            if (sharedIt->second.size() >= 2) params.push_back(sharedIt->second[1]);
                        }
                    }
                    const bool includeGlobalPhase = massDependenceConfig_.UseGlobalPhaseForWave(l_value, waveName);
                    if (includeGlobalPhase) {
                        const auto phaseIt = massDepPars.find(GlobalPhaseKeyForWave(waveName));
                        const double globalPhase = (phaseIt != massDepPars.end() && !phaseIt->second.empty()) ? phaseIt->second[0] : 0.0;
                        params.push_back(globalPhase);
                    }

                    double magnitude = massDepFuncs_[funcIndex].GetPWMagnitude(mass_bin, params, includeGlobalPhase);
                    double phase = massDepFuncs_[funcIndex].GetPWPhase(mass_bin, params, includeGlobalPhase);
                    
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
    int MassDependentFitter::GetFunctionIndexForWave(const TString& waveName) const {
        const auto it = waveToFunctionIndex_.find(waveName);
        if (it != waveToFunctionIndex_.end()) {
            return it->second;
        }
        
        LogMDWarn("LOOKUP", Form("wave=%s reason=unknown_wave_name", waveName.Data()));
        return -1;
    }

    void MassDependentFitter::ParameterManager::AddMassIndependentParameters(int seed) {
        randomSeed = seed;
        
        auto parNames = MassDependentFitter::GetParNames();
        auto& config = fitter.GetMassDependenceConfig();
        auto& hConfig = fitter.GetMomentsConfig();
        auto& massBins = fitter.GetMassBins();
        
        TRandom3 rng(seed);
        for (const double massBin : massBins) {
            for (const TString& parName : parNames) {
                // Parse parameter name to extract l value
                std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
                if (parts->GetEntries() < 3) continue;
                
                TString l_str = ((TObjString*)parts->At(1))->GetString();
                int l_value = l_str.Atoi();
                
                // Check if this parameter is needed for the H(L,M)s configuration  
                if (!fitter.ParameterNeededForMoments(l_value, hConfig)) {
                    LogMIInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), l_value));
                    continue;
                }
                
                // Only add if this l value should be mass independent
                if (!config.IsMassIndependent(l_value)) continue;

                TString name = Form("MI_%1.6f_%s", massBin, parName.Data());
                
                // Check if parameter already exists - skip if already added
                if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                    LogMIInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                    continue;
                }
                
                double initialValue = name.Contains("phi") ? 0.0 : rng.Uniform(-200, 200);

                parsList[name] = initialValue;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                totalNpars++;

                LogMIInfo("ADD", Form("param=%s l=%d index=%d value=%g",
                    name.Data(), l_value, totalNpars - 1, initialValue));
            }
        }
    }

    void MassDependentFitter::ParameterManager::AddFixedMassIndependentParametersForL(
        const std::vector<int>& fixedL,
        const std::map<int, double>& fixedValues,
        bool magnitudeOnly) {
        
        auto parNames = MassDependentFitter::GetParNames();
        auto& hConfig = fitter.GetMomentsConfig();
        auto& massBins = fitter.GetMassBins();
        
        for (const double massBin : massBins) {
            for (const TString& parName : parNames) {
                std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
                if (parts->GetEntries() < 3) continue;
                
                TString l_str = ((TObjString*)parts->At(1))->GetString();
                int l_value = l_str.Atoi();
                
                // Check if this L should be fixed
                if (std::find(fixedL.begin(), fixedL.end(), l_value) == fixedL.end()) {
                    continue;  // This L is not in the fixed list
                }
                
                // Check if this parameter is needed for the H(L,M)s configuration  
                if (!fitter.ParameterNeededForMoments(l_value, hConfig)) {
                    LogMIInfo("SKIP", Form("param=%s reason=fixed_not_needed_for_selected_moments l=%d", parName.Data(), l_value));
                    continue;
                }
                
                bool isPhase = parName.Contains("phi");
                if (magnitudeOnly && isPhase) continue;  // Skip phases if only fixing magnitudes
                
                TString name = Form("MI_%1.6f_%s", massBin, parName.Data());
                auto fixedIt = fixedValues.find(l_value);
                double fixedValue = (fixedIt != fixedValues.end()) ? fixedIt->second : 0.0;
                
                parsList[name] = fixedValue;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                fixedParNames.push_back(name);  // Mark as fixed
                
                LogMIInfo("ADD", Form("param=%s l=%d state=FIXED index=%d value=%g",
                    name.Data(), l_value, totalNpars, fixedValue));
                totalNpars++;
            }
        }
    }

    void MassDependentFitter::ParameterManager::AddFixedMassIndependentParametersForL(
        const std::vector<std::string>& lReflectivities,
        const std::map<std::string, double>& fixedValues,
        bool magnitudeOnly) {
        
        auto parNames = MassDependentFitter::GetParNames();
        auto& hConfig = fitter.GetMomentsConfig();
        auto& massBins = fitter.GetMassBins();
        
        for (const std::string& lReflectivity : lReflectivities) {
            // Parse L and reflectivity from string like "1+" or "2-"
            if (lReflectivity.empty()) {
                LogMIWarn("CONFIG", "invalid_l_reflectivity reason=empty_string");
                continue;
            }
            
            char reflectivity = lReflectivity.back();  // Get last character (+ or -)
            std::string l_str = lReflectivity.substr(0, lReflectivity.length() - 1);  // Get L value part
            
            if (reflectivity != '+' && reflectivity != '-') {
                LogMIWarn("CONFIG", Form("invalid_reflectivity value=%c", reflectivity));
                continue;
            }
            
            int l_value;
            try {
                l_value = std::stoi(l_str);
            } catch (const std::exception& e) {
                LogMIWarn("CONFIG", Form("invalid_L_value value=%s token=%s", l_str.c_str(), lReflectivity.c_str()));
                continue;
            }
            
            // Check if this L is needed for the H(L,M)s configuration  
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) {
                LogMIInfo("SKIP", Form("l=%d reason=fixed_not_needed_for_selected_moments", l_value));
                continue;
            }
            
            // Determine prefix based on reflectivity (a for +, b for -)
            std::string prefix = (reflectivity == '+') ? "a" : "b";
            
            LogMIInfo("INIT", Form("fixed_params l=%d reflectivity=%c prefix=%s", l_value, reflectivity, prefix.c_str()));
            
            for (const double massBin : massBins) {
                for (const TString& parName : parNames) {
                    // Check if this parameter matches our L and reflectivity
                    std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
                    if (parts->GetEntries() < 3) continue;
                    
                    TString par_prefix = ((TObjString*)parts->At(0))->GetString();
                    TString par_l_str = ((TObjString*)parts->At(1))->GetString();
                    int par_l_value = par_l_str.Atoi();
                    
                    // Check if this parameter matches our L value and reflectivity
                    if (par_l_value != l_value) continue;
                    
                    bool isCorrectReflectivity = false;
                    if (reflectivity == '+' && (par_prefix == "a" || par_prefix == "aphi")) {
                        isCorrectReflectivity = true;
                    } else if (reflectivity == '-' && (par_prefix == "b" || par_prefix == "bphi")) {
                        isCorrectReflectivity = true;
                    }
                    
                    if (!isCorrectReflectivity) continue;
                    
                    bool isPhase = parName.Contains("phi");
                    if (magnitudeOnly && isPhase) continue;  // Skip phases if only fixing magnitudes
                    
                    TString name = Form("MI_%1.6f_%s", massBin, parName.Data());
                    auto fixedIt = fixedValues.find(parName.Data());
                    double fixedValue = (fixedIt != fixedValues.end()) ? fixedIt->second : 0.0;
                    
                    parsList[name] = fixedValue;
                    nameToIndex[name] = totalNpars;
                    parIndexNames.push_back(name);
                    fixedParNames.push_back(name);  // Mark as fixed
                    
                    LogMIInfo("ADD", Form("param=%s l=%d reflectivity=%c state=FIXED index=%d value=%g",
                        name.Data(), l_value, reflectivity, totalNpars, fixedValue));
                    totalNpars++;
                }
            }
        }
    }

    void MassDependentFitter::ParameterManager::AddMassIndependentParametersForL(
        const std::vector<int>& targetL,
        const TString& resultTreeFile) {
        
        auto parNames = MassDependentFitter::GetParNames();
        auto& massBins = fitter.GetMassBins();
        
        // Load mass-independent parameters from result tree file
        TFile file(resultTreeFile, "READ");
        if (!file.IsOpen() || file.IsZombie()) {
            LogMIWarn("FILE", Form("failed_to_open path=%s", resultTreeFile.Data()));
            return;
        }
        
        TTree* tree = dynamic_cast<TTree*>(file.Get("result"));
        if (!tree) {
            LogMIWarn("FILE", Form("missing_tree tree=result path=%s", resultTreeFile.Data()));
            return;
        }

        LogMIInfo("FILE", Form("loading_from_file path=%s", resultTreeFile.Data()));

        // Convert targetL to set for faster lookup
        std::set<int> targetLSet(targetL.begin(), targetL.end());

        // Set up branches to read mass_bin and parameter values
        double mass_bin_from_file = 0.0;
        tree->SetBranchAddress("mass_bin", &mass_bin_from_file);
        
        // Create maps to store parameter values for each parameter name
        std::map<TString, double> paramValues;
        for (const TString& parName : parNames) {
            // Parse parameter name to extract l value
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int l_value = l_str.Atoi();
            
            // Check if this l value is in our target list
            if (targetLSet.find(l_value) == targetLSet.end()) {
                continue;
            }
            
            // Set up branch address for this parameter
            TBranch* branch = tree->GetBranch(parName);
            if (!branch) {
                LogMIWarn("FILE", Form("missing_branch branch=%s", parName.Data()));
                continue;
            }
            
            paramValues[parName] = 0.0;
            tree->SetBranchAddress(parName, &paramValues[parName]);
        }

        // Read each entry (one per mass bin) and extract parameter values
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t entry = 0; entry < nEntries; entry++) {
            tree->GetEntry(entry);
            
            // Check if this mass bin is in our target list
            bool foundMassBin = false;
            for (double targetMassBin : massBins) {
                if (TMath::Abs(mass_bin_from_file - targetMassBin) < 1e-6) {
                    foundMassBin = true;
                    break;
                }
            }
            
            if (!foundMassBin) continue;
            
            // Add parameters for this mass bin
            for (const auto& [parName, value] : paramValues) {
                // Create mass-independent parameter name with mass bin
                TString name = Form("MI_%1.6f_%s", mass_bin_from_file, parName.Data());
                
                // Check if parameter already exists - skip if already added
                if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                    LogMIInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                    continue;
                }
                
                parsList[name] = value;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                
                // Parse L value for logging
                std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
                TString l_str = ((TObjString*)parts->At(1))->GetString();
                int l_value = l_str.Atoi();
                
                LogMIInfo("ADD", Form("param=%s l=%d mass_bin=%g index=%d value=%g",
                    name.Data(), l_value, mass_bin_from_file, totalNpars, value));
                totalNpars++;
            }
        }
        
        file.Close();
    }

    void MassDependentFitter::ParameterManager::AddMassIndependentParametersForL(
        const std::vector<int>& targetL,
        const std::vector<double>& massBins,
        const TString& resultTreeFile) {
        
        auto parNames = MassDependentFitter::GetParNames();
        
        // Load mass-independent parameters from result tree file
        TFile file(resultTreeFile, "READ");
        if (!file.IsOpen() || file.IsZombie()) {
            LogMIWarn("FILE", Form("failed_to_open path=%s", resultTreeFile.Data()));
            return;
        }
        
        TTree* tree = dynamic_cast<TTree*>(file.Get("result"));
        if (!tree) {
            LogMIWarn("FILE", Form("missing_tree tree=result path=%s", resultTreeFile.Data()));
            return;
        }

        LogMIInfo("FILE", Form("loading_from_file path=%s selected_mass_bins=%zu", resultTreeFile.Data(), massBins.size()));

        // Convert targetL to set for faster lookup
        std::set<int> targetLSet(targetL.begin(), targetL.end());
        
        // Convert massBins to set for faster lookup
        std::set<double> targetMassBinSet;
        for (double mb : massBins) {
            targetMassBinSet.insert(mb);
        }

        // Set up branches to read mass_bin and parameter values
        double mass_bin_from_file = 0.0;
        tree->SetBranchAddress("mass_bin", &mass_bin_from_file);
        
        // Create maps to store parameter values for each parameter name
        std::map<TString, double> paramValues;
        for (const TString& parName : parNames) {
            // Parse parameter name to extract l value
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int l_value = l_str.Atoi();
            
            // Check if this l value is in our target list
            if (targetLSet.find(l_value) == targetLSet.end()) {
                continue;
            }
            
            // Set up branch address for this parameter
            TBranch* branch = tree->GetBranch(parName);
            if (!branch) {
                LogMIWarn("FILE", Form("missing_branch branch=%s", parName.Data()));
                continue;
            }
            
            paramValues[parName] = 0.0;
            tree->SetBranchAddress(parName, &paramValues[parName]);
        }

        // Read each entry (one per mass bin) and extract parameter values
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t entry = 0; entry < nEntries; entry++) {
            tree->GetEntry(entry);
            
            // Check if this mass bin is in our target list
            bool foundMassBin = false;
            for (double targetMassBin : massBins) {
                if (TMath::Abs(mass_bin_from_file - targetMassBin) < 1e-6) {
                    foundMassBin = true;
                    break;
                }
            }
            
            if (!foundMassBin) continue;
            
            // Add parameters for this mass bin
            for (const auto& [parName, value] : paramValues) {
                // Create mass-independent parameter name with mass bin
                TString name = Form("MI_%1.6f_%s", mass_bin_from_file, parName.Data());
                
                // Check if parameter already exists - skip if already added
                if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                    LogMIInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                    continue;
                }
                
                parsList[name] = value;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                
                // Parse L value for logging
                std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
                TString l_str = ((TObjString*)parts->At(1))->GetString();
                int l_value = l_str.Atoi();
                
                LogMIInfo("ADD", Form("param=%s l=%d mass_bin=%g index=%d value=%g",
                    name.Data(), l_value, mass_bin_from_file, totalNpars, value));
                totalNpars++;
            }
        }
        
        file.Close();
    }

    void MassDependentFitter::ParameterManager::AddMassDependentParameters(
        int seed) {
        
        auto parNames = MassDependentFitter::GetParNames();
        auto& config = fitter.GetMassDependenceConfig();
        auto& hConfig = fitter.GetMomentsConfig();
        
        randomSeed = seed;

        for (const auto& [l_value, waves] : config.massDependentWaves) {
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) continue;
            for (const TString& configuredWave : waves) {
                if (!MassDependenceConfig::HasPhaseTag(configuredWave)) continue;

                const TString waveName = MassDependenceConfig::StripPhaseTag(configuredWave);
                TString phaseParName = Form("MD_%s_globalphase", waveName.Data());
                if (parsList.find(phaseParName) != parsList.end() || nameToIndex.find(phaseParName) != nameToIndex.end()) {
                    continue;
                }

                parsList[phaseParName] = 0.0;
                nameToIndex[phaseParName] = totalNpars;
                parIndexNames.push_back(phaseParName);
                totalNpars++;

                std::cout << "Added mass-dependent parameter for l=" << l_value
                    << ": " << phaseParName
                    << " (global phase, index: " << totalNpars-1 << ") = 0" << std::endl;
            }
        }

        for (const TString& parName : parNames) {
            // Parse parameter name to extract l value
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int l_value = l_str.Atoi();
            
            // Check if this parameter is needed for the H(L,M)s configuration
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) {
                std::cout << "Skipping parameter " << parName << " l=" << l_value << std::endl;
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
                TString name = Form("MD_%s_%s_k", parName.Data(), waveName.Data());
                // parsList[name] = 10.0;  // Default initial value
                parsList[name] = TRandom3(seed).Uniform(-50, 50);  // Random initial value
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                totalNpars++;
                
                std::cout << "Added mass-dependent parameter for l=" << l_value 
                         << ": " << name << " (index: " << totalNpars-1 << ")" << std::endl;
            }
        }
    }


    void MassDependentFitter::ParameterManager::AddMassDependentParameters(const TString filePath,
                                                                                const bool isFixed) {
        auto parNames = MassDependentFitter::GetParNames();
        auto& config = fitter.GetMassDependenceConfig();
        auto& hConfig = fitter.GetMomentsConfig();
        
        // Load mass-dependent parameters from file
        TFile file(filePath, "READ");
        if (!file.IsOpen() || file.IsZombie()) {
            std::cerr << "Error opening file: " << filePath << std::endl;
            return;
        }
        TTree* tree = dynamic_cast<TTree*>(file.Get("mass_dependent_params"));
        if (!tree) {
            std::cerr << "Error: Tree 'mass_dependent_params' not found in file: " << filePath << std::endl;
            return;
        }

        tree->GetEntry(0);  // Load first entry to get parameter names
        cout << "Loading mass dependent parameter initial values from file: " << filePath << std::endl;

        for (const auto& [l_value, waves] : config.massDependentWaves) {
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) continue;
            for (const TString& configuredWave : waves) {
                if (!MassDependenceConfig::HasPhaseTag(configuredWave)) continue;

                const TString waveName = MassDependenceConfig::StripPhaseTag(configuredWave);
                TString phaseParName = Form("MD_%s_globalphase", waveName.Data());
                if (parsList.find(phaseParName) != parsList.end() || nameToIndex.find(phaseParName) != nameToIndex.end()) {
                    continue;
                }

                double value = 0.0;
                if (TLeaf* leaf = tree->GetLeaf(phaseParName)) {
                    value = leaf->GetValue();
                }

                parsList[phaseParName] = value;
                nameToIndex[phaseParName] = totalNpars;
                parIndexNames.push_back(phaseParName);
                if (isFixed) {
                    fixedParNames.push_back(phaseParName);
                }
                totalNpars++;

                std::cout << "Added mass-dependent parameter for l=" << l_value
                    << ": " << phaseParName
                    << " (global phase, " << (isFixed ? "FIXED" : "FREE")
                    << ", index: " << totalNpars-1 << ") = " << value << std::endl;
            }
        }

        for (const TString& parName : parNames) {
            // Parse parameter name to extract l value
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int l_value = l_str.Atoi();
            
            // Check if this parameter is needed for the H(L,M)s configuration
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) {
                std::cout << "Skipping parameter " << parName << " (L=" << l_value << " not needed for H(L,M)s)" << std::endl;
                continue;
            }
            
            // Check if this l value should be mass dependent
            if (!config.IsMassDependent(l_value)) continue;
            
            // Skip phase parameters for now
            if (parName.Contains("phi")) {
                continue;
            }

            // Add parameters for each wave that contributes to this l value
            std::vector<TString> waveNames = config.GetWavesForL(l_value);
            for (const TString& waveName : waveNames) {
                TString name = Form("MD_%s_%s_k", parName.Data(), waveName.Data());
                // parsList[name] = 10.0;  // Default initial value
                parsList[name] = tree->GetLeaf(name)->GetValue();
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                if (isFixed) {
                    fixedParNames.push_back(name);
                    cout << "Added fixed mass-dependent parameter for L=" << l_value 
                         << ": " << name << " (index: " << totalNpars << ") = " << parsList[name] << std::endl;
                }
                else {
                    std::cout << "Added free mass-dependent parameter for L=" << l_value 
                         << ": " << name << " (index: " << totalNpars << ") = " << parsList[name] << std::endl;
                }
                totalNpars++;
                
                
            }
        }
    }

    double MassDependentFitter::EvalChi2() const {
        double total_chi2 = 0.0;
        
        for (const auto& [mass_bin, pars] : pars_) {
            // Use current parameter values to evaluate chi2
            double chi2_for_bin = EvaluateChi2ForMassBin(mass_bin, const_cast<ParameterHelper&>(pars));
            total_chi2 += chi2_for_bin;
            
            LogFitInfo("CHI2", Form("mass_bin=%g chi2=%g", mass_bin, chi2_for_bin));
        }
        
        // Update the cached chi2 value
        const_cast<MassDependentFitter*>(this)->lastChi2_ = total_chi2;
        
        LogFitInfo("CHI2", Form("total=%g", total_chi2));
        return total_chi2;
    }

    void MassDependentFitter::PrintIncludedMoments() const {
        LogFitInfo("MOMENTS", "configuration_begin");
        
        if (hMomentsConfig_.includeAll) {
            LogFitInfo("MOMENTS", "include_all=1");
            return;
        }
        
        if (!hMomentsConfig_.includedL.empty()) {
            TString includedValues;
            for (size_t i = 0; i < hMomentsConfig_.includedL.size(); ++i) {
                includedValues += Form("%d", hMomentsConfig_.includedL[i]);
                if (i < hMomentsConfig_.includedL.size() - 1) includedValues += ",";
            }
            LogFitInfo("MOMENTS", Form("included_L=%s", includedValues.Data()));
        }
        
        if (!hMomentsConfig_.excludedL.empty()) {
            TString excludedValues;
            for (size_t i = 0; i < hMomentsConfig_.excludedL.size(); ++i) {
                excludedValues += Form("%d", hMomentsConfig_.excludedL[i]);
                if (i < hMomentsConfig_.excludedL.size() - 1) excludedValues += ",";
            }
            LogFitInfo("MOMENTS", Form("excluded_L=%s", excludedValues.Data()));
        }
        
        // Show which specific moments will be included
        for (int L = 0; L <= 4; ++L) {
            bool shouldInclude = hMomentsConfig_.ShouldIncludeL(L);
            LogFitInfo("MOMENTS", Form("L=%d include=%d", L, shouldInclude));
        }
    }

    int MassDependentFitter::ExtractLFromEquationName(const TString& eqnName) const {
        // Equation name format: "H_alpha_L_M"
        // std::cout << "Extracting L from equation name: " << eqnName << std::endl;
        if (!eqnName.BeginsWith("H_")) {
            return -1; // Not an H(L,M) equation
        }
        
        std::unique_ptr<TObjArray> parts(eqnName.Tokenize("_"));
        if (parts->GetEntries() >= 3) {
            // L is at index 2
            TString L_str = ((TObjString*)parts->At(2))->GetString();
            int L_value = L_str.Atoi();
            // std::cout << "  Parsed L=" << L_value << " from equation " << eqnName << std::endl;
            return L_value;
        }
        
        // std::cout << "  Could not parse L value from equation " << eqnName << std::endl;
        return -1; // Could not parse L value
    }

    bool MassDependentFitter::ParameterNeededForMoments(int l, const MomentsConfig& hConfig) const {
        // Check what H(L,M) L values are included
        std::set<int> includedMomentL;
        
        // Determine which H(L,M) L values are actually included
        for (int H_L = 0; H_L <= 6; ++H_L) {
            if (hConfig.ShouldIncludeL(H_L)) {
                includedMomentL.insert(H_L);
            }
        }
        
        // Rule 1: If only L=4 H(L,M)s are included -> need only l=2 amplitude parameters
        if (includedMomentL.size() == 1 && includedMomentL.count(4) == 1) {
            return (l == 2);
        }
        
        // Rule 2: If L=2 and L=4 H(L,M)s are included -> need l=0 and l=2 amplitude parameters
        if (includedMomentL.size() == 2 && 
            includedMomentL.count(2) == 1 && includedMomentL.count(4) == 1) {
            return (l == 0 || l == 2);
        }
        
        // Rule 3: If all L moments are included -> need all parameters l=0,1,2 (assuming lmax=2)
        if (hConfig.includeAll) {
            return (l <= 2);  // Include l=0,1,2 (assuming lmax=2)
        }
        
        // Default: if none of the specific cases match, include all parameters
        return (l <= 2);
    }

    void MassDependentFitter::ParameterManager::AddMassDependentParametersForL(
        const std::vector<int>& targetL,
        const int seed,
        const bool magnitudeOnly,
        const bool isFixed,
        const bool yieldOnly) {

        auto parNames = MassDependentFitter::GetParNames();
        auto& config = fitter.GetMassDependenceConfig();
        auto& hConfig = fitter.GetMomentsConfig();
        
        randomSeed = seed;
        TRandom3 rng(seed);


        for (const TString& parName : parNames) {
            // Parse parameter name to extract l value
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int l_value = l_str.Atoi();
            
            // Check if this l value is in the target list
            if (std::find(targetL.begin(), targetL.end(), l_value) == targetL.end()) {
                LogMDInfo("SKIP", Form("param=%s reason=targetL_mismatch l=%d", parName.Data(), l_value));
                continue;
            }
            
            // Check if this parameter is needed for the H(L,M)s configuration  
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) {
                LogMDInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), l_value));
                continue;
            }
            
            // Only add if this l value should be mass dependent according to config
            if (!config.IsMassDependent(l_value)) {
                LogMDInfo("SKIP", Form("param=%s reason=l_not_mass_dependent l=%d", parName.Data(), l_value));
                continue;
            }

            // Skip phase parameters if magnitudeOnly is true OR yieldOnly is true
            if ((magnitudeOnly || yieldOnly) && parName.Contains("phi")) {
                LogMDInfo("SKIP", Form("param=%s reason=phase_disabled magnitudeOnly=%d yieldOnly=%d",
                    parName.Data(), magnitudeOnly, yieldOnly));
                continue;
            }

            // Get waves for this l value
            std::vector<TString> waves = config.GetWavesForL(l_value);
            if (waves.empty()) {
                LogMDWarn("SKIP", Form("param=%s reason=no_waves_defined_for_l l=%d", parName.Data(), l_value));
                continue;
            }

            for (const TString& waveName : waves) {
                // Get function index for this wave
                int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                if (funcIndex < 0 || funcIndex >= static_cast<int>(fitter.massDepFuncs_.size())) {
                    LogMDWarn("INIT", Form("wave=%s reason=missing_function_index using_defaults=1", waveName.Data()));
                }

                // Check if this is a conformal polynomial (FuncType=4)
                bool isConformal = false;
                int polyOrder = 2;
                if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                    isConformal = (fitter.massDepFuncs_[funcIndex].GetFuncType() == 4);
                    if (isConformal) {
                        polyOrder = fitter.massDepFuncs_[funcIndex].GetPolyOrder();
                    }
                }

                // Determine which parameter types to add
                std::vector<TString> paramTypes;
                if (isConformal) {
                    // Conformal polynomial: add re_i and im_i for i=0 to polyOrder
                    // No k parameter needed - re_0 and im_0 provide the overall scale
                    for (int i = 0; i <= polyOrder; ++i) {
                        paramTypes.push_back(Form("re_%d", i));
                        paramTypes.push_back(Form("im_%d", i));
                    }
                    LogMDInfo("MODE", Form("wave=%s model=conformal order=%d n_params=%zu",
                        waveName.Data(), polyOrder, paramTypes.size()));
                } else if (yieldOnly) {
                    paramTypes = {"k"};  // Only coupling parameters
                    LogMDInfo("MODE", Form("wave=%s yieldOnly=1 params=k", waveName.Data()));
                } else {
                    paramTypes = {"k", "M", "width"};  // All parameters
                }

                for (const TString& paramType : paramTypes) {
                    TString name;
                    if ((paramType == "M" || paramType == "width") && !isConformal) {
                        name = Form("MD_%s_%s", SharedResonanceKeyForWave(waveName).Data(), paramType.Data());
                    } else {
                        name = Form("MD_%s_%s_%s", parName.Data(), waveName.Data(), paramType.Data());
                    }
                    
                    // Check if parameter already exists - skip if already added
                    if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                        LogMDInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                        continue;
                    }

                    double initialValue;
                    if (isConformal) {
                        // Conformal polynomial parameters
                        if (paramType.BeginsWith("re_") || paramType.BeginsWith("im_")) {
                            // Extract coefficient index
                            TString indexStr = paramType;
                            indexStr.Remove(0, 3);
                            int coeffIndex = indexStr.Atoi();
                            
                            // Scale down higher order terms
                            double scale = TMath::Power(0.5, coeffIndex);
                            initialValue = rng.Uniform(-50, 50) * scale;
                        } else {
                            initialValue = 0.0;
                        }
                    } else if (paramType == "k") {
                        initialValue = rng.Uniform(-30, 30);  // Random for coupling strength
                    } else if (paramType == "M") {
                        // Get mass from MassDependentFunction
                        if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                            initialValue = fitter.massDepFuncs_[funcIndex].GetResonanceMass();
                            LogMDInfo("INIT", Form("param=%s source=function_%d field=mass value=%g",
                                name.Data(), funcIndex, initialValue));
                        } else {
                            initialValue = 1.5;  // Fallback default
                            LogMDInfo("INIT", Form("param=%s source=default field=mass value=%g", name.Data(), initialValue));
                        }
                    } else {  // width
                        // Get width from MassDependentFunction
                        if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                            initialValue = fitter.massDepFuncs_[funcIndex].GetResonanceWidth();
                            LogMDInfo("INIT", Form("param=%s source=function_%d field=width value=%g",
                                name.Data(), funcIndex, initialValue));
                        } else {
                            initialValue = 0.1;  // Fallback default
                            LogMDInfo("INIT", Form("param=%s source=default field=width value=%g", name.Data(), initialValue));
                        }
                    }

                    parsList[name] = initialValue;
                    nameToIndex[name] = totalNpars;
                    parIndexNames.push_back(name);
                    
                    if (isFixed) {
                        fixedParNames.push_back(name);
                        LogMDInfo("ADD", Form("scope=L%d param=%s state=FIXED index=%d value=%g",
                            l_value, name.Data(), totalNpars, initialValue));
                    } else {
                        LogMDInfo("ADD", Form("scope=L%d param=%s state=FREE index=%d value=%g",
                            l_value, name.Data(), totalNpars, initialValue));
                    }
                    
                    totalNpars++;
                }
            }
        }

        for (const int l_value : targetL) {
            if (!config.IsMassDependent(l_value)) continue;
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) continue;
            const auto wavesIt = config.massDependentWaves.find(l_value);
            if (wavesIt == config.massDependentWaves.end()) continue;
            for (const TString& configuredWave : wavesIt->second) {
                if (!MassDependenceConfig::HasPhaseTag(configuredWave)) continue;

                const TString waveName = MassDependenceConfig::StripPhaseTag(configuredWave);
                TString phaseParName = Form("MD_%s_globalphase", waveName.Data());
                if (parsList.find(phaseParName) != parsList.end() || nameToIndex.find(phaseParName) != nameToIndex.end()) {
                    continue;
                }

                parsList[phaseParName] = 0.0;
                nameToIndex[phaseParName] = totalNpars;
                parIndexNames.push_back(phaseParName);
                if (isFixed) {
                    fixedParNames.push_back(phaseParName);
                }
                totalNpars++;

                LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase state=%s index=%d value=0",
                    l_value, phaseParName.Data(), isFixed ? "FIXED" : "FREE", totalNpars - 1));
            }
        }
    }

    void MassDependentFitter::ParameterManager::AddMassDependentParametersForL(
        const std::vector<int>& targetL,
        const TString filePath,
        const bool magnitudeOnly,
        const bool isFixed,
        const bool yieldOnly) {

        auto parNames = MassDependentFitter::GetParNames();
        auto& config = fitter.GetMassDependenceConfig();
        auto& hConfig = fitter.GetMomentsConfig();


        std::unique_ptr<TFile> file(TFile::Open(filePath, "READ"));
        if (!file || file->IsZombie()) {
            LogMDWarn("FILE", Form("failed_to_open path=%s", filePath.Data()));
            return;
        }

        TTree* tree = nullptr;
        file->GetObject("mass_dependent_params", tree);
        if (!tree) {
            LogMDWarn("FILE", Form("missing_tree tree=mass_dependent_params path=%s", filePath.Data()));
            file->Close();
            return;
        }

        LogMDInfo("FILE", Form("loading_from_file path=%s", filePath.Data()));
        
        // Read the tree entry once at the beginning
        if (tree->GetEntries() > 0) {
            tree->GetEntry(0);
        }

        for (const TString& parName : parNames) {
            // Parse parameter name to extract l value
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int l_value = l_str.Atoi();
            
            // Check if this l value is in the target list
            if (std::find(targetL.begin(), targetL.end(), l_value) == targetL.end()) {
                LogMDInfo("SKIP", Form("param=%s reason=targetL_mismatch l=%d", parName.Data(), l_value));
                continue;
            }
            
            // Check if this parameter is needed for the H(L,M)s configuration  
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) {
                LogMDInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), l_value));
                continue;
            }
            
            // Only add if this l value should be mass dependent according to config
            if (!config.IsMassDependent(l_value)) {
                LogMDInfo("SKIP", Form("param=%s reason=l_not_mass_dependent l=%d", parName.Data(), l_value));
                continue;
            }

            // Skip phase parameters if magnitudeOnly is true OR yieldOnly is true
            if ((magnitudeOnly || yieldOnly) && parName.Contains("phi")) {
                LogMDInfo("SKIP", Form("param=%s reason=phase_disabled magnitudeOnly=%d yieldOnly=%d",
                    parName.Data(), magnitudeOnly, yieldOnly));
                continue;
            }

            // Get waves for this l value
            std::vector<TString> waves = config.GetWavesForL(l_value);
            if (waves.empty()) {
                LogMDWarn("SKIP", Form("param=%s reason=no_waves_defined_for_l l=%d", parName.Data(), l_value));
                continue;
            }

            for (const TString& waveName : waves) {
                // Get function index for this wave
                int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                
                // Determine which parameter types to add based on yieldOnly
                std::vector<TString> paramTypes;
                if (yieldOnly) {
                    paramTypes = {"k"};  // Only coupling parameters
                    LogMDInfo("MODE", Form("wave=%s yieldOnly=1 params=k", waveName.Data()));
                } else {
                    paramTypes = {"k", "M", "width"};  // All parameters
                }
                
                for (const TString& paramType : paramTypes) {
                    TString name;
                    if (paramType == "M" || paramType == "width") {
                        name = Form("MD_%s_%s", SharedResonanceKeyForWave(waveName).Data(), paramType.Data());
                    } else {
                        name = Form("MD_%s_%s_%s", parName.Data(), waveName.Data(), paramType.Data());
                    }
                    
                    // Check if parameter already exists - skip if already added
                    if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                        LogMDInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                        continue;
                    }

                    // Try to read the parameter value from the tree using TLeaf
                    double value = 0.0;
                    TLeaf* leaf = tree->GetLeaf(name);
                    if (!leaf && (paramType == "M" || paramType == "width")) {
                        TString legacyName = Form("MD_%s_%s_%s", parName.Data(), waveName.Data(), paramType.Data());
                        leaf = tree->GetLeaf(legacyName);
                    }
                    if (leaf) {
                        value = leaf->GetValue();
                        LogMDInfo("LOAD", Form("param=%s source=file value=%g", name.Data(), value));
                    } else {
                        // Parameter not found in file, use values from MassDependentFunction or defaults
                        if (paramType == "k") {
                            value = 0.0;  // Default coupling
                            LogMDInfo("INIT", Form("param=%s source=default field=k value=%g", name.Data(), value));
                        } else if (paramType == "M") {
                            // Get mass from MassDependentFunction
                            if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                                value = fitter.massDepFuncs_[funcIndex].GetResonanceMass();
                                LogMDInfo("INIT", Form("param=%s source=function_%d field=mass value=%g",
                                    name.Data(), funcIndex, value));
                            } else {
                                value = 1.5;  // Fallback default
                                LogMDInfo("INIT", Form("param=%s source=default field=mass value=%g", name.Data(), value));
                            }
                        } else {  // width
                            // Get width from MassDependentFunction
                            if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                                value = fitter.massDepFuncs_[funcIndex].GetResonanceWidth();
                                LogMDInfo("INIT", Form("param=%s source=function_%d field=width value=%g",
                                    name.Data(), funcIndex, value));
                            } else {
                                value = 0.1;  // Fallback default
                                LogMDInfo("INIT", Form("param=%s source=default field=width value=%g", name.Data(), value));
                            }
                        }
                    }

                    parsList[name] = value;
                    nameToIndex[name] = totalNpars;
                    parIndexNames.push_back(name);
                    
                    if (isFixed) {
                        fixedParNames.push_back(name);
                        LogMDInfo("ADD", Form("param=%s state=FIXED index=%d value=%g", name.Data(), totalNpars, value));
                    } else {
                        LogMDInfo("ADD", Form("param=%s state=FREE index=%d value=%g", name.Data(), totalNpars, value));
                    }

                    totalNpars++;
                }
            }
        }

        for (const int l_value : targetL) {
            if (!config.IsMassDependent(l_value)) continue;
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) continue;
            const auto wavesIt = config.massDependentWaves.find(l_value);
            if (wavesIt == config.massDependentWaves.end()) continue;
            for (const TString& configuredWave : wavesIt->second) {
                if (!MassDependenceConfig::HasPhaseTag(configuredWave)) continue;

                const TString waveName = MassDependenceConfig::StripPhaseTag(configuredWave);
                TString phaseParName = Form("MD_%s_globalphase", waveName.Data());
                if (parsList.find(phaseParName) != parsList.end() || nameToIndex.find(phaseParName) != nameToIndex.end()) {
                    continue;
                }

                double value = 0.0;
                if (TLeaf* leaf = tree->GetLeaf(phaseParName)) {
                    value = leaf->GetValue();
                    LogMDInfo("LOAD", Form("param=%s source=file value=%g", phaseParName.Data(), value));
                } else {
                    LogMDInfo("INIT", Form("param=%s source=default field=globalphase value=0", phaseParName.Data()));
                }

                parsList[phaseParName] = value;
                nameToIndex[phaseParName] = totalNpars;
                parIndexNames.push_back(phaseParName);
                if (isFixed) {
                    fixedParNames.push_back(phaseParName);
                }
                totalNpars++;

                LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase state=%s index=%d value=%g",
                    l_value, phaseParName.Data(), isFixed ? "FIXED" : "FREE", totalNpars - 1, value));
            }
        }
        
        file->Close();
    }

    void MassDependentFitter::ParameterManager::AddMassDependentParametersForL(
        const std::vector<TString>& targetL,
        const int seed,
        const bool magnitudeOnly,
        const bool isFixed,
        const bool yieldOnly) {

        auto parNames = MassDependentFitter::GetParNames();
        auto& config = fitter.GetMassDependenceConfig();
        auto& hConfig = fitter.GetMomentsConfig();
        
        randomSeed = seed;
        TRandom3 rng(seed);

        std::set<int> targetLValues;
        for (const TString& amplitudeName : targetL) {
            std::unique_ptr<TObjArray> ampParts(amplitudeName.Tokenize("_"));
            if (ampParts->GetEntries() < 3) continue;
            TString lStr = ((TObjString*)ampParts->At(1))->GetString();
            targetLValues.insert(lStr.Atoi());
        }
        
        for (const TString& parName : parNames) {
            // Parse parameter name to extract l and m values
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString refl_str = ((TObjString*)parts->At(0))->GetString();
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            TString m_str = ((TObjString*)parts->At(2))->GetString();
            
            int l_value = l_str.Atoi();
            
            // Construct amplitude name from parameter name (e.g., "a_1_-1" from "a_1_-1")
            TString amplitudeName = refl_str + "_" + l_str + "_" + m_str;
            
            // Check if this amplitude is in the target list
            if (std::find(targetL.begin(), targetL.end(), amplitudeName) == targetL.end()) {
                LogMDInfo("SKIP", Form("param=%s reason=targetAmplitude_mismatch amplitude=%s", parName.Data(), amplitudeName.Data()));
                continue;
            }
            
            // Check if this parameter is needed for the H(L,M)s configuration  
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) {
                LogMDInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), l_value));
                continue;
            }
            
            // Only add if this l value should be mass dependent according to config
            if (!config.IsMassDependent(l_value)) {
                LogMDInfo("SKIP", Form("param=%s reason=l_not_mass_dependent l=%d", parName.Data(), l_value));
                continue;
            }

            // Skip phase parameters if magnitudeOnly is true OR yieldOnly is true
            if ((magnitudeOnly || yieldOnly) && parName.Contains("phi")) {
                LogMDInfo("SKIP", Form("param=%s reason=phase_disabled magnitudeOnly=%d yieldOnly=%d",
                    parName.Data(), magnitudeOnly, yieldOnly));
                continue;
            }

            // Get waves for this l value
            std::vector<TString> waves = config.GetWavesForL(l_value);
            if (waves.empty()) {
                LogMDWarn("SKIP", Form("param=%s reason=no_waves_defined_for_l l=%d", parName.Data(), l_value));
                continue;
            }

            for (const TString& waveName : waves) {
                // Get function index for this wave
                int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                if (funcIndex < 0 || funcIndex >= static_cast<int>(fitter.massDepFuncs_.size())) {
                    LogMDWarn("INIT", Form("wave=%s reason=missing_function_index using_defaults=1", waveName.Data()));
                }

                // Determine which parameter types to add based on yieldOnly
                std::vector<TString> paramTypes;
                if (yieldOnly) {
                    paramTypes = {"k"};  // Only coupling parameters
                    LogMDInfo("MODE", Form("wave=%s yieldOnly=1 params=k", waveName.Data()));
                } else {
                    paramTypes = {"k", "M", "width"};  // All parameters
                }

                for (const TString& paramType : paramTypes) {
                    TString name;
                    if (paramType == "M" || paramType == "width") {
                        name = Form("MD_%s_%s", SharedResonanceKeyForWave(waveName).Data(), paramType.Data());
                    } else {
                        name = Form("MD_%s_%s_%s", parName.Data(), waveName.Data(), paramType.Data());
                    }
                    
                    // Check if parameter already exists - skip if already added
                    if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                        LogMDInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                        continue;
                    }

                    double initialValue;
                    if (paramType == "k") {
                        initialValue = rng.Uniform(-30, 30);  // Random for coupling strength
                    } else if (paramType == "M") {
                        // Get mass from MassDependentFunction
                        if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                            initialValue = fitter.massDepFuncs_[funcIndex].GetResonanceMass();
                            LogMDInfo("INIT", Form("param=%s source=function_%d field=mass value=%g",
                                name.Data(), funcIndex, initialValue));
                        } else {
                            initialValue = 1.5;  // Fallback default
                            LogMDInfo("INIT", Form("param=%s source=default field=mass value=%g", name.Data(), initialValue));
                        }
                    } else {  // width
                        // Get width from MassDependentFunction
                        if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                            initialValue = fitter.massDepFuncs_[funcIndex].GetResonanceWidth();
                            LogMDInfo("INIT", Form("param=%s source=function_%d field=width value=%g",
                                name.Data(), funcIndex, initialValue));
                        } else {
                            initialValue = 0.1;  // Fallback default
                            LogMDInfo("INIT", Form("param=%s source=default field=width value=%g", name.Data(), initialValue));
                        }
                    }

                    parsList[name] = initialValue;
                    nameToIndex[name] = totalNpars;
                    parIndexNames.push_back(name);
                    
                    if (isFixed) {
                        fixedParNames.push_back(name);
                        LogMDInfo("ADD", Form("scope=%s param=%s state=FIXED index=%d value=%g",
                            amplitudeName.Data(), name.Data(), totalNpars, initialValue));
                    } else {
                        LogMDInfo("ADD", Form("scope=%s param=%s state=FREE index=%d value=%g",
                            amplitudeName.Data(), name.Data(), totalNpars, initialValue));
                    }
                    
                    totalNpars++;
                }
            }
        }

        for (const int l_value : targetLValues) {
            if (!config.IsMassDependent(l_value)) continue;
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) continue;
            const auto wavesIt = config.massDependentWaves.find(l_value);
            if (wavesIt == config.massDependentWaves.end()) continue;
            for (const TString& configuredWave : wavesIt->second) {
                if (!MassDependenceConfig::HasPhaseTag(configuredWave)) continue;

                const TString waveName = MassDependenceConfig::StripPhaseTag(configuredWave);
                TString phaseParName = Form("MD_%s_globalphase", waveName.Data());
                if (parsList.find(phaseParName) != parsList.end() || nameToIndex.find(phaseParName) != nameToIndex.end()) {
                    continue;
                }

                parsList[phaseParName] = 0.0;
                nameToIndex[phaseParName] = totalNpars;
                parIndexNames.push_back(phaseParName);
                if (isFixed) {
                    fixedParNames.push_back(phaseParName);
                }
                totalNpars++;

                LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase state=%s index=%d value=0",
                    l_value, phaseParName.Data(), isFixed ? "FIXED" : "FREE", totalNpars - 1));
            }
        }
    }

    void MassDependentFitter::ParameterManager::AddMassIndependentParametersForL(
        const std::vector<int>& targetL,
        const std::vector<double>& massBins,
        const std::map<std::string, std::pair<double,double>>& init_values) {
        
        LogMIInfo("INIT", Form("source=init_values targetL_count=%zu mass_bins=%zu init_entries=%zu",
            targetL.size(), massBins.size(), init_values.size()));
        
        // Iterate through mass bins
        for (size_t mb_idx = 0; mb_idx < massBins.size(); ++mb_idx) {
            double mass_bin = massBins[mb_idx];
            TString massBinStr = Form("%.6f", mass_bin);
            
            // Iterate through target L values
            for (size_t l_idx = 0; l_idx < targetL.size(); ++l_idx) {
                int l_value = targetL[l_idx];
                TString l_str = Form("%d", l_value);
                
                // Iterate through reflectivities (a, b) and M values (0 to L)
                for (const std::string& refl : {"a", "b"}) {
                    for (int m = 0; m <= l_value; ++m) {
                        // Create the key to look up in init_values map
                        std::string mag_key = refl + "_" + std::to_string(l_value) + "_" + std::to_string(m);
                        
                        // Check if this parameter exists in init_values
                        auto it = init_values.find(mag_key);
                        if (it == init_values.end()) {
                            LogMIWarn("SKIP", Form("param_key=%s reason=missing_in_init_values", mag_key.c_str()));
                            continue;
                        }
                        
                        double magnitude = it->second.first;
                        double phase = it->second.second;
                        
                        // Create magnitude parameter name
                        TString mag_param_name = Form("MI_%s_%s_%d_%d", massBinStr.Data(), refl.c_str(), l_value, m);
                        
                        // Check if parameter already exists
                        if (parsList.find(mag_param_name) != parsList.end() || 
                            nameToIndex.find(mag_param_name) != nameToIndex.end()) {
                            LogMIInfo("SKIP", Form("param=%s reason=already_exists", mag_param_name.Data()));
                            continue;
                        }
                        
                        parsList[mag_param_name] = magnitude;
                        nameToIndex[mag_param_name] = totalNpars;
                        parIndexNames.push_back(mag_param_name);
                        
                        LogMIInfo("ADD", Form("param=%s kind=magnitude index=%d value=%g",
                            mag_param_name.Data(), totalNpars, magnitude));
                        totalNpars++;
                        
                        // Create phase parameter name
                        TString phase_param_name = Form("MI_%s_%sphi_%d_%d", massBinStr.Data(), refl.c_str(), l_value, m);
                        
                        // Check if parameter already exists
                        if (parsList.find(phase_param_name) != parsList.end() || 
                            nameToIndex.find(phase_param_name) != nameToIndex.end()) {
                            LogMIInfo("SKIP", Form("param=%s reason=already_exists", phase_param_name.Data()));
                            continue;
                        }
                        
                        parsList[phase_param_name] = phase;
                        nameToIndex[phase_param_name] = totalNpars;
                        parIndexNames.push_back(phase_param_name);
                        
                        LogMIInfo("ADD", Form("param=%s kind=phase index=%d value=%g",
                            phase_param_name.Data(), totalNpars, phase));
                        totalNpars++;
                    }
                }
            }
        }
        
        LogMIInfo("SUMMARY", Form("total_parameters=%d", totalNpars));
    }

    void MassDependentFitter::ParameterManager::AddMassIndependentParametersForPhase(
        const std::vector<TString>& targetVariables) {

        auto parNames = MassDependentFitter::GetParNames();
        auto& massBins = fitter.GetMassBins();

        for (const double massBin : massBins) {
            for (const TString& parName : parNames) {
                // Check if this parameter name matches any of the target variables
                bool isTargetVariable = false;
                for (const TString& targetVar : targetVariables) {
                    if (parName == targetVar) {
                        isTargetVariable = true;
                        break;
                    }
                }
                
                if (!isTargetVariable) {
                    LogMIInfo("SKIP", Form("param=%s reason=not_in_target_variables", parName.Data()));
                    continue;
                }

                TString name = Form("MI_%1.6f_%s", massBin, parName.Data());
                
                // Check if parameter already exists - skip if already added
                if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                    LogMIInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                    continue;
                }
                
                // Initialize with zero as specified
                double initialValue = 0.0;

                parsList[name] = initialValue;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                totalNpars++;

                LogMIInfo("ADD", Form("param=%s kind=phase index=%d value=%g",
                    name.Data(), totalNpars - 1, initialValue));
            }
        }
    }

    void MassDependentFitter::ParameterManager::AddMassIndependentParametersForPhase(
        const TString& resultTreeFile,
        const std::vector<TString>& targetVariables) {

        auto parNames = MassDependentFitter::GetParNames();
        auto& massBins = fitter.GetMassBins();

        // Open the result file and read the result tree
        std::unique_ptr<TFile> file(TFile::Open(resultTreeFile, "READ"));
        if (!file || file->IsZombie()) {
            LogMIWarn("FILE", Form("failed_to_open path=%s", resultTreeFile.Data()));
            return;
        }

        TTree* tree = nullptr;
        file->GetObject("result", tree);
        if (!tree) {
            LogMIWarn("FILE", Form("missing_tree tree=result path=%s", resultTreeFile.Data()));
            file->Close();
            return;
        }

        LogMIInfo("FILE", Form("loading_from_file path=%s", resultTreeFile.Data()));

        // Set up branch for mass_bin
        double mass_bin_from_file = 0.0;
        tree->SetBranchAddress("mass_bin", &mass_bin_from_file);

        // Create a map to store parameter values for each parameter name
        std::map<TString, double> paramValues;
        
        // Set up branches for all target variables
        for (const TString& targetVar : targetVariables) {
            TBranch* branch = tree->GetBranch(targetVar);
            if (branch) {
                paramValues[targetVar] = 0.0;
                tree->SetBranchAddress(targetVar, &paramValues[targetVar]);
            } else {
                LogMIWarn("FILE", Form("missing_branch branch=%s", targetVar.Data()));
            }
        }

        // Read each entry (one per mass bin) and extract parameter values
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t entry = 0; entry < nEntries; entry++) {
            tree->GetEntry(entry);
            
            // Check if this mass bin is in our target list
            bool foundMassBin = false;
            for (double targetMassBin : massBins) {
                if (TMath::Abs(mass_bin_from_file - targetMassBin) < 1e-6) {
                    foundMassBin = true;
                    break;
                }
            }
            
            if (!foundMassBin) continue;
            
            // Add phase parameters for this mass bin
            for (const TString& parName : parNames) {
                // Check if this parameter name matches any of the target variables
                bool isTargetVariable = false;
                for (const TString& targetVar : targetVariables) {
                    if (parName == targetVar) {
                        isTargetVariable = true;
                        break;
                    }
                }
                
                if (!isTargetVariable) {
                    continue;
                }

                TString name = Form("MI_%1.6f_%s", mass_bin_from_file, parName.Data());
                
                // Check if parameter already exists - skip if already added
                if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                    LogMIInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                    continue;
                }
                
                // Get the value from the tree
                double value = 0.0;
                auto it = paramValues.find(parName);
                if (it != paramValues.end()) {
                    value = it->second;
                } else {
                    LogMIWarn("LOAD", Form("param=%s reason=missing_loaded_value using_default=0", parName.Data()));
                }

                parsList[name] = value;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                totalNpars++;

                LogMIInfo("ADD", Form("mass_bin=%g param=%s kind=phase index=%d value=%g",
                    mass_bin_from_file, name.Data(), totalNpars - 1, value));
            }
        }
        
        file->Close();
    }

    void MassDependentFitter::ParameterManager::AddMassDependentParametersForL(
        const std::vector<TString>& targetL,
        const TString filePath,
        const bool magnitudeOnly,
        const bool isFixed,
        const bool yieldOnly) {

        auto parNames = MassDependentFitter::GetParNames();
        auto& config = fitter.GetMassDependenceConfig();
        auto& hConfig = fitter.GetMomentsConfig();
        
        std::unique_ptr<TFile> file(TFile::Open(filePath, "READ"));
        if (!file || file->IsZombie()) {
            LogMDWarn("FILE", Form("failed_to_open path=%s", filePath.Data()));
            return;
        }

        TTree* tree = nullptr;
        file->GetObject("mass_dependent_params", tree);
        if (!tree) {
            LogMDWarn("FILE", Form("missing_tree tree=mass_dependent_params path=%s", filePath.Data()));
            file->Close();
            return;
        }

        LogMDInfo("FILE", Form("loading_from_file path=%s", filePath.Data()));
        
        // Read the tree entry once at the beginning
        if (tree->GetEntries() > 0) {
            tree->GetEntry(0);
        }

        std::set<int> targetLValues;
        for (const TString& amplitudeName : targetL) {
            std::unique_ptr<TObjArray> ampParts(amplitudeName.Tokenize("_"));
            if (ampParts->GetEntries() < 3) continue;
            TString lStr = ((TObjString*)ampParts->At(1))->GetString();
            targetLValues.insert(lStr.Atoi());
        }

        for (const TString& parName : parNames) {
            // Parse parameter name to extract l and m values
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString refl_str = ((TObjString*)parts->At(0))->GetString();
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            TString m_str = ((TObjString*)parts->At(2))->GetString();
            
            int l_value = l_str.Atoi();
            
            // Construct amplitude name from parameter name (e.g., "a_1_-1" from "a_1_-1")
            TString amplitudeName = refl_str + "_" + l_str + "_" + m_str;
            
            // Check if this amplitude is in the target list
            if (std::find(targetL.begin(), targetL.end(), amplitudeName) == targetL.end()) {
                LogMDInfo("SKIP", Form("param=%s reason=targetAmplitude_mismatch amplitude=%s", parName.Data(), amplitudeName.Data()));
                continue;
            }
            
            // Check if this parameter is needed for the H(L,M)s configuration  
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) {
                LogMDInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), l_value));
                continue;
            }
            
            // Only add if this l value should be mass dependent according to config
            if (!config.IsMassDependent(l_value)) {
                LogMDInfo("SKIP", Form("param=%s reason=l_not_mass_dependent l=%d", parName.Data(), l_value));
                continue;
            }

            // Skip phase parameters if magnitudeOnly is true OR yieldOnly is true
            if ((magnitudeOnly || yieldOnly) && parName.Contains("phi")) {
                LogMDInfo("SKIP", Form("param=%s reason=phase_disabled magnitudeOnly=%d yieldOnly=%d",
                    parName.Data(), magnitudeOnly, yieldOnly));
                continue;
            }

            // Get waves for this l value
            std::vector<TString> waves = config.GetWavesForL(l_value);
            if (waves.empty()) {
                LogMDWarn("SKIP", Form("param=%s reason=no_waves_defined_for_l l=%d", parName.Data(), l_value));
                continue;
            }

            for (const TString& waveName : waves) {
                // Get function index for this wave
                int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                
                // Determine which parameter types to add based on yieldOnly
                std::vector<TString> paramTypes;
                if (yieldOnly) {
                    paramTypes = {"k"};  // Only coupling parameters
                    LogMDInfo("MODE", Form("wave=%s yieldOnly=1 params=k", waveName.Data()));
                } else {
                    paramTypes = {"k", "M", "width"};  // All parameters
                }
                
                for (const TString& paramType : paramTypes) {
                    TString name;
                    if (paramType == "M" || paramType == "width") {
                        name = Form("MD_%s_%s", SharedResonanceKeyForWave(waveName).Data(), paramType.Data());
                    } else {
                        name = Form("MD_%s_%s_%s", parName.Data(), waveName.Data(), paramType.Data());
                    }
                    
                    // Check if parameter already exists - skip if already added
                    if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                        LogMDInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                        continue;
                    }

                    // Try to read the parameter value from the tree using TLeaf
                    double value = 0.0;
                    TLeaf* leaf = tree->GetLeaf(name);
                    if (!leaf && (paramType == "M" || paramType == "width")) {
                        TString legacyName = Form("MD_%s_%s_%s", parName.Data(), waveName.Data(), paramType.Data());
                        leaf = tree->GetLeaf(legacyName);
                    }
                    if (leaf) {
                        value = leaf->GetValue();
                        LogMDInfo("LOAD", Form("param=%s source=file value=%g", name.Data(), value));
                    } else {
                        // Parameter not found in file, use values from MassDependentFunction or defaults
                        if (paramType == "k") {
                            value = 0.0;  // Default coupling
                            LogMDInfo("INIT", Form("param=%s source=default field=k value=%g", name.Data(), value));
                        } else if (paramType == "M") {
                            // Get mass from MassDependentFunction
                            if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                                value = fitter.massDepFuncs_[funcIndex].GetResonanceMass();
                                LogMDInfo("INIT", Form("param=%s source=function_%d field=mass value=%g",
                                    name.Data(), funcIndex, value));
                            } else {
                                value = 1.5;  // Fallback default
                                LogMDInfo("INIT", Form("param=%s source=default field=mass value=%g", name.Data(), value));
                            }
                        } else {  // width
                            // Get width from MassDependentFunction
                            if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                                value = fitter.massDepFuncs_[funcIndex].GetResonanceWidth();
                                LogMDInfo("INIT", Form("param=%s source=function_%d field=width value=%g",
                                    name.Data(), funcIndex, value));
                            } else {
                                value = 0.1;  // Fallback default
                                LogMDInfo("INIT", Form("param=%s source=default field=width value=%g", name.Data(), value));
                            }
                        }
                    }

                    parsList[name] = value;
                    nameToIndex[name] = totalNpars;
                    parIndexNames.push_back(name);
                    
                    if (isFixed) {
                        fixedParNames.push_back(name);
                        LogMDInfo("ADD", Form("param=%s state=FIXED index=%d value=%g", name.Data(), totalNpars, value));
                    } else {
                        LogMDInfo("ADD", Form("param=%s state=FREE index=%d value=%g", name.Data(), totalNpars, value));
                    }

                    totalNpars++;
                }
            }
        }

        for (const int l_value : targetLValues) {
            if (!config.IsMassDependent(l_value)) continue;
            if (!fitter.ParameterNeededForMoments(l_value, hConfig)) continue;
            const auto wavesIt = config.massDependentWaves.find(l_value);
            if (wavesIt == config.massDependentWaves.end()) continue;
            for (const TString& configuredWave : wavesIt->second) {
                if (!MassDependenceConfig::HasPhaseTag(configuredWave)) continue;

                const TString waveName = MassDependenceConfig::StripPhaseTag(configuredWave);
                TString phaseParName = Form("MD_%s_globalphase", waveName.Data());
                if (parsList.find(phaseParName) != parsList.end() || nameToIndex.find(phaseParName) != nameToIndex.end()) {
                    continue;
                }

                double value = 0.0;

                parsList[phaseParName] = value;
                nameToIndex[phaseParName] = totalNpars;
                parIndexNames.push_back(phaseParName);
                if (isFixed) {
                    fixedParNames.push_back(phaseParName);
                }
                totalNpars++;

                LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase state=%s index=%d value=%g",
                    l_value, phaseParName.Data(), isFixed ? "FIXED" : "FREE", totalNpars - 1, value));
            }
        }
        
        file->Close();
    }

    std::map<std::string, std::pair<double, double>> MassDependentFitter::ParameterManager::GetAmplitudeValuesAtMassBins(
        const std::vector<int>& targetL,
        double massBinCenter) {
        
        std::map<std::string, std::pair<double, double>> amplitudes;  // key: "a_l_m" or "b_l_m", value: {magnitude, phase}
        
        // Format mass bin as string with 6 decimal places
        TString massBinStr = Form("%.6f", massBinCenter);
        
        // Now compute amplitudes for each L value
        for (int l_value : targetL) {
            // Iterate through reflectivities (a, b)
            for (const std::string& refl : {"a", "b"}) {
                // Iterate through M values (0 to L)
                for (int m = 0; m <= l_value; ++m) {
                    // Construct the key: "a_0_0", "a_0_1", "b_2_0", etc.
                    std::string key = refl + "_" + std::to_string(l_value) + "_" + std::to_string(m);
                    
                    // Look up magnitude parameter: MI_<mass_bin>_<key>
                    TString mag_param_name = Form("MI_%s_%s", massBinStr.Data(), key.c_str());
                    double magnitude = 0.0;
                    
                    auto mag_it = parsList.find(mag_param_name);
                    if (mag_it != parsList.end()) {
                        magnitude = mag_it->second;
                        LogMIInfo("LOOKUP", Form("param=%s found=1 value=%g", mag_param_name.Data(), magnitude));
                    } else {
                        LogMIWarn("LOOKUP", Form("param=%s found=0 default=0", mag_param_name.Data()));
                    }
                    
                    // Look up phase parameter: MI_<mass_bin>_<refl>phi_<l>_<m>
                    TString phase_param_name = Form("MI_%s_%sphi_%d_%d", massBinStr.Data(), refl.c_str(), l_value, m);
                    double phase = 0.0;
                    
                    auto phase_it = parsList.find(phase_param_name);
                    if (phase_it != parsList.end()) {
                        phase = phase_it->second;
                        LogMIInfo("LOOKUP", Form("param=%s found=1 value=%g", phase_param_name.Data(), phase));
                    } else {
                        LogMIWarn("LOOKUP", Form("param=%s found=0 default=0", phase_param_name.Data()));
                    }
                    
                    // Store the amplitude
                    amplitudes[key] = {magnitude, phase};
                    
                    LogMIInfo("AMPL", Form("key=%s mass_bin=%g magnitude=%g phase=%g",
                        key.c_str(), massBinCenter, magnitude, phase));
                }
            }
        }
        
        return amplitudes;
    }

}
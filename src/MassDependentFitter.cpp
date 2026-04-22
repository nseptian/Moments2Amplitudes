#include "MassDependentFitter.h"
#include "TObjString.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

namespace m2pw {

    namespace {
        TString ModelTypeName(int funcType) {
            switch (funcType) {
                case 0:
                    return "breitwigner";
                case 1:
                    return "two_breitwigner";
                case 2:
                    return "flatte";
                case 3:
                    return "polynomial";
                case 4:
                    return "flatte_plus_polynomial";
                default:
                    return Form("unknown_%d", funcType);
            }
        }

        void LogMDInfo(const TString& stage, const TString& message) {
            if (stage == "SKIP") return;
            std::cout << "[MD][" << stage << "] " << message << std::endl;
        }

        void LogMDWarn(const TString& stage, const TString& message) {
            std::cerr << "[MD][" << stage << "] " << message << std::endl;
        }

        void LogMIInfo(const TString& stage, const TString& message) {
            if (stage == "SKIP") return;
            std::cout << "[MI][" << stage << "] " << message << std::endl;
        }

        void LogMIWarn(const TString& stage, const TString& message) {
            std::cerr << "[MI][" << stage << "] " << message << std::endl;
        }

        void LogFitInfo(const TString& stage, const TString& message) {
            if (stage == "SKIP") return;
            std::cout << "[FIT][" << stage << "] " << message << std::endl;
        }

        void LogFitWarn(const TString& stage, const TString& message) {
            std::cerr << "[FIT][" << stage << "] " << message << std::endl;
        }

        TString GlobalPhaseKeyForWave(const TString& waveName) {
            return waveName;
        }

        TString GlobalPhaseKeyForWave(const TString& waveName, const TString& reflTag) {
            if (reflTag.Length() == 0) {
                return GlobalPhaseKeyForWave(waveName);
            }

            return Form("shared%s_%s", reflTag.Data(), waveName.Data());
        }

        TString ReflectivityTagFromParameterName(const TString& parName) {
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 1) {
                return TString();
            }

            const TString refl = ((TObjString*)parts->At(0))->GetString();
            if (refl == "a" || refl == "b") {
                return refl;
            }

            return TString();
        }

        TString SharedResonanceKeyForWave(const TString& waveName) {
            return Form("shared_%s", waveName.Data());
        }

        TString MassDependentBaseKey(const TString& parName) {
            const TString prefix = "MD_";
            if (!parName.BeginsWith(prefix)) {
                return parName;
            }

            if (parName.Contains("_re_")) {
                const int suffixPos = parName.Index("_re_");
                return parName(prefix.Length(), suffixPos - prefix.Length());
            }

            if (parName.Contains("_im_")) {
                const int suffixPos = parName.Index("_im_");
                return parName(prefix.Length(), suffixPos - prefix.Length());
            }

            const std::vector<TString> suffixes = {
                "_globalphase",
                "_m_threshold",
                "_m_expansion",
                "_g_etapi",
                "_g_KK",
                "_Mass1",
                "_Width1",
                "_Mass2",
                "_Width2",
                "_Mass",
                "_Width",
                "_k1",
                "_k2",
                "_k"
            };

            for (const auto& suffix : suffixes) {
                if (parName.EndsWith(suffix)) {
                    return parName(prefix.Length(), parName.Length() - prefix.Length() - suffix.Length());
                }
            }

            return parName(prefix.Length(), parName.Length() - prefix.Length());
        }

        bool BuildWaveParameters(
            const TString& baseName,
            const TString& waveName,
            int funcType,
            const std::map<TString, std::vector<double>>& massDepPars,
            std::vector<double>& params) {

            TString key = Form("%s_%s", baseName.Data(), waveName.Data());
            const auto it = massDepPars.find(key);
            if (it == massDepPars.end() || it->second.empty()) {
                return false;
            }

            params = it->second;

            if (funcType == 4) {
                if (params.size() < 5) {
                    return false;
                }
                TString sharedKey = SharedResonanceKeyForWave(waveName);
                auto sharedIt = massDepPars.find(sharedKey);
                if (sharedIt == massDepPars.end()) {
                    return false;
                }
                params.push_back(sharedIt->second[0]);
                return true;
            }

            if (funcType == 2) {
                if (params.size() < 3) {
                    return false;
                }

                TString sharedKey = SharedResonanceKeyForWave(waveName);
                auto sharedIt = massDepPars.find(sharedKey);
                if (sharedIt == massDepPars.end()) {
                    return false;
                }

                params.push_back(sharedIt->second[0]);
                return true;
            }

            if (funcType == 1) {
                if (params.size() < 2) {
                    return false;
                }

                auto sharedIt = massDepPars.find(SharedResonanceKeyForWave(waveName));
                if (sharedIt == massDepPars.end() || sharedIt->second.size() < 4) {
                    return false;
                }

                const double M1 = sharedIt->second[0];
                const double width1 = sharedIt->second[1];
                const double M2 = sharedIt->second[2];
                const double width2 = sharedIt->second[3];
                const double k1 = params[0];
                const double k2 = params[1];
                params = {k1, M1, width1, k2, M2, width2};
                return true;
            }

            if (params.size() == 1) {
                auto sharedIt = massDepPars.find(SharedResonanceKeyForWave(waveName));
                if (sharedIt == massDepPars.end() || sharedIt->second.size() < 2) {
                    return false;
                }

                params.push_back(sharedIt->second[0]);
                params.push_back(sharedIt->second[1]);
            }

            return true;
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

        BuildEquationLCache();
    }

    void MassDependentFitter::BuildEquationLCache() {
        equationLCache_.clear();
        for (const auto& [mass_bin, equations] : eqns_) {
            std::vector<int> lValues;
            lValues.reserve(equations.size());
            for (const auto& eqn : equations) {
                lValues.push_back(ExtractLFromEquationName(eqn.GetName()));
            }
            equationLCache_[mass_bin] = std::move(lValues);
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
        
        // Collect all unique wave model configs from the configuration.
        std::map<TString, MassDependenceConfig::WaveModelConfig> uniqueWaveConfigs;
        for (const auto& [ell_value, waveConfigs] : massDependenceConfig_.massDependentWaves) {
                for (const auto& waveConfig : waveConfigs) {
                    auto [it, inserted] = uniqueWaveConfigs.emplace(waveConfig.waveName, waveConfig);
                    if (!inserted) {
                        const bool sameModel = (it->second.modelType == waveConfig.modelType);
                        const bool sameOrder = (it->second.polyOrder == waveConfig.polyOrder);
                        if (!sameModel || !sameOrder) {
                            LogFitWarn("MDFUNC", Form("wave=%s reason=conflicting_model_config keep_first=1",
                                                        waveConfig.waveName.Data()));
                        }
                    }
                }
        }

        // Initialize functions for each unique wave config
        int funcIndex = 0;
        waveToFunctionIndex_.clear();
        
        for (const auto& [waveName, waveConfig] : uniqueWaveConfigs) {
            const int funcType = static_cast<int>(waveConfig.modelType);
            const int polyOrder = waveConfig.polyOrder;

            massDepFuncs_.emplace_back(static_cast<int>(massBins_.size()),
                                       massBins_.front(),
                                       bin_width,
                                       funcType,
                                       polyOrder);
            waveToFunctionIndex_[waveName] = funcIndex;
            if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Polynomial)) {
                LogFitInfo("MDFUNC", Form("index=%d wave=%s model=%s poly_order=%d",
                                           funcIndex,
                                           waveName.Data(),
                                           ModelTypeName(funcType).Data(),
                                           polyOrder));
            } else {
                LogFitInfo("MDFUNC", Form("index=%d wave=%s model=%s",
                                           funcIndex,
                                           waveName.Data(),
                                           ModelTypeName(funcType).Data()));
            }
            funcIndex++;
        }
        
        LogFitInfo("MDFUNC", Form("initialized_count=%zu", massDepFuncs_.size()));
    }

    void MassDependentFitter::SetMassDependenceConfig(const MassDependenceConfig& config) {
        massDependenceConfig_ = config;
     }

    void MassDependentFitter::SetMomentsConfig(const MomentsConfig& config) {
        hMomentsConfig_ = config;
     }

    void MassDependentFitter::SetMinimizerPrintLevel(int level) {
        minimizerPrintLevel_ = level;
            LogFitInfo("MIN", Form("set_print_level=%d", minimizerPrintLevel_));
    }

    void MassDependentFitter::SetMaxFunctionCalls(int maxCalls) {
        if (maxCalls <= 0) {
            LogFitWarn("MIN", Form("set_max_function_calls_failed value=%d reason=non_positive", maxCalls));
            return;
        }

        maxFunctionCalls_ = maxCalls;
        LogFitInfo("MIN", Form("set_max_function_calls=%d", maxFunctionCalls_));
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

        auto& equations = eqns_.at(massBin);
        const auto cacheIt = equationLCache_.find(massBin);

        for (size_t idx = 0; idx < equations.size(); ++idx) {
            auto& eqn = equations[idx];
            int L_value = -1;
            if (cacheIt != equationLCache_.end() && idx < cacheIt->second.size()) {
                L_value = cacheIt->second[idx];
            } else {
                // Fallback when cache is missing/out-of-sync.
                L_value = ExtractLFromEquationName(eqn.GetName());
            }
            
            // Apply H(L,M)s filtering
            bool shouldInclude = true;
            if (L_value >= 0) { // This is an H(L,M) equation
                shouldInclude = hMomentsConfig_.ShouldIncludeL(L_value);
            }
            
            if (shouldInclude) {
                eqn.SetNeedsRecalc();
                double eqn_chi2 = eqn.DoEvalSq(current_vals);
                chi2 += eqn_chi2;
            }
        }

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
            std::map<TString, std::complex<double>> amplitudeCache;

            auto getCachedAmplitude = [&](const TString& base_name, int ell_value) -> const std::complex<double>& {
                const TString cacheKey = MassDependentBaseKey(base_name);
                const auto cacheIt = amplitudeCache.find(cacheKey);
                if (cacheIt != amplitudeCache.end()) {
                    return cacheIt->second;
                }

                auto [insertedIt, _] = amplitudeCache.emplace(
                    cacheKey,
                    GetCombinedAmplitude(mass_bin, base_name, massDepPars, ell_value));
                return insertedIt->second;
            };
            
            // Set parameter values
            for (int j = 0; j < pars.Nvars(); ++j) {
                const TString par_name = pars.GetParName(j);
                const auto param_info = ParseParameterName(par_name);
                
                if (!param_info) continue;

                // l value already parsed by ParameterInfo
                int ell_value = param_info->l_str.Atoi();
                
                double value = 0.0;
                if (!param_info->isPhase) {
                    // Handle magnitude parameters
                    if (const auto* waveConfigs = massDependenceConfig_.GetWaveModel(ell_value);
                        waveConfigs && !waveConfigs->empty()) {
                            value = std::abs(getCachedAmplitude(par_name, ell_value));
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
                    if (!hasPhaseInMassIndepPars) {
                        if (const auto* waveConfigs = massDependenceConfig_.GetWaveModel(ell_value);
                            waveConfigs && !waveConfigs->empty()) {
                            value = std::arg(getCachedAmplitude(base_name, ell_value));
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
        std::vector<double> md_par_errors(nMDPars);
        int idx = 0;
        for (const auto& [parName, value] : paramManager.parsList) {
            // double currentValue = value;
            if (parName.BeginsWith("MD_")) {  // Only mass-dependent parameters
                LogFitInfo("RESULT", Form("tree=mass_dependent_params param=%s value=%g", parName.Data(), value));
                md_par_vals[idx] = value;
                
                // Get error from paramManager using nameToIndex
                double error = 0.0;
                auto it = paramManager.nameToIndex.find(parName);
                if (it != paramManager.nameToIndex.end() && it->second < static_cast<int>(paramManager.parErrors.size())) {
                    error = paramManager.parErrors[it->second];
                    LogFitInfo("RESULT", Form("tree=mass_dependent_params param=%s error=%g", parName.Data(), error));
                }
                md_par_errors[idx] = error;
                
                mdTree->Branch(parName, &md_par_vals[idx]);
                mdTree->Branch(parName + "_error", &md_par_errors[idx]);
                idx++;
            }
        }

        mdTree->Fill();

        
        auto chi2_tree = std::make_unique<TTree>("chi2", "chi2");
        double chi2 = lastChi2_;
        chi2_tree->Branch("chi2", &chi2);
        int seed_val = paramManager.randomSeed;
        bool isValid = minimizerIsValid_;
        bool isValidError = minimizerHasValidErrors_;
        int status = minimizerStatus_;

        LogFitInfo("MIN", Form("seed=%d is_valid=%d is_valid_error=%d status=%d", seed_val, isValid, isValidError, status));

        chi2_tree->Branch("isValid", &isValid);
        chi2_tree->Branch("isValidError", &isValidError);
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

    void MassDependentFitter::ParameterManager::PrintParameters() const {
        LogMDInfo("CONFIG", "registered_parameter_list_begin");
        for (int i = 0; i < totalNpars; ++i) {
            const TString& parName = parIndexNames[i];
            const auto valueIt = parsList.find(parName);
            const auto limitIt = parameterLimits.find(parName);
            const bool isFixed = std::find(fixedParNames.begin(), fixedParNames.end(), parName) != fixedParNames.end();

            const double initialValue = (valueIt != parsList.end()) ? valueIt->second : 0.0;
            TString limitText = "none";
            if (limitIt != parameterLimits.end()) {
                limitText = Form("[%g, %g]", limitIt->second.first, limitIt->second.second);
            }

            LogMDInfo("CONFIG", Form("index=%d name=%s initial=%g limit=%s status=%s",
                i,
                parName.Data(),
                initialValue,
                limitText.Data(),
                isFixed ? "fixed" : "free"));
        }
        LogMDInfo("CONFIG", "registered_parameter_list_end");
    }

    bool MassDependentFitter::ParameterManager::SetParameterInitialValue(const TString& parName, double value) {
        auto printRegisteredParameters = [this]() {
            LogMDWarn("CONFIG", "registered_parameter_list_begin");
            for (const auto& name : parIndexNames) {
                LogMDWarn("CONFIG", Form("registered_parameter=%s", name.Data()));
            }
            LogMDWarn("CONFIG", "registered_parameter_list_end");
        };

        auto it = parsList.find(parName);
        if (it == parsList.end()) {
            LogMDWarn("CONFIG", "Parameter is not found in the parameter list!");
            LogMDWarn("CONFIG", Form("set_initial_failed param=%s reason=not_found", parName.Data()));
            printRegisteredParameters();
            throw std::runtime_error("Parameter is not found in the parameter list!");
        }

        const double oldValue = it->second;
        it->second = value;
        LogMDInfo("CONFIG", Form("set_initial_ok param=%s old=%g new=%g", parName.Data(), oldValue, value));
        return true;
    }

    bool MassDependentFitter::ParameterManager::SetInitialValue(const TString& parName, double value, bool isFixed) {
        auto printRegisteredParameters = [this]() {
            LogMDWarn("CONFIG", "registered_parameter_list_begin");
            for (const auto& name : parIndexNames) {
                LogMDWarn("CONFIG", Form("registered_parameter=%s", name.Data()));
            }
            LogMDWarn("CONFIG", "registered_parameter_list_end");
        };

        auto indexIt = nameToIndex.find(parName);
        if (indexIt == nameToIndex.end() || indexIt->second < 0 || indexIt->second >= static_cast<int>(parIndexNames.size())) {
            LogMDWarn("CONFIG", "Parameter is not found in the parameter list!");
            LogMDWarn("CONFIG", Form("set_initial_failed param=%s reason=not_found_in_added_parameter_list", parName.Data()));
            printRegisteredParameters();
            throw std::runtime_error("Parameter is not found in the parameter list!");
        }

        if (parIndexNames[indexIt->second] != parName) {
            LogMDWarn("CONFIG", Form("set_initial_failed param=%s reason=index_mapping_mismatch", parName.Data()));
            printRegisteredParameters();
            throw std::runtime_error(TString::Format(
                "Parameter '%s' index mapping mismatch.",
                parName.Data()).Data());
        }

        auto valueIt = parsList.find(parName);
        if (valueIt == parsList.end()) {
            LogMDWarn("CONFIG", Form("set_initial_failed param=%s reason=missing_in_value_map", parName.Data()));
            printRegisteredParameters();
            throw std::runtime_error(TString::Format(
                "Parameter '%s' is registered but missing in value map.",
                parName.Data()).Data());
        }

        const double oldValue = valueIt->second;
        valueIt->second = value;
        LogMDInfo("CONFIG", Form("set_initial_update param=%s old=%g new=%g",
            parName.Data(), oldValue, value));

        const auto fixedIt = std::find(fixedParNames.begin(), fixedParNames.end(), parName);
        if (isFixed) {
            if (fixedIt == fixedParNames.end()) {
                fixedParNames.push_back(parName);
            }
            LogMDInfo("CONFIG", Form("set_initial_fix param=%s fixed=1", parName.Data()));
        } else {
            if (fixedIt != fixedParNames.end()) {
                fixedParNames.erase(fixedIt);
            }
            LogMDInfo("CONFIG", Form("set_initial_fix param=%s fixed=0", parName.Data()));
        }

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

    bool MassDependentFitter::ParameterManager::ConfigureBreitWignerDefaultsForWave(
        const TString& waveName,
        double mass,
        double width) {

        if (mass <= 0.0 || width <= 0.0) {
            LogMDWarn("CONFIG", Form("bw_defaults_failed wave=%s reason=invalid_values mass=%g width=%g",
                waveName.Data(), mass, width));
            return false;
        }

        bool updated = false;
        const TString sharedMassName = Form("MD_%s_Mass", SharedResonanceKeyForWave(waveName).Data());
        const TString sharedWidthName = Form("MD_%s_Width", SharedResonanceKeyForWave(waveName).Data());

        auto itM = parsList.find(sharedMassName);
        if (itM != parsList.end()) {
            itM->second = mass;
            updated = true;
        }

        auto itW = parsList.find(sharedWidthName);
        if (itW != parsList.end()) {
            itW->second = width;
            updated = true;
        }

        if (!updated) {
            LogMDWarn("CONFIG", Form("bw_defaults_failed wave=%s reason=shared_parameters_not_found", waveName.Data()));
            return false;
        }

        LogMDInfo("CONFIG", Form("bw_defaults_ok wave=%s mass=%g width=%g", waveName.Data(), mass, width));
        return true;
    }

    int MassDependentFitter::ParameterManager::ConfigureBreitWignerDefaults(
        int ellValue,
        double mass,
        double width) {

        const std::vector<TString> waves = fitter.GetMassDependenceConfig().GetWaves(ellValue);
        int successCount = 0;
        for (const TString& wave : waves) {
            if (ConfigureBreitWignerDefaultsForWave(wave, mass, width)) {
                successCount++;
            }
        }

        LogMDInfo("CONFIG", Form("bw_defaults_for_ell ell=%d attempted=%zu success=%d", ellValue, waves.size(), successCount));
        return successCount;
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

        TString parKey = MassDependentBaseKey(parName);
        massDepPars_[parKey].push_back(value);
    }

    // Static utility method
    std::vector<TString> MassDependentFitter::GetParNames() {
        std::vector<TString> parNames;
        for (int ell = 0; ell <= 2; ++ell) {
            for (int m = -ell; m <= ell; ++m) {
                parNames.push_back(Form("a_%d_%d", ell, m));
                parNames.push_back(Form("b_%d_%d", ell, m));
                parNames.push_back(Form("aphi_%d_%d", ell, m));
                parNames.push_back(Form("bphi_%d_%d", ell, m));
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
            minimizerIsValid_ = false;
            minimizerHasValidErrors_ = false;
            minimizerStatus_ = -1;
            return;
        }

        minimizer->SetPrintLevel(minimizerPrintLevel_);
        LogFitInfo("MIN", Form("minuit_print_level=%d", minimizerPrintLevel_));
        minimizer->SetMaxFunctionCalls(maxFunctionCalls_);
        LogFitInfo("MIN", Form("max_function_calls=%d", maxFunctionCalls_));

        Chi2Function chi2Function(*this, paramManager.parIndexNames);
        ROOT::Math::Functor functor(chi2Function, paramManager.totalNpars);
        minimizer->SetFunction(functor);

        for (int i = 0; i < paramManager.totalNpars; ++i) {
            const TString& parName = paramManager.parIndexNames[i];
            const bool fixed = std::find(paramManager.fixedParNames.begin(),
                                         paramManager.fixedParNames.end(),
                                         parName) != paramManager.fixedParNames.end();
            const double value = initialValues[i];

            // cout << "Setting parameter " << parName << " (index " << i << ") to initial value " << value 
            //      << (fixed ? " [fixed]" : "") << endl;

            if (fixed) {
                minimizer->SetFixedVariable(i, parName.Data(), value);
            } else {
                auto limitIt = paramManager.parameterLimits.find(parName);
                if (limitIt != paramManager.parameterLimits.end()) {
                    minimizer->SetLimitedVariable(i,
                                                  parName.Data(),
                                                  value,
                                                  0.01,
                                                  limitIt->second.first,
                                                  limitIt->second.second);
                } else {
                    minimizer->SetVariable(i, parName.Data(), value, 0.01);
                }
            }
        }

        const bool minimizerConverged = minimizer->Minimize();

        // isValid should represent minimization success, not covariance validity.
        minimizerIsValid_ = minimizerConverged && (minimizer->Status() == 0);
        minimizerHasValidErrors_ = minimizer->IsValidError();
        minimizerStatus_ = minimizer->Status();
        lastChi2_ = minimizer->MinValue();

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

    // Helper method for amplitude magnitude evaluation (single or coherent waves)
    std::complex<double> MassDependentFitter::GetCombinedAmplitude(
        double mass_bin,
        const TString& par_name,
        const std::map<TString, std::vector<double>>& massDepPars,
        int ell_value) const {

        std::complex<double> total_amplitude(0.0, 0.0);
        const std::vector<TString> waveNames = massDependenceConfig_.GetWaves(ell_value);

        for (const TString& waveName : waveNames) {
            const TString key = Form("%s_%s", par_name.Data(), waveName.Data());
            const auto it = massDepPars.find(key);
            if (it == massDepPars.end() || it->second.empty()) {
                continue;
            }

            const int funcIndex = GetFunctionIndexForWave(waveName);
            if (funcIndex < 0 || funcIndex >= static_cast<int>(massDepFuncs_.size())) {
                continue;
            }

            std::vector<double> params;
            const int funcType = massDepFuncs_[funcIndex].GetFuncType();
            if (!BuildWaveParameters(par_name, waveName, funcType, massDepPars, params)) {
                continue;
            }

            const TString reflTag = ReflectivityTagFromParameterName(par_name);
            auto phaseIt = massDepPars.find(GlobalPhaseKeyForWave(waveName, reflTag));
            if (phaseIt == massDepPars.end() && reflTag.Length() > 0) {
                phaseIt = massDepPars.find(GlobalPhaseKeyForWave(waveName));
            }

            const bool includeGlobalPhase = (phaseIt != massDepPars.end() && !phaseIt->second.empty());
            if (includeGlobalPhase) {
                params.push_back(phaseIt->second[0]);
            }

            total_amplitude += massDepFuncs_[funcIndex].GetAmplitude(mass_bin, params, includeGlobalPhase);
        }

        return total_amplitude;
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

    void MassDependentFitter::ParameterManager::AddMassIndependentParameter(
        const std::vector<int>& targetELL,
        const int seed) {
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
                int ell_value = l_str.Atoi();

                const bool targetRequested =
                    std::find(targetELL.begin(), targetELL.end(), ell_value) != targetELL.end();
                const bool forceFixedZero = targetRequested && ell_value == 1;

                if (!targetRequested) {
                    continue;
                }

                const bool neededForMoments = fitter.ParameterNeededForMoments(ell_value, hConfig);
                if (!neededForMoments && !forceFixedZero) {
                    continue;
                }
                
                // Check if this parameter is needed for the H(L,M)s configuration  
                if (!neededForMoments) {
                    LogMIInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), ell_value));
                }
                
                // Only add if this l value should be mass independent
                if (!config.IsMassIndependent(ell_value)) continue;

                TString name = Form("MI_%1.6f_%s", massBin, parName.Data());
                
                // Check if parameter already exists - skip if already added
                if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                    LogMIInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                    continue;
                }
                
                double initialValue = forceFixedZero ? 0.0 : (name.Contains("phi") ? 0.0 : rng.Uniform(-200, 200));

                parsList[name] = initialValue;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                if (forceFixedZero) {
                    fixedParNames.push_back(name);
                }
                totalNpars++;

                LogMIInfo("ADD", Form("param=%s l=%d index=%d value=%g",
                    name.Data(), ell_value, totalNpars - 1, initialValue));
            }
        }
    }

    void MassDependentFitter::ParameterManager::AddMassIndependentParameter(
        const std::vector<int>& targetELL,
        const TString filePath) {
        auto parNames = MassDependentFitter::GetParNames();
        auto& massBins = fitter.GetMassBins();
        
        // Load mass-independent parameters from result tree file
        TFile file(filePath, "READ");
        if (!file.IsOpen() || file.IsZombie()) {
            LogMIWarn("FILE", Form("failed_to_open path=%s", filePath.Data()));
            return;
        }
        
        TTree* tree = dynamic_cast<TTree*>(file.Get("result"));
        if (!tree) {
            LogMIWarn("FILE", Form("missing_tree tree=result path=%s", filePath.Data()));
            return;
        }

        LogMIInfo("FILE", Form("loading_from_file path=%s selected_mass_bins=%zu", filePath.Data(), massBins.size()));

        // Convert targetELL to set for faster lookup
        std::set<int> targetELLSet(targetELL.begin(), targetELL.end());
        
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
            int ell_value = l_str.Atoi();
            
            // Check if this l value is in our target list
            if (targetELLSet.find(ell_value) == targetELLSet.end()) {
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
                int ell_value = l_str.Atoi();
                
                LogMIInfo("ADD", Form("param=%s l=%d mass_bin=%g index=%d value=%g",
                    name.Data(), ell_value, mass_bin_from_file, totalNpars, value));
                totalNpars++;
            }
        }
        
        file.Close();
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
        lastChi2_ = total_chi2;
        
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

    bool MassDependentFitter::ParameterNeededForMoments(int ell, const MomentsConfig& hConfig) const {
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
            return (ell == 2);
        }
        
        // Rule 2: If L=2 and L=4 H(L,M)s are included -> need l=0 and l=2 amplitude parameters
        if (includedMomentL.size() == 2 && 
            includedMomentL.count(2) == 1 && includedMomentL.count(4) == 1) {
            return (ell == 0 || ell == 2);
        }
        
        // Rule 3: If all L moments are included -> need all parameters l=0,1,2 (assuming lmax=2)
        if (hConfig.includeAll) {
            return (ell <= 2);  // Include l=0,1,2 (assuming lmax=2)
        }
        
        // Default: if none of the specific cases match, include all parameters
        return (ell <= 2);
    }

    void MassDependentFitter::ParameterManager::AddMassDependentParameter(
        const std::vector<int>& targetELL,
        const int seed,
        const int option) {

        const bool includePhaseParameters = false;

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
            int ell_value = l_str.Atoi();
            
            // Check if this l value is in the target list
            if (std::find(targetELL.begin(), targetELL.end(), ell_value) == targetELL.end()) {
                LogMDInfo("SKIP", Form("param=%s reason=targetELL_mismatch ell=%d", parName.Data(), ell_value));
                continue;
            }
            
            // Check if this parameter is needed for the H(L,M)s configuration  
            if (!fitter.ParameterNeededForMoments(ell_value, hConfig)) {
                LogMDInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), ell_value));
                continue;
            }
            
            // Only add if this l value should be mass dependent according to config
            if (!config.IsMassDependent(ell_value)) {
                LogMDInfo("SKIP", Form("param=%s reason=l_not_mass_dependent l=%d", parName.Data(), ell_value));
                continue;
            }

            if (!includePhaseParameters && parName.Contains("phi")) {
                LogMDInfo("SKIP", Form("param=%s reason=phase_disabled option=%d", parName.Data(), option));
                continue;
            }

            // Get waves for this l value
            std::vector<TString> waves = config.GetWaves(ell_value);
            if (waves.empty()) {
                LogMDWarn("SKIP", Form("param=%s reason=no_waves_defined_for_l l=%d", parName.Data(), ell_value));
                continue;
            }

            for (const TString& waveName : waves) {
                // Get function index for this wave
                int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                if (funcIndex < 0 || funcIndex >= static_cast<int>(fitter.massDepFuncs_.size())) {
                    LogMDWarn("INIT", Form("wave=%s reason=missing_function_index using_defaults=1", waveName.Data()));
                }

                // Check if this wave uses conformal polynomial coefficients.
                bool isConformal = false;
                bool isFlattePlusConformal = false;
                int polyOrder = 2;
                int funcType = -1;
                if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                    funcType = fitter.massDepFuncs_[funcIndex].GetFuncType();
                    isConformal = (funcType == static_cast<int>(MassDependenceConfig::ModelType::Polynomial));
                    isFlattePlusConformal =
                        (funcType == static_cast<int>(MassDependenceConfig::ModelType::FlattePlusPolynomial));
                    if (isConformal || isFlattePlusConformal) {
                        polyOrder = fitter.massDepFuncs_[funcIndex].GetPolyOrder();
                    }
                }

                // Determine which parameter types to add
                std::vector<TString> paramTypes;
                if (isConformal) {
                    // Conformal polynomial: add re_i and im_i for i=0 to polyOrder
                    // Parameters: m_threshold, m_expansion, re_i, im_i
                    paramTypes = {"m_threshold", "m_expansion"};
                    for (int i = 0; i <= polyOrder; ++i) {
                        paramTypes.push_back(Form("re_%d", i));
                        paramTypes.push_back(Form("im_%d", i));
                    }
                    LogMDInfo("MODE", Form("wave=%s model=conformal order=%d n_params=%zu",
                        waveName.Data(), polyOrder, paramTypes.size()));
                } else if (isFlattePlusConformal) {
                    // Flatte + conformal polynomial:
                    // [k, g_etapi, g_KK, m_threshold, m_expansion, re_0, im_0, ..., re_N, im_N, Mass]
                    paramTypes = {"k", "g_etapi", "g_KK", "m_threshold", "m_expansion"};
                    for (int i = 0; i <= polyOrder; ++i) {
                        paramTypes.push_back(Form("re_%d", i));
                        paramTypes.push_back(Form("im_%d", i));
                    }
                    paramTypes.push_back("Mass");
                    LogMDInfo("MODE", Form("wave=%s model=flatte_plus_conformal order=%d n_params=%zu",
                        waveName.Data(), polyOrder, paramTypes.size()));
                } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::TwoBreitWigner)) {
                    paramTypes = {"k1", "Mass1", "Width1", "k2", "Mass2", "Width2"};
                } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Flatte)) {
                    paramTypes = {"k", "g_etapi", "g_KK", "Mass"};
                }
                else {
                    paramTypes = {"k", "Mass", "Width"};  // All parameters
                }

                for (const TString& paramType : paramTypes) {
                    TString name;
                    const bool isSharedResonanceParam = paramType.BeginsWith("Mass") || paramType.BeginsWith("Width");
                    if (isSharedResonanceParam && !isConformal) {
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
                    if (isConformal || isFlattePlusConformal) {
                        // Conformal polynomial parameters
                        if (paramType.BeginsWith("re_") || paramType.BeginsWith("im_")) {
                            // Extract coefficient index
                            TString indexStr = paramType;
                            indexStr.Remove(0, 3);
                            int coeffIndex = indexStr.Atoi();
                            
                            // Scale down higher order terms
                            double scale = TMath::Power(0.5, coeffIndex);
                            initialValue = rng.Uniform(-50, 50) * scale;
                        } else if (paramType == "k") {
                            initialValue = rng.Uniform(-30, 30);
                        } else if (paramType == "g_etapi") {
                            initialValue = 0.5;
                        } else if (paramType == "g_KK") {
                            initialValue = 0.3;
                        } else if (paramType == "m_threshold") {
                            initialValue = 0.683;
                        } else if (paramType == "m_expansion") {
                            initialValue = 0.98;
                        } else {
                            initialValue = 0.0;
                        }
                    } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Flatte) && paramType == "g_etapi") {
                        initialValue = 0.5;
                    } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Flatte) && paramType == "g_KK") {
                        initialValue = 0.3;
                    } else if (paramType.BeginsWith("k")) {
                        initialValue = rng.Uniform(-30, 30);  // Random for coupling strength
                    } else if (paramType == "Mass1") {
                        initialValue = 0.0;
                        LogMDInfo("INIT", Form("param=%s source=default field=Mass1 value=%g", name.Data(), initialValue));
                    } else if (paramType == "Width1") {
                        initialValue = 0.0;
                        LogMDInfo("INIT", Form("param=%s source=default field=Width1 value=%g", name.Data(), initialValue));
                    } else if (paramType == "Mass2") {
                        initialValue = 0.0;
                        LogMDInfo("INIT", Form("param=%s source=default field=Mass2 value=%g", name.Data(), initialValue));
                    } else if (paramType == "Width2") {
                        initialValue = 0.0;
                        LogMDInfo("INIT", Form("param=%s source=default field=Width2 value=%g", name.Data(), initialValue));
                    } else if (paramType == "Mass") {
                        // Use model-specific default mass if file/config value is absent.
                        if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Flatte)) {
                            initialValue = 1.001;
                        } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Polynomial) ||
                                   funcType == static_cast<int>(MassDependenceConfig::ModelType::FlattePlusPolynomial)) {
                            initialValue = 0.98;
                        } else {
                            initialValue = 0.0;
                        }
                        LogMDInfo("INIT", Form("param=%s source=default field=Mass value=%g", name.Data(), initialValue));
                    } else {  // Width
                        // Use model-specific default width/threshold scale if file/config value is absent.
                        if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Flatte)) {
                            initialValue = 0.075;
                        } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Polynomial) ||
                                   funcType == static_cast<int>(MassDependenceConfig::ModelType::FlattePlusPolynomial)) {
                            initialValue = 0.683;
                        } else {
                            initialValue = 0.0;
                        }
                        LogMDInfo("INIT", Form("param=%s source=default field=Width value=%g", name.Data(), initialValue));
                    }

                    parsList[name] = initialValue;
                    nameToIndex[name] = totalNpars;
                    parIndexNames.push_back(name);

                    LogMDInfo("ADD", Form("scope=L%d param=%s state=FREE index=%d value=%g",
                        ell_value, name.Data(), totalNpars, initialValue));
                    
                    totalNpars++;
                }
            }
        }

        if (option == 1 || option == 2) {
            for (const int ell_value : targetELL) {
                if (!config.IsMassDependent(ell_value)) continue;
                if (!fitter.ParameterNeededForMoments(ell_value, hConfig)) continue;

                const std::vector<TString> waves = config.GetWaves(ell_value);
                for (const TString& waveName : waves) {
                    if (option == 1) {
                        TString phaseParName = Form("MD_%s_globalphase", waveName.Data());
                        if (parsList.find(phaseParName) != parsList.end() || nameToIndex.find(phaseParName) != nameToIndex.end()) {
                            continue;
                        }

                        parsList[phaseParName] = 0.0;
                        nameToIndex[phaseParName] = totalNpars;
                        parIndexNames.push_back(phaseParName);
                        SetParameterLimits(phaseParName, -TMath::Pi(), TMath::Pi());
                        totalNpars++;

                        LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase option=%d state=FREE index=%d value=0",
                            ell_value, phaseParName.Data(), option, totalNpars - 1));
                        continue;
                    }

                    for (const TString& reflTag : {TString("a"), TString("b")}) {
                        TString phaseParName = Form("MD_shared%s_%s_globalphase", reflTag.Data(), waveName.Data());
                        if (parsList.find(phaseParName) != parsList.end() || nameToIndex.find(phaseParName) != nameToIndex.end()) {
                            continue;
                        }

                        parsList[phaseParName] = 0.0;
                        nameToIndex[phaseParName] = totalNpars;
                        parIndexNames.push_back(phaseParName);
                        SetParameterLimits(phaseParName, -TMath::Pi(), TMath::Pi());
                        totalNpars++;

                        LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase reflectivity=%s option=%d state=FREE index=%d value=0",
                            ell_value, phaseParName.Data(), reflTag.Data(), option, totalNpars - 1));
                    }
                }
            }
        }
    }

    void MassDependentFitter::ParameterManager::AddMassDependentParameter(
        const std::vector<int>& targetELL,
        const TString filePath,
        const int option) {

        const bool includePhaseParameters = false;

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
            int ell_value = l_str.Atoi();
            
            // Check if this l value is in the target list
            if (std::find(targetELL.begin(), targetELL.end(), ell_value) == targetELL.end()) {
                LogMDInfo("SKIP", Form("param=%s reason=targetELL_mismatch ell=%d", parName.Data(), ell_value));
                continue;
            }
            
            // Check if this parameter is needed for the H(L,M)s configuration  
            if (!fitter.ParameterNeededForMoments(ell_value, hConfig)) {
                LogMDInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), ell_value));
                continue;
            }
            
            // Only add if this l value should be mass dependent according to config
            if (!config.IsMassDependent(ell_value)) {
                LogMDInfo("SKIP", Form("param=%s reason=l_not_mass_dependent l=%d", parName.Data(), ell_value));
                continue;
            }

            if (!includePhaseParameters && parName.Contains("phi")) {
                LogMDInfo("SKIP", Form("param=%s reason=phase_disabled option=%d", parName.Data(), option));
                continue;
            }

            // Get waves for this l value
            std::vector<TString> waves = config.GetWaves(ell_value);
            if (waves.empty()) {
                LogMDWarn("SKIP", Form("param=%s reason=no_waves_defined_for_l l=%d", parName.Data(), ell_value));
                continue;
            }

            for (const TString& waveName : waves) {
                // Get function index for this wave
                int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                int funcType = (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size()))
                    ? fitter.massDepFuncs_[funcIndex].GetFuncType()
                    : -1;
                
                // Determine which parameter types to add from option.
                std::vector<TString> paramTypes;
                if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Polynomial)) {
                    int polyOrder = fitter.massDepFuncs_[funcIndex].GetPolyOrder();
                    paramTypes = {"m_threshold", "m_expansion"};
                    for (int i = 0; i <= polyOrder; ++i) {
                        paramTypes.push_back(Form("re_%d", i));
                        paramTypes.push_back(Form("im_%d", i));
                    }
                } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::FlattePlusPolynomial)) {
                    int polyOrder = fitter.massDepFuncs_[funcIndex].GetPolyOrder();
                    paramTypes = {"k", "g_etapi", "g_KK", "m_threshold", "m_expansion"};
                    for (int i = 0; i <= polyOrder; ++i) {
                        paramTypes.push_back(Form("re_%d", i));
                        paramTypes.push_back(Form("im_%d", i));
                    }
                    paramTypes.push_back("Mass");
                } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Flatte)) {
                    paramTypes = {"k", "g_etapi", "g_KK", "Mass"};
                } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::TwoBreitWigner)) {
                    paramTypes = {"k1", "Mass1", "Width1", "k2", "Mass2", "Width2"};
                } else {
                    paramTypes = {"k", "Mass", "Width"};  // All parameters
                }
                
                for (const TString& paramType : paramTypes) {
                    TString name;
                    const bool isSharedResonanceParam = paramType.BeginsWith("Mass") || paramType.BeginsWith("Width");
                    if (isSharedResonanceParam && 
                        funcType != static_cast<int>(MassDependenceConfig::ModelType::Polynomial)) {
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
                    if (!leaf && isSharedResonanceParam) {
                        TString legacyParamType = paramType;
                        legacyParamType.ReplaceAll("Mass", "M");
                        legacyParamType.ReplaceAll("Width", "width");
                        TString legacyName = Form("MD_%s_%s_%s", parName.Data(), waveName.Data(), legacyParamType.Data());
                        leaf = tree->GetLeaf(legacyName);
                    }
                    if (leaf) {
                        value = leaf->GetValue();
                        LogMDInfo("LOAD", Form("param=%s source=file value=%g", name.Data(), value));
                    } else {
                        // Parameter not found in file, use values from MassDependentFunction or defaults
                        if (paramType == "k") {
                            value = 1.0;
                            LogMDInfo("INIT", Form("param=%s source=default field=k value=%g", name.Data(), value));
                        } else if (paramType == "g_etapi") {
                            value = 0.5;
                            LogMDInfo("INIT", Form("param=%s source=default field=g_etapi value=%g", name.Data(), value));
                        } else if (paramType == "g_KK") {
                            value = 0.3;
                            LogMDInfo("INIT", Form("param=%s source=default field=g_KK value=%g", name.Data(), value));
                        } else if (paramType == "m_threshold") {
                            value = 0.683;
                            LogMDInfo("INIT", Form("param=%s source=default field=m_threshold value=%g", name.Data(), value));
                        } else if (paramType == "m_expansion") {
                            value = 0.98;
                            LogMDInfo("INIT", Form("param=%s source=default field=m_expansion value=%g", name.Data(), value));
                        } else if (paramType.BeginsWith("re_") || paramType.BeginsWith("im_")) {
                            value = 0.0;
                            LogMDInfo("INIT", Form("param=%s source=default field=conformal_coeff value=%g", name.Data(), value));
                        } else if (paramType.BeginsWith("k")) {
                            value = 0.0;  // Default coupling
                            LogMDInfo("INIT", Form("param=%s source=default field=k value=%g", name.Data(), value));
                        } else if (paramType == "Mass1") {
                            value = 0.0;
                            LogMDInfo("INIT", Form("param=%s source=default field=Mass1 value=%g", name.Data(), value));
                        } else if (paramType == "Width1") {
                            value = 0.0;
                            LogMDInfo("INIT", Form("param=%s source=default field=Width1 value=%g", name.Data(), value));
                        } else if (paramType == "Mass2") {
                            value = 0.0;
                            LogMDInfo("INIT", Form("param=%s source=default field=Mass2 value=%g", name.Data(), value));
                        } else if (paramType == "Width2") {
                            value = 0.0;
                            LogMDInfo("INIT", Form("param=%s source=default field=Width2 value=%g", name.Data(), value));
                        } else if (paramType == "Mass") {
                            // Use model-specific default mass if file/config value is absent.
                            if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Flatte)) {
                                value = 1.001;
                            } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Polynomial) ||
                                       funcType == static_cast<int>(MassDependenceConfig::ModelType::FlattePlusPolynomial)) {
                                value = 0.98;
                            } else {
                                value = 0.0;
                            }
                            LogMDInfo("INIT", Form("param=%s source=default field=Mass value=%g", name.Data(), value));
                        } else {  // Width
                            // Use model-specific default width/threshold scale if file/config value is absent.
                            if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Flatte)) {
                                value = 0.075;
                            } else if (funcType == static_cast<int>(MassDependenceConfig::ModelType::Polynomial) ||
                                       funcType == static_cast<int>(MassDependenceConfig::ModelType::FlattePlusPolynomial)) {
                                value = 0.683;
                            } else {
                                value = 0.0;
                            }
                            LogMDInfo("INIT", Form("param=%s source=default field=Width value=%g", name.Data(), value));
                        }
                    }

                    parsList[name] = value;
                    nameToIndex[name] = totalNpars;
                    parIndexNames.push_back(name);

                    LogMDInfo("ADD", Form("param=%s state=FREE index=%d value=%g", name.Data(), totalNpars, value));

                    totalNpars++;
                }
            }
        }

        if (option == 1 || option == 2) {
            for (const int ell_value : targetELL) {
                if (!config.IsMassDependent(ell_value)) continue;
                if (!fitter.ParameterNeededForMoments(ell_value, hConfig)) continue;

                const std::vector<TString> waves = config.GetWaves(ell_value);
                for (const TString& waveName : waves) {
                    if (option == 1) {
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
                        SetParameterLimits(phaseParName, -TMath::Pi(), TMath::Pi());
                        totalNpars++;

                        LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase option=%d state=FREE index=%d value=%g",
                            ell_value, phaseParName.Data(), option, totalNpars - 1, value));
                        continue;
                    }

                    for (const TString& reflTag : {TString("a"), TString("b")}) {
                        TString phaseParName = Form("MD_shared%s_%s_globalphase", reflTag.Data(), waveName.Data());
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
                        SetParameterLimits(phaseParName, -TMath::Pi(), TMath::Pi());
                        totalNpars++;

                        LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase reflectivity=%s option=%d state=FREE index=%d value=%g",
                            ell_value, phaseParName.Data(), reflTag.Data(), option, totalNpars - 1, value));
                    }
                }
            }
        }
        
        file->Close();
    }

    std::map<std::string, std::pair<double, double>> MassDependentFitter::ParameterManager::GetAmplitudeValuesAtMassBins(
        const std::vector<int>& targetELL,
        double massBinCenter) {
        
        std::map<std::string, std::pair<double, double>> amplitudes;  // key: "a_l_m" or "b_l_m", value: {magnitude, phase}
        
        // Format mass bin as string with 6 decimal places
        TString massBinStr = Form("%.6f", massBinCenter);
        
        // Now compute amplitudes for each L value
        for (int ell_value : targetELL) {
            // Iterate through reflectivities (a, b)
            for (const std::string& refl : {"a", "b"}) {
                // Iterate through M values (0 to L)
                for (int m = 0; m <= ell_value; ++m) {
                    // Construct the key: "a_0_0", "a_0_1", "b_2_0", etc.
                    std::string key = refl + "_" + std::to_string(ell_value) + "_" + std::to_string(m);
                    
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
                    TString phase_param_name = Form("MI_%s_%sphi_%d_%d", massBinStr.Data(), refl.c_str(), ell_value, m);
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
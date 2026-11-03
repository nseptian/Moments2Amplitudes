#ifndef PARAMETER_MANAGER_CPP
#define PARAMETER_MANAGER_CPP

#include "MassDependentFitter.h"
#include "FitInternalUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TMath.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <set>

namespace m2pw {

    bool MassDependentFitter::ParameterManager::SetGaussianConstraint(const TString& parName, double mean, double sigma) {
        if (sigma <= 0.0) {
            LogMDWarn("CONFIG", Form("set_gauss_constraint_failed param=%s reason=non_positive_sigma sigma=%g", parName.Data(), sigma));
            return false;
        }

        auto it = parsList.find(parName);
        if (it == parsList.end()) {
            LogMDWarn("CONFIG", Form("set_gauss_constraint_failed param=%s reason=not_found", parName.Data()));
            return false;
        }

        gaussianConstraints[parName] = {mean, sigma};
        LogMDInfo("CONFIG", Form("set_gauss_constraint_ok param=%s mean=%g sigma=%g", parName.Data(), mean, sigma));
        return true;
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

            TString constraintText = "none";
            const auto gIt = gaussianConstraints.find(parName);
            if (gIt != gaussianConstraints.end()) {
                constraintText = Form("gauss(mean=%g,sigma=%g)", gIt->second.first, gIt->second.second);
            }

            LogMDInfo("CONFIG", Form("index=%d name=%s initial=%g limit=%s constraint=%s status=%s",
                i,
                parName.Data(),
                initialValue,
                limitText.Data(),
                constraintText.Data(),
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
        
        for (int i = 0; i < totalNpars; ++i) {
            parErrors[i] = errors[i];
        }
        
        for (int i = 0; i < totalNpars; ++i) {
            const TString& parName = parIndexNames[i];
            
            if (parName.BeginsWith("MI_")) {
                int firstUnderscore = parName.Index("_", 3);
                if (firstUnderscore > 0) {
                    TString massBinStr = parName(3, firstUnderscore - 3);
                    double massBin = massBinStr.Atof();
                    TString actualParName = parName(firstUnderscore + 1, parName.Length() - firstUnderscore - 1);
                    
                    auto parsIt = fitter.pars_.find(massBin);
                    if (parsIt != fitter.pars_.end()) {
                        parsIt->second.SetCurrentError(actualParName, errors[i]);
                    }
                }
            }
        }
    }

    void MassDependentFitter::ParameterManager::AddMassIndependentParameter(
        const std::vector<int>& targetELL,
        const int seed) {
        randomSeed = seed;
        
        auto parNames = MassDependentFitter::GetParNames();
        auto& config = fitter.GetMassDependenceConfig();
        auto& hConfig = fitter.GetMomentsConfig();
        auto& massBins = fitter.GetMassBins();
        
        rng.SetSeed(seed);
        for (const double massBin : massBins) {
            for (const TString& parName : parNames) {
                std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
                if (parts->GetEntries() < 3) continue;
                
                TString l_str = ((TObjString*)parts->At(1))->GetString();
                int ell_value = l_str.Atoi();

                const bool targetRequested =
                    std::find(targetELL.begin(), targetELL.end(), ell_value) != targetELL.end();

                if (!targetRequested) {
                    continue;
                }

                WaveInitRule initRule;
                auto ruleIt = waveInitRules.find(ell_value);
                if (ruleIt != waveInitRules.end()) {
                    initRule = ruleIt->second;
                } else {
                    if (ell_value == 1) {
                        initRule.strategy = WaveInitRule::Strategy::FixedZero;
                        initRule.isFixed = true;
                    } else {
                        initRule.strategy = WaveInitRule::Strategy::RandomPositiveMagnitudeZeroPhase;
                        initRule.minVal = 0.0;
                        initRule.maxVal = 200.0;
                    }
                }

                const bool isFixed = initRule.isFixed || (initRule.strategy == WaveInitRule::Strategy::FixedZero) || (initRule.strategy == WaveInitRule::Strategy::FixedValue);

                const bool neededForMoments = fitter.ParameterNeededForMoments(ell_value, hConfig);
                if (!neededForMoments && !isFixed) {
                    continue;
                }
                
                if (!neededForMoments) {
                    LogMIInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), ell_value));
                }
                
                if (!config.IsMassIndependent(ell_value)) continue;

                TString name = Form("MI_%1.6f_%s", massBin, parName.Data());
                
                if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                    LogMIInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                    continue;
                }
                
                const bool isPhase = name.Contains("phi");
                double initialValue = 0.0;
                switch (initRule.strategy) {
                    case WaveInitRule::Strategy::FixedZero:
                        initialValue = 0.0;
                        break;
                    case WaveInitRule::Strategy::FixedValue:
                        initialValue = initRule.fixedValue;
                        break;
                    case WaveInitRule::Strategy::RandomPositiveMagnitudeZeroPhase:
                        initialValue = isPhase ? 0.0 : rng.Uniform(initRule.minVal, initRule.maxVal);
                        break;
                    case WaveInitRule::Strategy::RandomFullPhase:
                        initialValue = isPhase ? rng.Uniform(-TMath::Pi(), TMath::Pi()) : rng.Uniform(initRule.minVal, initRule.maxVal);
                        break;
                    default:
                        initialValue = isPhase ? 0.0 : rng.Uniform(initRule.minVal, initRule.maxVal);
                        break;
                }

                parsList[name] = initialValue;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                if (isFixed) {
                    fixedParNames.push_back(name);
                }
                totalNpars++;

                LogMIInfo("ADD", Form("param=%s l=%d index=%d value=%g fixed=%d",
                    name.Data(), ell_value, totalNpars - 1, initialValue, isFixed ? 1 : 0));
            }
        }
    }

    void MassDependentFitter::ParameterManager::AddMassIndependentParameter(
        const std::vector<int>& targetELL,
        const TString filePath) {
        auto parNames = MassDependentFitter::GetParNames();
        auto& massBins = fitter.GetMassBins();
        
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

        std::set<int> targetELLSet(targetELL.begin(), targetELL.end());
        
        double mass_bin_from_file = 0.0;
        tree->SetBranchAddress("mass_bin", &mass_bin_from_file);
        
        std::map<TString, double> paramValues;
        for (const TString& parName : parNames) {
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int ell_value = l_str.Atoi();
            
            if (targetELLSet.find(ell_value) == targetELLSet.end()) {
                continue;
            }
            
            TBranch* branch = tree->GetBranch(parName);
            if (!branch) {
                LogMIWarn("FILE", Form("missing_branch branch=%s", parName.Data()));
                continue;
            }
            
            paramValues[parName] = 0.0;
            tree->SetBranchAddress(parName, &paramValues[parName]);
        }

        Long64_t nEntries = tree->GetEntries();
        for (Long64_t entry = 0; entry < nEntries; entry++) {
            tree->GetEntry(entry);
            
            bool foundMassBin = false;
            for (double targetMassBin : massBins) {
                if (TMath::Abs(mass_bin_from_file - targetMassBin) < 1e-6) {
                    foundMassBin = true;
                    break;
                }
            }
            
            if (!foundMassBin) continue;
            
            for (const auto& [parName, value] : paramValues) {
                TString name = Form("MI_%1.6f_%s", mass_bin_from_file, parName.Data());
                if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                    LogMIInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                    continue;
                }
                
                parsList[name] = value;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                
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

    std::vector<TString> MassDependentFitter::ParameterManager::BuildParamTypesForWave(int funcType, int funcIndex) const {
        std::vector<TString> paramTypes;
        using MT = MassDependentFitter::MassDependenceConfig::ModelType;

        auto getPolyOrder = [&](int idx) -> int {
            if (idx >= 0 && idx < static_cast<int>(fitter.massDepFuncs_.size())) {
                return fitter.massDepFuncs_[idx].GetPolyOrder();
            }
            return -1;
        };

        if (funcType == static_cast<int>(MT::Polynomial)) {
            const int polyOrder = getPolyOrder(funcIndex);
            if (polyOrder >= 0) {
                paramTypes = {"m_threshold", "m_expansion"};
                for (int i = 0; i <= polyOrder; ++i) {
                    paramTypes.push_back(Form("re_%d", i));
                    paramTypes.push_back(Form("im_%d", i));
                }
            }
        } else if (funcType == static_cast<int>(MT::FlattePlusPolynomial)) {
            const int polyOrder = getPolyOrder(funcIndex);
            paramTypes = {"k", "Mass", "g_etapi", "g_KK", "g_etaprimepi"};
            if (polyOrder >= 0) {
                paramTypes.push_back("m_threshold");
                paramTypes.push_back("m_expansion");
                for (int i = 0; i <= polyOrder; ++i) {
                    paramTypes.push_back(Form("re_%d", i));
                    paramTypes.push_back(Form("im_%d", i));
                }
            }
        } else if (funcType == static_cast<int>(MT::BreitWignerPlusPolynomial)) {
            const int polyOrder = getPolyOrder(funcIndex);
            paramTypes = {"k", "Mass", "Width"};
            if (polyOrder >= 0) {
                paramTypes.push_back("m_threshold");
                paramTypes.push_back("m_expansion");
                for (int i = 0; i <= polyOrder; ++i) {
                    paramTypes.push_back(Form("re_%d", i));
                    paramTypes.push_back(Form("im_%d", i));
                }
            }
        } else if (funcType == static_cast<int>(MT::TwoBreitWignerPlusPolynomial)) {
            const int polyOrder = getPolyOrder(funcIndex);
            paramTypes = {"k1", "Mass1", "Width1", "k2", "Mass2", "Width2"};
            if (polyOrder >= 0) {
                paramTypes.push_back("m_threshold");
                paramTypes.push_back("m_expansion");
                for (int i = 0; i <= polyOrder; ++i) {
                    paramTypes.push_back(Form("re_%d", i));
                    paramTypes.push_back(Form("im_%d", i));
                }
            }
        } else if (funcType == static_cast<int>(MT::FlatteBreitWignerPlusPolynomial) || funcType == 4) {
            const int polyOrder = getPolyOrder(funcIndex);
            paramTypes = {"k_fl", "MassFlatte", "g_etapi", "g_KK", "g_etaprimepi", "k_bw", "MassBW", "WidthBW"};
            if (polyOrder >= 0) {
                paramTypes.push_back("m_threshold");
                paramTypes.push_back("m_expansion");
                for (int i = 0; i <= polyOrder; ++i) {
                    paramTypes.push_back(Form("re_%d", i));
                    paramTypes.push_back(Form("im_%d", i));
                }
            }
        } else if (funcType == static_cast<int>(MT::CustomFormula) || funcType == 5) {
            for (const auto& [ell, wConfigs] : fitter.GetMassDependenceConfig().massDependentWaves) {
                for (const auto& wc : wConfigs) {
                    if (fitter.GetFunctionIndexForWave(wc.waveName) == funcIndex && !wc.parameterNames.empty()) {
                        paramTypes = wc.parameterNames;
                        return paramTypes;
                    }
                }
            }
            if (paramTypes.empty()) {
                paramTypes = {"k", "Mass", "Width"};
            }
        }

        return paramTypes;
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
        rng.SetSeed(seed);

        for (const TString& parName : parNames) {
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int ell_value = l_str.Atoi();
            
            if (std::find(targetELL.begin(), targetELL.end(), ell_value) == targetELL.end()) {
                LogMDInfo("SKIP", Form("param=%s reason=targetELL_mismatch ell=%d", parName.Data(), ell_value));
                continue;
            }
            
            if (!fitter.ParameterNeededForMoments(ell_value, hConfig)) {
                LogMDInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), ell_value));
                continue;
            }
            
            if (!config.IsMassDependent(ell_value)) {
                LogMDInfo("SKIP", Form("param=%s reason=l_not_mass_dependent l=%d", parName.Data(), ell_value));
                continue;
            }

            if (!includePhaseParameters && parName.Contains("phi")) {
                LogMDInfo("SKIP", Form("param=%s reason=phase_disabled option=%d", parName.Data(), option));
                continue;
            }

            std::vector<TString> waves = config.GetWaves(ell_value);
            if (waves.empty()) {
                LogMDWarn("SKIP", Form("param=%s reason=no_waves_defined_for_l l=%d", parName.Data(), ell_value));
                continue;
            }

            for (const TString& waveName : waves) {
                int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                if (funcIndex < 0 || funcIndex >= static_cast<int>(fitter.massDepFuncs_.size())) {
                    LogMDWarn("INIT", Form("wave=%s reason=missing_function_index using_defaults=1", waveName.Data()));
                }

                int funcType = -1;
                if (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size())) {
                    funcType = fitter.massDepFuncs_[funcIndex].GetFuncType();
                }
                LogMDInfo("INIT", Form("wave=%s funcType=%d model=%s",
                    waveName.Data(), funcType, ModelTypeName(funcType).Data()));

                std::vector<TString> paramTypes = BuildParamTypesForWave(funcType, funcIndex);

                for (const TString& paramType : paramTypes) {
                    TString name;
                    const bool isSharedResonanceParam = paramType.BeginsWith("Mass") || paramType.BeginsWith("Width") || paramType.BeginsWith("g_");
                    if (isSharedResonanceParam) {
                        name = Form("MD_%s_%s", SharedResonanceKeyForWave(waveName).Data(), paramType.Data());
                    } else {
                        name = Form("MD_%s_%s_%s", parName.Data(), waveName.Data(), paramType.Data());
                    }
                    
                    if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                        LogMDInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                        continue;
                    }

                    double initialValue;
                    if (paramType.BeginsWith("re_") || paramType.BeginsWith("im_")) {
                        TString indexStr = paramType;
                        indexStr.Remove(0, 3);
                        int coeffIndex = indexStr.Atoi();
                        double scale = TMath::Power(0.5, coeffIndex);
                        initialValue = rng.Uniform(-50, 50) * scale;
                        LogMDInfo("INIT", Form("param=%s source=seed field=%s value=%g", name.Data(), paramType.Data(), initialValue));
                    } else if (paramType == "g_etapi") {
                        initialValue = 0.5;
                        LogMDInfo("INIT", Form("param=%s source=seed field=%s value=%g", name.Data(), paramType.Data(), initialValue));
                    } else if (paramType == "g_KK") {
                        initialValue = -99.0;
                        LogMDInfo("INIT", Form("param=%s source=seed field=%s value=%g", name.Data(), paramType.Data(), initialValue));
                    } else if (paramType == "g_etaprimepi") {
                        initialValue = -99.0;
                        LogMDInfo("INIT", Form("param=%s source=seed field=%s value=%g", name.Data(), paramType.Data(), initialValue));
                    } else if (paramType.BeginsWith("k")) {
                        initialValue = rng.Uniform(-30, 30);
                        LogMDInfo("INIT", Form("param=%s source=seed field=%s value=%g", name.Data(), paramType.Data(), initialValue));
                    } else { 
                        initialValue = 0.0;
                        LogMDInfo("INIT", Form("param=%s source=default field=%s value=%g", name.Data(), paramType.Data(), initialValue));
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
                        const int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                        const int funcType = (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size()))
                            ? fitter.massDepFuncs_[funcIndex].GetFuncType()
                            : -1;

                        std::vector<TString> phaseNames;
                        if (funcType == static_cast<int>(MassDependenceConfig::ModelType::FlatteBreitWignerPlusPolynomial) || funcType == 4) {
                            phaseNames = {
                                Form("MD_%s_globalphase_flatte", waveName.Data()),
                                Form("MD_%s_globalphase_bw", waveName.Data())
                            };
                        } else {
                            phaseNames = {Form("MD_%s_globalphase", waveName.Data())};
                        }

                        for (const TString& phaseParName : phaseNames) {
                            if (parsList.find(phaseParName) != parsList.end() || nameToIndex.find(phaseParName) != nameToIndex.end()) {
                                continue;
                            }

                            parsList[phaseParName] = 0.0;
                            nameToIndex[phaseParName] = totalNpars;
                            parIndexNames.push_back(phaseParName);
                            totalNpars++;

                            LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase option=%d state=FREE index=%d value=0",
                                ell_value, phaseParName.Data(), option, totalNpars - 1));
                        }
                        continue;
                    }

                    const int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                    const int funcType = (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size()))
                        ? fitter.massDepFuncs_[funcIndex].GetFuncType()
                        : -1;

                    for (const TString& reflTag : {TString("a"), TString("b")}) {
                        std::vector<TString> phaseNames;
                        if (funcType == static_cast<int>(MassDependenceConfig::ModelType::FlatteBreitWignerPlusPolynomial) || funcType == 4) {
                            phaseNames = {
                                Form("MD_shared%s_%s_globalphase_flatte", reflTag.Data(), waveName.Data()),
                                Form("MD_shared%s_%s_globalphase_bw", reflTag.Data(), waveName.Data())
                            };
                        } else {
                            phaseNames = {Form("MD_shared%s_%s_globalphase", reflTag.Data(), waveName.Data())};
                        }

                        for (const TString& phaseParName : phaseNames) {
                            if (parsList.find(phaseParName) != parsList.end() || nameToIndex.find(phaseParName) != nameToIndex.end()) {
                                continue;
                            }

                            parsList[phaseParName] = 0.0;
                            nameToIndex[phaseParName] = totalNpars;
                            parIndexNames.push_back(phaseParName);
                            totalNpars++;

                            LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase reflectivity=%s option=%d state=FREE index=%d value=0",
                                ell_value, phaseParName.Data(), reflTag.Data(), option, totalNpars - 1));
                        }
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
        
        if (tree->GetEntries() > 0) {
            tree->GetEntry(0);
        }

        for (const TString& parName : parNames) {
            std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
            if (parts->GetEntries() < 3) continue;
            
            TString l_str = ((TObjString*)parts->At(1))->GetString();
            int ell_value = l_str.Atoi();
            
            if (std::find(targetELL.begin(), targetELL.end(), ell_value) == targetELL.end()) {
                LogMDInfo("SKIP", Form("param=%s reason=targetELL_mismatch ell=%d", parName.Data(), ell_value));
                continue;
            }
            
            if (!fitter.ParameterNeededForMoments(ell_value, hConfig)) {
                LogMDInfo("SKIP", Form("param=%s reason=not_needed_for_selected_moments l=%d", parName.Data(), ell_value));
                continue;
            }
            
            if (!config.IsMassDependent(ell_value)) {
                LogMDInfo("SKIP", Form("param=%s reason=l_not_mass_dependent l=%d", parName.Data(), ell_value));
                continue;
            }

            if (!includePhaseParameters && parName.Contains("phi")) {
                LogMDInfo("SKIP", Form("param=%s reason=phase_disabled option=%d", parName.Data(), option));
                continue;
            }

            std::vector<TString> waves = config.GetWaves(ell_value);
            if (waves.empty()) {
                LogMDWarn("SKIP", Form("param=%s reason=no_waves_defined_for_l l=%d", parName.Data(), ell_value));
                continue;
            }

            for (const TString& waveName : waves) {
                int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                int funcType = (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size()))
                    ? fitter.massDepFuncs_[funcIndex].GetFuncType()
                    : -1;
                
                std::vector<TString> paramTypes = BuildParamTypesForWave(funcType, funcIndex);
                
                for (const TString& paramType : paramTypes) {
                    TString name;
                    const bool isSharedResonanceParam = paramType.BeginsWith("Mass") || paramType.BeginsWith("Width") || paramType.BeginsWith("g_");
                    if (isSharedResonanceParam) {
                        name = Form("MD_%s_%s", SharedResonanceKeyForWave(waveName).Data(), paramType.Data());
                    } else {
                        name = Form("MD_%s_%s_%s", parName.Data(), waveName.Data(), paramType.Data());
                    }
                    
                    if (parsList.find(name) != parsList.end() || nameToIndex.find(name) != nameToIndex.end()) {
                        LogMDInfo("SKIP", Form("param=%s reason=already_exists", name.Data()));
                        continue;
                    }

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
                        std::cout << "Parameter " << name << " not found in file. Using default value based on paramType." << std::endl;
                        if (paramType == "k") {
                            value = 1.0;
                            LogMDInfo("INIT", Form("param=%s source=default field=k value=%g", name.Data(), value));
                        } else if (paramType == "g_etapi") {
                            value = 0.5;
                            LogMDInfo("INIT", Form("param=%s source=default field=g_etapi value=%g", name.Data(), value));
                        } else if (paramType == "g_KK") {
                            value = -99.0;
                            LogMDInfo("INIT", Form("param=%s source=default field=g_KK value=%g", name.Data(), value));
                        } else if (paramType == "g_etaprimepi") {
                            value = -99.0;
                            LogMDInfo("INIT", Form("param=%s source=default field=g_etaprimepi value=%g", name.Data(), value));
                        } else if (paramType == "m_threshold") {
                            value = 0.683;
                            LogMDInfo("INIT", Form("param=%s source=default field=m_threshold value=%g", name.Data(), value));
                        } else if (paramType == "m_expansion") {
                            value = 1.2;
                            LogMDInfo("INIT", Form("param=%s source=default field=m_expansion value=%g", name.Data(), value));
                        } else { 
                            value = 0.0;
                            LogMDInfo("INIT", Form("param=%s source=default field=%s value=%g", name.Data(), paramType.Data(), value));
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
                        const int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                        const int funcType = (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size()))
                            ? fitter.massDepFuncs_[funcIndex].GetFuncType()
                            : -1;

                        std::vector<TString> phaseNames;
                        if (funcType == static_cast<int>(MassDependenceConfig::ModelType::FlatteBreitWignerPlusPolynomial) || funcType == 4) {
                            phaseNames = {
                                Form("MD_%s_globalphase_flatte", waveName.Data()),
                                Form("MD_%s_globalphase_bw", waveName.Data())
                            };
                        } else {
                            phaseNames = {Form("MD_%s_globalphase", waveName.Data())};
                        }

                        for (const TString& phaseParName : phaseNames) {
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
                            totalNpars++;

                            LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase option=%d state=FREE index=%d value=%g",
                                ell_value, phaseParName.Data(), option, totalNpars - 1, value));
                        }
                        continue;
                    }

                    const int funcIndex = fitter.GetFunctionIndexForWave(waveName);
                    const int funcType = (funcIndex >= 0 && funcIndex < static_cast<int>(fitter.massDepFuncs_.size()))
                        ? fitter.massDepFuncs_[funcIndex].GetFuncType()
                        : -1;

                    for (const TString& reflTag : {TString("a"), TString("b")}) {
                        std::vector<TString> phaseNames;
                        if (funcType == static_cast<int>(MassDependenceConfig::ModelType::FlatteBreitWignerPlusPolynomial) || funcType == 4) {
                            phaseNames = {
                                Form("MD_shared%s_%s_globalphase_flatte", reflTag.Data(), waveName.Data()),
                                Form("MD_shared%s_%s_globalphase_bw", reflTag.Data(), waveName.Data())
                            };
                        } else {
                            phaseNames = {Form("MD_shared%s_%s_globalphase", reflTag.Data(), waveName.Data())};
                        }

                        for (const TString& phaseParName : phaseNames) {
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
                            totalNpars++;

                            LogMDInfo("ADD", Form("scope=L%d param=%s role=globalphase reflectivity=%s option=%d state=FREE index=%d value=%g",
                                ell_value, phaseParName.Data(), reflTag.Data(), option, totalNpars - 1, value));
                        }
                    }
                }
            }
        }
        
        file->Close();
    }

    std::map<std::string, std::pair<double, double>> MassDependentFitter::ParameterManager::GetAmplitudeValuesAtMassBins(
        const std::vector<int>& targetELL,
        double massBinCenter) {
        
        std::map<std::string, std::pair<double, double>> amplitudes;
        TString massBinStr = Form("%.6f", massBinCenter);
        
        for (int ell_value : targetELL) {
            for (const std::string& refl : {"a", "b"}) {
                for (int m = 0; m <= ell_value; ++m) {
                    std::string key = refl + "_" + std::to_string(ell_value) + "_" + std::to_string(m);
                    
                    TString mag_param_name = Form("MI_%s_%s", massBinStr.Data(), key.c_str());
                    double magnitude = 0.0;
                    
                    auto mag_it = parsList.find(mag_param_name);
                    if (mag_it != parsList.end()) {
                        magnitude = mag_it->second;
                        LogMIInfo("LOOKUP", Form("param=%s found=1 value=%g", mag_param_name.Data(), magnitude));
                    } else {
                        LogMIWarn("LOOKUP", Form("param=%s found=0 default=0", mag_param_name.Data()));
                    }
                    
                    TString phase_param_name = Form("MI_%s_%sphi_%d_%d", massBinStr.Data(), refl.c_str(), ell_value, m);
                    double phase = 0.0;
                    
                    auto phase_it = parsList.find(phase_param_name);
                    if (phase_it != parsList.end()) {
                        phase = phase_it->second;
                        LogMIInfo("LOOKUP", Form("param=%s found=1 value=%g", phase_param_name.Data(), phase));
                    } else {
                        LogMIWarn("LOOKUP", Form("param=%s found=0 default=0", phase_param_name.Data()));
                    }
                    
                    amplitudes[key] = {magnitude, phase};
                    
                    LogMIInfo("AMPL", Form("key=%s mass_bin=%g magnitude=%g phase=%g",
                        key.c_str(), massBinCenter, magnitude, phase));
                }
            }
        }
        
        return amplitudes;
    }

    void MassDependentFitter::MakeResultTree(const ParameterManager& paramManager, const TString& fileName) const {
        if (pars_.empty()) {
            LogFitWarn("RESULT", "reason=no_parameters_to_save");
            return;
        }

        TFile* file = TFile::Open(fileName, "RECREATE");
        if (!file || file->IsZombie()) {
            if (file) delete file;
            LogFitWarn("RESULT", Form("failed_to_create_file path=%s", fileName.Data()));
            return;
        }

        TTree* tree = new TTree("result", "result");
        const int nMassBins = static_cast<int>(massBins_.size());
        const int nPars = pars_.begin()->second.Nvars();

        std::vector<double> par_vals(nPars);
        double mass_bin_center = 0.0;

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

        std::map<TString, double> hMomentBranchValues;
        for (const TString& momentName : allMomentNames) {
            hMomentBranchValues[momentName] = 0.0;
            tree->Branch(momentName, &hMomentBranchValues[momentName]);
        }
        
        LogFitInfo("RESULT", "filling_tree=result content=amplitude_and_moments");

        for (const auto& [mass_bin, pars] : pars_) {
            mass_bin_center = mass_bin;

            for (int i = 0; i < nPars; ++i) {
                const TString par_name = pars.GetParName(i);
                par_vals[i] = pars.GetCurrentVal(par_name);
            }

            for (const TString& momentName : allMomentNames) {
                hMomentBranchValues[momentName] = 0.0;

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

        LogFitInfo("RESULT", "filling_tree=mass_dependent_params");

        TTree* mdTree = new TTree("mass_dependent_params", "Mass Dependent Parameters");
        
        const size_t nMDPars = paramManager.parsList.size();
        std::vector<double> md_par_vals(nMDPars, 0.0);
        std::vector<double> md_par_errors(nMDPars, 0.0);
        size_t idx = 0;
        for (const auto& [parName, value] : paramManager.parsList) {
            if (parName.BeginsWith("MD_")) {
                LogFitInfo("RESULT", Form("tree=mass_dependent_params param=%s value=%g", parName.Data(), value));
                md_par_vals[idx] = value;
                
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

        TTree* allParamsTree = new TTree("all_params", "All Fitted Parameters");
        const size_t nTotalPars = paramManager.parIndexNames.size();
        std::vector<double> all_par_vals(nTotalPars, 0.0);
        std::vector<double> all_par_errors(nTotalPars, 0.0);
        for (size_t i = 0; i < nTotalPars; ++i) {
            const TString& pName = paramManager.parIndexNames[i];
            auto valIt = paramManager.parsList.find(pName);
            all_par_vals[i] = (valIt != paramManager.parsList.end()) ? valIt->second : 0.0;
            double error = 0.0;
            auto it = paramManager.nameToIndex.find(pName);
            if (it != paramManager.nameToIndex.end() && it->second < static_cast<int>(paramManager.parErrors.size())) {
                error = paramManager.parErrors[it->second];
            }
            all_par_errors[i] = error;
            allParamsTree->Branch(pName, &all_par_vals[i]);
            allParamsTree->Branch(pName + "_error", &all_par_errors[i]);
        }
        allParamsTree->Fill();

        TTree* chi2_tree = new TTree("chi2", "chi2");
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
        file->Close();
        delete file;
        LogFitInfo("RESULT", Form("saved path=%s", fileName.Data()));
    }

} // namespace m2pw

#endif // PARAMETER_MANAGER_CPP


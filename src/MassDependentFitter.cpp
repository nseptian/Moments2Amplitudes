#include "MassDependentFitter.h"
#include "FitInternalUtils.h"
#include "TObjString.h"
#include "ROOT/TThreadExecutor.hxx"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <numeric>
#include <thread>

namespace m2pw {


    MassDependentFitter::MassDependentFitter(const Setup& setup, 
                                           std::vector<double> mass_bins, 
                                           double l_max, double noise, 
                                           std::vector<TString> noUse,
                                           const MassDependenceConfig& massDepConfig,
                                           const MomentsConfig& hMomentsConfig) 
        : massBins_(std::move(mass_bins)), massDependenceConfig_(massDepConfig), hMomentsConfig_(hMomentsConfig) {
        
        InitializeMassBins(massBins_);
        InitializeMDFunctions();
        
        const auto validEqn = [&noUse](const TString& name) {
            if (name == "normalise") return false;
            return std::none_of(noUse.begin(), noUse.end(), 
                               [&name](const TString& match) { return name.Contains(match); });
        };

        for (const double& mass_bin : massBins_) {
            pars_[mass_bin] = ParameterHelper(setup);
            
            for (const auto& form : setup.ParameterFormulas()) {
                if (!validEqn(form->GetName())) continue;
                
                if (auto* var = dynamic_cast<RooFormulaVar*>(form)) {
                    if (IsConstraint(var->GetName())) {
                        eqns_[mass_bin].emplace_back(var, &pars_[mass_bin], noise);
                    } else {
                        others_[mass_bin].emplace_back(var, &pars_[mass_bin]);
                    }
                }
            }

            for (auto& eqn : eqns_[mass_bin]) {
                eqn.FindDependencies(others_[mass_bin]);
            }

            for (auto& eqn : eqns_[mass_bin]) {
                eqn.FindL();
                eqn.DoEval(pars_[mass_bin].CurrentVals());
            }
            for (auto& eqn : others_[mass_bin]) {
                eqn.DoEval(pars_[mass_bin].CurrentVals());
            }
        }

        BuildEquationLCache();
        BuildParameterInfoCache();
        BuildMassDependenceLookupCache();
        BuildParameterEvalPlanCache();
        BuildIndexVectors();
    }

    void MassDependentFitter::BuildIndexVectors() {
        massBinToIndex_.clear();
        eqnsVec_.clear();
        parsVec_.clear();
        eqnsVec_.reserve(massBins_.size());
        parsVec_.reserve(massBins_.size());
        for (size_t b = 0; b < massBins_.size(); ++b) {
            const double mass_bin = massBins_[b];
            massBinToIndex_[mass_bin] = b;
            eqnsVec_.push_back(&eqns_[mass_bin]);
            parsVec_.push_back(&pars_[mass_bin]);
        }
    }

    void MassDependentFitter::BuildEquationLCache() {
        equationLCache_.clear();
        equationLCacheVec_.clear();
        equationLCacheVec_.resize(massBins_.size());
        for (size_t b = 0; b < massBins_.size(); ++b) {
            const double mass_bin = massBins_[b];
            auto it = eqns_.find(mass_bin);
            if (it == eqns_.end()) continue;
            const auto& equations = it->second;
            std::vector<int> lValues;
            lValues.reserve(equations.size());
            for (const auto& eqn : equations) {
                lValues.push_back(ExtractLFromEquationName(eqn.GetName()));
            }
            equationLCache_[mass_bin] = lValues;
            equationLCacheVec_[b] = std::move(lValues);
        }
    }

    void MassDependentFitter::BuildParameterInfoCache() {
        parameterInfoCache_.clear();
        for (const auto& [mass_bin, pars] : pars_) {
            std::vector<ParameterInfo> infos;
            infos.reserve(pars.Nvars());

            for (int j = 0; j < pars.Nvars(); ++j) {
                infos.emplace_back(pars.GetParName(j));
            }

            parameterInfoCache_[mass_bin] = std::move(infos);
        }
    }

    void MassDependentFitter::BuildMassDependenceLookupCache() {
        massDependentEllCache_.clear();
        massDependentWaveNamesCache_.clear();

        for (const auto& [ell_value, waveConfigs] : massDependenceConfig_.massDependentWaves) {
            massDependentEllCache_[ell_value] = true;

            std::vector<TString> waveNames;
            waveNames.reserve(waveConfigs.size());
            for (const auto& waveConfig : waveConfigs) {
                waveNames.push_back(waveConfig.waveName);
            }

            massDependentWaveNamesCache_[ell_value] = std::move(waveNames);
        }
    }

    void MassDependentFitter::BuildParameterEvalPlanCache() {
        parameterEvalPlanCache_.clear();

        for (const auto& [mass_bin, pars] : pars_) {
            std::vector<ParameterEvalPlan> plans;
            plans.reserve(pars.Nvars());

            const auto infoCacheIt = parameterInfoCache_.find(mass_bin);
            for (int j = 0; j < pars.Nvars(); ++j) {
                const ParameterInfo* info = nullptr;
                std::unique_ptr<ParameterInfo> fallbackInfo;

                if (infoCacheIt != parameterInfoCache_.end() && j < static_cast<int>(infoCacheIt->second.size())) {
                    info = &infoCacheIt->second[j];
                } else {
                    fallbackInfo = ParseParameterName(pars.GetParName(j));
                    if (!fallbackInfo) {
                        continue;
                    }
                    info = fallbackInfo.get();
                }

                ParameterEvalPlan plan;
                plan.par_name = info->name;
                plan.base_name = info->base_name;
                plan.ell_value = info->ell_value;
                plan.isPhase = info->isPhase;
                plan.isMassDependent = IsMassDependentEll(info->ell_value);
                plan.isSharedResonance = (plan.par_name.BeginsWith("Mass") || plan.par_name.BeginsWith("Width") || plan.par_name.BeginsWith("g_"));
                
                plans.push_back(std::move(plan));
            }

            parameterEvalPlanCache_[mass_bin] = std::move(plans);
        }
    }

    bool MassDependentFitter::IsMassDependentEll(int ell_value) const {
        return massDependentEllCache_.find(ell_value) != massDependentEllCache_.end();
    }

    const std::vector<TString>& MassDependentFitter::GetMassDependentWaveNames(int ell_value) const {
        static const std::vector<TString> empty;
        const auto it = massDependentWaveNamesCache_.find(ell_value);
        if (it == massDependentWaveNamesCache_.end()) {
            return empty;
        }

        return it->second;
    }

    void MassDependentFitter::InitializeMassBins(const std::vector<double>& mass_bins) {
        if (mass_bins.empty()) {
            throw std::invalid_argument("Mass bins cannot be empty");
        }
    }

    void MassDependentFitter::InitializeMDFunctions() {
        if (massBins_.size() < 2) {
            throw std::invalid_argument("At least two mass bins required");
        }
        
        const double bin_width = massBins_[1] - massBins_[0];
        massDepFuncs_.clear();
        
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
            if (!waveConfig.formula.IsNull()) {
                massDepFuncs_.back().SetFormulaString(waveConfig.formula.Data());
            }
            if (!waveConfig.parameterNames.empty()) {
                massDepFuncs_.back().SetNCustomParams(static_cast<int>(waveConfig.parameterNames.size()));
            }
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
        BuildMassDependenceLookupCache();
        BuildParameterEvalPlanCache();
    }

    void MassDependentFitter::SetWaveCustomFormula(int ell, const TString& waveName,
                                                  std::function<std::complex<double>(double, const std::vector<double>&)> formula,
                                                  const TString& formulaStr) {
        int funcIndex = GetFunctionIndexForWave(waveName);
        if (funcIndex >= 0 && funcIndex < static_cast<int>(massDepFuncs_.size())) {
            massDepFuncs_[funcIndex].SetCustomFormula(formula, formulaStr.Data());
            LogFitInfo("FORMULA", Form("custom_formula_assigned ell=%d wave=%s formula=%s",
                ell, waveName.Data(), formulaStr.Data()));
        } else {
            LogFitWarn("FORMULA", Form("custom_formula_failed wave=%s not_found", waveName.Data()));
        }
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

    void MassDependentFitter::SetMinimizerStrategy(int strategy) {
        if (strategy < 0 || strategy > 2) {
            LogFitWarn("MIN", Form("set_strategy_failed value=%d reason=out_of_range_expected_0_to_2", strategy));
            return;
        }

        minimizerStrategy_ = strategy;
        LogFitInfo("MIN", Form("set_strategy=%d", minimizerStrategy_));
    }

    void MassDependentFitter::SetMinimizerTolerance(double tolerance) {
        if (tolerance <= 0.0) {
            LogFitWarn("MIN", Form("set_tolerance_failed value=%g reason=non_positive", tolerance));
            return;
        }

        minimizerTolerance_ = tolerance;
        LogFitInfo("MIN", Form("set_tolerance=%g", minimizerTolerance_));
    }

    void MassDependentFitter::SetMaxIterations(int maxIterations) {
        if (maxIterations <= 0) {
            LogFitWarn("MIN", Form("set_max_iterations_failed value=%d reason=non_positive", maxIterations));
            return;
        }

        maxIterations_ = maxIterations;
        LogFitInfo("MIN", Form("set_max_iterations=%d", maxIterations_));
    }

    void MassDependentFitter::SetNumThreads(int nthreads) {
        if (nthreads <= 0) {
            LogFitWarn("PARALLEL", Form("set_num_threads_non_positive value=%d action=use_auto_detect", nthreads));
            numThreads_ = 0;
        } else {
            numThreads_ = nthreads;
            LogFitInfo("PARALLEL", Form("set_num_threads=%d", numThreads_));
        }
    }

    unsigned int MassDependentFitter::NDim() const {
        if (pars_.empty()) return 0;
        
        unsigned int n = 0;
        for (const auto& [mass_bin, pars] : pars_) {
            n += pars.Nvars();
        }
        return n;
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

    void MassDependentFitter::PrintEquations(const TString opt, const double mass_bin_center) const {
        const auto it = eqns_.find(mass_bin_center);
        if (it != eqns_.end()) {
            for (const auto& eqn : it->second) {
                eqn.Print(opt);
            }
        } else {
            LogFitWarn("EQUATION", Form("mass_bin=%g reason=no_equations_found", mass_bin_center));
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
        minimizer->SetMaxIterations(maxIterations_);
        LogFitInfo("MIN", Form("max_iterations=%d", maxIterations_));
        minimizer->SetStrategy(minimizerStrategy_);
        LogFitInfo("MIN", Form("strategy=%d", minimizerStrategy_));
        minimizer->SetTolerance(minimizerTolerance_);
        LogFitInfo("MIN", Form("tolerance=%g", minimizerTolerance_));

        Chi2Function chi2Function(*this, paramManager.parIndexNames, paramManager);
        ROOT::Math::Functor functor(chi2Function, paramManager.totalNpars);
        minimizer->SetFunction(functor);

        for (int i = 0; i < paramManager.totalNpars; ++i) {
            const TString& parName = paramManager.parIndexNames[i];
            const bool fixed = std::find(paramManager.fixedParNames.begin(),
                                         paramManager.fixedParNames.end(),
                                         parName) != paramManager.fixedParNames.end();
            const double value = initialValues[i];

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

        inMinimization_ = true;
        const bool minimizerConverged = minimizer->Minimize();
        inMinimization_ = false;

        minimizerIsValid_ = minimizerConverged && (minimizer->Status() == 0);
        minimizerHasValidErrors_ = minimizer->IsValidError();
        minimizerStatus_ = minimizer->Status();
        lastChi2_ = minimizer->MinValue();

        minimizer->PrintResults();

        const double* parValues = minimizer->X();

        for (int i = 0; i < paramManager.totalNpars; ++i) {
            const TString& parName = paramManager.parIndexNames[i];
            paramManager.parsList[parName] = parValues[i];
        }

        chi2Function(parValues);
        paramManager.StoreResults(minimizer.get());
    }

    std::complex<double> MassDependentFitter::GetCombinedAmplitude(
        double mass_bin,
        const TString& par_name,
        const std::map<TString, std::vector<double>>& massDepPars,
        int ell_value) const {

        std::complex<double> total_amplitude(0.0, 0.0);
        const std::vector<TString>& waveNames = GetMassDependentWaveNames(ell_value);

        for (const TString& waveName : waveNames) {
            const TString key = Form("%s_%s", MassDependentBaseKey(par_name).Data(), waveName.Data());
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
            bool includeGlobalPhase = false;
            if (funcType == static_cast<int>(MassDependenceConfig::ModelType::FlatteBreitWignerPlusPolynomial) || funcType == 4) {
                double phaseFlatte = 0.0;
                double phaseBW = 0.0;
                bool foundAny = false;

                auto flatteIt = massDepPars.find(GlobalPhaseComponentKeyForWave(waveName, reflTag, "globalphase_flatte"));
                if (flatteIt == massDepPars.end() && reflTag.Length() > 0) {
                    flatteIt = massDepPars.find(GlobalPhaseComponentKeyForWave(waveName, "globalphase_flatte"));
                }
                if (flatteIt != massDepPars.end() && !flatteIt->second.empty()) {
                    phaseFlatte = flatteIt->second[0];
                    foundAny = true;
                }

                auto bwIt = massDepPars.find(GlobalPhaseComponentKeyForWave(waveName, reflTag, "globalphase_bw"));
                if (bwIt == massDepPars.end() && reflTag.Length() > 0) {
                    bwIt = massDepPars.find(GlobalPhaseComponentKeyForWave(waveName, "globalphase_bw"));
                }
                if (bwIt != massDepPars.end() && !bwIt->second.empty()) {
                    phaseBW = bwIt->second[0];
                    foundAny = true;
                }

                includeGlobalPhase = foundAny;
                if (includeGlobalPhase) {
                    params.push_back(phaseFlatte);
                    params.push_back(phaseBW);
                }
            } else {
                auto phaseIt = massDepPars.find(GlobalPhaseKeyForWave(waveName, reflTag));
                if (phaseIt == massDepPars.end() && reflTag.Length() > 0) {
                    phaseIt = massDepPars.find(GlobalPhaseKeyForWave(waveName));
                }

                includeGlobalPhase = (phaseIt != massDepPars.end() && !phaseIt->second.empty());
                if (includeGlobalPhase) {
                    params.push_back(phaseIt->second[0]);
                }
            }
            total_amplitude += massDepFuncs_[funcIndex].GetAmplitude(mass_bin, params, includeGlobalPhase);
        }

        return total_amplitude;
    }

    int MassDependentFitter::GetFunctionIndexForWave(const TString& waveName) const {
        const auto it = waveToFunctionIndex_.find(waveName);
        if (it != waveToFunctionIndex_.end()) {
            return it->second;
        }
        
        LogMDWarn("LOOKUP", Form("wave=%s reason=unknown_wave_name", waveName.Data()));
        return -1;
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
        
        for (int L = 0; L <= 4; ++L) {
            bool shouldInclude = hMomentsConfig_.ShouldIncludeL(L);
            LogFitInfo("MOMENTS", Form("L=%d include=%d", L, shouldInclude));
        }
    }

    void MassDependentFitter::PrintCovarianceMatrix(double massBin, bool printInverse) const {
        auto it = binCovarianceCache_.find(massBin);
        if (it == binCovarianceCache_.end()) {
            for (auto fit = binCovarianceCache_.begin(); fit != binCovarianceCache_.end(); ++fit) {
                if (std::abs(fit->first - massBin) < 1e-4) {
                    it = fit;
                    break;
                }
            }
        }
        if (it == binCovarianceCache_.end()) {
            std::cerr << "Mass bin " << massBin << " not found in covariance cache." << std::endl;
            return;
        }

        const auto& cache = it->second;
        const size_t K = cache.momentNames.size();
        std::cout << "\n=== Fitter Covariance Matrix for Mass Bin: " << it->first << " GeV (valid=" << cache.isValid << ") ===" << std::endl;
        std::cout << "Included moments (" << K << "): [";
        for (size_t i = 0; i < K; ++i) {
            std::cout << cache.momentNames[i] << (i + 1 < K ? ", " : "");
        }
        std::cout << "]" << std::endl;
        if (cache.isValid && K > 0) {
            if (printInverse) {
                std::cout << "Inverse Covariance Matrix (W = V^-1):" << std::endl;
                cache.invCovMatrix.Print();
            } else {
                std::cout << "Covariance Matrix (V):" << std::endl;
                cache.covMatrix.Print();
            }
        } else {
            std::cout << "No valid covariance matrix for this bin." << std::endl;
        }
    }

    void MassDependentFitter::PrintCovarianceMatrices(bool printInverse) const {
        std::cout << "\n========================================================" << std::endl;
        std::cout << " Fitter Covariance Matrices for " << binCovarianceCache_.size() << " Mass Bins" << std::endl;
        std::cout << "========================================================" << std::endl;
        for (const auto& [mass, cache] : binCovarianceCache_) {
            PrintCovarianceMatrix(mass, printInverse);
        }
    }

    int MassDependentFitter::ExtractLFromEquationName(const TString& eqnName) const {
        if (!eqnName.BeginsWith("H_")) {
            return -1;
        }
        
        std::unique_ptr<TObjArray> parts(eqnName.Tokenize("_"));
        if (parts->GetEntries() >= 3) {
            TString L_str = ((TObjString*)parts->At(2))->GetString();
            return L_str.Atoi();
        }
        
        return -1;
    }

    bool MassDependentFitter::ParameterNeededForMoments(int ell, const MomentsConfig& hConfig) const {
        std::set<int> includedMomentL;
        
        for (int H_L = 0; H_L <= 6; ++H_L) {
            if (hConfig.ShouldIncludeL(H_L)) {
                includedMomentL.insert(H_L);
            }
        }
        
        if (includedMomentL.size() == 1 && includedMomentL.count(4) == 1) {
            return (ell == 2);
        }
        
        if (includedMomentL.size() == 2 && 
            includedMomentL.count(2) == 1 && includedMomentL.count(4) == 1) {
            return (ell == 0 || ell == 2);
        }
        
        if (hConfig.includeAll) {
            return (ell <= 2);
        }
        
        return (ell <= 2);
    }

} // namespace m2pw

// Modular subcomponent implementations
#include "ParameterManager.cpp"
#include "Chi2Engine.cpp"
#include "MCMCFitter.cpp"

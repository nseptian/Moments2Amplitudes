#ifndef CHI2_ENGINE_CPP
#define CHI2_ENGINE_CPP

#include "MassDependentFitter.h"
#include "FitInternalUtils.h"
#include "ROOT/TThreadExecutor.hxx"
#include <TError.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <thread>

namespace m2pw {

    double MassDependentFitter::EvaluateChi2ForMassBin(size_t binIdx, ParameterHelper& pars) const {
        if (binIdx >= eqnsVec_.size() || !eqnsVec_[binIdx]) return 0.0;
        const auto current_vals = pars.CurrentVals();
        auto& equations = *eqnsVec_[binIdx];

        if (useCovarianceMatrix_ && binIdx < binCovarianceCacheVec_.size()) {
            const auto& cache = binCovarianceCacheVec_[binIdx];
            if (cache.isValid) {
                const size_t K = cache.eqnIndices.size();
                Double_t* delta = cache.deltaVec.GetMatrixArray();
                for (size_t i = 0; i < K; ++i) {
                    auto& eqn = equations[cache.eqnIndices[i]];
                    eqn.SetNeedsRecalc();
                    delta[i] = eqn.DoEvalResidual(current_vals);
                }
                return cache.invCovMatrix.Similarity(cache.deltaVec);
            }
        }

        double chi2 = 0.0;
        const auto& lValues = (binIdx < equationLCacheVec_.size()) ? equationLCacheVec_[binIdx] : std::vector<int>{};

        for (size_t idx = 0; idx < equations.size(); ++idx) {
            auto& eqn = equations[idx];
            int L_value = (idx < lValues.size()) ? lValues[idx] : ExtractLFromEquationName(eqn.GetName());
            
            bool shouldInclude = true;
            if (L_value >= 0) {
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

    double MassDependentFitter::EvaluateChi2ForMassBin(double massBin, ParameterHelper& pars) const {
        auto it = massBinToIndex_.find(massBin);
        if (it != massBinToIndex_.end()) {
            return EvaluateChi2ForMassBin(it->second, pars);
        }

        const auto current_vals = pars.CurrentVals();
        auto eqnsIt = eqns_.find(massBin);
        if (eqnsIt == eqns_.end()) return 0.0;
        auto& equations = eqnsIt->second;

        if (useCovarianceMatrix_) {
            auto covIt = binCovarianceCache_.find(massBin);
            if (covIt != binCovarianceCache_.end() && covIt->second.isValid) {
                const auto& cache = covIt->second;
                const size_t K = cache.eqnIndices.size();
                Double_t* delta = cache.deltaVec.GetMatrixArray();
                for (size_t i = 0; i < K; ++i) {
                    auto& eqn = equations[cache.eqnIndices[i]];
                    eqn.SetNeedsRecalc();
                    delta[i] = eqn.DoEvalResidual(current_vals);
                }
                return cache.invCovMatrix.Similarity(cache.deltaVec);
            }
        }

        double chi2 = 0.0;
        const auto cacheIt = equationLCache_.find(massBin);

        for (size_t idx = 0; idx < equations.size(); ++idx) {
            auto& eqn = equations[idx];
            int L_value = -1;
            if (cacheIt != equationLCache_.end() && idx < cacheIt->second.size()) {
                L_value = cacheIt->second[idx];
            } else {
                L_value = ExtractLFromEquationName(eqn.GetName());
            }
            
            bool shouldInclude = true;
            if (L_value >= 0) {
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

    double MassDependentFitter::DoEval(
        const std::map<TString, std::vector<double>>& massDepPars,
        const std::map<TString, std::map<TString, std::vector<double>>>& massIndepPars) {
        
        double total_chi2 = 0.0;
        std::map<TString, double> sharedResonanceValues;
        
        if (!massBins_.empty()) {
            const double mass_bin = massBins_[0];
            const TString mass_bin_str = TString::Format("%1.2f", mass_bin);
            auto& pars_first = pars_.at(mass_bin);
            const auto planCacheIt = parameterEvalPlanCache_.find(mass_bin);
            const auto massIndepBinIt = massIndepPars.find(mass_bin_str);
            const auto* massIndepBin = (massIndepBinIt != massIndepPars.end()) ? &massIndepBinIt->second : nullptr;
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

            for (int j = 0; j < pars_first.Nvars(); ++j) {
                const ParameterEvalPlan* param_plan = nullptr;
                ParameterEvalPlan fallbackPlan;
                if (planCacheIt != parameterEvalPlanCache_.end() && j < static_cast<int>(planCacheIt->second.size())) {
                    param_plan = &planCacheIt->second[j];
                } else {
                    auto fallbackInfo = ParseParameterName(pars_first.GetParName(j));
                    if (!fallbackInfo) continue;
                    fallbackPlan.par_name = fallbackInfo->name;
                    fallbackPlan.base_name = fallbackInfo->base_name;
                    fallbackPlan.ell_value = fallbackInfo->ell_value;
                    fallbackPlan.isPhase = fallbackInfo->isPhase;
                    fallbackPlan.isMassDependent = IsMassDependentEll(fallbackInfo->ell_value);
                    fallbackPlan.isSharedResonance = (fallbackPlan.par_name.BeginsWith("Mass") || fallbackPlan.par_name.BeginsWith("Width") || fallbackPlan.par_name.BeginsWith("g_"));
                    param_plan = &fallbackPlan;
                }

                if (!param_plan->isSharedResonance) continue;

                const TString& par_name = param_plan->par_name;
                const int ell_value = param_plan->ell_value;
                const bool isMassDependent = param_plan->isMassDependent;
                
                double value = 0.0;
                if (!param_plan->isPhase) {
                    if (isMassDependent) {
                        value = std::abs(getCachedAmplitude(par_name, ell_value));
                    } else {
                        if (massIndepBin) {
                            const auto par_it = massIndepBin->find(par_name);
                            if (par_it != massIndepBin->end() && !par_it->second.empty()) {
                                value = par_it->second[0];
                            }
                        }
                    }
                } else {
                    const TString& base_name = param_plan->base_name;
                    bool hasPhaseInMassIndepPars = false;
                    if (massIndepBin) {
                        const auto par_it = massIndepBin->find(par_name);
                        if (par_it != massIndepBin->end() && !par_it->second.empty()) {
                            hasPhaseInMassIndepPars = true;
                            value = par_it->second[0];
                        }
                    }
                    if (!hasPhaseInMassIndepPars && isMassDependent) {
                        value = std::arg(getCachedAmplitude(base_name, ell_value));
                    }
                }
                
                sharedResonanceValues[par_name] = value;
            }
        }
        
        auto evaluateBin = [&](size_t i) -> double {
            const double mass_bin = massBins_[i];
            const TString mass_bin_str = TString::Format("%1.2f", mass_bin);
            auto& pars = pars_.at(mass_bin);
            const auto planCacheIt = parameterEvalPlanCache_.find(mass_bin);
            const auto massIndepBinIt = massIndepPars.find(mass_bin_str);
            const auto* massIndepBin = (massIndepBinIt != massIndepPars.end()) ? &massIndepBinIt->second : nullptr;
            
            std::map<TString, std::complex<double>> amplitudeCache;

            auto getCachedAmplitude_local = [&](const TString& base_name, int ell_value) -> const std::complex<double>& {
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
            
            for (int j = 0; j < pars.Nvars(); ++j) {
                const ParameterEvalPlan* param_plan = nullptr;
                ParameterEvalPlan fallbackPlan;
                if (planCacheIt != parameterEvalPlanCache_.end() && j < static_cast<int>(planCacheIt->second.size())) {
                    param_plan = &planCacheIt->second[j];
                } else {
                    auto fallbackInfo = ParseParameterName(pars.GetParName(j));
                    if (!fallbackInfo) continue;

                    fallbackPlan.par_name = fallbackInfo->name;
                    fallbackPlan.base_name = fallbackInfo->base_name;
                    fallbackPlan.ell_value = fallbackInfo->ell_value;
                    fallbackPlan.isPhase = fallbackInfo->isPhase;
                    fallbackPlan.isMassDependent = IsMassDependentEll(fallbackInfo->ell_value);
                    fallbackPlan.isSharedResonance = (fallbackPlan.par_name.BeginsWith("Mass") || fallbackPlan.par_name.BeginsWith("Width") || fallbackPlan.par_name.BeginsWith("g_"));
                    param_plan = &fallbackPlan;
                }

                const TString& par_name = param_plan->par_name;
                const int ell_value = param_plan->ell_value;
                const bool isMassDependent = param_plan->isMassDependent;
                const bool isSharedResonance = param_plan->isSharedResonance;
                
                double value = 0.0;
                
                if (isSharedResonance) {
                    const auto it = sharedResonanceValues.find(par_name);
                    if (it != sharedResonanceValues.end()) {
                        value = it->second;
                    }
                } else if (!param_plan->isPhase) {
                    if (isMassDependent) {
                        value = std::abs(getCachedAmplitude_local(par_name, ell_value));
                    } else {
                        if (massIndepBin) {
                            const auto par_it = massIndepBin->find(par_name);
                            if (par_it != massIndepBin->end() && !par_it->second.empty()) {
                                value = par_it->second[0];
                            }
                        }
                    }
                } else {
                    const TString& base_name = param_plan->base_name;
                    bool hasPhaseInMassIndepPars = false;
                    if (massIndepBin) {
                        const auto par_it = massIndepBin->find(par_name);
                        if (par_it != massIndepBin->end() && !par_it->second.empty()) {
                            hasPhaseInMassIndepPars = true;
                            value = par_it->second[0];
                        }
                    }
                    if (!hasPhaseInMassIndepPars && isMassDependent) {
                        value = std::arg(getCachedAmplitude_local(base_name, ell_value));
                    }
                }
                
                pars.SetCurrentVal(par_name, value);
            }
            
            return EvaluateChi2ForMassBin(mass_bin, pars);
        };
        
        const size_t nBins = massBins_.size();
        const bool useParallel = enableParallelization_ && nBins >= 8 && numThreads_ != 1;

        if (useParallel) {
            const unsigned int hw = std::thread::hardware_concurrency();
            size_t requestedWorkers = (numThreads_ > 0) ? static_cast<size_t>(numThreads_) : ((hw > 0) ? static_cast<size_t>(hw) : 2);
            const size_t nWorkers = std::min(nBins, requestedWorkers);
            if (nWorkers > 1) {
                ROOT::TThreadExecutor executor(static_cast<UInt_t>(nWorkers));
                std::vector<size_t> workerIds(nWorkers);
                std::iota(workerIds.begin(), workerIds.end(), 0);

                auto evaluateWorker = [&](size_t workerId) -> double {
                    double localChi2 = 0.0;
                    for (size_t binIdx = workerId; binIdx < nBins; binIdx += nWorkers) {
                        localChi2 += evaluateBin(binIdx);
                    }
                    return localChi2;
                };

                const std::vector<double> partialChi2 = executor.Map(evaluateWorker, workerIds);
                total_chi2 = std::accumulate(partialChi2.begin(), partialChi2.end(), 0.0);
            } else {
                for (size_t i = 0; i < nBins; ++i) {
                    total_chi2 += evaluateBin(i);
                }
            }
        } else {
            for (size_t i = 0; i < nBins; ++i) {
                total_chi2 += evaluateBin(i);
            }
        }
        
        lastChi2_ = total_chi2;
        return total_chi2;
    }

    void MassDependentFitter::SetEquationValues(MassDependentMoments& moms) {
        binCovarianceCache_.clear();
        binCovarianceCacheVec_.clear();
        binCovarianceCacheVec_.resize(massBins_.size());

        for (size_t b = 0; b < massBins_.size(); ++b) {
            const double mass_bin = massBins_[b];
            auto eqnsIt = eqns_.find(mass_bin);
            if (eqnsIt == eqns_.end()) continue;
            auto& equations = eqnsIt->second;

            MomentHelper moments = moms.GetMoments(mass_bin);
            for (auto& eqn : equations) {
                eqn.FindL();
                eqn.SetEquationValue(moments.GetUnnormalizedVal(eqn.GetName()),
                                   moments.GetUnnormalizedError(eqn.GetName()));
                eqn.DoEval(pars_[mass_bin].CurrentVals());
            }

            BinCovarianceCache cache;
            const auto cacheIt = equationLCache_.find(mass_bin);
            for (size_t idx = 0; idx < equations.size(); ++idx) {
                auto& eqn = equations[idx];
                int L_value = -1;
                if (cacheIt != equationLCache_.end() && idx < cacheIt->second.size()) {
                    L_value = cacheIt->second[idx];
                } else {
                    L_value = ExtractLFromEquationName(eqn.GetName());
                }
                if (L_value < 0 || hMomentsConfig_.ShouldIncludeL(L_value)) {
                    // Ensure moment exists in dataset covariance matrix or has valid measurements
                    if (useCovarianceMatrix_ && moments.HasCovarianceMatrix()) {
                        if (moments.GetCovarianceIndex(eqn.GetName()) < 0) {
                            continue;
                        }
                    } else {
                        if (moments.GetUnnormalizedError(eqn.GetName()) <= 0.0 &&
                            moments.GetUnnormalizedVal(eqn.GetName()) == 0.0) {
                            continue;
                        }
                    }
                    cache.eqnIndices.push_back(idx);
                    cache.momentNames.push_back(eqn.GetName());
                }
            }

            const size_t K = cache.momentNames.size();
            if (K > 0) {
                cache.covMatrix.ResizeTo(K, K);
                cache.invCovMatrix.ResizeTo(K, K);
                cache.deltaVec.ResizeTo(K);

                if (!useCovarianceMatrix_) {
                    cache.covMatrix.Zero();
                    cache.invCovMatrix.Zero();
                    for (size_t i = 0; i < K; ++i) {
                        double err = moments.GetUnnormalizedError(cache.momentNames[i]);
                        double var = (err > 0.0) ? err * err : 1.0;
                        cache.covMatrix(i, i) = var;
                        cache.invCovMatrix(i, i) = 1.0 / var;
                    }
                    cache.isValid = true;
                } else {
                    TMatrixDSym V(K);
                    for (size_t i = 0; i < K; ++i) {
                        for (size_t j = i; j < K; ++j) {
                            double cov_val = moments.GetCovariance(cache.momentNames[i], cache.momentNames[j]);
                            V(i, j) = cov_val;
                            V(j, i) = cov_val;
                        }
                    }
                    cache.covMatrix = V;

                    // Temporarily silence internal ROOT TDecompLU errors during the invertibility check
                    int prevIgnore = gErrorIgnoreLevel;
                    gErrorIgnoreLevel = kBreak;
                    TMatrixDSym invV(V);
                    Double_t det = 0.0;
                    invV.Invert(&det);
                    gErrorIgnoreLevel = prevIgnore;

                    if (invV.IsValid() && det != 0.0 && !std::isnan(det) && !std::isinf(det)) {
                        cache.invCovMatrix = invV;
                        cache.isValid = true;
                    } else {
                        LogFitWarn("CHI2", Form("mass_bin=%g covariance_matrix_singular (det=%g) fallback_to_diagonal", mass_bin, det));
                        cache.invCovMatrix.Zero();
                        for (size_t i = 0; i < K; ++i) {
                            double var = V(i, i);
                            if (var <= 0.0) {
                                double err = moments.GetUnnormalizedError(cache.momentNames[i]);
                                var = (err > 0.0) ? err * err : 1.0;
                            }
                            cache.invCovMatrix(i, i) = 1.0 / var;
                        }
                        cache.isValid = true;
                    }
                }
            }
            binCovarianceCacheVec_[b] = cache;
            binCovarianceCache_[mass_bin] = std::move(cache);
        }
    }

    MassDependentFitter::Chi2Function::Chi2Function(
        MassDependentFitter& fitter, 
        const std::vector<TString>& parIndexNames,
        const ParameterManager& paramManager)
        : fitter_(fitter), parIndexNames_(parIndexNames), paramManager_(&paramManager) {
        BuildPrecompiledPlan();
    }

    void MassDependentFitter::Chi2Function::BuildPrecompiledPlan() {
        fastConstraints_.clear();
        amplitudes_.clear();
        binRouting_.clear();
        preassembledWaveParams_.clear();

        if (!paramManager_) {
            precompiled_ = false;
            return;
        }

        for (const auto& kv : paramManager_->gaussianConstraints) {
            const TString& pname = kv.first;
            auto it = paramManager_->nameToIndex.find(pname);
            if (it != paramManager_->nameToIndex.end() && kv.second.second > 0.0) {
                FastConstraint fc;
                fc.minuitIdx = it->second;
                fc.mean = kv.second.first;
                fc.invSigmaSq = 1.0 / (kv.second.second * kv.second.second);
                fastConstraints_.push_back(fc);
            }
        }

        std::map<TString, std::map<TString, std::vector<int>>> massIndepParIndices;
        std::map<TString, std::vector<int>> massDepParIndices;

        for (size_t i = 0; i < parIndexNames_.size(); ++i) {
            const TString& parName = parIndexNames_[i];
            if (parName.BeginsWith("MI_")) {
                int firstUnderscore = parName.Index("_", 3);
                int secondUnderscore = parName.Index("_", firstUnderscore + 1);
                int thirdUnderscore = parName.Index("_", secondUnderscore + 1);
                if (firstUnderscore > 0 && secondUnderscore > 0 && thirdUnderscore > 0) {
                    TString parKey = parName(firstUnderscore + 1, thirdUnderscore + 2 - firstUnderscore);
                    TString massBinStr = parName(3, firstUnderscore - 7);
                    massIndepParIndices[massBinStr][parKey].push_back(static_cast<int>(i));
                }
            } else if (parName.BeginsWith("MD_")) {
                TString parKey = MassDependentBaseKey(parName);
                massDepParIndices[parKey].push_back(static_cast<int>(i));
            }
        }

        std::map<std::pair<TString, int>, int> ampNameToIdx;

        for (size_t b = 0; b < fitter_.massBins_.size(); ++b) {
            const double mass_bin = fitter_.massBins_[b];
            const auto planCacheIt = fitter_.parameterEvalPlanCache_.find(mass_bin);
            if (planCacheIt == fitter_.parameterEvalPlanCache_.end()) continue;

            for (const auto& plan : planCacheIt->second) {
                if (!plan.isMassDependent) continue;

                const TString& baseName = plan.base_name;
                const int ellValue = plan.ell_value;
                auto key = std::make_pair(baseName, ellValue);

                if (ampNameToIdx.find(key) == ampNameToIdx.end()) {
                    int newAmpIdx = static_cast<int>(amplitudes_.size());
                    ampNameToIdx[key] = newAmpIdx;

                    PrecompiledAmplitude pAmp;
                    pAmp.baseName = baseName;
                    pAmp.ellValue = ellValue;

                    const std::vector<TString>& waveNames = fitter_.GetMassDependentWaveNames(ellValue);
                    for (const TString& waveName : waveNames) {
                        int funcIndex = fitter_.GetFunctionIndexForWave(waveName);
                        if (funcIndex < 0 || funcIndex >= static_cast<int>(fitter_.massDepFuncs_.size())) {
                            continue;
                        }

                        int funcType = fitter_.massDepFuncs_[funcIndex].GetFuncType();
                        std::vector<int> minuitIndices;
                        if (!BuildWaveParametersGeneric(baseName, waveName, funcType, massDepParIndices, minuitIndices)) {
                            continue;
                        }

                        const TString reflTag = ReflectivityTagFromParameterName(baseName);
                        bool includeGlobalPhase = false;

                        if (funcType == static_cast<int>(MassDependenceConfig::ModelType::FlatteBreitWignerPlusPolynomial) || funcType == 4) {
                            int phaseFlatteIdx = -1;
                            int phaseBWIdx = -1;
                            bool foundAny = false;

                            auto flatteIt = massDepParIndices.find(GlobalPhaseComponentKeyForWave(waveName, reflTag, "globalphase_flatte"));
                            if (flatteIt == massDepParIndices.end() && reflTag.Length() > 0) {
                                flatteIt = massDepParIndices.find(GlobalPhaseComponentKeyForWave(waveName, "globalphase_flatte"));
                            }
                            if (flatteIt != massDepParIndices.end() && !flatteIt->second.empty()) {
                                phaseFlatteIdx = flatteIt->second[0];
                                foundAny = true;
                            }

                            auto bwIt = massDepParIndices.find(GlobalPhaseComponentKeyForWave(waveName, reflTag, "globalphase_bw"));
                            if (bwIt == massDepParIndices.end() && reflTag.Length() > 0) {
                                bwIt = massDepParIndices.find(GlobalPhaseComponentKeyForWave(waveName, "globalphase_bw"));
                            }
                            if (bwIt != massDepParIndices.end() && !bwIt->second.empty()) {
                                phaseBWIdx = bwIt->second[0];
                                foundAny = true;
                            }

                            includeGlobalPhase = foundAny;
                            if (includeGlobalPhase) {
                                minuitIndices.push_back(phaseFlatteIdx);
                                minuitIndices.push_back(phaseBWIdx);
                            }
                        } else {
                            auto phaseIt = massDepParIndices.find(GlobalPhaseKeyForWave(waveName, reflTag));
                            if (phaseIt == massDepParIndices.end() && reflTag.Length() > 0) {
                                phaseIt = massDepParIndices.find(GlobalPhaseKeyForWave(waveName));
                            }

                            includeGlobalPhase = (phaseIt != massDepParIndices.end() && !phaseIt->second.empty());
                            if (includeGlobalPhase) {
                                minuitIndices.push_back(phaseIt->second[0]);
                            }
                        }

                        PrecompiledWave pWave;
                        pWave.funcIndex = funcIndex;
                        pWave.funcType = funcType;
                        pWave.includeGlobalPhase = includeGlobalPhase;
                        pWave.minuitIndices = std::move(minuitIndices);
                        pAmp.waves.push_back(std::move(pWave));
                    }

                    amplitudes_.push_back(std::move(pAmp));
                }
            }
        }

        binRouting_.clear();
        binRouting_.reserve(fitter_.massBins_.size());

        for (size_t b = 0; b < fitter_.massBins_.size(); ++b) {
            const double mass_bin = fitter_.massBins_[b];
            const TString mass_bin_str = TString::Format("%1.2f", mass_bin);
            const auto& pars = *fitter_.parsVec_[b];
            const auto planCacheIt = fitter_.parameterEvalPlanCache_.find(mass_bin);

            std::vector<BinParamRouting> binRoutings;
            binRoutings.reserve(pars.Nvars());

            for (int j = 0; j < pars.Nvars(); ++j) {
                const ParameterEvalPlan* plan = nullptr;
                if (planCacheIt != fitter_.parameterEvalPlanCache_.end() && j < static_cast<int>(planCacheIt->second.size())) {
                    plan = &planCacheIt->second[j];
                }
                if (!plan) continue;

                BinParamRouting routing;
                routing.varIdx = j;

                const TString& par_name = plan->par_name;
                const TString& base_name = plan->base_name;
                const int ell_value = plan->ell_value;
                const bool isMassDependent = plan->isMassDependent;

                bool foundMI = false;
                const auto binMIIt = massIndepParIndices.find(mass_bin_str);
                if (binMIIt != massIndepParIndices.end()) {
                    const auto parIt = binMIIt->second.find(par_name);
                    if (parIt != binMIIt->second.end() && !parIt->second.empty()) {
                        routing.source = ParamRoutingSource::DirectMinuit;
                        routing.minuitIdx = parIt->second[0];
                        foundMI = true;
                    }
                }

                if (!foundMI) {
                    if (isMassDependent) {
                        auto itAmp = ampNameToIdx.find({base_name, ell_value});
                        if (itAmp != ampNameToIdx.end()) {
                            routing.ampIdx = itAmp->second;
                            if (!plan->isPhase) {
                                routing.source = ParamRoutingSource::AmplitudeMagnitude;
                            } else {
                                routing.source = ParamRoutingSource::AmplitudePhase;
                            }
                        } else {
                            routing.source = ParamRoutingSource::ConstantValue;
                            routing.constVal = 0.0;
                        }
                    } else {
                        routing.source = ParamRoutingSource::ConstantValue;
                        routing.constVal = 0.0;
                    }
                }

                binRoutings.push_back(routing);
            }

            binRouting_.push_back(std::move(binRoutings));
        }

        preassembledWaveParams_.resize(amplitudes_.size());
        for (size_t a = 0; a < amplitudes_.size(); ++a) {
            preassembledWaveParams_[a].resize(amplitudes_[a].waves.size());
            for (size_t w = 0; w < amplitudes_[a].waves.size(); ++w) {
                preassembledWaveParams_[a][w].resize(amplitudes_[a].waves[w].minuitIndices.size(), 0.0);
            }
        }

        precompiled_ = true;
        LogFitInfo("OPT", Form("precompiled_plan_ready nBins=%zu nAmps=%zu nConstraints=%zu",
                               fitter_.massBins_.size(), amplitudes_.size(), fastConstraints_.size()));
    }

    double MassDependentFitter::Chi2Function::EvalFast(const double* mass_dep_pars) const {
        double penalty = 0.0;
        for (const auto& gc : fastConstraints_) {
            const double diff = mass_dep_pars[gc.minuitIdx] - gc.mean;
            penalty += diff * diff * gc.invSigmaSq;
        }

        for (size_t a = 0; a < amplitudes_.size(); ++a) {
            const auto& ampMeta = amplitudes_[a];
            for (size_t w = 0; w < ampMeta.waves.size(); ++w) {
                const auto& wave = ampMeta.waves[w];
                auto& pVec = preassembledWaveParams_[a][w];
                for (size_t p = 0; p < wave.minuitIndices.size(); ++p) {
                    const int idx = wave.minuitIndices[p];
                    pVec[p] = (idx >= 0) ? mass_dep_pars[idx] : 0.0;
                }
            }
        }

        const size_t nBins = fitter_.massBins_.size();
        double total_chi2 = 0.0;

        auto evaluateBinFast = [&](size_t i) -> double {
            const double mass_bin = fitter_.massBins_[i];
            ParameterHelper& pars = *fitter_.parsVec_[i];

            std::vector<std::complex<double>> binAmps(amplitudes_.size(), std::complex<double>(0.0, 0.0));
            for (size_t a = 0; a < amplitudes_.size(); ++a) {
                const auto& ampMeta = amplitudes_[a];
                std::complex<double> totalAmp(0.0, 0.0);
                for (size_t w = 0; w < ampMeta.waves.size(); ++w) {
                    const auto& wave = ampMeta.waves[w];
                    totalAmp += fitter_.massDepFuncs_[wave.funcIndex].GetAmplitude(
                        mass_bin, preassembledWaveParams_[a][w], wave.includeGlobalPhase);
                }
                binAmps[a] = totalAmp;
            }

            const auto& routingList = binRouting_[i];
            for (const auto& r : routingList) {
                double val = 0.0;
                switch (r.source) {
                    case ParamRoutingSource::DirectMinuit:
                        val = mass_dep_pars[r.minuitIdx];
                        break;
                    case ParamRoutingSource::AmplitudeMagnitude:
                        val = std::abs(binAmps[r.ampIdx]);
                        break;
                    case ParamRoutingSource::AmplitudePhase:
                        val = std::arg(binAmps[r.ampIdx]);
                        break;
                    case ParamRoutingSource::ConstantValue:
                        val = r.constVal;
                        break;
                }
                pars.SetCurrentVal(r.varIdx, val);
            }

            return fitter_.EvaluateChi2ForMassBin(i, pars);
        };

        const bool useParallel = fitter_.enableParallelization_ && nBins >= 8 && fitter_.numThreads_ != 1;
        if (useParallel) {
            const unsigned int hw = std::thread::hardware_concurrency();
            size_t requestedWorkers = (fitter_.numThreads_ > 0) ? static_cast<size_t>(fitter_.numThreads_) : ((hw > 0) ? static_cast<size_t>(hw) : 2);
            const size_t nWorkers = std::min(nBins, requestedWorkers);
            if (nWorkers > 1) {
                ROOT::TThreadExecutor executor(static_cast<UInt_t>(nWorkers));
                std::vector<size_t> workerIds(nWorkers);
                std::iota(workerIds.begin(), workerIds.end(), 0);

                auto evaluateWorker = [&](size_t workerId) -> double {
                    double localChi2 = 0.0;
                    for (size_t binIdx = workerId; binIdx < nBins; binIdx += nWorkers) {
                        localChi2 += evaluateBinFast(binIdx);
                    }
                    return localChi2;
                };

                const std::vector<double> partialChi2 = executor.Map(evaluateWorker, workerIds);
                total_chi2 = std::accumulate(partialChi2.begin(), partialChi2.end(), 0.0);
            } else {
                for (size_t i = 0; i < nBins; ++i) {
                    total_chi2 += evaluateBinFast(i);
                }
            }
        } else {
            for (size_t i = 0; i < nBins; ++i) {
                total_chi2 += evaluateBinFast(i);
            }
        }

        fitter_.lastChi2_ = total_chi2;
        return total_chi2 + penalty;
    }

    double MassDependentFitter::Chi2Function::operator()(const double* mass_dep_pars) const {
        if (precompiled_) {
            return EvalFast(mass_dep_pars);
        }

        MapParameters(mass_dep_pars);
        double chi2 = fitter_.DoEval(massDepPars_, massIndepPars_);

        double penalty = 0.0;
        if (paramManager_) {
            for (const auto& kv : paramManager_->gaussianConstraints) {
                const TString& pname = kv.first;
                const auto it = currentValues_.find(pname);
                if (it == currentValues_.end()) continue;
                const double mean = kv.second.first;
                const double sigma = kv.second.second;
                if (sigma <= 0) continue;
                const double diff = it->second - mean;
                penalty += (diff * diff) / (sigma * sigma);
            }
        }

        return chi2 + penalty;
    }

    void MassDependentFitter::Chi2Function::MapParameters(const double* mass_dep_pars) const {
        massDepPars_.clear();
        massIndepPars_.clear();
        currentValues_.clear();
        
        for (size_t i = 0; i < parIndexNames_.size(); ++i) {
            const TString& parName = parIndexNames_[i];
            currentValues_[parName] = mass_dep_pars[i];
            
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

    double MassDependentFitter::EvalChi2() const {
        double total_chi2 = 0.0;
        
        for (const auto& [mass_bin, pars] : pars_) {
            double chi2_for_bin = EvaluateChi2ForMassBin(mass_bin, const_cast<ParameterHelper&>(pars));
            total_chi2 += chi2_for_bin;
            LogFitInfo("CHI2", Form("mass_bin=%g chi2=%g", mass_bin, chi2_for_bin));
        }
        
        lastChi2_ = total_chi2;
        LogFitInfo("CHI2", Form("total=%g", total_chi2));
        return total_chi2;
    }

} // namespace m2pw

#endif // CHI2_ENGINE_CPP


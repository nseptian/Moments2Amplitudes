#ifndef FIT_INTERNAL_UTILS_H
#define FIT_INTERNAL_UTILS_H

#include "MassDependentFitter.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>

namespace m2pw {

    using MassDependenceConfig = MassDependentFitter::MassDependenceConfig;

    namespace internal {

        inline TString ModelTypeName(int funcType) {
            switch (funcType) {
                case 0:
                    return "polynomial";
                case 1:
                    return "flatte_plus_polynomial";
                case 2:
                    return "breitwigner_plus_polynomial";
                case 3:
                    return "two_breitwigner_plus_polynomial";
                case 4:
                    return "flatte_breitwigner_plus_polynomial";
                case 5:
                    return "custom_formula";
                default:
                    return Form("unknown_%d", funcType);
            }
        }

        inline void LogMDInfo(const TString& stage, const TString& message) {
            if (stage == "SKIP") return;
            std::cout << "[MD][" << stage << "] " << message << std::endl;
        }

        inline void LogMDWarn(const TString& stage, const TString& message) {
            std::cerr << "[MD][" << stage << "] " << message << std::endl;
        }

        inline void LogMIInfo(const TString& stage, const TString& message) {
            if (stage == "SKIP") return;
            std::cout << "[MI][" << stage << "] " << message << std::endl;
        }

        inline void LogMIWarn(const TString& stage, const TString& message) {
            std::cerr << "[MI][" << stage << "] " << message << std::endl;
        }

        inline void LogFitInfo(const TString& stage, const TString& message) {
            if (stage == "SKIP") return;
            std::cout << "[FIT][" << stage << "] " << message << std::endl;
        }

        inline void LogFitWarn(const TString& stage, const TString& message) {
            std::cerr << "[FIT][" << stage << "] " << message << std::endl;
        }

        inline TString GlobalPhaseKeyForWave(const TString& waveName) {
            return waveName;
        }

        inline TString GlobalPhaseKeyForWave(const TString& waveName, const TString& reflTag) {
            if (reflTag.Length() == 0) {
                return GlobalPhaseKeyForWave(waveName);
            }
            return Form("shared%s_%s", reflTag.Data(), waveName.Data());
        }

        inline TString GlobalPhaseComponentKeyForWave(const TString& waveName, const TString& component) {
            return Form("%s_%s", waveName.Data(), component.Data());
        }

        inline TString GlobalPhaseComponentKeyForWave(const TString& waveName, const TString& reflTag, const TString& component) {
            if (reflTag.Length() == 0) {
                return GlobalPhaseComponentKeyForWave(waveName, component);
            }
            return Form("shared%s_%s_%s", reflTag.Data(), waveName.Data(), component.Data());
        }

        inline TString ReflectivityTagFromParameterName(const TString& parName) {
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

        inline TString SharedResonanceKeyForWave(const TString& waveName) {
            return Form("shared_%s", waveName.Data());
        }

        inline TString MassDependentBaseKey(const TString& parName) {
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
                "_globalphase_flatte",
                "_globalphase_bw",
                "_globalphase",
                "_m_threshold",
                "_m_expansion",
                "_g_etapi",
                "_g_KK",
                "_g_etaprimepi",
                "_MassFlatte",
                "_MassBW",
                "_WidthBW",
                "_Mass1",
                "_Width1",
                "_Mass2",
                "_Width2",
                "_Mass",
                "_Width",
                "_k_fl",
                "_k_bw",
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

        template<typename T>
        inline bool BuildWaveParametersGeneric(
            const TString& baseName,
            const TString& waveName,
            int funcType,
            const std::map<TString, std::vector<T>>& massDepPars,
            std::vector<T>& params) {

            TString key = Form("%s_%s", baseName.Data(), waveName.Data());
            const auto it = massDepPars.find(key);
            if (it == massDepPars.end() || it->second.empty()) {
                return false;
            }

            params = it->second;

            if (funcType == static_cast<int>(MassDependenceConfig::ModelType::FlattePlusPolynomial)) {
                if (params.empty()) {
                    return false;
                }

                auto sharedIt = massDepPars.find(SharedResonanceKeyForWave(waveName));
                if (sharedIt == massDepPars.end() || sharedIt->second.size() < 4) {
                    return false;
                }
                
                std::vector<T> assembled;
                assembled.reserve(params.size() + 4);
                assembled.push_back(params[0]);         // k
                assembled.push_back(sharedIt->second[0]); // Mass
                assembled.push_back(sharedIt->second[1]); // g_etapi
                assembled.push_back(sharedIt->second[2]); // g_KK
                assembled.push_back(sharedIt->second[3]); // g_etaprimepi

                assembled.insert(assembled.end(), params.begin() + 1, params.end());
                params.swap(assembled);
                return true;
            }

            if (funcType == static_cast<int>(MassDependenceConfig::ModelType::BreitWignerPlusPolynomial)) {
                if (params.empty()) {
                    return false;
                }

                auto sharedIt = massDepPars.find(SharedResonanceKeyForWave(waveName));
                if (sharedIt == massDepPars.end() || sharedIt->second.size() < 2) {
                    return false;
                }

                std::vector<T> assembled;
                assembled.reserve(params.size() + 2);
                assembled.push_back(params[0]);         // k
                assembled.push_back(sharedIt->second[0]); // Mass
                assembled.push_back(sharedIt->second[1]); // Width
                assembled.insert(assembled.end(), params.begin() + 1, params.end());
                params.swap(assembled);
                return true;
            }

            if (funcType == static_cast<int>(MassDependenceConfig::ModelType::TwoBreitWignerPlusPolynomial)) {
                if (params.size() < 2) {
                    return false;
                }

                auto sharedIt = massDepPars.find(SharedResonanceKeyForWave(waveName));
                if (sharedIt == massDepPars.end() || sharedIt->second.size() < 4) {
                    return false;
                }

                std::vector<T> assembled;
                assembled.reserve(params.size() + 4);
                assembled.push_back(params[0]);           // k1
                assembled.push_back(sharedIt->second[0]); // Mass1
                assembled.push_back(sharedIt->second[1]); // Width1
                assembled.push_back(params[1]);           // k2
                assembled.push_back(sharedIt->second[2]); // Mass2
                assembled.push_back(sharedIt->second[3]); // Width2
                assembled.insert(assembled.end(), params.begin() + 2, params.end());
                params.swap(assembled);
                return true;
            }

            if (funcType == static_cast<int>(MassDependenceConfig::ModelType::FlatteBreitWignerPlusPolynomial) || funcType == 4) {
                if (params.size() < 2) {
                    return false;
                }

                auto sharedIt = massDepPars.find(SharedResonanceKeyForWave(waveName));
                if (sharedIt == massDepPars.end() || sharedIt->second.size() < 6) {
                    return false;
                }

                std::vector<T> assembled;
                assembled.reserve(params.size() + 6);
                assembled.push_back(params[0]);           // k_fl
                assembled.push_back(params[1]);           // k_bw
                assembled.push_back(sharedIt->second[0]); // g_etapi
                assembled.push_back(sharedIt->second[1]); // g_KK
                assembled.push_back(sharedIt->second[2]); // g_etaprimepi
                assembled.push_back(sharedIt->second[3]); // MassFlatte
                assembled.push_back(sharedIt->second[4]); // MassBW
                assembled.push_back(sharedIt->second[5]); // WidthBW
                assembled.insert(assembled.end(), params.begin() + 2, params.end());
                params.swap(assembled);
                return true;
            }

            if (funcType != static_cast<int>(MassDependenceConfig::ModelType::CustomFormula) && funcType != 5) {
                return true;
            }

            return true;
        }

        inline bool BuildWaveParameters(
            const TString& baseName,
            const TString& waveName,
            int funcType,
            const std::map<TString, std::vector<double>>& massDepPars,
            std::vector<double>& params) {
            return BuildWaveParametersGeneric(baseName, waveName, funcType, massDepPars, params);
        }

    } // namespace internal

    // Bring internal symbols into m2pw namespace for backward compatibility
    using namespace internal;

} // namespace m2pw

#endif // FIT_INTERNAL_UTILS_H

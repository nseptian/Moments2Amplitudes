#pragma once

#include "ParameterHelper.h"
#include "MomentHelper.h"
#include "MassDependentMoments.h"
#include "Equation.h"
#include "Setup.h"
#include <TString.h>
#include <Fit/Fitter.h>
#include "MassDependentFunction.h"
#include <TObjString.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <string>
#include <TRandom3.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TStopwatch.h>

namespace m2pw {

    using equation_t = std::vector<Equation>;
    using HS::FIT::Setup;

    class MassDependentFitter {
    public:
        /**
         * @brief Configuration for mass dependence of different l values
         */
        struct MassDependenceConfig {
            std::map<int, std::vector<TString>> massDependentWaves;  // l -> list of wave names
            std::vector<int> massIndependentL;                       // l values that are mass independent
            
            // Constructor with default behavior
            MassDependenceConfig() {
                massDependentWaves[2] = {"a2_1320"};  // Default: only l=2 with single wave
                massIndependentL = {0, 1};            // l=0,1 mass independent
            }
            
            // Check if a given l value is mass dependent
            bool IsMassDependent(int l) const {
                return massDependentWaves.find(l) != massDependentWaves.end();
            }
            
            // Check if a given l value is mass independent
            bool IsMassIndependent(int l) const {
                return std::find(massIndependentL.begin(), massIndependentL.end(), l) != massIndependentL.end();
            }

            static bool HasPhaseTag(const TString& waveName) {
                return waveName.EndsWith("_phase");
            }

            static TString StripPhaseTag(const TString& waveName) {
                if (!HasPhaseTag(waveName)) {
                    return waveName;
                }

                TString stripped = waveName;
                stripped.Remove(stripped.Length() - 6);
                return stripped;
            }

            bool UseGlobalPhaseForWave(int l, const TString& waveName) const {
                auto it = massDependentWaves.find(l);
                if (it == massDependentWaves.end()) {
                    return false;
                }

                for (const TString& configuredWave : it->second) {
                    if (StripPhaseTag(configuredWave) == waveName && HasPhaseTag(configuredWave)) {
                        return true;
                    }
                }

                return false;
            }

            bool IsGlobalPhaseEnabledForL(int l) const {
                auto it = massDependentWaves.find(l);
                if (it == massDependentWaves.end()) {
                    return false;
                }

                return std::any_of(it->second.begin(), it->second.end(),
                                   [](const TString& waveName) { return HasPhaseTag(waveName); });
            }
            
            // Get wave names for a given l value
            std::vector<TString> GetWavesForL(int l) const {
                auto it = massDependentWaves.find(l);
                if (it == massDependentWaves.end()) {
                    return {};
                }

                std::vector<TString> normalizedWaves;
                normalizedWaves.reserve(it->second.size());
                for (const TString& waveName : it->second) {
                    normalizedWaves.push_back(StripPhaseTag(waveName));
                }
                return normalizedWaves;
            }
        };

        /**
         * @brief Configuration for selecting which H moments to include in chi2
         */
        struct MomentsConfig {
            std::vector<int> includedL;     // L values of moments to include in chi2
            std::vector<int> excludedL;     // L values of moments to exclude from chi2
            bool includeAll;                // If true, include all moments (ignore other settings)
            
            // Constructor with default behavior (include all)
            MomentsConfig() : includeAll(true) {}
            
            // Constructor to include specific L values
            MomentsConfig(const std::vector<int>& included) 
                : includedL(included), includeAll(false) {}
            
            // Check if a given L value should be included in chi2
            bool ShouldIncludeL(int L) const {
                if (includeAll) return true;
                
                // Check exclusion list first
                if (std::find(excludedL.begin(), excludedL.end(), L) != excludedL.end()) {
                    return false;
                }
                
                // If inclusion list is empty, include all except excluded
                if (includedL.empty()) return true;
                
                // Check inclusion list
                return std::find(includedL.begin(), includedL.end(), L) != includedL.end();
            }
        };

        // Parameter mapping utilities
        struct ParameterInfo {
            TString name;
            TString refl_str;
            TString l_str;
            TString m_str;
            bool isPhase;
            
            ParameterInfo(const TString& parName);
        };
        
        // Constructor with both mass dependence and moments configuration
        MassDependentFitter(const Setup& setup, 
                                   std::vector<double> mass_bins, 
                                   double l_max, 
                                   double noise, 
                                   std::vector<TString> noUse,
                                   const MassDependenceConfig& massDepConfig,
                                   const MomentsConfig& hMomentsConfig);
        
        ~MassDependentFitter() = default;

        unsigned int NDim() const;
        
        /**
         * @brief Evaluate chi2 with current parameter values
         * @return Current chi2 value
         */
        double EvalChi2() const;
        
        /**
         * @brief Print which moments are included in chi2 calculation
         */
        void PrintIncludedMoments() const;

        void PrintEquations(const TString opt = "", const double mass_bin_center = 0.0) const;
        void SetEquationValues(MassDependentMoments& moms);
        void PrintParNameIndices() const;
        void PrintParCurrentVals(double mass_bin_center) const;     

        double DoEval(const std::map<TString, std::vector<double>>& massDepPars, 
                     const std::map<TString, std::map<TString, std::vector<double>>>& massIndepPars);
        
        // Getters
        const std::vector<double>& GetMassBins() const { return massBins_; }
        double GetLastChi2() const { return lastChi2_; }
        
        // Method to set/update configuration
        void SetMassDependenceConfig(const MassDependenceConfig& config);
        void SetMomentsConfig(const MomentsConfig& config);
        
        // Get current configuration
        const MassDependenceConfig& GetMassDependenceConfig() const { return massDependenceConfig_; }
        const MomentsConfig& GetMomentsConfig() const { return hMomentsConfig_; }

        // Static factory methods for common configurations
        static MassDependenceConfig CreateDefaultConfig() {
            return MassDependenceConfig(); // l=2 mass dependent with a2_1320, l=0,1 mass independent
        }

        static MassDependenceConfig CreateDWaveTwoBreitWignerConfig() {
            MassDependenceConfig config;
            // Single D-wave model with two resonances: a2(1320) + a2(1700).
            // Only one global phase is introduced and it rotates the a2(1700) term.
            config.massDependentWaves[2] = {"a2_1320_a2_1700_phase"};
            config.massIndependentL = {0, 1};
            return config;
        }
        
        static MassDependenceConfig CreateCustomConfig(const std::map<int, std::vector<TString>>& massDep,
                                                      const std::vector<int>& massIndep) {
            MassDependenceConfig config;
            config.massDependentWaves = massDep;
            config.massIndependentL = massIndep;
            return config;
        }

        // Static factory methods for H moment configurations
        static MomentsConfig CreateL4OnlyConfig() {
            return MomentsConfig({4});
        }
        
        static MomentsConfig CreateL2L4OnlyConfig() {
            return MomentsConfig({2, 4});
        }

        static MomentsConfig CreateL0L2L4OnlyConfig() {
            return MomentsConfig({0, 2, 4});
        }
        
        static MomentsConfig CreateIncludeAllConfig() {
            return MomentsConfig(); // Include all by default
        }

        // Parameter management and minimization
        struct ParameterManager {
            std::vector<TString> parIndexNames;
            std::map<TString, double> parsList;
            std::vector<TString> fixedParNames;  // For fixed parameters
            std::map<TString, int> nameToIndex;
            std::map<TString, std::pair<double, double>> parameterLimits;  // Optional per-parameter limits
            std::vector<double> parErrors;  // Store parameter errors
            int totalNpars = 0;
            int randomSeed = -1;
            MassDependentFitter& fitter;  // Reference to parent fitter

            // Constructor
            ParameterManager(MassDependentFitter& f) : fitter(f) {}

            // Add mass-independent parameters with random initialization
            void AddMassIndependentParameters(int seed);

            // Add mass-independent parameters for specific L values from a result tree
            void AddMassIndependentParametersForL(const std::vector<int>& targetL,
                                                 const TString& resultTreeFile);

            void AddMassIndependentParametersForL(const std::vector<int>& targetL,
                                                 const std::vector<double>& massBins,
                                                const TString& resultTreeFile);

            void AddMassIndependentParametersForL(const std::vector<int>& targetL,
                                                 const std::vector<double>& massBins,
                                                const std::map<std::string, std::pair<double,double>>& init_values);
            
            void AddMassIndependentParametersForPhase(const std::vector<TString>& targetVariables);

            void AddMassIndependentParametersForPhase(const TString& resultTreeFile,
                                                 const std::vector<TString>& targetVariables);

            // Add fixed mass-independent parameters for specific L values
            void AddFixedMassIndependentParametersForL(const std::vector<int>& fixedL,
                                                     const std::map<int, double>& fixedValues,
                                                     bool magnitudeOnly = true);

            // Add fixed mass-independent parameters for specific L and reflectivity (e.g., "1+", "2-")
            void AddFixedMassIndependentParametersForL(const std::vector<std::string>& lReflectivities,
                                                     const std::map<std::string, double>& fixedValues,
                                                     bool magnitudeOnly = true);

            void AddMassDependentParameters(int seed);

            // Add parameters with values from a file
            void AddMassDependentParameters(const TString filePath,
                                          const bool isFixed);

            // Add mass dependent parameters for specific waves with random initialization
            void AddMassDependentParametersForL(const std::vector<int>& targetL,
                                          const int seed,
                                          const bool magnitudeOnly,
                                          const bool isFixed,
                                          const bool yieldOnly = false);

            // Add mass dependent parameters for specific waves from a result tree
            void AddMassDependentParametersForL(const std::vector<int>& targetL,
                                          const TString filePath,
                                          const bool magnitudeOnly,
                                          const bool isFixed,
                                          const bool yieldOnly = false);

            void AddMassDependentParametersForL(const std::vector<TString>& targetL,
                                          const int seed,
                                          const bool magnitudeOnly,
                                          const bool isFixed,
                                          const bool yieldOnly = false);
            
            void AddMassDependentParametersForL(const std::vector<TString>& targetL,
                                          const TString filePath,
                                          const bool magnitudeOnly,
                                          const bool isFixed,
                                          const bool yieldOnly = false);

            bool SetParameterInitialValue(const TString& parName, double value);
            bool SetInitialValue(const TString& parName, double value, bool isFixed);
            bool SetParameterLimits(const TString& parName, double lower, double upper);
            void SetParameterInitialValues(const std::map<TString, double>& values);
            void SetParameterLimitsBatch(const std::map<TString, std::pair<double, double>>& limits);

            std::map<std::string, std::pair<double,double>> GetAmplitudeValuesAtMassBins(
                const std::vector<int>& targetL,
                double massBinCenter);

            std::vector<double> GetValues() const;
            
            // Store results from minimizer including errors
            void StoreResults(const ROOT::Math::Minimizer* minimizer);
        };

        void MakeResultTree(const ParameterManager& paramManager, const TString& fileName) const;

        // Chi2 function class for minimization
        class Chi2Function {
        private:
            MassDependentFitter& fitter_;
            const std::vector<TString>& parIndexNames_;
            mutable std::map<TString, std::vector<double>> massDepPars_;
            mutable std::map<TString, std::map<TString, std::vector<double>>> massIndepPars_;

        public:
            Chi2Function(MassDependentFitter& fitter, 
                        const std::vector<TString>& parIndexNames)
                : fitter_(fitter), parIndexNames_(parIndexNames) {}

            double operator()(const double* mass_dep_pars) const;

        private:
            void MapParameters(const double* mass_dep_pars) const;
            void MapMassIndependentParameter(const TString& parName, double value) const;
            void MapMassDependentParameter(const TString& parName, double value) const;
        };

        // Utility methods
        static std::vector<TString> GetParNames();
        void MinimizeChi2(ParameterManager& paramManager);

private:
    std::map<double, ParameterHelper> pars_;
    std::map<double, equation_t> eqns_;
    std::vector<MassDependentFunction> massDepFuncs_;
    std::vector<double> massBins_;
    std::map<TString, int> parNameToIndex_;
    MassDependenceConfig massDependenceConfig_;
    MomentsConfig hMomentsConfig_;  // Configuration for moment selection
    std::map<TString, int> waveToFunctionIndex_;  // Maps wave names to function indices

    // Caching for performance
    mutable double lastChi2_ = 0.0;
    mutable std::vector<double> cachedParams_;
    mutable bool cacheValid_ = false;

    // Minimizer status and validity
    bool minimizerIsValid_ = false;
    int minimizerStatus_ = -1;

    // Helper methods
    void InitializeMassBins(const std::vector<double>& mass_bins);
    void InitializeMDFunctions();

    std::unique_ptr<ParameterInfo> ParseParameterName(const TString& parName) const;

    // Chi2 calculation
    double EvaluateChi2ForMassBin(double massBin, ParameterHelper& pars) const;

    // Helper methods for amplitude-level evaluation (single or multiple waves)
    double GetAmplitudeMagnitude(double mass_bin,
                                const TString& par_name,
                                const std::map<TString, std::vector<double>>& massDepPars,
                                int l_value) const;

    double GetAmplitudePhase(double mass_bin,
                            const TString& base_name,
                            const std::map<TString, std::vector<double>>& massDepPars,
                            int l_value) const;

    // Helper method to map wave names to function indices
    int GetFunctionIndexForWave(const TString& waveName) const;

    // Helper method to extract L value from equation name (format: "H_alpha_L_M")
    int ExtractLFromEquationName(const TString& eqnName) const;

    // Helper method to determine if parameters for a given L are needed based on H moments config
    bool ParameterNeededForMoments(int l, const MomentsConfig& hConfig) const;
    };

    // Implementation of inline methods
    inline MassDependentFitter::ParameterInfo::ParameterInfo(const TString& parName) {
        name = parName;
        std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
        if (parts->GetEntries() >= 3) {
            refl_str = ((TObjString*)parts->At(0))->GetString();
            l_str = ((TObjString*)parts->At(1))->GetString();
            m_str = ((TObjString*)parts->At(2))->GetString();
            isPhase = (refl_str == "aphi" || refl_str == "bphi");
        }
    }
   
}
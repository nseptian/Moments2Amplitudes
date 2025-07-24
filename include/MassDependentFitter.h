#pragma once

#include "ParameterHelper.h"
#include "MomentHelper.h"
#include "MassDependentMoments.h"
#include "Equation.h"
#include "Setup.h"
#include <TString.h>
#include <Math/IParamFunction.h>
#include <Fit/Fitter.h>
#include "MassDependentFunction.h"
#include <TObjString.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <map>
#include <set>
#include <memory>
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
            
            // Get wave names for a given l value
            std::vector<TString> GetWavesForL(int l) const {
                auto it = massDependentWaves.find(l);
                return (it != massDependentWaves.end()) ? it->second : std::vector<TString>{};
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
        
        // Constructor with both mass dependence and H moments configuration
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
         * @brief Print which H moments are included in chi2 calculation
         */
        void PrintIncludedMoments() const;

        // Custom methods
        void PrintEquations(const TString opt = "", const double mass_bin_center = 0.0) const;
        void SetEquationValues(MassDependentMoments& moms);
        void PrintParNameIndices() const;
        void PrintParCurrentVals(double mass_bin_center) const;
        
        // Optimized evaluation methods
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
        
        static MassDependenceConfig CreateMultipleA2Config() {
            MassDependenceConfig config;
            config.massDependentWaves[2] = {"a2_1320", "a2_1700"};
            config.massIndependentL = {0, 1};
            return config;
        }
        
        static MassDependenceConfig CreateL0L2Config() {
            MassDependenceConfig config;
            config.massDependentWaves[0] = {"a0_980"};
            config.massDependentWaves[2] = {"a2_1320", "a2_1700"};
            config.massIndependentL = {1};
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
        
        static MomentsConfig CreateIncludeAllConfig() {
            return MomentsConfig(); // Include all by default
        }

        // Parameter management and minimization
        struct ParameterManager {
            std::vector<TString> parIndexNames;
            std::map<TString, double> parsList;
            std::vector<TString> fixedParNames;  // For fixed parameters
            std::map<TString, int> nameToIndex;
            int totalNpars = 0;
            int randomSeed = -1;

            // Add mass-independent parameters with random initialization
            void AddMassIndependentParameters(const std::vector<double>& massBins, 
                                            const std::vector<TString>& parNames, 
                                            int seed,
                                            const MassDependenceConfig& config,
                                            const MomentsConfig& hConfig,
                                            const MassDependentFitter& fitter);
            void AddMassDependentParameters(const std::vector<TString>& parNames, 
                                          int seed,
                                          const MassDependenceConfig& config,
                                          const MomentsConfig& hConfig,
                                          const MassDependentFitter& fitter);

            // Add parameters with values from a file
            void AddMassDependentParameters(const std::vector<TString>& parNames,
                                          const TString filePath,
                                          const MassDependenceConfig& config,
                                          const MomentsConfig& hConfig,
                                          const MassDependentFitter& fitter,
                                          const bool isFixed);

            std::vector<double> GetValues() const;
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
    MomentsConfig hMomentsConfig_;  // Configuration for H moment selection
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

    // Helper methods for single wave evaluation (no coherent sum)
    double GetSingleWaveMagnitude(double mass_bin, 
                                const TString& par_name,
                                const std::map<TString, std::vector<double>>& massDepPars,
                                const TString& waveName) const;

    double GetSingleWavePhase(double mass_bin, 
                            const TString& base_name,
                            const std::map<TString, std::vector<double>>& massDepPars,
                            const TString& waveName) const;

    // Helper methods for coherent amplitude combination (multiple waves)
    double GetCoherentMagnitudeForL(double mass_bin, 
                                   const TString& par_name,
                                   const std::map<TString, std::vector<double>>& massDepPars,
                                   int l_value) const;

    double GetCoherentPhaseForL(double mass_bin, 
                               const TString& base_name,
                               const std::map<TString, std::vector<double>>& massDepPars,
                               int l_value) const;

    // Helper method to map wave names to function indices
    int GetFunctionIndexForWave(const TString& waveName) const;

    // Helper method to extract L value from equation name (e.g., "H_0_0" -> 0)
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
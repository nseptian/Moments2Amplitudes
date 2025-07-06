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
#include <memory>
#include <TRandom3.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TStopwatch.h>

namespace m2pw {
    
    // Configuration constants
    namespace Config {
        constexpr double DEFAULT_CHI2_TOLERANCE = 1e-6;
        constexpr size_t MAX_CACHE_SIZE = 1000;
        const TString BW_NAMES[3] = {"k", "M", "width"};
    }

    using equation_t = std::vector<Equation>;
    using HS::FIT::Setup;

    class MassDependentEquationSolver {
    public:
        // Parameter mapping utilities
        struct ParameterInfo {
            TString name;
            TString refl_str;
            TString l_str;
            TString m_str;
            bool isPhase;
            
            ParameterInfo(const TString& parName);
        };

        MassDependentEquationSolver(const Setup& setup, std::vector<double> mass_bins, 
                                   double l_max, double noise, std::vector<TString> noUse);
        ~MassDependentEquationSolver() = default;

        // Core evaluation methods
        double DoEval(const double* mass_dep_pars);
        unsigned int NDim() const;

        // Custom methods
        void PrintEquations(const TString opt = "", const double mass_bin_center = 0.0) const;
        void SetEquationValues(MassDependentMoments& moms);
        void PrintParNameIndices() const;
        void PrintParCurrentVals(double mass_bin_center) const;
        
        // Optimized evaluation methods
        double DoEval(const std::map<TString, std::vector<double>>& massDepPars, 
                     const std::map<TString, std::map<TString, std::vector<double>>>& massIndepPars);
        
        // Memory efficient file operations
        void MakeResultTree(const TString& fileName) const;
        void MakeResultTree(const TString& fileName, int seed) const;
        
        // Getters
        const std::vector<double>& GetMassBins() const { return massBins_; }
        double GetLastChi2() const { return lastChi2_; }
        
        // Parameter management and minimization
        struct ParameterManager {
            std::vector<TString> parIndexNames;
            std::map<TString, double> parsList;
            std::map<TString, int> nameToIndex;
            int totalNpars = 0;

            void AddMassIndependentParameters(const std::vector<double>& massBins, 
                                            const std::vector<TString>& parNames, 
                                            int seed = 0);
            void AddMassDependentParameters(const std::vector<TString>& parNames, 
                                          int seed = 0);
            std::vector<double> GetInitialValues() const;
        };

        // Chi2 function class for minimization
        class Chi2Function {
        private:
            MassDependentEquationSolver& solver_;
            const std::vector<TString>& parIndexNames_;
            mutable std::map<TString, std::vector<double>> massDepPars_;
            mutable std::map<TString, std::map<TString, std::vector<double>>> massIndepPars_;

        public:
            Chi2Function(MassDependentEquationSolver& solver, 
                        const std::vector<TString>& parIndexNames)
                : solver_(solver), parIndexNames_(parIndexNames) {}

            double operator()(const double* mass_dep_pars) const;

        private:
            void MapParameters(const double* mass_dep_pars) const;
            void MapMassIndependentParameter(const TString& parName, double value) const;
            void MapMassDependentParameter(const TString& parName, double value) const;
        };

        // Utility methods
        static std::vector<TString> GetParNames();
        void MinimizeChi2(int seed = 0);
        void MinimizeChi2(const ParameterManager& paramManager, int seed = 0);

    private:
        // Core data structures - using consistent naming
        std::map<double, ParameterHelper> pars_;
        std::map<double, equation_t> eqns_;
        std::vector<MassDependentFunction> massDepFuncs_;
        std::vector<double> massBins_;
        std::map<TString, int> parNameToIndex_;
        
        // Caching for performance
        mutable double lastChi2_ = 0.0;
        mutable std::vector<double> cachedParams_;
        mutable bool cacheValid_ = false;
        
        // Helper methods
        void InitializeMassBins(const std::vector<double>& mass_bins);
        void InitializeFunctions(double l_max);
        void ClearCache() const;
        bool IsParameterCached(const double* params) const;
        
        std::unique_ptr<ParameterInfo> ParseParameterName(const TString& parName) const;
        
        // Chi2 calculation
        double EvaluateChi2ForMassBin(double massBin, ParameterHelper& pars) const;

    };

    // Implementation of inline methods
    inline MassDependentEquationSolver::ParameterInfo::ParameterInfo(const TString& parName) {
        name = parName;
        std::unique_ptr<TObjArray> parts(parName.Tokenize("_"));
        if (parts->GetEntries() >= 3) {
            refl_str = ((TObjString*)parts->At(0))->GetString();
            l_str = ((TObjString*)parts->At(1))->GetString();
            m_str = ((TObjString*)parts->At(2))->GetString();
            isPhase = (refl_str == "aphi" || refl_str == "bphi");
        }
    }

    inline void MassDependentEquationSolver::ClearCache() const {
        cacheValid_ = false;
        cachedParams_.clear();
    }

    inline bool MassDependentEquationSolver::IsParameterCached(const double* params) const {
        if (!cacheValid_ || cachedParams_.empty()) return false;
        
        const size_t nParams = NDim();
        if (cachedParams_.size() != nParams) return false;
        
        for (size_t i = 0; i < nParams; ++i) {
            if (std::abs(cachedParams_[i] - params[i]) > Config::DEFAULT_CHI2_TOLERANCE) {
                return false;
            }
        }
        return true;
    }

}

// extern MassDependentEquationSolver *geqn_solver;

// extern "C" {
//   double cpp_eval(const double * params, size_t d=0) {
//     Double_t result= geqn_solver->DoEval(params);
//     return result;
//   }
// };
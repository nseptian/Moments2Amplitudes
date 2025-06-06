#pragma once

#include "ParameterHelper.h"
#include "MomentHelper.h"
#include "MassDependentMoments.h"
#include "Equation.h"
#include "Setup.h"
#include <TString.h>
#include <Math/IParamFunction.h>
#include <Fit/Fitter.h>
#include <vector>
#include <map>

namespace m2pw {

    using equation_t = std::vector<Equation>;
    using HS::FIT::Setup;

    class MassDependentEquationSolver : public ROOT::Math::IMultiGenFunction {
    public:
        MassDependentEquationSolver(const Setup& setup, std::vector<Double_t> mass_bins, Double_t noise, std::vector<TString> noUse);
        ~MassDependentEquationSolver() override;

        ROOT::Math::IMultiGenFunction* Clone() const override {
            return nullptr; // Not implemented
        }

        void PrintEquations(const TString opt = "", const double mass_bin_center = 0.0) const;
        void SetEquationValues(MassDependentMoments& moms);

        double DoEval(const double* mass_dep_pars) const override;
        unsigned int NDim() const override {
            if (_pars.empty()) return 0;
            unsigned int n = 0;
            for (const auto& [mass_bin, pars] : _pars) {
                n += pars.Nvars();
            }
            // cout << "MassDependentEquationSolver::NDim() : Total number of parameters = " << n << endl;
            return n;
        }

    private:

        std::map<double, ParameterHelper> _pars;
        std::map<double, equation_t> _eqns;
        ROOT::Fit::Fitter _fitter;
        mutable double _chi2 = 0;
    };
}

// extern MassDependentEquationSolver *geqn_solver;

// extern "C" {
//   double cpp_eval(const double * params, size_t d=0) {
//     Double_t result= geqn_solver->DoEval(params);
//     return result;
//   }
// };
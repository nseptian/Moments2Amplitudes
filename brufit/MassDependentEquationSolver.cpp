#include "MassDependentEquationSolver.h"

namespace m2pw {

    MassDependentEquationSolver::MassDependentEquationSolver(const Setup& setup, vector<Double_t> mass_bins, Double_t noise, vector<TString> noUse){
        // Initialize equations for each mass bin
        const auto validEqn = [&noUse](const TString& name) {
            if (name.BeginsWith("a_")) return false; // Skip amplitude equations
            if (name == "normalise") return false; // Skip a_0_0
            for(const auto& match : noUse) {
               if(name.Contains(match)) return false;
            }
            return true;
        };

        for (const double& mass_bin : mass_bins) {
            // cout << endl << "=========================";
            // cout << "MassDependentEquationSolver::MassDependentEquationSolver() : massbin = " << mass_bin << endl;
            _pars[mass_bin] = ParameterHelper(setup);
            for (const auto& form : setup.ParameterFormulas()) {
                if (!validEqn(form->GetName())) continue; // Skip if not valid equation
                RooFormulaVar* var = dynamic_cast<RooFormulaVar*>(form);
                _eqns[mass_bin].push_back(Equation(var, &_pars[mass_bin], noise));
            }

            for (auto& eqn : _eqns[mass_bin]) {
                // eqn.FindDependencies(_eqns[mass_bin]);
                eqn.FindL();
                // cout << "MassDependentEquationSolver::MassDependentEquationSolver() : Evaluating equation: " << eqn.GetName() << endl;
                const auto& par_vals = _pars[mass_bin].CurrentVals();
                // cout << "MassDependentEquationSolver::MassDependentEquationSolver() : Current values: " << *par_vals << endl;
                eqn.DoEval(_pars[mass_bin].CurrentVals());
            }
            // cout << "Mass bin " << mass_bin << " has " << _eqns[mass_bin].size() << " equations." << endl;
        }

        auto getAllCurrentVals = [this]() -> double* {
            size_t total = 0;
            for (const auto& [mass_bin, helper] : _pars)
                total += helper.Nvars();

            double* arr = new double[total];
            size_t idx = 0;
            for (const auto& [mass_bin, helper] : _pars) {
                const double* vals = helper.CurrentVals();
                size_t n_vars = helper.Nvars();
                for (size_t i = 0; i < n_vars; ++i) {
                    arr[idx++] = vals[i];
                    cout << "idx = " << idx << endl;    
                }
            }
            return arr;
        };

        double* all_current_vals = getAllCurrentVals();

        _fitter.Config().SetParamsSettings(NDim(), all_current_vals);
        for (const auto& [mass_bin, helper] : _pars) {
            for (unsigned int i = 0; i < helper.Nvars(); ++i) {
                _fitter.Config().ParSettings(i).SetLimits(helper.Min(i), helper.Max(i));
            }
        }

        // if(geqn_solver==nullptr)geqn_solver=this;
        
    }


    MassDependentEquationSolver::~MassDependentEquationSolver() = default;

    void MassDependentEquationSolver::PrintEquations(const TString opt, const double mass_bin_center) const {
        // Print equations for specific mass bin center
        cout << "MassDependentEquationSolver::PrintEquations() : mass_bin_center = " << mass_bin_center << endl;
        auto it = _eqns.find(mass_bin_center);
        if (it != _eqns.end()) {
            for (const auto& eqn : it->second) {
                eqn.Print(opt);
            }
        } else {
            cout << "No equations found for mass bin center: " << mass_bin_center << endl;
        }
    }

    void MassDependentEquationSolver::SetEquationValues(MassDependentMoments& moms) {
        // Set the equation values based on the mass-dependent moments
        for (auto& [mass_bin, equations] : _eqns) {
            for (auto& eqn : equations) {
                eqn.FindL();
                // Set the moment to the value given in MassDependentMoments
                eqn.SetEquationValue(moms.GetMoments(mass_bin).GetVal(eqn.GetName()));
                eqn.DoEval(_pars[mass_bin].CurrentVals());
            }
        }
    }

    pair<double,double> BW(double* mass_dep_pars, double mass_bin_center) {
        // Example of a Breit-Wigner function evaluation
        double k = mass_dep_pars[0]; // Magnitude
        double M = mass_dep_pars[1]; // Mass
        double width = mass_dep_pars[2]; // Width

        double a = mass_bin_center * mass_bin_center - M * M;
        double b = M * width;

        double Re_BW = a * k / (a * a + b * b);
        double Im_BW = -b * k / (a * a + b * b);

        double magnitude = TMath::Sqrt(Re_BW * Re_BW + Im_BW * Im_BW);
        double phase = TMath::ATan2(Im_BW, Re_BW);
        return make_pair(magnitude, phase);
    }

    double MassDependentEquationSolver::DoEval(const double* mass_dep_pars) const {

        return 0.0; // Placeholder return value
    }

    // void MassDependentEquationSolver::SetEquationValues(MassDependentMoments& moms) {
    //     // Set the equation values based on the mass-dependent moments
    //     for (const auto& [mass_bin, equations] : _eqns) {
    //         for (auto& eqn : equations) {
    //             eqn.FindL();
    //             eqn.SetEquationValue(moms.GetMoments(mass_bin).GetVal(eqn.GetName()));
    //             eqn.DoEval(_pars[mass_bin].CurrentVals());
    //         }
    //     }
    // }
}
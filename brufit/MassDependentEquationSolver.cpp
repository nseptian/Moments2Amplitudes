#include "MassDependentEquationSolver.h"
#include "TObjString.h"

namespace m2pw {

    MassDependentEquationSolver::MassDependentEquationSolver(const Setup& setup, vector<double> mass_bins, double l_max, double noise, vector<TString> noUse): _massBins(mass_bins) {
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
        
        _massDepFuncs.emplace_back(mass_bins.size(), mass_bins.front(), mass_bins[1] - mass_bins[0], 0);

        int index = 0;
        for (int l=0; l <= l_max; ++l) {
            for (int m=-l; m <= l; ++m) {
                TString par_name = Form("a_%d_%d", l, m);
                _parNameToIndex[par_name] = index++;
            }
        }
        for (int l=0; l <= l_max; ++l) {
            for (int m=-l; m <= l; ++m) {
                TString par_name = Form("b_%d_%d", l, m);
                _parNameToIndex[par_name] = index++;
            }
        }

        // auto GetMomentCurrentVals = [this]() -> double* {
        //     size_t total = 0;
        //     for (const auto& [mass_bin, helper] : _pars)
        //         total += helper.Nvars();

        //     double* arr = new double[total];
        //     size_t idx = 0;
        //     for (const auto& [mass_bin, helper] : _pars) {
        //         const double* vals = helper.CurrentVals();
        //         size_t n_vars = helper.Nvars();
        //         for (size_t i = 0; i < n_vars; ++i) {
        //             arr[idx++] = vals[i];
        //             cout << "idx = " << idx << endl;    
        //         }
        //     }
        //     return arr;
        // };

        // double* Mome = GetMomentCurrentVals();

        // _fitter.Config().SetParamsSettings(NDim(), all_current_vals);
        // for (const auto& [mass_bin, helper] : _pars) {
        //     for (unsigned int i = 0; i < helper.Nvars(); ++i) {
        //         _fitter.Config().ParSettings(i).SetLimits(helper.Min(i), helper.Max(i));
        //     }
        // }

        // if(geqn_solver==nullptr)geqn_solver=this;
        
    }

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
                eqn.SetEquationValue(moms.GetMoments(mass_bin).GetUnnormalizedVal(eqn.GetName()),moms.GetMoments(mass_bin).GetUnnormalizedError(eqn.GetName()));
                eqn.DoEval(_pars[mass_bin].CurrentVals());
            }
        }
    }

    // void MassDependentEquationSolver::GetUnnormalizedMomentsError(MassDependentMoments& moms) {
    //     // Set the unnormalized moments error based on the mass-dependent moments
    //     for (auto& [mass_bin, equations] : _eqns) {
    //         for (auto& eqn : equations) {
    //             eqn.FindL();
    //             // Set the unnormalized moment error to the value given in MassDependentMoments
    //             eqn.SetUnnormalizedMomentsError(moms.GetMoments(mass_bin).GetUnnormalizedMomentsError(eqn.GetName()));
    //         }
    //     }
    // }

    double MassDependentEquationSolver::DoEval(const double* mass_dep_pars){
        double total_chi2 = 0.0;
        for (int i = 0; i < _massBins.size(); ++i) {
            // Set the current values for each mass bin
            // cout << endl << "Mass bin: " << _massBins[i] << endl;
            auto& par = _pars.at(_massBins[i]);
            for (int j = 0; j < par.Nvars(); ++j) {
                TString par_name = par.GetParName(j);
                auto parts = par_name.Tokenize("_");
                TString refl_str = ((TObjString*)parts->At(0))->GetString();
                TString l_str = ((TObjString*)parts->At(1))->GetString();
                TString m_str = ((TObjString*)parts->At(2))->GetString();

                if (refl_str == "a" || refl_str == "b") {
                    int par_index = _parNameToIndex.at(par_name);
                    double mass_dep_par = mass_dep_pars[_parNameToIndex.at(par_name)];
                    
                    // if l == 0 magnitude = exp, l==2 BW magnitude
                    double magnitude = 0.0;
                    // cout << "l_str = " << l_str << endl; 
                    if (l_str == "0" || l_str == "1") {
                        par.SetCurrentVal(par_name, mass_dep_par);
                    }
                    if (l_str == "2") {
                        magnitude = _massDepFuncs[0].GetPWMagnitude(_massBins[i], mass_dep_par);
                        // cout << "MassDependentEquationSolver::DoEval() : l_str = " << l_str << ", mass_bin = " << _massBins[i] << ", par_name = " << par_name << ", magnitude = " << magnitude << ", mass_dep_par = " << mass_dep_par << endl;
                        par.SetCurrentVal(par_name, magnitude);
                    }
                    
                    // cout << par_indqqex << " Magnitude for " << par_name << " = " << magnitude << endl;
                    
                }
                else {
                    refl_str = refl_str.ReplaceAll("phi","");
                    double mass_dep_par = mass_dep_pars[_parNameToIndex.at(TString::Format("%s_%s_%s", refl_str.Data(), l_str.Data(), m_str.Data()))];
                    double phase = 0.0;
                    if (l_str=="2") phase = _massDepFuncs[0].GetPWPhase(_massBins[i], mass_dep_par);
                    par.SetCurrentVal(par_name, phase);
                    // cout << "Phase for " << par_name << q" = " << phase << endl;
                }
            }

            // print the current parameter values for the mass bin
            // cout << "Current parameter values for mass bin " << _massBins[i] << ": ";
            // for (const auto& [name, val] : _pars.at(_massBins[i]).CurrentValsMap()) {
            //     cout << name << " = " << val << ", ";
            // }

            for (int j = 0; j < _eqns.at(_massBins[i]).size(); ++j) {
                auto& eqn = _eqns.at(_massBins[i])[j];
                // cout << "Evaluating equation: " << eqn.GetName() << endl;
                // cout << "Parameter values for equation " << eqn.GetName() << ": ";
                // eqn.Print("v");
                auto current_vals = _pars.at(_massBins[i]).CurrentVals();
                // cout << "Current values: " << endl;
                // for (int k = 0; k < _pars.at(_massBins[i]).Nvars(); ++k) {
                //     cout << _pars.at(_massBins[i]).GetParName(k) << " = " << current_vals[k] << endl;
                // }
                // TString formulaStr = eqn.GetOrigFormula();
                // cout << endl << "Formula for equation " << eqn.GetName() << ": " << formulaStr << endl;
                // RooFormulaVar* formulaVar = dynamic_cast<RooFormulaVar*>(eqn.GetRooFormulaVar());
                // cout << "RooFormulaVar for equation " << eqn.GetName() << ": " << formulaVar->expression() << endl;
                // cout << "chi = " << eqn.DoEval(current_vals) << endl;
                // cout << "Equation value for " << eqn.GetName() << ": " << eqn_val << endl;
                // cout << "Equation value: " << eqn_val << endl;
                // total_chi2 += eqn_val * eqn_val; // Assuming chi2 is the sum of squares of equation values
                eqn.SetNeedsRecalc();
                total_chi2 += eqn.DoEvalSq(current_vals);

                // get moments error from _moments.UnnormalizedmomentsError
                // cout << "Equation value for " << eqn.GetName() << ": " << chi2 << endl;
                // eqn.Print("v");
            }
            
        }

        // for (auto& [mass_bin, equations] : _eqns) {
        //     _pars[mass_bin].
        //     // Evaluate each equation for the given mass-dependent parameters
            
        // }
        return total_chi2; // Placeholder return value
    }

    void MassDependentEquationSolver::PrintParCurrentVals(double mass_bin_center) const {
        // Print the current parameter values for a specific mass bin center
        auto it = _pars.find(mass_bin_center);
        if (it != _pars.end()) {
            cout << "Current parameter values for mass bin center " << mass_bin_center << ":" << endl;
            it->second.Print("v");
        } else {
            cout << "No parameters found for mass bin center: " << mass_bin_center << endl;
        }
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
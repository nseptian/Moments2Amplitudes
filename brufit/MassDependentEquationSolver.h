#pragma once

#include "ParameterHelper.h"
#include "MomentHelper.h"
#include "MassDependentMoments.h"
#include "Equation.h"
#include "Setup.h"
#include <TString.h>
#include <Math/IParamFunction.h>
#include <Fit/Fitter.h>
#include <MassDependentFunction.h>
#include <TObjString.h>
#include <vector>
#include <map>

namespace m2pw {

    using equation_t = std::vector<Equation>;
    using HS::FIT::Setup;

    class MassDependentEquationSolver {
    public:
        MassDependentEquationSolver(const Setup& setup, std::vector<double> mass_bins, double l_max, double noise, std::vector<TString> noUse);
        ~MassDependentEquationSolver() = default;

        void PrintEquations(const TString opt = "", const double mass_bin_center = 0.0) const;
        void SetEquationValues(MassDependentMoments& moms);
        void PrintParNameIndices() const {
            for (const auto& [name, index] : _parNameToIndex) {
                std::cout << "Parameter: " << name << ", Index: " << index << std::endl;
            }
        }

        double DoEval(const double* mass_dep_pars);
        void PrintParCurrentVals(double mass_bin_center) const;

        double DoEval(const map<TString,vector<double>> massDepPars, const map<TString, map<TString, vector<double>>> massIndepPars) {
            double total_chi2 = 0.0;
            for (int i = 0; i < _massBins.size(); ++i) {
                TString mass_bin_str = TString::Format("%1.2f",_massBins[i]);
                // cout << "MassDependentEquationSolver::DoEval() : mass_bin_str = " << mass_bin_str << endl;
                auto& pars = _pars.at(_massBins[i]);
                // Set the current values for each mass bin
                // cout << endl << "Mass bin: " << mass_bin << endl;
                for (int j = 0; j < pars.Nvars(); ++j) {
                    TString par_name = pars.GetParName(j);
                    std::unique_ptr<TObjArray> parts(par_name.Tokenize("_"));
                    TString refl_str = ((TObjString*)parts->At(0))->GetString();
                    TString l_str = ((TObjString*)parts->At(1))->GetString();
                    TString m_str = ((TObjString*)parts->At(2))->GetString();

                    // cout << "MassDependentEquationSolver::DoEval() : l_str = " << l_str << ", mass_bin = " << mass_bin_str << ", par_name = " << par_name << endl;
                    if (!(refl_str == "aphi" || refl_str == "bphi")) {
                        double magnitude = 0.0;
                        if (l_str == "2") {
                            magnitude = _massDepFuncs[0].GetPWMagnitude(_massBins[i], massDepPars.at(par_name));
                            // cout << "par_name " << par_name << ", magnitude = " << magnitude << endl;
                        }
                        else {
                            // cout << "MassDependentEquationSolver::DoEval() : Using massIndepPars for " << par_name << endl;
                            magnitude = massIndepPars.at(mass_bin_str).at(par_name)[0];
                            // cout << "par_name " << par_name << ", magnitude = " << magnitude << endl;
                            // cout << "MassDependentEquationSolver::DoEval() : magnitude = " << magnitude << endl;
                        }
                        pars.SetCurrentVal(par_name, magnitude);
                    }
                    else {
                        // cout << "MassDependentEquationSolver::DoEval() : l_str = " << l_str << ", mass_bin = " << mass_bin_str << ", par_name = " << par_name << endl;
                        // refl_str = refl_str.ReplaceAll("phi", "");
                        double phase = 0.0;
                        // cout << "refl_str = " << refl_str << ", l_str = " << l_str << ", m_str = " << m_str << endl;
                        // cout << (l_str=="2") << endl;
                        if (l_str == "2") {
                            phase = _massDepFuncs[0].GetPWPhase(_massBins[i], massDepPars.at(par_name.ReplaceAll("phi", "")));
                            par_name = pars.GetParName(j);
                        }
                        else {
                            phase = massIndepPars.at(mass_bin_str).at(par_name)[0];
                            // cout << "MassDependentEquationSolver::DoEval() : refl_str = " << refl_str << " mass_bin = " << mass_bin_str << ", par_name = " << par_name << ", phase = " << phase << endl;
                        }
                        // cout << "MassDependentEquationSolver::DoEval() : refl_str = " << refl_str << " par_name = " << par_name << ", phase = " << phase << endl;
                        pars.SetCurrentVal(par_name, phase);
                    }
                }

                for (int j = 0; j < _eqns.at(_massBins[i]).size(); ++j) {
                    auto& eqn = _eqns.at(_massBins[i])[j];
                    eqn.SetNeedsRecalc();
                    total_chi2 += eqn.DoEvalSq(pars.CurrentVals());
                }
            }
            _chi2 = total_chi2;
            return total_chi2;
        }

        void MakeResultTree(TString fileName) {
            TFile *f = new TFile(fileName, "RECREATE");
            TTree *tree = new TTree("result", "result");
            const int NMassBins = _massBins.size();
            const int NPars = _pars.begin()->second.Nvars();
            double par_vals[NPars];

            double mass_bin_center;
            tree->Branch("mass_bin", &mass_bin_center);
            for (int i = 0; i < NPars; ++i) {
                TString par_name = _pars.begin()->second.GetParName(i);
                tree->Branch(par_name, &par_vals[i]);
            }

            cout << "Filling the result tree with parameter values..." << endl;

            for (const auto& [mass_bin, pars] : _pars) {
                mass_bin_center = mass_bin;
                for (int i = 0; i < NPars; ++i) {
                    par_vals[i] = pars.GetCurrentVal(pars.GetParName(i));
                    cout << "Mass bin: " << mass_bin << ", Parameter: " << pars.GetParName(i) << ", Value: " << par_vals[i] << endl;
                }
                tree->Fill();
            }

            double chi2 = 0.0;
            // write chi2 to root file, not in ttree
            TTree *chi2_tree = new TTree("chi2", "chi2");
            chi2_tree->Branch("chi2", &chi2);
            chi2 = 
            chi2_tree->Fill();

            f->Write();
            cout << "Result tree saved to " << fileName << endl;
            f->Close();
            
        }

        unsigned int NDim() const {
            if (_pars.empty()) return 0;
            unsigned int n = 0;
            for (const auto& [mass_bin, pars] : _pars) {
                n += pars.Nvars();
            }
            // cout << "MassDependentEquationSolver::NDim() : Total number of parameters = " << n << endl;
            return n;
        }

        vector<double> GetMassBins() const {
            return _massBins;
        }


    private:

        std::map<double, ParameterHelper> _pars;
        std::map<double, equation_t> _eqns;
        std::vector<MassDependentFunction> _massDepFuncs;
        std::vector<double> _massBins; 
        std::map<TString, int> _parNameToIndex;
        // std::map<TString, vector<double>> _massDepPars;
        // std::map<double, std::map<TString, vector<double>>> _massIndepPars;
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
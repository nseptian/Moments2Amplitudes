#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include "TStopwatch.h"
#include <memory>
#include <vector>
#include <map>

// Configuration constants
namespace Config {
    const Double_t FIRST_MASS_BIN_CENTER = 0.82; // in GeV
    const Double_t MASS_BIN_WIDTH = 0.04;
    const Int_t N_MASS_BINS = 3;
    const TString FIT_RESULTS_DIR = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t080100_m010200_MCMCN6000BI1000S08WCOV_R6.34/";
    const TString FIT_RESULTS_FILENAME = "ResultsBruMcmcCovariance.root";
    const TString BW_NAMES[3] = {"k", "M", "width"};
}

struct ParameterManager {
    std::vector<TString> parIndexNames;
    std::map<TString, double> parsList;
    std::map<TString, int> nameToIndex;
    int totalNpars = 0;

    void AddMassIndependentParameters(const std::vector<double>& massBins, 
                                    const std::vector<TString>& parNames, 
                                    int seed) {
        TRandom3 rng(seed);
        for (const double massBin : massBins) {
            for (const TString& parName : parNames) {
                if (parName.Contains("2")) continue;

                TString name = Form("MI_%1.6f_%s", massBin, parName.Data());
                double initialValue = name.Contains("_1_") ? 0.0 : rng.Uniform(-100, 100);
                
                parsList[name] = initialValue;
                nameToIndex[name] = totalNpars;
                parIndexNames.push_back(name);
                totalNpars++;
            }
        }
    }

    void AddMassDependentParameters(const std::vector<TString>& parNames, int seed) {
        for (const TString& parName : parNames) {
            if (!(parName.Contains("2")) || parName.Contains("phi_0") || 
                parName.Contains("phi_1") || parName.Contains("phi_2")) continue;

            TString name = Form("MD_%s_%s", parName.Data(), Config::BW_NAMES[0].Data());
            parsList[name] = 10.0;
            nameToIndex[name] = totalNpars;
            parIndexNames.push_back(name);
            totalNpars++;
        }
    }

    std::vector<double> GetInitialValues() const {
        std::vector<double> values(totalNpars);
        for (int i = 0; i < totalNpars; ++i) {
            values[i] = parsList.at(parIndexNames[i]);
        }
        return values;
    }
};

// Forward declarations
std::vector<TString> GetParNames();
class Chi2Function;

// Chi2Function class definition
class Chi2Function {
private:
    m2pw::MassDependentFitter& fitter_;
    const std::vector<TString>& parIndexNames_;
    mutable std::map<TString, std::vector<double>> massDepPars_;
    mutable std::map<TString, std::map<TString, std::vector<double>>> massIndepPars_;

public:
    Chi2Function(m2pw::MassDependentFitter& fitter, 
                 const std::vector<TString>& parIndexNames)
        : fitter_(fitter), parIndexNames_(parIndexNames) {}

    double operator()(const double* mass_dep_pars) const {
        // Map parameters and use the optimized DoEval
        MapParameters(mass_dep_pars);
        return fitter_.DoEval(massDepPars_, massIndepPars_);
    }

private:
    void MapParameters(const double* mass_dep_pars) const {
        massDepPars_.clear();
        massIndepPars_.clear();
        
        for (size_t i = 0; i < parIndexNames_.size(); ++i) {
            const TString& parName = parIndexNames_[i];
            
            if (parName.BeginsWith("MI_")) {
                MapMassIndependentParameter(parName, mass_dep_pars[i]);
            } else if (parName.BeginsWith("MD_")) {
                MapMassDependentParameter(parName, mass_dep_pars[i]);
            }
        }
    }

    void MapMassIndependentParameter(const TString& parName, double value) const {
        int firstUnderscore = parName.Index("_", 3);
        int secondUnderscore = parName.Index("_", firstUnderscore + 1);
        int thirdUnderscore = parName.Index("_", secondUnderscore + 1);
        
        TString parKey = parName(firstUnderscore + 1, thirdUnderscore + 2 - firstUnderscore);
        TString massBinStr = parName(3, firstUnderscore - 7);
        
        massIndepPars_[massBinStr][parKey].push_back(value);
    }

    void MapMassDependentParameter(const TString& parName, double value) const {
        int lastUnderscore = parName.Last('_');
        TString parKey = parName(3, lastUnderscore - 3);
        massDepPars_[parKey].push_back(value);
    }
};

std::vector<TString> GetParNames() {
    std::vector<TString> parNames;
    for (int l = 0; l <= 2; ++l) {
        for (int m = -l; m <= l; ++m) {
            parNames.push_back(Form("a_%d_%d", l, m));
            parNames.push_back(Form("b_%d_%d", l, m));
            parNames.push_back(Form("aphi_%d_%d", l, m));
            parNames.push_back(Form("bphi_%d_%d", l, m));
        }
    }
    return parNames;
}

void MinimizeChi2(m2pw::MassDependentFitter& fitter, int seed = 0) {
    // Setup parameter manager
    ParameterManager paramManager;
    std::vector<TString> parNames = GetParNames();
    std::vector<double> massBins = fitter.GetMassBins();
    
    paramManager.AddMassIndependentParameters(massBins, parNames, seed);
    paramManager.AddMassDependentParameters(parNames, seed);
    
    std::cout << "Total number of parameters: " << paramManager.totalNpars << std::endl;

    // Get initial parameter values
    std::vector<double> initialValues = paramManager.GetInitialValues();

    // Setup minimizer with Chi2Function wrapper
    ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
    Chi2Function chi2Function(fitter, paramManager.parIndexNames);
    ROOT::Math::Functor functor(chi2Function, paramManager.totalNpars);
    
    TStopwatch timer;
    timer.Start();

    std::cout << "Chi2 before minimization: " << chi2Function(initialValues.data()) << std::endl;

    minimizer.SetFunction(functor);
    
    // Set parameters
    for (int i = 0; i < paramManager.totalNpars; ++i) {
        const TString& parName = paramManager.parIndexNames[i];
        minimizer.SetVariable(i, parName.Data(), initialValues[i], 1.0);
    }

    minimizer.SetPrintLevel(2);
    bool isValid = minimizer.Minimize();

    double minChi2 = minimizer.MinValue();
    std::cout << "Chi2 after minimization: " << minChi2 << std::endl;

    timer.Stop();
    std::cout << "Minimization time: " << timer.RealTime() << " seconds" << std::endl;

    minimizer.PrintResults();

    if (isValid) {
        fitter.MakeResultTree(Form("MassDependentResults_%d.root", seed));
    }
    
    fitter.PrintEquations("v", Config::FIRST_MASS_BIN_CENTER);
    fitter.PrintParCurrentVals(Config::FIRST_MASS_BIN_CENTER);
}


MassDependentMoments LoadMomentsData() {
    MassDependentMoments massDepMoments;
    
    for (Int_t i = 0; i < Config::N_MASS_BINS; ++i) {
        Double_t massValue = Config::FIRST_MASS_BIN_CENTER + i * Config::MASS_BIN_WIDTH;
        const TString massBinDirName = Form("Mpi0eta%1.6f_/", massValue);
        const TString filePath = Config::FIT_RESULTS_DIR + massBinDirName + Config::FIT_RESULTS_FILENAME;

        std::cout << "Loading moments from file: " << filePath << std::endl;

        std::unique_ptr<TFile> file(TFile::Open(filePath));
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filePath << std::endl;
            continue;
        }

        // Use stack allocation instead of new
        MomentHelper moments;
        moments.Set(filePath, 1, 1);
        massDepMoments.SetMoments(massValue, moments);
    }
    
    return massDepMoments;
}

std::vector<double> GenerateMassBins() {
    std::vector<double> massBins;
    massBins.reserve(Config::N_MASS_BINS);
    
    for (Int_t i = 0; i < Config::N_MASS_BINS; ++i) {
        Double_t massValue = Config::FIRST_MASS_BIN_CENTER + i * Config::MASS_BIN_WIDTH;
        massBins.push_back(massValue);
    }
    
    return massBins;
}

void TestMassDependentMoments(int seed = 0) {
    auto& setup = ConfigureAmpsNoValues(2, 2, 2); // (Lmax, MMax, Nref)
    
    // Load moments data
    MassDependentMoments massDepMoments = LoadMomentsData();
    
    // Generate mass bins
    std::vector<double> massBins = GenerateMassBins();
    
    // Setup fitter
    m2pw::MassDependentFitter fitter{setup, massBins, 2, 0.0, {"H_3"}};
    fitter.SetEquationValues(massDepMoments);
    
    // Run minimization with new interface
    fitter.MinimizeChi2(seed);
}
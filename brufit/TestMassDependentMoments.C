#include "ConfigureAmpsNoValues.C"
#include "MassDependentMoments.h"
#include "MassDependentFunction.h"
#include "TStopwatch.h"

// pair<double magnitude, double phase> MassDependentFCN(double *mass_dep_pars, double massBinCenter){
//     double k = mass_dep_pars[0];
//     double M = mass_dep_pars[1];
//     double width = mass_dep_pars[2];

//     double a = massBinCenter*massBinCenter - M*M;
//     double b = M*width;

//     double Re_BW = a*k / (a*a + b*b);
//     double Im_BW = -b*k / (a*a + b*b);

//     double magnitude = TMath::Sqrt(Re_BW*Re_BW + Im_BW*Im_BW);
//     double phase = TMath::ATan2(Im_BW, Re_BW);
//     return make_pair(magnitude, phase);
// }

// double EvalChi2MassBin(double *mass_dep_pars, double massBinCenter, const MassDependentMoments& massDepMoments){
//     double chi2 = 0.0;
//     auto moments = massDepMoments.GetMoments(massBinCenter);

//     m2pw::EquationSolver solver{setup, 0.0, {"H_3"}};
//     solver.SetEquationValues(moments);

//     return solver.DoEvalSq()

// }

vector<TString> GetParNames() {
    vector<TString> parNames;
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


void MinimizeChi2(m2pw::MassDependentEquationSolver& solver, int seed = 0) {
    ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
    vector<TString> parNames = GetParNames();
    map<TString, vector<double>> MassDepPars;
    map<double, map<TString, vector<double>>> MassIndepPars;

    vector<double> MassBins = solver.GetMassBins();

    vector<TString> parIndexNames;

    int total_npars = 0;

    // const TString AmpName[2] = {"Magnitude", "Phase"};
    const TString BWName[3] = {"k", "M", "width"};

    map<TString, double> pars_list;
    map<TString, int> NameToIndex;

    for (const double massBin: MassBins) {
        for (const TString& parName: parNames) {
            // MassIndepPars[massBin][parName] = {0.0, 0.0};
            // total_npars+=2;
            if (parName.Contains("2")) continue;

            // TString name = Form("MI_%1.6f_%s_%s", massBin, parName.Data(), AmpName[0].Data());
            TString name = Form("MI_%1.6f_%s", massBin, parName.Data());
            pars_list[name] = TRandom3(seed).Uniform(-100,100);
            NameToIndex[name] = total_npars;
            total_npars++;
            parIndexNames.push_back(name);
            // name = Form("MI_%1.6f_%s_%s", massBin, parName.Data(), AmpName[1].Data());
            // parIndexNames.push_back(name);
            // pars_list[name] = TRandom3(seed).Uniform(-TMath::Pi(), TMath::Pi());
            // NameToIndex[name] = total_npars;
            // total_npars++;
        }
    }
    for (const TString& parName : parNames) {
        // MassDepPars[parName] = {0.0,0.13182,0.1111};
        
        if (!(parName.Contains("2")) || parName.Contains("phi_0") || parName.Contains("phi_1")) continue;
        if (parName.Contains("phi_2")) continue;
        TString name = Form("MD_%s_%s", parName.Data(), BWName[0].Data());
        

        pars_list[name] = TRandom3(seed).Uniform(-10,10);
        NameToIndex[name] = total_npars;
        total_npars++;
        parIndexNames.push_back(name);
        // name = Form("MD_%s_%s", parName.Data(), BWName[1].Data());
        // pars_list[name] = 1.3182;
        // NameToIndex[name] = total_npars;
        // total_npars++;
        // parIndexNames.push_back(name);
        // name = Form("MD_%s_%s", parName.Data(), BWName[2].Data());
        // pars_list[name] = 0.1111;
        // NameToIndex[name] = total_npars;
        // total_npars++;
        // parIndexNames.push_back(name);
        // total_npars+=3;
    }

    double* mass_dep_pars = new double[total_npars];

    for (unsigned int i= 0; i < parIndexNames.size(); ++i) {
        mass_dep_pars[i] = pars_list[parIndexNames[i]];
    }

    cout << "Total number of parameters: " << total_npars << endl;
    cout << "Parameter names:" << endl;
    for (const auto& name : parIndexNames) {
        cout << name << endl;
    }

    map<TString, vector<double>> massDepPars;
    map<TString, map<TString, vector<double>>> massIndepPars;

    ROOT::Math::Functor f([&solver,&parIndexNames,&massDepPars,&massIndepPars](const double* mass_dep_pars){
        massDepPars.clear();
        massIndepPars.clear();

        // TODO: map *mass_dep_pars to massDepPars and massIndepPars
        // cout << "Mapping mass_dep_pars to massDepPars and massIndepPars..." << endl;
        // cout << "parIndexNames.size() = " << parIndexNames.size() << endl;
        for (unsigned int i = 0; i < parIndexNames.size(); ++i) {
            TString parName = parIndexNames[i];
            // cout << i << ". ";
            if (parName.BeginsWith("MI_")) {
                double value = mass_dep_pars[i];
                // Find the position of the first underscore after "MI_"
                int firstUnderscore = parName.Index("_", 3); // after "MI"
                int secondUnderscore = parName.Index("_", firstUnderscore + 1);
                int thirdUnderscore = parName.Index("_", secondUnderscore + 1);

                // cout << firstUnderscore << ", " << secondUnderscore << ", " << thirdUnderscore << ", " << fourthUnderscore << endl;
                        
                // Parameter name is between firstUnderscore+1 and thirdUnderscore-1
                // cout << firstUnderscore << ", " << secondUnderscore << ", " << thirdUnderscore << ", " << fourthUnderscore << endl;
                TString parKey = parName(firstUnderscore + 1, thirdUnderscore + 2 - firstUnderscore);
                        
                TString massBinStr = parName(3, firstUnderscore - 7);
                // double massBin = massBinStr.Atof();
                        
                // cout << "parName = " << parName << ", massBin = " << massBinStr << ", parKey = " << parKey << ", value = " << value << endl;
                massIndepPars[massBinStr][parKey].push_back(value);
            } else if (parName.BeginsWith("MD_")) {
                // Find the last underscore (before "k", "M", or "width")
                int lastUnderscore = parName.Last('_');
                // Parameter name is between "MD_" (index 3) and last underscore
                TString parKey = parName(3, lastUnderscore - 3);
                // cout << "parName = " << parName << ", parKey = " << parKey << ", value = " << mass_dep_pars[i] << endl;
                massDepPars[parKey].push_back(mass_dep_pars[i]);
            }
        }
        // cout << "Finished mapping parameters." << endl;
        return solver.DoEval(massDepPars, massIndepPars);
        // cout << "Finished evaluating equations." << endl;

        // cout << massIndepPars[1.38]["a_0_0"][1] << endl;
        // return 0.0;

    }, total_npars);

    // ROOT::Math::Functor f([&solver](const map<TString,vector<double>> massDepPars, const map<double, map<TString, vector<double>>> massIndepPars) {
    //     return solver.DoEval(massDepPars, massIndepPars);
    // }, total_npars); 

    TStopwatch timer;
    timer.Start();

    cout << "Chi2 before minimization: " << f(mass_dep_pars) << endl;

    minimizer.SetFunction(f);

    for (unsigned int i = 0; i < parIndexNames.size(); ++i) {
        TString parName = parIndexNames[i];
        minimizer.SetVariable(i, parName.Data(), mass_dep_pars[i], 1);
        // cout << "Setting variable: " << parName << " with initial value: " << mass_dep_pars[i] << endl;
    }

    // Set initial random values for the parameters
    // const int npar = 18; // Example number of parameters
    // double initial_values[npar];
    // for (int i = 0; i < npars; ++i) {
    //     initial_values[i] = TRandom3(0).Uniform(100,1000);
    //     minimizer.SetVariable(i, Form("par_%d", i), initial_values[i], 100);
    // }

    // minimizer.SetStrategy(2); // Set strategy for minimization
    minimizer.SetPrintLevel(2); // Set print level for minimization output

    bool isValid = minimizer.Minimize();

    const double* xs = minimizer.X();

    double minChi2 = minimizer.MinValue();

    cout << "Chi2 after minimization: " << minChi2 << endl;

    timer.Stop();

    cout << "Minimization time: " << timer.RealTime() << " seconds" << endl;

    minimizer.PrintResults();

    if (isValid) solver.MakeResultTree(Form("MassDependentResults_%d.root",seed));
    // solver.MakeResultTree(Form("MassDependentResults_%d.root",seed));

    // std::cout << "Minimum chi2: " << minChi2 << std::endl;
    // std::cout << "Best-fit parameters:" << std::endl;
    // std::vector<double> best_fit(xs, xs + npars);
    // for (int i = 0; i < npars; ++i) {
    //     std::cout << "Parameter " << i << ": " << best_fit[i] << std::endl;
    // }
    solver.PrintEquations("v", 1.14);
    solver.PrintParCurrentVals(1.14);
    return;
}


void TestMassDependentMoments(int seed = 0){
    auto& setup = ConfigureAmpsNoValues(2,2,2); // (Lmax,MMax,Nref)
    // setup.SetParVal("aphi_0_0",0,kTRUE); //fix S real
    // setup.SetParVal("bphi_0_0",0,kTRUE);
    // setup.SetParVal("a_2_2",0,kTRUE);
    
    // P-waves to 0
    // setup.SetParVal("a_1_0",0,kTRUE);
    // setup.SetParVal("a_1_-1",0,kTRUE);
    // setup.SetParVal("a_1_1",0,kTRUE);
    // setup.SetParVal("b_1_0",0,kTRUE);
    // setup.SetParVal("b_1_-1",0,kTRUE);
    // setup.SetParVal("b_1_1",0,kTRUE);
    // setup.SetParVal("aphi_1_0",0,kTRUE);
    // setup.SetParVal("bphi_1_0",0,kTRUE);
    // setup.SetParVal("aphi_1_1",0,kTRUE);
    // setup.SetParVal("bphi_1_1",0,kTRUE);
    // setup.SetParVal("aphi_1_-1",0,kTRUE);
    // setup.SetParVal("bphi_1_-1",0,kTRUE);


    // auto& formulas = setup.ParameterFormulas();

    MassDependentMoments massDepMoments;

    
    const Double_t firstMassBinCenter = 0.82; // in GeV
    const Double_t massBinWidth = 0.04;
    const Int_t nMassBins = 20;
    const TString fitResultsDir = "/d/home/septian/EtaPi0Analysis/run_merged/fitMoment_GlueX1_2019_11_t010100_m010200_MCMCN6000BI1000S08WCOV_R6.34/";
    const TString fitResultsFilename = "ResultsBruMcmcCovariance.root";

    TFile *f[nMassBins];

    // map<double, vector<Equation>> massDepMomentsEqs;
    vector<double> massBinCenters;

    for (Int_t i = 0; i < nMassBins; ++i) {
        // cout << "Mass bin " << i << ": ";
        Double_t massValue = firstMassBinCenter + i * massBinWidth;
        massBinCenters.push_back(massValue);
        // cout << "Mass value = " << massValue << endl;

        const TString massBinDirName = Form("Mpi0eta%1.6f_/", massValue);
        // cout << "Mass bin directory name: " << massBinDirName << endl;
        
        const TString filePath = fitResultsDir + massBinDirName + fitResultsFilename;

        cout << "Loading moments from file: " << filePath << endl;

        f[i] = TFile::Open(filePath);
        if (!f[i] || f[i]->IsZombie()) {
            cerr << "Error opening file: " << filePath << endl;
            continue;
        }

        MomentHelper* moments = new MomentHelper();
        
        moments->Set(filePath, 1, 1);

        // cout << "Moments loaded for mass value: " << massValue << endl;
        // moments->PrintVals();
        // cout << endl;

        
        massDepMoments.SetMoments(massValue, *moments);
        // cout << "Moments set for mass value: " << massValue << endl;

    }
    // Print all moments for each mass value
    // cout << "Mass-dependent moments:" << endl;
    // massDepMoments.PrintMoments();

    m2pw::MassDependentEquationSolver solver{setup, massBinCenters, 2, 0.0, {"H_3"}};
    // solver.PrintEquations("v", firstMassBinCenter);
    solver.SetEquationValues(massDepMoments);
    // solver.PrintEquations("v", firstMassBinCenter);

    // double mass_dep_pars[19] = {100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
    //                             200, 210, 220, 230, 240, 250, 260, 270, 3}; // Example values for mass-dependent parameters
    // solver.PrintParNameIndices();
    // solver.DoEval(mass_dep_pars);
    // solver.PrintEquations("v", firstMassBinCenter);
    // std::vector<vector<double>> bestChi2Params;
    // for (int i = 0; i< 10; i++){
    //     std::vector<double> bestChi2 = MinimizeChi2(solver, 11);
    //     bestChi2Params.push_back(bestChi2);
    // }

    // double mass_dep_pars[18];
    MinimizeChi2(solver, seed);
    

    // for (int i = 0; i< best_fit.size(); i++){
    //     mass_dep_pars[i] = best_fit[i];
    // }
    // cout << "Chi2 after minimization = " << solver.DoEval(mass_dep_pars) << endl;
    // solver.PrintEquations("v", firstMassBinCenter);
    // solver.PrintParCurrentVals(firstMassBinCenter);
    // // solver.PrintParNameIndices();

    // calculate sum of pars squared
    

    // Example of how to evaluate chi2 for a specific mass bin
    // double massBinCenter = 0.82; // Example mass bin center
    // double mass_dep_pars[] = {1000};

    // Double_t sum_of_amp_sq = 0.0;

    // m2pw::MassDependentFunction massDepFunction{18, firstMassBinCenter, massBinWidth, 0};
    // for (int i = 0; i < 18; ++i) {
    //     double mag = massDepFunction.GetPWMagnitude(firstMassBinCenter, mass_dep_pars[i]);
    //     sum_of_amp_sq += 2*mag*mag;
    //     cout << "Magnitude for parameter " << i << ": " << mag << endl;
    // }
    // cout << "Sum of amplitudes squared: " << sum_of_amp_sq << endl;
    // double mag = massDepFunction.GetPWMagnitude(massBinCenter, mass_dep_pars);
    // double phase = massDepFunction.GetPWPhase(massBinCenter, mass_dep_pars);
    // cout << "Magnitude: " << mag << ", Phase: " << phase << endl;
}
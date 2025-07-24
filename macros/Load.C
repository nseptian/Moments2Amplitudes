#include <TROOT.h>

void Load(){
  gSystem->AddIncludePath("-I../include");
  gROOT->ProcessLine(".L ../src/ParameterHelper.cpp+");
  gROOT->ProcessLine(".L ../src/Equation.cpp+");
  gROOT->ProcessLine(".L ../src/EquationSolver.cpp+");
  // gROOT->ProcessLine(".L ../include/MomentHelper.h");
  // gROOT->ProcessLine(".L ../include/MassDependentMoments.h");
  gROOT->ProcessLine(".L ../src/MassDependentFitter.cpp+");
  // gROOT->ProcessLine(".L ../include/MassDependentFunction.h");
}

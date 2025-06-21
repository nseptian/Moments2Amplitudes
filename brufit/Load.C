#include <TROOT.h>

void Load(){
  gROOT->ProcessLine(".L ParameterHelper.cpp+");
  gROOT->ProcessLine(".L Equation.cpp+");
  gROOT->ProcessLine(".L EquationSolver.cpp+");
  gROOT->ProcessLine(".L MomentHelper.h+");
  gROOT->ProcessLine(".L MassDependentMoments.h+");
  gROOT->ProcessLine(".L MassDependentEquationSolver.cpp+");
  gROOT->ProcessLine(".L MassDependentFunction.h");
}

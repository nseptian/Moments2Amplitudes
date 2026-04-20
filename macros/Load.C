#include <iostream>

void CleanAclicArtifacts(const TString& srcDir) {
  gSystem->Exec(Form("rm -f %s/*ACLiC* %s/*.so %s/*.pcm %s/*.rootmap %s/*.d",
                     srcDir.Data(),
                     srcDir.Data(),
                     srcDir.Data(),
                     srcDir.Data(),
                     srcDir.Data()));
}

void Load(){
  // Prefer the run-time project directory, then fall back to the legacy env var.
  TString projectDir = "/d/grid17/septian/Moments2Amplitudes_dev/";
  // cout << "Using project directory: " << projectDir << endl;
  TString srcDir = projectDir + "/src";

  // CleanAclicArtifacts(srcDir);
  
  // // Load source files with absolute paths
  gROOT->ProcessLine(Form(".L %s/ParameterHelper.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/Equation.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/EquationSolver.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/MomentHelper.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/MassDependentFitter.cpp+", srcDir.Data()));
}
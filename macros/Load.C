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
  // Ensure BruFit environment and classes (e.g. RooHSSphHarmonic) are loaded
  if (!TClass::GetClass("RooHSSphHarmonic")) {
    TString brufitDir = gSystem->Getenv("BRUFIT");
    if (brufitDir.IsNull()) {
      if (gSystem->AccessPathName("/d/grid17/septian/brufit_dev/brufit/macros/LoadBru.C") == 0) {
        brufitDir = "/d/grid17/septian/brufit_dev/brufit";
      } else if (gSystem->AccessPathName("/d/home/septian/brufit/macros/LoadBru.C") == 0) {
        brufitDir = "/d/home/septian/brufit";
      }
    }
    if (!brufitDir.IsNull() && gSystem->AccessPathName(brufitDir + "/macros/LoadBru.C") == 0) {
      gROOT->ProcessLine(Form(".x %s/macros/LoadBru.C", brufitDir.Data()));
    }
  }

  // Determine repository root directory
  TString projectDir = gSystem->Getenv("MOMS2AMPS");
  if (projectDir.IsNull()) {
    TString currentDir = gSystem->WorkingDirectory();
    if (gSystem->AccessPathName(currentDir + "/src/MassDependentFitter.cpp") == 0) {
      projectDir = currentDir;
    } else if (gSystem->AccessPathName(currentDir + "/../src/MassDependentFitter.cpp") == 0) {
      projectDir = gSystem->DirName(currentDir);
    } else {
      projectDir = "/d/home/septian/Moments2Amplitudes";
    }
  }
  TString srcDir = projectDir + "/src";

  // CleanAclicArtifacts(srcDir);
  
  // Suppress duplicate TClassTable warnings while loading ACLiC libraries
  int prevLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;

  // Load source files with absolute paths
  gROOT->ProcessLine(Form(".L %s/ParameterHelper.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/Equation.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/EquationSolver.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/MomentHelper.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/MassDependentFitter.cpp+", srcDir.Data()));

  gErrorIgnoreLevel = prevLevel;
}
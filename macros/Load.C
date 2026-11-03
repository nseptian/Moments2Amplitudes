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
  if (TClass::GetClass("m2pw::MassDependentFitter")) return;

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
    if (!brufitDir.IsNull()) {
      gSystem->Setenv("BRUFIT", brufitDir.Data());
      if (gSystem->AccessPathName(brufitDir + "/macros/LoadBru.C") == 0) {
        gROOT->ProcessLine(Form(".x %s/macros/LoadBru.C", brufitDir.Data()));
      }
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
  TString incDir = projectDir + "/include";
  TString macDir = projectDir + "/macros";

  // Register include paths for Cling interpreter and ACLiC compilation
  gInterpreter->AddIncludePath(incDir.Data());
  gInterpreter->AddIncludePath(macDir.Data());
  gSystem->AddIncludePath(Form("-I%s", incDir.Data()));
  gSystem->AddIncludePath(Form("-I%s", macDir.Data()));

  TString brufitDirEnv = gSystem->Getenv("BRUFIT");
  if (!brufitDirEnv.IsNull()) {
    gInterpreter->AddIncludePath(brufitDirEnv + "/include");
    gInterpreter->AddIncludePath(brufitDirEnv + "/core");
    gSystem->AddIncludePath(Form("-I%s/include", brufitDirEnv.Data()));
    gSystem->AddIncludePath(Form("-I%s/core", brufitDirEnv.Data()));
  }

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
#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TClass.h>

void Compile() {
  TString projectDir = gSystem->Getenv("MOMS2AMPS");
  if (projectDir.IsNull()) {
    TString currentDir = gSystem->WorkingDirectory();
    if (gSystem->AccessPathName(currentDir + "/src/MassDependentFitter.cpp") == 0) {
      projectDir = currentDir;
    } else {
      projectDir = "/d/home/septian/Moments2Amplitudes";
    }
  }

  TString srcDir = projectDir + "/src";
  TString incDir = projectDir + "/include";
  TString macDir = projectDir + "/macros";

  // Ensure BruFit environment and classes are loaded
  TString brufitDir = gSystem->Getenv("BRUFIT");
  if (brufitDir.IsNull()) {
    if (gSystem->AccessPathName("/d/grid17/septian/brufit_dev/brufit_negative/macros/LoadBru.C") == 0) {
      brufitDir = "/d/grid17/septian/brufit_dev/brufit_negative";
    }
  }
  if (!brufitDir.IsNull()) {
    gSystem->Setenv("BRUFIT", brufitDir.Data());
    if (gSystem->AccessPathName(brufitDir + "/macros/LoadBru.C") == 0) {
      gROOT->ProcessLine(Form(".x %s/macros/LoadBru.C", brufitDir.Data()));
    }
  }

  // Register include paths
  gInterpreter->AddIncludePath(incDir.Data());
  gInterpreter->AddIncludePath(macDir.Data());
  gSystem->AddIncludePath(Form("-I%s", incDir.Data()));
  gSystem->AddIncludePath(Form("-I%s", macDir.Data()));

  if (!brufitDir.IsNull()) {
    gInterpreter->AddIncludePath(brufitDir + "/install/include");
    gInterpreter->AddIncludePath(brufitDir + "/include");
    gInterpreter->AddIncludePath(brufitDir + "/core");
    gInterpreter->AddIncludePath(brufitDir + "/utility");
    gSystem->AddIncludePath(Form("-I%s/install/include", brufitDir.Data()));
    gSystem->AddIncludePath(Form("-I%s/include", brufitDir.Data()));
    gSystem->AddIncludePath(Form("-I%s/core", brufitDir.Data()));
    gSystem->AddIncludePath(Form("-I%s/utility", brufitDir.Data()));

    gSystem->AddDynamicPath(brufitDir + "/install/lib64");
    gSystem->AddDynamicPath(brufitDir + "/install/lib");
    gSystem->AddDynamicPath(brufitDir + "/build/lib");
    gSystem->AddDynamicPath(brufitDir + "/lib");
  }

  const std::vector<TString> targets = {
    "ParameterHelper.cpp",
    "Equation.cpp",
    "EquationSolver.cpp",
    "MomentHelper.cpp",
    "MassDependentFitter.cpp"
  };

  std::cout << "[Compile.C] Starting sequential ACLiC compilation of " << targets.size() << " libraries..." << std::endl;

  for (const auto& target : targets) {
    TString fullPath = srcDir + "/" + target;
    std::cout << "[Compile.C] Compiling " << target << "++ ..." << std::endl;
    Long_t ret = gROOT->ProcessLine(Form(".L %s++", fullPath.Data()));
    if (ret != 0) {
      std::cerr << "[Compile.C] ERROR: Failed to compile " << target << " (return code " << ret << ")" << std::endl;
      gSystem->Exit(1);
      return;
    }
  }

  std::cout << "[Compile.C] All libraries successfully compiled!" << std::endl;
}

void Load(){
  // Get project directory from environment variable
  TString projectDir = gSystem->Getenv("MOMS2AMPS");
  if (projectDir == "") {
    std::cerr << "ERROR: MOMS2AMPS environment variable not set!" << std::endl;
    return;
  }
  
  TString srcDir = projectDir + "/src";
  
  // // Load source files with absolute paths
  gROOT->ProcessLine(Form(".L %s/ParameterHelper.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/Equation.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/EquationSolver.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/MomentHelper.cpp+", srcDir.Data()));
  gROOT->ProcessLine(Form(".L %s/MassDependentFitter.cpp+", srcDir.Data()));
}
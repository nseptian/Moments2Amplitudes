// Example: Setting up conformal polynomial parameters from scratch
//
// This shows how to initialize conformal polynomial parameters
// when you don't have previous fit results

#include "MassDependentFitter.h"
#include <iostream>
#include <fstream>

void SetupConformalPolynomialParameters() {
    
    std::cout << "========================================" << std::endl;
    std::cout << "Setting Up Conformal Polynomial Parameters" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Choose polynomial order
    int polynomial_order = 2;  // Start with order 2 (quadratic)
    
    std::cout << "\nPolynomial order: " << polynomial_order << std::endl;
    std::cout << "Number of complex coefficients: " << (polynomial_order + 1) << std::endl;
    std::cout << "Total parameters: " << 1 + 2*(polynomial_order + 1) << " (k + real/imag pairs)" << std::endl;
    
    // Wave name
    TString wave_name = "conformal_S0";
    
    // Initialize parameters with reasonable starting values
    // For S-wave, typically start with real values around the data scale
    
    std::vector<double> initial_params;
    
    // params[0] = k (not used, set to 1.0)
    initial_params.push_back(1.0);
    
    // Constant term (z^0): set to typical amplitude scale
    double a0_init = 100.0;   // Real part
    double b0_init = 50.0;    // Imag part
    initial_params.push_back(a0_init);
    initial_params.push_back(b0_init);
    std::cout << "\nz^0 coefficient: (" << a0_init << ", " << b0_init << "i)" << std::endl;
    
    // Higher order terms: start with smaller values
    for (int i = 1; i <= polynomial_order; ++i) {
        double scale_factor = TMath::Power(0.5, i);  // Decrease with order
        double a_i = a0_init * scale_factor * (i % 2 == 0 ? 1 : -1);  // Alternate signs
        double b_i = b0_init * scale_factor * (i % 2 == 0 ? 1 : -1);
        
        initial_params.push_back(a_i);
        initial_params.push_back(b_i);
        std::cout << "z^" << i << " coefficient: (" << a_i << ", " << b_i << "i)" << std::endl;
    }
    
    // Create parameter names
    std::cout << "\n========================================" << std::endl;
    std::cout << "Parameter Names for Configuration" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::vector<TString> param_names;
    
    // k parameter (coupling)
    param_names.push_back(Form("MD_k_%s_k", wave_name.Data()));
    std::cout << param_names.back() << " = " << initial_params[0] << std::endl;
    
    // z^0 coefficients (use _m suffix for order 0)
    param_names.push_back(Form("MD_a_0_m_%s_k", wave_name.Data()));
    param_names.push_back(Form("MD_b_0_m_%s_k", wave_name.Data()));
    std::cout << param_names[param_names.size()-2] << " = " << initial_params[1] << std::endl;
    std::cout << param_names[param_names.size()-1] << " = " << initial_params[2] << std::endl;
    
    // z^i coefficients (i >= 1)
    for (int i = 1; i <= polynomial_order; ++i) {
        // For order i, use _{i-1} suffix
        int suffix = i - 1;
        param_names.push_back(Form("MD_a_0_%d_%s_k", suffix, wave_name.Data()));
        param_names.push_back(Form("MD_b_0_%d_%s_k", suffix, wave_name.Data()));
        
        int a_idx = 1 + 2*i;
        int b_idx = 1 + 2*i + 1;
        std::cout << param_names[param_names.size()-2] << " = " << initial_params[a_idx] << std::endl;
        std::cout << param_names[param_names.size()-1] << " = " << initial_params[b_idx] << std::endl;
    }
    
    // Save to a ROOT file for later use
    std::cout << "\n========================================" << std::endl;
    std::cout << "Saving Initial Parameters" << std::endl;
    std::cout << "========================================" << std::endl;
    
    TString output_file = "conformal_initial_params.root";
    TFile* f = new TFile(output_file, "RECREATE");
    
    // Create a tree to store parameters
    TTree* tree = new TTree("initial_params", "Initial conformal polynomial parameters");
    
    TString param_name;
    double param_value;
    
    tree->Branch("param_name", &param_name);
    tree->Branch("param_value", &param_value);
    
    for (size_t i = 0; i < param_names.size(); ++i) {
        param_name = param_names[i];
        param_value = initial_params[i];
        tree->Fill();
    }
    
    tree->Write();
    f->Close();
    
    std::cout << "\nParameters saved to: " << output_file << std::endl;
    std::cout << "Number of parameters: " << param_names.size() << std::endl;
    
    // Also save as a text file for reference
    std::string txt_file = "conformal_initial_params.txt";
    std::ofstream txt(txt_file);
    txt << "# Initial parameters for conformal polynomial (order " << polynomial_order << ")\n";
    txt << "# Wave: " << wave_name << "\n";
    txt << "# Format: parameter_name = value\n\n";
    
    for (size_t i = 0; i < param_names.size(); ++i) {
        txt << param_names[i] << " = " << initial_params[i] << "\n";
    }
    txt.close();
    
    std::cout << "Parameters also saved to: " << txt_file << std::endl;
    
    // Show usage example
    std::cout << "\n========================================" << std::endl;
    std::cout << "Usage in Your Fit Script" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nIn your fit configuration:\n" << std::endl;
    std::cout << "// 1. Configure wave" << std::endl;
    std::cout << "mdConfig.massDependentWaves[0] = {\"" << wave_name << "\"};" << std::endl;
    std::cout << "\n// 2. Load initial parameters" << std::endl;
    std::cout << "fitter.AddMassDependentParametersForL({0}, \"" << output_file << "\");" << std::endl;
    std::cout << "\n// OR manually add parameters:" << std::endl;
    std::cout << "std::vector<double> params = {";
    for (size_t i = 0; i < initial_params.size(); ++i) {
        std::cout << initial_params[i];
        if (i < initial_params.size() - 1) std::cout << ", ";
    }
    std::cout << "};" << std::endl;
    std::cout << "std::vector<TString> names = {";
    for (size_t i = 0; i < param_names.size(); ++i) {
        std::cout << "\"" << param_names[i] << "\"";
        if (i < param_names.size() - 1) std::cout << ", ";
    }
    std::cout << "};" << std::endl;
    std::cout << "// Add to fitter..." << std::endl;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Setup Complete!" << std::endl;
    std::cout << "========================================" << std::endl;
}

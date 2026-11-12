// Test configuration for conformal polynomial S-wave parameterization
//
// This demonstrates how to use the conformal polynomial (FuncType = 4) 
// for mass-dependent S-wave amplitudes instead of Flatté
//
// Usage: root -l TestConformalPolynomial.C

#include "MassDependentFitter.h"
#include "MassDependentMoments.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>

void TestConformalPolynomial() {
    
    std::cout << "========================================" << std::endl;
    std::cout << "Testing Conformal Polynomial Amplitude" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Define mass bins
    std::vector<double> mass_bins;
    for (double m = 0.82; m <= 1.58; m += 0.04) {
        mass_bins.push_back(m);
    }
    
    std::cout << "\nMass bins: " << mass_bins.size() << " bins from " 
              << mass_bins.front() << " to " << mass_bins.back() << " GeV" << std::endl;
    
    // Configure mass dependence with conformal polynomial for S-wave
    m2pw::MassDependentFitter::MassDependenceConfig mdConfig;
    
    // S-wave (L=0) using conformal polynomial
    mdConfig.massDependentWaves[0] = {"conformal_S0"};  // wave name contains "conformal"
    
    // D-wave (L=2) using standard resonances
    mdConfig.massDependentWaves[2] = {"a2_1320", "a2_1700"};
    
    std::cout << "\nMass-dependent wave configuration:" << std::endl;
    std::cout << "  L=0 (S-wave): conformal_S0 (conformal polynomial, default order 2)" << std::endl;
    std::cout << "  L=2 (D-wave): a2_1320, a2_1700 (Breit-Wigner)" << std::endl;
    std::cout << "\nNote: To specify order explicitly, use wave name like:" << std::endl;
    std::cout << "  \"conformal_S0_order3\" for order-3 polynomial" << std::endl;
    std::cout << "  \"conformal_S0_order1\" for order-1 polynomial" << std::endl;
    
    // Test conformal mapping parameters
    std::cout << "\n========================================" << std::endl;
    std::cout << "Conformal Mapping Parameters" << std::endl;
    std::cout << "========================================" << std::endl;
    
    double m_eta = 0.54786;   // GeV
    double m_pi0 = 0.13498;   // GeV
    double m_threshold = m_eta + m_pi0;
    double s_min = m_threshold * m_threshold;
    
    double m_expansion = 0.98;  // a0(980) mass
    double s_0 = m_expansion * m_expansion;
    
    std::cout << "\nηπ⁰ system:" << std::endl;
    std::cout << "  Threshold: m_thr = " << m_threshold << " GeV" << std::endl;
    std::cout << "  s_min = " << s_min << " GeV²" << std::endl;
    std::cout << "\nExpansion point: a₀(980)" << std::endl;
    std::cout << "  m_0 = " << m_expansion << " GeV" << std::endl;
    std::cout << "  s_0 = " << s_0 << " GeV²" << std::endl;
    
    // Compute conformal variables for mass range
    std::cout << "\n========================================" << std::endl;
    std::cout << "Conformal Variable z(m)" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "z = (√s - √(s₀-s_min)) / (√s + √(s₀-s_min))" << std::endl;
    
    std::cout << "\nMass [GeV]    s [GeV²]     z" << std::endl;
    std::cout << "------------ ------------ -----------" << std::endl;
    
    double sqrt_s0_minus_smin = TMath::Sqrt(s_0 - s_min);
    
    for (size_t i = 0; i < mass_bins.size(); i += 5) {  // Print every 5th bin
        double m = mass_bins[i];
        double s = m * m;
        double z = (m - sqrt_s0_minus_smin) / (m + sqrt_s0_minus_smin);
        
        std::cout << std::fixed << std::setprecision(3) 
                  << m << "        " << s << "      " << z << std::endl;
    }
    
    // Test amplitude evaluation with example parameters
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Amplitude Evaluation" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Create a test function
    m2pw::MassDependentFunction conformal_func(
        static_cast<int>(mass_bins.size()),
        mass_bins.front(),
        mass_bins[1] - mass_bins[0],
        4  // FuncType = 4 for conformal polynomial
    );
    
    // Example parameters for order-2 polynomial
    // Format: [re_0, im_0, re_1, im_1, re_2, im_2]
    std::vector<double> test_params = {
        100.0,    // re_0 (real part of constant term - overall scale)
        50.0,     // im_0 (imag part of constant term)
        -20.0,    // re_1 (real part of z^1 term)
        30.0,     // im_1 (imag part of z^1 term)
        5.0,      // re_2 (real part of z^2 term)
        -10.0     // im_2 (imag part of z^2 term)
    };
    
    std::cout << "\nTest parameters (order 2 polynomial):" << std::endl;
    std::cout << "  Coefficient 0: (" << test_params[0] << ", " << test_params[1] << "i)" << std::endl;
    std::cout << "  Coefficient 1: (" << test_params[2] << ", " << test_params[3] << "i)" << std::endl;
    std::cout << "  Coefficient 2: (" << test_params[4] << ", " << test_params[5] << "i)" << std::endl;
    
    std::cout << "\nAmplitude T(m) = ∑ (a_i + i*b_i) z^i" << std::endl;
    std::cout << "\nMass [GeV]    |T(m)|     Phase[deg]" << std::endl;
    std::cout << "------------ ---------- -----------" << std::endl;
    
    for (size_t i = 0; i < mass_bins.size(); i += 5) {
        double m = mass_bins[i];
        double mag = conformal_func.GetPWMagnitude(m, test_params);
        double phase = conformal_func.GetPWPhase(m, test_params) * 180.0 / TMath::Pi();
        
        std::cout << std::fixed << std::setprecision(3)
                  << m << "        " << std::setprecision(2) << mag 
                  << "       " << std::setprecision(1) << phase << std::endl;
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Usage in Mass-Dependent Fit" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nTo use conformal polynomial in your fit:" << std::endl;
    std::cout << "\n1. Configure mass-dependent waves:" << std::endl;
    std::cout << "   mdConfig.massDependentWaves[0] = {\"conformal_S0\"};  // S-wave" << std::endl;
    std::cout << "\n2. Add mass-dependent parameters for each L:" << std::endl;
    std::cout << "   fitter.AddMassDependentParametersForL({0}, \"parameter_file.root\");" << std::endl;
    std::cout << "\n3. Parameter naming convention:" << std::endl;
    std::cout << "   MD_a_0_m_conformal_S0_re_0  -> re_0 (real part of z^0, overall scale)" << std::endl;
    std::cout << "   MD_a_0_m_conformal_S0_im_0  -> im_0 (imag part of z^0)" << std::endl;
    std::cout << "   MD_a_0_m_conformal_S0_re_1  -> re_1 (real part of z^1)" << std::endl;
    std::cout << "   MD_a_0_m_conformal_S0_im_1  -> im_1 (imag part of z^1)" << std::endl;
    std::cout << "   etc." << std::endl;
    std::cout << "\n4. Default conformal parameters (from FuncType 4):" << std::endl;
    std::cout << "   s_min = 0.466 GeV² (ηπ⁰ threshold)" << std::endl;
    std::cout << "   s_0 = 0.96 GeV² (a₀(980) resonance)" << std::endl;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Complete!" << std::endl;
    std::cout << "========================================" << std::endl;
}

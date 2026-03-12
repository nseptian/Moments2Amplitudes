#pragma once

#include "TMath.h"
#include <complex>
#include <vector>

using namespace std;

namespace m2pw{

  class MassDependentFunction{
    public:
        MassDependentFunction() = default;
        MassDependentFunction(const MassDependentFunction&) = default;
        MassDependentFunction(MassDependentFunction&&) = default;
        MassDependentFunction& operator=(const MassDependentFunction&) = default;
        MassDependentFunction& operator=(MassDependentFunction&&) = default;
        ~MassDependentFunction() = default;

        MassDependentFunction(int NMassBins, double FirstMassBinCenter, double MassBinWidth, int FuncType, int PolyOrder = 2) : 
            _NMassBins(NMassBins),
            _FirstMassBinCenter(FirstMassBinCenter),
            _MassBinWidth(MassBinWidth),
            _FuncType(FuncType),
            _PolyOrder(PolyOrder) {}

        // Single parameter version
        double GetPWMagnitude(double mass, const double params) const {
            return GetPWMagnitude(mass, vector<double>{params}, false);
        }

        double GetPWPhase(double mass, const double params) const {
            return GetPWPhase(mass, vector<double>{params}, false);
        }

        // Vector parameter version (main implementation)
        double GetPWMagnitude(double mass, const vector<double>& params, bool include_globalphase = false) const {
            complex<double> amp = GetAmplitude(mass, params, include_globalphase);
            return abs(amp);
        }

        double GetPWPhase(double mass, const vector<double>& params, bool include_globalphase = false) const {
            complex<double> amp = GetAmplitude(mass, params, include_globalphase);
            return arg(amp);
        }

        int GetNMassBins() const { return _NMassBins; }
        double GetFirstMassBinCenter() const { return _FirstMassBinCenter; }
        double GetMassBinWidth() const { return _MassBinWidth; }
        int GetFuncType() const { return _FuncType; }
        int GetPolyOrder() const { return _PolyOrder; }
        
        // Get number of parameters needed for this function type
        int GetNParameters() const {
            if (_FuncType == 4) {  // Conformal polynomial
                return 2*(_PolyOrder + 1);  // (re_i, im_i) for i=0 to N
            } else if (_FuncType == 5) {  // Flatté + conformal polynomial
                return 3 + 2*(_PolyOrder + 1);  // scale, g_etapi, g_KK + polynomial coeffs
            } else if (_FuncType == 2) {  // Flatté
                return 3;  // g1, g2, M0
            } else {  // Breit-Wigner
                return 3;  // k, M, width
            }
        }

        // Get resonance mass based on function type
        double GetResonanceMass() const {
            if (_FuncType == 0) return 1.3182;      // a2(1320)
            else if (_FuncType == 1) return 1.706;  // a2(1700)
            else if (_FuncType == 2) return 1.001;  // a0(980) bare mass
            else if (_FuncType == 3) return 1.354;  // pi1(1400)
            else if (_FuncType == 4) return 0.98;   // Conformal polynomial expansion point (a0(980))
            else if (_FuncType == 5) return 0.98;   // Flatté + conformal poly (a0(980))
            else return 1.5;  // Default
        }

        // Get resonance width based on function type
        double GetResonanceWidth() const {
            if (_FuncType == 0) return 0.107;       // a2(1320)
            else if (_FuncType == 1) return 0.380;  // a2(1700)
            else if (_FuncType == 2) return 0.075;  // a0(980) - approximate width
            else if (_FuncType == 3) return 0.330;  // pi1(1400)
            else if (_FuncType == 4) return 0.683;  // Conformal polynomial threshold (η + π⁰)
            else if (_FuncType == 5) return 0.683;  // Flatté + conformal poly threshold (η + π⁰)
            else return 0.15;  // Default
        }

    private:
        int _NMassBins = 0;
        double _FirstMassBinCenter = 0.0;
        double _MassBinWidth = 0.0;
        int _FuncType = 0;
        int _PolyOrder = 2;  // Polynomial order for conformal polynomial (FuncType=4)

        // Core amplitude calculation
        complex<double> GetAmplitude(double mass, const vector<double>& params, bool include_globalphase) const {
            if (params.empty()) {
                return complex<double>(0.0, 0.0);
            }

            // Check mass bin bounds
            int bin = static_cast<int>((mass - _FirstMassBinCenter) / _MassBinWidth);
            if (bin < 0 || bin >= _NMassBins) {
                return complex<double>(0.0, 0.0);
            }

            vector<double> coreParams = params;
            double phase = 0.0;
            if (include_globalphase && params.size() > 1) {
                phase = params.back();
                coreParams.pop_back();
            }

            if (coreParams.empty()) {
                return complex<double>(0.0, 0.0);
            }

            double k = coreParams[0];  // Coupling parameter always first
            
            switch (_FuncType) {
                case 0:  // a2(1320) - Breit-Wigner
                case 1:  // a2(1700) - Breit-Wigner  
                case 3:  // pi1(1400) - Breit-Wigner
                    return GetBreitWignerAmplitude(mass, k, coreParams)*exp(complex<double>(0, phase));
                    
                case 2:  // a0(980) - Flatté
                    return GetFlatteAmplitude(mass, k, coreParams)*exp(complex<double>(0, phase));
                    
                case 4:  // Conformal polynomial
                    return GetConformalPolynomialAmplitude(mass, coreParams)*exp(complex<double>(0, phase));
                    
                case 5:  // Flatté + conformal polynomial
                    return GetFlattePlusConformalAmplitude(mass, coreParams)*exp(complex<double>(0, phase));
                    
                default:
                    return complex<double>(0.0, 0.0);
            }
        }

        // Unified Breit-Wigner calculation for cases 0, 1, 3
        complex<double> GetBreitWignerAmplitude(double mass, double k, const vector<double>& params) const {
            // Get mass and width - either from params or defaults
            double M = (params.size() >= 2) ? params[1] : GetResonanceMass();
            double width = (params.size() >= 3) ? params[2] : GetResonanceWidth();
            
            double s = mass * mass;
            double M_sq = M * M;
            
            // Breit-Wigner: k / (M² - s - iMΓ)
            complex<double> denominator(M_sq - s, -M * width);
            return k / denominator;
        }

        // Flatté amplitude for case 2
        complex<double> GetFlatteAmplitude(double mass, double k, const vector<double>& params) const {
            double g1 = k;  // coupling to ηπ
            double g2 = (params.size() >= 2) ? params[1] : 0.340 * g1;  // coupling to KK̄
            double M0 = (params.size() >= 3) ? params[2] : GetResonanceMass();  // bare mass
            
            double s = mass * mass;
            double M0_sq = M0 * M0;
            
            // Phase space factors
            double rho_etapi = GetPhaseSpaceFactor(mass, 0.547853, 0.13957); // η-π
            double rho_KK = GetPhaseSpaceFactor(mass, 0.493677, 0.493677);   // K-K
            
            // Flatté denominator: M₀² - s - i(g₁²ρ₁ + g₂²ρ₂)
            complex<double> denominator(M0_sq - s, -(g1*g1*rho_etapi + g2*g2*rho_KK));
            return g1 / denominator;
        }
        
        // Phase space factor: ρ = 2k/√s
        double GetPhaseSpaceFactor(double mass, double m1, double m2) const {
            double s = mass * mass;
            double threshold_sq = (m1 + m2) * (m1 + m2);
            
            if (s <= threshold_sq) return 0.0;
            
            // Breakup momentum: k = √[(s-(m₁+m₂)²)(s-(m₁-m₂)²)]/(2√s)
            double diff_sq = (m1 - m2) * (m1 - m2);
            double k = TMath::Sqrt((s - threshold_sq) * (s - diff_sq)) / (2.0 * mass);
            
            return 2.0 * k / mass;
        }

        // Conformal polynomial amplitude for case 4
        complex<double> GetConformalPolynomialAmplitude(double mass, const vector<double>& params) const {
            // Parameters structure:
            // params[0] = re_0 (real part of z^0 coefficient)
            // params[1] = im_0 (imag part of z^0 coefficient)
            // params[2] = re_1 (real part of z^1 coefficient)
            // params[3] = im_1 (imag part of z^1 coefficient)
            // ... and so on for higher orders
            
            // Default conformal parameters
            double m_threshold = GetResonanceWidth();  // threshold mass (η + π⁰)
            double s_min = m_threshold * m_threshold;  // threshold in s
            double m_expansion = GetResonanceMass();   // expansion point mass (a0(980))
            double s_0 = m_expansion * m_expansion;    // expansion point in s
            
            // Compute conformal variable z
            double s = mass * mass;
            double sqrt_s = mass;
            double sqrt_s0_minus_smin = TMath::Sqrt(s_0 - s_min);
            double z = (sqrt_s - sqrt_s0_minus_smin) / (sqrt_s + sqrt_s0_minus_smin);
            
            // Build polynomial: sum_{i=0}^{N} (re_i + i*im_i) * z^i
            complex<double> amplitude(0.0, 0.0);
            
            if (params.size() < 2) {
                // Minimum: need at least re_0 and im_0
                return complex<double>(0.0, 0.0);
            }
            
            // Determine polynomial order from number of parameters
            // params layout: [re_0, im_0, re_1, im_1, re_2, im_2, ...]
            // For order N, we need 2*(N+1) parameters
            int order = params.size() / 2 - 1;
            
            double z_power = 1.0;
            for (int i = 0; i <= order; ++i) {
                size_t re_idx = 2*i;      // index for re_i
                size_t im_idx = 2*i + 1;  // index for im_i
                
                if (re_idx < params.size() && im_idx < params.size()) {
                    complex<double> coeff(params[re_idx], params[im_idx]);
                    amplitude += coeff * z_power;
                    z_power *= z;
                }
            }
            
            return amplitude;
        }

        // Flatté + conformal polynomial amplitude for case 5
        complex<double> GetFlattePlusConformalAmplitude(double mass, const vector<double>& params) const {
            // Parameters structure:
            // params[0] = scale (scale factor for Flatté amplitude)
            // params[1] = g_etapi (coupling to ηπ channel)
            // params[2] = g_KK (coupling to KK̄ channel)
            // params[3] = re_0 (real part of z^0 coefficient)
            // params[4] = im_0 (imag part of z^0 coefficient)
            // params[5] = re_1 (real part of z^1 coefficient)
            // params[6] = im_1 (imag part of z^1 coefficient)
            // ... and so on for higher orders
            
            if (params.size() < 5) {  // Minimum: scale, g_etapi, g_KK, re_0, im_0
                return complex<double>(0.0, 0.0);
            }
            
            // Extract Flatté parameters
            double scale = params[0];
            double g_etapi = params[1];
            double g_KK = params[2];
            double m0 = GetResonanceMass();  // Fixed at 0.98 GeV for a0(980)
            
            // Compute Flatté amplitude
            double s = mass * mass;
            double M0_sq = m0 * m0;
            
            // Phase space factors with proper analytic continuation below threshold
            double m_etapi_threshold = 0.683;  // η + π⁰ threshold in GeV
            double m_KK_threshold = 0.987354;  // K + K̄ threshold in GeV
            
            // η-π channel: m_η = 0.547853, m_π = 0.13957
            complex<double> rho_etapi = GetComplexPhaseSpaceFactor(mass, 0.547853, 0.13957);
            // K-K̄ channel: m_K = 0.493677
            complex<double> rho_KK = GetComplexPhaseSpaceFactor(mass, 0.493677, 0.493677);
            
            // Flatté denominator: M₀² - s - i*m₀*(g_ηπ²*ρ_ηπ + g_KK²*ρ_KK)
            complex<double> width_term = m0 * (g_etapi * g_etapi * rho_etapi + g_KK * g_KK * rho_KK);
            complex<double> denominator(M0_sq - s, 0.0);
            denominator -= complex<double>(0.0, 1.0) * width_term;
            
            complex<double> A_flatte = scale / denominator;
            
            // Compute conformal polynomial background
            double m_threshold = GetResonanceWidth();  // threshold mass (η + π⁰)
            double s_min = m_threshold * m_threshold;
            double m_expansion = GetResonanceMass();   // expansion point (a0(980))
            double s_0 = m_expansion * m_expansion;
            
            // Conformal variable z
            double sqrt_s = mass;
            double sqrt_s0_minus_smin = TMath::Sqrt(s_0 - s_min);
            double z = (sqrt_s - sqrt_s0_minus_smin) / (sqrt_s + sqrt_s0_minus_smin);
            
            // Build polynomial from remaining parameters
            complex<double> A_poly(0.0, 0.0);
            int poly_start = 3;  // Polynomial coefficients start at index 3
            int n_poly_params = params.size() - poly_start;
            
            if (n_poly_params >= 2) {
                int order = n_poly_params / 2 - 1;
                double z_power = 1.0;
                
                for (int i = 0; i <= order; ++i) {
                    size_t re_idx = poly_start + 2*i;
                    size_t im_idx = poly_start + 2*i + 1;
                    
                    if (re_idx < params.size() && im_idx < params.size()) {
                        complex<double> coeff(params[re_idx], params[im_idx]);
                        A_poly += coeff * z_power;
                        z_power *= z;
                    }
                }
            }
            
            // Total amplitude: A(s) = scale * A_Flatté(s) + A_poly(z)
            return A_flatte + A_poly;
        }
        
        // Complex phase space factor with analytic continuation below threshold
        complex<double> GetComplexPhaseSpaceFactor(double mass, double m1, double m2) const {
            double s = mass * mass;
            double threshold = m1 + m2;
            double threshold_sq = threshold * threshold;
            
            if (s >= threshold_sq) {
                // Above threshold: real and positive
                double diff_sq = (m1 - m2) * (m1 - m2);
                double k = TMath::Sqrt((s - threshold_sq) * (s - diff_sq)) / (2.0 * mass);
                double rho = 2.0 * k / mass;
                return complex<double>(rho, 0.0);
            } else {
                // Below threshold: imaginary (analytic continuation)
                double diff_sq = (m1 - m2) * (m1 - m2);
                double k_sq = (threshold_sq - s) * (s - diff_sq) / (4.0 * s);
                
                if (k_sq > 0) {
                    double k = TMath::Sqrt(k_sq);
                    double rho = 2.0 * k / mass;
                    return complex<double>(0.0, rho);  // Pure imaginary
                } else {
                    return complex<double>(0.0, 0.0);
                }
            }
        }
  };
}
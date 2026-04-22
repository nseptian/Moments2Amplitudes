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
            if (_FuncType == 3) {  // Polynomial
                return 2 + 2*(_PolyOrder + 1);  // m_threshold, m_expansion, (re_i, im_i) for i=0 to N
            } else if (_FuncType == 4) {  // Flatte + polynomial
                return 6 + 2*(_PolyOrder + 1);  // k, g_etapi, g_KK, m_threshold, m_expansion, Mass + polynomial coeffs
            } else if (_FuncType == 2) {  // Flatté
                return 4;
            } else if (_FuncType == 1) {  // Coherent sum of two Breit-Wigners
                return 6; // k1, M1, width1, k2, M2, width2 (global phase handled separately)
            } else {  // Breit-Wigner
                return 3;  // k, M, width
            }
        }

    private:
        static constexpr double kDefaultFlatteMass = 1.001;
        static constexpr double kPolynomialThresholdMass = 0.683;
        static constexpr double kPolynomialExpansionMass = 0.98;

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

            switch (_FuncType) {
                case 0:  // Breit-Wigner
                    {
                    double k = coreParams[0];
                    return GetBreitWignerAmplitude(mass, k, coreParams)*exp(complex<double>(0, phase));
                    }
                    
                case 2:  // a0(980) - Flatté
                    {
                    double k = coreParams[0];
                    return GetFlatteAmplitude(mass, k, coreParams)*exp(complex<double>(0, phase));
                    }
                    
                case 3:  // Polynomial
                    {
                        if (coreParams.size() < 4) {
                            return complex<double>(0.0, 0.0);
                        }
                        const double m_threshold = coreParams[0];
                        const double m_expansion = coreParams[1];
                        vector<double> polyParams(coreParams.begin() + 2, coreParams.end());
                        return GetConformalPolynomialAmplitude(mass, polyParams, m_threshold, m_expansion) * exp(complex<double>(0, phase));
                    }
                    
                case 4:  // Flatte + polynomial
                    {
                        const size_t requiredParams = static_cast<size_t>(GetNParameters());
                        if (coreParams.size() < requiredParams) {
                            return complex<double>(0.0, 0.0);
                        }

                        // Parameter layout for FuncType=4:
                        // [0]=k, [1]=g_etapi, [2]=g_KK, [3]=m_threshold, [4]=m_expansion, [5...end-1]=conformal polynomial coeffs, [end]=Mass
                        const double k = coreParams[0];
                        const double g_etapi = coreParams[1];
                        const double g_KK = coreParams[2];
                        const double m_threshold = coreParams[3];
                        const double m_expansion = coreParams[4];

                        const double s = mass * mass;
                        const double m0 = (coreParams.size() > 5) ? coreParams.back() : 0.98;  // Use passed Mass or default
                        const double m0_poly_sq = m_expansion * m_expansion;
                        const complex<double> rho_etapi = GetComplexPhaseSpaceFactor(mass, 0.547853, 0.13957);
                        const complex<double> rho_KK = GetComplexPhaseSpaceFactor(mass, 0.493677, 0.493677);
                        const complex<double> width_term = m0 * (g_etapi * g_etapi * rho_etapi + g_KK * g_KK * rho_KK);
                        complex<double> denominator(m0_poly_sq - s, 0.0);
                        denominator -= complex<double>(0.0, 1.0) * width_term;

                        complex<double> flatteAmp = 1.0 / denominator;
                        vector<double> polyParams(coreParams.begin() + 5, coreParams.end() - 1);  // Exclude Mass at end

                        complex<double> polyAmp = GetConformalPolynomialAmplitude(mass, polyParams, m_threshold, m_expansion);
                        complex<double> flatteTerm = k * flatteAmp * exp(complex<double>(0, phase));
                        return flatteTerm + polyAmp;
                    }

                case 1: // Coherent sum of two Breit-Wigners
                    {
                        if (coreParams.size() < 6) {
                            return complex<double>(0.0, 0.0);
                        }

                        const double k1 = coreParams[0];
                        const double M1 = coreParams[1];
                        const double w1 = coreParams[2];

                        const double k2 = coreParams[3];
                        const double M2 = coreParams[4];
                        const double w2 = coreParams[5];

                        const double s = mass * mass;
                        const complex<double> denom1(M1 * M1 - s, -M1 * w1);
                        const complex<double> denom2(M2 * M2 - s, -M2 * w2);

                        // a2(1320) phase is fixed to 0; optional global phase rotates only a2(1700).
                        complex<double> amp_1320 = k1 / denom1;
                        complex<double> amp_1700 = (k2 / denom2) * exp(complex<double>(0, phase));
                        return amp_1320 + amp_1700;
                    }
                    
                default:
                    return complex<double>(0.0, 0.0);
            }
        }

        // Unified Breit-Wigner calculation
        complex<double> GetBreitWignerAmplitude(double mass, double k, const vector<double>& params) const {
            // Require explicit mass and width from parameters.
            if (params.size() < 3) {
                return complex<double>(0.0, 0.0);
            }

            double M = params[1];
            double width = params[2];
            
            double s = mass * mass;
            double M_sq = M * M;
            
            // Breit-Wigner: k / (M² - s - iMΓ)
            complex<double> denominator(M_sq - s, -M * width);
            return k / denominator;
        }

        // Flatté amplitude for case 2
        complex<double> GetFlatteAmplitude(double mass, double k, const vector<double>& params) const {
            double g_etapi = params[1];
            double g_KK = params[2];
            double M0 = params[3];
            
            double s = mass * mass;
            double M0_sq = M0 * M0;
            
            // Phase space factors
            double rho_etapi = GetPhaseSpaceFactor(mass, 0.547853, 0.13957); // η-π
            double rho_KK = GetPhaseSpaceFactor(mass, 0.493677, 0.493677);   // K-K
            
            // Flatté denominator: M₀² - s - i(g₁²ρ₁ + g₂²ρ₂)
            complex<double> denominator(M0_sq - s, -(g_etapi*g_etapi*rho_etapi + g_KK*g_KK*rho_KK));
            return k / denominator;
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

        // Polynomial amplitude (conformal mapping basis)
        complex<double> GetConformalPolynomialAmplitude(
            double mass,
            const vector<double>& params,
            double thresholdMass = kPolynomialThresholdMass,
            double expansionMass = kPolynomialExpansionMass) const {
            // Parameters structure:
            // params[0] = re_0 (real part of z^0 coefficient)
            // params[1] = im_0 (imag part of z^0 coefficient)
            // params[2] = re_1 (real part of z^1 coefficient)
            // params[3] = im_1 (imag part of z^1 coefficient)
            // ... and so on for higher orders
            
            // Default conformal parameters
            double m_threshold = thresholdMass;  // threshold mass (η + π⁰)
            double s_min = m_threshold * m_threshold;  // threshold in s
            double m_expansion = expansionMass;   // expansion point mass
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
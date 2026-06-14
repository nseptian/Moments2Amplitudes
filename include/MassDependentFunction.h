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
        complex<double> GetAmplitude(double mass, const vector<double>& params, bool include_globalphase = false) const {
            return ComputeAmplitude(mass, params, include_globalphase);
        }

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
        // If _PolyOrder == -1 the polynomial block is disabled and contributes 0 parameters.
        int GetNParameters() const {
            int base = 0;
            switch (_FuncType) {
                case 0: // Polynomial (pure conformal polynomial)
                    base = 0;
                    break;
                case 1: // Flatté + polynomial
                    base = 5; // k,  Mass, g_etapi, g_KK, g_etaprimepi
                    break;
                case 2: // Breit-Wigner + polynomial
                    base = 3; // k, M, width
                    break;
                case 3: // Two Breit-Wigner + polynomial
                    base = 6; // k1, M1, width1, k2, M2, width2
                    break;
                case 4: // Flatte + Breit-Wigner + polynomial
                    base = 8; // k_fl, g_etapi, g_KK, g_etaprimepi, M_fl, k_bw, M_bw, width_bw
                    break;
            }

            int polyBlock = 0;
            if (_PolyOrder >= 0) {
                // polynomial block layout: m_threshold, m_expansion, then 2*(order+1) real/imag coefficients
                polyBlock = 2 + 2*(_PolyOrder + 1);
            }

            if (_FuncType == 0) {
                return polyBlock; // pure polynomial function
            }

            return base + polyBlock;
        }

    private:
        // static constexpr double kDefaultFlatteMass = 1.001;
        // static constexpr double kPolynomialThresholdMass = 0.683;
        // static constexpr double kPolynomialExpansionMass = 0.98;

        int _NMassBins = 0;
        double _FirstMassBinCenter = 0.0;
        double _MassBinWidth = 0.0;
        int _FuncType = 0;
        int _PolyOrder = 2;  // Polynomial order for conformal polynomial (used when ModelType indicates polynomial)

        // Core amplitude calculation
        complex<double> ComputeAmplitude(double mass, const vector<double>& params, bool include_globalphase) const {
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
            double phase_flatte = 0.0;
            double phase_bw = 0.0;
            if (include_globalphase && params.size() > 1) {
                if (_FuncType == 4 && params.size() > 2) {
                    // case-4 uses two independent phases appended in order: [phase_flatte, phase_bw]
                    phase_bw = coreParams.back();
                    coreParams.pop_back();
                    phase_flatte = coreParams.back();
                    coreParams.pop_back();
                } else {
                    phase = coreParams.back();
                    coreParams.pop_back();
                }
            }

            if (coreParams.empty()) {
                return complex<double>(0.0, 0.0);
            }
            // Determine whether a polynomial block is present in the parameters.
            const bool hasPoly = (_PolyOrder >= 0);
            const int polyBlockSize = hasPoly ? (2 + 2*(_PolyOrder + 1)) : 0; // threshold, expansion, coeffs

            // cout << "Computing amplitude for mass = " << mass << ", funcType = " << _FuncType 
            //      << ", hasPoly = " << hasPoly << ", polyOrder = " << _PolyOrder 
            //      << ", total params = " << params.size() << ", coreParams = " << coreParams.size() << endl;

            switch (_FuncType) {
                case 0:  // Polynomial (pure conformal polynomial)
                    {
                        // Pure polynomial function: expect [m_threshold, m_expansion, coeffs...]
                        if (coreParams.size() < 4) {
                            return complex<double>(0.0, 0.0);
                        }
                        const double m_threshold = coreParams[0];
                        const double m_expansion = coreParams[1];
                        vector<double> polyParams(coreParams.begin() + 2, coreParams.end());
                        return GetConformalPolynomialAmplitude(mass, polyParams, m_threshold, m_expansion) * exp(complex<double>(0, phase));
                    }

                case 1: // Flatté + optional polynomial
                    {
                        // Flatté base layout without polynomial: [k, Mass, g_etapi, g_KK, g_etaprimepi]
                        const int baseReq = 5; // k + tail (Mass,g_etapi,g_KK,g_etaprimepi)
                        // cout << "FuncType 1: baseReq = " << baseReq << ", coreParams size = " << coreParams.size() << endl;

                        if (coreParams.size() < static_cast<size_t>(baseReq)) return complex<double>(0.0, 0.0);

                        double k = coreParams[0];
                        vector<double> flatteParams(coreParams.begin(), coreParams.begin() + baseReq);
                        // cout << "Flatte params: ";
                        // for (size_t i = 0; i < flatteParams.size(); ++i) {
                        //     cout << flatteParams[i] << " ";
                        // }
                        // cout << endl;
                        complex<double> flatteAmp = GetFlatteAmplitude(mass, k, flatteParams) * exp(complex<double>(0, phase));
                        if (!hasPoly) return flatteAmp;

                        // If polynomial block present, it is expected appended after base params
                        if (coreParams.size() < static_cast<size_t>(baseReq + polyBlockSize)) return flatteAmp;
                        double m_threshold = coreParams[baseReq];
                        double m_expansion = coreParams[baseReq + 1];
                        vector<double> polyParams(coreParams.begin() + baseReq + 2, coreParams.begin() + baseReq + polyBlockSize);
                        complex<double> polyAmp = GetConformalPolynomialAmplitude(mass, polyParams, m_threshold, m_expansion);
                        return flatteAmp + polyAmp;
                    }

                case 2:  // Breit-Wigner + optional polynomial
                    {
                        // Base BW expects [k, M, width] then optional poly block appended
                        // cout << "FuncType 2: baseReq = " << 3 << ", coreParams size = " << coreParams.size() << endl;

                        const int baseReq = 3;
                        if (coreParams.size() < static_cast<size_t>(baseReq)) return complex<double>(0.0, 0.0);
                        double k = coreParams[0];
                        vector<double> bwParams(coreParams.begin(), coreParams.begin() + baseReq);
                        complex<double> bwAmp = GetBreitWignerAmplitude(mass, k, bwParams) * exp(complex<double>(0, phase));

                        if (!hasPoly) return bwAmp;

                        // Extract polynomial block appended after base params
                        if (coreParams.size() < static_cast<size_t>(baseReq + polyBlockSize)) return bwAmp;
                        double m_threshold = coreParams[baseReq];
                        double m_expansion = coreParams[baseReq + 1];
                        vector<double> polyParams(coreParams.begin() + baseReq + 2, coreParams.begin() + baseReq + polyBlockSize);
                        complex<double> polyAmp = GetConformalPolynomialAmplitude(mass, polyParams, m_threshold, m_expansion);
                        return bwAmp + polyAmp;
                    }

                case 3: // Coherent sum of two Breit-Wigners + optional polynomial
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

                        // resonance 1 phase is fixed to 0; optional global phase rotates only resonance 2.
                        complex<double> amp_1 = k1 / denom1;
                        complex<double> amp_2 = (k2 / denom2) * exp(complex<double>(0, phase));
                        complex<double> baseAmp = amp_1 + amp_2;

                        if (!hasPoly) return baseAmp;

                        // If polynomial block present, it is expected appended after base params
                        const size_t baseReq = 6;
                        if (coreParams.size() < static_cast<size_t>(baseReq + polyBlockSize)) return baseAmp;
                        double m_threshold = coreParams[baseReq];
                        double m_expansion = coreParams[baseReq + 1];
                        vector<double> polyParams(coreParams.begin() + baseReq + 2, coreParams.begin() + baseReq + polyBlockSize);
                        complex<double> polyAmp = GetConformalPolynomialAmplitude(mass, polyParams, m_threshold, m_expansion);
                        return baseAmp + polyAmp;
                    }

                case 4: // Coherent sum of Flatté + Breit-Wigner + optional polynomial
                    {
                        // Flatté part expects [k, g_etapi, g_KK, g_etaprimepi, Mass], BW part expects [k, M, width], then optional poly block appended
                        const int flatteReq = 5;
                        const int bwReq = 3;
                        if (coreParams.size() < static_cast<size_t>(flatteReq + bwReq)) return complex<double>(0.0, 0.0);

                        double k_fl = coreParams[0];
                        vector<double> flatteParams(coreParams.begin(), coreParams.begin() + flatteReq);
                        complex<double> flatteAmp = GetFlatteAmplitude(mass, k_fl, flatteParams) * exp(complex<double>(0, phase_flatte));

                        double k_bw = coreParams[flatteReq];
                        vector<double> bwParams(coreParams.begin() + flatteReq, coreParams.begin() + flatteReq + bwReq);
                        complex<double> bwAmp = GetBreitWignerAmplitude(mass, k_bw, bwParams) * exp(complex<double>(0, phase_bw));

                        complex<double> baseAmp = flatteAmp + bwAmp;

                        // If polynomial block present, it is expected appended after base params
                        const size_t baseReq = flatteReq + bwReq;
                        if (coreParams.size() < static_cast<size_t>(baseReq + polyBlockSize)) return baseAmp;
                        double m_threshold = coreParams[baseReq];
                        double m_expansion = coreParams[baseReq + 1];
                        vector<double> polyParams(coreParams.begin() + baseReq + 2, coreParams.begin() + baseReq + polyBlockSize);
                        complex<double> polyAmp = GetConformalPolynomialAmplitude(mass, polyParams, m_threshold, m_expansion);
                        return baseAmp + polyAmp;
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
            if (params.size() < 5) {
                return complex<double>(0.0, 0.0);
            }

            double M0 = params[1];
            double g_etapi = params[2];
            double g_KK = params[3];
            double g_etaprimepi = params[4];

            double g_etapi_sqr = g_etapi * g_etapi;
            double g_KK_sqr = g_KK * g_KK;
            double g_etaprimepi_sqr = g_etaprimepi * g_etaprimepi;
            if (g_KK < 0.0) {
                const double g_etapi_per_g_KK_sqr = 1.05; // (g_etapi/g_KK)² from CBAR
                g_KK_sqr = g_etapi_sqr * g_etapi_per_g_KK_sqr;
            }
            if (g_etaprimepi < 0.0) {
                const double g_etaprimepi_per_g_etapi_sqr = 0.772; // (g_etaprimepi/g_etapi)² from CBAR
                g_etaprimepi_sqr = g_etapi_sqr * g_etaprimepi_per_g_etapi_sqr;
            }
            
            double s = mass * mass;
            double M0_sq = M0 * M0;
            
            // Phase space factors
            double rho_etapi = GetPhaseSpaceFactor(mass, 0.547853, 0.13957); // η-π
            double rho_KK = GetPhaseSpaceFactor(mass, 0.493677, 0.493677);   // K-K
            double rho_etaprimepi = GetPhaseSpaceFactor(mass, 0.95778, 0.13957); // η'-π
            
            // Flatté denominator: M₀² - s - i(g₁²ρ₁ + g₂²ρ₂)
            complex<double> denominator(M0_sq - s, -(g_etapi_sqr*rho_etapi + g_KK_sqr*rho_KK + g_etaprimepi_sqr*rho_etaprimepi));
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
            double thresholdMass,
            double expansionMass) const {
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
            
            // Compute canonical conformal variable z using threshold-shifted sqrt(s - s_min)
            double s = mass * mass;
            double sqrt_s_minus_smin = 0.0;
            if (s > s_min) {
                sqrt_s_minus_smin = TMath::Sqrt(s - s_min);
            }
            double sqrt_s0_minus_smin = TMath::Sqrt(s_0 - s_min);
            double z = (sqrt_s_minus_smin - sqrt_s0_minus_smin) /
                       (sqrt_s_minus_smin + sqrt_s0_minus_smin);
            
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
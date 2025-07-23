#pragma once

#include "TMath.h"
#include <complex.h>

namespace m2pw{

  class MassDependentFunction{
    public:
        MassDependentFunction() = default;
        MassDependentFunction(const MassDependentFunction&) = default;
        MassDependentFunction(MassDependentFunction&&) = default;
        MassDependentFunction& operator=(const MassDependentFunction&) = default;
        MassDependentFunction& operator=(MassDependentFunction&&) = default;
        ~MassDependentFunction() = default;

        MassDependentFunction(int NMassBins, double FirstMassBinCenter, double MassBinWidth, int FuncType) : 
            _NMassBins(NMassBins),
            _FirstMassBinCenter(FirstMassBinCenter),
            _MassBinWidth(MassBinWidth),
            _FuncType(FuncType) {}

        double GetPWMagnitude(double mass, const double params) const {
            int bin = static_cast<int>((mass - _FirstMassBinCenter) / _MassBinWidth);
            if (bin < 0 || bin >= _NMassBins) {
                return 0.0; // Out of bounds
            }
            if (_FuncType == 0) { // Breit-Wigner for a2(1320)
                double k = params;
                double M = 1.3182;
                double width = 0.107;
                double a = mass * mass - M * M;
                double b = M * width;
                double Re_BW = a * k / (a * a + b * b);
                double Im_BW = -b * k / (a * a + b * b);
                double magnitude = TMath::Sqrt(Re_BW * Re_BW + Im_BW * Im_BW);

                return magnitude;
            } else if (_FuncType == 1) { // Breit-Wigner for a2(1700)
                double k = params;
                double M = 1.706;
                double width = 0.380;
                double a = mass * mass - M * M;
                double b = M * width;
                double Re_BW = a * k / (a * a + b * b);
                double Im_BW = -b * k / (a * a + b * b);
                double magnitude = TMath::Sqrt(Re_BW * Re_BW + Im_BW * Im_BW);

                return magnitude;
             } else if (_FuncType == 2) { // Flatté for a0(980) in ηπ
                double g1 = params; // coupling to ηπ
                double g2 = 0.340 * g1; // coupling to KK̄ (typical ratio for a0(980))
                double M0 = 1.001; // bare mass in GeV
                
                complex<double> amp = GetFlatteAmplitude(mass, M0, g1, g2);
                return abs(amp);
            }
            return 0.0; // Default case
        }

        double GetPWMagnitude(double mass, const vector<double> params) const {
            int bin = static_cast<int>((mass - _FirstMassBinCenter) / _MassBinWidth);
            if (bin < 0 || bin >= _NMassBins) {
                return 0.0; // Out of bounds
            }
            if (_FuncType == 0) { // Breit-Wigner for a2(1320)
                double k = params[0];
                double M = params[1];
                double width = params[2];
                if (params.size()==3) {
                    M = params[1];
                    width = params[2];
                }
                double a = mass * mass - M * M;
                double b = M * width;
                double Re_BW = a * k / (a * a + b * b);
                double Im_BW = -b * k / (a * a + b * b);
                double magnitude = TMath::Sqrt(Re_BW * Re_BW + Im_BW * Im_BW);

                return magnitude;
            } else if (_FuncType == 1) { // Breit-Wigner for a2(1700)
                double k = params[0];
                double M = 1.706;
                double width = 0.380;
                if (params.size()==3) {
                    M = params[1];
                    width = params[2];
                }
                double a = mass * mass - M * M;
                double b = M * width;
                double Re_BW = a * k / (a * a + b * b);
                double Im_BW = -b * k / (a * a + b * b);
                double magnitude = TMath::Sqrt(Re_BW * Re_BW + Im_BW * Im_BW);

                return magnitude;
            } else if (_FuncType == 2) { // Flatté for a0(980) in ηπ
                double g1 = params[0]; // coupling to ηπ
                double g2 = params.size() > 1 ? params[1] : 0.340 * g1; // coupling to KK̄
                double M0 = params.size() > 2 ? params[2] : 1.001; // bare mass
                
                complex<double> amp = GetFlatteAmplitude(mass, M0, g1, g2);
                return abs(amp);
            }
            return 0.0; // Default case
        }

        double GetPWPhase(double mass, const double params) const {
            int bin = static_cast<int>((mass - _FirstMassBinCenter) / _MassBinWidth);
            if (bin < 0 || bin >= _NMassBins) {
                return 0.0; // Out of bounds
            }
            if (_FuncType == 0) { // Breit-Wigner for a2(1320)
                double k = params;
                double M = 1.3182;
                double width = 0.107;
                double a = mass * mass - M * M;
                double b = M * width;
                double Re_BW = a * k / (a * a + b * b);
                double Im_BW = -b * k / (a * a + b * b);
                double phase = TMath::ATan2(Im_BW, Re_BW);
                return phase;
            } else if (_FuncType == 1) { // Breit-Wigner for a2(1700)
                double k = params;
                double M = 1.706;
                double width = 0.380;
                double a = mass * mass - M * M;
                double b = M * width;
                double Re_BW = a * k / (a * a + b * b);
                double Im_BW = -b * k / (a * a + b * b);
                double phase = TMath::ATan2(Im_BW, Re_BW);
                return phase;
            } else if (_FuncType == 2) { // Flatté for a0(980) in ηπ
                double g1 = params;
                double g2 = 0.340 * g1;
                double M0 = 1.001;
                
                complex<double> amp = GetFlatteAmplitude(mass, M0, g1, g2);
                return arg(amp);
            }
            return 0.0; // Default case
        }

        double GetPWPhase(double mass, const vector<double> params) const {
            int bin = static_cast<int>((mass - _FirstMassBinCenter) / _MassBinWidth);
            if (bin < 0 || bin >= _NMassBins) {
                return 0.0; // Out of bounds
            }
            if (_FuncType == 0) { // Breit-Wigner for a2(1320)
                double k = params[0];
                double M = params[1];
                double width = params[2];
                if (params.size()==3) {
                    M = params[1];
                    width = params[2];
                }
                double a = mass * mass - M * M;
                double b = M * width;
                double Re_BW = a * k / (a * a + b * b);
                double Im_BW = -b * k / (a * a + b * b);
                double phase = TMath::ATan2(Im_BW, Re_BW);
                return phase;
            } else if (_FuncType == 1) { // Breit-Wigner for a2(1700)
                double k = params[0];
                double M = 1.706;
                double width = 0.380;
                if (params.size()==3) {
                    M = params[1];
                    width = params[2];
                }
                double a = mass * mass - M * M;
                double b = M * width;
                double Re_BW = a * k / (a * a + b * b);
                double Im_BW = -b * k / (a * a + b * b);
                double phase = TMath::ATan2(Im_BW, Re_BW);
                return phase;
            } else if (_FuncType == 2) { // Flatté for a0(980) in ηπ
                double g1 = params[0];
                double g2 = params.size() > 1 ? params[1] : 0.340 * g1;
                double M0 = params.size() > 2 ? params[2] : 1.001;
                
                complex<double> amp = GetFlatteAmplitude(mass, M0, g1, g2);
                return arg(amp);
            }
            return 0.0; // Default case
        }

        int GetNMassBins() const { return _NMassBins; }
        double GetFirstMassBinCenter() const { return _FirstMassBinCenter; }
        double GetMassBinWidth() const { return _MassBinWidth; }
        int GetFuncType() const { return _FuncType; }

    private:
        int _NMassBins = 0;
        double _FirstMassBinCenter = 0.0;
        double _MassBinWidth = 0.0;
        int _FuncType = 0; // 0 for Breit-Wigner, 1 for kMatrix, etc

        // Flatté amplitude calculation for a0(980)
        complex<double> GetFlatteAmplitude(double mass, double M0, double g1, double g2) const {
            double s = mass * mass;
            double M0_sq = M0 * M0;
            
            // Phase space factors
            double rho_etapi = GetPhaseSpaceFactor(mass, 0.547853, 0.13957); // η=547.853 MeV, π=139.57 MeV
            double rho_KK = GetPhaseSpaceFactor(mass, 0.493677, 0.493677);  // K±=493.677 MeV
            
            // Flatté denominator: M₀² - s - i(g₁²ρ₁ + g₂²ρ₂)
            complex<double> denominator = complex<double>(M0_sq - s, -(g1*g1*rho_etapi + g2*g2*rho_KK));
            
            // Return g₁/denominator for ηπ channel
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
  };
}
#pragma once

#include "MomentHelper.h"
#include <map>
#include <stdexcept>

class MassDependentMoments {
    public:
        MassDependentMoments() = default;

        MassDependentMoments(Double_t NMassBins, Double_t MassBinWidth, Double_t FirstMassBinCenter,
                          const TString& FitResultsDir, const TString& FitResultsFilename, const TString& MassBinName)
            : _NMassBins(NMassBins), _MassBinWidth(MassBinWidth), _FirstMassBinCenter(FirstMassBinCenter),
              _FitResultsDir(FitResultsDir), _FitResultsFilename(FitResultsFilename), _MassBinName(MassBinName) {

            for (Double_t i = 0; i < _NMassBins; ++i) {
                Double_t massValue = _FirstMassBinCenter + i * _MassBinWidth;
                _MassBins.push_back(massValue);

                const TString filePath = Form("%s/%s%1.6f_/", _FitResultsDir.Data(), _MassBinName.Data(), massValue);
                std::cout << "Loading moments from: " << filePath << std::endl;
                MomentHelper moments;
                moments.Set(filePath, 1, 1);
                _moments[massValue] = moments;
            }
        }

        MassDependentMoments(Double_t NMassBins, Double_t MassBinWidth, Double_t FirstMassBinCenter,
                          const TString& FitResultsDir, const TString& FitResultsFilename, const TString& MassBinName, const TString& VarBinName)
            : _NMassBins(NMassBins), _MassBinWidth(MassBinWidth), _FirstMassBinCenter(FirstMassBinCenter),
              _FitResultsDir(FitResultsDir), _FitResultsFilename(FitResultsFilename), _MassBinName(MassBinName), _VarBinName(VarBinName) {

            for (Double_t i = 0; i < _NMassBins; ++i) {
                Double_t massValue = _FirstMassBinCenter + i * _MassBinWidth;
                _MassBins.push_back(massValue);

                const TString filePath = Form("%s/%s%1.6f_%s_/", _FitResultsDir.Data(), _MassBinName.Data(), massValue, _VarBinName.Data());
                std::cout << "Loading moments from: " << filePath << std::endl;
                MomentHelper moments;
                moments.Set(filePath, 1, 1);
                _moments[massValue] = moments;
            }
        }

        // Function to set the mass-dependent moments
        void SetMoments(const Double_t massValue, MomentHelper& moments) {
            _moments[massValue] = moments;
        }

        // Function to get the mass-dependent moments
        MomentHelper GetMoments(const Double_t massValue) const {
            auto it = _moments.find(massValue);
            if (it != _moments.end()) {
                return it->second;
            } else {
                throw std::runtime_error("Mass value not found in moments map.");
            }
        }

        void PrintMoments() {
            std::cout << "Mass-dependent moments size: " << _moments.size() << std::endl;
            for (auto it = _moments.begin(); it != _moments.end(); ++it) {
                std::cout << std::endl;
                std::cout << "Mass: " << it->first << std::endl;
                it->second.PrintVals();
            }
        }
        
        size_t GetNMassBins() {
            return _moments.size();
        }

        const std::vector<Double_t>& GetMassBins() const {
            return _MassBins;
        }

    private:
        std::map<Double_t, MomentHelper> _moments;
        Double_t _NMassBins;
        Double_t _MassBinWidth;
        Double_t _FirstMassBinCenter;
        TString _FitResultsDir;
        TString _FitResultsFilename;
        TString _MassBinName;
        TString _VarBinName;
        std::vector<Double_t> _MassBins;
};
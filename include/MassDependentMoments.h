#pragma once

#include "MomentHelper.h"
#include <map>
#include <stdexcept>
#include <TGraphErrors.h>
#include <TKey.h>
#include <cmath>
#include <algorithm>

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
                const TString fullFilePath = filePath + _FitResultsFilename;
                
                // Check if file exists
                TFile* testFile = TFile::Open(fullFilePath);
                if (!testFile || testFile->IsZombie()) {
                    if (testFile) delete testFile;
                    throw std::runtime_error(Form("File not found: %s", fullFilePath.Data()));
                }
                delete testFile;
                
                std::cout << "Loading moments from: " << filePath << std::endl;
                MomentHelper moments;
                moments.Set(filePath, _FitResultsFilename, 1, 1);
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
                const TString fullFilePath = filePath + _FitResultsFilename;
                
                // Check if file exists
                TFile* testFile = TFile::Open(fullFilePath);
                if (!testFile || testFile->IsZombie()) {
                    if (testFile) delete testFile;
                    throw std::runtime_error(Form("File not found: %s", fullFilePath.Data()));
                }
                delete testFile;
                
                std::cout << "Loading moments from: " << filePath << std::endl;
                MomentHelper moments;
                moments.Set(filePath, _FitResultsFilename, 1, 1);
                _moments[massValue] = moments;
            }
        }

        /**
         * @brief Load moments directly from a list of explicit mass bins
         */
        MassDependentMoments(const std::vector<Double_t>& massBins,
                             const TString& FitResultsDir,
                             const TString& FitResultsFilename,
                             const TString& MassBinName,
                             const TString& VarBinName = "")
            : _NMassBins(massBins.size()),
              _FitResultsDir(FitResultsDir),
              _FitResultsFilename(FitResultsFilename),
              _MassBinName(MassBinName),
              _VarBinName(VarBinName),
              _MassBins(massBins) {

            if (_MassBins.size() >= 2) {
                _MassBinWidth = _MassBins[1] - _MassBins[0];
                _FirstMassBinCenter = _MassBins[0];
            } else if (!_MassBins.empty()) {
                _MassBinWidth = 0.0;
                _FirstMassBinCenter = _MassBins[0];
            }

            for (const Double_t massValue : _MassBins) {
                TString filePath;
                if (_VarBinName.IsNull() || _VarBinName.EqualTo("")) {
                    filePath = Form("%s/%s%1.6f_/", _FitResultsDir.Data(), _MassBinName.Data(), massValue);
                } else {
                    filePath = Form("%s/%s%1.6f_%s_/", _FitResultsDir.Data(), _MassBinName.Data(), massValue, _VarBinName.Data());
                }
                const TString fullFilePath = filePath + _FitResultsFilename;

                TFile* testFile = TFile::Open(fullFilePath);
                if (!testFile || testFile->IsZombie()) {
                    if (testFile) delete testFile;
                    throw std::runtime_error(Form("File not found: %s", fullFilePath.Data()));
                }
                delete testFile;

                std::cout << "Loading moments from: " << filePath << std::endl;
                MomentHelper moments;
                moments.Set(filePath, _FitResultsFilename, 1, 1);
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
            }
            for (const auto& [mass, moms] : _moments) {
                if (std::abs(mass - massValue) < 1e-5) {
                    return moms;
                }
            }
            throw std::runtime_error(Form("Mass value %g not found in moments map.", massValue));
        }

        void PrintMoments() {
            std::cout << "Mass-dependent moments size: " << _moments.size() << std::endl;
            for (auto it = _moments.begin(); it != _moments.end(); ++it) {
                std::cout << std::endl;
                std::cout << "Mass: " << it->first << std::endl;
                it->second.PrintVals();
            }
        }

        /**
         * @brief Print the estimated covariance matrix for a specific mass bin
         * @param massValue Center of the mass bin in GeV
         */
        void PrintCovarianceMatrix(const Double_t massValue) const {
            auto it = _moments.find(massValue);
            if (it != _moments.end()) {
                std::cout << "\n=== Estimated Covariance Matrix for Mass Bin: " << massValue << " GeV ===" << std::endl;
                it->second.PrintCovarianceMatrix();
                return;
            }
            for (const auto& [mass, moms] : _moments) {
                if (std::abs(mass - massValue) < 1e-4) {
                    std::cout << "\n=== Estimated Covariance Matrix for Mass Bin: " << mass << " GeV ===" << std::endl;
                    moms.PrintCovarianceMatrix();
                    return;
                }
            }
            std::cerr << "Mass bin " << massValue << " not found." << std::endl;
        }

        /**
         * @brief Print the estimated covariance matrix for every mass bin
         */
        void PrintCovarianceMatrices() const {
            std::cout << "\n========================================================" << std::endl;
            std::cout << " Estimated Covariance Matrices for All " << _moments.size() << " Mass Bins" << std::endl;
            std::cout << "========================================================" << std::endl;
            for (const auto& [mass, moms] : _moments) {
                std::cout << "\n--- Mass Bin: " << mass << " GeV ---" << std::endl;
                moms.PrintCovarianceMatrix();
            }
        }
        
        size_t GetNMassBins() {
            return _moments.size();
        }

        const std::vector<Double_t>& GetMassBins() const {
            return _MassBins;
        }

        /**
         * @brief Check whether MCMC samples are available across loaded mass bins
         */
        bool HasMCMCSamples() const {
            return GetMCMCEntries() > 0;
        }

        /**
         * @brief Get minimum number of MCMC entries across all mass bins
         */
        Long64_t GetMCMCEntries() const {
            if (_moments.empty()) return 0;
            Long64_t minEntries = -1;
            for (const auto& [mass, momHelper] : _moments) {
                Long64_t n = momHelper.GetMCMCEntries();
                if (minEntries < 0 || n < minEntries) {
                    minEntries = n;
                }
            }
            return (minEntries < 0) ? 0 : minEntries;
        }

        /**
         * @brief Load moment values for a specific MCMC entry into all mass bins
         * @param entryIndex MCMC tree entry index
         */
        void LoadMCMCEntry(Long64_t entryIndex) {
            for (auto& [mass, momHelper] : _moments) {
                momHelper.LoadMCMCEntry(entryIndex);
            }
        }

        /**
         * @brief Get sample weight for a specific MCMC entry
         * @param entryIndex MCMC tree entry index
         * @return Weight value, or 1.0 if not available / out of range
         */
        Double_t GetMCMCWeight(Long64_t entryIndex) const {
            if (_moments.empty()) return 1.0;
            return _moments.begin()->second.GetMCMCWeight(entryIndex);
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
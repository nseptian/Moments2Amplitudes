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
                moments.Set(filePath, 1, 1);
                _moments[massValue] = moments;
            }
        }

        // =========================================================================
        // OPTION 2 (Commented out for now - working on MCMCTree first):
        // Load unnormalized moments as a function of mass from a ROOT file containing TGraphErrors
        // (e.g., UnnormalizedMoments_vs_Mass.root).
        // =========================================================================
        /*
        MassDependentMoments(const TString& momentsVsMassFilePath,
                             const TString& FitResultsDir = "",
                             const TString& FitResultsFilename = "ResultsBruMcmcFixedCovariance.root",
                             const TString& MassBinName = "Mpi0eta",
                             const TString& VarBinName = "")
            : _FitResultsDir(FitResultsDir), _FitResultsFilename(FitResultsFilename),
              _MassBinName(MassBinName), _VarBinName(VarBinName) {
            
            TFile* file = TFile::Open(momentsVsMassFilePath);
            if (!file || file->IsZombie()) {
                if (file) delete file;
                throw std::runtime_error(Form("Cannot open moments file: %s", momentsVsMassFilePath.Data()));
            }

            // Map: mass bin center -> (moment name -> pair(val, error))
            std::map<Double_t, std::map<TString, std::pair<Double_t, Double_t>>> massBinMoments;

            TIter next(file->GetListOfKeys());
            TKey* key = nullptr;
            while ((key = (TKey*)next())) {
                TString keyName = key->GetName();
                if (!keyName.BeginsWith("gr_unnormalized_moments_") || !keyName.EndsWith("_vs_mass")) {
                    continue;
                }
                TString middle = keyName;
                middle.ReplaceAll("gr_unnormalized_moments_", "");
                middle.ReplaceAll("_vs_mass", "");
                TString momentName = "H_" + middle;

                TGraphErrors* gr = dynamic_cast<TGraphErrors*>(key->ReadObj());
                if (!gr) continue;

                const int nPoints = gr->GetN();
                for (int p = 0; p < nPoints; ++p) {
                    Double_t mass = 0.0;
                    Double_t val = 0.0;
                    gr->GetPoint(p, mass, val);
                    Double_t err = gr->GetErrorY(p);
                    Double_t roundedMass = std::round(mass * 1e6) / 1e6;
                    massBinMoments[roundedMass][momentName] = {val, err};
                }
                delete gr;
            }
            delete file;

            if (massBinMoments.empty()) {
                throw std::runtime_error(Form("No unnormalized moments found in: %s", momentsVsMassFilePath.Data()));
            }

            _NMassBins = massBinMoments.size();
            for (const auto& [mass, momMap] : massBinMoments) {
                _MassBins.push_back(mass);
                MomentHelper moments;
                for (const auto& [momName, valErr] : momMap) {
                    moments.SetUnnormalized(momName, valErr.first, valErr.second);
                }

                bool loadedCov = false;
                if (!_FitResultsDir.IsNull() && !_FitResultsDir.EqualTo("")) {
                    TString filePath;
                    if (_VarBinName.IsNull() || _VarBinName.EqualTo("")) {
                        filePath = Form("%s/%s%1.6f_/", _FitResultsDir.Data(), _MassBinName.Data(), mass);
                    } else {
                        filePath = Form("%s/%s%1.6f_%s_/", _FitResultsDir.Data(), _MassBinName.Data(), mass, _VarBinName.Data());
                    }
                    const TString fullFilePath = filePath + _FitResultsFilename;

                    TFile* testFile = TFile::Open(fullFilePath);
                    if (testFile && !testFile->IsZombie()) {
                        if (testFile->Get("PosteriorCovariance")) {
                            delete testFile;
                            moments.SetFilePath(filePath);
                            loadedCov = moments.LoadCovarianceMatrix(fullFilePath, "PosteriorCovariance");
                        } else if (testFile->Get("MCMCTree")) {
                            delete testFile;
                            moments.Set(filePath, 1, 1);
                            // Re-apply graph values and errors to ensure they are preserved
                            for (const auto& [momName, valErr] : momMap) {
                                moments.SetUnnormalized(momName, valErr.first, valErr.second);
                            }
                            loadedCov = moments.HasCovarianceMatrix();
                        } else {
                            delete testFile;
                        }
                    } else if (testFile) {
                        delete testFile;
                    }
                }

                if (!loadedCov) {
                    moments.SetDiagonalCovariance();
                }

                _moments[mass] = moments;
            }

            if (_MassBins.size() >= 2) {
                _MassBinWidth = _MassBins[1] - _MassBins[0];
                _FirstMassBinCenter = _MassBins[0];
            } else if (!_MassBins.empty()) {
                _MassBinWidth = 0.0;
                _FirstMassBinCenter = _MassBins[0];
            }
        }
        */

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
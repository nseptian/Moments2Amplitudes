#pragma once

#include "CrossSection.h"
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooArgList.h>
#include <TString.h>
#include <TFile.h>
#include <TMatrixDSym.h>
#include <TSystem.h>
#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <memory>

/**
 * @brief Class to load and manage moments from fit results
 * 
 * This class provides functionality to load moments from BruFit results,
 * including both MCMC and direct fit results. It supports both normalized
 * and unnormalized moments with error calculations.
 */
class MomentHelper {
public:
  
  MomentHelper() = default;

  // Copy constructor - safely resizes ROOT TMatrixDSym
  MomentHelper(const MomentHelper& other)
    : _filepath(other._filepath),
      _moments(other._moments),
      _moments_err(other._moments_err),
      _unnormalized_moments(other._unnormalized_moments),
      _unnormalized_moments_err(other._unnormalized_moments_err),
      _cov_moment_names(other._cov_moment_names),
      _cov_name_to_index(other._cov_name_to_index),
      _has_cov(other._has_cov) {
    if (other._has_cov && other._unnormalized_cov.GetNrows() > 0) {
      _unnormalized_cov.ResizeTo(other._unnormalized_cov.GetNrows(), other._unnormalized_cov.GetNcols());
      _unnormalized_cov = other._unnormalized_cov;
    }
  }

  // Copy assignment operator - safely resizes ROOT TMatrixDSym before assignment
  MomentHelper& operator=(const MomentHelper& other) {
    if (this != &other) {
      _filepath = other._filepath;
      _moments = other._moments;
      _moments_err = other._moments_err;
      _unnormalized_moments = other._unnormalized_moments;
      _unnormalized_moments_err = other._unnormalized_moments_err;
      _cov_moment_names = other._cov_moment_names;
      _cov_name_to_index = other._cov_name_to_index;
      _has_cov = other._has_cov;
      if (other._has_cov && other._unnormalized_cov.GetNrows() > 0) {
        _unnormalized_cov.ResizeTo(other._unnormalized_cov.GetNrows(), other._unnormalized_cov.GetNcols());
        _unnormalized_cov = other._unnormalized_cov;
      } else {
        _unnormalized_cov.Clear();
      }
    }
    return *this;
  }

  // Move constructor
  MomentHelper(MomentHelper&& other) noexcept
    : _filepath(std::move(other._filepath)),
      _moments(std::move(other._moments)),
      _moments_err(std::move(other._moments_err)),
      _unnormalized_moments(std::move(other._unnormalized_moments)),
      _unnormalized_moments_err(std::move(other._unnormalized_moments_err)),
      _cov_moment_names(std::move(other._cov_moment_names)),
      _cov_name_to_index(std::move(other._cov_name_to_index)),
      _has_cov(other._has_cov) {
    if (other._has_cov && other._unnormalized_cov.GetNrows() > 0) {
      _unnormalized_cov.ResizeTo(other._unnormalized_cov.GetNrows(), other._unnormalized_cov.GetNcols());
      _unnormalized_cov = other._unnormalized_cov;
    }
  }

  // Move assignment operator
  MomentHelper& operator=(MomentHelper&& other) noexcept {
    if (this != &other) {
      _filepath = std::move(other._filepath);
      _moments = std::move(other._moments);
      _moments_err = std::move(other._moments_err);
      _unnormalized_moments = std::move(other._unnormalized_moments);
      _unnormalized_moments_err = std::move(other._unnormalized_moments_err);
      _cov_moment_names = std::move(other._cov_moment_names);
      _cov_name_to_index = std::move(other._cov_name_to_index);
      _has_cov = other._has_cov;
      if (other._has_cov && other._unnormalized_cov.GetNrows() > 0) {
        _unnormalized_cov.ResizeTo(other._unnormalized_cov.GetNrows(), other._unnormalized_cov.GetNcols());
        _unnormalized_cov = other._unnormalized_cov;
      } else {
        _unnormalized_cov.Clear();
      }
    }
    return *this;
  }

  /**
   * @brief Set a moment value directly
   * @param mom Moment name
   * @param val Moment value
   */
  void Set(const TString& mom, Double_t val) {
    _moments[mom] = val;
  }

  /**
   * @brief Load moments from a BruFit result file
   * @param filename Path to the result file
   * @param BruFitResult Flag to load BruFit results (default: 1)
   * @param MCMCResult Flag to load MCMC results (default: 1)
   */
  void Set(const TString& filepath, const Int_t BruFitResult = 1, const Int_t MCMCResult = 1){
    _filepath = filepath;
    if (!BruFitResult) return;
    
    if (MCMCResult) {
      LoadMCMCResults(_filepath + "/ResultsBruMcmcCovariance.root");
    } else {
      LoadDirectResults(_filepath + "/ResultsBruFit.root");
    }
  }

  /**
   * @brief Set the directory path for result files
   * @param filepath Path to directory containing result files
   */
  void SetFilePath(const TString& filepath) {
    _filepath = filepath;
  }

  /**
   * @brief Get the directory path for result files
   * @return Directory path string
   */
  const TString& GetFilePath() const {
    return _filepath;
  }

  /**
   * @brief Set moments from a RooArgSet
   * @param pars RooArgSet containing moment parameters
   */
  void Set(const RooArgSet* pars) {
    if (!pars) return;
    
    for (auto* mom : *pars) {
      RooRealVar* rmom = dynamic_cast<RooRealVar*>(mom);
      if (!rmom) continue;
      if (mom->GetName()[0] != 'H') continue;
      _moments[TString(mom->GetName())] = rmom->getVal();
    }
  }

  /**
   * @brief Get normalized moment value
   * @param name Moment name
   * @return Moment value, or 0.0 if not found
   */
  Double_t GetVal(const TString& name) const {
    auto it = _moments.find(name);
    return (it != _moments.end()) ? it->second : 0.0;
  }

  /**
   * @brief Get unnormalized moment value
   * @param name Moment name
   * @return Unnormalized moment value, or 0.0 if not found
   */
  Double_t GetUnnormalizedVal(const TString& name) const {
    auto it = _unnormalized_moments.find(name);
    return (it != _unnormalized_moments.end()) ? it->second : 0.0;
  }

  /**
   * @brief Get unnormalized moment error
   * @param name Moment name
   * @return Unnormalized moment error, or 0.0 if not found
   */
  Double_t GetUnnormalizedError(const TString& name) const {
    auto it = _unnormalized_moments_err.find(name);
    return (it != _unnormalized_moments_err.end()) ? it->second : 0.0;
  }

  /**
   * @brief Get normalized moment error
   * @param name Moment name
   * @return Moment error, or 0.0 if not found
   */
  Double_t GetError(const TString& name) const {
    auto it = _moments_err.find(name);
    return (it != _moments_err.end()) ? it->second : 0.0;
  }

  /**
   * @brief Set an unnormalized moment value and error directly
   * @param mom Moment name
   * @param val Moment value
   * @param error Moment error (default: 0.0)
   */
  void SetUnnormalized(const TString& mom, Double_t val, Double_t error = 0.0) {
    _unnormalized_moments[mom] = val;
    _unnormalized_moments_err[mom] = error;
  }

  /**
   * @brief Check whether a covariance matrix is available
   */
  bool HasCovarianceMatrix() const { return _has_cov; }

  /**
   * @brief Get reference to the unnormalized covariance matrix
   */
  const TMatrixDSym& GetCovarianceMatrix() const { return _unnormalized_cov; }

  /**
   * @brief Get list of moment names corresponding to covariance matrix rows/columns
   */
  const std::vector<TString>& GetCovarianceMomentNames() const { return _cov_moment_names; }

  /**
   * @brief Get index of a moment in the covariance matrix
   * @param name Moment name
   * @return Index if found, or -1
   */
  int GetCovarianceIndex(const TString& name) const {
    auto it = _cov_name_to_index.find(name);
    return (it != _cov_name_to_index.end()) ? it->second : -1;
  }

  /**
   * @brief Get covariance between two moments
   * @param name1 First moment name
   * @param name2 Second moment name
   * @return Covariance value, or variance if name1 == name2 and no cov matrix, or 0.0
   */
  Double_t GetCovariance(const TString& name1, const TString& name2) const {
    if (!_has_cov) {
      if (name1 == name2) {
        Double_t err = GetUnnormalizedError(name1);
        return err * err;
      }
      return 0.0;
    }
    auto it1 = _cov_name_to_index.find(name1);
    auto it2 = _cov_name_to_index.find(name2);
    if (it1 != _cov_name_to_index.end() && it2 != _cov_name_to_index.end()) {
      return _unnormalized_cov(it1->second, it2->second);
    }
    if (name1 == name2) {
      Double_t err = GetUnnormalizedError(name1);
      return err * err;
    }
    return 0.0;
  }

  /**
   * @brief Build a diagonal covariance matrix from existing unnormalized errors
   */
  void SetDiagonalCovariance() {
    _cov_moment_names.clear();
    _cov_name_to_index.clear();
    for (const auto& mom : _unnormalized_moments_err) {
      _cov_name_to_index[mom.first] = _cov_moment_names.size();
      _cov_moment_names.push_back(mom.first);
    }
    const int K = _cov_moment_names.size();
    if (K > 0) {
      _unnormalized_cov.ResizeTo(K, K);
      _unnormalized_cov.Zero();
      for (int i = 0; i < K; ++i) {
        Double_t err = _unnormalized_moments_err[_cov_moment_names[i]];
        _unnormalized_cov(i, i) = (err > 0.0) ? err * err : 1.0;
      }
      _has_cov = true;
    }
  }

  /**
   * @brief Set covariance matrix directly
   */
  void SetCovarianceMatrix(const std::vector<TString>& names, const TMatrixDSym& cov) {
    _cov_moment_names = names;
    _cov_name_to_index.clear();
    for (size_t i = 0; i < names.size(); ++i) {
      _cov_name_to_index[names[i]] = i;
    }
    _unnormalized_cov.ResizeTo(cov.GetNrows(), cov.GetNcols());
    _unnormalized_cov = cov;
    _has_cov = true;
  }

  /**
   * @brief Print the covariance matrix and associated moment names
   */
  void PrintCovarianceMatrix() const {
    if (!_has_cov) {
      std::cout << "No covariance matrix available." << std::endl;
      return;
    }
    const int K = _cov_moment_names.size();
    std::cout << "Covariance matrix (" << K << "x" << K << "):" << std::endl;
    std::cout << "Moments: [";
    for (int i = 0; i < K; ++i) {
      std::cout << _cov_moment_names[i] << (i + 1 < K ? ", " : "");
    }
    std::cout << "]" << std::endl;
    _unnormalized_cov.Print();
  }

  // =========================================================================
  // OPTION 2 (Commented out for now - working on MCMCTree first):
  // Load covariance matrix from Results ROOT file (e.g. PosteriorCovariance)
  // and convert to unnormalized moments using error propagation.
  // =========================================================================
  /*
  bool LoadCovarianceMatrix(const TString& filePath, const TString& matrixName = "PosteriorCovariance") {
    TFile* file = TFile::Open(filePath);
    if (!file || file->IsZombie()) {
      if (file) delete file;
      return false;
    }

    auto cov = file->Get<TMatrixDSym>(matrixName);
    if (!cov) {
      delete file;
      return false;
    }

    auto dataset = file->Get<RooDataSet>("FinalParameters");
    if (!dataset) {
      delete file;
      return false;
    }

    const RooArgSet* vars = dataset->get();
    if (!vars) {
      delete file;
      return false;
    }

    TString dirPath = gSystem->DirName(filePath);
    _filepath = dirPath;
    Double_t acceptance = 1.0;
    TFile* csFile = TFile::Open(dirPath + "/ResultsCrossSection.root");
    if (csFile && !csFile->IsZombie()) {
      HS::FIT::CrossSection* cs = (HS::FIT::CrossSection*)csFile->Get("cs");
      if (cs) {
        acceptance = cs->GetAcceptance();
      }
      delete csFile;
    } else {
      if (csFile) delete csFile;
      acceptance = GetAcceptance();
    }
    if (acceptance <= 0.0) acceptance = 1.0;

    int nRows = cov->GetNrows();
    std::vector<TString> varNames(nRows);
    int yldIdx = -1;
    Double_t yldVal = 1.0;
    Double_t yldVar = 0.0;

    int idx = 0;
    for (auto* obj : *vars) {
      auto var = dynamic_cast<RooRealVar*>(obj);
      if (!var || idx >= nRows) continue;
      TString vname = var->GetName();
      varNames[idx] = vname;
      if (vname == "Yld_Moments") {
        yldIdx = idx;
        yldVal = var->getVal();
        yldVar = (*cov)(idx, idx);
      }
      idx++;
    }

    std::vector<TString> momNames;
    std::vector<int> momOrigIdx;

    for (int i = 0; i < nRows; ++i) {
      if (varNames[i].BeginsWith("H_") && varNames[i] != "H_0_0_0") {
        momNames.push_back(varNames[i]);
        momOrigIdx.push_back(i);
      }
    }
    int h000IdxInMatrix = -1;
    if (yldIdx >= 0) {
      momNames.push_back("H_0_0_0");
      h000IdxInMatrix = static_cast<int>(momNames.size()) - 1;
    }

    const int K = momNames.size();
    if (K == 0) {
      delete file;
      return false;
    }

    _unnormalized_cov.ResizeTo(K, K);
    _cov_moment_names = momNames;
    _cov_name_to_index.clear();
    for (int i = 0; i < K; ++i) {
      _cov_name_to_index[momNames[i]] = i;
    }

    const Double_t invAcc2 = 1.0 / (acceptance * acceptance);

    for (int i = 0; i < K; ++i) {
      bool iIsH000 = (i == h000IdxInMatrix);
      TString name_i = momNames[i];
      Double_t mean_i = iIsH000 ? 2.0 : GetVal(name_i);
      if (!iIsH000 && mean_i == 0.0 && yldVal > 0.0) {
        mean_i = acceptance * GetUnnormalizedVal(name_i) / yldVal;
      }
      int orig_i = iIsH000 ? yldIdx : momOrigIdx[i];

      for (int k = i; k < K; ++k) {
        bool kIsH000 = (k == h000IdxInMatrix);
        TString name_k = momNames[k];
        Double_t mean_k = kIsH000 ? 2.0 : GetVal(name_k);
        if (!kIsH000 && mean_k == 0.0 && yldVal > 0.0) {
          mean_k = acceptance * GetUnnormalizedVal(name_k) / yldVal;
        }
        int orig_k = kIsH000 ? yldIdx : momOrigIdx[k];

        Double_t cov_ik = 0.0;
        if (iIsH000 && kIsH000) {
          cov_ik = 4.0 * invAcc2 * yldVar;
        } else if (iIsH000) {
          Double_t cov_Y_Hk = (yldIdx >= 0) ? (*cov)(yldIdx, orig_k) : 0.0;
          cov_ik = 2.0 * invAcc2 * (yldVal * cov_Y_Hk + mean_k * yldVar);
        } else if (kIsH000) {
          Double_t cov_Y_Hi = (yldIdx >= 0) ? (*cov)(yldIdx, orig_i) : 0.0;
          cov_ik = 2.0 * invAcc2 * (yldVal * cov_Y_Hi + mean_i * yldVar);
        } else {
          Double_t cov_Hi_Hk = (*cov)(orig_i, orig_k);
          Double_t cov_Y_Hi = (yldIdx >= 0) ? (*cov)(yldIdx, orig_i) : 0.0;
          Double_t cov_Y_Hk = (yldIdx >= 0) ? (*cov)(yldIdx, orig_k) : 0.0;
          cov_ik = invAcc2 * (yldVal * yldVal * cov_Hi_Hk +
                              yldVal * mean_i * cov_Y_Hk +
                              yldVal * mean_k * cov_Y_Hi +
                              mean_i * mean_k * yldVar);
        }

        _unnormalized_cov(i, k) = cov_ik;
        _unnormalized_cov(k, i) = cov_ik;
      }
    }

    _has_cov = true;
    delete file;
    return true;
  }
  */

  /**
   * @brief Print moment values with optional prefix filter
   * @param prefix Prefix to filter moments (default: empty = print all)
   */
  void PrintVals(const TString& prefix = "") const {
    for (const auto& mom : _moments) {
      if (mom.first.BeginsWith(prefix)) {
        Double_t error = GetError(mom.first);
        std::cout << mom.first << " = " << mom.second << " +- " << error << std::endl;
      }
    }
  }

private:
  
  /**
   * @brief Load MCMC results from file
   * @param filename Path to the result file
   */
  void LoadMCMCResults(const TString& filename) {
    auto file = std::unique_ptr<TFile>(TFile::Open(filename));
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
    }
    
    auto MCMC_tree = file->Get<TTree>("MCMCTree");
    if (!MCMC_tree) {
      std::cerr << "Error: MCMCTree not found in file " << filename << std::endl;
      return;
    }
    
    ProcessMCMCTree(MCMC_tree);
    
    // Set the normalization moment
    _moments["H_0_0_0"] = 2.0;
  }

  Double_t GetAcceptance() const {
    TFile *file = TFile::Open(_filepath + "/ResultsCrossSection.root");
    if (!file || !file->IsOpen()) {
      std::cerr << "Error: Cannot open file " << _filepath + "/ResultsCrossSection.root" << std::endl;
      return 1.0; // Default acceptance if file cannot be opened
    }
    HS::FIT::CrossSection* cross_section = (HS::FIT::CrossSection*)file->Get("cs");
    if (!cross_section) {
      std::cerr << "Error: CrossSection not found in file " << _filepath + "/ResultsCrossSection.root" << std::endl;
      return 1.0; // Default acceptance if cross section is not found
    }
    return cross_section->GetAcceptance();
  }

  /**
   * @brief Load direct fit results from file
   * @param filename Path to the result file
   */
  void LoadDirectResults(const TString& filename) {
    auto file = std::unique_ptr<TFile>(TFile::Open(filename));
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
    }
    
    auto dataset = file->Get<RooDataSet>("FinalParameters");
    if (!dataset) {
      std::cerr << "Error: FinalParameters dataset not found in file " << filename << std::endl;
      return;
    }
    
    Set(dataset->get());
    Set("H_0_0_0", 2.00000); // Normalization following BruFit
  }

  /**
   * @brief Process MCMC tree to extract moments
   * @param MCMC_tree Pointer to the MCMC tree
   */
  void ProcessMCMCTree(TTree* MCMC_tree) {
    Double_t Yld_Moments = 0.0;
    MCMC_tree->SetBranchAddress("Yld_Moments", &Yld_Moments);

    const Int_t nBranches = MCMC_tree->GetListOfBranches()->GetEntries();
    std::vector<TString> hBranches;

    for (Int_t i = 0; i < nBranches; ++i) {
      TBranch* var = static_cast<TBranch*>(MCMC_tree->GetListOfBranches()->At(i));
      TString var_name = var->GetName();
      
      if (var_name.BeginsWith("H_") && var_name != "H_0_0_0") {
        ProcessMomentBranch(MCMC_tree, var_name, Yld_Moments);
        hBranches.push_back(var_name);
      } else if (var_name.BeginsWith("Yld_")) {
        ProcessYieldBranch(MCMC_tree, var_name);
      }
    }

    ComputeMCMCCovariance(MCMC_tree, hBranches);
  }

  /**
   * @brief Compute sample covariance matrix of unnormalized moments from MCMCTree
   * @param tree Pointer to MCMCTree
   * @param hBranches List of H_ branch names
   */
  void ComputeMCMCCovariance(TTree* tree, const std::vector<TString>& hBranches) {
    const Long64_t nEntries = tree->GetEntries();
    if (nEntries <= 1) return;

    tree->ResetBranchAddresses();

    Double_t yld = 0.0;
    if (tree->GetBranch("Yld_Moments") == nullptr) return;
    tree->SetBranchAddress("Yld_Moments", &yld);

    Double_t acceptance = GetAcceptance();
    if (acceptance <= 0.0) acceptance = 1.0;

    std::vector<TString> momNames = hBranches;
    momNames.push_back("H_0_0_0");
    const int K = momNames.size();

    std::vector<Double_t> branchVals(hBranches.size(), 0.0);
    for (size_t i = 0; i < hBranches.size(); ++i) {
      tree->SetBranchAddress(hBranches[i], &branchVals[i]);
    }

    std::vector<std::vector<Double_t>> samples(K, std::vector<Double_t>(nEntries, 0.0));
    std::vector<Double_t> means(K, 0.0);

    for (Long64_t j = 0; j < nEntries; ++j) {
      tree->GetEntry(j);
      for (size_t i = 0; i < hBranches.size(); ++i) {
        Double_t unnorm = (yld * branchVals[i]) / acceptance;
        samples[i][j] = unnorm;
        means[i] += unnorm;
      }
      Double_t unnorm0 = (2.0 * yld) / acceptance;
      samples[K - 1][j] = unnorm0;
      means[K - 1] += unnorm0;
    }

    for (int i = 0; i < K; ++i) {
      means[i] /= nEntries;
    }

    _cov_moment_names = momNames;
    _cov_name_to_index.clear();
    for (int i = 0; i < K; ++i) {
      _cov_name_to_index[momNames[i]] = i;
    }

    _unnormalized_cov.ResizeTo(K, K);
    for (int i = 0; i < K; ++i) {
      for (int k = i; k < K; ++k) {
        Double_t sum_cov = 0.0;
        for (Long64_t j = 0; j < nEntries; ++j) {
          sum_cov += (samples[i][j] - means[i]) * (samples[k][j] - means[k]);
        }
        Double_t c = sum_cov / (nEntries - 1);
        _unnormalized_cov(i, k) = c;
        _unnormalized_cov(k, i) = c;
      }
    }
    _has_cov = true;
    tree->ResetBranchAddresses();
  }

  /**
   * @brief Process a moment branch from MCMC tree
   * @param tree Pointer to the MCMC tree
   * @param var_name Name of the variable
   * @param Yld_Moments Yield moments value
   */
  void ProcessMomentBranch(TTree* tree, const TString& var_name, Double_t& Yld_Moments) {
    Double_t val;
    tree->SetBranchAddress(var_name, &val);
    
    Double_t mean = 0.0;
    Double_t sigma = 0.0;
    Double_t unnormalized_mean = 0.0;
    Double_t unnormalized_sigma = 0.0;
    
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t j = 0; j < nEntries; ++j) {
      tree->GetEntry(j);
      mean += val;
      sigma += val * val;
      unnormalized_mean += Yld_Moments * val;
      unnormalized_sigma += (Yld_Moments * val) * (Yld_Moments * val);
    }
    
    // Calculate normalized moments
    mean /= nEntries;
    sigma = TMath::Sqrt(sigma / nEntries - mean * mean);
    _moments[var_name] = mean;
    _moments_err[var_name] = sigma;
    
    // Calculate unnormalized moments
    unnormalized_mean /= nEntries;
    Double_t acceptance = GetAcceptance();
    // cout << "Acceptance: " << acceptance << endl;
    if (acceptance <= 0.0) {
      std::cerr << "Warning: Acceptance is zero or negative, setting to 1.0" << std::endl;
      acceptance = 1.0; // Default to 1.0 if acceptance is invalid
    }
    // unnormalized_mean /= acceptance;
    // cout << "Unnormalized mean for " << var_name << ": " << unnormalized_mean << endl;
    unnormalized_sigma = TMath::Sqrt(unnormalized_sigma / nEntries - unnormalized_mean * unnormalized_mean);
    // cout << "Acceptance: " << acceptance << ", Unnormalized sigma for " << var_name << ": " << unnormalized_sigma << endl;
    _unnormalized_moments[var_name] = unnormalized_mean/acceptance;
    _unnormalized_moments_err[var_name] = unnormalized_sigma/acceptance;
    
    // TODO: include efficiency correction to the unnormalized moments
  }

  /**
   * @brief Process a yield branch from MCMC tree
   * @param tree Pointer to the MCMC tree
   * @param var_name Name of the variable
   */
  void ProcessYieldBranch(TTree* tree, const TString& var_name) {
    Double_t val;
    tree->SetBranchAddress(var_name, &val);
    
    Double_t mean = 0.0;
    Double_t sigma = 0.0;
    
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t j = 0; j < nEntries; ++j) {
      tree->GetEntry(j);
      mean += 2 * val;
      sigma += 4 * val * val;
    }
    
    mean /= nEntries;
    sigma = TMath::Sqrt(sigma / nEntries - mean * mean);
    Double_t acceptance = GetAcceptance();
    if (acceptance <= 0.0) {
      std::cerr << "Warning: Acceptance is zero or negative, setting to 1.0" << std::endl;
      acceptance = 1.0; // Default to 1.0 if acceptance is invalid
    }
    _unnormalized_moments["H_0_0_0"] = mean/acceptance;
    _unnormalized_moments_err["H_0_0_0"] = sigma/acceptance;
  }
  
  // Member variables
  TString _filepath;
  std::map<TString, Double_t> _moments;                    ///< Normalized moments
  std::map<TString, Double_t> _moments_err;                ///< Normalized moment errors
  std::map<TString, Double_t> _unnormalized_moments;       ///< Unnormalized moments
  std::map<TString, Double_t> _unnormalized_moments_err;   ///< Unnormalized moment errors
  TMatrixDSym _unnormalized_cov;                           ///< Unnormalized moments covariance matrix
  std::vector<TString> _cov_moment_names;                  ///< Names of moments in covariance matrix
  std::map<TString, int> _cov_name_to_index;               ///< Index lookup in covariance matrix
  bool _has_cov = false;                                   ///< Flag indicating if covariance matrix is available
};

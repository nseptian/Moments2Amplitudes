#pragma once

#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooArgList.h>
#include <TString.h>
#include <TFile.h>
#include <map>
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
  void Set(const TString& filename, const Int_t BruFitResult = 1, const Int_t MCMCResult = 1) {
    if (!BruFitResult) return;
    
    if (MCMCResult) {
      LoadMCMCResults(filename);
    } else {
      LoadDirectResults(filename);
    }
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
    for (Int_t i = 0; i < nBranches; ++i) {
      TBranch* var = static_cast<TBranch*>(MCMC_tree->GetListOfBranches()->At(i));
      TString var_name = var->GetName();
      
      if (var_name.BeginsWith("H_")) {
        ProcessMomentBranch(MCMC_tree, var_name, Yld_Moments);
      } else if (var_name.BeginsWith("Yld_")) {
        ProcessYieldBranch(MCMC_tree, var_name);
      }
    }
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
    unnormalized_sigma = TMath::Sqrt(unnormalized_sigma / nEntries - unnormalized_mean * unnormalized_mean);
    _unnormalized_moments[var_name] = unnormalized_mean;
    _unnormalized_moments_err[var_name] = unnormalized_sigma;
    
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
    _unnormalized_moments["H_0_0_0"] = mean;
    _unnormalized_moments_err["H_0_0_0"] = sigma;
  }
  
  // Member variables
  std::map<TString, Double_t> _moments;                    ///< Normalized moments
  std::map<TString, Double_t> _moments_err;                ///< Normalized moment errors
  std::map<TString, Double_t> _unnormalized_moments;       ///< Unnormalized moments
  std::map<TString, Double_t> _unnormalized_moments_err;   ///< Unnormalized moment errors
};

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

// Class to load moments from fit results
class MomentHelper{

 public:
  
  MomentHelper()=default;


  void Set(TString mom,Double_t val){
    _moments[mom]= val;
  }

  void Set(const TString& filename, const Int_t BruFitResult=1, const Int_t MCMCResult=1){
    if (BruFitResult) {
      if (MCMCResult) {
        auto file = std::unique_ptr<TFile>(TFile::Open(filename));
        auto MCMC_tree = file->Get<TTree>("MCMCTree");
        
        Double_t Yld_Moments = 0.0;
        MCMC_tree->SetBranchAddress("Yld_Moments", &Yld_Moments);

        for (size_t i = 0; i < MCMC_tree->GetListOfBranches()->GetEntries(); ++i) {
          TBranch* var = (TBranch*)MCMC_tree->GetListOfBranches()->At(i);
          TString var_name = var->GetName();
          if (var_name.BeginsWith("H_")) {
            Double_t val;
            MCMC_tree->SetBranchAddress(var_name, &val);
            Double_t mean=0.0;
            Double_t sigma=0.0;
            Double_t unnormalized_mean=0.0;
            Double_t unnormalized_sigma=0.0;

            for (size_t i = 0; i < MCMC_tree->GetEntries(); ++i) {
              MCMC_tree->GetEntry(i);
              mean += val;
              sigma += val*val;
              unnormalized_mean += Yld_Moments*val /2.0;
              unnormalized_sigma += (Yld_Moments*val/2.0)*(Yld_Moments*val/2.0);
            }
            sigma = TMath::Sqrt(sigma/MCMC_tree->GetEntries());
            _moments_err[var_name] = sigma;
            mean /= MCMC_tree->GetEntries();
            _moments[var_name] = mean;
            unnormalized_sigma = TMath::Sqrt(unnormalized_sigma/MCMC_tree->GetEntries());
            _unnormalized_moments_err[var_name] = unnormalized_sigma;
            unnormalized_mean /= MCMC_tree->GetEntries();
            _unnormalized_moments[var_name] = unnormalized_mean;
            // to do: include efficiency correction to the unnormalized moments
          }
          else if (var_name.BeginsWith("Yld_")) {
            Double_t val;
            MCMC_tree->SetBranchAddress(var_name, &val);
            Double_t mean=0.0;
            Double_t sigma=0.0;

            for (size_t i = 0; i < MCMC_tree->GetEntries(); ++i) {
              MCMC_tree->GetEntry(i);
              mean += val;
              sigma += val*val;
            }
            sigma = TMath::Sqrt(sigma/MCMC_tree->GetEntries());
            _unnormalized_moments_err["H_0_0_0"] = sigma;
            mean /= MCMC_tree->GetEntries();
            _unnormalized_moments["H_0_0_0"] = mean;
          }
        }
        // Set the normalization moment
        _moments["H_0_0_0"] = 2.0;
      }
      else {
        auto file = std::unique_ptr<TFile>(TFile::Open(filename));
        auto fileMoments = file->Get<RooDataSet>("FinalParameters")->get();
        Set(fileMoments);
        Set("H_0_0_0", 2.00000); //Normalization following BruFit
      }
    }     
  }

  void Set(const RooArgSet* pars){
    for(auto* mom:*pars){
      RooRealVar* rmom= dynamic_cast<RooRealVar*>(mom);
      if(rmom==nullptr) continue;
      if(mom->GetName()[0]!='H') continue;
      _moments[TString(mom->GetName())]=rmom->getVal();
    }
  }

  Double_t GetVal(const TString& name){
    if(_moments.find(name)==_moments.end()) return 0.0;
    return _moments[name];
  }

  void PrintVals(const TString& prefix=""){
    for(const auto& mom:_moments){
      if(mom.first.BeginsWith(prefix)){
        std::cout<<mom.first<<" = "<<mom.second<<" +- "<<_moments_err[mom.first]<<std::endl;
      }
    }
  }
  
 private:
  
  std::map<TString,Double_t> _moments;
  std::map<TString,Double_t> _moments_err;

  std::map<TString,Double_t> _unnormalized_moments;
  std::map<TString,Double_t> _unnormalized_moments_err;
  
};

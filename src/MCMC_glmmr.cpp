#include <RcppArmadillo.h>
#include "MCMC_glmmr.h"

//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat glmmr_sgd_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                            Rcpp::List SgdObj_List, Rcpp::List Para_List,
                            bool Interactive) {

  //Convet Rcpp::Lists to C++ structs
  datobj DatObj = ConvertDatObj(DatObj_List);
  hypara HyPara = ConvertHyPara(HyPara_List);
  sgdobj SgdObj = ConvertSgdObj(SgdObj_List);
  para Para = ConvertPara(Para_List);

  //Set objects to be used in the for loop
  int n_epochs = SgdObj.n_epochs;
  int n = DatObj.n;
  int S = SgdObj.S;
  arma::Col<int> Seqn = DatObj.Seqn;
  arma::colvec Probn = DatObj.Probn;
  int n_omega = DatObj.n_omega;
  std::pair<para, sgdobj> Update;
  arma::vec WhichBurnInProgress = SgdObj.WhichBurnInProgress;
  arma::vec WhichBurnInProgressInt = SgdObj.WhichBurnInProgressInt;
  
  //User output
  BeginBurnInProgress(SgdObj, Interactive);
  
  //Loop over epochs
  arma::mat Omega_out(n_omega, n_epochs);
  for (arma::uword e = 1; e < (n_epochs + 1); e++) {

    // Rcpp::Rcout << std::fixed << Para.Beta << Para.Sigma << std::endl;
    
    //Check for user interrupt every 1 epochs
    if (e % 1 == 0) Rcpp::checkUserInterrupt();
    
    //Sample mini-batches
    arma::colvec samps = sampleRcpp(Seqn, S, true, Probn);

    //Calculate gradient contribution for each subject
    arma::colvec grad_likelihood(n_omega, arma::fill::zeros);
    for (arma::uword s = 0; s < S; s++) {
      
      //Sample random effects
      int id = samps(s);
      arma::mat Gamma_i = SampleGamma(id, DatObj, HyPara, SgdObj, Para);
      
      //Compute likelihood contribution for group i
      grad_likelihood += ComputeGradienti(id, Gamma_i, DatObj, HyPara, SgdObj, Para);
      
    //End loop over mini-batch samples 
    }
    
    //Calculate gradient contribution from the prior          
    arma::colvec grad_prior = ComputeGradientPrior(DatObj, HyPara, Para);
  
    //Final gradient computation
    arma::colvec grad = grad_prior + (n / S) * grad_likelihood;

    //Take a step in parameter space and update all parameters 
    Update = UpdateOmega(grad, DatObj, HyPara, SgdObj, Para);
    Para = Update.first;
    SgdObj = Update.second;

    //Update burn-in progress bar
    if (Interactive) if (std::find(WhichBurnInProgress.begin(), WhichBurnInProgress.end(), e) != WhichBurnInProgress.end())
      UpdateBurnInBar(e, SgdObj);
    if (!Interactive) if (std::find(WhichBurnInProgressInt.begin(), WhichBurnInProgressInt.end(), e) != WhichBurnInProgressInt.end())
      UpdateBurnInBarInt(e, SgdObj);
    
    //Save output
    Omega_out.col(e - 1) = Para.Omega;
    
  //End loop over epochs
  }

  //Return raw samples
  return Omega_out;

//End MCMC sampler function
}

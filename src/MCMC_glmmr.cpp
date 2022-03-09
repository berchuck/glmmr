#include <RcppArmadillo.h>
#include "MCMC_glmmr.h"

//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List glmmr_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                      Rcpp::List TuningObj_List, Rcpp::List Para_List,
                      bool Interactive) {

  //Convert Rcpp::Lists to C++ structs
  datobj DatObj = ConvertDatObj(DatObj_List);
  hypara HyPara = ConvertHyPara(HyPara_List);
  tuning TuningObj = ConvertTuningObj(TuningObj_List);
  para Para = ConvertPara(Para_List);

  //Set objects to be used in the for loop
  int NEpochs = TuningObj.NEpochs;
  int NTotal = TuningObj.NTotal;
  int NUnits = DatObj.NUnits;
  int S = TuningObj.S;
  int NOmega = DatObj.NOmega;
  int NKeep = TuningObj.NKeep;
  arma::Col<int> SeqNUnits = DatObj.SeqNUnits;
  arma::colvec ProbNUnits = DatObj.ProbNUnits;
  arma::vec WhichKeep = TuningObj.WhichKeep;
  arma::vec WhichMAPProgress = TuningObj.WhichMAPProgress;
  arma::vec WhichMAPProgressInt = TuningObj.WhichMAPProgressInt;
  arma::vec WhichSamplerProgress = TuningObj.WhichSamplerProgress;
  std::pair<para, tuning> Update;
  
  //User output
  BeginMAPProgress(TuningObj, Interactive);
  
  //Initialize objects
  int Id;
  arma::mat OmegaMat(NOmega, NKeep), GammaI(DatObj.Q, TuningObj.R);
  arma::colvec Samps(S), GradLikelihood(NOmega), GradPrior(NOmega), Grad(NOmega);
  arma::colvec OmegaMAP(NOmega);
  
  //Loop over epochs
  for (arma::uword e = 1; e < (NTotal + 1); e++) {

    //Check for user interrupt every 1 epochs
    if (e % 1 == 0) Rcpp::checkUserInterrupt();
    
    //Sample mini-batches
    Samps = sampleRcpp(SeqNUnits, S, true, ProbNUnits);

    //Calculate gradient contribution for each subject
    GradLikelihood.zeros();
    for (arma::uword s = 0; s < S; s++) {

      //Sample random effects
      Id = Samps(s);
      GammaI = SampleGamma(Id, DatObj, TuningObj, Para);

      //Compute likelihood contribution for group i
      GradLikelihood += ComputeGradientI(Id, GammaI, DatObj, TuningObj, Para);

    //End loop over mini-batch samples
    }
    
    //Calculate gradient contribution from the prior          
    GradPrior = ComputeGradientPrior(DatObj, HyPara, Para);
  
    //Final gradient computation
    Grad = GradPrior + (NUnits / S) * GradLikelihood;

    //Update omega
    Update = UpdateOmega(e, Grad, DatObj, TuningObj, Para);
    Para = Update.first;
    TuningObj = Update.second;
    
    //Update other parameters
    Para = UpdatePara(DatObj, Para);
    
    //Save parameters
    if (e == NEpochs) OmegaMAP = Para.Omega;
    if (std::find(WhichKeep.begin(), WhichKeep.end(), e) != WhichKeep.end())
      OmegaMat.cols(find(e == WhichKeep)) = Para.Omega;
      
    //Update MAP progress bar
    if (Interactive) if (std::find(WhichMAPProgress.begin(), WhichMAPProgress.end(), e) != WhichMAPProgress.end())
      UpdateMAPBar(e, TuningObj);
    if (!Interactive) if (std::find(WhichMAPProgressInt.begin(), WhichMAPProgressInt.end(), e) != WhichMAPProgressInt.end())
      UpdateMAPBarInt(e, TuningObj);
    
    //Update SGMCMC progress bar
    if (e == NEpochs) Rcpp::Rcout << std::fixed << "\nSampler progress:  0%.. ";
    if (std::find(WhichSamplerProgress.begin(), WhichSamplerProgress.end(), e) != WhichSamplerProgress.end())
      SamplerProgress(e, TuningObj);
    
  //End loop over epochs
  }

  //Return raw samples
  return Rcpp::List::create(Rcpp::Named("map") = OmegaMAP,
                            Rcpp::Named("samples") = OmegaMat);
  
//End MCMC sampler function
}

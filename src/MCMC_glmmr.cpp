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
  int NTune = TuningObj.NTune;
  int NUnits = DatObj.NUnits;
  int S = TuningObj.S;
  int NOmega = DatObj.NOmega;
  int NKeep = TuningObj.NKeep;
  int AlgorithmInd = DatObj.AlgorithmInd;
  int Timer = DatObj.Timer;
  arma::Col<int> SeqNUnits = DatObj.SeqNUnits;
  arma::colvec ProbNUnits = DatObj.ProbNUnits;
  arma::vec WhichKeep = TuningObj.WhichKeep;
  arma::vec WhichMAPProgress = TuningObj.WhichMAPProgress;
  arma::vec WhichMAPProgressInt = TuningObj.WhichMAPProgressInt;
  arma::vec WhichSamplerProgress = TuningObj.WhichSamplerProgress;
  arma::vec WhichSamplerProgressInt = TuningObj.WhichSamplerProgressInt;
  arma::vec WhichTuneProgress = TuningObj.WhichTuneProgress;
  arma::vec WhichTuneProgressInt = TuningObj.WhichTuneProgressInt;
  arma::vec WhichPilotAdapt = TuningObj.WhichPilotAdapt;
  std::pair<para, tuning> Update;
  
  //User output
  BeginMAPProgress(TuningObj, Interactive, AlgorithmInd);
  
  //Initialize objects
  int Id;
  arma::mat OmegaMat(NOmega, NKeep), GammaI(DatObj.Q, TuningObj.R);
  arma::colvec GradLikelihood(NOmega), GradPrior(NOmega), Grad(NOmega);
  arma::colvec OmegaMAP(NOmega), TimerOut(NKeep);
  arma::uvec Samps(S);
  arma::wall_clock timer;
  timer.tic();
  double secs;  

  //Loop over epochs
  for (arma::uword e = 0; e < NTotal; e++) {

    //Check for user interrupt every 1 epochs
    if (e % 1 == 0) Rcpp::checkUserInterrupt();
    
    //Stochastic gradient algorithms
    if (AlgorithmInd < 3) {
    
      //Sample mini-batches
      // arma::colvec Samps = sampleRcpp(SeqNUnits, S, true, ProbNUnits);
      Samps = arma::randperm(NUnits, S);
      
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
    
      //Update Tune progress bar
      if (AlgorithmInd > 0) { // sgld and sgld_corrected
        if (e == NEpochs) BeginTuneProgress(TuningObj, Interactive);
        if (Interactive) if (std::find(WhichTuneProgress.begin(), WhichTuneProgress.end(), e) != WhichTuneProgress.end())
          UpdateTuneBar(e, TuningObj);
        if (!Interactive) if (std::find(WhichTuneProgressInt.begin(), WhichTuneProgressInt.end(), e) != WhichTuneProgressInt.end())
          UpdateTuneBarInt(e, TuningObj);
      }
    
    //End if statement over gradient algorithms
    }
    
    //Gibbs Algorithm
    if (AlgorithmInd == 3) {
      
      //Update Omega
      Para = Sampleomega(DatObj, Para);
      
      //Update Beta
      Para = SampleBeta(DatObj, Para);
      
      //Update Gamma
      Para = SampleGammaGibbs(DatObj, Para);
      
      //Update L
      Update = UpdateL(DatObj, HyPara, TuningObj, Para);
      Para = Update.first;
      TuningObj = Update.second;

      //Update D
      Update = UpdateD(DatObj, HyPara, TuningObj, Para);
      Para = Update.first;
      TuningObj = Update.second;
      
    }
    
    //Save parameters
    if (AlgorithmInd == 0) { // sgd
      OmegaMat.col(e) = Para.Omega;
    }
    if (e == (NEpochs - 2)) OmegaMAP = Para.Omega;
    if (AlgorithmInd > 0) { // sgld, sgld_corrected, gibbs
      if (std::find(WhichKeep.begin(), WhichKeep.end(), e) != WhichKeep.end())
        OmegaMat.cols(find(e == WhichKeep)) = Para.Omega;
    } 
    
    //Pilot adaptation
    if (std::find(WhichPilotAdapt.begin(), WhichPilotAdapt.end(), e) != WhichPilotAdapt.end())
      TuningObj = PilotAdaptation(DatObj, TuningObj);
    
    //Update MAP progress bar
    if (Interactive) if (std::find(WhichMAPProgress.begin(), WhichMAPProgress.end(), e) != WhichMAPProgress.end())
      UpdateMAPBar(e, TuningObj, AlgorithmInd);
    if (!Interactive) if (std::find(WhichMAPProgressInt.begin(), WhichMAPProgressInt.end(), e) != WhichMAPProgressInt.end())
      UpdateMAPBarInt(e, TuningObj);
  
    //Update SGMCMC progress bar
    if (AlgorithmInd > 0) { // sgld and sgld_corrected
      if (e == (NEpochs + NTune)) BeginSamplerProgress(TuningObj, Interactive);
      if (Interactive) if (std::find(WhichSamplerProgress.begin(), WhichSamplerProgress.end(), e) != WhichSamplerProgress.end())
        UpdateSamplerBar(e, TuningObj);
      if (!Interactive) if (std::find(WhichSamplerProgressInt.begin(), WhichSamplerProgressInt.end(), e) != WhichSamplerProgressInt.end())
        UpdateSamplerBarInt(e, TuningObj);
    }
    
    //Compute SGLD Correction
    if (AlgorithmInd == 2 & e == (NEpochs - 1)) {
      Update = ComputeSGLDCorrection(DatObj, TuningObj, Para, Interactive);
      Para = Update.first;
      TuningObj = Update.second;
    }
    
    //Check to see if the sampler should break
    if (Timer > 0) {
      secs = timer.toc();
      TimerOut(e) = secs;
      if (secs > Timer) {
        break;
      }
    }
    
  //End loop over epochs
  }
  
  //Output Metropolis object for summary
  Rcpp::List Metropolis;
  if (AlgorithmInd == 3) Metropolis = OutputMetrObj(TuningObj);

  //Return raw samples
  return Rcpp::List::create(Rcpp::Named("map") = OmegaMAP, 
                           Rcpp::Named("samples") = OmegaMat,
                           Rcpp::Named("metropolis") = Metropolis,
                           Rcpp::Named("timer") = TimerOut,
                           Rcpp::Named("epsilon") = TuningObj.EpsilonSGLD,
                           Rcpp::Named("counter") = TuningObj.Counter,
                           Rcpp::Named("epsilontunecounter") = TuningObj.EpsilonTuneCounter);

//End MCMC sampler function
}

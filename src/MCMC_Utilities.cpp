#include <RcppArmadillo.h>
#include "MCMC_glmmr.h"

//Initiate burn-in progress bar--------------------------------------------------------------------------
void BeginMAPProgress(tuning TuningObj, bool Interactive) {

  //Set MCMC object
  int BarLength = TuningObj.BarLength;

  //Initialize burn-in bar
  if (Interactive) {
    Rcpp::Rcout << std::fixed << "Identifying MAP:  |";
    for (int i = 0; i < BarLength - 1; i++) Rcpp::Rcout << std::fixed << " ";
    Rcpp::Rcout << std::fixed <<  "|" << std::fixed;
  }
  if (!Interactive) {
    Rcpp::Rcout << std::fixed << "Identifying MAP:  0%..  ";
  }

}



//Initiate burn-in progress bar-------------------------------------------------------------------------------------
void SamplerProgress(int e, tuning TuningObj) {

  //Set MCMC object
  int NSims = TuningObj.NSims;
  int NEpochs = TuningObj.NEpochs;

  //Add a new percentage
  Rcpp::Rcout.precision(0);
  if (e < NSims + NEpochs) Rcpp::Rcout << std::fixed << 100 * (e - NEpochs) / NSims << "%.. ";
  if (e == NSims + NEpochs) Rcpp::Rcout << std::fixed << 100 * (e - NEpochs) / NSims << "%!";

}



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateMAPBarInt(int e, tuning TuningObj) {

  //Set MCMC object
  arma::vec WhichMAPProgressInt = TuningObj.WhichMAPProgressInt;
  arma::uvec NewStarBoolean = find(e == WhichMAPProgressInt);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);

  //Add percentage to submited job mode
  Rcpp::Rcout.precision(0);
  if (NewStar == 0) Rcpp::Rcout << std::fixed << "10%.. ";
  if (NewStar == 1) Rcpp::Rcout << std::fixed << "20%.. ";
  if (NewStar == 2) Rcpp::Rcout << std::fixed << "30%.. ";
  if (NewStar == 3) Rcpp::Rcout << std::fixed << "40%.. ";
  if (NewStar == 4) Rcpp::Rcout << std::fixed << "50%.. ";
  if (NewStar == 5) Rcpp::Rcout << std::fixed << "60%.. ";
  if (NewStar == 6) Rcpp::Rcout << std::fixed << "70%.. ";
  if (NewStar == 7) Rcpp::Rcout << std::fixed << "80%.. ";
  if (NewStar == 8) Rcpp::Rcout << std::fixed << "90%.. ";
  if (NewStar == 9) Rcpp::Rcout << std::fixed << "100%!";

}



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateMAPBar(int e, tuning TuningObj) {

  //Set MCMC object
  arma::vec WhichMAPProgress = TuningObj.WhichMAPProgress;
  int BarLength = TuningObj.BarLength;

  //Add a new star
  arma::uvec NewStarBoolean = find(e == WhichMAPProgress);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);
  Rcpp::Rcout << std::fixed << "\rIdentifying MAP:  |";
  for (int i = 0; i < NewStar; i++) Rcpp::Rcout << std::fixed << "*";
  for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
  Rcpp::Rcout << std::fixed << "|";

}
#include <RcppArmadillo.h>
#include "MCMC_glmmr.h"

//Initiate burn-in progress bar--------------------------------------------------------------------------
void BeginMAPProgress(tuning TuningObj, bool Interactive, int AlgorithmInd) {

  //Set MCMC object
  int BarLength = TuningObj.BarLength;

  //Initialize burn-in bar
  if (Interactive) {
    if (AlgorithmInd < 3) Rcpp::Rcout << std::fixed << "Identifying MAP:  |";
    if (AlgorithmInd == 3) Rcpp::Rcout << std::fixed << "Burnin Progress:  |";
    for (int i = 0; i < BarLength - 1; i++) Rcpp::Rcout << std::fixed << " ";
    Rcpp::Rcout << std::fixed <<  "|" << std::fixed;
  }
  if (!Interactive) {
    if (AlgorithmInd < 3) Rcpp::Rcout << std::fixed << "Identifying MAP:  0%..  ";
    if (AlgorithmInd == 3) Rcpp::Rcout << std::fixed << "Burnin Progress:  0%..  ";
  }

}



//Initiate burn-in progress bar--------------------------------------------------------------------------
void BeginSGLDProgress(tuning TuningObj, bool Interactive) {
  
  //Set MCMC object
  int BarLength = TuningObj.BarLength;
  
  //Initialize burn-in bar
  if (Interactive) {
    Rcpp::Rcout << std::fixed << "\nSGLD Correction:  |";
    for (int i = 0; i < BarLength - 1; i++) Rcpp::Rcout << std::fixed << " ";
    Rcpp::Rcout << std::fixed <<  "|" << std::fixed;
  }
  if (!Interactive) {
    Rcpp::Rcout << std::fixed << "\nSGLD Correction:  0%..  ";
  }
  
}



//Initiate burn-in progress bar--------------------------------------------------------------------------
void BeginTuneProgress(tuning TuningObj, bool Interactive) {
  
  //Set MCMC object
  int BarLength = TuningObj.BarLength;
  
  //Initialize burn-in bar
  if (Interactive) {
    Rcpp::Rcout << std::fixed << "\nTuning Epsilon :  |";
    for (int i = 0; i < BarLength - 1; i++) Rcpp::Rcout << std::fixed << " ";
    Rcpp::Rcout << std::fixed <<  "|" << std::fixed;
  }
  if (!Interactive) {
    Rcpp::Rcout << std::fixed << "\nTuning Epsilon :  0%..  ";
  }
  
}



//Initiate burn-in progress bar--------------------------------------------------------------------------
void BeginSamplerProgress(tuning TuningObj, bool Interactive) {
  
  //Set MCMC object
  int BarLength = TuningObj.BarLength;
  
  //Initialize burn-in bar
  if (Interactive) {
    Rcpp::Rcout << std::fixed << "\nSample Progress:  |";
    for (int i = 0; i < BarLength - 1; i++) Rcpp::Rcout << std::fixed << " ";
    Rcpp::Rcout << std::fixed <<  "|" << std::fixed;
  }
  if (!Interactive) {
    Rcpp::Rcout << std::fixed << "\nSample Progress:  0%..  ";
  }
  
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
void UpdateSGLDBarInt(int e, tuning TuningObj) {
  
  //Set MCMC object
  arma::vec WhichSGLDProgressInt = TuningObj.WhichSGLDProgressInt;
  arma::uvec NewStarBoolean = find(e == WhichSGLDProgressInt);
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
void UpdateTuneBarInt(int e, tuning TuningObj) {
  
  //Set MCMC object
  arma::vec WhichTuneProgressInt = TuningObj.WhichTuneProgressInt;
  arma::uvec NewStarBoolean = find(e == WhichTuneProgressInt);
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
void UpdateSamplerBarInt(int e, tuning TuningObj) {
  
  //Set MCMC object
  arma::vec WhichSamplerProgressInt = TuningObj.WhichSamplerProgressInt;
  arma::uvec NewStarBoolean = find(e == WhichSamplerProgressInt);
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
void UpdateMAPBar(int e, tuning TuningObj, int AlgorithmInd) {

  //Set MCMC object
  arma::vec WhichMAPProgress = TuningObj.WhichMAPProgress;
  int BarLength = TuningObj.BarLength;

  //Add a new star
  arma::uvec NewStarBoolean = find(e == WhichMAPProgress);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);
  if (AlgorithmInd < 3) Rcpp::Rcout << std::fixed << "\rIdentifying MAP:  |";
  if (AlgorithmInd == 3) Rcpp::Rcout << std::fixed << "\rBurnin Progress:  |";
  for (int i = 0; i < NewStar; i++) Rcpp::Rcout << std::fixed << "*";
  for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
  Rcpp::Rcout << std::fixed << "|";

}



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateSGLDBar(int e, tuning TuningObj) {
  
  //Set MCMC object
  arma::vec WhichSGLDProgress = TuningObj.WhichSGLDProgress;
  int BarLength = TuningObj.BarLength;
  
  //Add a new star
  arma::uvec NewStarBoolean = find(e == WhichSGLDProgress);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);
  Rcpp::Rcout << std::fixed << "\rSGLD Correction:  |";
  for (int i = 0; i < NewStar; i++) Rcpp::Rcout << std::fixed << "*";
  for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
  Rcpp::Rcout << std::fixed << "|";
  
}



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateTuneBar(int e, tuning TuningObj) {
  
  //Set MCMC object
  arma::vec WhichTuneProgress = TuningObj.WhichTuneProgress;
  int BarLength = TuningObj.BarLength;
  
  //Add a new star
  arma::uvec NewStarBoolean = find(e == WhichTuneProgress);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);
  Rcpp::Rcout << std::fixed << "\rTuning Epsilon :  |";
  for (int i = 0; i < NewStar; i++) Rcpp::Rcout << std::fixed << "*";
  for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
  Rcpp::Rcout << std::fixed << "|";
  
}



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateSamplerBar(int e, tuning TuningObj) {
  
  //Set MCMC object
  arma::vec WhichSamplerProgress = TuningObj.WhichSamplerProgress;
  int BarLength = TuningObj.BarLength;
  
  //Add a new star
  arma::uvec NewStarBoolean = find(e == WhichSamplerProgress);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);
  Rcpp::Rcout << std::fixed << "\rSample Progress:  |";
  for (int i = 0; i < NewStar; i++) Rcpp::Rcout << std::fixed << "*";
  for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
  Rcpp::Rcout << std::fixed << "|";
  
}



//Function to pilot adapt tuning parameter--------------------------------------------------------------
double PilotAdaptFunc(double TuningParameter, double AcceptancePct) {
  
  //Adjust tuning parameter using scaling based on size of acceptance rate
  if (AcceptancePct >= 0.90) TuningParameter *= 1.3;
  if ( (AcceptancePct >= 0.75 ) & (AcceptancePct < 0.90 ) ) TuningParameter *= 1.2;
  if ( (AcceptancePct >= 0.45 ) & (AcceptancePct < 0.75 ) ) TuningParameter *= 1.1;
  if ( (AcceptancePct <= 0.25 ) & (AcceptancePct > 0.15 ) ) TuningParameter *= 0.9;
  if ( (AcceptancePct <= 0.15 ) & (AcceptancePct > 0.10 ) ) TuningParameter *= 0.8;
  if (AcceptancePct <= 0.10) TuningParameter *= 0.7;
  return TuningParameter;
  
}



//Function for implementing pilot adaptation in MCMC sampler--------------------------------------------
tuning PilotAdaptation(datobj DatObj, tuning TuningObj) {
  
  //Set data objects
  int Q = DatObj.Q;
  int NL = DatObj.NL;
  
  //Set Metropolis objects
  arma::vec MetropL = TuningObj.MetropL;
  arma::vec AcceptanceL = TuningObj.AcceptanceL;
  arma::vec MetropD = TuningObj.MetropD;
  arma::vec AcceptanceD = TuningObj.AcceptanceD;

  //Set MCMC objects
  int PilotAdaptDenominator = TuningObj.PilotAdaptDenominator;
  
  //Get acceptance percentages
  arma::vec PctL = AcceptanceL / double(PilotAdaptDenominator);
  arma::vec PctD = AcceptanceD / double(PilotAdaptDenominator);

  //Update Tuning Parameter
  for (int i = 0; i < NL; i++) MetropL(i) = PilotAdaptFunc(MetropL(i), PctL(i));
  for (int i = 0; i < Q; i++) MetropD(i) = PilotAdaptFunc(MetropD(i), PctD(i));
  TuningObj.MetropL = MetropL;
  TuningObj.MetropD = MetropD;

  //Zero the acceptance counters
  AcceptanceL.zeros();
  AcceptanceD.zeros();
  TuningObj.AcceptanceL = AcceptanceL;
  TuningObj.AcceptanceD = AcceptanceD;
  return TuningObj;
  
}



//Output Metropolis object for summary-------------------------------------------------------------------
Rcpp::List OutputMetrObj(tuning TuningObj) {
  
  Rcpp::List Out = Rcpp::List::create(Rcpp::Named("AcceptanceL") = TuningObj.AcceptanceL,
                                      Rcpp::Named("MetropL") = TuningObj.MetropL,
                                      Rcpp::Named("AcceptanceD") = TuningObj.AcceptanceD,
                                      Rcpp::Named("MetropD") = TuningObj.MetropD);
  return Out;
  
}

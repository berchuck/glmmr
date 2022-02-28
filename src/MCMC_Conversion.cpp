#include <RcppArmadillo.h>
#include "MCMC_glmmr.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobj ConvertDatObj(Rcpp::List DatObj_List) {

  //Set objects from List
  arma::colvec Y = DatObj_List["Y"];
  arma::mat X = DatObj_List["X"];
  arma::mat Z = DatObj_List["Z"];
  arma::Col<int> group = DatObj_List["group"];
  arma::Col<int> group2 = DatObj_List["group2"];
  int N = DatObj_List["N"];
  int n = DatObj_List["n"];
  int p = DatObj_List["p"];
  int q = DatObj_List["q"];
  int n_L = DatObj_List["n_L"];
  int FamilyInd = DatObj_List["FamilyInd"];
  arma::mat EyeQ = DatObj_List["EyeQ"];
  arma::Col<int> Seqn = arma::linspace<arma::Col<int>>(1, n, n);
  arma::colvec Probn = arma::ones<arma::colvec>(n) / n;
  int n_omega = p + q + n_L;
  
  //Convert to C++ struct
  datobj DatObj;
  DatObj.Y = Y;
  DatObj.X = X;
  DatObj.Z = Z;
  DatObj.group = group;
  DatObj.group2 = group2;
  DatObj.N = N;
  DatObj.n = n;
  DatObj.p = p;
  DatObj.q = q;
  DatObj.n_L = n_L;
  DatObj.FamilyInd = FamilyInd;
  DatObj.EyeQ = EyeQ;
  DatObj.Seqn = Seqn;
  DatObj.Probn = Probn;
  DatObj.n_omega = n_omega;
  return DatObj;

}



//Function to convert Rcpp::List HyPara to a custom C++ struct hypara--------------------------------------------------
hypara ConvertHyPara(Rcpp::List HyPara_List) {
  
  //Set objects from List
  double Eta = HyPara_List["Eta"];
  double Nu = HyPara_List["Nu"];
  
  //Convert to C++ struct
  hypara HyPara;
  HyPara.Eta = Eta;
  HyPara.Nu = Nu;
  return HyPara;

}


  
//Function to convert Rcpp::List SgdObj to a custom C++ struct sgdobj-----------------------------------------------
sgdobj ConvertSgdObj(Rcpp::List SgdObj_List) {
  
  //Set objects from List
  double Epsilon = SgdObj_List["Epsilon"];
  arma::vec M_nadam = SgdObj_List["M_nadam"];
  arma::vec N_nadam = SgdObj_List["N_nadam"];
  double Mu_nadam = SgdObj_List["Mu_nadam"];
  double Alpha_nadam = SgdObj_List["Alpha_nadam"];
  double Nu_nadam = SgdObj_List["Nu_nadam"];
  double S = SgdObj_List["S"];
  double n_epochs = SgdObj_List["n_epochs"];
  double R = SgdObj_List["R"];
  int BarLength = SgdObj_List["BarLength"];
  arma::vec WhichBurnInProgressInt = SgdObj_List["WhichBurnInProgressInt"];
  arma::vec WhichBurnInProgress = SgdObj_List["WhichBurnInProgress"];
  
  //Convert to C++ struct
  sgdobj SgdObj;
  SgdObj.Epsilon = Epsilon;
  SgdObj.M_nadam = M_nadam;
  SgdObj.N_nadam = N_nadam;
  SgdObj.Mu_nadam = Mu_nadam;
  SgdObj.Alpha_nadam = Alpha_nadam;
  SgdObj.Nu_nadam = Nu_nadam;
  SgdObj.S = S;
  SgdObj.n_epochs = n_epochs;
  SgdObj.R = R;
  SgdObj.BarLength = BarLength;
  SgdObj.WhichBurnInProgressInt = WhichBurnInProgressInt;
  SgdObj.WhichBurnInProgress = WhichBurnInProgress;
  return SgdObj;
  
}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
para ConvertPara(Rcpp::List Para_List) {

  //Set objects from List
  arma::colvec Beta = Para_List["Beta"];
  arma::colvec l = Para_List["l"];
  arma::colvec d = Para_List["d"];
  arma::colvec Omega = Para_List["Omega"];
  arma::mat z = Para_List["z"];
  arma::mat L = Para_List["L"];
  arma::vec Loverz = Para_List["Loverz"];
  arma::mat D = Para_List["D"];
  arma::mat Upsilon = Para_List["Upsilon"];
  arma::mat Sigma = Para_List["Sigma"];
  arma::mat LInv = Para_List["LInv"];
  arma::mat tLInv = Para_List["tLInv"];
  arma::mat UpsilonInv = Para_List["UpsilonInv"];
  arma::mat DInv = Para_List["DInv"];
  arma::mat SigmaInv = Para_List["SigmaInv"];
  arma::mat gradLz = Para_List["gradLz"];
  arma::mat gradzl = Para_List["gradzl"];
  arma::mat gradLl = Para_List["gradLl"];
    
  //Convert to C++ struct
  para Para;
  Para.Beta = Beta;
  Para.l = l;
  Para.d = d;
  Para.Omega = Omega;
  Para.z = z;
  Para.L = L;
  Para.Loverz = Loverz;
  Para.D = D;
  Para.Upsilon = Upsilon;
  Para.Sigma = Sigma;
  Para.LInv = LInv;
  Para.tLInv = tLInv;
  Para.UpsilonInv = UpsilonInv;
  Para.DInv = DInv;
  Para.SigmaInv = SigmaInv;
  Para.gradLz = gradLz;
  Para.gradzl = gradzl;
  Para.gradLl = gradLl;
  return Para;
}
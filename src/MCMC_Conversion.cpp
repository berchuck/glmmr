#include <RcppArmadillo.h>
#include "MCMC_glmmr.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobj ConvertDatObj(Rcpp::List DatObj_List) {

  //Set objects from List
  arma::colvec Y = DatObj_List["Y"];
  arma::mat X = DatObj_List["X"];
  // arma::umat Zlocations = DatObj_List["Zlocations"];
  // arma::vec Zvalues = DatObj_List["Zvalues"];
  Rcpp::List ZList = DatObj_List["ZList"];
  arma::Col<int> Group = DatObj_List["Group"];
  arma::Col<int> Group2 = DatObj_List["Group2"];
  int N = DatObj_List["N"];
  int NOmega = DatObj_List["NOmega"];
  int NUnits = DatObj_List["NUnits"];
  int P = DatObj_List["P"];
  int Q = DatObj_List["Q"];
  int NL = DatObj_List["NL"];
  int FamilyInd = DatObj_List["FamilyInd"];
  int AlgorithmInd = DatObj_List["AlgorithmInd"];
  arma::mat EyeQ = DatObj_List["EyeQ"];
  arma::mat EyeNOmega = DatObj_List["EyeNOmega"];
  arma::Col<int> SeqNUnits = arma::linspace<arma::Col<int>>(1, NUnits, NUnits);
  arma::colvec ProbNUnits = arma::ones<arma::colvec>(NUnits) / NUnits;
  
  //Construct Z
  arma::field<arma::mat> Z(NUnits);
  for (arma::uword i = 0; i < NUnits; i++) {
   arma::mat z_i = ZList[i];
    Z(i) = z_i; 
  }
  // arma::sp_mat Z(Zlocations, Zvalues, N, Q * NUnits);
  
  //Convert to C++ struct
  datobj DatObj;
  DatObj.Y = Y;
  DatObj.X = X;
  DatObj.Z = Z;
  DatObj.Group = Group;
  DatObj.Group2 = Group2;
  DatObj.N = N;
  DatObj.NUnits = NUnits;
  DatObj.P = P;
  DatObj.Q = Q;
  DatObj.NL = NL;
  DatObj.FamilyInd = FamilyInd;
  DatObj.AlgorithmInd = AlgorithmInd;
  DatObj.EyeQ = EyeQ;
  DatObj.EyeNOmega = EyeNOmega;
  DatObj.SeqNUnits = SeqNUnits;
  DatObj.ProbNUnits = ProbNUnits;
  DatObj.NOmega = NOmega;
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


  
//Function to convert Rcpp::List TuningObj to a custom C++ struct tuning-----------------------------------------------
tuning ConvertTuningObj(Rcpp::List TuningObj_List) {
  
  //Set objects from List
  double EpsilonNADAM = TuningObj_List["EpsilonNADAM"];
  double MuNADAM = TuningObj_List["MuNADAM"];
  double AlphaNADAM = TuningObj_List["AlphaNADAM"];
  double NuNADAM = TuningObj_List["NuNADAM"];
  arma::vec MNADAM = TuningObj_List["MNADAM"];
  arma::vec NNADAM = TuningObj_List["NNADAM"];
  int S = TuningObj_List["S"];
  int S_SGLD = TuningObj_List["S_SGLD"];
  int NEpochs = TuningObj_List["NEpochs"];
  int R = TuningObj_List["R"];
  double EpsilonSGLD = TuningObj_List["EpsilonSGLD"];
  int NSims = TuningObj_List["NSims"];
  int NKeep = TuningObj_List["NKeep"];
  int NTotal = TuningObj_List["NTotal"];
  int BarLength = TuningObj_List["BarLength"];
  arma::vec WhichKeep = TuningObj_List["WhichKeep"];
  arma::vec WhichMAPProgress = TuningObj_List["WhichMAPProgress"];
  arma::vec WhichMAPProgressInt = TuningObj_List["WhichMAPProgressInt"];
  arma::vec WhichSamplerProgress = TuningObj_List["WhichSamplerProgress"];
  arma::vec WhichSamplerProgressInt = TuningObj_List["WhichSamplerProgressInt"];
  arma::vec WhichSGLDProgress = TuningObj_List["WhichSGLDProgress"];
  arma::vec WhichSGLDProgressInt = TuningObj_List["WhichSGLDProgressInt"];
  
  //Convert to C++ struct
  tuning TuningObj;
  TuningObj.EpsilonNADAM = EpsilonNADAM;
  TuningObj.MuNADAM = MuNADAM;
  TuningObj.AlphaNADAM = AlphaNADAM;
  TuningObj.NuNADAM = NuNADAM;
  TuningObj.MNADAM = MNADAM;
  TuningObj.NNADAM = NNADAM;
  TuningObj.S = S;
  TuningObj.S_SGLD = S_SGLD;
  TuningObj.NEpochs = NEpochs;
  TuningObj.R = R;
  TuningObj.EpsilonSGLD = EpsilonSGLD;
  TuningObj.NSims = NSims;
  TuningObj.NKeep = NKeep;
  TuningObj.NTotal = NTotal;
  TuningObj.BarLength = BarLength;
  TuningObj.WhichKeep = WhichKeep;
  TuningObj.WhichMAPProgress = WhichMAPProgress;
  TuningObj.WhichMAPProgressInt = WhichMAPProgressInt;
  TuningObj.WhichSamplerProgress = WhichSamplerProgress;
  TuningObj.WhichSamplerProgressInt = WhichSamplerProgressInt;
  TuningObj.WhichSGLDProgress = WhichSGLDProgress;
  TuningObj.WhichSGLDProgressInt = WhichSGLDProgressInt;
  return TuningObj;
  
}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
para ConvertPara(Rcpp::List Para_List) {

  //Set objects from List
  arma::colvec Beta = Para_List["Beta"];
  arma::colvec l = Para_List["l"];
  arma::colvec d = Para_List["d"];
  arma::colvec Omega = Para_List["Omega"];
  arma::mat Z = Para_List["Z"];
  arma::mat L = Para_List["L"];
  arma::vec LoverZ = Para_List["LoverZ"];
  arma::mat D = Para_List["D"];
  arma::mat Upsilon = Para_List["Upsilon"];
  arma::mat Sigma = Para_List["Sigma"];
  arma::mat LInv = Para_List["LInv"];
  arma::mat tLInv = Para_List["tLInv"];
  arma::mat UpsilonInv = Para_List["UpsilonInv"];
  arma::mat DInv = Para_List["DInv"];
  arma::mat SigmaInv = Para_List["SigmaInv"];
  arma::mat GradLZ = Para_List["GradLZ"];
  arma::mat GradZl = Para_List["GradZl"];
  arma::mat GradLl = Para_List["GradLl"];
  arma::mat SigmaPrime = Para_List["SigmaPrime"];
    
  //Convert to C++ struct
  para Para;
  Para.Beta = Beta;
  Para.l = l;
  Para.d = d;
  Para.Omega = Omega;
  Para.Z = Z;
  Para.L = L;
  Para.LoverZ = LoverZ;
  Para.D = D;
  Para.Upsilon = Upsilon;
  Para.Sigma = Sigma;
  Para.LInv = LInv;
  Para.tLInv = tLInv;
  Para.UpsilonInv = UpsilonInv;
  Para.DInv = DInv;
  Para.SigmaInv = SigmaInv;
  Para.GradLZ = GradLZ;
  Para.GradZl = GradZl;
  Para.GradLl = GradLl;
  Para.SigmaPrime = SigmaPrime;
  return Para;
}
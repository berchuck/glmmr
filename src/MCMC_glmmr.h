//
//  functions used in glmmr package
//

#ifndef __glmmr__
#define __glmmr__

//GLMMR Function
Rcpp::List glmmr_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                      Rcpp::List TuningObj_List, Rcpp::List Para_List,
                      bool Interactive);

//STRUCT DEFINITIONS
struct datobj {
  arma::colvec Y;
  arma::mat X;
  arma::mat Z;
  arma::Col<int> Group;
  arma::Col<int> Group2;
  int N;
  int NUnits;
  int P;
  int Q;
  int NL;
  int FamilyInd;
  int AlgorithmInd;
  arma::mat EyeQ;
  arma::mat EyeNOmega;
  arma::Col<int> SeqNUnits;
  arma::colvec ProbNUnits;
  int NOmega;
};
struct hypara {
  double Eta;
  double Nu;
};
struct tuning {
  double EpsilonNADAM;
  double MuNADAM;
  double AlphaNADAM;
  double NuNADAM;
  arma::vec MNADAM;
  arma::vec NNADAM;
  int S;
  int NEpochs;
  int R;
  double EpsilonSGLD;
  int NSims;
  int NTotal;
  int NKeep;
  int BarLength;
  arma::vec WhichKeep;
  arma::vec WhichMAPProgress;
  arma::vec WhichMAPProgressInt;
  arma::vec WhichSamplerProgress;
  arma::vec WhichSamplerProgressInt;
  arma::vec WhichSGLDProgress;
  arma::vec WhichSGLDProgressInt;
};
struct para {
  arma::colvec Beta;
  arma::colvec l;
  arma::colvec d;
  arma::colvec Omega;
  arma::mat Z;
  arma::mat L;
  arma::vec LoverZ;
  arma::mat D;
  arma::mat Upsilon;
  arma::mat Sigma;
  arma::mat LInv;
  arma::mat tLInv;
  arma::mat UpsilonInv;
  arma::mat DInv;
  arma::mat SigmaInv;
  arma::mat GradLZ;
  arma::mat GradZl;
  arma::mat GradLl;
  arma::mat SigmaPrime;
};

//DISTRIBUTION FUNCTIONS
arma::vec rnormRcpp(int n, double mean, double sd);
arma::vec sampleRcpp(arma::Col<int> const& x, int size, bool replace, arma::vec const& prob);
double rtnormRcppMSM(double mean, double sd, double lower, double upper);
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
arma::colvec rmvnormRcppRobust(arma::colvec const& Mu, arma::mat const& Sigma);
double pnormRcpp(double q);
double lpnormRcpp(double q);
double UpperpnormRcpp(double q);
double lUpperpnormRcpp(double q);
double rigammaRcpp(double Alpha, double Theta);
double rgammaRcpp(double Alpha, double Theta);
arma::mat rwishRcpp(double n, arma::mat const& V);
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti);
double randuRcpp();
double rtnormRcpp(double mean, double sd, bool Above);
arma::vec rtnormRcppMSM(int N, arma::vec const& mean, arma::vec const& sd, double lower, double upper);
arma::vec pgRcpp(arma::vec const& b, arma::vec const& c);

//CONVERSION FUNCTIONS
datobj ConvertDatObj(Rcpp::List DatObj_List);
hypara ConvertHyPara(Rcpp::List HyPara_List);
tuning ConvertTuningObj(Rcpp::List TuningObj_List);
para ConvertPara(Rcpp::List Para_List);

//UTILITY FUNCTIONS
para ComputeSGLDCorrection(datobj DatObj, tuning TuningObj, para Para, bool Interactive);
para UpdatePara(datobj DatObj, para Para);
std::pair<para, tuning> UpdateOmega(int e, arma::colvec const& Grad, datobj DatObj, tuning TuningObj, para Para);
arma::rowvec get_grad_1_L(arma::mat const& L, int q);
arma::mat get_grad_w_L(arma::mat const& L, arma::colvec const& v, arma::colvec const& w, int n_L, int q);
arma::mat ComputeSigmaHatI(int Id, datobj DatObj, tuning TuningObj, para Para);
arma::mat ComputeSGLDCorrection(datobj DatObj, tuning TuningObj, para Para);
arma::colvec ComputeGradientPrior(datobj DatObj, hypara HyPara, para Para);
arma::colvec ComputeGradientI(int Id, arma::mat const& Gamma_i, datobj DatObj, tuning TuningObj, para Para);
arma::mat SampleGamma(int id, datobj DatObj, tuning TuningObj, para Para);
arma::colvec vecLT(arma::mat const& x);
arma::mat GetL(arma::mat const& Z, int Q);
arma::mat GetZ(arma::vec const& l, int Q);
arma::mat CholInv(arma::mat const& Cov);
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol);

//PROGRESS BARS
void BeginMAPProgress(tuning TuningObj, bool Interactive);
void UpdateMAPBar(int e, tuning TuningObj);
void UpdateMAPBarInt(int e, tuning TuningObj);
void BeginSamplerProgress(tuning TuningObj, bool Interactive);
void UpdateSamplerBarInt(int e, tuning TuningObj);
void UpdateSamplerBar(int e, tuning TuningObj);
void BeginSGLDProgress(tuning TuningObj, bool Interactive);
void UpdateSGLDBarInt(int e, tuning TuningObj);
void UpdateSGLDBar(int e, tuning TuningObj);

#endif // __glmmr__

//
//  functions used in glmmr package
//

#ifndef __glmmr__
#define __glmmr__

//SGD Function
arma::mat glmmr_sgd_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                         Rcpp::List SgdObj_List, Rcpp::List Para_List,
                         bool Interactive);

//STRUCT DEFINITIONS
struct datobj {
  arma::colvec Y;
  arma::mat X;
  arma::mat Z;
  arma::Col<int> group;
  arma::Col<int> group2;
  int N;
  int n;
  int p;
  int q;
  int n_L;
  int FamilyInd;
  arma::mat EyeQ;
  arma::Col<int> Seqn;
  arma::colvec Probn;
  int n_omega;
};
struct hypara {
  double Eta;
  double Nu;
};
struct sgdobj {
  double Epsilon;
  arma::vec M_nadam;
  arma::vec N_nadam;
  double Mu_nadam;
  double Alpha_nadam;
  double Nu_nadam;
  double S;
  double n_epochs;
  double R;
  int BarLength;
  arma::vec WhichBurnInProgressInt;
  arma::vec WhichBurnInProgress;
};
struct para {
  arma::colvec Beta;
  arma::colvec l;
  arma::colvec d;
  arma::colvec Omega;
  arma::mat z;
  arma::mat L;
  arma::vec Loverz;
  arma::mat D;
  arma::mat Upsilon;
  arma::mat Sigma;
  arma::mat LInv;
  arma::mat tLInv;
  arma::mat UpsilonInv;
  arma::mat DInv;
  arma::mat SigmaInv;
  arma::mat gradLz;
  arma::mat gradzl;
  arma::mat gradLl;
};

//DISTRIBUTION FUNCTIONS
arma::vec rnormRcpp(int n, double mean, double sd);
arma::vec sampleRcpp(arma::Col<int> const& x, int size, bool replace, arma::vec const& prob);
double rtnormRcppMSM(double mean, double sd, double lower, double upper);
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
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
sgdobj ConvertSgdObj(Rcpp::List SgdObj_List);
para ConvertPara(Rcpp::List Para_List);

//UTILITY FUNCTIONS
std::pair<para, sgdobj> UpdateOmega(arma::colvec const& grad, datobj DatObj, hypara HyPara, sgdobj SgdObj, para Para);
arma::colvec ComputeGradientPrior(datobj DatObj, hypara HyPara, para Para);
arma::colvec ComputeGradienti(int id, arma::mat const& Gamma_i, datobj DatObj, hypara HyPara, sgdobj SgdObj, para Para);
arma::mat SampleGamma(int id, datobj DatObj, hypara HyPara, sgdobj SgdObj, para Para);
arma::colvec vecLT(arma::mat const& x);
arma::mat GetL(arma::mat const& z, int q);
arma::mat GetZ(arma::vec const& l, int q);
arma::mat CholInv(arma::mat const& Cov);
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol);
void UpdateBurnInBar(int e, sgdobj SgdObj);
void UpdateBurnInBarInt(int e, sgdobj SgdObj);
void BeginBurnInProgress(sgdobj SgdObj, bool Interactive);

#endif // __glmmr__

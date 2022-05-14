#include <RcppArmadillo.h>
#include "MCMC_glmmr.h"

//Function to sample latent polya-gamma process using Gibbs sampling step------------------------------------------------
para Sampleomega(datobj DatObj, para Para) {

  //Set data objects
  int NUnits = DatObj.NUnits;
  arma::mat X = DatObj.X;
  arma::field<arma::mat> Z = DatObj.Z;
  arma::Col<int> Group = DatObj.Group;
  
  //Set parameters
  arma::colvec omega = Para.omega;
  arma::mat Gamma = Para.Gamma;
  arma::colvec Beta = Para.Beta;
  
  //Loop over unit index
  for (arma::uword i = 0; i < NUnits; i++) {
    
    //Prepare data
    arma::uvec indeces_row = arma::find(Group == i);
    arma::mat x_i = X.rows(indeces_row);
    arma::mat z_i = Z(i);
    arma::colvec theta_i = x_i * Beta + z_i * Gamma.col(i);
    int n_i = x_i.n_rows;
    
    //Sample latent Variable from full conditional
    omega(indeces_row) = pgRcpp(arma::ones(n_i), theta_i);
    
  //End loop over units 
  }
  
  //Return parameters
  Para.omega = omega;
  return Para;
  
}



//Function to sample random effect------------------------------------------------
para SampleBeta(datobj DatObj, para Para) {
  
  //Set data objects
  int NUnits = DatObj.NUnits;
  int P = DatObj.P;
  arma::colvec Y = DatObj.Y;
  arma::mat X = DatObj.X;
  arma::field<arma::mat> Z = DatObj.Z;
  arma::Col<int> Group = DatObj.Group;

  //Set parameters
  arma::colvec omega = Para.omega;
  arma::mat Gamma = Para.Gamma;

  //Compute moment sums
  arma::mat sum1(P, P, arma::fill::zeros);
  arma::colvec sum2(P, arma::fill::zeros);
  for (arma::uword i = 0; i < NUnits; i++) {
    arma::uvec indeces_row = arma::find(Group == i);
    arma::mat x_i = X.rows(indeces_row);
    arma::mat z_i = Z(i);
    arma::colvec y_i = Y(indeces_row);
    arma::colvec omega_i = omega(indeces_row);
    arma::colvec ystar = ((y_i - 0.5) / omega_i);
    arma::mat D_omega_i = arma::diagmat(omega_i);
    arma::mat tx_i_D_omega_i = arma::trans(x_i) * D_omega_i;
    sum1 += tx_i_D_omega_i * x_i;
    sum2 += tx_i_D_omega_i * (ystar - z_i * Gamma.col(i));
  }
    
  //Compute moments
  arma::mat cov_beta(P, P);
  bool not_singular = arma::inv_sympd(cov_beta, sum1);
  if (!not_singular) Rcpp::stop("Beta covariance singular");
  arma::colvec mean_beta = cov_beta * sum2;
    
  //Sample beta
  // Beta = rmvnormRcpp(1, mean_beta, cov_beta);
  arma::colvec Beta = rmvnormRcppRobust(mean_beta, cov_beta);
  
  //Return Beta
  Para.Beta = Beta;
  Para.Omega = arma::join_cols(Beta, Para.l, Para.d);
  return Para;
  
}



//Function to sample random effectsp------------------------------------------------
para SampleGammaGibbs(datobj DatObj, para Para) {
  
  //Set data objects
  int NUnits = DatObj.NUnits;
  int Q = DatObj.Q;
  arma::colvec Y = DatObj.Y;
  arma::mat X = DatObj.X;
  arma::field<arma::mat> Z = DatObj.Z;
  arma::Col<int> Group = DatObj.Group;
  arma::Col<int> Group2 = DatObj.Group2;
  
  //Set parameters
  arma::colvec omega = Para.omega;
  arma::mat Gamma = Para.Gamma;
  arma::colvec Beta = Para.Beta;
  arma::mat SigmaInv = Para.SigmaInv;
  
  //Loop over unit index
  for (arma::uword i = 0; i < NUnits; i++) {
    
    //Prepare data
    arma::uvec indeces_row = arma::find(Group == i);
    arma::mat x_i = X.rows(indeces_row);
    arma::colvec x_i_beta = x_i * Beta;
    arma::mat z_i = Z(i);
    arma::colvec y_i = Y(indeces_row);
    arma::colvec omega_i = omega(indeces_row);

    //Compute moments
    arma::mat D_omega_i = arma::diagmat(omega_i);
    arma::mat tz_i_D_omega_i = arma::trans(z_i) * D_omega_i;
    // arma::mat cov_gamma = CholInv(tz_i_D_omega_i * z_i + SigmaInv);
    arma::mat cov_gamma(Q, Q);
    bool not_singular = arma::inv_sympd(cov_gamma, tz_i_D_omega_i * z_i + SigmaInv);
    if (!not_singular) Rcpp::stop("Gamma covariance singular");
    arma::colvec ystar = ((y_i - 0.5) / omega_i);
    arma::colvec mean_gamma = cov_gamma * (tz_i_D_omega_i * (ystar - x_i_beta));
    
    //Sample gamma
    // Gamma.col(r) = rmvnormRcpp(1, mean_gamma, cov_gamma);
    Gamma.col(i) = rmvnormRcppRobust(mean_gamma, cov_gamma);
    
    //End loop over R  
  }
  
  //Return gamma
  Para.Gamma = Gamma;
  return Para;
  
}



//Sample L using a Metropolis step------------------------------------------------
std::pair<para, tuning> UpdateL(datobj DatObj, hypara HyPara, tuning TuningObj, para Para) {
  
  //Set data objects
  int NL = DatObj.NL;
  arma::mat EyeQ = DatObj.EyeQ;
  int Q = DatObj.Q;
  int NUnits = DatObj.NUnits;
  
  //Set hyperparameters
  double Eta = HyPara.Eta;
  
  //Set Metropolis Tuning Objects
  arma::vec MetropL = TuningObj.MetropL;
  arma::vec AcceptanceL = TuningObj.AcceptanceL;
  
  //Set parameter objects
  arma::mat LInv = Para.LInv;
  arma::mat DInv = Para.DInv;
  arma::vec l = Para.l;
  arma::mat Gamma = Para.Gamma;
  arma::mat L = Para.L;
  arma::mat Z = Para.Z;
  
  //Loop over visits
  for (arma::uword k = 0; k < NL; k++) {
  
    //Numerical fix for when L becomes singular
    arma::vec lProposal = l;
    arma::mat LInvProposal(Q, Q);
    arma::mat ZProposal;
    arma::mat LProposal;
    double lkProposal;
    bool not_singular = false;
    while (!not_singular) {
      
      //Sample proposal
      lkProposal = arma::as_scalar(rnormRcpp(1, l(k), sqrt(MetropL(k))));
      lProposal(k) = lkProposal;
      
      //Update L
      ZProposal = GetZ(lProposal, Q);
      LProposal = GetL(ZProposal, Q);
      // Rcpp::Rcout << std::fixed << L.diag() << arma::zeros(Q) << std::endl;
      not_singular = arma::solve(LInvProposal, arma::trimatl(LProposal), EyeQ);

    }
    
    //Random Effect Component (Rooti is the cholesky of the inverse of Sigma)
    arma::mat RootiProposal = LInvProposal * DInv;
    arma::mat Rooti = LInv * DInv;
    double Component1A = 0;
    double Component1B = 0;
    for (arma::uword i = 0; i < NUnits; i++) {
      Component1A += lndMvn(Gamma.col(i), arma::zeros(Q), RootiProposal);
      Component1B += lndMvn(Gamma.col(i), arma::zeros(Q), Rooti);
    }
    double Component1 = Component1A - Component1B;
    
    //Prior components
    double Component21A = 0;
    double Component21B = 0;
    double Component22A = 0;
    double Component22B = 0;
    double sum21A;
    double sum21B;
    double sum22A;
    double sum22B;
    for (arma::uword j = 1; j < Q; j++) {
      sum21A = 0;
      sum21B = 0;
      for (arma::uword r = 0; r < j; r++) {
        sum21A += LProposal(j, r) * LProposal(j, r);
        sum21B += L(j, r) * L(j, r);
      }
      Component21A += ((Q - (j + 1) + 2 * Eta - 2) / 2) * log(1 - sum21A);
      Component21B += ((Q - (j + 1) + 2 * Eta - 2) / 2) * log(1 - sum21B);
    }
    if (Q > 2) {
      for (arma::uword j = 1; j < Q; j++) {
        for (arma::uword i = (j + 1); i < Q; i++) {
          sum22A = 0;
          sum22B = 0;
          for (arma::uword r = 0; r < (j - 1); r++) {
            sum22A += LProposal(i, r) * LProposal(i, r);
            sum22B += L(i, r) * L(i, r);
          }
          Component22A += 0.5 * log(1 - sum22A);
          Component22B += 0.5 * log(1 - sum22B);
        }
      }
    }
    double Component23A = -2 * arma::as_scalar(arma::sum(arma::log(cosh(lProposal))));
    double Component23B = -2 * arma::as_scalar(arma::sum(arma::log(cosh(l))));
    double Component2A = Component21A + Component22A + Component23A;
    double Component2B = Component21B + Component22B + Component23B;
    double Component2 = Component2A - Component2B;
    
    //Log acceptance ratio
    double LogR = Component1 + Component2;
    
    //Metropolis update
    double RandU = randuRcpp();
    if (log(RandU) < LogR) {
      
      //Keep count of acceptances
      AcceptanceL(k)++;
      
      //Update parameters output
      l = lProposal;
      Z = ZProposal;
      L = LProposal;
      LInv = LInvProposal;

    }
    
  //End loop over L entries
  }
  
  //Update Metropolis object
  TuningObj.AcceptanceL = AcceptanceL;
  
  //Update Para objects
  Para.l = l;
  Para.Z = Z;
  Para.L = L;
  Para.LoverZ = vecLT(L / Z);
  Para.Upsilon = L * arma::trans(L);
  Para.Sigma = Para.D * Para.Upsilon * Para.D;
  Para.LInv = LInv;
  Para.tLInv = arma::trans(LInv);
  Para.UpsilonInv = Para.tLInv * LInv;
  Para.SigmaInv = Para.DInv * Para.UpsilonInv * Para.DInv;
  Para.GradLZ = arma::diagmat(Para.LoverZ);
  Para.GradZl = arma::diagmat(arma::pow(arma::cosh(l), -2));
  Para.GradLl = Para.GradLZ * Para.GradZl;
  Para.Omega = arma::join_cols(Para.Beta, Para.l, Para.d);
  
  //Return final object
  return std::pair<para, tuning>(Para, TuningObj);
  
}



//Sample D using a Metropolis step------------------------------------------------
std::pair<para, tuning> UpdateD(datobj DatObj, hypara HyPara, tuning TuningObj, para Para) {
  
  //Set data objects
  arma::mat EyeQ = DatObj.EyeQ;
  int Q = DatObj.Q;
  int NUnits = DatObj.NUnits;
  
  //Set hyperparameter objects
  double Nu = HyPara.Nu;
  
  //Set Metropolis Tuning Objects
  arma::vec MetropD = TuningObj.MetropD;
  arma::vec AcceptanceD = TuningObj.AcceptanceD;
  
  //Set parameter objects
  arma::vec d = Para.d;
  arma::mat LInv = Para.LInv;
  arma::mat Gamma = Para.Gamma;
  arma::mat DInv = Para.DInv;
  
  //Loop over visits
  for (int k = 0; k < Q; k++) {
    
    //Numerical fix for when d becomes too large
    arma::vec dProposal = d;
    arma::mat DInvProposal;
    double dkProposal;
    double sigmakProposal = arma::datum::inf;
    while (!arma::is_finite(sigmakProposal)) {
  
      //Sample proposal
      dkProposal = arma::as_scalar(rnormRcpp(1, d(k), sqrt(MetropD(k))));
      dProposal(k) = dkProposal;
      sigmakProposal = exp(dkProposal);
      
      //Update DInv
      DInvProposal = arma::diagmat(arma::exp(-dProposal));
      
    }
    
    //Random Effect Component (Rooti is the cholesky of the inverse of Sigma)
    arma::mat RootiProposal = LInv * DInvProposal;
    arma::mat Rooti = LInv * DInv;
    double Component1A = 0;
    double Component1B = 0;
    for (arma::uword i = 0; i < NUnits; i++) {
      Component1A += lndMvn(Gamma.col(i), arma::zeros(Q), RootiProposal);
      Component1B += lndMvn(Gamma.col(i), arma::zeros(Q), Rooti);
    }
    double Component1 = Component1A - Component1B;
    
    //Prior components
    double Component2A = dkProposal - 0.5 * (Nu + 1) * log(1 + exp(2 * dkProposal) / Nu);
    double Component2B = d(k) - 0.5 * (Nu + 1) * log(1 + exp(2 * d(k)) / Nu);
    double Component2 = Component2A - Component2B;
    
    // Rcpp::Rcout << std::fixed << Component1 << " " << Component2 << std::endl;
    
    //Log acceptance ratio
    double LogR = Component1 + Component2;
    
    //Metropolis update
    double RandU = randuRcpp();
    if (log(RandU) < LogR) {
      
      //Keep count of acceptances
      AcceptanceD(k)++;
      
      //Update parameters output
      d = dProposal;
      DInv = DInvProposal;
      
    }
    
    //End loop over L entries
  }
  
  //Update Metropolis object
  TuningObj.AcceptanceD = AcceptanceD;
  
  //Update Para objects
  Para.d = d;
  Para.D = arma::diagmat(arma::exp(d));
  Para.Sigma = Para.D * Para.Upsilon * Para.D;
  Para.DInv = DInv;
  Para.SigmaInv = Para.DInv * Para.UpsilonInv * Para.DInv;
  Para.Omega = arma::join_cols(Para.Beta, Para.l, Para.d);
  
  
  //Return final object
  return std::pair<para, tuning>(Para, TuningObj);
  
}



//Function to compute the SGLD correction
para ComputeSGLDCorrection(datobj DatObj, tuning TuningObj, para Para, bool Interactive) {
  
  //Set data objects
  int NUnits = DatObj.NUnits;
  int NOmega = DatObj.NOmega;
  arma::mat EyeNOmega = DatObj.EyeNOmega;
  
  //Set tuning objects
  arma::vec WhichSGLDProgress = TuningObj.WhichSGLDProgress;
  arma::vec WhichSGLDProgressInt = TuningObj.WhichSGLDProgressInt;
  double EpsilonSGLD = TuningObj.EpsilonSGLD;
  int S_SGLD = TuningObj.S_SGLD;
  
  //User output
  BeginSGLDProgress(TuningObj, Interactive);
  
  //Initialize objects
  arma::mat SigmaSum(NOmega, NOmega, arma::fill::zeros);
  arma::uvec Samps = arma::randperm(NUnits, S_SGLD);
  
  //Loop over all units
  for (arma::uword i = 0; i < S_SGLD; i++) {
    SigmaSum += ComputeSigmaHatI(Samps(i), DatObj, TuningObj, Para);
    if (Interactive) if (std::find(WhichSGLDProgress.begin(), WhichSGLDProgress.end(), i) != WhichSGLDProgress.end())
      UpdateSGLDBar(i, TuningObj);
    if (!Interactive) if (std::find(WhichSGLDProgressInt.begin(), WhichSGLDProgressInt.end(), i) != WhichSGLDProgressInt.end())
      UpdateSGLDBarInt(i, TuningObj);
  }
  arma::mat SigmaPrime = arma::chol(SigmaSum) / sqrt(S_SGLD);
  arma::mat Sigma = sqrt(2) * EyeNOmega - sqrt(EpsilonSGLD) * SigmaPrime;
  Para.SigmaPrime = EpsilonSGLD * Sigma * Sigma.t();
  return Para;
  
}



//Function to compute Monte Carlo variance for unit i
arma::mat ComputeSigmaHatI(int Id, datobj DatObj, tuning TuningObj, para Para) {
  
  //Begin by computing gamma and mu hat for subject i
  arma::mat GammaI = SampleGamma(Id, DatObj, TuningObj, Para);
  arma::colvec MuHatI = ComputeGradientI(Id, GammaI, DatObj, TuningObj, Para);
  
  //Set data objects
  int P = DatObj.P;
  int Q = DatObj.Q;
  int NL = DatObj.NL;
  int NOmega = DatObj.NOmega;
  arma::colvec Y = DatObj.Y;
  arma::mat X = DatObj.X;
  arma::field<arma::mat> Z = DatObj.Z;
  arma::Col<int> Group = DatObj.Group;
  arma::Col<int> Group2 = DatObj.Group2;
  
  //Set tuning objects
  int R = TuningObj.R;
  
  //Set parameters
  arma::colvec Beta = Para.Beta;
  arma::mat DInv = Para.DInv;
  arma::mat LInv = Para.LInv;
  arma::mat GradLl = Para.GradLl;
  arma::mat L = Para.L;
  arma::mat UpsilonInv = Para.UpsilonInv;
  arma::colvec d = Para.d;
  
  //Prepare data
  arma::uvec indeces_row = arma::find(Group == Id);
  arma::uvec indeces_col = arma::find(Group2 == Id);
  arma::mat x_i = X.rows(indeces_row);
  arma::colvec x_i_beta = x_i * Beta;
  // arma::mat z_i = Z.submat(indeces_row, indeces_col);
  arma::mat z_i = Z(Id);
  arma::colvec y_i = Y(indeces_row);
  
  //Compute GLMM mean parameter
  arma::mat theta_i = arma::repmat(x_i_beta, 1, R) + z_i * GammaI;
  arma::mat fit_i = arma::exp(theta_i);
  bool any_inf = fit_i.has_inf();
  arma::mat pi_i = fit_i / (1 + fit_i);
  if (any_inf) {
    arma::uvec where_inf = arma::find_nonfinite(fit_i); 
    for (arma::uword i = 0; i < where_inf.size(); i++) {
      arma::uword index = where_inf(i);
      if (fit_i(index) == arma::datum::inf) pi_i(index) = 1;
      if (fit_i(index) == -arma::datum::inf) pi_i(index) = 0;
    }
  }
  
  //Initialize objects
  arma::colvec grad_beta_likelihood(P, arma::fill::zeros);
  arma::colvec grad_L_likelihood(NL, arma::fill::zeros); 
  arma::colvec grad_D_likelihood(Q, arma::fill::zeros); 
  arma::colvec grad_re_d(Q, arma::fill::zeros);
  arma::colvec grad_Omega_r;
  arma::mat SigmaHatOmega(NOmega, NOmega, arma::fill::zeros);
  
  //Compute gradients by looping over R
  for (arma::uword r = 1; r < R; r++) {
    arma::colvec gamma_ir = GammaI.col(r);
    grad_beta_likelihood = arma::trans(arma::trans(y_i - pi_i.col(r)) * x_i);
    arma::colvec v = DInv * gamma_ir;
    arma::colvec w = LInv * v;
    grad_L_likelihood = (get_grad_1_L(L, Q) - arma::trans(w) * get_grad_w_L(L, v, w, NL, Q)) * GradLl;
    arma::mat Gamma_r = arma::diagmat(gamma_ir);
    arma::mat M = Gamma_r * UpsilonInv * Gamma_r;
    for (arma::uword k = 0; k < Q; k++) {
      double sum1 = 0;
      for (arma::uword h = 0; h < Q; h++) sum1 += M(h, k) * exp(-d(h));
      grad_re_d(k) = -1 + exp(-d(k)) * sum1;
    }
    grad_D_likelihood = grad_re_d;
    grad_Omega_r = arma::join_cols(grad_beta_likelihood, grad_L_likelihood, grad_D_likelihood);
    SigmaHatOmega += (grad_Omega_r - MuHatI) * arma::trans(grad_Omega_r - MuHatI);
  }
  
  //Final likelihood contribution
  SigmaHatOmega /= ((R - 1) * (R - 1));
  return SigmaHatOmega;
  
//End function to compute variance of the Monte Carlo integral       
} 



//Function to update parameters using NADAM
para UpdatePara(datobj DatObj, para Para) {
  
  //Set data objects
  int P = DatObj.P;
  int NL = DatObj.NL;
  int Q = DatObj.Q;
  arma::mat EyeQ = DatObj.EyeQ;
  
  //Set parameters
  arma::colvec Omega = Para.Omega;
  
  //Update other parameters
  arma::colvec Beta = Omega(arma::span(0, P - 1));
  arma::colvec l = Omega(arma::span(P, P + NL - 1));
  arma::colvec d = Omega(arma::span(P + NL, P + NL + Q - 1));
  
  //Update L
  arma::mat Z = GetZ(l, Q);
  arma::mat L = GetL(Z, Q);
  // Rcpp::Rcout << std::fixed << L.diag() << arma::zeros(Q) << std::endl;
  bool close_to_singular = any(L.diag() == 0);
  // bool close_to_singular = approx_equal(L.diag(), arma::zeros(Q), "absdiff", 0.0001);
  if (close_to_singular) Rcpp::stop("L is almost singular. Consider decreasing EpsilonSGLD.");
  arma::mat LInv(Q, Q);
  bool not_singular = arma::solve(LInv, arma::trimatl(L), EyeQ);
  if (!not_singular) Rcpp::stop("L is singular. Decrease EpsilonSGLD.");
  
  //Save parameters
  Para.Omega = Omega;
  Para.Beta = Beta;
  Para.l = l;
  Para.d = d;
  Para.Z = Z;
  Para.L = L;
  Para.LoverZ = vecLT(L / Z);
  Para.D = arma::diagmat(arma::exp(d));
  Para.Upsilon = L * arma::trans(L);
  Para.Sigma = Para.D * Para.Upsilon * Para.D;
  Para.LInv = LInv;
  Para.tLInv = arma::trans(LInv);
  Para.UpsilonInv = Para.tLInv * LInv;
  Para.DInv = arma::diagmat(arma::exp(-d));
  Para.SigmaInv = Para.DInv * Para.UpsilonInv * Para.DInv;
  Para.GradLZ = arma::diagmat(Para.LoverZ);
  Para.GradZl = arma::diagmat(arma::pow(arma::cosh(l), -2));
  Para.GradLl = Para.GradLZ * Para.GradZl;
  
  //Return updated object
  return Para;

}



//Function to update Omega
std::pair<para, tuning> UpdateOmega(int e, arma::colvec const& Grad, datobj DatObj, tuning TuningObj, para Para) {
  
  //Set data objects
  int NOmega = DatObj.NOmega;
  int AlgorithmInd = DatObj.AlgorithmInd;

  //Set tuning objects
  double EpsilonNADAM = TuningObj.EpsilonNADAM;
  double MuNADAM = TuningObj.MuNADAM;
  double NuNADAM = TuningObj.NuNADAM;
  double AlphaNADAM = TuningObj.AlphaNADAM;
  arma::vec MNADAM = TuningObj.MNADAM;
  arma::vec NNADAM = TuningObj.NNADAM;
  double EpsilonSGLD = TuningObj.EpsilonSGLD;
  int NEpochs = TuningObj.NEpochs;
  
  //Set parameters
  arma::colvec Omega = Para.Omega;
  
  //Update omega using NADAM
  if (e < NEpochs) {
    double mhat, nhat;
    for (arma::uword i = 0; i < NOmega; i++) {
      MNADAM(i) = MuNADAM * MNADAM(i) + (1 - MuNADAM) * Grad(i);
      NNADAM(i) = NuNADAM * NNADAM(i) + (1 - NuNADAM) * pow(Grad(i), 2);
      mhat = (MuNADAM * MNADAM(i) / (1 - MuNADAM)) + ((1 - MuNADAM) * Grad(i) / (1 - MuNADAM));
      nhat = NuNADAM * NNADAM(i) / (1 - NuNADAM);
      Omega(i) = Omega(i) + AlphaNADAM / (sqrt(nhat) + EpsilonNADAM) * mhat;
    }
  }
  
  //Update omega using SGLD with or without the correction
  if (AlgorithmInd > 0) {
    if (e >= NEpochs) {
      if (AlgorithmInd == 1) Omega += (0.5 * EpsilonSGLD * Grad + rnormRcpp(NOmega, 0, sqrt(EpsilonSGLD))); // SGLD from Welling et al. 2011
      if (AlgorithmInd == 2) Omega += EpsilonSGLD * Grad + rmvnormRcpp(1, arma::zeros(NOmega), Para.SigmaPrime); // SGLD with correction
    }
  }
  //Update parameter object
  Para.Omega = Omega;
  
  //Update tuning object
  if (e < NEpochs) {
    TuningObj.MNADAM = MNADAM;
    TuningObj.NNADAM = NNADAM;
  }
  
  //Return updated Omega
  return std::pair<para, tuning>(Para, TuningObj);
  
}



//Function to compute gradient computation for prior
arma::colvec ComputeGradientPrior(datobj DatObj, hypara HyPara, para Para) {
  
  //Set data objects
  int Q = DatObj.Q;
  int P = DatObj.P;

  //Set hyperparameters
  double Eta = HyPara.Eta;
  double Nu = HyPara.Nu;
  
  //Set parameters
  arma::mat L = Para.L;
  arma::mat GradLl = Para.GradLl;
  arma::mat Z = Para.Z;
  arma::colvec d = Para.d;
  
  //Prior contribution for Beta
  arma::colvec grad_beta_prior(P, arma::fill::zeros);
  
  //Prior contribution for l
  arma::mat grad_1_L(Q, Q, arma::fill::zeros), grad_2_L(Q, Q, arma::fill::zeros);
  for (arma::uword i = 0; i < Q; i++) {
    for (arma::uword j = 0; j < Q; j++) {
      if (i > j) {
        double sum1 = 0;
        for (arma::uword h = 0; h < i; h++) sum1 += L(i, h) * L(i, h);
        grad_1_L(i, j) = -((Q - (i + 1) + 2 * Eta - 2) * L(i, j)) / (1 - sum1);
        if (i > 1) {
          if (j < (i - 1)) {
            double sum2 = 0;
            for (arma::uword r = j; r < (i - 1); r++) {
              double sum3 = 0;
              for (arma::uword k = 0; k < (r + 1); k++) sum3 += L(i, k) * L(i, k);
              sum2 += 1 / (1 - sum3);
            }
            grad_2_L(i, j) = -L(i, j) * sum2;
          }
        }
      }
    }
  }
  arma::colvec grad_1_l = vecLT(grad_1_L) * GradLl;
  arma::colvec grad_2_l = vecLT(grad_2_L) * GradLl;
  arma::colvec grad_3_l = -2 * arma::diagmat(vecLT(Z));
  arma::colvec grad_l_prior = grad_1_l + grad_2_l + grad_3_l;
  
  //Prior contribution for d
  arma::colvec grad_d_prior(Q);
  for (arma::uword k = 0; k < Q; k++) grad_d_prior(k) = (Nu - Nu * exp(2 * d(k))) / (Nu + exp(2 * d(k)));

  //Final prior gradient contribution
  arma::colvec grad_prior = arma::join_cols(grad_beta_prior, grad_l_prior, grad_d_prior);
  return grad_prior;
  
}



//Gradient computation
arma::rowvec get_grad_1_L(arma::mat const& L, int q) {
  arma::mat grad_1_L(q, q, arma::fill::zeros);
  for (arma::uword i = 0; i < q; i++) {
    for (arma::uword j = 0; j < q; j++) {
      if (i > j) {
        double sum1 = 0;
        for (arma::uword k = 0; k < i; k++) sum1 += L(i, k) * L(i, k);
        grad_1_L(i, j) = L(i, j) / (1 - sum1);
      }
    }
  }
  return arma::trans(vecLT(grad_1_L));
}



//Gradient computation
arma::mat get_grad_w_L(arma::mat const& L, arma::colvec const& v, arma::colvec const& w, int n_L, int q) {
  arma::cube grad_w_l_k(q, q, q);
  for (arma::uword k = 0; k < q; k++) {
    for (arma::uword i = 0; i < q; i++) {
      for (arma::uword j = 0; j < q; j++) {
        if (i <= j) grad_w_l_k(k, i, j) = 0;
        if (i > j) {
          if (k == 0) grad_w_l_k(k, i, j) = 0;
          if (k > 0) {
            if (i == k) {
              double sum1 = 0, sum2 = 0;
              for (arma::uword h = 0; h < i; h++) sum1 += L(i, h) * w(h);
              for (arma::uword h = 0; h < i; h++) sum2 += L(i, h) * L(i, h);
              double component1 = v(i) - sum1;
              double component2 = pow(1 - sum2, -1.5);
              double component3 = w(j) * pow(1 - sum2, -0.5);
              grad_w_l_k(k, i, j) = L(i, j) * component1 * component2 - component3;
            }
            if (i > k) grad_w_l_k(k, i, j) = 0;
            if (i < k) {
              int c_diff = k - i;
              double sum3 = 0, sum4 = 0;
              for (arma::uword h = 0; h < (i + c_diff); h++) sum3 += L(i + c_diff, h) * L(i + c_diff, h);
              for (arma::uword r = i; r < (i + c_diff); r++) sum4 += L(i + c_diff, r) * grad_w_l_k(r, i, j);
              double component1 = pow(1 - sum3, -0.5);
              double component2 = -sum4;
              grad_w_l_k(k, i, j) = component1 * component2;
            }
            
          }
        }
      }
    }
  }
  arma::mat out(n_L, q);
  for (arma::uword k = 0; k < q; k++) {
    out.col(k) = vecLT(grad_w_l_k.row(k));
  }
  return arma::trans(out);
}



//Function to compute likelihood gradient component for unit i
arma::colvec ComputeGradientI(int Id, arma::mat const& GammaI, datobj DatObj, tuning TuningObj, para Para) {

  //Set data objects
  int P = DatObj.P;
  int Q = DatObj.Q;
  int NL = DatObj.NL;
  arma::colvec Y = DatObj.Y;
  arma::mat X = DatObj.X;
  arma::field<arma::mat> Z = DatObj.Z;
  arma::Col<int> Group = DatObj.Group;
  arma::Col<int> Group2 = DatObj.Group2;

  //Set tuning objects
  int R = TuningObj.R;
  
  //Set parameters
  arma::colvec Beta = Para.Beta;
  arma::mat DInv = Para.DInv;
  arma::mat LInv = Para.LInv;
  arma::mat GradLl = Para.GradLl;
  arma::mat L = Para.L;
  arma::mat UpsilonInv = Para.UpsilonInv;
  arma::colvec d = Para.d;
  
  //Prepare data
  arma::uvec indeces_row = arma::find(Group == Id);
  arma::uvec indeces_col = arma::find(Group2 == Id);
  arma::mat x_i = X.rows(indeces_row);
  arma::colvec x_i_beta = x_i * Beta;
  // arma::mat z_i = Z.submat(indeces_row, indeces_col);
  arma::mat z_i = Z(Id);
  arma::colvec y_i = Y(indeces_row);
  
  //Compute GLMM mean parameter
  arma::mat theta_i = arma::repmat(x_i_beta, 1, R) + z_i * GammaI;
  arma::mat fit_i = arma::exp(theta_i);
  bool any_inf = fit_i.has_inf();
  arma::mat pi_i = fit_i / (1 + fit_i);
  if (any_inf) {
    arma::uvec where_inf = arma::find_nonfinite(fit_i); 
    for (arma::uword i = 0; i < where_inf.size(); i++) {
      arma::uword index = where_inf(i);
      if (fit_i(index) == arma::datum::inf) pi_i(index) = 1;
      if (fit_i(index) == -arma::datum::inf) pi_i(index) = 0;
    }
  }
  
  //Initialize objects
  arma::colvec grad_beta_likelihood(P, arma::fill::zeros);
  arma::colvec grad_L_likelihood(NL, arma::fill::zeros); 
  arma::colvec grad_D_likelihood(Q, arma::fill::zeros); 
  arma::colvec grad_re_d(Q, arma::fill::zeros); 
    
  //Compute gradients by looping over R
  for (arma::uword r = 1; r < R; r++) {
    arma::colvec gamma_ir = GammaI.col(r);
    grad_beta_likelihood += arma::trans(arma::trans(y_i - pi_i.col(r)) * x_i);
    arma::colvec v = DInv * gamma_ir;
    arma::colvec w = LInv * v;
    grad_L_likelihood += (get_grad_1_L(L, Q) - arma::trans(w) * get_grad_w_L(L, v, w, NL, Q)) * GradLl;
    arma::mat Gamma_r = arma::diagmat(gamma_ir);
    arma::mat M = Gamma_r * UpsilonInv * Gamma_r;
    for (arma::uword k = 0; k < Q; k++) {
      double sum1 = 0;
      for (arma::uword h = 0; h < Q; h++) sum1 += M(h, k) * exp(-d(h));
      grad_re_d(k) = -1 + exp(-d(k)) * sum1;
    }
    grad_D_likelihood += grad_re_d;
  }
  
  //Final likelihood contribution
  grad_beta_likelihood /= (R - 1);
  grad_L_likelihood /= (R - 1);
  grad_D_likelihood /= (R - 1);
  arma::colvec grad_likelihood_i = arma::join_cols(grad_beta_likelihood, grad_L_likelihood, grad_D_likelihood);
  return grad_likelihood_i;
  
//End function to compute gradient component       
} 



//Function to sample random effects
arma::mat SampleGamma(int Id, datobj DatObj, tuning TuningObj, para Para) {

  //Set data objects
  int Q = DatObj.Q;
  arma::colvec Y = DatObj.Y;
  arma::mat X = DatObj.X;
  arma::field<arma::mat> Z = DatObj.Z;
  arma::Col<int> Group = DatObj.Group;
  arma::Col<int> Group2 = DatObj.Group2;
  
  //Set tuning objects
  int R = TuningObj.R;
  
  //Set parameters
  arma::colvec Beta = Para.Beta;
  arma::mat SigmaInv = Para.SigmaInv;
  
  //Prepare data
  arma::uvec indeces_row = arma::find(Group == Id);
  arma::uvec indeces_col = arma::find(Group2 == Id);
  arma::mat x_i = X.rows(indeces_row);
  arma::colvec x_i_beta = x_i * Beta;
  // arma::mat z_i = Z.submat(indeces_row, indeces_col);
  arma::mat z_i = Z(Id);
  arma::colvec y_i = Y(indeces_row);
  int n_i = x_i.n_rows;
  
  //Initialize output object
  arma::mat Gamma(Q, R, arma::fill::zeros);
  
  //Loop over R
  for (arma::uword r = 1; r < R; r++) {
    
    // Sample omega
    arma::colvec omega_i = pgRcpp(arma::ones(n_i), x_i_beta + z_i * Gamma.col(r - 1));
    
    //Compute moments
    arma::mat D_omega_i = arma::diagmat(omega_i);
    arma::mat tz_i_D_omega_i = arma::trans(z_i) * D_omega_i;
    // arma::mat cov_gamma = CholInv(tz_i_D_omega_i * z_i + SigmaInv);
    arma::mat cov_gamma(Q, Q);
    bool not_singular = arma::inv_sympd(cov_gamma, tz_i_D_omega_i * z_i + SigmaInv);
    if (!not_singular) Rcpp::stop("Gamma covariance singular");
    arma::colvec ystar = ((y_i - 0.5) / omega_i);
    arma::colvec mean_gamma = cov_gamma * (tz_i_D_omega_i * (ystar - x_i_beta));
    
    //Sample gamma
    // Gamma.col(r) = rmvnormRcpp(1, mean_gamma, cov_gamma);
    Gamma.col(r) = rmvnormRcppRobust(mean_gamma, cov_gamma);
    
  //End loop over R  
  }
  
  //Return gamma
  return Gamma;
  
}



//Function to get z matrix-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetZ(arma::vec const& l, int Q) {
  arma::mat z(Q, Q, arma::fill::zeros);
  arma::uvec upper_indices = arma::trimatu_ind(arma::size(z), 1);
  z.elem(upper_indices) += arma::tanh(l);
  return z.t();
}



//Function to get L matrix-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetL(arma::mat const& Z, int Q) {
  arma::mat L(Q, Q, arma::fill::zeros);
  for (arma::uword i = 0; i < Q; i++) {
    for (arma::uword j = 0; j < Q; j++) {
      if (j == 0) {
        if (i == 0) L(i, j) = 1;
        if (i > j) L(i, j) = Z(i, j);
      }
      if (j > 0) {
        arma::rowvec Li = L.row(i);
        if (i > j) L(i, j) = Z(i, j) * sqrt(1 - arma::as_scalar(arma::sum(arma::pow(Li.subvec(0, j - 1), 2))));         
        if (i == j) L(i, j) = sqrt(1 - arma::as_scalar(arma::sum(arma::pow(Li.subvec(0, j - 1), 2))));         
      }
    }
  }
  return L;
}

  
  
//Function that vectorizes the lower triangular components of a matrix row-wise
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec vecLT(arma::mat const& x) {
  arma::mat B = x.t();
  arma::uvec upper_indices = arma::trimatu_ind(arma::size(x), 1);
  return arma::vectorise(B.elem(upper_indices));
}  
  
//Matrix inverse using cholesky decomoposition for covariances-------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CholInv(arma::mat const& Cov) {
  return arma::inv_sympd(Cov);
}



//Matrix inverse of for 3x3 matrix-----------------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat Inv3(arma::mat const& A) {
  arma::mat result = arma::mat(3, 3);
  double determinant = A(0, 0) * ( A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2) ) - A(0, 1) * ( A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0) ) + A(0, 2) * ( A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0) );
  double invdet = 1 / determinant;
  result(0, 0) =  ( A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2) ) * invdet;
  result(1, 0) = -( A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1) ) * invdet;
  result(2, 0) =  ( A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1) ) * invdet;
  result(0, 1) = -( A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0) ) * invdet;
  result(1, 1) =  ( A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0) ) * invdet;
  result(2, 1) = -( A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2) ) * invdet;
  result(0, 2) =  ( A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1) ) * invdet;
  result(1, 2) = -( A(0, 0) * A(2, 1) - A(2, 0) * A(0, 1) ) * invdet;
  result(2, 2) =  ( A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1) ) * invdet;
  return result;
}



//Matrix inverse of for 2x2 matrix-----------------------------------------------------------------------
arma::mat Inv2(arma::mat const& A) {
  arma::mat result = arma::mat(2, 2);
  double determinant = ( A(0, 0) * A(1, 1) ) - ( A(0, 1) * A(0, 1) );
  double invdet = 1 / determinant;
  result(0, 0) =  A(1, 1) * invdet;
  result(0, 1) = -A(1, 0) * invdet;
  result(1, 0) = -A(0, 1) * invdet;
  result(1, 1) =  A(0, 0) * invdet;
  return result;
}



//Function for making an upper diagonal matrix symmetric-------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat makeSymm(arma::mat const& A) {
  return arma::symmatu(A);
}



//Function that checks numerical equality of two objects against a tolerance-----------------------------
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}

#include <RcppArmadillo.h>
#include "MCMC_glmmr.h"

//Function to update parameters using NADAM
std::pair<para, sgdobj> UpdateOmega(arma::colvec const& grad, datobj DatObj, hypara HyPara, sgdobj SgdObj, para Para) {

  //Set data objects
  int p = DatObj.p;
  int n_L = DatObj.n_L;
  int q = DatObj.q;
  int n_omega = DatObj.n_omega;
  arma::mat EyeQ = DatObj.EyeQ;
  
  //Set hyperparameters

  //Set sgd objects
  arma::vec M_nadam = SgdObj.M_nadam;
  arma::vec N_nadam = SgdObj.N_nadam;
  double Nu_nadam = SgdObj.Nu_nadam;
  double Mu_nadam = SgdObj.Mu_nadam;
  double Alpha_nadam = SgdObj.Alpha_nadam;
  double Epsilon = SgdObj.Epsilon;

  //Set parameters
  arma::colvec Omega = Para.Omega;
  
  //Update omega using NADAM
  double mhat, nhat;
  for (arma::uword i = 0; i < n_omega; i++) {
    M_nadam(i) = Mu_nadam * M_nadam(i) + (1 - Mu_nadam) * grad(i);
    N_nadam(i) = Nu_nadam * N_nadam(i) + (1 - Nu_nadam) * pow(grad(i), 2);
    mhat = (Mu_nadam * M_nadam(i) / (1 - Mu_nadam)) + ((1 - Mu_nadam) * grad(i) / (1 - Mu_nadam));
    nhat = Nu_nadam * N_nadam(i) / (1 - Nu_nadam);
    Omega(i) = Omega(i) + Alpha_nadam / (sqrt(nhat) + Epsilon) * mhat;
  }

  //Update other parameters
  arma::colvec Beta = Omega(arma::span(0, p - 1));
  arma::colvec l = Omega(arma::span(p, p + n_L - 1));
  arma::colvec d = Omega(arma::span(p + n_L, p + n_L + q - 1));
  
  //Save parameters
  Para.Omega = Omega;
  Para.Beta = Beta;
  Para.l = l;
  Para.d = d;
  Para.z = GetZ(l, q);
  Para.L = GetL(Para.z, q);
  Para.Loverz = vecLT(Para.L / Para.z);
  Para.D = arma::diagmat(arma::exp(d));
  Para.Upsilon = Para.L * arma::trans(Para.L);
  Para.Sigma = Para.D * Para.Upsilon * Para.D;
  Para.LInv = arma::solve(arma::trimatl(Para.L), EyeQ);
  Para.tLInv = arma::trans(Para.LInv);
  Para.UpsilonInv = Para.tLInv * Para.LInv;
  Para.DInv = arma::diagmat(arma::exp(-d));
  Para.SigmaInv = Para.DInv * Para.UpsilonInv * Para.DInv;
  Para.gradLz = arma::diagmat(Para.Loverz);
  Para.gradzl = arma::diagmat(arma::pow(arma::cosh(l), -2));
  Para.gradLl = Para.gradLz * Para.gradzl;
  
  //Update sgd object
  SgdObj.M_nadam = M_nadam;
  SgdObj.N_nadam = N_nadam;
  
  //Return updated objects
  return std::pair<para, sgdobj>(Para, SgdObj);
  
}



//Function to compute gradient computation for prior
arma::colvec ComputeGradientPrior(datobj DatObj, hypara HyPara, para Para) {
  
  //Set data objects
  int q = DatObj.q;
  int p = DatObj.p;

  //Set hyperparameters
  double Eta = HyPara.Eta;
  double Nu = HyPara.Nu;
  
  //Set parameters
  arma::mat L = Para.L;
  arma::mat gradLl = Para.gradLl;
  arma::mat z = Para.z;
  arma::colvec d = Para.d;
  
  //Prior contribution for Beta
  arma::colvec grad_beta_prior(p, arma::fill::zeros);
  
  //Prior contribution for l
  arma::mat grad_1_L(q, q, arma::fill::zeros), grad_2_L(q, q, arma::fill::zeros);
  for (arma::uword i = 0; i < q; i++) {
    for (arma::uword j = 0; j < q; j++) {
      if (i > j) {
        double sum1 = 0;
        for (arma::uword h = 0; h < i; h++) sum1 += L(i, h) * L(i, h);
        grad_1_L(i, j) = -((q - (i + 1) + 2 * Eta - 2) * L(i, j)) / (1 - sum1);
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
  arma::colvec grad_1_l = vecLT(grad_1_L) * gradLl;
  arma::colvec grad_2_l = vecLT(grad_2_L) * gradLl;
  arma::colvec grad_3_l = -2 * arma::diagmat(vecLT(z));
  arma::colvec grad_l_prior = grad_1_l + grad_2_l + grad_3_l;
  
  //Prior contribution for d
  arma::colvec grad_d_prior(q);
  for (arma::uword k = 0; k < q; k++) grad_d_prior(k) = (Nu - Nu * exp(2 * d(k))) / (Nu + exp(2 * d(k)));

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
arma::colvec ComputeGradienti(int id, arma::mat const& Gamma_i, datobj DatObj, hypara HyPara, sgdobj SgdObj, para Para) {

  //Set data objects
  int p = DatObj.p;
  int q = DatObj.q;
  int n_L = DatObj.n_L;
  arma::colvec Y = DatObj.Y;
  arma::mat X = DatObj.X;
  arma::mat Z = DatObj.Z;
  arma::Col<int> group = DatObj.group;
  arma::Col<int> group2 = DatObj.group2;
  
  //Set hyperparameters
  
  //Set sgd objects
  int R = SgdObj.R;
  
  //Set parameters
  arma::colvec Beta = Para.Beta;
  arma::mat SigmaInv = Para.SigmaInv;
  arma::mat DInv = Para.DInv;
  arma::mat LInv = Para.LInv;
  arma::mat gradLl = Para.gradLl;
  arma::mat L = Para.L;
  arma::mat UpsilonInv = Para.UpsilonInv;
  arma::colvec d = Para.d;
  
  //Prepare data
  arma::uvec indeces_row = arma::find(group == id);
  arma::uvec indeces_col = arma::find(group2 == id);
  arma::mat x_i = X.rows(indeces_row);
  arma::colvec x_i_beta = x_i * Beta;
  arma::mat z_i = Z.submat(indeces_row, indeces_col);
  arma::colvec y_i = Y(indeces_row);
  
  //Compute GLMM mean parameter
  arma::mat theta_i = arma::repmat(x_i_beta, 1, R) + z_i * Gamma_i;
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
  arma::colvec grad_beta_likelihood(p, arma::fill::zeros);
  arma::colvec grad_L_likelihood(n_L, arma::fill::zeros); 
  arma::colvec grad_D_likelihood(q, arma::fill::zeros); 
  arma::colvec grad_re_d(q, arma::fill::zeros); 
    
  //Compute gradients by looping over R
  for (arma::uword r = 1; r < R; r++) {
    arma::colvec gamma_ir = Gamma_i.col(r);
    grad_beta_likelihood += arma::trans(arma::trans(y_i - pi_i.col(r)) * x_i);
    arma::colvec v = DInv * gamma_ir;
    arma::colvec w = LInv * v;
    grad_L_likelihood += (get_grad_1_L(L, q) - arma::trans(w) * get_grad_w_L(L, v, w, n_L, q)) * gradLl;
    arma::mat Gamma_r = arma::diagmat(gamma_ir);
    arma::mat M = Gamma_r * UpsilonInv * Gamma_r;
    for (arma::uword k = 0; k < q; k++) {
      double sum1 = 0;
      for (arma::uword h = 0; h < q; h++) sum1 += M(h, k) * exp(-d(h));
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
arma::mat SampleGamma(int id, datobj DatObj, hypara HyPara, sgdobj SgdObj, para Para) {

  //Set data objects
  int q = DatObj.q;
  arma::colvec Y = DatObj.Y;
  arma::mat X = DatObj.X;
  arma::mat Z = DatObj.Z;
  arma::Col<int> group = DatObj.group;
  arma::Col<int> group2 = DatObj.group2;
  
  //Set hyperparameters

  //Set sgd objects
  int R = SgdObj.R;
  
  //Set parameters
  arma::colvec Beta = Para.Beta;
  arma::mat SigmaInv = Para.SigmaInv;
  
  //Prepare data
  arma::uvec indeces_row = arma::find(group == id);
  arma::uvec indeces_col = arma::find(group2 == id);
  arma::mat x_i = X.rows(indeces_row);
  arma::colvec x_i_beta = x_i * Beta;
  arma::mat z_i = Z.submat(indeces_row, indeces_col);
  arma::colvec y_i = Y(indeces_row);
  int n_i = x_i.n_rows;
  
  //Initialize output object
  arma::mat Gamma(q, R, arma::fill::zeros);
  
  //Loop over R
  for (arma::uword r = 1; r < R; r++) {
    
    // Sample omega
    arma::colvec omega_i = pgRcpp(arma::ones(n_i), x_i_beta + z_i * Gamma.col(r - 1));
    
    //Compute moments
    arma::mat D_omega_i = arma::diagmat(omega_i);
    arma::mat tz_i_D_omega_i = arma::trans(z_i) * D_omega_i;
    arma::mat cov_gamma = CholInv(tz_i_D_omega_i * z_i + SigmaInv);
    arma::colvec ystar = ((y_i - 0.5) / omega_i);
    arma::colvec mean_gamma = cov_gamma * (tz_i_D_omega_i * (ystar - x_i_beta));
    
    //Sample gamma
    Gamma.col(r) = rmvnormRcpp(1, mean_gamma, cov_gamma);
  
  //End loop over R  
  }
  
  //Return gamma
  return Gamma;
  
}



//Function to get z matrix-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetZ(arma::vec const& l, int q) {
  arma::mat z(q, q, arma::fill::zeros);
  arma::uvec upper_indices = arma::trimatu_ind(arma::size(z), 1);
  z.elem(upper_indices) += arma::tanh(l);
  return z.t();
}



//Function to get L matrix-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetL(arma::mat const& z, int q) {
  arma::mat L(q, q, arma::fill::zeros);
  for (arma::uword i = 0; i < q; i++) {
    for (arma::uword j = 0; j < q; j++) {
      if (j == 0) {
        if (i == 0) L(i, j) = 1;
        if (i > j) L(i, j) = z(i, j);
      }
      if (j > 0) {
        arma::rowvec Li = L.row(i);
        if (i > j) L(i, j) = z(i, j) * sqrt(1 - arma::as_scalar(arma::sum(arma::pow(Li.subvec(0, j - 1), 2))));         
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

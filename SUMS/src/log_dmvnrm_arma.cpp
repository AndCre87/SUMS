// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "log_dmvnrm_arma.h"

// double log_dmvnrm_arma(arma::rowvec x, arma::rowvec mu, arma::mat Sigma){
//   int p = Sigma.n_cols;
//   double out, rootsum, constants, other_terms;
//   arma::mat root = arma::inv(arma::chol(Sigma));
//   rootsum = arma::sum(log(root.diag()));
//   constants = - (double)p/2.0 * M_LN_2PI;
//   other_terms = rootsum + constants;
// 
//   arma::rowvec z;
//   z = (x - mu) * root;
//   out = other_terms - 0.5 * arma::dot(z,z);
// 
//   return out;
// }


//This function takes the precision!
// [[Rcpp::export]]
double log_dmvnrm_arma(arma::rowvec x, arma::rowvec mu, arma::mat Omega){
  double out = 0.0, sign = 1.0, log_det_Omega = 0.0, p = Omega.n_rows;
  arma::log_det(log_det_Omega, sign, Omega);

  out = - 0.5 * ( p * M_LN_2PI - log_det_Omega + as_scalar((x - mu) * Omega * (x.t() - mu.t())) );
  
  return(out);
}


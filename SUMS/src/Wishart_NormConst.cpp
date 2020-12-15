// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "Wishart_NormConst.h"

double Wishart_NormConst(double nu, arma::mat Psi){
  
  int p = Psi.n_rows, j;
  double log_I, log_det_Psi, sign = 1;
  
  arma::log_det(log_det_Psi, sign, Psi);
  
  log_I = (p * (nu + p - 1) / 2) * log(2) + p * (p - 1) / 4 * log(M_PI);
  log_I += - ((nu + p - 1)/2) * log_det_Psi;
  for(j = 0; j < p; j++){
    log_I += lgamma( (nu + p - 1 - j) / 2); 
  }
  
  return log_I; 
}



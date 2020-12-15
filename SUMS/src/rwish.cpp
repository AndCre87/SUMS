// This function is a wrapper for arma::wishrnd( Psi, nu )

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "rwish.h"

// [[Rcpp::export]]
arma::mat rwish( double nu, arma::mat Psi ){
  
  arma::mat Omega = arma::wishrnd( Psi, nu );
  
  return Omega;
}


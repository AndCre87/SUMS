// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "log_dgamma_arma.h"

double log_dgamma_arma(double x, double a, double b){ 

  double out = 0.0;
  out = a * log(b) - lgamma(a) + (a - 1) * log(x) - b * x;

  return(out);
}
#ifndef LOG_DMVNRM_ARMA_H
#define LOG_DMVNRM_ARMA_H

#include <RcppArmadillo.h>

double log_dmvnrm_arma(arma::rowvec x, arma::rowvec mu, arma::mat Omega);

#endif
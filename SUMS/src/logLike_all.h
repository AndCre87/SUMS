#ifndef LOGLIKE_ALL_H
#define LOGLIKE_ALL_H

#include <RcppArmadillo.h>

arma::vec logLike_all(arma::uvec setC, std::vector<std::vector<arma::vec>> Y, std::vector<std::vector<arma::vec>> H, arma::mat phi_star, arma::mat n_times, std::vector<std::vector<arma::vec>> epsilon, std::vector<arma::mat> X, std::vector<std::vector<arma::mat>> Z, arma::field<arma::mat>  beta, arma::field<arma::mat> gamma, arma::vec n_rates_cum, arma::vec dY, arma::vec dH, bool impute_missing_Y, bool impute_missing_H);
  
#endif
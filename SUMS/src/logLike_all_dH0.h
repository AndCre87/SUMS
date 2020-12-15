#ifndef LOGLIKE_ALL_DH0_H
#define LOGLIKE_ALL_DH0_H

#include <RcppArmadillo.h>

arma::vec logLike_all_dH0(arma::uvec setC, std::vector<std::vector<arma::vec>> Y, arma::mat phi_star, arma::mat n_times, std::vector<std::vector<arma::vec>> epsilon, std::vector<arma::mat> X, std::vector<std::vector<arma::mat>> Z, arma::field<arma::mat> beta, arma::field<arma::mat> gamma, arma::vec n_rates_cum, arma::vec dY, bool impute_missing);

#endif
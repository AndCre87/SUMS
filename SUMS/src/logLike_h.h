#ifndef LOGLIKE_H_H
#define LOGLIKE_H_H

#include <RcppArmadillo.h>

double logLike_h(arma::uvec setC, std::vector<std::vector<arma::vec>> Y, arma::vec s, arma::mat phi_star, arma::vec n_times, std::vector<std::vector<arma::vec>> epsilon, arma::mat X, std::vector<arma::mat> Z, arma::mat beta, arma::mat gamma, arma::vec n_rates_cum, arma::vec dY, int ind_h, bool impute_missing_Y);

#endif
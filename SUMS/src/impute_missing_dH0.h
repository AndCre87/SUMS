#ifndef IMPUTE_MISSING_DH0_H
#define IMPUTE_MISSING_DH0_H

#include <RcppArmadillo.h>

std::vector<std::vector<arma::vec>> impute_missing_dH0(std::vector<std::vector<arma::vec>> Y, arma::mat is_NA_Y, arma::vec c, arma::mat phi_star, arma::mat n_times, std::vector<std::vector<arma::vec>> epsilon, std::vector<arma::mat> X, std::vector<std::vector<arma::mat>> Z, arma::field<arma::mat> beta, arma::field<arma::mat> gamma, arma::vec n_rates_cum, arma::vec dY, bool impute_missing);

#endif
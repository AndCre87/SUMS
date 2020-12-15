#ifndef IMPUTE_MISSING_H
#define IMPUTE_MISSING_H

#include <RcppArmadillo.h>

typedef std::tuple< std::vector<std::vector<arma::vec>>, std::vector<std::vector<arma::vec>> > impute_missing_t;

impute_missing_t impute_missing(std::vector<std::vector<arma::vec>> Y, arma::mat is_NA_Y, std::vector<std::vector<arma::vec>> H, arma::mat is_NA_H, arma::vec c, arma::mat phi_star, arma::mat n_times, std::vector<std::vector<arma::vec>> epsilon, std::vector<arma::mat> X, std::vector<std::vector<arma::mat>> Z, arma::field<arma::mat> beta, arma::field<arma::mat> gamma, arma::vec n_rates_cum, arma::vec dY, arma::vec dH, bool impute_missing_Y, bool impute_missing_H);
  
#endif
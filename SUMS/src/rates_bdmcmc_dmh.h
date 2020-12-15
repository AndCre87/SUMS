#ifndef RATES_BDMCMC_DMH_H
#define RATES_BDMCMC_DMH_H

#include <RcppArmadillo.h>

arma::vec rates_bdmcmc_dmh( double log_ratio_d, bool size_based, double a_d, double b_d, double size_G0, arma::mat G0, arma::mat Psi, arma::mat Psi_star,
                                     arma::mat Omega, arma::mat Omega_tilde, double nu, double nu_star, arma::vec n_states_cum );

#endif
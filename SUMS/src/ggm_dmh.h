#ifndef GGM_DMH_H
#define GGM_DMH_H

#include <RcppArmadillo.h>

typedef std::tuple<arma::mat, arma::mat, arma::mat, double, double> ggm_DMH_t;
  
  ggm_DMH_t ggm_DMH( arma::mat G_mat, arma::mat G0_mat, double d, bool size_based_prior, double a_d, double b_d, double size_G0,
              arma::mat Ts, arma::mat Ti, arma::mat Omega, arma::vec n_states_cum, double threshold,
              double nu, double nu_star, arma::mat Psi, arma::mat Psi_star, int n_edges);

#endif
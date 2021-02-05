// This function and the ones called are modifications of the originals presented in (Mohammadi et al. 2015)

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <tuple>
#include "rgwish_c.h"
#include "rates_bdmcmc_dmh.h"
#include "ggm_dmh.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// birth-death MCMC for Gaussian Graphical models  
// Based on Double Metropolis-Hastings
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

ggm_DMH_t ggm_DMH( arma::mat G, arma::mat G0, double eta, bool size_based_prior, double a_eta, double b_eta, double size_G0,
              arma::mat Ts, arma::mat Ti, arma::mat Omega, arma::vec n_states_cum, double threshold,
              double nu, double nu_star, arma::mat Psi, arma::mat Psi_star, int n_edges){
  
  int p0 = G0.n_cols, p_tot = G.n_cols;
  int qp0 = p0 * (p0-1.0)/2.0;
  
  arma::vec log_rates(qp0,arma::fill::zeros);
  arma::mat Omega_tilde(p_tot,p_tot);
  
  // Log-prior ratio as for death rate
  //Can be extended to edge/nodes-dependent probabilities!
  double log_ratio_eta = 0.0;
  if(!size_based_prior){//Fill with log ratio of prior part for the graph, when proposing a deletion (numerator)
    log_ratio_eta = log(1.0 - eta) - log(eta);
  }

  for(int i_n = 0; i_n < n_edges; i_n ++){
    // - - - STEP 1: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - -|		
    // sampling from Omega and sigma for double Metropolis-Hastings
    arma::mat Omega_tilde = rgwish_c( nu, Ti, G, threshold);		

    log_rates = rates_bdmcmc_dmh( log_ratio_eta, size_based_prior, a_eta, b_eta, size_G0, G0, Psi, Psi_star, Omega, Omega_tilde, nu, nu_star, n_states_cum );

    // Selecting an edge based on birth and death rates
    arma::vec rates = exp(log_rates-max(log_rates));
    rates = rates/sum(rates);
    
    double aux_runif = arma::randu();
    arma::vec rates_cum = arma::cumsum(rates);
    arma::uvec hh_vec = arma::find(rates_cum >= aux_runif, 1, "first");
    int edge_index = hh_vec(0);
    
    int selected_edge_h = edge_index;
    int counter = 0;
    
    while(selected_edge_h >= 0){
      counter ++;
      selected_edge_h -= counter;
    }
    selected_edge_h += counter;
    int selected_edge_k = counter;

    // Updating G0 (graph) based on selected edge
    G0(selected_edge_h, selected_edge_k) = 1.0 - G0(selected_edge_h, selected_edge_k);
    G0(selected_edge_k, selected_edge_h) = G0(selected_edge_h, selected_edge_k);

    //Compute graph size (no need for loop)
    size_G0 += 2.0 * G0(selected_edge_h, selected_edge_k) - 1.0;

    // Applying the update to graph G
    arma::uvec ind_G_h = arma::regspace<arma::uvec>(n_states_cum(selected_edge_h), n_states_cum(selected_edge_h+1)-1);
    arma::uvec ind_G_k = arma::regspace<arma::uvec>(n_states_cum(selected_edge_k), n_states_cum(selected_edge_k+1)-1);

    G(ind_G_h, ind_G_k).fill(G0(selected_edge_h, selected_edge_k));
    G(ind_G_k, ind_G_h).fill(G0(selected_edge_h, selected_edge_k));

    // - - -- STEP 2: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
    Omega = rgwish_c( nu_star, Ts, G, threshold);
  }
  
  double sum_weights = sum(exp(log_rates));
  
  return ggm_DMH_t(Omega, G, G0, size_G0, sum_weights);
}

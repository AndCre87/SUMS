// This function and the ones called are modifications of the originals presented in (Mohammadi et al. 2015)

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <omp.h>
#include "Wishart_NormConst.h"
#include "rates_bdmcmc_dmh.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Parallel Computation for birth-death rates for double BD-MCMC algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

arma::vec rates_bdmcmc_dmh( double log_ratio_eta, bool size_based_prior, double a_eta, double b_eta, double size_G0, arma::mat G0, arma::mat Psi, arma::mat Psi_star,
                                arma::mat Omega, arma::mat Omega_tilde, double nu, double nu_star, arma::vec n_states_cum ){
  
  int p0 = G0.n_cols, p_tot = Omega.n_cols;
  
  double qp0 = p0*(p0-1.0)/2.0, sign = 1.0, log_det_Omega_h= 0.0, log_det_Omega_tilde_h = 0.0;
  arma::vec log_rates(qp0);
  
  arma::uvec ind_G_aux = arma::regspace<arma::uvec>(0, p_tot-1);
  
  // #pragma omp parallel for
  for(int h = 0; h < (p0-1); h++ ){			
    
    arma::uvec ind_G_h = arma::regspace<arma::uvec>(n_states_cum(h), n_states_cum(h+1)-1);
    int ph = ind_G_h.n_elem;
    
    arma::mat Psi_hh = Psi.submat(ind_G_h, ind_G_h);
    arma::mat Psi_star_hh = Psi_star.submat(ind_G_h, ind_G_h);
    
    arma::mat Omega_hh = Omega.submat(ind_G_h, ind_G_h);
    arma::mat Omega_tilde_hh = Omega_tilde.submat(ind_G_h, ind_G_h);
    
    for(int k = h+1; k < p0; k++ ){
      
      arma::uvec ind_G_k = arma::regspace<arma::uvec>(n_states_cum(k), n_states_cum(k+1)-1);
      int pk = ind_G_k.n_elem;
      
      arma::mat Psi_kk = Psi.submat(ind_G_k, ind_G_k);
      arma::mat Psi_star_kk = Psi_star.submat(ind_G_k, ind_G_k);
      
      arma::uvec ind_G_hk = join_cols(ind_G_h,ind_G_k);
      
      arma::mat Psi_hkhk = Psi.submat(ind_G_hk, ind_G_hk);
      arma::mat Psi_star_hkhk = Psi_star.submat(ind_G_hk, ind_G_hk);
      
      arma::uvec ind_G_nok = ind_G_aux;
      ind_G_nok.shed_rows(ind_G_k(0), ind_G_k(pk-1));
      arma::uvec ind_G_noh = ind_G_aux;
      ind_G_noh.shed_rows(ind_G_h(0), ind_G_h(ph-1));
      arma::uvec ind_G_nohk = intersect(ind_G_noh,ind_G_nok);
      
      arma::mat Psi_hnohk = Psi.submat(ind_G_h, ind_G_nohk);
      arma::mat Psi_nohknohk = Psi.submat(ind_G_nohk, ind_G_nohk);
      arma::mat Psi_star_hnohk = Psi_star.submat(ind_G_h, ind_G_nohk);
      arma::mat Psi_star_nohknohk = Psi_star.submat(ind_G_nohk, ind_G_nohk);  
      arma::mat Psi_hhc = Psi_hh - Psi_hnohk * arma::inv_sympd(Psi_nohknohk) * Psi_hnohk.t();
      arma::mat Psi_star_hhc = Psi_star_hh - Psi_star_hnohk * arma::inv_sympd(Psi_star_nohknohk) * Psi_star_hnohk.t();
      
      
      //Compute log-normalizing constants from Wishart distributions
      double log_I_star_kk = Wishart_NormConst(nu_star, Psi_star_kk);
      double log_I_star_hkhk = Wishart_NormConst(nu_star, Psi_star_hkhk);
      double log_I_star_hhc = Wishart_NormConst(nu_star + ph, Psi_star_hhc);
      
      double log_I_kk = Wishart_NormConst(nu, Psi_kk);
      double log_I_hkhk = Wishart_NormConst(nu, Psi_hkhk);
      double log_I_hhc = Wishart_NormConst(nu + ph, Psi_hhc);
      
      //Compute Omega_a
      arma::mat Omega_knok = Omega.submat(ind_G_k, ind_G_nok);
      arma::mat Omega_noknok = Omega.submat(ind_G_nok, ind_G_nok);
      arma::mat Omega_a = Omega;
      Omega_a.submat(ind_G_k, ind_G_k) = Omega_knok * arma::inv_sympd(Omega_noknok) * Omega_knok.t();
      if(G0(h,k) == 1){
        arma::mat zeros_aux1(ph,pk,arma::fill::zeros);
        Omega_a.submat(ind_G_h, ind_G_k) = zeros_aux1;
        Omega_a.submat(ind_G_k, ind_G_h) = zeros_aux1.t();
      }
      
      arma::mat Omega_tilde_knok = Omega_tilde.submat(ind_G_k, ind_G_nok);
      arma::mat Omega_tilde_noknok = Omega_tilde.submat(ind_G_nok, ind_G_nok);
      arma::mat Omega_tilde_a = Omega_tilde;
      Omega_tilde_a.submat(ind_G_k, ind_G_k) = Omega_tilde_knok * arma::inv_sympd(Omega_tilde_noknok) * Omega_tilde_knok.t();
      if(G0(h,k) == 1){
        arma::mat zeros_aux2(ph,pk,arma::fill::zeros);
        Omega_tilde_a.submat(ind_G_h, ind_G_k) = zeros_aux2;
        Omega_tilde_a.submat(ind_G_k, ind_G_h) = zeros_aux2.t();
      }
      
      //Compute Omega_b
      arma::mat Omega_hknohk = Omega.submat(ind_G_hk, ind_G_nohk);
      arma::mat Omega_nohknohk = Omega.submat(ind_G_nohk, ind_G_nohk);
      arma::mat Omega_b = Omega;
      Omega_b.submat(ind_G_hk, ind_G_hk) = Omega_hknohk * arma::inv_sympd(Omega_nohknohk) * Omega_hknohk.t();
      arma::mat Omega_b_hnohc = Omega_b.submat(ind_G_h, ind_G_h);
      
      arma::mat Omega_tilde_hknohk = Omega_tilde.submat(ind_G_hk, ind_G_nohk);
      arma::mat Omega_tilde_nohknohk = Omega_tilde.submat(ind_G_nohk, ind_G_nohk);
      arma::mat Omega_tilde_b = Omega_tilde;
      Omega_tilde_b.submat(ind_G_hk, ind_G_hk) = Omega_tilde_hknohk * arma::inv_sympd(Omega_tilde_nohknohk) * Omega_tilde_hknohk.t();      
      arma::mat Omega_tilde_b_hnohc = Omega_tilde_b.submat(ind_G_h, ind_G_h);
      
      //Compute log-H functions
      
      arma::log_det(log_det_Omega_h, sign, Omega_hh - Omega_b_hnohc);
      double log_H = - 0.5 * ( arma::trace(Psi_star * (Omega_a - Omega_b)) - arma::trace(Psi_star_hhc * (Omega_hh - Omega_b_hnohc)) + log_det_Omega_h);
      log_H += log_I_star_kk + log_I_star_hhc - log_I_star_hkhk;
      
      arma::log_det(log_det_Omega_tilde_h, sign, Omega_tilde_hh - Omega_tilde_b_hnohc);
      double log_H_tilde = - 0.5 * (arma::trace(Psi * (Omega_tilde_a - Omega_tilde_b)) - arma::trace(Psi_hhc * (Omega_tilde_hh - Omega_tilde_b_hnohc)) + log_det_Omega_tilde_h);
      log_H_tilde += log_I_kk + log_I_hhc - log_I_hkhk;
      
      double log_rate = log_H - log_H_tilde;

      if(G0(h,k) == 0){// It is a birth event
        log_rate = - log_rate;
        //Prior on graph
        if(size_based_prior){        
          log_rate += log(a_eta + size_G0) - log(b_eta + qp0 - size_G0 - 1.0);
        }else{
          log_rate += - log_ratio_eta;
        }
      }else{// It is a death event
        //Prior on graph
        if(size_based_prior){
          log_rate += log(b_eta + qp0 - size_G0) - log(a_eta + size_G0 - 1.0);
        }else{
          log_rate += log_ratio_eta;
        }
      }

      //Fill upper-tri part!
      log_rates(k*(k-1)/2 + h) = log_rate;
    }
  }
  return log_rates;
}


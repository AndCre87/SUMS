//Gibbs sampler for estimating Covariance and Precision matrix in Gaussian graphical model

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <tuple>
#include "ggm_dmh.h"
#include "rgwish_c.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
List DisProgrCov_Gibbs(List MCMC_input){
  
  // Data (here the phi's of the model)
  arma::mat phi = MCMC_input["data"];
  arma::vec n_rates = MCMC_input["n_rates"];
  arma::vec n_rates_cum = arma::cumsum(n_rates);

  // Extract sizes
  int N = phi.n_rows;
  int p_tot = phi.n_cols;
  int p0 = n_rates.n_elem - 1;
  
  // Extract hyperparameters
  List Param_list = MCMC_input["Param_list"];
  double nu = Param_list["nu"];
  // double nu_star = nu + N + 1; //Update conditionally to mu
  double nu_star = nu + N; //Marginal update
  arma::mat Psi = Param_list["Psi"];
  
  // Initialize graphs G0 and G as empty
  arma::mat G0(p0,p0,arma::fill::zeros), G(p_tot,p_tot,arma::fill::zeros);
  // But recall that some edges are forced in graph G!
  for(int j = 0; j < p0; j++){
    arma::uvec ind_G_j = arma::regspace<arma::uvec>(n_rates_cum(j), n_rates_cum(j+1)-1);
    for(int j1 = 0; j1 < (n_rates(j+1) - 1); j1++){
      for(int j2 = j1 + 1; j2 < n_rates(j+1); j2++){
        G(ind_G_j(j1),ind_G_j(j2)) = 1.0;
        G(ind_G_j(j2),ind_G_j(j1)) = G(ind_G_j(j1),ind_G_j(j2));
      }
    }
  }
  
  // Hyperpriors
  arma::rowvec mu = Param_list["mu"], m_mu = Param_list["m_mu"];
  arma::mat Sigma_m_mu(p_tot,p_tot,arma::fill::eye), Prec_m_mu = arma::inv_sympd(Sigma_m_mu);
  
  bool update_mu = Param_list["update_mu"];
  
  bool update_m_mu = Param_list["update_m_mu"];
  if(update_m_mu){
    Sigma_m_mu = as<arma::mat>(Param_list["Sigma_m_mu"]);
  }
  
  double a_k0 = 0.0, b_k0 = 0.0, k0 = Param_list["k0"];
  bool update_k0 = Param_list["update_k0"];
  if(update_k0){
    a_k0 = Param_list["a_k0"];
    b_k0 = Param_list["b_k0"];
  }
  
  //We have different prior options for Graph
  double n_edges0 = p0 * (p0 - 1) / 2, size_G0 = 0.0;
  double a_d = 0.0, b_d = 0.0, d = 0.0, a_lambda = 0.0, d_accept = 0.0, d_count = 0.0, s_d = 0.0, lambda = 0.0, lambda_accept = 0.0, lambda_count = 0.0, s_lambda = 0.0, mu_d = 0.0;
  bool size_based_prior = Param_list["size_based_prior"], update_d = false, d_beta_prior = false;
  if(size_based_prior){
    a_d = Param_list["a_d"], b_d = Param_list["b_d"];
    d = 0.5; // Initialize
  }else{
    update_d = Param_list["update_d"];
    if(update_d){
      d_beta_prior = Param_list["d_beta_prior"];
      if(d_beta_prior || size_based_prior){
        a_d = Param_list["a_d"], b_d = Param_list["b_d"];
        d = 0.5; // Initialize
      }else{
        a_lambda = Param_list["a_lambda"];
        lambda = 1.0; // Initialize
        d = 0.5;
        s_d = 0.01;
        s_lambda = 0.01;
      }
    }else{
      d = Param_list["d"];
    }
  }
  
  // Initialize matrices and vectors
  arma::mat Omega(p_tot,p_tot,arma::fill::eye), Ti = arma::chol( arma::inv_sympd( Psi ) );
  double sum_weights;
  
  // Iterations after first burn-in
  List Alg_list = MCMC_input["Alg_list"];
  double n_burn1 = Alg_list["n_burn1"], n_burn2 = Alg_list["n_burn2"], thin = Alg_list["thin"], n_save = Alg_list["n_save"], threshold = Alg_list["threshold"];
  int n_tot = n_burn1 + n_burn2 + thin*n_save, n_edges = Alg_list["n_edges"];
  
  // Adaptation
  NumericVector ADAPT(4);
  ADAPT(0) = n_burn1; //"burn-in" for adaptation
  ADAPT(1) = 0.7; //exponent for adaptive step
  ADAPT(2) = 0.234; //reference acceptance rate
  ADAPT(3) = 0.001; //for multivariate updates
  
  
  // Output lists
  List ggm_result;
  arma::cube Omega_out(p_tot, p_tot, n_save, arma::fill::zeros), G_out(p_tot, p_tot, n_save, arma::fill::zeros), G0_out(p0, p0, n_save, arma::fill::zeros);
  
  arma::mat mu_out(n_save, p_tot, arma::fill::zeros), m_mu_out(n_save, p_tot, arma::fill::zeros);
  
  arma::vec d_out(n_save, arma::fill::zeros), size_G0_out(n_save, arma::fill::zeros), sum_weights_out(n_save,arma::fill::zeros), k0_out(n_save, arma::fill::zeros), lambda_out(n_save,arma::fill::zeros), mu_d_out(n_save,arma::fill::zeros);
  
  Progress progr(n_tot, true);
  
  for(int g = 0; g < n_tot; g++){
    
    ///////////////////////////
    // Update G0, G, Omega_G //
    ///////////////////////////

    //Marginal update
    arma::rowvec phi_mean = sum(phi) / N;
    arma::mat Psi_star = Psi + k0 * N / (k0 + N) * (m_mu.t() - phi_mean.t()) * (m_mu - phi_mean);
    for(int i = 0; i < N; i++){
      Psi_star = Psi_star + (phi.row(i).t()  - phi_mean.t()) * (phi.row(i)  - phi_mean);
    }
    
    arma::mat Ts = arma::chol( arma::inv_sympd( Psi_star ) );
    
    std::tie(Omega, G, G0, size_G0, sum_weights) = ggm_DMH(G, G0, d, size_based_prior, a_d, b_d, size_G0, Ts, Ti, Omega, n_rates_cum, threshold, nu, nu_star, Psi, Psi_star, n_edges);

    

    
    
    
    
    ////////////////////////////
    // Update Omega_G (extra) //
    ////////////////////////////
    
    Omega = rgwish_c( nu_star, Ts, G, threshold );
    
    
    
    
    
    ///////////////
    // Update mu //
    ///////////////
    
    // Conditionally independent prior on mu
    if(update_mu){
      arma::rowvec sum_of_phi = sum(phi);
      arma::mat var_mu = arma::inv_sympd((N + k0) * Omega);
      arma::vec mean_mu = (k0 * m_mu.t() + sum_of_phi.t())/(N + k0);
      mu = arma::mvnrnd(mean_mu, var_mu).t();
    }
    
    
    
    
    /////////////////
    // Update m_mu //
    /////////////////
    
    //Additional prior on mean of mu
    if(update_m_mu){
      arma::mat var_m_mu = arma::inv_sympd(Prec_m_mu + k0 * Omega);
      arma::vec mean_m_mu = var_m_mu * k0 * Omega * mu.t();
      m_mu = arma::mvnrnd(mean_m_mu, var_m_mu).t();
    }
    
    
    
    
    ///////////////
    // Update k0 //
    ///////////////
    
    //Additional prior on k0
    if(update_k0){
      double k0_aux = b_k0 + .5 * as_scalar((mu - m_mu) * Omega * (mu.t() - m_mu.t()));
      k0 = arma::randg( arma::distr_param( a_k0 + p_tot/2.0, 1.0/k0_aux ));
    }
    
    
    
    
    //////////////
    // Update d //
    //////////////
    if(!size_based_prior){
      if(update_d){
        if(d_beta_prior){
          d = R::rbeta(a_d + size_G0, b_d + n_edges0 - size_G0);
        }else{
          //Update d \sim logit-Normal(0,lambda^2)
          double d_new = log(d/(1.0 - d)) + arma::randn() * sqrt(s_d);
          d_new = 1.0/(1.0 + exp(-d_new));
          
          double log_ratio_d = size_G0 * (log(d_new) - log(d)) + (n_edges0 - size_G0) * (log(1 - d_new) - log(1 - d)) - .5/pow(lambda,2) * (pow(log(d_new) - log(1 - d_new) - mu_d,2) - pow(log(d) - log(1 - d) - mu_d,2));
          
          double accept_d = 1.0;
          if( arma::is_finite(log_ratio_d) ){
            if(log_ratio_d < 0){
              accept_d = exp(log_ratio_d);
            }
          }else{
            accept_d = 0.0;
          }
          
          d_accept += accept_d;
          d_count ++;
          
          if( arma::randu() < accept_d ){
            d = d_new;
          }
          
          s_d = s_d + pow(g+1,-ADAPT(1))*(accept_d - ADAPT(2));
          if(s_d > exp(50)){
            s_d = exp(50);
          }else{
            if(s_d < exp(-50)){
              s_d = exp(-50);
            }
          }
          
          
          
          
          //Update mu_d (univariate conjugate normal)
          mu_d = (log(d) - log(1.0 - d))/2 + arma::randn() * sqrt(pow(lambda,2.0)/2.0);
          
          
          
          
          //Update lambda \sim Cauchy(a)
          double lambda_new = lambda * exp(arma::randn() * sqrt(s_lambda));
          
          double log_ratio_lambda = log(a_lambda + pow(lambda,2.0)) - log(a_lambda + pow(lambda_new,2.0)) - 0.5 * pow(log(d) - log(1.0 - d) - mu_d,2.0) * (1.0 / pow(lambda_new,2.0) - 1.0 / pow(lambda,2.0));
          
          double accept_lambda = 1.0;
          if( arma::is_finite(log_ratio_lambda) ){
            if(log_ratio_lambda < 0){
              accept_lambda = exp(log_ratio_lambda);
            }
          }else{
            accept_lambda = 0.0;
          }
          
          lambda_accept += accept_lambda;
          lambda_count ++;
          
          if( arma::randu() < accept_lambda ){
            lambda = lambda_new;
          }
          
          s_lambda = s_lambda + pow(g+1,-ADAPT(1))*(accept_lambda - ADAPT(2));
          if(s_lambda > exp(50)){
            s_lambda = exp(50);
          }else{
            if(s_lambda < exp(-50)){
              s_lambda = exp(-50);
            }
          }
          
        }
      }
    }
    
    
    
    
    
    //Save output for this iteration
    if( (g + 1 > (n_burn1 + n_burn2)) & (((g + 1 - n_burn1 - n_burn2) / thin - floor((g + 1 - n_burn1 - n_burn2) / thin)) == 0 )){
      
      int iter = (g + 1 - n_burn1 - n_burn2)/thin - 1;
      
      Omega_out.slice(iter) = Omega;
      G_out.slice(iter) = G;
      G0_out.slice(iter) = G0;
      
      mu_out.row(iter) = mu;
      m_mu_out.row(iter) = m_mu;
      
      d_out(iter) = d;
      lambda_out(iter) = lambda;
      mu_d_out(iter) = mu_d;
      size_G0_out(iter) = size_G0;
      sum_weights_out(iter) = sum_weights;
      k0_out(iter) = k0;
    }
    
    //Progress bar increment
    progr.increment(); 
  }
  
  if(update_d){
    if(!d_beta_prior){
      Rcout << "d a.r. = " << d_accept/d_count << "\n";
      Rcout << "lambda a.r. = " << lambda_accept/lambda_count << "\n";
    }
  }
  
  return List::create(Named("Omega_out") = Omega_out, Named("G_out") = G_out, Named("G0_out") = G0_out, Named("sizeG0_out") = size_G0_out, Named("sum_weights_out") = sum_weights_out, Named("mu_out") = mu_out, Named("m_mu_out") = m_mu_out, Named("d_out") = d_out, Named("k0_out") = k0_out, Named("mu_d_out") = mu_d_out, Named("lambda_out") = lambda_out);
}


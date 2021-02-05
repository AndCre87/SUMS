// Main Gibbs sampler function for Disease Progression model (parametric)

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <tuple>
// #include <omp.h>
#include "impute_missing.h"
#include "ggm_dmh.h"
#include "logLike_all.h"
#include "logLike_h.h"
#include "rgwish_c.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
List DisProgrGibbs_Param(List MCMC_input){
  
  // Data
  std::vector<std::vector<arma::vec>> Y = MCMC_input["Y_obs"], H = MCMC_input["H_obs"], epsilon = MCMC_input["epsilon"];
  
  //Missing values at time0 that require imputation
  arma::mat is_NA_Y = as<arma::mat>(MCMC_input["is_NA_Y"]);
  arma::mat is_NA_H = as<arma::mat>(MCMC_input["is_NA_H"]);
  //Introduce a bool variable saying if we are conditioning to the initial value or not
  //This is linked to the imputation of missing values for us
  //But we could still define a model for t0 with values observed at all times/processes
  bool impute_missing_Y = false, impute_missing_H = false;
  if(arma::accu(is_NA_Y) > 0){
    impute_missing_Y = true;
  }
  if(arma::accu(is_NA_H) > 0){
    impute_missing_H = true;
  }
  
  
  //If they are NULL, they are replaced by zeroes so they won't contribute int he likelihood term
  bool is_null_Z = MCMC_input["is_null_Z"];
  std::vector<std::vector<arma::mat>> Z = as<std::vector<std::vector<arma::mat>>>(MCMC_input["Z"]);
  
  bool is_null_X = MCMC_input["is_null_X"];
  std::vector<arma::mat> X = as<std::vector<arma::mat>>(MCMC_input["X"]);
  
  
  arma::vec n_rates = MCMC_input["n_rates"], n_rates_cum = arma::cumsum(n_rates);
  arma::mat n_times = MCMC_input["n_times"];
  
  // Extract sizes
  int N = Y.size();
  std::vector<arma::vec> Hi = H[0], Yi = Y[0];
  int pY = Yi.size();
  int pH = Hi.size();
  int p0 = pY + pH;
  
  arma::vec g(pY,arma::fill::ones), q(pY,arma::fill::ones);
  for(int h = 0; h< pY; h++){
    
    if(!is_null_X){
      arma::mat X_h = X[h];
      g(h) = X_h.n_cols;
    }else{
      g(h) = 1; //Assumption when X is null
    }
    
    if(!is_null_Z){
      std::vector<arma::mat> Z_h = Z[h];
      arma::mat Z_h_i = Z_h[0];
      q(h) = Z_h_i.n_cols;
    }else{
      q(h) = 1; //Assumption when Z is null
    }
    
  }
  
  //Compute number of rates (dimension of phi_star)
  int p_tot = arma::accu(n_rates); //First is equal to zero
  
  //We assume general number of states
  arma::vec dY = MCMC_input["dY"];
  arma::vec dH = MCMC_input["dH"];
  
  // Extract hyperparameters
  List Param_list = MCMC_input["Param_list"];
  double nu = Param_list["nu"], nu_star = nu + N;
  arma::mat Psi = Param_list["Psi"];
  
  // Initialize graphs G0 and G as empty
  arma::mat G0(p0,p0,arma::fill::zeros), G(p_tot,p_tot,arma::fill::zeros);
  // But recall that some edges are forced in graph G!
  // #pragma omp parallel for
  for(int j = 0; j < p0; j++){
    arma::uvec ind_G_j = arma::regspace<arma::uvec>(n_rates_cum(j), n_rates_cum(j+1)-1);
    for(int j1 = 0; j1 < (n_rates(j+1) - 1); j1++){
      for(int j2 = j1 + 1; j2 < n_rates(j+1); j2++){
        G(ind_G_j(j1),ind_G_j(j2)) = 1;
        G(ind_G_j(j2),ind_G_j(j1)) = G(ind_G_j(j1),ind_G_j(j2));
      }
    }
  }
  
  // Initalize rates
  arma::mat S_phi_star(p_tot,p_tot,arma::fill::eye);
  S_phi_star = 0.01 * S_phi_star;
  
  // Hyperpriors
  arma::rowvec mu = Param_list["mu"], m_mu = Param_list["m_mu"];
  arma::mat Sigma_m_mu(p_tot,p_tot,arma::fill::eye);
  
  bool update_mu = Param_list["update_mu"];
  
  bool update_m_mu = Param_list["update_m_mu"];
  if(update_m_mu){
    Sigma_m_mu = as<arma::mat>(Param_list["Sigma_m_mu"]);
  }
  arma::mat Prec_m_mu = arma::inv_sympd(Sigma_m_mu);
  
  double a_k0 = 0.0, b_k0 = 0.0, k0 = Param_list["k0"];
  bool update_k0 = Param_list["update_k0"];
  if(update_k0){
    a_k0 = Param_list["a_k0"];
    b_k0 = Param_list["b_k0"];
  }
  
  //We have different prior options for Graph
  double n_edges0 = p0 * (p0 - 1) / 2, size_G0 = 0.0;
  double a_eta = 0.0, b_eta = 0.0, eta = 0.0, a_lambda = 0.0, eta_accept = 0.0, eta_count = 0.0, s_eta = 0.0, lambda = 0.0, lambda_accept = 0.0, lambda_count = 0.0, s_lambda = 0.0, mu_eta = 0.0;
  bool size_based_prior = Param_list["size_based_prior"], update_eta = false, eta_beta_prior = false;
  if(size_based_prior){
    a_eta = Param_list["a_eta"], b_eta = Param_list["b_eta"];
    eta = 0.5; // Initialize
  }else{
    update_eta = Param_list["update_eta"];
    if(update_eta){
      eta_beta_prior = Param_list["eta_beta_prior"];
      if(eta_beta_prior || size_based_prior){
        a_eta = Param_list["a_eta"], b_eta = Param_list["b_eta"];
        eta = 0.5; // Initialize
      }else{
        a_lambda = Param_list["a_lambda"];
        lambda = 1.0; // Initialize
        eta = 0.5;
        s_eta = 0.01;
        s_lambda = 0.01;
      }
    }else{
      eta = Param_list["eta"];
    }
  }
  
  // Covariates, beta and gamma coefficients
  arma::vec dim_beta = g % dY % (dY - 1.0);
  arma::vec s_d_beta = dim_beta, beta_accept(pY,arma::fill::zeros), beta_count(pY,arma::fill::zeros);
  arma::field<arma::mat> beta(pY);
  std::vector<arma::mat> mu_beta(pY), mu_beta_vec(pY), U_beta(pY), V_beta(pY), VU_beta(pY), VU_beta_inv(pY), eye_mat_beta(pY);
  std::vector<arma::mat> beta_vec(pY), sum_beta(pY), prod_beta(pY), S_beta(pY);
  
  if(!is_null_X){
    
    dim_beta = g % dY % (dY - 1.0);
    s_d_beta = pow(2.4,2)/dim_beta;
    for(int h = 0; h < pY; h++){
      arma::mat beta_aux(g(h),dY(h)*(dY(h) - 1.0),arma::fill::zeros);
      beta(h) = beta_aux;
    }
    
    mu_beta = as<std::vector<arma::mat>>(Param_list["mu_beta"]);
    U_beta = as<std::vector<arma::mat>>(Param_list["U_beta"]);
    V_beta = as<std::vector<arma::mat>>(Param_list["V_beta"]);
    for(int h = 0; h < pY; h++){
      //beta
      arma::mat mu_beta_aux = mu_beta[h];
      mu_beta_vec[h]  = vectorise(mu_beta_aux);
      
      arma::mat U_beta_aux = U_beta[h];
      arma::mat V_beta_aux = V_beta[h];
      arma::mat VU_beta_aux = arma::kron(V_beta_aux,U_beta_aux), VU_beta_inv_aux = arma::inv_sympd(VU_beta_aux);
      VU_beta[h] = VU_beta_aux;
      VU_beta_inv[h] = VU_beta_inv_aux;
      
      arma::mat eye_mat_beta_aux(dim_beta(h),dim_beta(h),arma::fill::eye);
      eye_mat_beta[h] = eye_mat_beta_aux;
    }
    
    //Adaptation for beta's and gamma's
    for(int h = 0; h < pY; h++){
      //beta
      arma::vec beta_vec_aux(dim_beta(h),arma::fill::zeros);
      beta_vec[h] = beta_vec_aux;
      sum_beta[h] = beta_vec_aux;
      arma::mat beta_prod_aux = beta_vec_aux * beta_vec_aux.t();
      prod_beta[h] = beta_prod_aux;
      
      arma::mat eye_mat_beta_aux(dim_beta(h),dim_beta(h),arma::fill::eye);
      S_beta[h] = 0.01 * eye_mat_beta_aux;
    }
  }else{//Initialise beta to vectors of zeros
    dim_beta = g % dY % (dY - 1.0);
    for(int h = 0; h < pY; h++){
      arma::mat beta_aux(1,dY(h)*(dY(h) - 1.0),arma::fill::zeros);
      beta(h) = beta_aux;
    }
  }
  
  arma::vec dim_gamma = q % dY % (dY - 1.0);
  arma::vec s_d_gamma = dim_gamma, gamma_accept(pY,arma::fill::zeros), gamma_count(pY,arma::fill::zeros);
  arma::field<arma::mat> gamma(pY); //Only these two are arma::field's! For saving them in output (List)
  std::vector<arma::mat> mu_gamma(pY), mu_gamma_vec(pY), U_gamma(pY), V_gamma(pY), VU_gamma(pY), VU_gamma_inv(pY), eye_mat_gamma(pY);
  std::vector<arma::mat> gamma_vec(pY), sum_gamma(pY), prod_gamma(pY), S_gamma(pY);
  
  if(!is_null_Z){
    
    dim_gamma = q % dY % (dY - 1.0);
    s_d_gamma = pow(2.4,2)/dim_gamma;
    for(int h = 0; h < pY; h++){
      arma::mat gamma_aux(q(h),dY(h)*(dY(h) - 1.0),arma::fill::zeros);
      gamma(h) = gamma_aux;
    }
    
    mu_gamma = as<std::vector<arma::mat>>(Param_list["mu_gamma"]);
    U_gamma = as<std::vector<arma::mat>>(Param_list["U_gamma"]);
    V_gamma = as<std::vector<arma::mat>>(Param_list["V_gamma"]);
    for(int h = 0; h < pY; h++){
      //gamma
      arma::mat mu_gamma_aux = mu_gamma[h];
      mu_gamma_vec[h]  = vectorise(mu_gamma_aux);
      
      arma::mat U_gamma_aux = U_gamma[h];
      arma::mat V_gamma_aux = V_gamma[h];
      arma::mat VU_gamma_aux = arma::kron(V_gamma_aux,U_gamma_aux), VU_gamma_inv_aux = arma::inv_sympd(VU_gamma_aux);
      VU_gamma[h] = VU_gamma_aux;
      VU_gamma_inv[h] = VU_gamma_inv_aux;
      
      arma::mat eye_mat_gamma_aux(dim_gamma(h),dim_gamma(h),arma::fill::eye);
      eye_mat_gamma[h] = eye_mat_gamma_aux;
    }
    
    //Adaptation for beta's and gamma's
    for(int h = 0; h < pY; h++){
      //gamma
      arma::vec gamma_vec_aux(dim_gamma(h),arma::fill::zeros);
      gamma_vec[h] = gamma_vec_aux;
      sum_gamma[h] = gamma_vec_aux;
      arma::mat gamma_prod_aux = gamma_vec_aux * gamma_vec_aux.t();
      prod_gamma[h] = gamma_prod_aux;
      
      arma::mat eye_mat_gamma_aux(dim_gamma(h),dim_gamma(h),arma::fill::eye);
      S_gamma[h] = 0.01 * eye_mat_gamma_aux;
    }
  }else{//Initialise gamma to vectors of zeros
    dim_gamma = q % dY % (dY - 1.0);
    for(int h = 0; h < pY; h++){
      arma::mat gamma_aux(1,dY(h)*(dY(h) - 1.0),arma::fill::zeros);
      gamma(h) = gamma_aux;
    }
  }
  
  // Initialize matrices and vectors
  arma::uvec setC_all = arma::regspace<arma::uvec>(0,N-1);
  arma::mat Omega(p_tot,p_tot,arma::fill::eye), Sigma = arma::inv_sympd(Omega), Ti = arma::chol( arma::inv_sympd( Psi ) );
  double sum_weights;
  
  
  //Initialize from P0
  arma::vec c = MCMC_input["c_init"]; //allocation variables
  arma::mat phi_star = mvnrnd(mu.t(), Sigma, N).t();
  
  // Iterations after first burn-in
  List Alg_list = MCMC_input["Alg_list"];
  double n_burn1 = Alg_list["n_burn1"], n_burn2 = Alg_list["n_burn2"], thin = Alg_list["thin"], n_save = Alg_list["n_save"], threshold = Alg_list["threshold"];
  int n_tot = n_burn1 + n_burn2 + thin*n_save, n_edges = Alg_list["n_edges"];
  
  // Adaptation
  arma::vec ADAPT(4);
  ADAPT(0) = n_burn1; //"burn-in" for adaptation
  ADAPT(1) = 0.7; //exponent for adaptive step
  ADAPT(2) = 0.234; //reference acceptance rate
  ADAPT(3) = 0.001; //for multivariate updates
  
  
  // Output lists
  List beta_out(n_save), gamma_out(n_save), Y_out(n_save), H_out(n_save);
  arma::cube phi_star_out(N,p_tot,n_save,arma::fill::zeros), Omega_out(p_tot,p_tot,n_save,arma::fill::zeros), G_out(p_tot,p_tot,n_save,arma::fill::zeros), G0_out(p0,p0,n_save,arma::fill::zeros);
  arma::mat mu_out(n_save, p_tot, arma::fill::zeros), m_mu_out(n_save, p_tot, arma::fill::zeros);
  arma::vec eta_out(n_save, arma::fill::zeros), size_G0_out(n_save, arma::fill::zeros), sum_weights_out(n_save,arma::fill::zeros), k0_out(n_save, arma::fill::zeros), lambda_out(n_save,arma::fill::zeros), mu_eta_out(n_save,arma::fill::zeros);
  
  // Main Gibbs
  Progress progr(n_tot, true);
  
  for(int it = 0; it < n_tot; it++){
    
    ///////////////////////////
    // Impute missing values //
    ///////////////////////////
    std::tie(Y, H) = impute_missing(Y, is_NA_Y, H, is_NA_H, c, phi_star, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H);
    
    
    
    ///////////////////////////
    // Update G0, G, Omega_G //
    ///////////////////////////
    
    //Marginal update
    arma::rowvec phi_star_mean = sum(phi_star) / N;
    arma::mat Psi_star = Psi + k0 * N / (k0 + N) * (m_mu.t() - phi_star_mean.t()) * (m_mu - phi_star_mean);
    for(int i = 0; i < N; i++){
      Psi_star += (phi_star.row(i).t()  - phi_star_mean.t()) * (phi_star.row(i)  - phi_star_mean);
    }
    
    //Cholesky gives upper tri here
    arma::mat Ts = arma::chol( arma::inv_sympd( Psi_star ) );
    
    std::tie(Omega, G, G0, size_G0, sum_weights) = ggm_DMH(G, G0, eta, size_based_prior, a_eta, b_eta, size_G0, Ts, Ti, Omega, n_rates_cum, threshold, nu, nu_star, Psi, Psi_star, n_edges);
    
    
    
    
    
    ////////////////////////////
    // Update Omega_G (extra) //
    ////////////////////////////
    
    Omega = rgwish_c( nu_star, Ts, G, threshold );
    
    
    
    
    
    ///////////////
    // Update mu //
    ///////////////
    
    // Conditionally independent prior on mu
    if(update_mu){
      arma::rowvec sum_of_phi_star = arma::sum(phi_star);
      arma::mat var_mu = arma::inv_sympd((N + k0) * Omega);
      arma::vec mean_mu = (k0 * m_mu.t() + sum_of_phi_star.t())/(N + k0);
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
      k0 = arma::randg( arma::distr_param( a_k0 + p_tot/2, 1.0/k0_aux ));
    }
    
    
    
    
    
    ////////////////
    // Update eta //
    ////////////////
    if(!size_based_prior){
      if(update_eta){
        if(eta_beta_prior){
          eta = R::rbeta(a_eta + size_G0, b_eta + n_edges0 - size_G0);
        }else{
          //Update eta \sim logit-Normal(0,lambda^2)
          double eta_new = log(eta/(1.0 - eta)) + arma::randn() * sqrt(s_eta);
          eta_new = 1.0/(1.0 + exp(-eta_new));
          
          double log_ratio_eta = size_G0 * (log(eta_new) - log(eta)) + (n_edges0 - size_G0) * (log(1 - eta_new) - log(1 - eta)) - .5/pow(lambda,2) * (pow(log(eta_new) - log(1 - eta_new) - mu_eta,2) - pow(log(eta) - log(1 - eta) - mu_eta,2));
          
          double accept_eta = 1.0;
          if( arma::is_finite(log_ratio_eta) ){
            if(log_ratio_eta < 0){
              accept_eta = exp(log_ratio_eta);
            }
          }else{
            accept_eta = 0.0;
          }
          
          eta_accept += accept_eta;
          eta_count ++;
          
          if( arma::randu() < accept_eta ){
            eta = eta_new;
          }
          
          s_eta += pow(it+1,-ADAPT(1))*(accept_eta - ADAPT(2));
          if(s_eta > exp(50)){
            s_eta = exp(50);
          }else{
            if(s_eta < exp(-50)){
              s_eta = exp(-50);
            }
          }
          
          
          
          
          //Update mu_eta (univariate conjugate normal)
          mu_eta = (log(eta) - log(1.0 - eta))/2 + arma::randn() * sqrt(pow(lambda,2.0)/2.0);
          
          
          
          
          //Update lambda \sim Cauchy(a)
          double lambda_new = lambda * exp(arma::randn() * sqrt(s_lambda));
          
          double log_ratio_lambda = log(a_lambda + pow(lambda,2.0)) - log(a_lambda + pow(lambda_new,2.0)) - 0.5 * pow(log(eta) - log(1.0 - eta) - mu_eta,2.0) * (1.0 / pow(lambda_new,2.0) - 1.0 / pow(lambda,2.0));
          
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
          
          s_lambda += pow(it+1,-ADAPT(1))*(accept_lambda - ADAPT(2));
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
    
    
    
    
    
    if(!is_null_X){
      /////////////////
      // Update beta // one for each process with regression coefficients
      /////////////////
      
      arma::vec randu_cache_beta(pY,arma::fill::randu);
      arma::field<arma::vec> randn_cache_beta(pY);
      for(int h = 0; h < pY; h++){
        arma::vec beta_vec_h = beta_vec[h];
        arma::mat S_beta_h = S_beta[h];
        randn_cache_beta(h) = arma::mvnrnd(beta_vec_h, S_beta_h);
      }
      
      // #pragma omp parallel for
      for(int h = 0; h < pY; h++){
        
        //Propose new value (mv-normal)
        arma::vec beta_vec_h = beta_vec[h];
        arma::vec beta_vec_new = randn_cache_beta(h);
        arma::mat beta_new(beta_vec_new);
        beta_new.reshape(g(h),dY(h)*(dY(h) - 1));
        
        //Computing MH ratio:
        arma::vec mu_beta_vec_h = mu_beta_vec[h];
        arma::mat VU_beta_inv_h = VU_beta_inv[h];
        
        //Prior (proposal is symmetric)
        double log_ratio_beta = - .5 * (arma::as_scalar((beta_vec_new.t() - mu_beta_vec_h.t()) * VU_beta_inv_h * (beta_vec_new - mu_beta_vec_h)) - arma::as_scalar((beta_vec_h.t() - mu_beta_vec_h.t()) * VU_beta_inv_h * (beta_vec_h - mu_beta_vec_h)));
        
        //Likelihood (only consider the contribution from h-th process)
        arma::mat beta_h = beta(h);
        arma::mat X_h = X[h];
        arma::mat gamma_h = gamma(h);
        std::vector<arma::mat> Z_h = Z[h];
        
        double num_loglike = 0.0, den_loglike = 0.0;
        
        //sum because the other elements (different from c(i)) will be zero
        num_loglike = logLike_h(setC_all, Y, c, phi_star, n_times.col(h), epsilon, X_h, Z_h, beta_new, gamma_h, n_rates_cum, dY, h, impute_missing_Y);
        den_loglike = logLike_h(setC_all, Y, c, phi_star, n_times.col(h), epsilon, X_h, Z_h, beta_h, gamma_h, n_rates_cum, dY, h, impute_missing_Y);
        
        log_ratio_beta += num_loglike - den_loglike;
        
        double accept_beta = 1.0;
        if( arma::is_finite(log_ratio_beta) ){
          if(log_ratio_beta < 0){
            accept_beta = exp(log_ratio_beta);
          }
        }else{
          accept_beta = 0.0;
        }
        
        beta_accept(h) += accept_beta;
        beta_count(h) ++;
        
        if( randu_cache_beta(h) < accept_beta ){
          beta_vec_h = beta_vec_new;
          beta_vec[h] = beta_vec_new;
          beta(h) = beta_new;
        }
        
        arma::vec sum_beta_h = sum_beta[h];
        sum_beta_h += beta_vec_h;
        sum_beta[h] = sum_beta_h;
        
        arma::mat prod_beta_h = prod_beta[h];
        prod_beta_h += beta_vec_h * beta_vec_h.t();
        prod_beta[h] = prod_beta_h;
        
        s_d_beta(h) += pow(it+1,-ADAPT(1))*(accept_beta - ADAPT(2));
        if(s_d_beta(h) > exp(50)){
          s_d_beta(h) = exp(50);
        }else{
          if(s_d_beta(h) < exp(-50)){
            s_d_beta(h) = exp(-50);
          }
        }
        if(it > (ADAPT(0) - 1)){
          arma::mat eye_mat_beta_h = eye_mat_beta[h];
          arma::mat S_beta_h = S_beta[h];
          S_beta_h = s_d_beta(h)/it * (prod_beta_h - sum_beta_h * sum_beta_h.t()/(it+1.0)) + s_d_beta(h) * pow(0.1,2.0) / dim_beta(h) * eye_mat_beta_h;
          // S_beta_h = s_d_beta(h)/it * (prod_beta_h - sum_beta_h * sum_beta_h.t()/(it+1.0)) + s_d_beta(h) * ADAPT(3) * eye_mat_beta_h;
          S_beta[h] = S_beta_h;
        }
      }
    }
    
    
    
    if(!is_null_Z){
      //////////////////
      // Update gamma // one for each process with regression coefficients
      //////////////////
      
      arma::vec randu_cache_gamma(pY,arma::fill::randu);
      arma::field<arma::vec> randn_cache_gamma(pY);
      for(int h = 0; h < pY; h++){
        arma::vec gamma_vec_h = gamma_vec[h];
        arma::mat S_gamma_h = S_gamma[h];
        randn_cache_gamma(h) = arma::mvnrnd(gamma_vec_h, S_gamma_h);
      }
      
      // #pragma omp parallel for
      for(int h = 0; h < pY; h++){
        //Propose new value (mv-normal)
        arma::vec gamma_vec_h = gamma_vec[h];
        arma::vec gamma_vec_new = randn_cache_gamma(h);
        arma::mat gamma_new(gamma_vec_new);
        gamma_new.reshape(q(h),dY(h)*(dY(h) - 1));
        
        //Computing MH ratio:
        arma::vec mu_gamma_vec_h = mu_gamma_vec[h];
        arma::mat VU_gamma_inv_h = VU_gamma_inv[h];
        
        //Prior (proposal is symmetric)
        double log_ratio_gamma = - .5 * (arma::as_scalar((gamma_vec_new.t() - mu_gamma_vec_h.t()) * VU_gamma_inv_h * (gamma_vec_new - mu_gamma_vec_h)) - arma::as_scalar((gamma_vec_h.t() - mu_gamma_vec_h.t()) * VU_gamma_inv_h * (gamma_vec_h - mu_gamma_vec_h)));
        
        //Likelihood (only consider the contribution from h-th process)
        arma::mat beta_h = beta(h);
        arma::mat X_h = X[h];
        arma::mat gamma_h = gamma(h);
        std::vector<arma::mat> Z_h = Z[h];
        
        double num_loglike = 0.0, den_loglike = 0.0;
        
        //sum because the other elements (different from c(i)) will be zero
        num_loglike = logLike_h(setC_all, Y, c, phi_star, n_times.col(h), epsilon, X_h, Z_h, beta_h, gamma_new, n_rates_cum, dY, h, impute_missing_Y);
        den_loglike = logLike_h(setC_all, Y, c, phi_star, n_times.col(h), epsilon, X_h, Z_h, beta_h, gamma_h, n_rates_cum, dY, h, impute_missing_Y);
        
        log_ratio_gamma += num_loglike - den_loglike;
        
        double accept_gamma = 1.0;
        if( arma::is_finite(log_ratio_gamma) ){
          if(log_ratio_gamma < 0){
            accept_gamma = exp(log_ratio_gamma);
          }
        }else{
          accept_gamma = 0.0;
        }
        
        gamma_accept(h) += accept_gamma;
        gamma_count(h) ++;
        
        if( randu_cache_gamma(h) < accept_gamma ){
          gamma_vec_h = gamma_vec_new;
          gamma_vec[h] = gamma_vec_new;
          gamma(h) = gamma_new;
        }
        
        arma::vec sum_gamma_h = sum_gamma[h];
        sum_gamma_h += gamma_vec_h;
        sum_gamma[h] = sum_gamma_h;
        
        arma::mat prod_gamma_h = prod_gamma[h];
        prod_gamma_h += gamma_vec_h * gamma_vec_h.t();
        prod_gamma[h] = prod_gamma_h;
        
        s_d_gamma(h) += pow(it+1,-ADAPT(1))*(accept_gamma - ADAPT(2));
        if(s_d_gamma(h) > exp(50)){
          s_d_gamma(h) = exp(50);
        }else{
          if(s_d_gamma(h) < exp(-50)){
            s_d_gamma(h) = exp(-50);
          }
        }
        if(it > (ADAPT(0) - 1)){
          arma::mat eye_mat_gamma_h = eye_mat_gamma[h];
          arma::mat S_gamma_h = S_gamma[h];
          S_gamma_h = s_d_gamma(h)/it * (prod_gamma_h - sum_gamma_h * sum_gamma_h.t()/(it+1.0)) + s_d_gamma(h) * pow(0.1,2.0) / dim_gamma(h) * eye_mat_gamma_h;
          // S_gamma_h = s_d_gamma(h)/it * (prod_gamma_h - sum_gamma_h * sum_gamma_h.t()/(it+1.0)) + s_d_gamma(h) * ADAPT(3) * eye_mat_gamma_h;
          S_gamma[h] = S_gamma_h;
        }
      }
    }
    
    
    
    
    //////////////////////
    // sample locations //
    //////////////////////
    
    arma::vec randu_cache_phi(N,arma::fill::randu);
    arma::mat randn_cache_phi(N,p_tot);
    for(int i = 0; i < N; i++){
      randn_cache_phi.row(i) = arma::mvnrnd(phi_star.row(i).t(), S_phi_star).t();
    }
    
    // #pragma omp parallel for
    for(int i = 0; i < N; i++){
      
      //Sample new phi_star_m
      arma::rowvec phi_star_new = randn_cache_phi.row(i);
      
      arma::uvec setC1(1);
      setC1(0) = i;
      
      //Compute likelihood for vector phi_star_j (consider all the processes now!)
      double num_loglike = 0.0, den_loglike = 0.0;
      
      num_loglike = arma::accu(logLike_all(setC1, Y, H, phi_star_new, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
      den_loglike = arma::accu(logLike_all(setC1, Y, H, phi_star.row(i), n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
      
      double log_ratio_phi = num_loglike - den_loglike;
      
      //Prior (symmetric proposal in phi_star!)
      arma::rowvec phi_star_new_minus_mu = phi_star_new - mu;
      arma::rowvec phi_star_minus_mu = phi_star.row(i) - mu;
      log_ratio_phi += - .5 * ( arma::as_scalar(phi_star_new_minus_mu * Omega * phi_star_new_minus_mu.t()) - arma::as_scalar(phi_star_minus_mu * Omega * phi_star_minus_mu.t()) );
      
      double accept_phi = 1.0;
      if( arma::is_finite(log_ratio_phi) ){
        if(log_ratio_phi < 0){
          accept_phi = exp(log_ratio_phi);
        }
      }else{
        accept_phi = 0.0;
      }
      
      if( randu_cache_phi(i) < accept_phi ){
        phi_star.row(i) = phi_star_new;
      }
    }
    
    
    //Save output for this iteration
    if( (it + 1 > (n_burn1 + n_burn2)) & (((it + 1 - n_burn1 - n_burn2) / thin - floor((it + 1 - n_burn1 - n_burn2) / thin)) == 0 )){
      
      int iter = (it + 1 - n_burn1 - n_burn2)/thin - 1;
      
      // Lists
      if(!is_null_X){
        beta_out[iter] = beta;
      }
      if(!is_null_Z){
        gamma_out[iter] = gamma;
      }
      if(impute_missing_Y){
        Y_out[iter] = Y;
      }
      if(impute_missing_H){
        H_out[iter] = H;
      }
      
      // Cubes
      phi_star_out.slice(iter) = phi_star;
      Omega_out.slice(iter) = Omega;
      G_out.slice(iter) = G;
      G0_out.slice(iter) = G0;
      
      // Mats
      mu_out.row(iter) = mu;
      m_mu_out.row(iter) = m_mu;
      
      // Vecs
      eta_out(iter) = eta;
      lambda_out(iter) = lambda;
      mu_eta_out(iter) = mu_eta;
      size_G0_out(iter) = size_G0;
      sum_weights_out(iter) = sum_weights;
      k0_out(iter) = k0;
    }
    
    //Progress bar increment
    progr.increment(); 
  }
  
  //Print acceptance rates
  if(!is_null_X){
    Rcout << "beta a.r. = " << beta_accept/beta_count << "\n";
  }
  if(!is_null_Z){
    Rcout << "gamma a.r. = " << gamma_accept/gamma_count << "\n";
  }
  if(update_eta){
    if(!eta_beta_prior){
      Rcout << "eta a.r. = " << eta_accept/eta_count << "\n";
      Rcout << "lambda a.r. = " << lambda_accept/lambda_count << "\n";
    }
  }
  
  List Graph_List = List::create(Named("Omega_out") = Omega_out, Named("G_out") = G_out, Named("G0_out") = G0_out, Named("sizeG0_out") = size_G0_out, Named("sum_weights_out") = sum_weights_out, Named("eta_out") = eta_out, Named("mu_eta_out") = mu_eta_out, Named("lambda_out") = lambda_out);
  
  List XZ_List;
  if(!is_null_X){
    if(!is_null_Z){
      XZ_List = List::create(Named("beta_out") = beta_out, Named("gamma_out") = gamma_out);
    }else{
      XZ_List = List::create(Named("beta_out") = beta_out);
    }
  }else{
    if(!is_null_Z){
      XZ_List = List::create(Named("gamma_out") = gamma_out);
    }//Else we are fitting a model without covariates...not available yet
  }
  
  if(impute_missing_Y){
    if(impute_missing_H){
      return List::create(Named("Graph_List") = Graph_List, Named("mu_out") = mu_out, Named("m_mu_out") = m_mu_out, Named("k0_out") = k0_out, Named("phi_star_out") = phi_star_out, Named("XZ_List") = XZ_List, Named("Y_out") = Y_out, Named("H_out") = H_out);
    }else{
      return List::create(Named("Graph_List") = Graph_List, Named("mu_out") = mu_out, Named("m_mu_out") = m_mu_out, Named("k0_out") = k0_out, Named("phi_star_out") = phi_star_out, Named("XZ_List") = XZ_List, Named("Y_out") = Y_out);
    }
  }else{
    if(impute_missing_H){
      return List::create(Named("Graph_List") = Graph_List, Named("mu_out") = mu_out, Named("m_mu_out") = m_mu_out, Named("k0_out") = k0_out, Named("phi_star_out") = phi_star_out, Named("XZ_List") = XZ_List, Named("H_out") = H_out);
    }else{
      return List::create(Named("Graph_List") = Graph_List, Named("mu_out") = mu_out, Named("m_mu_out") = m_mu_out, Named("k0_out") = k0_out, Named("phi_star_out") = phi_star_out, Named("XZ_List") = XZ_List);
    }
  }
}


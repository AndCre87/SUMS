// Main Gibbs sampler function for Disease Progression model (mixture)

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <omp.h>
#include <tuple>
#include "impute_missing.h"
#include "ggm_dmh.h"
#include "logLike_all.h"
#include "logLike_h.h"
#include "rgwish_c.h"
#include "log_dgamma_arma.h"
#include "log_dmvnrm_arma.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
List DisProgrGibbs_Mix(List MCMC_input){
  
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
  int p_tot = arma::accu(n_rates);
  
  //We assume general number of states
  arma::vec dY = MCMC_input["dY"];
  arma::vec dH = MCMC_input["dH"];
  
  
  // Extract hyperparameters
  List Param_list = MCMC_input["Param_list"];
  double nu = Param_list["nu"];
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
  double n_edges0 = p0 * (p0 - 1) / 2, size_G0 = 0;
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
  arma::field<arma::mat> gamma(pY); //Only these are arma::field's! For saving them in output (List)
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
  
  //Initialize BNP lists
  double a1 = 0.0, b1 = 0.0, gamma_S = 0.0, s_gamma_S = 0.01, gamma_S_accept = 0.0, gamma_S_count = 0.0;
  double a2 = 0.0, b2 = 0.0, Lambda = 0.0;
  bool update_gamma_S = false, update_Lambda = false;
  
  gamma_S = Param_list["gamma_S"];
  update_gamma_S = Param_list["update_gamma_S"];
  if(update_gamma_S){
    a1 = Param_list["a1"];
    b1 = Param_list["b1"];
  }
  
  Lambda = Param_list["Lambda"];
  update_Lambda = Param_list["update_Lambda"];
  if(update_Lambda){
    a2 = Param_list["a2"];
    b2 = Param_list["b2"];
  }
  
  
  //Algorithm specs
  List Alg_list = MCMC_input["Alg_list"];
  
  //We initialize based on s_init
  bool update_c = Alg_list["update_c"];
  arma::vec c = MCMC_input["c_init"], c_star = unique(c); //allocation variables
  int K_N = c_star.n_elem, M = K_N, M_na = M - K_N;
  double u = 1.0, T_sum = 1.0;
  arma::vec nj(M,arma::fill::zeros), S_m(M,arma::fill::ones); // cluster sizes and weights
  for(int i = 0; i < N; i++){
    nj(c(i)) ++;
  }
  //Initialize locations
  arma::mat phi_star = arma::mvnrnd(mu.t(), Sigma, M).t();
  
  //For Split & Merge
  bool SM_alg = Alg_list["SM_alg"];
  int N2 = N*(N - 1)/2;
  arma::mat lambda_dist_mat = MCMC_input["lambda_dist_mat"], eye_mat_SM(p_tot,p_tot,arma::fill::eye);
  arma::vec prob_ij(N2);
  prob_ij.fill(1.0/N2);
  arma::vec prob_ij_cum = arma::cumsum(prob_ij);
  double SM_accept = 0, s_SM = Alg_list["s_SM"];
  int SM_nSplit = 0, SM_nMerge = 0, SM_count = 0;
  bool use_Sigma = Alg_list["use_Sigma"], update_phij = Alg_list["update_phij"];
  
  // Iterations after first burn-in
  double n_burn1 = Alg_list["n_burn1"], n_burn2 = Alg_list["n_burn2"], thin = Alg_list["thin"], n_save = Alg_list["n_save"], threshold = Alg_list["threshold"];
  int n_tot = n_burn1 + n_burn2 + thin*n_save, n_edges = Alg_list["n_edges"], Gibbs_its = Alg_list["Gibbs_its"];
  
  //To monitor allocation probabilities
  int n_nan_ci = 0, n_inf_ci = 0;
  
  // Adaptation
  arma::vec ADAPT(4);
  ADAPT(0) = n_burn1; //"burn-in" for adaptation
  ADAPT(1) = 0.7; //exponent for adaptive step
  ADAPT(2) = 0.234; //reference acceptance rate
  ADAPT(3) = 0.001; //for multivariate updates
  
  //Binder
  arma::uvec Binder_index(1);
  double Binder_const = 0.5;
  arma::rowvec c_out_row(N), Binder_f(n_save), Binder_est(N);
  arma::mat pij(N,N), aux_mat(N,N);
  
  
  // Output lists
  List phi_star_out(n_save), S_m_out(n_save), beta_out(n_save), gamma_out(n_save), Y_out(n_save), H_out(n_save);
  arma::cube Omega_out(p_tot,p_tot,n_save,arma::fill::zeros), G_out(p_tot,p_tot,n_save,arma::fill::zeros), G0_out(p0,p0,n_save,arma::fill::zeros);
  arma::mat mu_out(n_save, p_tot, arma::fill::zeros), m_mu_out(n_save, p_tot, arma::fill::zeros), c_out(n_save,N,arma::fill::zeros);
  arma::vec d_out(n_save, arma::fill::zeros), size_G0_out(n_save, arma::fill::zeros), sum_weights_out(n_save,arma::fill::zeros), u_out(n_save, arma::fill::zeros), K_N_out(n_save,arma::fill::zeros), M_out(n_save,arma::fill::zeros), Lambda_out(n_save,arma::fill::zeros), gamma_S_out(n_save,arma::fill::zeros), k0_out(n_save, arma::fill::zeros), lambda_out(n_save,arma::fill::zeros), mu_d_out(n_save,arma::fill::zeros), Entropy_out(n_save, arma::fill::zeros);
  
  // Main Gibbs
  Progress progr(n_tot, true);
  
  for(int it = 0; it < n_tot; it++){
    // Rcout << "it = " << it << "\n";
    
    
    // Rcout << "Missing\n";
    ///////////////////////////
    // Impute missing values //
    ///////////////////////////
    std::tie(Y, H) = impute_missing(Y, is_NA_Y, H, is_NA_H, c, phi_star, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H);
    
    
    
    
    
    
    // Rcout << "Omega, G\n";
    ///////////////////////////
    // Update G0, G, Omega_G //
    ///////////////////////////
    
    //Marginal update
    arma::rowvec phi_star_mean = arma::sum(phi_star) / M;
    arma::mat Psi_star = Psi + k0 * M / (k0 + M) * (m_mu.t() - phi_star_mean.t()) * (m_mu - phi_star_mean);
    for(int m = 0; m < M; m++){
      Psi_star += (phi_star.row(m).t()  - phi_star_mean.t()) * (phi_star.row(m)  - phi_star_mean);
    }
    
    //Cholesky gives upper tri here
    arma::mat Ts = arma::chol( arma::inv_sympd( Psi_star ) );
    // nu_star = nu + M + 1;
    double nu_star = nu + M;
    
    std::tie(Omega, G, G0, size_G0, sum_weights) = ggm_DMH(G, G0, d, size_based_prior, a_d, b_d, size_G0, Ts, Ti, Omega, n_rates_cum, threshold, nu, nu_star, Psi, Psi_star, n_edges);
    
    
    
    
    
    
    
    // Rcout << "Omega\n";
    ////////////////////////////
    // Update Omega_G (extra) //
    ////////////////////////////
    
    Omega = rgwish_c( nu_star, Ts, G, threshold );
    Sigma = arma::inv_sympd(Omega);
    
    
    
    
    // Rcout << "mu\n";
    ///////////////
    // Update mu //
    ///////////////
    
    // Conditionally independent prior on mu
    if(update_mu){
      arma::rowvec sum_of_phi_star = arma::sum(phi_star);
      arma::mat var_mu = arma::inv_sympd((M + k0) * Omega);
      arma::vec mean_mu = (k0 * m_mu.t() + sum_of_phi_star.t())/(M + k0);
      mu = arma::mvnrnd(mean_mu, var_mu).t();
    }
    
    
    
    
    // Rcout << "m_mu\n";
    /////////////////
    // Update m_mu //
    /////////////////
    
    //Additional prior on mean of mu
    if(update_m_mu){
      arma::mat var_m_mu = arma::inv_sympd(Prec_m_mu + k0 * Omega);
      arma::vec mean_m_mu = var_m_mu * k0 * Omega * mu.t();
      m_mu = arma::mvnrnd(mean_m_mu, var_m_mu).t();
    }
    
    
    
    
    // Rcout << "k0\n";
    ///////////////
    // Update k0 //
    ///////////////
    
    //Additional prior on k0
    if(update_k0){
      double k0_aux = b_k0 + .5 * as_scalar((mu - m_mu) * Omega * (mu.t() - m_mu.t()));
      k0 = arma::randg( arma::distr_param( a_k0 + p_tot/2, 1/k0_aux ));
    }
    
    
    
    
    // Rcout << "d\n";
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
          
          s_d += pow(it+1,-ADAPT(1))*(accept_d - ADAPT(2));
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
      // Rcout << "beta\n";
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
      // Rcout << "gamma\n";
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
    
    
    
    // Rcout << "u\n";
    //////////////
    // Sample U //
    //////////////
    
    T_sum = sum(S_m);
    u = arma::randg( arma::distr_param( N/1.0, 1.0/T_sum ));
    
    
    
    if(update_c){
      // Rcout << "ci\n";
      ////////////////////////
      // Sample allocations //
      ////////////////////////
      
      arma::vec randu_cache_c(N,arma::fill::randu);
      
      // #pragma omp parallel for
      for(int i = 0; i < N; i++){
        
        arma::uvec setC1(1);
        setC1(0) = i;
        
        arma::vec loglike_vec;
        
        loglike_vec = logLike_all(setC1, Y, H, phi_star, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H);
        
        arma::vec prob_m(M);
        
        if(loglike_vec.has_nan()){//Likelihood = 0
          loglike_vec.replace(arma::datum::nan, - arma::datum::inf);
          n_nan_ci ++;
        }
        if(loglike_vec.has_inf()){//Replace with current max (not inf)
          loglike_vec.replace(arma::datum::inf, max(loglike_vec.elem( arma::find_finite(loglike_vec) )));
          n_inf_ci ++;
        }
        prob_m = log(S_m) + loglike_vec;
        
        prob_m = exp(prob_m - max(prob_m));
        prob_m = prob_m/sum(prob_m);
        
        arma::vec prob_m_cum = arma::cumsum(prob_m);
        arma::uvec hh_vec = arma::find(prob_m_cum >= randu_cache_c(i), 1, "first");
        c(i) = hh_vec(0);
      }
      c_star = arma::unique(c); //This returns a vector of unique elements in ascending order!
      K_N = c_star.n_elem;
      M_na = M - K_N;
      
      
      
      // Rcout << "Re-order\n";
      //////////////
      // Re-order //
      //////////////
      arma::vec c_aux = c;
      c_aux.fill(-1);
      
      // #pragma omp parallel for
      for(int m = 0; m < K_N; m ++){
        int c_m = c_star(m);
        
        arma::uvec setCm = arma::find(c == c_m);
        nj(m) = setCm.n_elem;
        
        for(int ind = 0; ind < nj(m); ind++){
          c_aux(setCm(ind)) = m;
        }
      }
      c = c_aux;
      c_star = arma::unique(c); //This returns a vector of unique elements in ascending order!
    }
    
    
    
    // Rcout << "Alloc\n";
    ///////////////////////////////////////////////
    // sample unnormalized weights and locations //
    ///////////////////////////////////////////////
    
    // Allocated components
    arma::vec randu_cache_phi(K_N,arma::fill::randu);
    arma::mat randn_cache_phi(K_N,p_tot);
    for(int m = 0; m < K_N; m++){
      S_m(m) = arma::randg( arma::distr_param( nj(m) + gamma_S, 1.0 / (1.0 + u) ));
      randn_cache_phi.row(m) = arma::mvnrnd(phi_star.row(m).t(), S_phi_star).t();
    }
    
    // #pragma omp parallel for
    for(int m = 0; m < K_N; m++){
      //Locations are not conjugate
      //Sample new phi_star_m
      arma::rowvec phi_star_m_curr = phi_star.row(m);
      arma::rowvec phi_star_m_new = randn_cache_phi.row(m);
      
      arma::uvec setCm = arma::find(c == m);
      
      //Compute likelihood for vector phi_star_j
      double num_loglike = 0.0, den_loglike = 0.0;
      
      num_loglike = arma::accu(logLike_all(setCm, Y, H, phi_star_m_new, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
      den_loglike = arma::accu(logLike_all(setCm, Y, H, phi_star_m_curr, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
      
      double log_ratio_alloc = num_loglike - den_loglike;
      
      //Prior (symmetric proposal in phi_star!)
      arma::rowvec phi_star_new_minus_mu = phi_star_m_new - mu;
      arma::rowvec phi_star_minus_mu = phi_star_m_curr - mu;
      log_ratio_alloc += - .5 * ( arma::as_scalar(phi_star_new_minus_mu * Omega * phi_star_new_minus_mu.t()) - arma::as_scalar(phi_star_minus_mu * Omega * phi_star_minus_mu.t()) );
      
      double accept_alloc = 1.0;
      if( arma::is_finite(log_ratio_alloc) ){
        if(log_ratio_alloc < 0){
          accept_alloc = exp(log_ratio_alloc);
        }
      }else{
        accept_alloc = 0.0;
      }
      
      if( randu_cache_phi(m) < accept_alloc ){
        phi_star.row(m) = phi_star_m_new;
      }
    }
    
    
    
    
    // Rcout << "non-Alloc\n";
    ///////////////////////////////////////////////////////
    // Sample M-K_N (number of non-allocated components) //
    ///////////////////////////////////////////////////////
    double lphi_u = - gamma_S * log(1.0 + u);
    double Lambda_u = exp(log(Lambda) + lphi_u);
    double M_na_new = R::rpois(Lambda_u);
    double pi_M_na = Lambda_u/(K_N + Lambda_u);
    if( arma::randu() < pi_M_na ){
      M_na_new ++;
    }
    
    int M_new = K_N + M_na_new;
    
    S_m.resize(M_new);
    phi_star.resize(M_new,p_tot);
    nj.resize(M_new);
    
    M_na = M_na_new;
    M = M_new;
    
    
    // Non-Allocated components
    for(int m = K_N; m < M; m ++){
      S_m(m) = arma::randg( arma::distr_param( gamma_S, 1/(1 + u) ));
      //Sample new phi_star_m from prior
      phi_star.row(m) = arma::mvnrnd(mu.t(), Sigma).t();
    }
    
    
    
    
    
    ////////////////////
    // Update gamma_S //
    ////////////////////
    if(update_gamma_S){
      // Rcout << "gamma_S\n";
      //Propose new value
      double gamma_S_new = gamma_S * exp(R::rnorm(0,1) * sqrt(s_gamma_S));
      
      //Prior and proposal part
      double log_ratio_gamma_S = a1 *(log(gamma_S_new) - log(gamma_S)) - b1 * (gamma_S_new - gamma_S);
      //Likelihood part
      log_ratio_gamma_S += log(Lambda / pow(1 + u, gamma_S_new) + K_N) - log(Lambda / pow(1 + u, gamma_S) + K_N);
      log_ratio_gamma_S += Lambda * (1 / pow(1 + u, gamma_S_new) - 1 / pow(1 + u, gamma_S));
      log_ratio_gamma_S += - K_N * log(1 + u) * (gamma_S_new - gamma_S);
      log_ratio_gamma_S += arma::accu(lgamma(gamma_S_new + nj.rows(0,K_N-1)) - lgamma(gamma_S + nj.rows(0,K_N-1))) - K_N * (lgamma(gamma_S_new) - lgamma(gamma_S));
      
      double accept_gamma_S = 1.0;
      if( arma::is_finite(log_ratio_gamma_S) ){
        if(log_ratio_gamma_S < 0){
          accept_gamma_S = exp(log_ratio_gamma_S);
        }
      }else{
        accept_gamma_S = 0.0;
      }
      
      gamma_S_accept += accept_gamma_S;
      gamma_S_count ++;
      
      if( arma::randu() < accept_gamma_S ){
        gamma_S = gamma_S_new;
      }
      
      s_gamma_S = s_gamma_S + pow(it+1,-ADAPT(1)) * (accept_gamma_S - ADAPT(2));
      if(s_gamma_S > exp(50)){
        s_gamma_S = exp(50);
      }else{
        if(s_gamma_S < exp(-50)){
          s_gamma_S = exp(-50);
        }
      }
    }
    
    
    
    
    ///////////////////
    // Update Lambda //
    ///////////////////
    if(update_Lambda){
      // Rcout << "Lambda\n";
      double lphi_u = - gamma_S * log(1.0 + u);
      double phi_b = 1 + b2 - exp(lphi_u);
      double a_K_N_minus1 = a2 + K_N - 1;
      
      double pi_Lambda_num = exp(lphi_u + log(a_K_N_minus1));
      double pi_Lambda_den = pi_Lambda_num + K_N * phi_b;
      double pi_Lambda = pi_Lambda_num / pi_Lambda_den;
      
      if(arma::randu() < pi_Lambda){
        Lambda = arma::randg( arma::distr_param( a_K_N_minus1 + 1, 1/phi_b ));
      }else{
        Lambda = arma::randg( arma::distr_param( a_K_N_minus1, 1/phi_b ));
      }
    }
    
    
    
    
    if(update_c){
      //////////////////////////
      // Split and Merge step //
      //////////////////////////
      if((SM_alg) & (it % Gibbs_its == 0)){
        
        // Rcout << "SM alg\n";
        
        //Select pair of indices at random from matrix disposition (upper-tri)
        arma::uvec index_vec = arma::find(prob_ij_cum >= arma::randu(), 1, "first");
        int index_SM = index_vec(0);
        
        int counter = 0;
        while(index_SM >= 0){
          counter = counter + 1;
          index_SM -= counter;
        }
        int i = index_SM + counter;
        int j = counter;
        // Rcout << "(i,j) = (" << i << "," << j << ")\n";
        
        //Extract current allocation variables
        int c_i = c(i);
        int c_j = c(j);
        // Rcout << "(c_i,c_j) = (" << c_i << "," << c_j << ")\n";
        
        //Intialisation of new vectors
        arma::vec c_new = c;
        arma::vec nj_new = nj;
        arma::vec S_m_new = S_m;
        arma::mat phi_star_new = phi_star;
        
        double log_ratio_SM = 0.0;
        
        bool is_Split = true;
        
        if(c_i == c_j){//Split
          
          // Rcout << "Split\n";
          
          //Current cluster label and values
          int c_c = c_j;
          
          double S_m_c = S_m(c_c);
          arma::rowvec phi_star_c = phi_star.row(c_c);
          
          //Find indices excepts (i,j)
          arma::uvec setC_minusij = arma::find(c == c_c);
          setC_minusij.shed_rows(arma::find(setC_minusij == i)); //remove i-th element
          setC_minusij.shed_rows(arma::find(setC_minusij == j)); //remove j-th element
          int n_setC_minusij = setC_minusij.n_elem;
          
          //Allocate new component
          //New label 
          int c_i_new = K_N;
          c_new(i) = c_i_new;
          int c_j_new = c_c;
          c_new(j) = c_j_new;
          
          //Start by putting one element in each cluster
          nj_new.insert_rows(c_i_new,1);
          nj_new(c_i_new) = 1.0;
          nj_new(c_j_new) = 1.0;
          
          //Propose allocations
          //Exp of minus the distance from the "centres" (yi, yj) of the new clusters / 2.0
          //Distances calculated outside using Optimal Matching Analysis
          for(int ind =  0; ind < n_setC_minusij; ind++){
            
            int i2 = setC_minusij(ind);
            
            arma::vec prob_2(2,arma::fill::zeros);
            prob_2(0) =  - lambda_dist_mat(i2,i) / 2.0;
            prob_2(1) =  - lambda_dist_mat(i2,j) / 2.0;
            
            prob_2 = exp(prob_2-max(prob_2));
            prob_2 = prob_2/sum(prob_2);
            
            arma::vec prob_2_cum = arma::cumsum(prob_2);
            arma::uvec hh_vec = arma::find(prob_2_cum >= arma::randu(), 1, "first");
            int hh = hh_vec(0);
            
            c_new(i2) = c_i_new * (1 - hh) + c_j_new * hh;
            nj_new(c_new(i2))++;
            
            log_ratio_SM -= log(prob_2(hh));
          }
          // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
          
          
          // New component for unnormalised weigths and for latent variables
          S_m_new.insert_rows(c_i_new,1);
          phi_star_new.insert_rows(c_i_new, 1);
          
          if(update_phij){
            for(int m = 0; m < 2; m++){
              int c_m = c_i_new * (1 - m) + c_j_new * m;
              
              arma::uvec setCm = arma::find(c_new == c_m);
              double nj_m = setCm.n_elem;
              
              S_m_new(c_m) = arma::randg( arma::distr_param( gamma_S + nj_m, 1.0 / (1.0 + u) ));
              log_ratio_SM -= log_dgamma_arma(S_m_new(c_m), gamma_S + nj_m, 1.0 + u);
              
              // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
              
              //Compute mean vector for proposal: estimator of log-transition rates in the cluster
              arma::rowvec phi_prop_mean(p_tot,arma::fill::zeros);
              
              //Response processes
              for(int h = 0; h < pY; h ++){
                // Rcout << "h = " << h << "\n";
                
                int dY_h = dY[h];
                arma::mat trans_mat(dY_h,dY_h,arma::fill::zeros);
                arma::vec sojourn_times(dY_h,arma::fill::zeros);
                
                for(int ind = 0; ind < nj_m; ind++){
                  
                  int subj = setCm(ind);
                  std::vector<arma::vec> Y_i = Y[subj];
                  arma::vec Y_i_h = Y_i[h];
                  
                  //Number of times
                  int n_times_i_h = n_times(subj,h);
                  
                  //Time intervals
                  std::vector<arma::vec> epsilon_i = epsilon[subj];
                  arma::vec epsilon_i_h = epsilon_i[h];
                  
                  //Value of response at time 0
                  double Yihj_minus_1 = Y_i_h(0);
                  
                  for(int jt = 1; jt < n_times_i_h; jt++){
                    
                    //Value of response at time j
                    double Yihj = Y_i_h(jt);
                    
                    if(Yihj_minus_1 != Yihj){//We do not consider transitions to the same state
                      trans_mat(Yihj_minus_1,Yihj) ++;
                    }
                    sojourn_times(Yihj_minus_1) += epsilon_i_h(jt-1);
                  }
                }
                
                //Fill-in components of mean vector for proposal
                arma::uvec ind_dY = arma::regspace<arma::uvec>(n_rates_cum(h),n_rates_cum(h+1)-1);
                int count_Y = 0;
                for(int i2 = 0; i2< dY_h; i2 ++){
                  for(int j2 = 0; j2 < dY_h; j2++){
                    if(i2 != j2){
                      if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                        phi_prop_mean(ind_dY(count_Y)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                      }else{
                        phi_prop_mean(ind_dY(count_Y)) = mu(ind_dY(count_Y));
                      }
                      count_Y++;
                    }
                  }
                }
              }
              
              //Covariate processes
              for(int l = 0; l < pH; l ++){
                // Rcout << "l = " << l << "\n";
                
                int dH_l = dH[l];
                arma::mat trans_mat(dH_l,dH_l,arma::fill::zeros);
                arma::vec sojourn_times(dH_l,arma::fill::zeros);
                
                for(int ind = 0; ind < nj_m; ind++){
                  
                  int subj = setCm(ind);
                  std::vector<arma::vec> H_i = H[subj];
                  arma::vec H_i_l = H_i[l];
                  
                  //Number of times
                  int n_times_i_l = n_times(subj,pY+l);
                  
                  //Time intervals
                  std::vector<arma::vec> epsilon_i = epsilon[subj];
                  arma::vec epsilon_i_l = epsilon_i[pY+l];
                  
                  //Value of response at time 0
                  double Hilj_minus_1 = H_i_l(0);
                  
                  for(int jt = 1; jt < n_times_i_l; jt++){
                    
                    //Value of response at time j
                    double Hilj = H_i_l(jt);
                    
                    if(Hilj_minus_1 != Hilj){//We do not consider transitions to the same state
                      trans_mat(Hilj_minus_1,Hilj) ++;
                    }
                    sojourn_times(Hilj_minus_1) += epsilon_i_l(jt-1);
                  }
                }
                
                //Fill-in components of mean vector for proposal
                arma::uvec ind_dH = arma::regspace<arma::uvec>(n_rates_cum(pY+l),n_rates_cum(pY+l+1)-1);
                int count_H = 0;
                for(int i2 = 0; i2< dH_l; i2 ++){
                  for(int j2 = 0; j2 < dH_l; j2++){
                    if(i2 != j2){
                      if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                        phi_prop_mean(ind_dH(count_H)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                      }else{
                        phi_prop_mean(ind_dH(count_H)) = mu(ind_dH(count_H));
                      }
                      count_H++;
                    }
                  }
                }
              }
              
              
              // Rcout << "phi_prop_mean = " << phi_prop_mean << "\n";
              
              //Sample new phi_star_m using a Gaussian random-walk centred on the estimator in the cluster and with covariance:
              if(use_Sigma){//Sigma
                phi_star_new.row(c_m) = arma::mvnrnd(phi_prop_mean.t(), Sigma).t();
                log_ratio_SM -= log_dmvnrm_arma(phi_star_new.row(c_m), phi_prop_mean, Omega);
                // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
              }else{//or identity with small diagonal values
                phi_star_new.row(c_m) = arma::mvnrnd(phi_prop_mean.t(), eye_mat_SM * s_SM).t();
                log_ratio_SM -= log_dmvnrm_arma(phi_star_new.row(c_m), phi_prop_mean, eye_mat_SM / s_SM);
                // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
              }
            }
          }else{ //Only consider changes in new cluster
            int c_m = c_i_new;
            
            arma::uvec setCm = arma::find(c_new == c_m);
            double nj_m = setCm.n_elem;
            
            S_m_new(c_m) = arma::randg( arma::distr_param( gamma_S + nj_m, 1.0 / (1.0 + u) ));
            log_ratio_SM -= log_dgamma_arma(S_m_new(c_m), gamma_S + nj_m, 1.0 + u);
            
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            
            //Compute mean vector for proposal: estimator of log-transition rates in the cluster
            arma::rowvec phi_prop_mean(p_tot,arma::fill::zeros);
            
            //Response processes
            for(int h = 0; h < pY; h ++){
              // Rcout << "h = " << h << "\n";
              
              int dY_h = dY[h];
              arma::mat trans_mat(dY_h,dY_h,arma::fill::zeros);
              arma::vec sojourn_times(dY_h,arma::fill::zeros);
              
              for(int ind = 0; ind < nj_m; ind++){
                
                int subj = setCm(ind);
                std::vector<arma::vec> Y_i = Y[subj];
                arma::vec Y_i_h = Y_i[h];
                
                //Number of times
                int n_times_i_h = n_times(subj,h);
                
                //Time intervals
                std::vector<arma::vec> epsilon_i = epsilon[subj];
                arma::vec epsilon_i_h = epsilon_i[h];
                
                //Value of response at time 0
                double Yihj_minus_1 = Y_i_h(0);
                
                for(int jt = 1; jt < n_times_i_h; jt++){
                  
                  //Value of response at time j
                  double Yihj = Y_i_h(jt);
                  
                  if(Yihj_minus_1 != Yihj){//We do not consider transitions to the same state
                    trans_mat(Yihj_minus_1,Yihj) ++;
                  }
                  sojourn_times(Yihj_minus_1) += epsilon_i_h(jt-1);
                }
              }
              
              //Fill-in components of mean vector for proposal
              arma::uvec ind_dY = arma::regspace<arma::uvec>(n_rates_cum(h),n_rates_cum(h+1)-1);
              int count_Y = 0;
              for(int i2 = 0; i2< dY_h; i2 ++){
                for(int j2 = 0; j2 < dY_h; j2++){
                  if(i2 != j2){
                    if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                      phi_prop_mean(ind_dY(count_Y)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                    }else{
                      phi_prop_mean(ind_dY(count_Y)) = mu(ind_dY(count_Y));
                    }
                    count_Y++;
                  }
                }
              }
            }
            
            //Covariate processes
            for(int l = 0; l < pH; l ++){
              // Rcout << "l = " << l << "\n";
              
              int dH_l = dH[l];
              arma::mat trans_mat(dH_l,dH_l,arma::fill::zeros);
              arma::vec sojourn_times(dH_l,arma::fill::zeros);
              
              for(int ind = 0; ind < nj_m; ind++){
                
                int subj = setCm(ind);
                std::vector<arma::vec> H_i = H[subj];
                arma::vec H_i_l = H_i[l];
                
                //Number of times
                int n_times_i_l = n_times(subj,pY+l);
                
                //Time intervals
                std::vector<arma::vec> epsilon_i = epsilon[subj];
                arma::vec epsilon_i_l = epsilon_i[pY+l];
                
                //Value of response at time 0
                double Hilj_minus_1 = H_i_l(0);
                
                for(int jt = 1; jt < n_times_i_l; jt++){
                  
                  //Value of response at time j
                  double Hilj = H_i_l(jt);
                  
                  if(Hilj_minus_1 != Hilj){//We do not consider transitions to the same state
                    trans_mat(Hilj_minus_1,Hilj) ++;
                  }
                  sojourn_times(Hilj_minus_1) += epsilon_i_l(jt-1);
                }
              }
              
              //Fill-in components of mean vector for proposal
              arma::uvec ind_dH = arma::regspace<arma::uvec>(n_rates_cum(pY+l),n_rates_cum(pY+l+1)-1);
              int count_H = 0;
              for(int i2 = 0; i2< dH_l; i2 ++){
                for(int j2 = 0; j2 < dH_l; j2++){
                  if(i2 != j2){
                    if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                      phi_prop_mean(ind_dH(count_H)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                    }else{
                      phi_prop_mean(ind_dH(count_H)) = mu(ind_dH(count_H));
                    }
                    count_H++;
                  }
                }
              }
            }          
            
            // Rcout << "phi_prop_mean = " << phi_prop_mean << "\n";
            
            //Sample new phi_star_m using a Gaussian random-walk centred on the estimator in the cluster and with covariance:
            if(use_Sigma){//Sigma
              phi_star_new.row(c_m) = arma::mvnrnd(phi_prop_mean.t(), Sigma).t();
              log_ratio_SM -= log_dmvnrm_arma(phi_star_new.row(c_m), phi_prop_mean, Omega);
              // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            }else{//or identity with small diagonal values
              phi_star_new.row(c_m) = arma::mvnrnd(phi_prop_mean.t(), eye_mat_SM * s_SM).t();
              log_ratio_SM -= log_dmvnrm_arma(phi_star_new.row(c_m), phi_prop_mean, eye_mat_SM / s_SM);
              // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            }
            
          }
          
          //Numerator of the proposal is a reversed Merge step
          if(update_phij){//Only if we update the old label we need to consider this, otherwise prob = 1
            arma::uvec setCc = arma::find(c == c_c);
            double nj_c = setCc.n_elem;
            
            //Unnormalised weights
            log_ratio_SM += log_dgamma_arma(S_m_c, gamma_S + nj_c, 1.0 + u);
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            
            //Current phi_star_c centred on estimator of log-transition rates in the cluster
            arma::rowvec phi_prop_mean(p_tot,arma::fill::zeros);
            
            //Response processes
            for(int h = 0; h < pY; h ++){
              // Rcout << "h = " << h << "\n";
              
              int dY_h = dY[h];
              arma::mat trans_mat(dY_h,dY_h,arma::fill::zeros);
              arma::vec sojourn_times(dY_h,arma::fill::zeros);
              
              for(int ind = 0; ind < nj_c; ind++){
                
                int subj = setCc(ind);
                std::vector<arma::vec> Y_i = Y[subj];
                arma::vec Y_i_h = Y_i[h];
                
                //Number of times
                int n_times_i_h = n_times(subj,h);
                
                //Time intervals
                std::vector<arma::vec> epsilon_i = epsilon[subj];
                arma::vec epsilon_i_h = epsilon_i[h];
                
                //Value of response at time 0
                double Yihj_minus_1 = Y_i_h(0);
                
                for(int jt = 1; jt < n_times_i_h; jt++){
                  
                  //Value of response at time j
                  double Yihj = Y_i_h(jt);
                  
                  if(Yihj_minus_1 != Yihj){//We do not consider transitions to the same state
                    trans_mat(Yihj_minus_1,Yihj) ++;
                  }
                  sojourn_times(Yihj_minus_1) += epsilon_i_h(jt-1);
                }
              }
              
              //Fill-in components of mean vector for proposal
              arma::uvec ind_dY = arma::regspace<arma::uvec>(n_rates_cum(h),n_rates_cum(h+1)-1);
              int count_Y = 0;
              for(int i2 = 0; i2< dY_h; i2 ++){
                for(int j2 = 0; j2 < dY_h; j2++){
                  if(i2 != j2){
                    if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                      phi_prop_mean(ind_dY(count_Y)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                    }else{
                      phi_prop_mean(ind_dY(count_Y)) = mu(ind_dY(count_Y));
                    }
                    count_Y++;
                  }
                }
              }
            }
            
            //Covariate processes
            for(int l = 0; l < pH; l ++){
              // Rcout << "l = " << l << "\n";
              
              int dH_l = dH[l];
              arma::mat trans_mat(dH_l,dH_l,arma::fill::zeros);
              arma::vec sojourn_times(dH_l,arma::fill::zeros);
              
              for(int ind = 0; ind < nj_c; ind++){
                
                int subj = setCc(ind);
                std::vector<arma::vec> H_i = H[subj];
                arma::vec H_i_l = H_i[l];
                
                //Number of times
                int n_times_i_l = n_times(subj,pY+l);
                
                //Time intervals
                std::vector<arma::vec> epsilon_i = epsilon[subj];
                arma::vec epsilon_i_l = epsilon_i[pY+l];
                
                //Value of response at time 0
                double Hilj_minus_1 = H_i_l(0);
                
                for(int jt = 1; jt < n_times_i_l; jt++){
                  
                  //Value of response at time j
                  double Hilj = H_i_l(jt);
                  
                  if(Hilj_minus_1 != Hilj){//We do not consider transitions to the same state
                    trans_mat(Hilj_minus_1,Hilj) ++;
                  }
                  sojourn_times(Hilj_minus_1) += epsilon_i_l(jt-1);
                }
              }
              
              //Fill-in components of mean vector for proposal
              arma::uvec ind_dH = arma::regspace<arma::uvec>(n_rates_cum(pY+l),n_rates_cum(pY+l+1)-1);
              int count_H = 0;
              for(int i2 = 0; i2< dH_l; i2 ++){
                for(int j2 = 0; j2 < dH_l; j2++){
                  if(i2 != j2){
                    if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                      phi_prop_mean(ind_dH(count_H)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                    }else{
                      phi_prop_mean(ind_dH(count_H)) = mu(ind_dH(count_H));
                    }
                    count_H++;
                  }
                }
              }
            }
            
            // Rcout << "phi_star_c = " << phi_star_c << "\n";
            // Rcout << "phi_prop_mean = " << phi_prop_mean << "\n";
            
            //Current phi_star_c proposed from a Gaussian centred on the estimator in the cluster and with covariance:
            if(use_Sigma){//Sigma
              log_ratio_SM += log_dmvnrm_arma(phi_star_c, phi_prop_mean, Omega);
              // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            }else{//or identity with small diagonal values
              log_ratio_SM += log_dmvnrm_arma(phi_star_c, phi_prop_mean, eye_mat_SM / s_SM);
              // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            }
          }
          
          
          
          //Prior part with added new component
          if(update_phij){
            log_ratio_SM += - u * (S_m_new(c_i_new) + S_m_new(c_j_new) - S_m_c) + nj_new(c_i_new) * log(S_m_new(c_i_new)) + nj_new(c_j_new) * log(S_m_new(c_j_new)) - nj(c_c) * log(S_m_c) + log(Lambda) - log(M);
            log_ratio_SM += log_dgamma_arma(S_m_new(c_i_new), gamma_S, 1.0) + log_dgamma_arma(S_m_new(c_j_new), gamma_S, 1.0) - log_dgamma_arma(S_m_c, gamma_S, 1.0);
            log_ratio_SM += log_dmvnrm_arma(phi_star_new.row(c_i_new), mu, Omega) + log_dmvnrm_arma(phi_star_new.row(c_j_new), mu, Omega) - log_dmvnrm_arma(phi_star_c, mu, Omega);
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
          }else{
            log_ratio_SM += - u * S_m_new(c_i_new) + nj_new(c_i_new) * (log(S_m_new(c_i_new)) - log(S_m_c)) + log(Lambda) - log(M);
            log_ratio_SM += log_dgamma_arma(S_m_new(c_i_new), gamma_S, 1.0);
            log_ratio_SM += log_dmvnrm_arma(phi_star_new.row(c_i_new), mu, Omega);
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
          }
          
          
          if(update_phij){
            //Loglike for new set Ci
            arma::uvec setC_i_new = arma::find(c_new == c_i_new);
            log_ratio_SM += arma::accu(logLike_all(setC_i_new, Y, H, phi_star_new.row(c_i_new), n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            
            //Loglike for new set Cj
            arma::uvec setC_j_new = arma::find(c_new == c_j_new);
            log_ratio_SM += arma::accu(logLike_all(setC_j_new, Y, H, phi_star_new.row(c_j_new), n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            
            //Loglike for old set C
            arma::uvec setC_ij = arma::find(c == c_c);
            log_ratio_SM -= arma::accu(logLike_all(setC_ij, Y, H, phi_star_c, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
          }else{
            //Loglike part
            arma::uvec setC_i_new = arma::find(c_new == c_i_new);
            log_ratio_SM += arma::accu(logLike_all(setC_i_new, Y, H, phi_star_new.row(c_i_new), n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            log_ratio_SM -= arma::accu(logLike_all(setC_i_new, Y, H, phi_star_c, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
          }
          
        }else{//Merge
          // Rcout << "Merge\n";
          
          is_Split = false;
          
          //Current cluster label and values
          double S_m_i = S_m(c_i);
          double S_m_j = S_m(c_j);
          arma::rowvec phi_star_i = phi_star.row(c_i);
          arma::rowvec phi_star_j = phi_star.row(c_j);
          
          //Elements in set Ci now have same allocation as set Cj
          for(int i2 = 0; i2 < N; i2++){
            if(c_new(i2) == c_i){
              c_new(i2) = c_j;
            }
          }
          
          //Update allocation variable after removing one cluster
          for(int i2 = 0; i2 < N; i2 ++){
            if(c_new(i2) > c_i){
              c_new(i2) = c_new(i2) - 1;
            }
          }
          
          //New label is the same as cj but now we have removed ci, so it could be c_j - 1!
          int c_c_new = c_new(j);
          
          //We merge data in cluster c_c_new after removing c_i-th component
          nj_new.shed_row(c_i);
          nj_new(c_c_new) = nj(c_i) + nj(c_j);
          
          // Propose unnormalised weight and latent variables
          S_m_new.shed_row(c_i);
          phi_star_new.shed_row(c_i);
          
          if(update_phij){//Consider changes only if udpate old cluster label, otherwise prob = 1
            arma::uvec setC_new = arma::find(c_new == c_c_new);
            double nj_m_new = setC_new.n_elem;
            
            S_m_new(c_c_new) = arma::randg( arma::distr_param( gamma_S + nj_m_new, 1.0 / (1.0 + u) ));
            log_ratio_SM -= log_dgamma_arma(S_m_new(c_c_new), gamma_S + nj_m_new, 1.0 + u);
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            
            //Compute mean vector for proposal: estimator of log-transition rates in the cluster
            arma::rowvec phi_prop_mean(p_tot,arma::fill::zeros);
            
            //Response processes
            for(int h = 0; h < pY; h ++){
              // Rcout << "h = " << h << "\n";
              
              int dY_h = dY[h];
              arma::mat trans_mat(dY_h,dY_h,arma::fill::zeros);
              arma::vec sojourn_times(dY_h,arma::fill::zeros);
              
              for(int ind = 0; ind < nj_m_new; ind++){
                
                int subj = setC_new(ind);
                std::vector<arma::vec> Y_i = Y[subj];
                arma::vec Y_i_h = Y_i[h];
                
                //Number of times
                int n_times_i_h = n_times(subj,h);
                
                //Time intervals
                std::vector<arma::vec> epsilon_i = epsilon[subj];
                arma::vec epsilon_i_h = epsilon_i[h];
                
                //Value of response at time 0
                double Yihj_minus_1 = Y_i_h(0);
                
                for(int jt = 1; jt < n_times_i_h; jt++){
                  
                  //Value of response at time j
                  double Yihj = Y_i_h(jt);
                  
                  if(Yihj_minus_1 != Yihj){//We do not consider transitions to the same state
                    trans_mat(Yihj_minus_1,Yihj) ++;
                  }
                  sojourn_times(Yihj_minus_1) += epsilon_i_h(jt-1);
                }
              }
              
              //Fill-in components of mean vector for proposal
              arma::uvec ind_dY = arma::regspace<arma::uvec>(n_rates_cum(h),n_rates_cum(h+1)-1);
              int count_Y = 0;
              for(int i2 = 0; i2< dY_h; i2 ++){
                for(int j2 = 0; j2 < dY_h; j2++){
                  if(i2 != j2){
                    if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                      phi_prop_mean(ind_dY(count_Y)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                    }else{
                      phi_prop_mean(ind_dY(count_Y)) = mu(ind_dY(count_Y));
                    }
                    count_Y++;
                  }
                }
              }
            }
            
            //Covariate processes
            for(int l = 0; l < pH; l ++){
              // Rcout << "l = " << l << "\n";
              
              int dH_l = dH[l];
              arma::mat trans_mat(dH_l,dH_l,arma::fill::zeros);
              arma::vec sojourn_times(dH_l,arma::fill::zeros);
              
              for(int ind = 0; ind < nj_m_new; ind++){
                
                int subj = setC_new(ind);
                std::vector<arma::vec> H_i = H[subj];
                arma::vec H_i_l = H_i[l];
                
                //Number of times
                int n_times_i_l = n_times(subj,pY+l);
                
                //Time intervals
                std::vector<arma::vec> epsilon_i = epsilon[subj];
                arma::vec epsilon_i_l = epsilon_i[pY+l];
                
                //Value of response at time 0
                double Hilj_minus_1 = H_i_l(0);
                
                for(int jt = 1; jt < n_times_i_l; jt++){
                  
                  //Value of response at time j
                  double Hilj = H_i_l(jt);
                  
                  if(Hilj_minus_1 != Hilj){//We do not consider transitions to the same state
                    trans_mat(Hilj_minus_1,Hilj) ++;
                  }
                  sojourn_times(Hilj_minus_1) += epsilon_i_l(jt-1);
                }
              }
              
              //Fill-in components of mean vector for proposal
              arma::uvec ind_dH = arma::regspace<arma::uvec>(n_rates_cum(pY+l),n_rates_cum(pY+l+1)-1);
              int count_H = 0;
              for(int i2 = 0; i2< dH_l; i2 ++){
                for(int j2 = 0; j2 < dH_l; j2++){
                  if(i2 != j2){
                    if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                      phi_prop_mean(ind_dH(count_H)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                    }else{
                      phi_prop_mean(ind_dH(count_H)) = mu(ind_dH(count_H));
                    }
                    count_H++;
                  }
                }
              }
            }
            
            // Rcout << "phi_prop_mean = " << phi_prop_mean << "\n";
            
            //Sample new phi_star_m using a Gaussian random-walk centred on the estimator in the cluster and with covariance:
            if(use_Sigma){//Sigma
              phi_star_new.row(c_c_new) = arma::mvnrnd(phi_prop_mean.t(), Sigma).t();
              log_ratio_SM -= log_dmvnrm_arma(phi_star_new.row(c_c_new), phi_prop_mean, Omega);
              // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            }else{//or identity with small diagonal values
              phi_star_new.row(c_c_new) = arma::mvnrnd(phi_prop_mean.t(), eye_mat_SM * s_SM).t();
              log_ratio_SM -= log_dmvnrm_arma(phi_star_new.row(c_c_new), phi_prop_mean, eye_mat_SM / s_SM);
              // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            }
          }
          
          //Numerator proposal is a reverse Split step
          //Find sets that are currently split
          arma::uvec setC_minusij = arma::find((c == c_i)||(c == c_j));
          setC_minusij.shed_rows(arma::find(setC_minusij == i)); //remove i-th element
          setC_minusij.shed_rows(arma::find(setC_minusij == j)); //remove j-th element
          int n_setC_minusij = setC_minusij.n_elem;
          
          // Contribution from allocation proposal
          for(int ind =  0; ind < n_setC_minusij; ind++){
            
            int i2 = setC_minusij(ind);
            
            arma::vec prob_2(2,arma::fill::zeros);
            prob_2(0) =  - lambda_dist_mat(i2,i) / 2.0;
            prob_2(1) =  - lambda_dist_mat(i2,j) / 2.0;
            
            prob_2 = exp(prob_2-max(prob_2));
            prob_2 = prob_2/sum(prob_2);
            
            if(c(i2) == c_i){
              log_ratio_SM += log(prob_2(0));
            }else{
              log_ratio_SM += log(prob_2(1));
            }
          }
          // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
          
          //Unnormalised weights and latent variables
          if(update_phij){
            for(int m = 0; m < 2; m++){
              int c_m = c_i * (1 - m) + c_j * m;
              
              arma::uvec setCm = arma::find(c == c_m);
              double nj_m = setCm.n_elem;
              
              log_ratio_SM += log_dgamma_arma(S_m(c_m), gamma_S + nj_m, 1.0 + u);
              // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
              
              //Compute mean vector for proposal: estimator of log-transition rates in the cluster
              arma::rowvec phi_prop_mean(p_tot,arma::fill::zeros);
              
              //Response processes
              for(int h = 0; h < pY; h ++){
                // Rcout << "h = " << h << "\n";
                
                int dY_h = dY[h];
                arma::mat trans_mat(dY_h,dY_h,arma::fill::zeros);
                arma::vec sojourn_times(dY_h,arma::fill::zeros);
                
                for(int ind = 0; ind < nj_m; ind++){
                  
                  int subj = setCm(ind);
                  std::vector<arma::vec> Y_i = Y[subj];
                  arma::vec Y_i_h = Y_i[h];
                  
                  //Number of times
                  int n_times_i_h = n_times(subj,h);
                  
                  //Time intervals
                  std::vector<arma::vec> epsilon_i = epsilon[subj];
                  arma::vec epsilon_i_h = epsilon_i[h];
                  
                  //Value of response at time 0
                  double Yihj_minus_1 = Y_i_h(0);
                  
                  for(int jt = 1; jt < n_times_i_h; jt++){
                    
                    //Value of response at time j
                    double Yihj = Y_i_h(jt);
                    
                    if(Yihj_minus_1 != Yihj){//We do not consider transitions to the same state
                      trans_mat(Yihj_minus_1,Yihj) ++;
                    }
                    sojourn_times(Yihj_minus_1) += epsilon_i_h(jt-1);
                  }
                }
                
                //Fill-in components fo mean vector for proposal
                arma::uvec ind_dY = arma::regspace<arma::uvec>(n_rates_cum(h),n_rates_cum(h+1)-1);
                int count_Y = 0;
                for(int i2 = 0; i2< dY_h; i2 ++){
                  for(int j2 = 0; j2 < dY_h; j2++){
                    if(i2 != j2){
                      if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                        phi_prop_mean(ind_dY(count_Y)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                      }else{
                        phi_prop_mean(ind_dY(count_Y)) = mu(ind_dY(count_Y));
                      }
                      count_Y++;
                    }
                  }
                }
              }
              
              //Covariate processes
              for(int l = 0; l < pH; l ++){
                // Rcout << "l = " << l << "\n";
                
                int dH_l = dH[l];
                arma::mat trans_mat(dH_l,dH_l,arma::fill::zeros);
                arma::vec sojourn_times(dH_l,arma::fill::zeros);
                
                for(int ind = 0; ind < nj_m; ind++){
                  
                  int subj = setCm(ind);
                  std::vector<arma::vec> H_i = H[subj];
                  arma::vec H_i_l = H_i[l];
                  
                  //Number of times
                  int n_times_i_l = n_times(subj,pY+l);
                  
                  //Time intervals
                  std::vector<arma::vec> epsilon_i = epsilon[subj];
                  arma::vec epsilon_i_l = epsilon_i[pY+l];
                  
                  //Value of response at time 0
                  double Hilj_minus_1 = H_i_l(0);
                  
                  for(int jt = 1; jt < n_times_i_l; jt++){
                    
                    //Value of response at time j
                    double Hilj = H_i_l(jt);
                    
                    if(Hilj_minus_1 != Hilj){//We do not consider transitions to the same state
                      trans_mat(Hilj_minus_1,Hilj) ++;
                    }
                    sojourn_times(Hilj_minus_1) += epsilon_i_l(jt-1);
                  }
                }
                
                //Fill-in components fo mean vector for proposal
                arma::uvec ind_dH = arma::regspace<arma::uvec>(n_rates_cum(pY+l),n_rates_cum(pY+l+1)-1);
                int count_H = 0;
                for(int i2 = 0; i2< dH_l; i2 ++){
                  for(int j2 = 0; j2 < dH_l; j2++){
                    if(i2 != j2){
                      if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                        phi_prop_mean(ind_dH(count_H)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                      }else{
                        phi_prop_mean(ind_dH(count_H)) = mu(ind_dH(count_H));
                      }
                      count_H++;
                    }
                  }
                }
              }
              
              // Rcout << "phi_star.row(c_m) = " << phi_star.row(c_m) << "\n";
              // Rcout << "phi_prop_mean = " << phi_prop_mean << "\n";
              
              //Sample new phi_star_m using a Gaussian random-walk centred on the estimator in the cluster and with covariance:
              if(use_Sigma){//Sigma
                log_ratio_SM += log_dmvnrm_arma(phi_star.row(c_m), phi_prop_mean, Omega);
                // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
              }else{//or identity with small diagonal values
                log_ratio_SM += log_dmvnrm_arma(phi_star.row(c_m), phi_prop_mean, eye_mat_SM / s_SM);
                // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
              }
            }
          }else{//Only consider "new" cluster c_i
            int c_m = c_i;
            
            arma::uvec setCm = arma::find(c == c_m);
            double nj_m = setCm.n_elem;
            
            log_ratio_SM += log_dgamma_arma(S_m(c_m), gamma_S + nj_m, 1.0 + u);
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            
            //Compute mean vector for proposal: estimator of log-transition rates in the cluster
            arma::rowvec phi_prop_mean(p_tot,arma::fill::zeros);
            
            //Response processes
            for(int h = 0; h < pY; h ++){
              // Rcout << "h = " << h << "\n";
              
              int dY_h = dY[h];
              arma::mat trans_mat(dY_h,dY_h,arma::fill::zeros);
              arma::vec sojourn_times(dY_h,arma::fill::zeros);
              
              for(int ind = 0; ind < nj_m; ind++){
                
                int subj = setCm(ind);
                std::vector<arma::vec> Y_i = Y[subj];
                arma::vec Y_i_h = Y_i[h];
                
                //Number of times
                int n_times_i_h = n_times(subj,h);
                
                //Time intervals
                std::vector<arma::vec> epsilon_i = epsilon[subj];
                arma::vec epsilon_i_h = epsilon_i[h];
                
                //Value of response at time 0
                double Yihj_minus_1 = Y_i_h(0);
                
                for(int jt = 1; jt < n_times_i_h; jt++){
                  
                  //Value of response at time j
                  double Yihj = Y_i_h(jt);
                  
                  if(Yihj_minus_1 != Yihj){//We do not consider transitions to the same state
                    trans_mat(Yihj_minus_1,Yihj) ++;
                  }
                  sojourn_times(Yihj_minus_1) += epsilon_i_h(jt-1);
                }
              }
              
              //Fill-in components fo mean vector for proposal
              arma::uvec ind_dY = arma::regspace<arma::uvec>(n_rates_cum(h),n_rates_cum(h+1)-1);
              int count_Y = 0;
              for(int i2 = 0; i2< dY_h; i2 ++){
                for(int j2 = 0; j2 < dY_h; j2++){
                  if(i2 != j2){
                    if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                      phi_prop_mean(ind_dY(count_Y)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                    }else{
                      phi_prop_mean(ind_dY(count_Y)) = mu(ind_dY(count_Y));
                    }
                    count_Y++;
                  }
                }
              }
            }
            
            //Covariate processes
            for(int l = 0; l < pH; l ++){
              // Rcout << "l = " << l << "\n";
              
              int dH_l = dH[l];
              arma::mat trans_mat(dH_l,dH_l,arma::fill::zeros);
              arma::vec sojourn_times(dH_l,arma::fill::zeros);
              
              for(int ind = 0; ind < nj_m; ind++){
                
                int subj = setCm(ind);
                std::vector<arma::vec> H_i = H[subj];
                arma::vec H_i_l = H_i[l];
                
                //Number of times
                int n_times_i_l = n_times(subj,pY+l);
                
                //Time intervals
                std::vector<arma::vec> epsilon_i = epsilon[subj];
                arma::vec epsilon_i_l = epsilon_i[pY+l];
                
                //Value of response at time 0
                double Hilj_minus_1 = H_i_l(0);
                
                for(int jt = 1; jt < n_times_i_l; jt++){
                  
                  //Value of response at time j
                  double Hilj = H_i_l(jt);
                  
                  if(Hilj_minus_1 != Hilj){//We do not consider transitions to the same state
                    trans_mat(Hilj_minus_1,Hilj) ++;
                  }
                  sojourn_times(Hilj_minus_1) += epsilon_i_l(jt-1);
                }
              }
              
              //Fill-in components fo mean vector for proposal
              arma::uvec ind_dH = arma::regspace<arma::uvec>(n_rates_cum(pY+l),n_rates_cum(pY+l+1)-1);
              int count_H = 0;
              for(int i2 = 0; i2< dH_l; i2 ++){
                for(int j2 = 0; j2 < dH_l; j2++){
                  if(i2 != j2){
                    if((trans_mat(i2,j2) > 0) & (sojourn_times(i2) > 0)){
                      phi_prop_mean(ind_dH(count_H)) = log(trans_mat(i2,j2)) - log(sojourn_times(i2));
                    }else{
                      phi_prop_mean(ind_dH(count_H)) = mu(ind_dH(count_H));
                    }
                    count_H++;
                  }
                }
              }
            }
            
            
            // Rcout << "phi_star.row(c_m) = " << phi_star.row(c_m) << "\n";
            // Rcout << "phi_prop_mean = " << phi_prop_mean << "\n";
            
            //Sample new phi_star_m using a Gaussian random-walk centred on the estimator in the cluster and with covariance:
            if(use_Sigma){//Sigma
              log_ratio_SM += log_dmvnrm_arma(phi_star.row(c_m), phi_prop_mean, Omega);
              // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            }else{//or identity with small diagonal values
              log_ratio_SM += log_dmvnrm_arma(phi_star.row(c_m), phi_prop_mean, eye_mat_SM / s_SM);
              // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            }
            
          }
          
          //Prior part extra component
          if(update_phij){
            log_ratio_SM += - u * (S_m_new(c_c_new) - S_m_i - S_m_j) + nj_new(c_c_new) * log(S_m_new(c_c_new)) - nj(c_i) * log(S_m_i) - nj(c_j) * log(S_m_j) + log(M-1.0) - log(Lambda);
            log_ratio_SM += log_dgamma_arma(S_m_new(c_c_new), gamma_S, 1.0) - log_dgamma_arma(S_m_i, gamma_S, 1.0) - log_dgamma_arma(S_m_j, gamma_S, 1.0);
            log_ratio_SM += log_dmvnrm_arma(phi_star_new.row(c_c_new), mu, Omega) - log_dmvnrm_arma(phi_star_i, mu, Omega) - log_dmvnrm_arma(phi_star_j, mu, Omega);
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
          }else{
            log_ratio_SM += - u * (- S_m_i) + nj(c_i) * (log(S_m_new(c_c_new)) - log(S_m_i)) + log(M-1.0) - log(Lambda);
            log_ratio_SM += - log_dgamma_arma(S_m_i, gamma_S, 1.0);
            log_ratio_SM += - log_dmvnrm_arma(phi_star_i, mu, Omega);
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
          }
          
          if(update_phij){
            //Loglike for new set C
            arma::uvec setC_ij_new = arma::find(c_new == c_c_new);
            log_ratio_SM += arma::accu(logLike_all(setC_ij_new, Y, H, phi_star_new.row(c_c_new), n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            
            //Loglike for old set Ci
            arma::uvec setC_i = arma::find(c == c_i);
            log_ratio_SM -= arma::accu(logLike_all(setC_i, Y, H, phi_star_i, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            
            //Loglike for old set Cj
            arma::uvec setC_j = arma::find(c == c_j);
            log_ratio_SM -= arma::accu(logLike_all(setC_j, Y, H, phi_star_j, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
          }else{
            //Loglike part
            arma::uvec setC_i = arma::find(c == c_i);
            log_ratio_SM += arma::accu(logLike_all(setC_i, Y, H, phi_star_new.row(c_c_new), n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
            log_ratio_SM -= arma::accu(logLike_all(setC_i, Y, H, phi_star_i, n_times, epsilon, X, Z, beta, gamma, n_rates_cum, dY, dH, impute_missing_Y, impute_missing_H));
            // Rcout << "log_ratio_SM = " << log_ratio_SM << "\n";
          }
          
        }//End Split/Merge calculations
        double accept_SM = 1.0;
        if( arma::is_finite(log_ratio_SM) ){
          if(log_ratio_SM < 0){
            accept_SM = exp(log_ratio_SM);
          }
        }else{
          accept_SM = 0.0;
        }
        
        SM_accept += accept_SM;
        SM_count ++;
        
        if( arma::randu<double>() < accept_SM ){
          // Rcout << "accept!\n";
          
          if(is_Split){
            SM_nSplit ++;
          }else{
            SM_nMerge ++;
          }
          
          c = c_new;
          nj = nj_new;
          S_m = S_m_new;
          phi_star = phi_star_new;
        }
        c_star = arma::unique(c); //This returns a vector of unique elements in ascending order!
        K_N = c_star.n_elem;
        M = K_N + M_na;
      }
    }
    
    
    
    //Save output for this iteration
    if( (it + 1 > (n_burn1 + n_burn2)) & (((it + 1 - n_burn1 - n_burn2) / thin - floor((it + 1 - n_burn1 - n_burn2) / thin)) == 0 )){
      
      int iter = (it + 1 - n_burn1 - n_burn2)/thin - 1;
      
      // Lists
      phi_star_out[iter] = phi_star;
      S_m_out[iter] = S_m;
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
      Omega_out.slice(iter) = Omega;
      G_out.slice(iter) = G;
      G0_out.slice(iter) = G0;
      
      // Mats
      mu_out.row(iter) = mu;
      m_mu_out.row(iter) = m_mu;
      c_out.row(iter) = c.t();
      
      // Vecs
      d_out(iter) = d;
      lambda_out(iter) = lambda;
      mu_d_out(iter) = mu_d;
      size_G0_out(iter) = size_G0;
      sum_weights_out(iter) = sum_weights;
      u_out(iter) = u;
      K_N_out(iter) = K_N;
      M_out(iter) = M;
      gamma_S_out(iter) = gamma_S;
      Lambda_out(iter) = Lambda;
      k0_out(iter) = k0;
      
      //Entropy
      for(int j = 0; j < K_N; j++){
        Entropy_out(iter) -= nj(j) / N * log(nj(j) / N);
      }
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
  if(update_gamma_S){
    Rcout << "gamma_S a.r. = " << gamma_S_accept/gamma_S_count << "\n";
  }
  if(update_d){
    if(!d_beta_prior){
      Rcout << "d a.r. = " << d_accept/d_count << "\n";
      Rcout << "lambda a.r. = " << lambda_accept/lambda_count << "\n";
    }
  }
  if(SM_alg){
    Rcout << "SM a.r. = " << SM_accept/SM_count << "\n";
    Rcout << "SM_nSplit = " << SM_nSplit << "\n";
    Rcout << "SM_nMerge = " << SM_nMerge << "\n";
  }
  Rcout << "n_nan_ci = " << n_nan_ci << "\n";
  Rcout << "n_inf_ci = " << n_inf_ci << "\n";
  
  
  List Binder_List;
  if(update_c){
    //Compute Binder estimate with equal costs and return it
    pij.zeros();
    for(int it = 0; it < n_save; it++){
      c_out_row = c_out.row(it);
      for(int i = 0; i < (N-1); i++){
        for(int j = (i+1); j < N; j++){
          pij(i,j) += (c_out_row(i) == c_out_row(j));
        }
      }
    }
    pij = pij/n_save; //Diagonal is zero but it is not a problem
    
    Binder_f.zeros();
    for(int it = 0; it < n_save; it++){
      c_out_row = c_out.row(it);
      for(int i = 0; i < (N-1); i++){
        for(int j = (i+1); j < N; j++){
          Binder_f(it) += (pij(i,j) - Binder_const) * (c_out_row(i) == c_out_row(j));
        }
      }
    }
    Binder_index = arma::find(Binder_f == max(Binder_f),1);
    Binder_est = c_out.row(arma::conv_to<int>::from(Binder_index));
    
    Binder_List = List::create(Named("pij") = pij, Named("Binder_f") = Binder_f, Named("Binder_est") = Binder_est, Named("Entropy_out") = Entropy_out);
  }else{
    Binder_List = List::create(Named("Entropy_out") = Entropy_out);
  }
  
  List Graph_List = List::create(Named("Omega_out") = Omega_out, Named("G_out") = G_out, Named("G0_out") = G0_out, Named("sizeG0_out") = size_G0_out, Named("sum_weights_out") = sum_weights_out, Named("d_out") = d_out, Named("mu_d_out") = mu_d_out, Named("lambda_out") = lambda_out);
  
  
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
  
  List SM_List = List::create(Named("SM_nSplit") = SM_nSplit, Named("SM_nMerge") = SM_nMerge, Named("SM_accept") = SM_accept, Named("SM_count") = SM_count);
  
  if(impute_missing_Y){
    if(impute_missing_H){
      return List::create(Named("Graph_List") = Graph_List, Named("mu_out") = mu_out, Named("m_mu_out") = m_mu_out, Named("k0_out") = k0_out, Named("phi_star_out") = phi_star_out, Named("S_m_out") = S_m_out, Named("c_out") = c_out, Named("u_out") = u_out, Named("K_N_out") = K_N_out, Named("M_out") = M_out, Named("gamma_S_out") = gamma_S_out, Named("Lambda_out") = Lambda_out, Named("XZ_List") = XZ_List, Named("SM_List") = SM_List, Named("Binder_List") = Binder_List, Named("Y_out") = Y_out, Named("H_out") = H_out);
    }else{
      return List::create(Named("Graph_List") = Graph_List, Named("mu_out") = mu_out, Named("m_mu_out") = m_mu_out, Named("k0_out") = k0_out, Named("phi_star_out") = phi_star_out, Named("S_m_out") = S_m_out, Named("c_out") = c_out, Named("u_out") = u_out, Named("K_N_out") = K_N_out, Named("M_out") = M_out, Named("gamma_S_out") = gamma_S_out, Named("Lambda_out") = Lambda_out, Named("XZ_List") = XZ_List, Named("SM_List") = SM_List, Named("Binder_List") = Binder_List, Named("Y_out") = Y_out);
    }
  }else{
    if(impute_missing_H){
      return List::create(Named("Graph_List") = Graph_List, Named("mu_out") = mu_out, Named("m_mu_out") = m_mu_out, Named("k0_out") = k0_out, Named("phi_star_out") = phi_star_out, Named("S_m_out") = S_m_out, Named("c_out") = c_out, Named("u_out") = u_out, Named("K_N_out") = K_N_out, Named("M_out") = M_out, Named("gamma_S_out") = gamma_S_out, Named("Lambda_out") = Lambda_out, Named("XZ_List") = XZ_List, Named("SM_List") = SM_List, Named("Binder_List") = Binder_List, Named("H_out") = H_out);
    }else{
      return List::create(Named("Graph_List") = Graph_List, Named("mu_out") = mu_out, Named("m_mu_out") = m_mu_out, Named("k0_out") = k0_out, Named("phi_star_out") = phi_star_out, Named("S_m_out") = S_m_out, Named("c_out") = c_out, Named("u_out") = u_out, Named("K_N_out") = K_N_out, Named("M_out") = M_out, Named("gamma_S_out") = gamma_S_out, Named("Lambda_out") = Lambda_out, Named("XZ_List") = XZ_List, Named("SM_List") = SM_List, Named("Binder_List") = Binder_List);
    }
  }
}

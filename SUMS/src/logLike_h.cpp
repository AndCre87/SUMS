//Likelihood contributions for updating beta_h and gamma_h (coefficients of covariates for h-th process)
//In this case we need the likelihood of each data point given the ci-th phi_star
//This means evaluating the log-likelihood of each subject and its associated phi_star value only for the h-th process involving covariates (Y^h)
//Than returning the sum of these log-likelihoods

// [[Rcpp::depends(RcppArmadillo)]]
// #define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
// #include <omp.h>
#include "logLike_h.h"
#include "expmat_Tesch1.h"

double logLike_h(arma::uvec setC, std::vector<std::vector<arma::vec>> Y, arma::vec c, arma::mat phi_star, arma::vec n_times_h, std::vector<std::vector<arma::vec>> epsilon, arma::mat X_h, std::vector<arma::mat> Z_h, arma::mat beta_h, arma::mat gamma_h, arma::vec n_rates_cum, arma::vec dY, int ind_h, bool impute_missing_Y){
  
  int n_ind = setC.n_elem;
  arma::vec loglike(n_ind,arma::fill::zeros);
  
  arma::uvec ind_dY = arma::regspace<arma::uvec>(n_rates_cum(ind_h),n_rates_cum(ind_h+1)-1);
  int dY_h = dY[ind_h]; //Number of states
  
  // #pragma omp parallel for
  for(int ind = 0; ind < n_ind; ind++){
    
    int i = setC(ind);
    
    arma::mat Z_h_i = Z_h[i];
    
    std::vector<arma::vec> epsilon_i = epsilon[i];
    arma::vec epsilon_i_h = epsilon_i[ind_h];
    
    int n_times_h_i = n_times_h(i);
    
    int m = c(i); //Component for i-th data point
    
    arma::rowvec phi_star_m = phi_star.row(m);
    
    arma::rowvec phi_star_tilde_m = phi_star_m(ind_dY).t();
    
    //Loglike from response process Y
    std::vector<arma::vec> Y_i = Y[i];
    arma::vec Y_i_h = Y_i[ind_h];
    
    if(impute_missing_Y){
      int j = 0;
      //We assume time-varying continuous covariate at time zero is observed
      arma::rowvec lambda_Yihj = exp( phi_star_tilde_m + X_h.row(i) * beta_h + Z_h_i.row(j) * gamma_h );
      
      //Value of response at time j
      double Yihj = Y_i_h(j);
      
      //Known expressions of stationary distribution
      //Solutions obtained from Matlab symbolic calculations
      if(dY_h == 2){
        loglike(ind) += log(lambda_Yihj(1 - Yihj))  - log(arma::sum(lambda_Yihj));
      }
      if(dY_h == 3){
        //Easier to populate the generator matrix first (no need for diagonal)
        arma::mat Q_Y(dY_h, dY_h,arma::fill::zeros);
        int count_dY = 0;
        for(int i2 = 0; i2< dY_h; i2 ++){
          for(int j2 = 0; j2 < dY_h; j2++){
            if(i2 != j2){
              Q_Y(i2,j2) = lambda_Yihj(count_dY);
              count_dY++;
            }
          }
        }
        
        double denom = (Q_Y(0,2)*Q_Y(1,0) + Q_Y(0,1)*Q_Y(1,2) + Q_Y(0,2)*Q_Y(1,2) + Q_Y(0,1)*Q_Y(2,0) + Q_Y(0,1)*Q_Y(2,1) + Q_Y(0,2)*Q_Y(2,1) + Q_Y(1,0)*Q_Y(2,0) + Q_Y(1,0)*Q_Y(2,1) + Q_Y(1,2)*Q_Y(2,0));
        
        if(Yihj == 0){
          loglike(ind) += log(Q_Y(1,0)*Q_Y(2,0) + Q_Y(1,0)*Q_Y(2,1) + Q_Y(1,2)*Q_Y(2,0)) - log(denom);
        }
        if(Yihj == 1){
          loglike(ind) += log(Q_Y(0,1)*Q_Y(2,0) + Q_Y(0,1)*Q_Y(2,1) + Q_Y(0,2)*Q_Y(2,1)) - log(denom);
        }
        if(Yihj == 2){
          loglike(ind) += log(Q_Y(0,2)*Q_Y(1,0) + Q_Y(0,1)*Q_Y(1,2) + Q_Y(0,2)*Q_Y(1,2)) - log(denom);
        }
        
      }
    }
    
    double Yihj_minus_1 = Y_i_h(0);
    
    for(int j = 1; j < n_times_h_i; j++){
      
      arma::rowvec lambda_Yihj = exp( phi_star_tilde_m + X_h.row(i) * beta_h + Z_h_i.row(j) * gamma_h );
      
      double eps_ij = epsilon_i_h(j-1);
      
      //Value of response at time j
      double Yihj = Y_i_h(j);
      
      if(dY_h > 2){
        arma::mat Q_matrix(dY_h,dY_h,arma::fill::zeros);
        int count_dY = 0;
        for(int i2 = 0; i2< dY_h; i2 ++){
          for(int j2 = 0; j2 < dY_h; j2++){
            if(i2 != j2){
              Q_matrix(i2,j2) = lambda_Yihj(count_dY);
              count_dY++;
            }
          }
          Q_matrix(i2,i2) = - arma::accu(Q_matrix.row(i2));
        } 
        
        // arma::mat P_matrix = arma::expmat(eps_ij * Q_matrix);
        arma::mat P_matrix = expmat_Tesch1(eps_ij * Q_matrix);
        
        //Select entry in the Pmatrix
        loglike(ind) += log(P_matrix(Yihj_minus_1,Yihj));
        
      }else{
        double lambda_sum_Y = arma::sum(lambda_Yihj);

        double prs_Y = lambda_Yihj(Yihj_minus_1) * ( 1 - exp( - (lambda_sum_Y * eps_ij) ) );

        loglike(ind) += (Yihj == Yihj_minus_1) * (log(lambda_sum_Y - prs_Y)) + (Yihj != Yihj_minus_1) * log(prs_Y) - log(lambda_sum_Y);
      }
      
      Yihj_minus_1 = Yihj;
    }
  }
  
  return arma::accu(loglike);
}
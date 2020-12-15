//Impute missing values at time 0 for a model where those states are distributed as their stationary distribution

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <tuple>
#include "expmat_Tesch1.h"
#include "impute_missing_dH0.h"

std::vector<std::vector<arma::vec>> impute_missing_dH0(std::vector<std::vector<arma::vec>> Y, arma::mat is_NA_Y, arma::vec c, arma::mat phi_star, arma::mat n_times, std::vector<std::vector<arma::vec>> epsilon, std::vector<arma::mat> X, std::vector<std::vector<arma::mat>> Z, arma::field<arma::mat> beta, arma::field<arma::mat> gamma, arma::vec n_rates_cum, arma::vec dY, bool impute_missing){
  
  int N = Y.size();
  
  std::vector<arma::vec> Y_i = Y[0];
  int pY = Y_i.size();
  
  //Response processes
  if(impute_missing){
    for(int i = 0; i < N; i++){
      std::vector<arma::vec> Y_i = Y[i];
      std::vector<arma::vec> epsilon_i = epsilon[i];
      
      int m = c(i);
      arma::rowvec phi_m = phi_star.row(m);
      
      for(int h = 0; h < pY; h ++){
        arma::vec Y_i_h = Y_i[h];
        
        if(is_NA_Y(i,h) == 1){//Missing time zero for this subject and this process
          
          //Time intervals
          arma::vec epsilon_i_h = epsilon_i[h];
          
          //Constant covariates
          arma::mat X_h = X[h];
          arma::mat beta_h = beta[h];
          
          //Time-varying covariates
          std::vector<arma::mat> Z_h = Z[h];
          arma::mat Z_h_i = Z_h[i];
          arma::mat gamma_h = gamma[h];
          
          //Process rates
          arma::uvec ind_dY = arma::regspace<arma::uvec>(n_rates_cum(h),n_rates_cum(h+1)-1);
          arma::rowvec phi_tilde_m = phi_m(ind_dY).t();
          
          //Sampling probabilities
          int dY_h = dY[h];
          arma::vec sampling_prob(dY_h,arma::fill::zeros);
          
          ///////////////////////////
          //Contribution from time 0
          int j = 0;
          
          //We assume time-varying continuous covariate at time zero is observed
          arma::rowvec lambda_Yihj = exp(phi_tilde_m + X_h.row(i) * beta_h + Z_h_i.row(j) * gamma_h);
          
          //Known expressions of stationary distribution
          //Solutions obtained from Matlab symbolic calculations
          if(dY_h == 2){
            sampling_prob(0) += log(lambda_Yihj(1))  - log(arma::sum(lambda_Yihj));
            sampling_prob(1) += log(lambda_Yihj(0))  - log(arma::sum(lambda_Yihj));
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
            
            //We will normalise again so no need for the denominator
            // double denom = (Q_Y(0,2)*Q_Y(1,0) + Q_Y(0,1)*Q_Y(1,2) + Q_Y(0,2)*Q_Y(1,2) + Q_Y(0,1)*Q_Y(2,0) + Q_Y(0,1)*Q_Y(2,1) + Q_Y(0,2)*Q_Y(2,1) + Q_Y(1,0)*Q_Y(2,0) + Q_Y(1,0)*Q_Y(2,1) + Q_Y(1,2)*Q_Y(2,0));
            
            sampling_prob(0) += log(Q_Y(1,0)*Q_Y(2,0) + Q_Y(1,0)*Q_Y(2,1) + Q_Y(1,2)*Q_Y(2,0)); // - log(denom);
            sampling_prob(1) += log(Q_Y(0,1)*Q_Y(2,0) + Q_Y(0,1)*Q_Y(2,1) + Q_Y(0,2)*Q_Y(2,1)); // - log(denom);
            sampling_prob(2) += log(Q_Y(0,2)*Q_Y(1,0) + Q_Y(0,1)*Q_Y(1,2) + Q_Y(0,2)*Q_Y(1,2)); // - log(denom);
          }
          
          /////////////////////////////////////////////////////
          //Contribution from time 1 (transition from t0 to t1)
          j = 1;
          
          double Yihj = Y_i_h(j);
          
          //We assume time-varying continuous covariate at time zero is observed
          lambda_Yihj = exp(phi_tilde_m + X_h.row(i) * beta_h + Z_h_i.row(j) * gamma_h);
          
          double eps_ij = epsilon_i_h(j-1);
          
          arma::mat Q_Y(dY_h, dY_h,arma::fill::zeros);
          int count_dY = 0;
          for(int i2 = 0; i2< dY_h; i2 ++){
            for(int j2 = 0; j2 < dY_h; j2++){
              if(i2 != j2){
                Q_Y(i2,j2) = lambda_Yihj(count_dY);
                count_dY++;
              }
            }
            Q_Y(i2,i2) = - arma::accu(Q_Y.row(i2));
          } 
          arma::mat P_matrix = expmat_Tesch1(eps_ij * Q_Y);
          
          sampling_prob += log(P_matrix.col(Yihj));
          
          //Fix possible NA's due to numerical overflow
          if(sampling_prob.has_nan()){
            sampling_prob.replace(arma::datum::nan, - arma::datum::inf);
          }
          
          //Normalise and sample
          sampling_prob = exp(sampling_prob-max(sampling_prob));
          sampling_prob = sampling_prob/arma::sum(sampling_prob);
          
          double aux_runif = arma::randu();
          arma::vec sampling_prob_cum = arma::cumsum(sampling_prob);
          arma::uvec hh_vec = arma::find(sampling_prob_cum >= aux_runif, 1, "first");
          Y_i_h(0) = hh_vec(0);
        }
        Y_i[h] = Y_i_h;
      }
      Y[i] = Y_i;
    }
  }
  
  return Y;
}
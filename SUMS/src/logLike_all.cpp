//Likelihood contribution for given values of phi_star given one or more data points
//Here we compute the likelihood of a set of subjects (maybe one or all of them) for EACH given value of phi_star (could be only one)
//For instance, in the parametric case (or in a singleton), we compute the loglike for one subject for each row of phi_star
//This involves all the processes, not only the ones with covariates!

// [[Rcpp::depends(RcppArmadillo)]]
// #define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
// #include <omp.h>
#include "logLike_all.h"
#include "expmat_Tesch1.h"

arma::vec logLike_all(arma::uvec setC, std::vector<std::vector<arma::vec>> Y, std::vector<std::vector<arma::vec>> H, arma::mat phi_star, arma::mat n_times, std::vector<std::vector<arma::vec>> epsilon, std::vector<arma::mat> X, std::vector<std::vector<arma::mat>> Z, arma::field<arma::mat> beta, arma::field<arma::mat> gamma, arma::vec n_rates_cum, arma::vec dY, arma::vec dH, bool impute_missing_Y, bool impute_missing_H){
  
  int n_ind = setC.n_elem, M = phi_star.n_rows;
  
  std::vector<arma::vec> Y_i = Y[0], H_i = H[0];
  int pY = Y_i.size(), pH = H_i.size();
  
  arma::vec loglike_vec(M,arma::fill::zeros);
  
  // #pragma omp parallel for
  for(int m = 0; m < M; m++){
    
    double loglike = 0.0;
    
    arma::rowvec phi_star_m = phi_star.row(m);
    
    //Loglike from response processes Y's
    for(int h = 0; h < pY; h ++){
      
      arma::uvec ind_dY = arma::regspace<arma::uvec>(n_rates_cum(h),n_rates_cum(h+1)-1);
      
      arma::rowvec phi_star_tilde_m = phi_star_m(ind_dY).t();
      
      std::vector<arma::mat> Z_h = Z[h];
      arma::mat X_h = X[h];
      
      arma::mat beta_h = beta[h];
      arma::mat gamma_h = gamma[h];
      int dY_h = dY[h];
      
      for(int ind = 0; ind < n_ind; ind++){
        
        int i = setC(ind);
        
        arma::mat Z_h_i = Z_h[i];
        
        std::vector<arma::vec> epsilon_i = epsilon[i];
        arma::vec epsilon_i_h = epsilon_i[h];
        
        int n_times_i_h = n_times(i,h);
        
        std::vector<arma::vec> Y_i = Y[i];
        arma::vec Y_i_h = Y_i[h];
        
        if(impute_missing_Y){
          int j = 0;
          //We assume time-varying continuous covariate at time zero is observed
          arma::rowvec lambda_Yihj = exp( phi_star_tilde_m + X_h.row(i) * beta_h + Z_h_i.row(j) * gamma_h );
          
          //Value of response at time j
          double Yihj = Y_i_h(j);
          
          //Known expressions of stationary distribution
          //Solutions obtained from Matlab symbolic calculations
          if(dY_h == 2){
            loglike += log(lambda_Yihj(1 - Yihj))  - log(arma::sum(lambda_Yihj));
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
              loglike += log(Q_Y(1,0)*Q_Y(2,0) + Q_Y(1,0)*Q_Y(2,1) + Q_Y(1,2)*Q_Y(2,0)) - log(denom);
            }
            if(Yihj == 1){
              loglike += log(Q_Y(0,1)*Q_Y(2,0) + Q_Y(0,1)*Q_Y(2,1) + Q_Y(0,2)*Q_Y(2,1)) - log(denom);
            }
            if(Yihj == 2){
              loglike += log(Q_Y(0,2)*Q_Y(1,0) + Q_Y(0,1)*Q_Y(1,2) + Q_Y(0,2)*Q_Y(1,2)) - log(denom);
            }
            
          }
        }
        
        double Yihj_minus_1 = Y_i_h(0);
        for(int j = 1; j < n_times_i_h; j++){
          
          arma::rowvec lambda_Yihj = exp( phi_star_tilde_m + X_h.row(i) * beta_h + Z_h_i.row(j) * gamma_h );
          
          double eps_ij = epsilon_i_h(j-1);
          
          //Value of response at time j
          double Yihj = Y_i_h(j);
          
          if(dY_h > 2){
            arma::mat Q_matrix(dY_h, dY_h,arma::fill::zeros);
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
            loglike += log(P_matrix(Yihj_minus_1,Yihj));
            
          }else{
            double lambda_sum_Y = arma::sum(lambda_Yihj);
            
            double prs_Y = lambda_Yihj(Yihj_minus_1) * ( 1 - exp( - (lambda_sum_Y * eps_ij) ) );
            
            loglike += (Yihj == Yihj_minus_1) * (log(lambda_sum_Y - prs_Y)) + (Yihj != Yihj_minus_1) * log(prs_Y) - log(lambda_sum_Y);
          }
          
          Yihj_minus_1 = Yihj;
        }
      }
    }
    
    
    //Loglike from symptom processes H's (convenient to reverse the for loops because we do not have beta and gamma lists)
    for(int ind = 0; ind < n_ind; ind++){
      
      int i = setC(ind);
      
      std::vector<arma::vec> epsilon_i = epsilon[i];
      
      std::vector<arma::vec> H_i = H[i];
      for(int l = 0; l < pH; l ++){
        
        int n_times_i_l = n_times(i,pY+l);
        
        arma::vec epsilon_i_l = epsilon_i[pY+l];
        
        arma::uvec ind_dH = arma::regspace<arma::uvec>(n_rates_cum(pY+l),n_rates_cum(pY+l+1)-1);
        int dH_l = dH[l];
        
        arma::rowvec lambda_Hl = exp(phi_star_m(ind_dH).t());
        
        arma::vec H_i_l = H_i[l];
        
        if(impute_missing_H){
          int j = 0;
          
          //Value of response at time j
          double Hilj = H_i_l(j);
          
          //Known expressions of stationary distribution
          //Solutions obtained from Matlab symbolic calculations
          if(dH_l == 2){
            loglike += log(lambda_Hl(1 - Hilj))  - log(arma::sum(lambda_Hl));
          }
          if(dH_l == 3){
            //Easier to populate the generator matrix first (no need for diagonal)
            arma::mat Q_H(dH_l, dH_l,arma::fill::zeros);
            int count_dH = 0;
            for(int i2 = 0; i2< dH_l; i2 ++){
              for(int j2 = 0; j2 < dH_l; j2++){
                if(i2 != j2){
                  Q_H(i2,j2) = lambda_Hl(count_dH);
                  count_dH++;
                }
              }
            }
            
            double denom = (Q_H(0,2)*Q_H(1,0) + Q_H(0,1)*Q_H(1,2) + Q_H(0,2)*Q_H(1,2) + Q_H(0,1)*Q_H(2,0) + Q_H(0,1)*Q_H(2,1) + Q_H(0,2)*Q_H(2,1) + Q_H(1,0)*Q_H(2,0) + Q_H(1,0)*Q_H(2,1) + Q_H(1,2)*Q_H(2,0));
            
            if(Hilj == 0){
              loglike += log(Q_H(1,0)*Q_H(2,0) + Q_H(1,0)*Q_H(2,1) + Q_H(1,2)*Q_H(2,0)) - log(denom);
            }
            if(Hilj == 1){
              loglike += log(Q_H(0,1)*Q_H(2,0) + Q_H(0,1)*Q_H(2,1) + Q_H(0,2)*Q_H(2,1)) - log(denom);
            }
            if(Hilj == 2){
              loglike += log(Q_H(0,2)*Q_H(1,0) + Q_H(0,1)*Q_H(1,2) + Q_H(0,2)*Q_H(1,2)) - log(denom);
            }
            
          }
        }
        
        double Hilj_minus_1 = H_i_l(0);
        
        for(int j = 1; j < n_times_i_l; j++){
          
          double eps_ij = epsilon_i_l(j-1);
          
          //Value of response at time j
          double Hilj = H_i_l(j);
          
          if(dH_l > 2){
            
            arma::mat Q_matrix(dH_l, dH_l, arma::fill::zeros);
            int count_dH = 0;
            for(int i2 = 0; i2< dH_l; i2 ++){
              for(int j2 = 0; j2 < dH_l; j2++){
                if(i2 != j2){
                  Q_matrix(i2,j2) = lambda_Hl(count_dH);
                  count_dH++;
                }
              }
              Q_matrix(i2,i2) = - arma::accu(Q_matrix.row(i2));
            } 
            
            arma::mat P_matrix = expmat_Tesch1(eps_ij * Q_matrix);
            
            //Select entry in the Pmatrix
            loglike += log(P_matrix(Hilj_minus_1,Hilj));
            
          }else{
            double lambda_sum_H = arma::sum(lambda_Hl);
            
            double prs_H = lambda_Hl(Hilj_minus_1) * ( 1 - exp( - (lambda_sum_H * eps_ij) ) );
            
            loglike += (Hilj == Hilj_minus_1) * (log(lambda_sum_H - prs_H)) + (Hilj != Hilj_minus_1) * log(prs_H) - log(lambda_sum_H);
          }
          Hilj_minus_1 = Hilj;
        }
      }
    }
    loglike_vec(m) = loglike;
  }
  
  return loglike_vec;
}
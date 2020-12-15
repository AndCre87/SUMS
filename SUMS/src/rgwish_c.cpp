// G-Wishart sampler based on Lenkoski and Dobra (2013)

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "rgwish_c.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// G is adjacency matrix which has zero in its diagonal // threshold = 1e-8 but can be changed
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

// [[Rcpp::export]]
arma::mat rgwish_c( double nu, arma::mat Ti, arma::mat G, double threshold ){
  
  //Recall that Ti = chol(inv(Psi))

  int p = G.n_cols;	
  
  arma::mat Omega = arma::wishrnd( Ti, nu + p - 1, Ti ), sigma_start = arma::inv_sympd(Omega), sigma = sigma_start;

  double mean_diff = 1.0;
  while( mean_diff > threshold ){

    arma::mat sigma_last = sigma;

    for(int i = 0; i < p; i++){

      // Count size of node
      int size_node = arma::accu(G.row(i));

      if( size_node > 0 ){
        arma::uvec N_i(size_node);
        arma::vec beta_star(p,arma::fill::zeros), sigma_start_N_i(size_node, arma::fill::zeros);
        int l = 0;
        for(int j = 0; j < p; j++ ){
          if( G(i,j) == 1 ){
            sigma_start_N_i(l) = sigma_start(i,j);
            N_i(l) = j;
            l++;
          }else{
            beta_star(j) = 0.0; 
          }
        }
        
        arma::mat sigma_N_i = sigma.submat(N_i, N_i);

        sigma_start_N_i = arma::inv_sympd(sigma_N_i) * sigma_start_N_i;

        for(int j = 0; j < size_node; j++ ){
          beta_star(N_i(j)) = sigma_start_N_i(j);
        }

        arma::vec sigma_start_i = sigma * beta_star;
        
        //Fill matrix by skipping position i
        for(int j = 0; j < i; j++ ){
          sigma(i,j) = sigma_start_i(j);
          sigma(j,i) = sigma_start_i(j);
        }
        
        for(int j = i + 1; j < p; j++ ){
          sigma(i,j) = sigma_start_i(j);
          sigma(j,i) = sigma_start_i(j);
        }
      }else{
        //Fill matrix by skipping position i
        for(int j = 0; j < i; j++ ){
          sigma(i,j) = 0.0;
          sigma(j,i) = 0.0;
        }
        
        for(int j = i + 1; j < p; j++ ){
          sigma(i,j) = 0.0;
          sigma(j,i) = 0.0;
        }
      }
    }

    arma::mat abs_diff_mat = arma::abs( sigma - sigma_last );
    mean_diff =  arma::accu(abs_diff_mat) / pow(p,2); //In L1 entrywise matrix norm (mean)
    // mean_diff =  arma::accu(abs_diff_mat); //In L1 entrywise matrix norm (theory?)
  }

  Omega = arma::inv_sympd(sigma);
  return Omega;
}


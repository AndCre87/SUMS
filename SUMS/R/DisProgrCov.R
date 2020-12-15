## Wrapper to call parametric mixture model for estimating Covariance matrix (mostly a test function)

DisProgrCov = function( data, Alg_list = Alg_list, Param_list = NULL ){
  
  if( thin < 1 ) stop( "Number of thinned iterations must be positive" )

  phi <- data$phi
  if( any( is.na( phi ) ) ) stop( "There are some missing values and I cannot proceed!!!" )	

  n_rates <- data$n_rates
  n_states <- data$n_states
  dY <- n_states[[1]]
  if( min(dY) < 1 ) stop( "Number of states in response processes must be positive" )
  dH <- n_states[[2]]
  if( min(dH) < 1 ) stop( "Number of states in symptom processes must be positive" )
  
  if(is.null(n_rates)) stop( "Provide number of states for each Markov process involved" )
  
  p0 <- (length(n_rates)-1)
  if( p0 < 3 ) stop( "Number of variables/nodes ('p0') must be more than 2" )

  if(is.null(Param_list)){
    stop("Provide some form of hyperparameter specification!")
  }
  
  # - -  main BDMCMC algorithms implemented in C++ - - - - - - - - - - - - - - - - - - - - - - - |
  MCMC_input <- list(data = phi, Alg_list = Alg_list, n_rates = n_rates, Param_list = Param_list)
  
  MCMC_output = DisProgrCov_Gibbs(MCMC_input)
  
  return( MCMC_output )   
}


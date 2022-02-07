## Wrapper to call the mixture model with some variations (parametric, finite mixture with random M, DP)

DisProgr = function( data, model_list, Alg_list, Param_list = NULL ){
  
  #Extract data from list (names are important!)
  Y_obs <- data$Y
  H_obs <- data$H
  n_times <- data$n_times
  epsilon <- data$epsilon
  
  #Is there any NA at time zero? (from outside wrapper we only leave these type of NA)
  N <- length(Y_obs)
  pY <- length(Y_obs[[1]])
  is_NA_Y <- matrix(0, N, pY)
  for(i in 1:N){
    for(h in 1:pY){
      if(is.na(Y_obs[[i]][[h]][1])){
        is_NA_Y[i,h] <- 1
      }
    }
  }
  
  is_NA_H <- NULL
  if(!is.null(H_obs)){
    pH <- length(H_obs[[1]])
    is_NA_H <- matrix(0, N, pH)
    for(i in 1:N){
      for(l in 1:pH){
        if(is.na(H_obs[[i]][[l]][1])){
          is_NA_H[i,l] <- 1
        }
      }
    }
  }
  
  
  #Fixed covariates might be absent
  X <- data$X
  is_null_X <- FALSE
  if(is.null(X)){
    is_null_X <- TRUE
    
    #Initialise to one covariate per process equal to zeros
    g <- rep(1,pY)
    X <- list()
    for(h in 1:pY){
      X[[h]] <- matrix(0, N, g[h])
    }
  }
  
  #Time-varying covariates might be absent
  Z <- data$Z
  is_null_Z <- FALSE
  if(is.null(Z)){
    is_null_Z <- TRUE
    
    #Initialise to one covariate per process equal to zeros
    q <- 1
    Z <- vector("list", length = pY)
    for(h in 1:pY){
      Z_h <- vector("list", length = N)
      for(i in 1:N){
        n_times_i <- n_times[i,h]
        times_i <- times[[i]][[h]]
        Z_h[[i]] <- matrix(0, n_times_i,1)
      }
      Z[[h]] <- Z_h
    }
  }
  
  n_rates <- data$n_rates
  n_states <- data$n_states
  dY <- n_states[[1]]
  if( min(dY) < 1 ) stop( "Number of states in response processes must be positive" )
  dH <- n_states[[2]]
  
  if( thin < 1 ) stop( "Number of thinned iteration must be positive" )
  
  if( any( is.na( data ) ) ) stop( "There are some missing values and I cannot proceed!!!" )	
  
  if(is.null(n_rates)) stop( "Provide number of states for each Markov process involved" )
  
  p0 <- (length(n_rates)-1)
  if( p0 < 2 ) print( "Number of variables/nodes ('p0') is less than 2. No graph will be inferred." )
  
  if(is.null(Param_list)){
    stop("Provide some form of hyperparameter specification!")
  }
  
  c_init <- model_list$c_init
  
  #Even if H and dH are NULL you can still pass them
  MCMC_input <- list(Y_obs = Y_obs, dY = dY, is_NA_Y = is_NA_Y, H_obs = H_obs, dH = dH, is_NA_H = is_NA_H, X = X, Z = Z, c_init = c_init, Alg_list = Alg_list, n_rates = n_rates,
                     n_times = n_times, epsilon = epsilon, Param_list = Param_list, is_null_X = is_null_X, is_null_Z = is_null_Z)
  
  #Select Gibbs algorithm
  is_parametric <- model_list$is_parametric
  
  if(is.null(dH)){
    if(is_parametric){  
      MCMC_output = DisProgrGibbs_Param_dH0(MCMC_input)
    }else{
      is_DP <- model_list$is_DP
      if(is_DP){
        MCMC_output = DisProgrGibbs_DP_dH0(MCMC_input)
      }else{
        MCMC_output = DisProgrGibbs_Mix_dH0(MCMC_input)
      }
    } 
  }else{
    if(is_parametric){  
      MCMC_output = DisProgrGibbs_Param(MCMC_input)
    }else{
      is_DP <- model_list$is_DP
      if(is_DP){
        MCMC_output = DisProgrGibbs_DP(MCMC_input)
      }else{
        MCMC_output = DisProgrGibbs_Mix(MCMC_input)
      }
    }
  }
  return( MCMC_output )   
}


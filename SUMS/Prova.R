#Simulations with modified BDgraph package

rm(list = ls())
set.seed(123)
library("BayesDiseaseProgr")
library("fields")

# #Check some matrix products
# p <- 5
# N <- 10
# Omega <- rwish(n = 1, 2*p, diag(p))
# 
# mu <- rmvnorm(n = 1, mean = rep(0,p), sigma = solve(Omega))
# 
# phi <- rmvnorm(n = N, mean = c(mu), sigma = solve(Omega))
# sum_phi <- colSums(phi)
# 
# exp1 <- 0
# for(i in 1:N){
#   exp1 <- exp1 + mu %*% Omega %*% phi[i,]
# }
# exp2 <- mu %*% Omega %*% sum_phi
# 
# 
# tr1 <- 0
# for(i in 1:N){
#   tr1 <- tr1 + sum(diag( phi[i,] %*% t( phi[i,]) %*% Omega))
# }
# tr2 <- sum(diag(sum_phi %*% t(sum_phi) %*% Omega))




# #Check normalizing constant in formulas
# set.seed(123)
# library("MCMCpack")
# p <- 2
# 
# b <- 5
# aux <- rmvnorm(p, sigma = diag(10))
# D <- aux %*% t(aux)
# D_inv <- solve(D)
# D112 <- D[1,1] - D[1,2]^2 / D[2,2]
# 
# A <- rwish(b + p - 1, D_inv)
# A112 <- A[1,1]
# #- A[1,2]^2/A[2,2]
# 
# #Check formulas
# dwish(A, b + p - 1, D_inv)
# (2^(2*(b+1)/2) * det(D)^(-(b+1)/2) * gamma(b/2) * gamma((b+1)/2) * pi^(1/2))^(-1) * det(A)^((b-2)/2) * exp(-0.5*sum(diag(D %*% A)))
# 
# 
# J_num <- dwish(A112, b+1, 1/D112)/dwish(A, b + p - 1, D_inv) * det(A)^((b-2)/2) * exp(-0.5*sum(diag(D %*% A)))
# 
# J_BDgraph <- (2*pi/D[2,2])^(.5)* (2^(b/2) * D[2,2]^(-b/2) * gamma(b/2)) *A[1,1]^((b-1)/2)*exp(-.5*A112*D112)


# #The functions I need are compiled functions, so they do not come with usual installation
# untar(download.packages(pkgs = "BDgraph", destdir = ".", type = "source")[,2])









#Number of patients
N = 150;
#Number of response processes
pY = 3;
#Number of covariate processes
pH = 7;
#Tot number of porcesses
p0 = pY+pH;

#Set graph structure(s)
d_simul <- 0.05
G0_simul <- matrix(0, p0, p0)
G0_simul[lower.tri(G0_simul)] <- c(runif(p0*(p0-1)/2) <= d_simul)
G0_simul <- G0_simul + t(G0_simul)

#Assume 2-state responses
dY <- c(2, 3, 2)
dH <- c(2, 3, 2, 4, 5, 3, 2)
n_states <- list(dY,dH)
n_rates <- 0
for(h in 1:pY){
  n_rates <- c(n_rates, dY[h] * (dY[h]-1))
}
for(l in 1:pH){
  n_rates <- c(n_rates, dH[l] * (dH[l]-1))
}
n_rates_cum <- cumsum(n_rates)
p_tot <- sum(n_rates)
G_simul <- matrix(0, nrow = p_tot, ncol = p_tot)
diag(G0_simul) <- rep(1,p0)
for(i in 1:p0){
  for(j in i:p0){
    if(G0_simul[i,j]){
      G_simul[c((n_rates_cum[i]+1) : n_rates_cum[i+1]), c((n_rates_cum[j]+1) : n_rates_cum[j+1])] <- 1
      G_simul[c((n_rates_cum[j]+1) : n_rates_cum[j+1]), c((n_rates_cum[i]+1) : n_rates_cum[i+1])] <- 1 #For symmetry
    }
  }
}
diag(G0_simul) <- rep(0,p0)
diag(G_simul) <- rep(0,p_tot)

#Prior specification
nu_simul <- 5
# Psi <- rmvnorm(p_tot, sigma = 0.1*diag(2*p_tot))
# Psi <- Psi %*% t(Psi)
Psi_simul <- diag(p_tot)

#Simulate data
mu_simul <- rep(2,p_tot)
#Simulate data
Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

phi <- rmvnorm(n = N, mean = c(mu_simul), sigma = solve(Omega_simul))

data <- list(phi = phi, n_states = n_states, n_rates = n_rates)

#Run BDMCMC algorithm for GGMs
n_burn1 <- 10
n_burn2 <- 10
n_save <- 10
thin <- 1
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8
SM_alg <- FALSE
only_SM <- FALSE
Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold, SM_alg = SM_alg)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

#Prior on mu
update_mu <- TRUE
if(update_mu){# We provide hyperparameters for mu
  m_mu <- rep(0,p_tot)
  mu <- m_mu
  mu_list <- list(update_mu = update_mu, mu = mu, m_mu = m_mu)
}else{# We fix the values of mu
  mu <- mu_simul
  mu_list <- list(update_mu = update_mu, mu = mu)
}
Param_list <- modifyList(mu_list,Param_list)

#Prior on m_mu
update_m_mu <- FALSE
if(update_mu){# We provide hyperparameters for m_mu
  Sigma_m_mu <- diag(p_tot)
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu, Sigma_m_mu = Sigma_m_mu)
}else{# We fix the values of m_mu
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu)
}
Param_list <- modifyList(m_mu_list,Param_list)

#Prior on k0
update_k0 <- TRUE
if(update_k0){# We provide hyperparameters for k0
  mm <- 1
  ss <- 10
  a_k0 <- mm^2/ss
  b_k0 <- mm/ss
  k0 <- mm
  k0_list <- list(update_k0 = update_k0, k0 = k0, a_k0 = a_k0, b_k0 = b_k0)
}else{# We fix the values of k0
  k0 <- k0_simul
  k0_list <- list(update_k0 = update_k0, k0 = k0)
}
Param_list <- modifyList(k0_list,Param_list)

#Prior on graph
size_based_prior <- FALSE
if(size_based_prior){
  #The name of the parameters is the same as for d, because it's like marginalizeing wrt d
  a_d <- 1
  b_d <- 1
  G_list <- list(size_based_prior = size_based_prior, a_d = a_d, b_d = b_d)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_d <- TRUE
  if(update_d){
    #Different prior options for d
    d_beta_prior <- FALSE
    if(d_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_d <- mm*(mm*(1 - mm)/ss - 1)
      b_d <- (mm*(1 - mm)^2/ss - (1 - mm))
      d <- mm
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_d = a_d, b_d = b_d, d = d)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      d <- (1 + exp(-mm))^(-1)
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_lambda = a_lambda, d = d)
    }
  }else{# We fix the values of d
    d <- 0.25
    d_list <- list(update_d = update_d, d = d)
  }
  Param_list <- modifyList(d_list,Param_list)
}
Param_list <- modifyList(G_list,Param_list)


# Main Gibbs sampler (only for covariance part, to test and for the future!)
n_rep <- 1
diff_time <- rep(0, n_rep)
for(i_rep in 1:n_rep){
  start_time <- Sys.time()
  sample.bdmcmc <- DisProgrCov( data = data, Alg_list = Alg_list, Param_list = Param_list)
  end_time <- Sys.time()
  diff_time[i_rep] <- end_time - start_time
}
plot(diff_time)
mean(diff_time)


####
#Some plots to check
####
Omega_out <- sample.bdmcmc$Omega_out
G_out <- sample.bdmcmc$G_out
G0_out <- sample.bdmcmc$G0_out

Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

layout(matrix(c(1,2),nrow = 1, ncol = 2))

col_val <- c(Omega_mean, Omega_simul)
col_map <- colorRampPalette(c("white", "red"))(99)

image(Omega_mean, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_mean, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))
image(Omega_simul, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_simul, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))

image(G_mean)
image(G_simul)

image(G0_mean)
image(G0_simul)



sum_weights_out <- 1/sample.bdmcmc$sum_weights_out
Omega_meanW <- matrix(0,p_tot,p_tot)
G_meanW <- matrix(0,p_tot,p_tot)
G0_meanW <- matrix(0,p0,p0)
for(it in 1:n_save){
  Omega_meanW = Omega_meanW + Omega_out[,,it] * sum_weights_out[it]
  G_meanW = G_meanW + G_out[,,it] * sum_weights_out[it]
  G0_meanW = G0_meanW + G0_out[,,it] * sum_weights_out[it]
}
Omega_meanW <- Omega_meanW / sum(sum_weights_out)
G_meanW <- G_meanW / sum(sum_weights_out)
G0_meanW <- G0_meanW / sum(sum_weights_out)

image(Omega_meanW)
image(Omega_simul)

image(G_meanW)
image(G_simul)

image(G0_meanW)
image(G0_simul)





mu_mean <- apply(sample.bdmcmc$mu_out, 2, mean)
plot(c(1:p_tot), mu_simul, col = "red", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~mu))
lines(c(1:p_tot), mu_mean, col = "blue")

m_mu_mean <- apply(sample.bdmcmc$m_mu_out, 2, mean)
plot(c(1:p_tot), m_mu_simul, col = "red", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~m[mu]))
lines(c(1:p_tot), m_mu_mean, col = "blue")

d_out <- sample.bdmcmc$d_out
hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20)
abline(v = d_simul, lwd = 3, col = "red")

k0_out <- sample.bdmcmc$k0_out
hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), breaks = 20)
abline(v = k0_simul, lwd = 3, col = "red")

sizeG0_out <- sample.bdmcmc$sizeG0_out
hist(sizeG0_out, col = "lightblue", main = bquote("Posterior of graph size"), breaks = 20)
abline(v = sum(G0_simul)/2, lwd = 3, col = "red")
plot(sizeG0_out, type = "l")
abline(h = sum(G0_simul)/2, lwd = 3, col = "red")



















#################################################################################################################
#################################################################################################################

## Here we do a fuller simulation from a DP model where we also have the Markov processes
## We will try to handle directly different time points for different processes and patients

rm(list = ls())
set.seed(123)
library("BayesDiseaseProgr")
library("ggplot2")
library("msm")
library("fields")


######
#Simulate sensible data
######

#Number of patients
N = 50;
#Number of response processes
pY = 2;
#Number of covariate processes
pH = 3;
#Tot number of porcesses
p0 = pY+pH;

#Set graph structure(s)
# d_simul <- 0.25
G0_simul <- matrix(0, p0, p0)
# G0_simul[lower.tri(G0_simul)] <- c(runif(p0*(p0-1)/2) <= d_simul)
G0_simul[1,2] <- 1
G0_simul[1,4] <- 1
G0_simul[2,3] <- 1
G0_simul[3,4] <- 1
G0_simul[4,5] <- 1
G0_simul[2,5] <- 1

G0_simul <- G0_simul + t(G0_simul)

#Assume 2-state responses
dY <- c(2, 2)
dH <- c(2, 2, 2)
n_states <- list(dY,dH)
n_rates <- 0
for(h in 1:pY){
  n_rates <- c(n_rates, dY[h] * (dY[h]-1))
}
for(l in 1:pH){
  n_rates <- c(n_rates, dH[l] * (dH[l]-1))
}
n_rates_cum <- cumsum(n_rates)
p_tot <- sum(n_rates)
G_simul <- matrix(0, nrow = p_tot, ncol = p_tot)
diag(G0_simul) <- rep(1,p0)
for(i in 1:p0){
  for(j in i:p0){
    if(G0_simul[i,j]){
      G_simul[c((n_rates_cum[i]+1) : n_rates_cum[i+1]), c((n_rates_cum[j]+1) : n_rates_cum[j+1])] <- 1
      G_simul[c((n_rates_cum[j]+1) : n_rates_cum[j+1]), c((n_rates_cum[i]+1) : n_rates_cum[i+1])] <- 1 #For symmetry
    }
  }
}
diag(G0_simul) <- rep(0,p0)
diag(G_simul) <- rep(0,p_tot)

#Prior specification
nu_simul <- 5
Psi_simul <- diag(p_tot)

#Simulate data
Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

K_N_simul <- 3
c_simul <- sample.int(K_N_simul, size = N, replace = TRUE, prob = c(1/3,1/6,1/2))
phi_simul <- rmvnorm(n = K_N_simul, mean = c(mu_simul), sigma = solve(Omega_simul))
lambda_simul <- exp(phi_simul)

#Number of visits per patient per process (for msm package)
n_times <- matrix(9 + rpois(N*p0, 1), N, p0)
#Times of visit (0 and Tmax included!)
T_max <- 30
times <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i,]
  times_i <- list()
  for(ip in 1:p0){
    times_i[[ip]] <- sort(c(0,runif(n_times_i[ip]-2, min = 0, max = T_max), T_max))
  }
  times[[i]] <- times_i
}

#Transition rates for Y are covariate- and time-dependent
#Transition rates for processes H are constant (and H are binary processes)
lambda_Y_tilde_simul <- vector("list", length = K_N_simul)
lambda_H_simul <- vector("list", length = K_N_simul)
for(j in 1:K_N_simul){
  lambda_Y_tilde_simul_j <- list()
  for(h in 1:pY){
    lambda_Y_tilde_simul_j[[h]] <- lambda_simul[j,(n_rates_cum[h]+1):n_rates_cum[h+1]]
  }
  lambda_Y_tilde_simul[[j]] <- lambda_Y_tilde_simul_j
  
  lambda_H_simul_j <- list()
  for(l in 1:pH){
    lambda_H_simul_j[[l]] <- lambda_simul[j,(n_rates_cum[pY+l]+1):n_rates_cum[pY+l+1]]
  }
  lambda_H_simul[[j]] <- lambda_H_simul_j
}

#Constant covariates, but they can differ in each process i,...,pY
g <- c(2,4) #Simulate age and sex without intercept (it's the lambda_tilde!)
X <- list()

aux <- matrix(1,N,g[1])
aux[,1] <- rnorm(N, mean = 0, sd = 1) #standardized
aux[,2] <- (runif(N) <= 0.5)
X[[1]] <- aux

X[[2]] <- matrix(rgamma(N*g[2], 1, 1), N, g[2])


#Coefficient beta (different for each rate and process 1,...,pY)
beta_simul <- vector("list", length = pY)
for(h in 1:pY){
  beta_simul[[h]] <- matrix(rnorm(g[h]*dY[h]*(dY[h] - 1), sd = 0.1),g[h],dY[h]*(dY[h] - 1),byrow = TRUE)
}

#Time-varying covariates (not random for msm package, no need for averaging)
#Also depends on the processes
#Here we put one equal to zero
q <- c(2,3)
Z_msm <- vector("list", length = pY)
Z <- vector("list", length = pY)

#Careful with dimension q!
Z_msm_h <- vector("list", length = N)
Z_h <- vector("list", length = N)
h = 1
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm(times_i/T_max, sd = sqrt(0.5))
  Z2 <- rnorm(cos(times_i/T_max*2*pi), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h

h = 2
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rep(0,n_times_i)
  Z2 <- rep(0,n_times_i)
  Z3 <- rep(0,n_times_i)
  Z_msm_h[[i]] <- cbind(Z1,Z2,Z3)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z3 <- (Z3[1:(n_times_i-1)] + Z3[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2,Z3)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h


#Coefficient gamma (different for each rate and process 1,...,pY)
gamma_simul <- vector("list", length = pY)
for(h in 1:pY){
  gamma_simul[[h]] <- matrix(rnorm(q[h]*dY[h]*(dY[h] - 1), sd = 0.1),q[h],dY[h]*(dY[h] - 1),byrow = TRUE)
}

#Constant part of transition rates (needed in sim.msm, the time-varying covariates are included in the function used for the simulation)
lambda_Y_msm <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_msm_i <- list()
  for(h in 1:pY){
    lambda_Y_msm_i[[h]] <- lambda_Y_tilde_simul[[c_simul[i]]][[h]] * exp( X[[h]][i,] %*% beta_simul[[h]] )
  }
  lambda_Y_msm[[i]] <- lambda_Y_msm_i
}

#Generate the data (...a bit tricky with msm...
#...but we can do it independently and at different times for each process and patient...)

#Response process
Y <- vector("list", length = N)
lambda_Y_hat <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i,]
  
  lambda_Y_msm_i <- lambda_Y_msm[[i]]
  
  Y_i <- list()
  lambda_Y_hat_i <- list()
  for(h in 1:pY){
    
    lambda_Y_msm_i_h <- lambda_Y_msm_i[[h]]
    Q_Y <- matrix(0, dY[h], dY[h])
    for(i2 in 1:dY[h]){
      Q_Y[i2,c(1:dY[h])[-i2]] <- lambda_Y_msm_i_h[(i2-1)*(dY[h]-1) + c(1:(dY[h]-1))]
      Q_Y[i2,i2] <- -sum(Q_Y[i2,])
    }
    
    #Here we can include the time-varying covariates via sim.msm
    times_i_h <- times_i[[h]]
    n_times_i_h <- n_times_i[h]
    
    Z_msm_i_h <- Z_msm[[h]][[i]]
    gamma_simul_h <- gamma_simul[[h]]
    
    Y_i_h_msm <- sim.msm(qmatrix = Q_Y, maxtime = T_max, covs = Z_msm_i_h, beta = gamma_simul_h, obstimes = times_i_h, start = 2, mintime = 0)
    # Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
    
    out_aux <- outer(times_i_h[2:(n_times_i_h-1)], Y_i_h_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    Y_i[[h]] <- Y_i_h_msm$states[c(1,ind-1,length(Y_i_h_msm$times))] - 1 #(0,1)
    
    #For estimation of the initial transition rates we do not need any special q matrix
    #It only has to indicate id there are some constrained zeros in the matrix to be estimated (not for us I think)
    #Problems arise when there are no or only one transitions (zeros)
    data_Y_i_df <- data.frame(states = Y_i[[h]]+1, time = times_i_h)
    lambda_aux <- crudeinits.msm(states ~ time, data = data_Y_i_df, qmatrix = Q_Y)
    #If zeroes we impose arbitrary value
    lambda_aux[lambda_aux == 0] <- 0.5
    lambda_aux_vec <- lambda_aux[1,1]
    for(i2 in 1:dY[h]){
      lambda_aux_vec <- c(lambda_aux_vec, lambda_aux[i2,c(1:dY[h])[-i2]])
    }
    lambda_aux_vec <- lambda_aux_vec[-1]
    lambda_Y_hat_i[[h]] <- lambda_aux_vec
  }
  Y[[i]] <- Y_i
  lambda_Y_hat[[i]] <- lambda_Y_hat_i
}


#Covariate processes
H <- vector("list", length = N)
lambda_H_hat <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i,]
  
  lambda_H_simul_i <- lambda_H_simul[[c_simul[i]]]
  
  H_i <-vector("list", length = pH)
  lambda_H_i_hat <-vector("list", length = pH)
  for(l in 1:pH){
    
    times_i_l <- times_i[[pY+l]]
    n_times_i_l <- n_times_i[pY+l]
    
    lambda_H_simul_i_l <- lambda_H_simul_i[[l]]
    #No covariates
    Q_H <- matrix(0, dH[l], dH[l])
    for(i2 in 1:dH[l]){
      Q_H[i2,c(1:dH[l])[-i2]] <- lambda_H_simul_i_l[(i2-1)*(dH[l]-1) + c(1:(dH[l]-1))]
      Q_H[i2,i2] <- -sum(Q_H[i2,])
    }
    H_i_l_msm <- sim.msm(Q_H, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i_l, start = 2, mintime = 0)
    
    #Exclude first and last time points
    out_aux <- outer(times_i_l[2:(n_times_i_l-1)], H_i_l_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    H_i[[l]] <- H_i_l_msm$states[c(1,ind-1,length(H_i_l_msm$times))]  - 1 #(0,1)
    
    #For estimation of the initial transition rates we do not need any special q matrix
    #It only has to indicate id there are some constrained zeros in the matrix to be estimated (not for us I think)
    #Problems arise when there are no or only one transitions (zeros)
    data_H_i_l_df <- data.frame(states = H_i[[l]]+1, time = times_i_l)
    lambda_aux <- crudeinits.msm(states ~ time, data = data_H_i_l_df, qmatrix = Q_H)
    lambda_aux[lambda_aux == 0] <- 0.5
    lambda_aux_vec <- lambda_aux[1,1]
    for(i2 in 1:dH[l]){
      lambda_aux_vec <- c(lambda_aux_vec, lambda_aux[i2,c(1:dH[l])[-i2]])
    }
    lambda_aux_vec <- lambda_aux_vec[-1]
    lambda_H_i_hat[[l]] <- lambda_aux_vec
  }
  H[[i]] <- H_i
  lambda_H_hat[[i]] <- lambda_H_i_hat
}

#Put the estimated lambdas in a list
lambda_hat_all <- matrix(0,N,p_tot)
for(i in 1:N){
  lambda_hat_i <- 0
  for(h in 1:pY){
    lambda_hat_i <- c(lambda_hat_i, lambda_Y_hat[[i]][[h]])
  }
  for(l in 1:pH){
    lambda_hat_i <- c(lambda_hat_i, lambda_H_hat[[i]][[l]])
  }
  lambda_hat_all[i,] <- lambda_hat_i[2:(p_tot+1)]
}


#Interval lengths
epsilon <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  epsilon_i <- vector("list", length = p0)
  for(ip in 1:p0){
    epsilon_i[[ip]] <- diff(times_i[[ip]])
  }
  epsilon[[i]] <- epsilon_i
}

data <- list(Y = Y, H = H, lambda_hat_all = lambda_hat_all, n_states = n_states, n_rates = n_rates, n_times = n_times, epsilon = epsilon, X = X, Z = Z)

## Now that we have the data, we need to initialize the parameters and priors
n_burn1 <- 100
n_burn2 <- 100
n_save <- 1000
thin <- 1
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8

Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

#Prior on mu
update_mu <- TRUE
if(update_mu){# We provide hyperparameters for mu
  m_mu <- rep(0,p_tot)
  mu <- m_mu
  mu_list <- list(update_mu = update_mu, mu = mu, m_mu = m_mu)
}else{# We fix the values of mu
  mu <- mu_simul
  mu_list <- list(update_mu = update_mu, mu = mu)
}
Param_list <- modifyList(mu_list,Param_list)

#Prior on m_mu
update_m_mu <- FALSE
if(update_mu){# We provide hyperparameters for m_mu
  Sigma_m_mu <- diag(p_tot)
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu, Sigma_m_mu = Sigma_m_mu)
}else{# We fix the values of m_mu
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu)
}
Param_list <- modifyList(m_mu_list,Param_list)

#Prior on k0
update_k0 <- TRUE
if(update_k0){# We provide hyperparameters for k0
  mm <- 1
  ss <- 1
  a_k0 <- mm^2/ss
  b_k0 <- mm/ss
  k0 <- mm
  k0_list <- list(update_k0 = update_k0, k0 = k0, a_k0 = a_k0, b_k0 = b_k0)
}else{# We fix the values of k0
  k0 <- k0_simul
  k0_list <- list(update_k0 = update_k0, k0 = k0)
}
Param_list <- modifyList(k0_list,Param_list)

#Prior on graph
size_based_prior <- TRUE
if(size_based_prior){
  #The name of the parameters is the same as for d, because it's like marginalizeing wrt d
  a_d <- 1
  b_d <- 1
  G_list <- list(size_based_prior = size_based_prior, a_d = a_d, b_d = b_d)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_d <- TRUE
  if(update_d){
    #Different prior options for d
    d_beta_prior <- TRUE
    if(d_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_d <- mm*(mm*(1 - mm)/ss - 1)
      b_d <- (mm*(1 - mm)^2/ss - (1 - mm))
      d <- mm
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_d = a_d, b_d = b_d, d = d)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      d <- (1 + exp(-mm))^(-1)
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_lambda = a_lambda, d = d)
    }
  }else{# We fix the values of d
    d <- 0.25
    d_list <- list(update_d = update_d, d = d)
  }
  Param_list <- modifyList(d_list,Param_list)
}
Param_list <- modifyList(G_list,Param_list)


#Prior on beta
mu_beta <- vector("list", length = pY)
U_beta <- vector("list", length = pY)
V_beta <- vector("list", length = pY)
for(h in 1:pY){
  mu_beta[[h]] <- matrix(0,g[h],dY[h]*(dY[h] - 1))
  U_beta[[h]] <- diag(g[h])
  V_beta[[h]] <- diag(dY[h]*(dY[h] - 1))
}
beta_list <- list(mu_beta = mu_beta, U_beta = U_beta, V_beta = V_beta)
Param_list <- modifyList(beta_list, Param_list)

#Prior on gamma
mu_gamma <- vector("list", length = pY)
U_gamma <- vector("list", length = pY)
V_gamma <- vector("list", length = pY)
for(h in 1:pY){
  mu_gamma[[h]] <- matrix(0,q[h],dY[h]*(dY[h] - 1))
  U_gamma[[h]] <- diag(q[h])
  V_gamma[[h]] <- diag(dY[h]*(dY[h] - 1))
}
gamma_list <- list(mu_gamma = mu_gamma, U_gamma = U_gamma, V_gamma = V_gamma)
Param_list <- modifyList(gamma_list, Param_list)

#A-priori expected number of components
Mm <- 5

is_DP <- TRUE
is_parametric <- FALSE
if(is_parametric){#We can fit a model for a given clustering
  c_init <- c(1:N)-1
  alpha_list <- list(update_alpha = FALSE, alpha = 0)
}else{
  if(is_DP){#We can fit a DP model
    #Prior on alpha
    # #Based on prior expected number of components
    alpha <- matrix(seq(0.001,2,length = 50000),50000,1)
    sum_alpha <- apply(alpha,1,function(alpha) sum(alpha/(alpha + c(1:N) - 1)))
    alpha_Mm <- alpha[which.min(abs(Mm - sum_alpha))]
    sum_alpha[which.min(abs(Mm - sum_alpha))]
    
    update_alpha <- TRUE
    if(update_alpha){# We provide hyperparameters for alpha
      mm <- alpha_Mm
      ss <- 1
      a_alpha <- mm^2/ss
      b_alpha <- mm/ss
      alpha <- mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha, a_alpha = a_alpha, b_alpha = b_alpha)
    }else{# We fix the values of alpha
      alpha <- alpha_Mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha)
    }
  }else{
    alpha_list <- list(update_alpha = FALSE, alpha = 0)
  }
  c_init <- rep(1,N)-1
}
Param_list <- modifyList(alpha_list,Param_list)

model_list <- list(is_parametric = is_parametric, is_DP = is_DP , c_init = c_init)

# Main Gibbs sampler
sample.bdmcmc <- DisProgr( data = data, model_list = model_list, Alg_list = Alg_list, Param_list = Param_list)




#Produce plots for report

Graph_list <- sample.bdmcmc$Graph_List

Omega_out <- Graph_list$Omega_out
G_out <- Graph_list$G_out
G0_out <- Graph_list$G0_out

Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

layout(matrix(c(1,2),nrow = 1, ncol = 2))

col_val <- c(Omega_mean, Omega_simul)
col_map <- colorRampPalette(c("white", "red"))(99)

image(Omega_mean, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_mean, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))
image(Omega_simul, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_simul, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))

image(G_mean)
image(G_simul)

image(G0_mean)
image(G0_simul)



sum_weights_out <- 1/Graph_list$sum_weights_out
Omega_meanW <- matrix(0,p_tot,p_tot)
G_meanW <- matrix(0,p_tot,p_tot)
G0_meanW <- matrix(0,p0,p0)
for(it in 1:n_save){
  Omega_meanW = Omega_meanW + Omega_out[,,it] * sum_weights_out[it]
  G_meanW = G_meanW + G_out[,,it] * sum_weights_out[it]
  G0_meanW = G0_meanW + G0_out[,,it] * sum_weights_out[it]
}
Omega_meanW <- Omega_meanW / sum(sum_weights_out)
G_meanW <- G_meanW / sum(sum_weights_out)
G0_meanW <- G0_meanW / sum(sum_weights_out)

image(Omega_meanW)
image(Omega_simul)

image(G_meanW)
image(G_simul)

image(G0_meanW)
image(G0_simul)

dev.off()

#d
d_out <- Graph_list$d_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20)
abline(v = d_simul, lwd = 3, col = "red")
plot(d_out, type = "l")

#mu
mu_out <- sample.bdmcmc$mu_out
matplot(t(mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~mu), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#m_mu
m_mu_out <- sample.bdmcmc$m_mu_out
matplot(t(m_mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~m[mu]), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(m_mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), m_mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(m_mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#k0
k0_out <- sample.bdmcmc$k0_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), breaks = 20)
abline(v = k0_simul, lwd = 3, col = "red")
plot(k0_out, type = "l")



#Phi
phi_star_out <- sample.bdmcmc$phi_star_out
c_out <- sample.bdmcmc$c_out
phi_mean <- matrix(0, N, p_tot)
for(it in 1:n_save){
  phi_mean <- phi_mean + phi_star_out[[it]][c_out[it,]+1,]
}
phi_mean <- phi_mean/n_save
matplot(c(1:p_tot), t(phi_mean), col = "grey", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~phi), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
matplot(c(1:p_tot), t(phi_simul), col = "red", lwd = 1.5, lty = 1, type = "l", add = TRUE)
legend("topright", legend = c("MCMC", "truth"), pch = 19, col = c("grey", "red"), bty = "n", cex = 2)


#Number of clusters
K_N_out <- rep(0,n_save)
for(it in 1:n_save){
  K_N_out[it] <- dim(phi_star_out[[it]])[1]
}
plot(table(K_N_out)/n_save, col = "blue", lwd = 2)


#alpha
alpha_out <- sample.bdmcmc$alpha_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(alpha_out, col = "lightblue", main = bquote("Posterior of "~alpha), breaks = 20)
plot(alpha_out, type = "l")



#beta
beta_out <- sample.bdmcmc$beta_out
for(h in 1:pY){
  beta_out_h <- array(0, dim = c(g[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    beta_out_h[,,it] <- beta_out[[it]][[h]]
  }
  
  beta_df <- data.frame(dim_g = rep(c(1:(g[h]*dY[h]*(dY[h] - 1))), n_save), beta_vec = c(beta_out_h))
  beta_true_df <- data.frame( x = c(1:(g[h]*dY[h]*(dY[h] - 1))), y = c(beta_simul[[h]]) ) 
  
  ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
    geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
    labs(x = "g x dY", y = "", title = bquote("Posterior of "~beta)) +
    theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
}


#gamma
gamma_out <- sample.bdmcmc$gamma_out
for(h in 1:pY){
  gamma_out_h <- array(0, dim = c(q[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    gamma_out_h[,,it] <- gamma_out[[it]][[h]]
  }
  
  gamma_df <- data.frame(dim_g = rep(c(1:(q[h]*dY[h]*(dY[h] - 1))), n_save), gamma_vec = c(gamma_out_h))
  gamma_true_df <- data.frame( x = c(1:(q[h]*dY[h]*(dY[h] - 1))), y = c(gamma_simul[[h]]) ) 
  
  ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
    geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
    labs(x = "g x dY", y = "", title = bquote("Posterior of "~gamma)) +
    theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
}


#To see how the transition rates are estimated I want to show the evolution of trajectories

phi_star_out <- sample.bdmcmc$phi_star_out
c_out <- sample.bdmcmc$c_out
beta_out <- sample.bdmcmc$beta_out
gamma_out <- sample.bdmcmc$gamma_out

#Reconstruct transition rates for plots
#It's parametric so phi_star_out is a matrix
lambda_Y_out <- vector("list", length = n_save)
lambda_H_out <- vector("list", length = n_save)
for(it in 1:n_save){
  phi_star_it_out <- phi_star_out[[it]]
  
  beta_out_it <- beta_out[[it]]
  gamma_out_it <- gamma_out[[it]]
  
  lambda_Yi_out <- vector("list", length = N)
  lambda_Hi_out <- vector("list", length = N)
  for(i in 1:N){
    n_times_i <- n_times[i,]
    lambda_Y_i <- vector("list", length = pY)
    for(h in 1:pY){
      ind_dY <- (n_rates_cum[h]+1):n_rates_cum[h+1]
      
      n_times_i_h <- n_times_i[h]
      lambda_Y_i_h <- matrix(0, n_times_i_h, dY[h]*(dY[h]-1))
      
      X_h <- X[[h]]
      Z_h <- Z[[h]]
      Z_h_i <- Z_h[[i]]
      beta_it_h_out <- beta_out_it[[h]]
      gamma_it_h_out <- gamma_out_it[[h]]
      
      lambda_Y_i[[h]] <- matrix(exp( phi_star_it_out[c_out[it,i]+1,ind_dY] + X_h[i,] %*% beta_it_h_out), nrow = n_times_i[h]-1, ncol = dY[h]*(dY[h]-1), byrow = TRUE) * exp( Z_h_i %*% gamma_it_h_out )
    }
    lambda_H_i <- vector("list", length = pH)
    for(l in 1:pH){
      ind_dH <- (n_rates_cum[pY+l]+1):n_rates_cum[pY+l+1]
      lambda_H_i[[l]] <- exp(phi_star_it_out[c_out[it,i]+1,ind_dH])
    }
    lambda_Yi_out[[i]] <- lambda_Y_i
    lambda_Hi_out[[i]] <- lambda_H_i
  }
  lambda_Y_out[[it]] <- lambda_Yi_out
  lambda_H_out[[it]] <- lambda_Hi_out
}

#Reconstruct simulated transition rates
lambda_Y_simul <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i,]
  lambda_Y_simul_i <- vector("list", length = pY)
  for(h in 1:pY){
    n_times_i_h <- n_times_i[h]
    X_h <- X[[h]]
    Z_h <- Z[[h]]
    Z_h_i <- Z_h[[i]]
    beta_simul_h <- beta_simul[[h]]
    gamma_simul_h <- gamma_simul[[h]]
    lambda_Y_simul_i[[h]] <- matrix(lambda_Y_tilde_simul[[c_simul[i]]][[h]] * exp( X_h[i,] %*% beta_simul_h ), n_times_i_h - 1, dY[h] * (dY[h] - 1), byrow = TRUE) * exp( Z_h_i %*% gamma_simul_h )
  }
  lambda_Y_simul[[i]] <- lambda_Y_simul_i
}


#Select patient
i_plot <- 28
Y_i <- Y[[i_plot]]
H_i <- H[[i_plot]]
times_i_plot <- times[[i_plot]]
n_times_i_plot <- n_times[i_plot,]


#Simulated transition probabilities
prs_simul_log <- vector("list", length(p0))

lambda_Y_simul_i_plot <- lambda_Y_simul[[i_plot]]
epsilon_i_plot <- epsilon[[i_plot]]
for(h in 1:pY){
  n_times_i_h_plot <- n_times_i_plot[h]
  prs_simul_log_h <- rep(0,n_times_i_h_plot)
  Y_i_h <- Y_i[[h]]
  lambda_Y_simul_i_h_plot <- lambda_Y_simul_i_plot[[h]]
  epsilon_i_h_plot <- epsilon_i_plot[[h]]
  for(j in 2:n_times_i_h_plot){
    Yihj <- Y_i_h[j-1]
    prs_simul_log_h[j] <- log(lambda_Y_simul_i_h_plot[j-1,Yihj+1]) - log(sum(lambda_Y_simul_i_h_plot[j-1,])) + log( 1 - exp(- ((sum(lambda_Y_simul_i_h_plot[j-1,])) * epsilon_i_h_plot[j-1])) )
  }
  prs_simul_log[[h]] <- prs_simul_log_h
}

lambda_H_simul_i_plot <- lambda_H_simul[[c_simul[i_plot]]]
for(l in 1:pH){
  n_times_i_l_plot <- n_times_i_plot[pY+l]
  prs_simul_log_l <- rep(0,n_times_i_l_plot)
  H_i_l <- H_i[[l]]
  lambda_H_simul_i_l_plot <- lambda_H_simul_i_plot[[l]]
  epsilon_i_l_plot <- epsilon_i_plot[[pY+l]]
  for(j in 2:n_times_i_l_plot){
    Hilj <- H_i_l[j-1]
    prs_simul_log_l[j] <- log(lambda_H_simul_i_l_plot[Hilj+1]) - log(sum(lambda_H_simul_i_l_plot)) + log( 1 - exp(- (sum(lambda_H_simul_i_l_plot) * epsilon_i_l_plot[j-1])) )
  }
  prs_simul_log[[pY+l]] <- prs_simul_log_l
}

#Estimated transition probabilities
prs_out_log <- vector("list", length = p0)

for(h in 1:pY){
  n_times_i_h_plot <- n_times_i_plot[h]
  prs_out_log_h <- matrix(0,n_times_i_h_plot,n_save)
  Y_i_h <- Y_i[[h]]
  epsilon_i_h_plot <- epsilon_i_plot[[h]]
  for(j in 2:n_times_i_h_plot){
    Yihj <- Y_i_h[j-1]
    for(it in 1:n_save){
      lambda_Y_it_i_plot_h_out <- lambda_Y_out[[it]][[i_plot]][[h]]
      prs_out_log_h[j,it] <- log(lambda_Y_it_i_plot_h_out[j-1,Yihj+1]) - log(sum(lambda_Y_it_i_plot_h_out[j-1,])) + log( 1 - exp(- (sum(lambda_Y_it_i_plot_h_out[j-1,]) * epsilon_i_h_plot[j-1])) )
    }
  }
  prs_out_log[[h]] <- prs_out_log_h
}

for(l in 1:pH){
  n_times_i_l_plot <- n_times_i_plot[pY+l]
  prs_out_log_l <- matrix(0,n_times_i_l_plot,n_save)
  H_i_l <- H_i[[l]]
  epsilon_i_l_plot <- epsilon_i_plot[[pY+l]]
  for(j in 2:n_times_i_l_plot){
    Hilj <- H_i_l[j-1]
    for(it in 1:n_save){
      lambda_H_it_i_plot_l_out <- lambda_H_out[[it]][[i_plot]][[l]]
      prs_out_log_l[j,it] <- log(lambda_H_it_i_plot_l_out[Hilj+1]) - log(sum(lambda_H_it_i_plot_l_out)) + log( 1 - exp(- (sum(lambda_H_it_i_plot_l_out) * epsilon_i_l_plot[j-1])) )
    }
  }
  prs_out_log[[pY+l]] <- prs_out_log_l
}


#Colours of processes
col_proc <- rgb(255/255,51/255,51/255)

scale_red1 <- seq(0, 200, length = pY)
scale_red2 <- seq(0, 255, length = pY)
for(h in 1:pY){
  col_proc <- c(col_proc, rgb(255/255,scale_red1[h]/255,scale_red2[h]/255))
}
scale_blue1 <- seq(20, 200, length = pH)
scale_blue2 <- seq(120, 255, length = pH)
for(l in 1:pH){
  col_proc <- c(col_proc, rgb(scale_blue1[l]/255,scale_blue2[l]/255,255/255))
}
col_proc <- col_proc[-1]

#Labels of processes
proc_lab <- c(paste0("Y", c(1:pY)), paste0("H", c(1:pH)))

#Credibility intervals
prs_CI_low <- vector("list", length = p0)
prs_CI_up <- vector("list", length = p0)
for(ip in 1:p0){
  prs_CI_low[[ip]] <- apply(prs_out_log[[ip]],1, quantile, probs = 0.025)
  prs_CI_up[[ip]] <- apply(prs_out_log[[ip]],1, quantile, probs = 0.975)
}
prs_CI_df <- data.frame(prs_CI_low = unlist(prs_CI_low), prs_CI_up = unlist(prs_CI_up), proc = unlist(mapply(rep,proc_lab,n_times_i_plot), use.names = FALSE), times_vec = round(unlist(times_i_plot),3) )
prs_simul_df <- data.frame(x_times = unlist(times_i_plot), prs_simul_vec = unlist(prs_simul_log), proc = unlist(mapply(rep,proc_lab,n_times_i_plot), use.names = FALSE), times_vec = round(unlist(times_i_plot),3) )
prs_df <- data.frame(x_times = unlist(mapply(rep,times_i_plot,n_save)), prs_vec = unlist(prs_out_log), proc = unlist(mapply(rep,proc_lab,n_times_i_plot * n_save), use.names = FALSE), times_vec = round(unlist(mapply(rep,times_i_plot,n_save)),3) )


ggplot(prs_df, aes(x = x_times, y = prs_vec, col = proc)) + 
  geom_line(data = prs_simul_df, aes(x = x_times, y = prs_simul_vec, col = proc), size = 1.25, linetype = 2) + 
  stat_summary(fun.y = 'mean', geom = 'line', size = 1.25) +
  stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2) + 
  labs(x = "time", y = bquote(log(p[t]^{rs})), title = bquote("Posterior of "~log(p[t]^{rs})~" for subject "~.(i_plot))) +
  theme(text = element_text(size=35), axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 15)) +
  facet_grid(rows = vars(proc), scales="free") +
  # the labels must match what you specified above
  scale_fill_manual(name = proc_lab, values = col_proc) +
  scale_color_manual(name = proc_lab, values = col_proc) +
  theme_bw()





#####Binder loss function equal costs
K_Binder <- 1/2
#For selected combinations of parameters
c_out <- sample.bdmcmc$c_out

pij <- matrix(0,N,N)
cij <- vector("list", length = n_save)
for(g in 1:n_save){
  cij[[g]] <- outer(c_out[g,], c_out[g,], "==")
  pij <- pij + cij[[g]]
}
pij <- pij/n_save


Binder_f <- rep(0,n_save)
for(g in 1:n_save){
  aux <- (pij - K_Binder) * as.matrix(cij[[g]])
  aux <-  aux[upper.tri(aux)]
  Binder_f[g] <- sum(aux)
}
Binder_ind <- which.max(Binder_f)
Binder_out <- c_out[Binder_ind,] + 1

graphics.off()
library("gplots")
Binder_sort <- sort(Binder_out, index.return = TRUE)
Binder_sort <- Binder_sort$ix

heatmap.2(pij[Binder_sort,Binder_sort], ColSideColor = rainbow(K_N_simul)[c_simul[Binder_sort]], col=rev(heat.colors(100)), trace = "none", key  = TRUE, Rowv = FALSE, Colv = FALSE, dendrogram = "none", main = "Posterior co-clustering probabilities", cex.main = 3, cex.lab = 2, cex.axis = 1.5)


















#################################################################################################################
#################################################################################################################

## Here we do a fuller simulation from an AntMAN Mixture model where we also have the Markov processes
## We will try to handle directly different time points for different processes and patients
## We will try also SM algorithm


rm(list = ls())
set.seed(123)
library("BayesDiseaseProgr")
library("ggplot2")
library("msm")
library("fields")


######
#Simulate sensible data
######

#Number of patients
N = 150;
#Number of response processes
pY = 2;
#Number of covariate processes
pH = 3;
#Tot number of porcesses
p0 = pY+pH;

#Set graph structure(s)
# d_simul <- 0.25
G0_simul <- matrix(0, p0, p0)
# G0_simul[lower.tri(G0_simul)] <- c(runif(p0*(p0-1)/2) <= d_simul)
G0_simul[1,2] <- 1
G0_simul[2,3] <- 1
G0_simul[3,4] <- 1
G0_simul[2,5] <- 1

G0_simul <- G0_simul + t(G0_simul)

#Assume 2-state responses
dY <- c(2, 2)
dH <- c(2, 2, 2)
n_states <- list(dY,dH)
n_rates <- 0
for(h in 1:pY){
  n_rates <- c(n_rates, dY[h] * (dY[h]-1))
}
for(l in 1:pH){
  n_rates <- c(n_rates, dH[l] * (dH[l]-1))
}
n_rates_cum <- cumsum(n_rates)
p_tot <- sum(n_rates)
G_simul <- matrix(0, nrow = p_tot, ncol = p_tot)
diag(G0_simul) <- rep(1,p0)
for(i in 1:p0){
  for(j in i:p0){
    if(G0_simul[i,j]){
      G_simul[c((n_rates_cum[i]+1) : n_rates_cum[i+1]), c((n_rates_cum[j]+1) : n_rates_cum[j+1])] <- 1
      G_simul[c((n_rates_cum[j]+1) : n_rates_cum[j+1]), c((n_rates_cum[i]+1) : n_rates_cum[i+1])] <- 1 #For symmetry
    }
  }
}
diag(G0_simul) <- rep(0,p0)
diag(G_simul) <- rep(0,p_tot)

#Prior specification
nu_simul <- p0 + 2
Psi_simul <- diag(p_tot)

#Simulate data
Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

K_N_simul <- 3
c_simul <- sample.int(K_N_simul, size = N, replace = TRUE, prob = c(1/3,1/6,1/2))
phi_simul <- rmvnorm(n = K_N_simul, mean = c(mu_simul), sigma = solve(Omega_simul))
lambda_simul <- exp(phi_simul)

#Number of visits per patient per process (for msm package)
n_times <- matrix(9 + rpois(N*p0, 1), N, p0)
#Times of visit (0 and Tmax included!)
T_max <- 30
times <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i,]
  times_i <- list()
  for(ip in 1:p0){
    times_i[[ip]] <- sort(c(0,runif(n_times_i[ip]-2, min = 0, max = T_max), T_max))
  }
  times[[i]] <- times_i
}

#Transition rates for Y are covariate- and time-dependent
#Transition rates for processes H are constant (and H are binary processes)
lambda_Y_tilde_simul <- vector("list", length = K_N_simul)
lambda_H_simul <- vector("list", length = K_N_simul)
for(j in 1:K_N_simul){
  lambda_Y_tilde_simul_j <- list()
  for(h in 1:pY){
    lambda_Y_tilde_simul_j[[h]] <- lambda_simul[j,(n_rates_cum[h]+1):n_rates_cum[h+1]]
  }
  lambda_Y_tilde_simul[[j]] <- lambda_Y_tilde_simul_j
  
  lambda_H_simul_j <- list()
  for(l in 1:pH){
    lambda_H_simul_j[[l]] <- lambda_simul[j,(n_rates_cum[pY+l]+1):n_rates_cum[pY+l+1]]
  }
  lambda_H_simul[[j]] <- lambda_H_simul_j
}

#Constant covariates, but they can differ in each process i,...,pY
g <- c(2,4) #Simulate age and sex without intercept (it's the lambda_tilde!)
X <- list()

aux <- matrix(1,N,g[1])
aux[,1] <- rnorm(N, mean = 0, sd = 1) #standardized
aux[,2] <- (runif(N) <= 0.5)
X[[1]] <- aux

X[[2]] <- matrix(rgamma(N*g[2], 1, 1), N, g[2])


#Coefficient beta (different for each rate and process 1,...,pY)
beta_simul <- vector("list", length = pY)
for(h in 1:pY){
  beta_simul[[h]] <- matrix(rnorm(g[h]*dY[h]*(dY[h] - 1), sd = 0.1),g[h],dY[h]*(dY[h] - 1),byrow = TRUE)
}

#Time-varying covariates (not random for msm package, no need for averaging)
#Also depends on the processes
#Here we put one equal to zero
q <- c(2,3)
Z_msm <- vector("list", length = pY)
Z <- vector("list", length = pY)

#Careful with dimension q!
Z_msm_h <- vector("list", length = N)
Z_h <- vector("list", length = N)
h = 1
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm(times_i/T_max, sd = sqrt(0.5))
  Z2 <- rnorm(cos(times_i/T_max*2*pi), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h

h = 2
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm(sqrt(times_i/T_max), sd = sqrt(0.5))
  Z2 <- rnorm((times_i/T_max)^2, sd = sqrt(0.5))
  Z3 <- rnorm((times_i/T_max)^3, sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2,Z3)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z3 <- (Z3[1:(n_times_i-1)] + Z3[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2,Z3)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h


#Coefficient gamma (different for each rate and process 1,...,pY)
gamma_simul <- vector("list", length = pY)
for(h in 1:pY){
  gamma_simul[[h]] <- matrix(rnorm(q[h]*dY[h]*(dY[h] - 1), sd = 0.1),q[h],dY[h]*(dY[h] - 1),byrow = TRUE)
}

#Constant part of transition rates (needed in sim.msm, the time-varying covariates are included in the function used for the simulation)
lambda_Y_msm <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_msm_i <- list()
  for(h in 1:pY){
    lambda_Y_msm_i[[h]] <- lambda_Y_tilde_simul[[c_simul[i]]][[h]] * exp( X[[h]][i,] %*% beta_simul[[h]] )
  }
  lambda_Y_msm[[i]] <- lambda_Y_msm_i
}

#Generate the data (...a bit tricky with msm...
#...but we can do it independently and at different times for each process and patient...)

#Response process
Y <- vector("list", length = N)
lambda_Y_hat <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i,]
  
  lambda_Y_msm_i <- lambda_Y_msm[[i]]
  
  Y_i <- list()
  lambda_Y_hat_i <- list()
  for(h in 1:pY){
    
    lambda_Y_msm_i_h <- lambda_Y_msm_i[[h]]
    Q_Y <- matrix(0, dY[h], dY[h])
    for(i2 in 1:dY[h]){
      Q_Y[i2,c(1:dY[h])[-i2]] <- lambda_Y_msm_i_h[(i2-1)*(dY[h]-1) + c(1:(dY[h]-1))]
      Q_Y[i2,i2] <- -sum(Q_Y[i2,])
    }
    
    #Here we can include the time-varying covariates via sim.msm
    times_i_h <- times_i[[h]]
    n_times_i_h <- n_times_i[h]
    
    Z_msm_i_h <- Z_msm[[h]][[i]]
    gamma_simul_h <- gamma_simul[[h]]
    
    Y_i_h_msm <- sim.msm(qmatrix = Q_Y, maxtime = T_max, covs = Z_msm_i_h, beta = gamma_simul_h, obstimes = times_i_h, start = 2, mintime = 0)
    # Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
    
    out_aux <- outer(times_i_h[2:(n_times_i_h-1)], Y_i_h_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    Y_i[[h]] <- Y_i_h_msm$states[c(1,ind-1,length(Y_i_h_msm$times))] - 1 #(0,1)
    
    #For estimation of the initial transition rates we do not need any special q matrix
    #It only has to indicate id there are some constrained zeros in the matrix to be estimated (not for us I think)
    #Problems arise when there are no or only one transitions (zeros)
    data_Y_i_df <- data.frame(states = Y_i[[h]]+1, time = times_i_h)
    lambda_aux <- crudeinits.msm(states ~ time, data = data_Y_i_df, qmatrix = Q_Y)
    #If zeroes we impose arbitrary value
    lambda_aux[lambda_aux == 0] <- 1
    lambda_aux_vec <- lambda_aux[1,1]
    for(i2 in 1:dY[h]){
      lambda_aux_vec <- c(lambda_aux_vec, lambda_aux[i2,c(1:dY[h])[-i2]])
    }
    lambda_aux_vec <- lambda_aux_vec[-1]
    lambda_Y_hat_i[[h]] <- lambda_aux_vec
  }
  Y[[i]] <- Y_i
  lambda_Y_hat[[i]] <- lambda_Y_hat_i
}


#Covariate processes
H <- vector("list", length = N)
lambda_H_hat <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i,]
  
  lambda_H_simul_i <- lambda_H_simul[[c_simul[i]]]
  
  H_i <-vector("list", length = pH)
  lambda_H_i_hat <-vector("list", length = pH)
  for(l in 1:pH){
    
    times_i_l <- times_i[[pY+l]]
    n_times_i_l <- n_times_i[pY+l]
    
    lambda_H_simul_i_l <- lambda_H_simul_i[[l]]
    #No covariates
    Q_H <- matrix(0, dH[l], dH[l])
    for(i2 in 1:dH[l]){
      Q_H[i2,c(1:dH[l])[-i2]] <- lambda_H_simul_i_l[(i2-1)*(dH[l]-1) + c(1:(dH[l]-1))]
      Q_H[i2,i2] <- -sum(Q_H[i2,])
    }
    H_i_l_msm <- sim.msm(Q_H, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i_l, start = 2, mintime = 0)
    
    #Exclude first and last time points
    out_aux <- outer(times_i_l[2:(n_times_i_l-1)], H_i_l_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    H_i[[l]] <- H_i_l_msm$states[c(1,ind-1,length(H_i_l_msm$times))]  - 1 #(0,1)
    
    #For estimation of the initial transition rates we do not need any special q matrix
    #It only has to indicate id there are some constrained zeros in the matrix to be estimated (not for us I think)
    #Problems arise when there are no or only one transitions (zeros)
    data_H_i_l_df <- data.frame(states = H_i[[l]]+1, time = times_i_l)
    lambda_aux <- crudeinits.msm(states ~ time, data = data_H_i_l_df, qmatrix = Q_H)
    lambda_aux[lambda_aux == 0] <- 1
    lambda_aux_vec <- lambda_aux[1,1]
    for(i2 in 1:dH[l]){
      lambda_aux_vec <- c(lambda_aux_vec, lambda_aux[i2,c(1:dH[l])[-i2]])
    }
    lambda_aux_vec <- lambda_aux_vec[-1]
    lambda_H_i_hat[[l]] <- lambda_aux_vec
  }
  H[[i]] <- H_i
  lambda_H_hat[[i]] <- lambda_H_i_hat
}

#Put the estimated lambdas in a list
lambda_hat_all <- matrix(0,N,p_tot)
for(i in 1:N){
  lambda_hat_i <- 0
  for(h in 1:pY){
    lambda_hat_i <- c(lambda_hat_i, lambda_Y_hat[[i]][[h]])
  }
  for(l in 1:pH){
    lambda_hat_i <- c(lambda_hat_i, lambda_H_hat[[i]][[l]])
  }
  lambda_hat_all[i,] <- lambda_hat_i[2:(p_tot+1)]
}


#Interval lengths
epsilon <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  epsilon_i <- vector("list", length = p0)
  for(ip in 1:p0){
    epsilon_i[[ip]] <- diff(times_i[[ip]])
  }
  epsilon[[i]] <- epsilon_i
}

data <- list(Y = Y, H = H, lambda_hat_all = lambda_hat_all, n_states = n_states, n_rates = n_rates, n_times = n_times, epsilon = epsilon, X = X, Z = Z)

## Now that we have the data, we need to initialize the parameters and priors
n_burn1 <- 10
n_burn2 <- 10
n_save <- 100
thin <- 1
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8

SM_alg <- TRUE

Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold, SM_alg = SM_alg, t_SM_split = 3, t_SM_merge = 3)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

#A-priori expected number of components
Mm <- 5

#Prior on Lambda
update_Lambda <- TRUE
if(update_Lambda){# We provide hyperparameters for Lambda
  #These are based on mean and var for M (number of components)
  mm <- 1
  ss <- 1 #Variance quite small?
  a2 <- mm^2/ss
  b2 <- mm/ss
  Lambda <- mm
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda, a2 = a2, b2 = b2)
}else{# We fix the values of Lambda
  Lambda <- 1
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda)
}
Param_list <- modifyList(Lambda_list,Param_list)

update_gamma_S <- TRUE
if(update_gamma_S){# We provide hyperparameters for gamma_S
  mm <- 1
  ss <- 1
  a1 <- mm^2/ss
  b1 <- mm/ss
  gamma_S <- mm
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S, a1 = a1, b1 = b1)
}else{# We fix the values of gamma_S
  gamma_S <- 1
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S)
}
Param_list <- modifyList(gamma_S_list,Param_list)

#Prior on mu
update_mu <- TRUE
if(update_mu){# We provide hyperparameters for mu
  m_mu <- rep(0,p_tot)
  mu <- m_mu
  mu_list <- list(update_mu = update_mu, mu = mu, m_mu = m_mu)
}else{# We fix the values of mu
  mu <- mu_simul
  mu_list <- list(update_mu = update_mu, mu = mu)
}
Param_list <- modifyList(mu_list,Param_list)

#Prior on m_mu
update_m_mu <- FALSE
if(update_mu){# We provide hyperparameters for m_mu
  Sigma_m_mu <- diag(p_tot)
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu, Sigma_m_mu = Sigma_m_mu)
}else{# We fix the values of m_mu
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu)
}
Param_list <- modifyList(m_mu_list,Param_list)

#Prior on k0
update_k0 <- TRUE
if(update_k0){# We provide hyperparameters for k0
  mm <- 1
  ss <- 1
  a_k0 <- mm^2/ss
  b_k0 <- mm/ss
  k0 <- mm
  k0_list <- list(update_k0 = update_k0, k0 = k0, a_k0 = a_k0, b_k0 = b_k0)
}else{# We fix the values of k0
  k0 <- k0_simul
  k0_list <- list(update_k0 = update_k0, k0 = k0)
}
Param_list <- modifyList(k0_list,Param_list)

#Prior on graph
size_based_prior <- TRUE
if(size_based_prior){
  #The name of the parameters is the same as for d, because it's like marginalizeing wrt d
  a_d <- 1
  b_d <- 1
  G_list <- list(size_based_prior = size_based_prior, a_d = a_d, b_d = b_d)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_d <- TRUE
  if(update_d){
    #Different prior options for d
    d_beta_prior <- FALSE
    if(d_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_d <- mm*(mm*(1 - mm)/ss - 1)
      b_d <- (mm*(1 - mm)^2/ss - (1 - mm))
      d <- mm
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_d = a_d, b_d = b_d, d = d)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      d <- (1 + exp(-mm))^(-1)
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_lambda = a_lambda, d = d)
    }
  }else{# We fix the values of d
    d <- 0.25
    d_list <- list(update_d = update_d, d = d)
  }
  Param_list <- modifyList(d_list,Param_list)
}
Param_list <- modifyList(G_list,Param_list)


#Prior on beta
mu_beta <- vector("list", length = pY)
U_beta <- vector("list", length = pY)
V_beta <- vector("list", length = pY)
for(h in 1:pY){
  mu_beta[[h]] <- matrix(0,g[h],dY[h]*(dY[h] - 1))
  U_beta[[h]] <- diag(g[h])
  V_beta[[h]] <- diag(dY[h]*(dY[h] - 1))
}
beta_list <- list(mu_beta = mu_beta, U_beta = U_beta, V_beta = V_beta)
Param_list <- modifyList(beta_list, Param_list)

#Prior on gamma
mu_gamma <- vector("list", length = pY)
U_gamma <- vector("list", length = pY)
V_gamma <- vector("list", length = pY)
for(h in 1:pY){
  mu_gamma[[h]] <- matrix(0,q[h],dY[h]*(dY[h] - 1))
  U_gamma[[h]] <- diag(q[h])
  V_gamma[[h]] <- diag(dY[h]*(dY[h] - 1))
}
gamma_list <- list(mu_gamma = mu_gamma, U_gamma = U_gamma, V_gamma = V_gamma)
Param_list <- modifyList(gamma_list, Param_list)


is_DP <- FALSE
is_parametric <- FALSE
if(is_parametric){#We can fit a model for a given clustering
  c_init <- c(1:N)-1
  alpha_list <- list(update_alpha = FALSE, alpha = 0)
}else{
  if(is_DP){#We can fit a DP model
    #Prior on alpha
    # #Based on prior expected number of components
    alpha <- matrix(seq(0.001,2,length = 50000),50000,1)
    sum_alpha <- apply(alpha,1,function(alpha) sum(alpha/(alpha + c(1:N) - 1)))
    alpha_Mm <- alpha[which.min(abs(Mm - sum_alpha))]
    sum_alpha[which.min(abs(Mm - sum_alpha))]
    
    update_alpha <- TRUE
    if(update_alpha){# We provide hyperparameters for alpha
      mm <- alpha_Mm
      ss <- 1
      a_alpha <- mm^2/ss
      b_alpha <- mm/ss
      alpha <- mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha, a_alpha = a_alpha, b_alpha = b_alpha)
    }else{# We fix the values of alpha
      alpha <- alpha_Mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha)
    }
  }else{
    alpha_list <- list(update_alpha = FALSE, alpha = 0)
  }
  c_init <- rep(1,N)-1
}
Param_list <- modifyList(alpha_list,Param_list)

model_list <- list(is_parametric = is_parametric, is_DP = is_DP , c_init = c_init)

# Main Gibbs sampler
sample.bdmcmc <- DisProgr( data = data, model_list = model_list, Alg_list = Alg_list, Param_list = Param_list)




#Produce plots for report

Graph_list <- sample.bdmcmc$Graph_List

Omega_out <- Graph_list$Omega_out
G_out <- Graph_list$G_out
G0_out <- Graph_list$G0_out

Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

layout(matrix(c(1,2),nrow = 1, ncol = 2))

col_val <- c(Omega_mean, Omega_simul)
col_map <- colorRampPalette(c("white", "red"))(99)

image(Omega_mean, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_mean, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))
image(Omega_simul, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_simul, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))

image(G_mean)
image(G_simul)

image(G0_mean)
image(G0_simul)



sum_weights_out <- 1/Graph_list$sum_weights_out
Omega_meanW <- matrix(0,p_tot,p_tot)
G_meanW <- matrix(0,p_tot,p_tot)
G0_meanW <- matrix(0,p0,p0)
for(it in 1:n_save){
  Omega_meanW = Omega_meanW + Omega_out[,,it] * sum_weights_out[it]
  G_meanW = G_meanW + G_out[,,it] * sum_weights_out[it]
  G0_meanW = G0_meanW + G0_out[,,it] * sum_weights_out[it]
}
Omega_meanW <- Omega_meanW / sum(sum_weights_out)
G_meanW <- G_meanW / sum(sum_weights_out)
G0_meanW <- G0_meanW / sum(sum_weights_out)

image(Omega_meanW)
image(Omega_simul)

image(G_meanW)
image(G_simul)

image(G0_meanW)
image(G0_simul)

dev.off()

#d
d_out <- Graph_list$d_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20)
abline(v = d_simul, lwd = 3, col = "red")
plot(d_out, type = "l")

#mu
mu_out <- sample.bdmcmc$mu_out
matplot(t(mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~mu), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#m_mu
m_mu_out <- sample.bdmcmc$m_mu_out
matplot(t(m_mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~m[mu]), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(m_mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), m_mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(m_mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#k0
k0_out <- sample.bdmcmc$k0_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), breaks = 20)
abline(v = k0_simul, lwd = 3, col = "red")
plot(k0_out, type = "l")



#Phi
phi_star_out <- sample.bdmcmc$phi_star_out
c_out <- sample.bdmcmc$c_out
phi_mean <- matrix(0, N, p_tot)
for(it in 1:n_save){
  phi_mean <- phi_mean + phi_star_out[[it]][c_out[it,]+1,]
}
phi_mean <- phi_mean/n_save
matplot(c(1:p_tot), t(phi_mean), col = "grey", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~phi), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
matplot(c(1:p_tot), t(phi_simul), col = "red", lwd = 1.5, lty = 1, type = "l", add = TRUE)
legend("topright", legend = c("MCMC", "truth"), pch = 19, col = c("grey", "red"), bty = "n", cex = 2)


#Number of clusters
K_N_out <- rep(0,n_save)
for(it in 1:n_save){
  K_N_out[it] <- dim(phi_star_out[[it]])[1]
}
plot(table(K_N_out)/n_save, col = "blue", lwd = 2)


#beta
beta_out <- sample.bdmcmc$beta_out
for(h in 1:pY){
  beta_out_h <- array(0, dim = c(g[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    beta_out_h[,,it] <- beta_out[[it]][[h]]
  }
  
  beta_df <- data.frame(dim_g = rep(c(1:(g[h]*dY[h]*(dY[h] - 1))), n_save), beta_vec = c(beta_out_h))
  beta_true_df <- data.frame( x = c(1:(g[h]*dY[h]*(dY[h] - 1))), y = c(beta_simul[[h]]) ) 
  
  ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
    geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
    labs(x = "g x dY", y = "", title = bquote("Posterior of "~beta)) +
    theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
}


#gamma
gamma_out <- sample.bdmcmc$gamma_out
for(h in 1:pY){
  gamma_out_h <- array(0, dim = c(q[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    gamma_out_h[,,it] <- gamma_out[[it]][[h]]
  }
  
  gamma_df <- data.frame(dim_g = rep(c(1:(q[h]*dY[h]*(dY[h] - 1))), n_save), gamma_vec = c(gamma_out_h))
  gamma_true_df <- data.frame( x = c(1:(q[h]*dY[h]*(dY[h] - 1))), y = c(gamma_simul[[h]]) ) 
  
  ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
    geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
    labs(x = "g x dY", y = "", title = bquote("Posterior of "~gamma)) +
    theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
}



#Lambda
Lambda_out <- sample.bdmcmc$Lambda_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(Lambda_out, col = "lightblue", main = bquote("Posterior of "~Lambda), breaks = 20)
plot(Lambda_out, type = "l")


#gamma_S
gamma_S_out <- sample.bdmcmc$gamma_S_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(gamma_S_out, col = "lightblue", main = bquote("Posterior of "~gamma[S]), breaks = 20)
plot(gamma_S_out, type = "l")



#To see how the transition rates are estimated I want to show the evolution of trajectories

phi_star_out <- sample.bdmcmc$phi_star_out
c_out <- sample.bdmcmc$c_out
beta_out <- sample.bdmcmc$beta_out
gamma_out <- sample.bdmcmc$gamma_out

#Reconstruct transition rates for plots
#It's parametric so phi_star_out is a matrix
lambda_Y_out <- vector("list", length = n_save)
lambda_H_out <- vector("list", length = n_save)
for(it in 1:n_save){
  phi_star_it_out <- phi_star_out[[it]]
  
  beta_out_it <- beta_out[[it]]
  gamma_out_it <- gamma_out[[it]]
  
  lambda_Yi_out <- vector("list", length = N)
  lambda_Hi_out <- vector("list", length = N)
  for(i in 1:N){
    n_times_i <- n_times[i,]
    lambda_Y_i <- vector("list", length = pY)
    for(h in 1:pY){
      ind_dY <- (n_rates_cum[h]+1):n_rates_cum[h+1]
      
      n_times_i_h <- n_times_i[h]
      lambda_Y_i_h <- matrix(0, n_times_i_h, dY[h]*(dY[h]-1))
      
      X_h <- X[[h]]
      Z_h <- Z[[h]]
      Z_h_i <- Z_h[[i]]
      beta_it_h_out <- beta_out_it[[h]]
      gamma_it_h_out <- gamma_out_it[[h]]
      
      lambda_Y_i[[h]] <- matrix(exp( phi_star_it_out[c_out[it,i]+1,ind_dY] + X_h[i,] %*% beta_it_h_out), nrow = n_times_i[h]-1, ncol = dY[h]*(dY[h]-1), byrow = TRUE) * exp( Z_h_i %*% gamma_it_h_out )
    }
    lambda_H_i <- vector("list", length = pH)
    for(l in 1:pH){
      ind_dH <- (n_rates_cum[pY+l]+1):n_rates_cum[pY+l+1]
      lambda_H_i[[l]] <- exp(phi_star_it_out[c_out[it,i]+1,ind_dH])
    }
    lambda_Yi_out[[i]] <- lambda_Y_i
    lambda_Hi_out[[i]] <- lambda_H_i
  }
  lambda_Y_out[[it]] <- lambda_Yi_out
  lambda_H_out[[it]] <- lambda_Hi_out
}

#Reconstruct simulated transition rates
lambda_Y_simul <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i,]
  lambda_Y_simul_i <- vector("list", length = pY)
  for(h in 1:pY){
    n_times_i_h <- n_times_i[h]
    X_h <- X[[h]]
    Z_h <- Z[[h]]
    Z_h_i <- Z_h[[i]]
    beta_simul_h <- beta_simul[[h]]
    gamma_simul_h <- gamma_simul[[h]]
    lambda_Y_simul_i[[h]] <- matrix(lambda_Y_tilde_simul[[c_simul[i]]][[h]] * exp( X_h[i,] %*% beta_simul_h ), n_times_i_h - 1, dY[h] * (dY[h] - 1), byrow = TRUE) * exp( Z_h_i %*% gamma_simul_h )
  }
  lambda_Y_simul[[i]] <- lambda_Y_simul_i
}


#Select patient
i_plot <- 28
Y_i <- Y[[i_plot]]
H_i <- H[[i_plot]]
times_i_plot <- times[[i_plot]]
n_times_i_plot <- n_times[i_plot,]


#Simulated transition probabilities
prs_simul_log <- vector("list", length(p0))

lambda_Y_simul_i_plot <- lambda_Y_simul[[i_plot]]
epsilon_i_plot <- epsilon[[i_plot]]
for(h in 1:pY){
  n_times_i_h_plot <- n_times_i_plot[h]
  prs_simul_log_h <- rep(0,n_times_i_h_plot)
  Y_i_h <- Y_i[[h]]
  lambda_Y_simul_i_h_plot <- lambda_Y_simul_i_plot[[h]]
  epsilon_i_h_plot <- epsilon_i_plot[[h]]
  for(j in 2:n_times_i_h_plot){
    Yihj <- Y_i_h[j-1]
    prs_simul_log_h[j] <- log(lambda_Y_simul_i_h_plot[j-1,Yihj+1]) - log(sum(lambda_Y_simul_i_h_plot[j-1,])) + log( 1 - exp(- ((sum(lambda_Y_simul_i_h_plot[j-1,])) * epsilon_i_h_plot[j-1])) )
  }
  prs_simul_log[[h]] <- prs_simul_log_h
}

lambda_H_simul_i_plot <- lambda_H_simul[[c_simul[i_plot]]]
for(l in 1:pH){
  n_times_i_l_plot <- n_times_i_plot[pY+l]
  prs_simul_log_l <- rep(0,n_times_i_l_plot)
  H_i_l <- H_i[[l]]
  lambda_H_simul_i_l_plot <- lambda_H_simul_i_plot[[l]]
  epsilon_i_l_plot <- epsilon_i_plot[[pY+l]]
  for(j in 2:n_times_i_l_plot){
    Hilj <- H_i_l[j-1]
    prs_simul_log_l[j] <- log(lambda_H_simul_i_l_plot[Hilj+1]) - log(sum(lambda_H_simul_i_l_plot)) + log( 1 - exp(- (sum(lambda_H_simul_i_l_plot) * epsilon_i_l_plot[j-1])) )
  }
  prs_simul_log[[pY+l]] <- prs_simul_log_l
}

#Estimated transition probabilities
prs_out_log <- vector("list", length = p0)

for(h in 1:pY){
  n_times_i_h_plot <- n_times_i_plot[h]
  prs_out_log_h <- matrix(0,n_times_i_h_plot,n_save)
  Y_i_h <- Y_i[[h]]
  epsilon_i_h_plot <- epsilon_i_plot[[h]]
  for(j in 2:n_times_i_h_plot){
    Yihj <- Y_i_h[j-1]
    for(it in 1:n_save){
      lambda_Y_it_i_plot_h_out <- lambda_Y_out[[it]][[i_plot]][[h]]
      prs_out_log_h[j,it] <- log(lambda_Y_it_i_plot_h_out[j-1,Yihj+1]) - log(sum(lambda_Y_it_i_plot_h_out[j-1,])) + log( 1 - exp(- (sum(lambda_Y_it_i_plot_h_out[j-1,]) * epsilon_i_h_plot[j-1])) )
    }
  }
  prs_out_log[[h]] <- prs_out_log_h
}

for(l in 1:pH){
  n_times_i_l_plot <- n_times_i_plot[pY+l]
  prs_out_log_l <- matrix(0,n_times_i_l_plot,n_save)
  H_i_l <- H_i[[l]]
  epsilon_i_l_plot <- epsilon_i_plot[[pY+l]]
  for(j in 2:n_times_i_l_plot){
    Hilj <- H_i_l[j-1]
    for(it in 1:n_save){
      lambda_H_it_i_plot_l_out <- lambda_H_out[[it]][[i_plot]][[l]]
      prs_out_log_l[j,it] <- log(lambda_H_it_i_plot_l_out[Hilj+1]) - log(sum(lambda_H_it_i_plot_l_out)) + log( 1 - exp(- (sum(lambda_H_it_i_plot_l_out) * epsilon_i_l_plot[j-1])) )
    }
  }
  prs_out_log[[pY+l]] <- prs_out_log_l
}


#Colours of processes
col_proc <- rgb(255/255,51/255,51/255)

scale_red1 <- seq(0, 200, length = pY)
scale_red2 <- seq(0, 255, length = pY)
for(h in 1:pY){
  col_proc <- c(col_proc, rgb(255/255,scale_red1[h]/255,scale_red2[h]/255))
}
scale_blue1 <- seq(20, 200, length = pH)
scale_blue2 <- seq(120, 255, length = pH)
for(l in 1:pH){
  col_proc <- c(col_proc, rgb(scale_blue1[l]/255,scale_blue2[l]/255,255/255))
}
col_proc <- col_proc[-1]

#Labels of processes
proc_lab <- c(paste0("Y", c(1:pY)), paste0("H", c(1:pH)))

#Credibility intervals
prs_CI_low <- vector("list", length = p0)
prs_CI_up <- vector("list", length = p0)
for(ip in 1:p0){
  prs_CI_low[[ip]] <- apply(prs_out_log[[ip]],1, quantile, probs = 0.025)
  prs_CI_up[[ip]] <- apply(prs_out_log[[ip]],1, quantile, probs = 0.975)
}
prs_CI_df <- data.frame(prs_CI_low = unlist(prs_CI_low), prs_CI_up = unlist(prs_CI_up), proc = unlist(mapply(rep,proc_lab,n_times_i_plot), use.names = FALSE), times_vec = round(unlist(times_i_plot),3) )
prs_simul_df <- data.frame(x_times = unlist(times_i_plot), prs_simul_vec = unlist(prs_simul_log), proc = unlist(mapply(rep,proc_lab,n_times_i_plot), use.names = FALSE), times_vec = round(unlist(times_i_plot),3) )
prs_df <- data.frame(x_times = unlist(mapply(rep,times_i_plot,n_save)), prs_vec = unlist(prs_out_log), proc = unlist(mapply(rep,proc_lab,n_times_i_plot * n_save), use.names = FALSE), times_vec = round(unlist(mapply(rep,times_i_plot,n_save)),3) )


ggplot(prs_df, aes(x = x_times, y = prs_vec, col = proc)) + 
  geom_line(data = prs_simul_df, aes(x = x_times, y = prs_simul_vec, col = proc), size = 1.25, linetype = 2) + 
  stat_summary(fun.y = 'mean', geom = 'line', size = 1.25) +
  stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2) + 
  labs(x = "time", y = bquote(log(p[t]^{rs})), title = bquote("Posterior of "~log(p[t]^{rs})~" for subject "~.(i_plot))) +
  theme(text = element_text(size=35), axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 15)) +
  facet_grid(rows = vars(proc), scales="free") +
  # the labels must match what you specified above
  scale_fill_manual(name = proc_lab, values = col_proc) +
  scale_color_manual(name = proc_lab, values = col_proc) +
  theme_bw()





#####Binder loss function equal costs
K_Binder <- 1/2
#For selected combinations of parameters
c_out <- sample.bdmcmc$c_out

pij <- matrix(0,N,N)
cij <- vector("list", length = n_save)
for(g in 1:n_save){
  cij[[g]] <- outer(c_out[g,], c_out[g,], "==")
  pij <- pij + cij[[g]]
}
pij <- pij/n_save


Binder_f <- rep(0,n_save)
for(g in 1:n_save){
  aux <- (pij - K_Binder) * as.matrix(cij[[g]])
  aux <-  aux[upper.tri(aux)]
  Binder_f[g] <- sum(aux)
}
Binder_ind <- which.max(Binder_f)
Binder_out <- c_out[Binder_ind,] + 1
  
graphics.off()
library("gplots")
Binder_sort <- sort(Binder_out, index.return = TRUE)
Binder_sort <- Binder_sort$ix

heatmap.2(pij[Binder_sort,Binder_sort], ColSideColor = rainbow(K_N_simul)[c_simul[Binder_sort]], col=rev(heat.colors(100)), trace = "none", key  = TRUE, Rowv = FALSE, Colv = FALSE, dendrogram = "none", main = "Posterior co-clustering probabilities", cex.main = 3, cex.lab = 2, cex.axis = 1.5)











#################################################################################################################
# dH = 0
#################################################################################################################

## Here we do a fuller simulation from a parametric model where we also have the Markov processes
## We will try to handle directly different time points for different processes and patients

rm(list = ls())
set.seed(123)
library("BayesDiseaseProgr")
library("ggplot2")
library("msm")
library("fields")


######
#Simulate sensible data
######

#Number of patients
N = 50;
#Number of response processes
pY = 5;
#Number of covariate processes
pH = 0;
#Tot number of porcesses
p0 = pY+pH;

#Set graph structure(s)
# d_simul <- 0.25
G0_simul <- matrix(0, p0, p0)
# G0_simul[lower.tri(G0_simul)] <- c(runif(p0*(p0-1)/2) <= d_simul)
G0_simul[1,2] <- 1
G0_simul[1,4] <- 1
G0_simul[2,3] <- 1
G0_simul[3,4] <- 1
G0_simul[4,5] <- 1
G0_simul[2,5] <- 1

G0_simul <- G0_simul + t(G0_simul)

#Assume 2-state responses
dY <- c(2, 2, 3, 4, 5)
dH <- NULL
n_states <- list(dY,dH)
n_rates <- 0
for(h in 1:pY){
  n_rates <- c(n_rates, dY[h] * (dY[h]-1))
}
n_rates_cum <- cumsum(n_rates)
p_tot <- sum(n_rates)
G_simul <- matrix(0, nrow = p_tot, ncol = p_tot)
diag(G0_simul) <- rep(1,p0)
for(i in 1:p0){
  for(j in i:p0){
    if(G0_simul[i,j]){
      G_simul[c((n_rates_cum[i]+1) : n_rates_cum[i+1]), c((n_rates_cum[j]+1) : n_rates_cum[j+1])] <- 1
      G_simul[c((n_rates_cum[j]+1) : n_rates_cum[j+1]), c((n_rates_cum[i]+1) : n_rates_cum[i+1])] <- 1 #For symmetry
    }
  }
}
diag(G0_simul) <- rep(0,p0)
diag(G_simul) <- rep(0,p_tot)

#Prior specification
nu_simul <- 5
Psi_simul <- diag(p_tot)

#Simulate data
Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

phi_simul <- rmvnorm(n = N, mean = c(mu_simul), sigma = solve(Omega_simul))
lambda_simul <- exp(phi_simul)

#Number of visits per patient per process (for msm package)
n_times <- matrix(9 + rpois(N*p0, 1), N, p0)
#Times of visit (0 and Tmax included!)
T_max <- 30
times <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i,]
  times_i <- list()
  for(ip in 1:p0){
    times_i[[ip]] <- sort(c(0,runif(n_times_i[ip]-2, min = 0, max = T_max), T_max))
  }
  times[[i]] <- times_i
}

#Transition rates for Y are covariate- and time-dependent
#Transition rates for processes H are constant (and H are binary processes)
lambda_Y_tilde_simul <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_tilde_simul_i <- list()
  for(h in 1:pY){
    lambda_Y_tilde_simul_i[[h]] <- lambda_simul[i,(n_rates_cum[h]+1):n_rates_cum[h+1]]
  }
  lambda_Y_tilde_simul[[i]] <- lambda_Y_tilde_simul_i
}

#Constant covariates, but they can differ in each process i,...,pY
g <- c(2, 4, 3, 2, 1) #Simulate age and sex without intercept (it's the lambda_tilde!)
X <- list()

aux <- matrix(1,N,g[1])
aux[,1] <- rnorm(N, mean = 0, sd = 1) #standardized
aux[,2] <- (runif(N) <= 0.5)
X[[1]] <- aux

X[[2]] <- matrix(rgamma(N*g[2], 1, 1), N, g[2])

X[[3]] <- cbind(runif(N) <= 0.25, rnorm(N, mean = -1, sd = 1), runif(N) <= 0.1)

aux <- matrix(1,N,g[1])
aux[,1] <- rnorm(N, mean = 2, sd = 1) #standardized
aux[,2] <- (runif(N) > 0.3)
X[[4]] <- aux
  
X[[5]] <- matrix(runif(N) <= 0.15)


#Coefficient beta (different for each rate and process 1,...,pY)
beta_simul <- vector("list", length = pY)
for(h in 1:pY){
  beta_simul[[h]] <- matrix(rnorm(g[h]*dY[h]*(dY[h] - 1), sd = 0.1),g[h],dY[h]*(dY[h] - 1),byrow = TRUE)
}

#Time-varying covariates (not random for msm package, no need for averaging)
#Also depends on the processes
#Here we put one equal to zero
q <- c(2,3,2,2,2)
Z_msm <- vector("list", length = pY)
Z <- vector("list", length = pY)

#Careful with dimension q!
Z_msm_h <- vector("list", length = N)
Z_h <- vector("list", length = N)
h = 1
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm(times_i/T_max, sd = sqrt(0.5))
  Z2 <- rnorm(cos(times_i/T_max*2*pi), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h

h = 2
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rep(0,n_times_i)
  Z2 <- rep(0,n_times_i)
  Z3 <- rep(0,n_times_i)
  Z_msm_h[[i]] <- cbind(Z1,Z2,Z3)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z3 <- (Z3[1:(n_times_i-1)] + Z3[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2,Z3)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h

h = 3
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rnorm((times_i/T_max)^2, sd = sqrt(0.5))
  Z2 <- rnorm((times_i/T_max)^3, sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h


h = 4
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rnorm((times_i/T_max)^(.2), sd = sqrt(0.5))
  Z2 <- rnorm((times_i/T_max)^(.3), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h


h = 5
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rnorm(exp(times_i/T_max), sd = sqrt(0.5))
  Z2 <- rnorm(sin(times_i/T_max), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h



#Coefficient gamma (different for each rate and process 1,...,pY)
gamma_simul <- vector("list", length = pY)
for(h in 1:pY){
  gamma_simul[[h]] <- matrix(rnorm(q[h]*dY[h]*(dY[h] - 1), sd = 0.1),q[h],dY[h]*(dY[h] - 1),byrow = TRUE)
}

#Constant part of transition rates (needed in sim.msm, the time-varying covariates are included in the function used for the simulation)
lambda_Y_msm <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_msm_i <- list()
  for(h in 1:pY){
    lambda_Y_msm_i[[h]] <- lambda_Y_tilde_simul[[i]][[h]] * exp( X[[h]][i,] %*% beta_simul[[h]] )
  }
  lambda_Y_msm[[i]] <- lambda_Y_msm_i
}

#Generate the data (...a bit tricky with msm...
#...but we can do it independently and at different times for each process and patient...)

#Response process
Y <- vector("list", length = N)
lambda_Y_hat <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i,]
  
  lambda_Y_msm_i <- lambda_Y_msm[[i]]
  
  Y_i <- list()
  lambda_Y_hat_i <- list()
  for(h in 1:pY){
    
    lambda_Y_msm_i_h <- lambda_Y_msm_i[[h]]
    Q_Y <- matrix(0, dY[h], dY[h])
    for(i2 in 1:dY[h]){
      Q_Y[i2,c(1:dY[h])[-i2]] <- lambda_Y_msm_i_h[(i2-1)*(dY[h]-1) + c(1:(dY[h]-1))]
      Q_Y[i2,i2] <- -sum(Q_Y[i2,])
    }
    
    #Here we can include the time-varying covariates via sim.msm
    times_i_h <- times_i[[h]]
    n_times_i_h <- n_times_i[h]
    
    Z_msm_i_h <- Z_msm[[h]][[i]]
    gamma_simul_h <- gamma_simul[[h]]
    
    Y_i_h_msm <- sim.msm(qmatrix = Q_Y, maxtime = T_max, covs = Z_msm_i_h, beta = gamma_simul_h, obstimes = times_i_h, start = 2, mintime = 0)
    # Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
    
    out_aux <- outer(times_i_h[2:(n_times_i_h-1)], Y_i_h_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    Y_i[[h]] <- Y_i_h_msm$states[c(1,ind-1,length(Y_i_h_msm$times))] - 1 #(0,1)
    
    #For estimation of the initial transition rates we do not need any special q matrix
    #It only has to indicate id there are some constrained zeros in the matrix to be estimated (not for us I think)
    #Problems arise when there are no or only one transitions (zeros)
    data_Y_i_df <- data.frame(states = Y_i[[h]]+1, time = times_i_h)
    lambda_aux <- crudeinits.msm(states ~ time, data = data_Y_i_df, qmatrix = Q_Y)
    #If zeroes we impose arbitrary value
    lambda_aux[lambda_aux == 0] <- 0.5
    lambda_aux_vec <- lambda_aux[1,1]
    for(i2 in 1:dY[h]){
      lambda_aux_vec <- c(lambda_aux_vec, lambda_aux[i2,c(1:dY[h])[-i2]])
    }
    lambda_aux_vec <- lambda_aux_vec[-1]
    lambda_Y_hat_i[[h]] <- lambda_aux_vec
  }
  Y[[i]] <- Y_i
  lambda_Y_hat[[i]] <- lambda_Y_hat_i
}


#Put the estimated lambdas in a list
lambda_hat_all <- matrix(0,N,p_tot)
for(i in 1:N){
  lambda_hat_i <- 0
  for(h in 1:pY){
    lambda_hat_i <- c(lambda_hat_i, lambda_Y_hat[[i]][[h]])
  }
  lambda_hat_all[i,] <- lambda_hat_i[2:(p_tot+1)]
}


#Interval lengths
epsilon <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  epsilon_i <- vector("list", length = p0)
  for(ip in 1:p0){
    epsilon_i[[ip]] <- diff(times_i[[ip]])
  }
  epsilon[[i]] <- epsilon_i
}

data <- list(Y = Y, lambda_hat_all = lambda_hat_all, n_states = n_states, n_rates = n_rates, n_times = n_times, epsilon = epsilon, X = X, Z = Z)

## Now that we have the data, we need to initialize the parameters and priors
n_burn1 <- 10
n_burn2 <- 10
n_save <- 100
thin <- 1
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8

Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

#Prior on mu
update_mu <- TRUE
if(update_mu){# We provide hyperparameters for mu
  m_mu <- rep(0,p_tot)
  mu <- m_mu
  mu_list <- list(update_mu = update_mu, mu = mu, m_mu = m_mu)
}else{# We fix the values of mu
  mu <- mu_simul
  mu_list <- list(update_mu = update_mu, mu = mu)
}
Param_list <- modifyList(mu_list,Param_list)

#Prior on m_mu
update_m_mu <- FALSE
if(update_mu){# We provide hyperparameters for m_mu
  Sigma_m_mu <- diag(p_tot)
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu, Sigma_m_mu = Sigma_m_mu)
}else{# We fix the values of m_mu
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu)
}
Param_list <- modifyList(m_mu_list,Param_list)

#Prior on k0
update_k0 <- TRUE
if(update_k0){# We provide hyperparameters for k0
  mm <- 1
  ss <- 1
  a_k0 <- mm^2/ss
  b_k0 <- mm/ss
  k0 <- mm
  k0_list <- list(update_k0 = update_k0, k0 = k0, a_k0 = a_k0, b_k0 = b_k0)
}else{# We fix the values of k0
  k0 <- k0_simul
  k0_list <- list(update_k0 = update_k0, k0 = k0)
}
Param_list <- modifyList(k0_list,Param_list)

#Prior on graph
size_based_prior <- TRUE
if(size_based_prior){
  #The name of the parameters is the same as for d, because it's like marginalizeing wrt d
  a_d <- 1
  b_d <- 1
  G_list <- list(size_based_prior = size_based_prior, a_d = a_d, b_d = b_d)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_d <- TRUE
  if(update_d){
    #Different prior options for d
    d_beta_prior <- TRUE
    if(d_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_d <- mm*(mm*(1 - mm)/ss - 1)
      b_d <- (mm*(1 - mm)^2/ss - (1 - mm))
      d <- mm
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_d = a_d, b_d = b_d, d = d)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      d <- (1 + exp(-mm))^(-1)
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_lambda = a_lambda, d = d)
    }
  }else{# We fix the values of d
    d <- 0.25
    d_list <- list(update_d = update_d, d = d)
  }
  Param_list <- modifyList(d_list,Param_list)
}
Param_list <- modifyList(G_list,Param_list)


#Prior on beta
mu_beta <- vector("list", length = pY)
U_beta <- vector("list", length = pY)
V_beta <- vector("list", length = pY)
for(h in 1:pY){
  mu_beta[[h]] <- matrix(0,g[h],dY[h]*(dY[h] - 1))
  U_beta[[h]] <- diag(g[h])
  V_beta[[h]] <- diag(dY[h]*(dY[h] - 1))
}
beta_list <- list(mu_beta = mu_beta, U_beta = U_beta, V_beta = V_beta)
Param_list <- modifyList(beta_list, Param_list)

#Prior on gamma
mu_gamma <- vector("list", length = pY)
U_gamma <- vector("list", length = pY)
V_gamma <- vector("list", length = pY)
for(h in 1:pY){
  mu_gamma[[h]] <- matrix(0,q[h],dY[h]*(dY[h] - 1))
  U_gamma[[h]] <- diag(q[h])
  V_gamma[[h]] <- diag(dY[h]*(dY[h] - 1))
}
gamma_list <- list(mu_gamma = mu_gamma, U_gamma = U_gamma, V_gamma = V_gamma)
Param_list <- modifyList(gamma_list, Param_list)

is_DP <- FALSE
is_parametric <- TRUE
if(is_parametric){#We can fit a model for a given clustering
  c_init <- c(1:N)-1
  alpha_list <- list(update_alpha = FALSE, alpha = 0)
}else{
  if(is_DP){#We can fit a DP model
    #Prior on alpha
    # #Based on prior expected number of components
    alpha <- matrix(seq(0.001,2,length = 50000),50000,1)
    sum_alpha <- apply(alpha,1,function(alpha) sum(alpha/(alpha + c(1:N) - 1)))
    alpha_Mm <- alpha[which.min(abs(Mm - sum_alpha))]
    sum_alpha[which.min(abs(Mm - sum_alpha))]
    
    update_alpha <- TRUE
    if(update_alpha){# We provide hyperparameters for alpha
      mm <- alpha_Mm
      ss <- 1
      a_alpha <- mm^2/ss
      b_alpha <- mm/ss
      alpha <- mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha, a_alpha = a_alpha, b_alpha = b_alpha)
    }else{# We fix the values of alpha
      alpha <- alpha_Mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha)
    }
  }else{
    alpha_list <- list(update_alpha = FALSE, alpha = 0)
  }
  c_init <- rep(1,N)-1
}
Param_list <- modifyList(alpha_list,Param_list)

model_list <- list(is_parametric = is_parametric, is_DP = is_DP , c_init = c_init)

# Main Gibbs sampler
sample.bdmcmc <- DisProgr( data = data, model_list = model_list, Alg_list = Alg_list, Param_list = Param_list)


# #Save output
# save.image("OUTPUT_DiseaseProgr_Parametric.RData")




#Produce plots for report

Graph_list <- sample.bdmcmc$Graph_List

Omega_out <- Graph_list$Omega_out
G_out <- Graph_list$G_out
G0_out <- Graph_list$G0_out

Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

layout(matrix(c(1,2),nrow = 1, ncol = 2))

col_val <- c(Omega_mean, Omega_simul)
col_map <- colorRampPalette(c("white", "red"))(99)

image(Omega_mean, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_mean, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))
image(Omega_simul, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_simul, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))

image(G_mean)
image(G_simul)

image(G0_mean)
image(G0_simul)



sum_weights_out <- 1/Graph_list$sum_weights_out
Omega_meanW <- matrix(0,p_tot,p_tot)
G_meanW <- matrix(0,p_tot,p_tot)
G0_meanW <- matrix(0,p0,p0)
for(it in 1:n_save){
  Omega_meanW = Omega_meanW + Omega_out[,,it] * sum_weights_out[it]
  G_meanW = G_meanW + G_out[,,it] * sum_weights_out[it]
  G0_meanW = G0_meanW + G0_out[,,it] * sum_weights_out[it]
}
Omega_meanW <- Omega_meanW / sum(sum_weights_out)
G_meanW <- G_meanW / sum(sum_weights_out)
G0_meanW <- G0_meanW / sum(sum_weights_out)

image(Omega_meanW)
image(Omega_simul)

image(G_meanW)
image(G_simul)

image(G0_meanW)
image(G0_simul)

dev.off()

#d
d_out <- Graph_list$d_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20)
abline(v = d_simul, lwd = 3, col = "red")
plot(d_out, type = "l")

#mu
mu_out <- sample.bdmcmc$mu_out
matplot(t(mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~mu), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#m_mu
m_mu_out <- sample.bdmcmc$m_mu_out
matplot(t(m_mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~m[mu]), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(m_mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), m_mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(m_mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#k0
k0_out <- sample.bdmcmc$k0_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), breaks = 20)
abline(v = k0_simul, lwd = 3, col = "red")
plot(k0_out, type = "l")



#Phi
phi_star_out <- sample.bdmcmc$phi_star_out
phi_mean <- matrix(0, N, p_tot)
for(it in 1:n_save){
  phi_mean <- phi_mean + phi_star_out[,,it]
}
phi_mean <- phi_mean/n_save
matplot(c(1:p_tot), t(phi_mean), col = "grey", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~phi), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
matplot(c(1:p_tot), t(phi_simul), col = "red", lwd = 1.5, lty = 1, type = "l", add = TRUE)
legend("topright", legend = c("MCMC", "truth"), pch = 19, col = c("grey", "red"), bty = "n", cex = 2)


#beta
beta_out <- sample.bdmcmc$beta_out
for(h in 1:pY){
  beta_out_h <- array(0, dim = c(g[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    beta_out_h[,,it] <- beta_out[[it]][[h]]
  }
  
  beta_df <- data.frame(dim_g = rep(c(1:(g[h]*dY[h]*(dY[h] - 1))), n_save), beta_vec = c(beta_out_h))
  beta_true_df <- data.frame( x = c(1:(g[h]*dY[h]*(dY[h] - 1))), y = c(beta_simul[[h]]) ) 
  
  ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
    geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
    labs(x = "g x dY", y = "", title = bquote("Posterior of "~beta)) +
    theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
}


#gamma
gamma_out <- sample.bdmcmc$gamma_out
for(h in 1:pY){
  gamma_out_h <- array(0, dim = c(q[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    gamma_out_h[,,it] <- gamma_out[[it]][[h]]
  }
  
  gamma_df <- data.frame(dim_g = rep(c(1:(q[h]*dY[h]*(dY[h] - 1))), n_save), gamma_vec = c(gamma_out_h))
  gamma_true_df <- data.frame( x = c(1:(q[h]*dY[h]*(dY[h] - 1))), y = c(gamma_simul[[h]]) ) 
  
  ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
    geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
    labs(x = "g x dY", y = "", title = bquote("Posterior of "~gamma)) +
    theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
}


#To see how the transition rates are estimated I want to show the evolution of trajectories

phi_star_out <- sample.bdmcmc$phi_star_out
beta_out <- sample.bdmcmc$beta_out
gamma_out <- sample.bdmcmc$gamma_out

#Reconstruct transition rates for plots
#It's parametric so phi_star_out is a matrix
lambda_Y_out <- vector("list", length = n_save)
for(it in 1:n_save){
  beta_out_it <- beta_out[[it]]
  gamma_out_it <- gamma_out[[it]]
  
  lambda_Yi_out <- vector("list", length = N)
  for(i in 1:N){
    n_times_i <- n_times[i,]
    lambda_Y_i <- vector("list", length = pY)
    for(h in 1:pY){
      ind_dY <- (n_rates_cum[h]+1):n_rates_cum[h+1]
      
      n_times_i_h <- n_times_i[h]
      lambda_Y_i_h <- matrix(0, n_times_i_h, dY[h]*(dY[h]-1))
      
      X_h <- X[[h]]
      Z_h <- Z[[h]]
      Z_h_i <- Z_h[[i]]
      beta_it_h_out <- beta_out_it[[h]]
      gamma_it_h_out <- gamma_out_it[[h]]
      
      lambda_Y_i[[h]] <- matrix(exp( phi_star_out[i,ind_dY,it] + X_h[i,] %*% beta_it_h_out), nrow = n_times_i[h]-1, ncol = dY[h]*(dY[h]-1), byrow = TRUE) * exp( Z_h_i %*% gamma_it_h_out )
    }
    lambda_Yi_out[[i]] <- lambda_Y_i
  }
  lambda_Y_out[[it]] <- lambda_Yi_out
}

#Reconstruct simulated transition rates
lambda_Y_simul <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i,]
  lambda_Y_simul_i <- vector("list", length = pY)
  for(h in 1:pY){
    n_times_i_h <- n_times_i[h]
    X_h <- X[[h]]
    Z_h <- Z[[h]]
    Z_h_i <- Z_h[[i]]
    beta_simul_h <- beta_simul[[h]]
    gamma_simul_h <- gamma_simul[[h]]
    lambda_Y_simul_i[[h]] <- matrix(lambda_Y_tilde_simul[[i]][[h]] * exp( X_h[i,] %*% beta_simul_h ), n_times_i_h - 1, dY[h] * (dY[h] - 1), byrow = TRUE) * exp( Z_h_i %*% gamma_simul_h )
  }
  lambda_Y_simul[[i]] <- lambda_Y_simul_i
}


#Select patient
i_plot <- 28
Y_i <- Y[[i_plot]]
times_i_plot <- times[[i_plot]]
n_times_i_plot <- n_times[i_plot,]


#Simulated transition probabilities
prs_simul_log <- vector("list", length(p0))

lambda_Y_simul_i_plot <- lambda_Y_simul[[i_plot]]
epsilon_i_plot <- epsilon[[i_plot]]
for(h in 1:pY){
  n_times_i_h_plot <- n_times_i_plot[h]
  prs_simul_log_h <- rep(0,n_times_i_h_plot)
  Y_i_h <- Y_i[[h]]
  lambda_Y_simul_i_h_plot <- lambda_Y_simul_i_plot[[h]]
  epsilon_i_h_plot <- epsilon_i_plot[[h]]
  for(j in 2:n_times_i_h_plot){
    Yihj <- Y_i_h[j-1]
    prs_simul_log_h[j] <- log(lambda_Y_simul_i_h_plot[j-1,Yihj+1]) - log(sum(lambda_Y_simul_i_h_plot[j-1,])) + log( 1 - exp(- ((sum(lambda_Y_simul_i_h_plot[j-1,])) * epsilon_i_h_plot[j-1])) )
  }
  prs_simul_log[[h]] <- prs_simul_log_h
}


#Estimated transition probabilities
prs_out_log <- vector("list", length = p0)

for(h in 1:pY){
  n_times_i_h_plot <- n_times_i_plot[h]
  prs_out_log_h <- matrix(0,n_times_i_h_plot,n_save)
  Y_i_h <- Y_i[[h]]
  epsilon_i_h_plot <- epsilon_i_plot[[h]]
  for(j in 2:n_times_i_h_plot){
    Yihj <- Y_i_h[j-1]
    for(it in 1:n_save){
      lambda_Y_it_i_plot_h_out <- lambda_Y_out[[it]][[i_plot]][[h]]
      prs_out_log_h[j,it] <- log(lambda_Y_it_i_plot_h_out[j-1,Yihj+1]) - log(sum(lambda_Y_it_i_plot_h_out[j-1,])) + log( 1 - exp(- (sum(lambda_Y_it_i_plot_h_out[j-1,]) * epsilon_i_h_plot[j-1])) )
    }
  }
  prs_out_log[[h]] <- prs_out_log_h
}



#Colours of processes
col_proc <- rgb(255/255,51/255,51/255)

scale_red1 <- seq(0, 120, length = pY)
scale_red2 <- seq(0, 255, length = pY)
for(h in 1:pY){
  col_proc <- c(col_proc, rgb(255/255,scale_red1[h]/255,scale_red2[h]/255))
}

#Labels of processes
proc_lab <- paste0("Y", c(1:pY))

#Credibility intervals
prs_CI_low <- vector("list", length = p0)
prs_CI_up <- vector("list", length = p0)
for(ip in 1:p0){
  prs_CI_low[[ip]] <- apply(prs_out_log[[ip]],1, quantile, probs = 0.025)
  prs_CI_up[[ip]] <- apply(prs_out_log[[ip]],1, quantile, probs = 0.975)
}
prs_CI_df <- data.frame(prs_CI_low = unlist(prs_CI_low), prs_CI_up = unlist(prs_CI_up), proc = unlist(mapply(rep,proc_lab,n_times_i_plot), use.names = FALSE), times_vec = round(unlist(times_i_plot),3) )
prs_simul_df <- data.frame(x_times = unlist(times_i_plot), prs_simul_vec = unlist(prs_simul_log), proc = unlist(mapply(rep,proc_lab,n_times_i_plot), use.names = FALSE), times_vec = round(unlist(times_i_plot),3) )
prs_df <- data.frame(x_times = unlist(mapply(rep,times_i_plot,n_save)), prs_vec = unlist(prs_out_log), proc = unlist(mapply(rep,proc_lab,n_times_i_plot * n_save), use.names = FALSE), times_vec = round(unlist(mapply(rep,times_i_plot,n_save)),3) )


ggplot(prs_df, aes(x = x_times, y = prs_vec, col = proc)) + 
  geom_line(data = prs_simul_df, aes(x = x_times, y = prs_simul_vec, col = proc), size = 1.25, linetype = 2) + 
  stat_summary(fun.y = 'mean', geom = 'line', size = 1.25) +
  stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2) + 
  labs(x = "time", y = bquote(log(p[t]^{rs})), title = bquote("Posterior of "~log(p[t]^{rs})~" for subject "~.(i_plot))) +
  theme(text = element_text(size=35), axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 15)) +
  facet_grid(rows = vars(proc), scales="free") +
  # the labels must match what you specified above
  scale_fill_manual(name = proc_lab, values = col_proc) +
  scale_color_manual(name = proc_lab, values = col_proc) +
  theme_bw()













#################################################################################################################
# dH = 0
#################################################################################################################

## Here we do a fuller simulation from a DP model where we also have the Markov processes
## We will try to handle directly different time points for different processes and patients

rm(list = ls())
set.seed(123)
library("BayesDiseaseProgr")
library("ggplot2")
library("msm")
library("fields")


######
#Simulate sensible data
######

#Number of patients
N = 50;
#Number of response processes
pY = 5;
#Number of covariate processes
pH = 0;
#Tot number of porcesses
p0 = pY+pH;

#Set graph structure(s)
# d_simul <- 0.25
G0_simul <- matrix(0, p0, p0)
# G0_simul[lower.tri(G0_simul)] <- c(runif(p0*(p0-1)/2) <= d_simul)
G0_simul[1,2] <- 1
G0_simul[1,4] <- 1
G0_simul[2,3] <- 1
G0_simul[3,4] <- 1
G0_simul[4,5] <- 1
G0_simul[2,5] <- 1

G0_simul <- G0_simul + t(G0_simul)

#Assume 2-state responses
dY <- c(2, 2, 3, 4, 5)
dH <- NULL
n_states <- list(dY,dH)
n_rates <- 0
for(h in 1:pY){
  n_rates <- c(n_rates, dY[h] * (dY[h]-1))
}
n_rates_cum <- cumsum(n_rates)
p_tot <- sum(n_rates)
G_simul <- matrix(0, nrow = p_tot, ncol = p_tot)
diag(G0_simul) <- rep(1,p0)
for(i in 1:p0){
  for(j in i:p0){
    if(G0_simul[i,j]){
      G_simul[c((n_rates_cum[i]+1) : n_rates_cum[i+1]), c((n_rates_cum[j]+1) : n_rates_cum[j+1])] <- 1
      G_simul[c((n_rates_cum[j]+1) : n_rates_cum[j+1]), c((n_rates_cum[i]+1) : n_rates_cum[i+1])] <- 1 #For symmetry
    }
  }
}
diag(G0_simul) <- rep(0,p0)
diag(G_simul) <- rep(0,p_tot)

#Prior specification
nu_simul <- 5
Psi_simul <- diag(p_tot)

#Simulate data
Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

K_N_simul <- 3
c_simul <- sample.int(K_N_simul, size = N, replace = TRUE, prob = c(1/3,1/6,1/2))
phi_simul <- rmvnorm(n = K_N_simul, mean = c(mu_simul), sigma = solve(Omega_simul))
lambda_simul <- exp(phi_simul)

#Number of visits per patient per process (for msm package)
n_times <- matrix(9 + rpois(N*p0, 1), N, p0)
#Times of visit (0 and Tmax included!)
T_max <- 30
times <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i,]
  times_i <- list()
  for(ip in 1:p0){
    times_i[[ip]] <- sort(c(0,runif(n_times_i[ip]-2, min = 0, max = T_max), T_max))
  }
  times[[i]] <- times_i
}

#Transition rates for Y are covariate- and time-dependent
lambda_Y_tilde_simul <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_tilde_simul_i <- list()
  for(h in 1:pY){
    lambda_Y_tilde_simul_i[[h]] <- lambda_simul[c_simul[i],(n_rates_cum[h]+1):n_rates_cum[h+1]]
  }
  lambda_Y_tilde_simul[[i]] <- lambda_Y_tilde_simul_i
}

#Constant covariates, but they can differ in each process i,...,pY
g <- c(2, 4, 3, 2, 1) #Simulate age and sex without intercept (it's the lambda_tilde!)
X <- list()

aux <- matrix(1,N,g[1])
aux[,1] <- rnorm(N, mean = 0, sd = 1) #standardized
aux[,2] <- (runif(N) <= 0.5)
X[[1]] <- aux

X[[2]] <- matrix(rgamma(N*g[2], 1, 1), N, g[2])

X[[3]] <- cbind(runif(N) <= 0.25, rnorm(N, mean = -1, sd = 1), runif(N) <= 0.1)

aux <- matrix(1,N,g[1])
aux[,1] <- rnorm(N, mean = 2, sd = 1) #standardized
aux[,2] <- (runif(N) > 0.3)
X[[4]] <- aux

X[[5]] <- matrix(runif(N) <= 0.15)


#Coefficient beta (different for each rate and process 1,...,pY)
beta_simul <- vector("list", length = pY)
for(h in 1:pY){
  beta_simul[[h]] <- matrix(rnorm(g[h]*dY[h]*(dY[h] - 1), sd = 0.1),g[h],dY[h]*(dY[h] - 1),byrow = TRUE)
}

#Time-varying covariates (not random for msm package, no need for averaging)
#Also depends on the processes
#Here we put one equal to zero
q <- c(2,3,2,2,2)
Z_msm <- vector("list", length = pY)
Z <- vector("list", length = pY)

#Careful with dimension q!
Z_msm_h <- vector("list", length = N)
Z_h <- vector("list", length = N)
h = 1
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm(times_i/T_max, sd = sqrt(0.5))
  Z2 <- rnorm(cos(times_i/T_max*2*pi), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h

h = 2
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rep(0,n_times_i)
  Z2 <- rep(0,n_times_i)
  Z3 <- rep(0,n_times_i)
  Z_msm_h[[i]] <- cbind(Z1,Z2,Z3)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z3 <- (Z3[1:(n_times_i-1)] + Z3[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2,Z3)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h

h = 3
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rnorm((times_i/T_max)^2, sd = sqrt(0.5))
  Z2 <- rnorm((times_i/T_max)^3, sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h


h = 4
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rnorm((times_i/T_max)^(.2), sd = sqrt(0.5))
  Z2 <- rnorm((times_i/T_max)^(.3), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h


h = 5
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rnorm(exp(times_i/T_max), sd = sqrt(0.5))
  Z2 <- rnorm(sin(times_i/T_max), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h



#Coefficient gamma (different for each rate and process 1,...,pY)
gamma_simul <- vector("list", length = pY)
for(h in 1:pY){
  gamma_simul[[h]] <- matrix(rnorm(q[h]*dY[h]*(dY[h] - 1), sd = 0.1),q[h],dY[h]*(dY[h] - 1),byrow = TRUE)
}

#Constant part of transition rates (needed in sim.msm, the time-varying covariates are included in the function used for the simulation)
lambda_Y_msm <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_msm_i <- list()
  for(h in 1:pY){
    lambda_Y_msm_i[[h]] <- lambda_Y_tilde_simul[[i]][[h]] * exp( X[[h]][i,] %*% beta_simul[[h]] )
  }
  lambda_Y_msm[[i]] <- lambda_Y_msm_i
}

#Generate the data (...a bit tricky with msm...
#...but we can do it independently and at different times for each process and patient...)

#Response process
Y <- vector("list", length = N)
lambda_Y_hat <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i,]
  
  lambda_Y_msm_i <- lambda_Y_msm[[i]]
  
  Y_i <- list()
  lambda_Y_hat_i <- list()
  for(h in 1:pY){
    
    lambda_Y_msm_i_h <- lambda_Y_msm_i[[h]]
    Q_Y <- matrix(0, dY[h], dY[h])
    for(i2 in 1:dY[h]){
      Q_Y[i2,c(1:dY[h])[-i2]] <- lambda_Y_msm_i_h[(i2-1)*(dY[h]-1) + c(1:(dY[h]-1))]
      Q_Y[i2,i2] <- -sum(Q_Y[i2,])
    }
    
    #Here we can include the time-varying covariates via sim.msm
    times_i_h <- times_i[[h]]
    n_times_i_h <- n_times_i[h]
    
    Z_msm_i_h <- Z_msm[[h]][[i]]
    gamma_simul_h <- gamma_simul[[h]]
    
    Y_i_h_msm <- sim.msm(qmatrix = Q_Y, maxtime = T_max, covs = Z_msm_i_h, beta = gamma_simul_h, obstimes = times_i_h, start = 2, mintime = 0)
    # Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
    
    out_aux <- outer(times_i_h[2:(n_times_i_h-1)], Y_i_h_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    Y_i[[h]] <- Y_i_h_msm$states[c(1,ind-1,length(Y_i_h_msm$times))] - 1 #(0,1)
    
    #For estimation of the initial transition rates we do not need any special q matrix
    #It only has to indicate id there are some constrained zeros in the matrix to be estimated (not for us I think)
    #Problems arise when there are no or only one transitions (zeros)
    data_Y_i_df <- data.frame(states = Y_i[[h]]+1, time = times_i_h)
    lambda_aux <- crudeinits.msm(states ~ time, data = data_Y_i_df, qmatrix = Q_Y)
    #If zeroes we impose arbitrary value
    lambda_aux[lambda_aux == 0] <- 0.5
    lambda_aux_vec <- lambda_aux[1,1]
    for(i2 in 1:dY[h]){
      lambda_aux_vec <- c(lambda_aux_vec, lambda_aux[i2,c(1:dY[h])[-i2]])
    }
    lambda_aux_vec <- lambda_aux_vec[-1]
    lambda_Y_hat_i[[h]] <- lambda_aux_vec
  }
  Y[[i]] <- Y_i
  lambda_Y_hat[[i]] <- lambda_Y_hat_i
}


#Put the estimated lambdas in a list
lambda_hat_all <- matrix(0,N,p_tot)
for(i in 1:N){
  lambda_hat_i <- 0
  for(h in 1:pY){
    lambda_hat_i <- c(lambda_hat_i, lambda_Y_hat[[i]][[h]])
  }
  lambda_hat_all[i,] <- lambda_hat_i[2:(p_tot+1)]
}


#Interval lengths
epsilon <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  epsilon_i <- vector("list", length = p0)
  for(ip in 1:p0){
    epsilon_i[[ip]] <- diff(times_i[[ip]])
  }
  epsilon[[i]] <- epsilon_i
}

data <- list(Y = Y, lambda_hat_all = lambda_hat_all, n_states = n_states, n_rates = n_rates, n_times = n_times, epsilon = epsilon, X = X, Z = Z)

## Now that we have the data, we need to initialize the parameters and priors
n_burn1 <- 10
n_burn2 <- 10
n_save <- 100
thin <- 1
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8

Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

#Prior on mu
update_mu <- TRUE
if(update_mu){# We provide hyperparameters for mu
  m_mu <- rep(0,p_tot)
  mu <- m_mu
  mu_list <- list(update_mu = update_mu, mu = mu, m_mu = m_mu)
}else{# We fix the values of mu
  mu <- mu_simul
  mu_list <- list(update_mu = update_mu, mu = mu)
}
Param_list <- modifyList(mu_list,Param_list)

#Prior on m_mu
update_m_mu <- FALSE
if(update_mu){# We provide hyperparameters for m_mu
  Sigma_m_mu <- diag(p_tot)
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu, Sigma_m_mu = Sigma_m_mu)
}else{# We fix the values of m_mu
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu)
}
Param_list <- modifyList(m_mu_list,Param_list)

#Prior on k0
update_k0 <- TRUE
if(update_k0){# We provide hyperparameters for k0
  mm <- 1
  ss <- 1
  a_k0 <- mm^2/ss
  b_k0 <- mm/ss
  k0 <- mm
  k0_list <- list(update_k0 = update_k0, k0 = k0, a_k0 = a_k0, b_k0 = b_k0)
}else{# We fix the values of k0
  k0 <- k0_simul
  k0_list <- list(update_k0 = update_k0, k0 = k0)
}
Param_list <- modifyList(k0_list,Param_list)

#Prior on graph
size_based_prior <- TRUE
if(size_based_prior){
  #The name of the parameters is the same as for d, because it's like marginalizeing wrt d
  a_d <- 1
  b_d <- 1
  G_list <- list(size_based_prior = size_based_prior, a_d = a_d, b_d = b_d)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_d <- TRUE
  if(update_d){
    #Different prior options for d
    d_beta_prior <- TRUE
    if(d_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_d <- mm*(mm*(1 - mm)/ss - 1)
      b_d <- (mm*(1 - mm)^2/ss - (1 - mm))
      d <- mm
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_d = a_d, b_d = b_d, d = d)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      d <- (1 + exp(-mm))^(-1)
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_lambda = a_lambda, d = d)
    }
  }else{# We fix the values of d
    d <- 0.25
    d_list <- list(update_d = update_d, d = d)
  }
  Param_list <- modifyList(d_list,Param_list)
}
Param_list <- modifyList(G_list,Param_list)


#Prior on beta
mu_beta <- vector("list", length = pY)
U_beta <- vector("list", length = pY)
V_beta <- vector("list", length = pY)
for(h in 1:pY){
  mu_beta[[h]] <- matrix(0,g[h],dY[h]*(dY[h] - 1))
  U_beta[[h]] <- diag(g[h])
  V_beta[[h]] <- diag(dY[h]*(dY[h] - 1))
}
beta_list <- list(mu_beta = mu_beta, U_beta = U_beta, V_beta = V_beta)
Param_list <- modifyList(beta_list, Param_list)

#Prior on gamma
mu_gamma <- vector("list", length = pY)
U_gamma <- vector("list", length = pY)
V_gamma <- vector("list", length = pY)
for(h in 1:pY){
  mu_gamma[[h]] <- matrix(0,q[h],dY[h]*(dY[h] - 1))
  U_gamma[[h]] <- diag(q[h])
  V_gamma[[h]] <- diag(dY[h]*(dY[h] - 1))
}
gamma_list <- list(mu_gamma = mu_gamma, U_gamma = U_gamma, V_gamma = V_gamma)
Param_list <- modifyList(gamma_list, Param_list)

#A-priori expected number of components
Mm <- 5

is_DP <- TRUE
is_parametric <- FALSE
if(is_parametric){#We can fit a model for a given clustering
  c_init <- c(1:N)-1
  alpha_list <- list(update_alpha = FALSE, alpha = 0)
}else{
  if(is_DP){#We can fit a DP model
    #Prior on alpha
    # #Based on prior expected number of components
    alpha <- matrix(seq(0.001,2,length = 50000),50000,1)
    sum_alpha <- apply(alpha,1,function(alpha) sum(alpha/(alpha + c(1:N) - 1)))
    alpha_Mm <- alpha[which.min(abs(Mm - sum_alpha))]
    sum_alpha[which.min(abs(Mm - sum_alpha))]
    
    update_alpha <- TRUE
    if(update_alpha){# We provide hyperparameters for alpha
      mm <- alpha_Mm
      ss <- 1
      a_alpha <- mm^2/ss
      b_alpha <- mm/ss
      alpha <- mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha, a_alpha = a_alpha, b_alpha = b_alpha)
    }else{# We fix the values of alpha
      alpha <- alpha_Mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha)
    }
  }else{
    alpha_list <- list(update_alpha = FALSE, alpha = 0)
  }
  c_init <- rep(1,N)-1
}
Param_list <- modifyList(alpha_list,Param_list)

model_list <- list(is_parametric = is_parametric, is_DP = is_DP , c_init = c_init)

# Main Gibbs sampler
sample.bdmcmc <- DisProgr( data = data, model_list = model_list, Alg_list = Alg_list, Param_list = Param_list)




#Produce plots for report

Graph_list <- sample.bdmcmc$Graph_List

Omega_out <- Graph_list$Omega_out
G_out <- Graph_list$G_out
G0_out <- Graph_list$G0_out

Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

layout(matrix(c(1,2),nrow = 1, ncol = 2))

col_val <- c(Omega_mean, Omega_simul)
col_map <- colorRampPalette(c("white", "red"))(99)

image(Omega_mean, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_mean, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))
image(Omega_simul, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_simul, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))

image(G_mean)
image(G_simul)

image(G0_mean)
image(G0_simul)



sum_weights_out <- 1/Graph_list$sum_weights_out
Omega_meanW <- matrix(0,p_tot,p_tot)
G_meanW <- matrix(0,p_tot,p_tot)
G0_meanW <- matrix(0,p0,p0)
for(it in 1:n_save){
  Omega_meanW = Omega_meanW + Omega_out[,,it] * sum_weights_out[it]
  G_meanW = G_meanW + G_out[,,it] * sum_weights_out[it]
  G0_meanW = G0_meanW + G0_out[,,it] * sum_weights_out[it]
}
Omega_meanW <- Omega_meanW / sum(sum_weights_out)
G_meanW <- G_meanW / sum(sum_weights_out)
G0_meanW <- G0_meanW / sum(sum_weights_out)

image(Omega_meanW)
image(Omega_simul)

image(G_meanW)
image(G_simul)

image(G0_meanW)
image(G0_simul)

dev.off()

#d
d_out <- Graph_list$d_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20)
abline(v = d_simul, lwd = 3, col = "red")
plot(d_out, type = "l")

#mu
mu_out <- sample.bdmcmc$mu_out
matplot(t(mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~mu), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#m_mu
m_mu_out <- sample.bdmcmc$m_mu_out
matplot(t(m_mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~m[mu]), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(m_mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), m_mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(m_mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#k0
k0_out <- sample.bdmcmc$k0_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), breaks = 20)
abline(v = k0_simul, lwd = 3, col = "red")
plot(k0_out, type = "l")



#Phi
phi_star_out <- sample.bdmcmc$phi_star_out
c_out <- sample.bdmcmc$c_out
phi_mean <- matrix(0, N, p_tot)
for(it in 1:n_save){
  phi_mean <- phi_mean + phi_star_out[[it]][c_out[it,]+1,]
}
phi_mean <- phi_mean/n_save
matplot(c(1:p_tot), t(phi_mean), col = "grey", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~phi), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
matplot(c(1:p_tot), t(phi_simul), col = "red", lwd = 1.5, lty = 1, type = "l", add = TRUE)
legend("topright", legend = c("MCMC", "truth"), pch = 19, col = c("grey", "red"), bty = "n", cex = 2)


#Number of clusters
K_N_out <- rep(0,n_save)
for(it in 1:n_save){
  K_N_out[it] <- dim(phi_star_out[[it]])[1]
}
plot(table(K_N_out)/n_save, col = "blue", lwd = 2)


#alpha
alpha_out <- sample.bdmcmc$alpha_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(alpha_out, col = "lightblue", main = bquote("Posterior of "~alpha), breaks = 20)
plot(alpha_out, type = "l")



#beta
beta_out <- sample.bdmcmc$beta_out
for(h in 1:pY){
  beta_out_h <- array(0, dim = c(g[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    beta_out_h[,,it] <- beta_out[[it]][[h]]
  }
  
  beta_df <- data.frame(dim_g = rep(c(1:(g[h]*dY[h]*(dY[h] - 1))), n_save), beta_vec = c(beta_out_h))
  beta_true_df <- data.frame( x = c(1:(g[h]*dY[h]*(dY[h] - 1))), y = c(beta_simul[[h]]) ) 
  
  ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
    geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
    labs(x = "g x dY", y = "", title = bquote("Posterior of "~beta)) +
    theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
}


#gamma
gamma_out <- sample.bdmcmc$gamma_out
for(h in 1:pY){
  gamma_out_h <- array(0, dim = c(q[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    gamma_out_h[,,it] <- gamma_out[[it]][[h]]
  }
  
  gamma_df <- data.frame(dim_g = rep(c(1:(q[h]*dY[h]*(dY[h] - 1))), n_save), gamma_vec = c(gamma_out_h))
  gamma_true_df <- data.frame( x = c(1:(q[h]*dY[h]*(dY[h] - 1))), y = c(gamma_simul[[h]]) ) 
  
  ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
    geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
    labs(x = "g x dY", y = "", title = bquote("Posterior of "~gamma)) +
    theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
}


#To see how the transition rates are estimated I want to show the evolution of trajectories

phi_star_out <- sample.bdmcmc$phi_star_out
c_out <- sample.bdmcmc$c_out
beta_out <- sample.bdmcmc$beta_out
gamma_out <- sample.bdmcmc$gamma_out

#Reconstruct transition rates for plots
#It's parametric so phi_star_out is a matrix
lambda_Y_out <- vector("list", length = n_save)
for(it in 1:n_save){
  phi_star_it_out <- phi_star_out[[it]]
  
  beta_out_it <- beta_out[[it]]
  gamma_out_it <- gamma_out[[it]]
  
  lambda_Yi_out <- vector("list", length = N)
  for(i in 1:N){
    n_times_i <- n_times[i,]
    lambda_Y_i <- vector("list", length = pY)
    for(h in 1:pY){
      ind_dY <- (n_rates_cum[h]+1):n_rates_cum[h+1]
      
      n_times_i_h <- n_times_i[h]
      lambda_Y_i_h <- matrix(0, n_times_i_h, dY[h]*(dY[h]-1))
      
      X_h <- X[[h]]
      Z_h <- Z[[h]]
      Z_h_i <- Z_h[[i]]
      beta_it_h_out <- beta_out_it[[h]]
      gamma_it_h_out <- gamma_out_it[[h]]
      
      lambda_Y_i[[h]] <- matrix(exp( phi_star_it_out[c_out[it,i]+1,ind_dY] + X_h[i,] %*% beta_it_h_out), nrow = n_times_i[h]-1, ncol = dY[h]*(dY[h]-1), byrow = TRUE) * exp( Z_h_i %*% gamma_it_h_out )
    }

    lambda_Yi_out[[i]] <- lambda_Y_i
  }
  lambda_Y_out[[it]] <- lambda_Yi_out
}

#Reconstruct simulated transition rates
lambda_Y_simul <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i,]
  lambda_Y_simul_i <- vector("list", length = pY)
  for(h in 1:pY){
    n_times_i_h <- n_times_i[h]
    X_h <- X[[h]]
    Z_h <- Z[[h]]
    Z_h_i <- Z_h[[i]]
    beta_simul_h <- beta_simul[[h]]
    gamma_simul_h <- gamma_simul[[h]]
    lambda_Y_simul_i[[h]] <- matrix(lambda_Y_tilde_simul[[c_simul[i]]][[h]] * exp( X_h[i,] %*% beta_simul_h ), n_times_i_h - 1, dY[h] * (dY[h] - 1), byrow = TRUE) * exp( Z_h_i %*% gamma_simul_h )
  }
  lambda_Y_simul[[i]] <- lambda_Y_simul_i
}


#Select patient
i_plot <- 28
Y_i <- Y[[i_plot]]
times_i_plot <- times[[i_plot]]
n_times_i_plot <- n_times[i_plot,]


#Simulated transition probabilities
prs_simul_log <- vector("list", length(p0))

lambda_Y_simul_i_plot <- lambda_Y_simul[[i_plot]]
epsilon_i_plot <- epsilon[[i_plot]]
for(h in 1:pY){
  n_times_i_h_plot <- n_times_i_plot[h]
  prs_simul_log_h <- rep(0,n_times_i_h_plot)
  Y_i_h <- Y_i[[h]]
  lambda_Y_simul_i_h_plot <- lambda_Y_simul_i_plot[[h]]
  epsilon_i_h_plot <- epsilon_i_plot[[h]]
  for(j in 2:n_times_i_h_plot){
    Yihj <- Y_i_h[j-1]
    prs_simul_log_h[j] <- log(lambda_Y_simul_i_h_plot[j-1,Yihj+1]) - log(sum(lambda_Y_simul_i_h_plot[j-1,])) + log( 1 - exp(- ((sum(lambda_Y_simul_i_h_plot[j-1,])) * epsilon_i_h_plot[j-1])) )
  }
  prs_simul_log[[h]] <- prs_simul_log_h
}

#Estimated transition probabilities
prs_out_log <- vector("list", length = p0)

for(h in 1:pY){
  n_times_i_h_plot <- n_times_i_plot[h]
  prs_out_log_h <- matrix(0,n_times_i_h_plot,n_save)
  Y_i_h <- Y_i[[h]]
  epsilon_i_h_plot <- epsilon_i_plot[[h]]
  for(j in 2:n_times_i_h_plot){
    Yihj <- Y_i_h[j-1]
    for(it in 1:n_save){
      lambda_Y_it_i_plot_h_out <- lambda_Y_out[[it]][[i_plot]][[h]]
      prs_out_log_h[j,it] <- log(lambda_Y_it_i_plot_h_out[j-1,Yihj+1]) - log(sum(lambda_Y_it_i_plot_h_out[j-1,])) + log( 1 - exp(- (sum(lambda_Y_it_i_plot_h_out[j-1,]) * epsilon_i_h_plot[j-1])) )
    }
  }
  prs_out_log[[h]] <- prs_out_log_h
}


#Colours of processes
col_proc <- rgb(255/255,51/255,51/255)

scale_red1 <- seq(0, 120, length = pY)
scale_red2 <- seq(0, 255, length = pY)
for(h in 1:pY){
  col_proc <- c(col_proc, rgb(255/255,scale_red1[h]/255,scale_red2[h]/255))
}

#Labels of processes
proc_lab <- paste0("Y", c(1:pY))

#Credibility intervals
prs_CI_low <- vector("list", length = p0)
prs_CI_up <- vector("list", length = p0)
for(ip in 1:p0){
  prs_CI_low[[ip]] <- apply(prs_out_log[[ip]],1, quantile, probs = 0.025)
  prs_CI_up[[ip]] <- apply(prs_out_log[[ip]],1, quantile, probs = 0.975)
}
prs_CI_df <- data.frame(prs_CI_low = unlist(prs_CI_low), prs_CI_up = unlist(prs_CI_up), proc = unlist(mapply(rep,proc_lab,n_times_i_plot), use.names = FALSE), times_vec = round(unlist(times_i_plot),3) )
prs_simul_df <- data.frame(x_times = unlist(times_i_plot), prs_simul_vec = unlist(prs_simul_log), proc = unlist(mapply(rep,proc_lab,n_times_i_plot), use.names = FALSE), times_vec = round(unlist(times_i_plot),3) )
prs_df <- data.frame(x_times = unlist(mapply(rep,times_i_plot,n_save)), prs_vec = unlist(prs_out_log), proc = unlist(mapply(rep,proc_lab,n_times_i_plot * n_save), use.names = FALSE), times_vec = round(unlist(mapply(rep,times_i_plot,n_save)),3) )


ggplot(prs_df, aes(x = x_times, y = prs_vec, col = proc)) + 
  geom_line(data = prs_simul_df, aes(x = x_times, y = prs_simul_vec, col = proc), size = 1.25, linetype = 2) + 
  stat_summary(fun.y = 'mean', geom = 'line', size = 1.25) +
  stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2) + 
  labs(x = "time", y = bquote(log(p[t]^{rs})), title = bquote("Posterior of "~log(p[t]^{rs})~" for subject "~.(i_plot))) +
  theme(text = element_text(size=35), axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 15)) +
  facet_grid(rows = vars(proc), scales="free") +
  # the labels must match what you specified above
  scale_fill_manual(name = proc_lab, values = col_proc) +
  scale_color_manual(name = proc_lab, values = col_proc) +
  theme_bw()





#####Binder loss function equal costs
K_Binder <- 1/2
#For selected combinations of parameters
c_out <- sample.bdmcmc$c_out

pij <- matrix(0,N,N)
cij <- vector("list", length = n_save)
for(g in 1:n_save){
  cij[[g]] <- outer(c_out[g,], c_out[g,], "==")
  pij <- pij + cij[[g]]
}
pij <- pij/n_save


Binder_f <- rep(0,n_save)
for(g in 1:n_save){
  aux <- (pij - K_Binder) * as.matrix(cij[[g]])
  aux <-  aux[upper.tri(aux)]
  Binder_f[g] <- sum(aux)
}
Binder_ind <- which.max(Binder_f)
Binder_out <- c_out[Binder_ind,] + 1

graphics.off()
library("gplots")
Binder_sort <- sort(Binder_out, index.return = TRUE)
Binder_sort <- Binder_sort$ix

heatmap.2(pij[Binder_sort,Binder_sort], ColSideColor = rainbow(K_N_simul)[c_simul[Binder_sort]], col=rev(heat.colors(100)), trace = "none", key  = TRUE, Rowv = FALSE, Colv = FALSE, dendrogram = "none", main = "Posterior co-clustering probabilities", cex.main = 3, cex.lab = 2, cex.axis = 1.5)














#################################################################################################################
# dH = 0
#################################################################################################################

## Here we do a fuller simulation from an AntMAN Mixture model where we also have the Markov processes
## We will try to handle directly different time points for different processes and patients
## We will try also SM algorithm


rm(list = ls())
set.seed(123)
library("BayesDiseaseProgr")
library("ggplot2")
library("msm")
library("fields")


######
#Simulate sensible data
######

#Number of patients
N = 50;
#Number of response processes
pY = 5;
#Number of covariate processes
pH = 0;
#Tot number of porcesses
p0 = pY+pH;

#Set graph structure(s)
# d_simul <- 0.25
G0_simul <- matrix(0, p0, p0)
# G0_simul[lower.tri(G0_simul)] <- c(runif(p0*(p0-1)/2) <= d_simul)
G0_simul[1,2] <- 1
G0_simul[1,4] <- 1
G0_simul[2,3] <- 1
G0_simul[3,4] <- 1
G0_simul[4,5] <- 1
G0_simul[2,5] <- 1

G0_simul <- G0_simul + t(G0_simul)

#Assume 2-state responses
dY <- c(2, 2, 3, 4, 5)
dH <- NULL
n_states <- list(dY,dH)
n_rates <- 0
for(h in 1:pY){
  n_rates <- c(n_rates, dY[h] * (dY[h]-1))
}
n_rates_cum <- cumsum(n_rates)
p_tot <- sum(n_rates)
G_simul <- matrix(0, nrow = p_tot, ncol = p_tot)
diag(G0_simul) <- rep(1,p0)
for(i in 1:p0){
  for(j in i:p0){
    if(G0_simul[i,j]){
      G_simul[c((n_rates_cum[i]+1) : n_rates_cum[i+1]), c((n_rates_cum[j]+1) : n_rates_cum[j+1])] <- 1
      G_simul[c((n_rates_cum[j]+1) : n_rates_cum[j+1]), c((n_rates_cum[i]+1) : n_rates_cum[i+1])] <- 1 #For symmetry
    }
  }
}
diag(G0_simul) <- rep(0,p0)
diag(G_simul) <- rep(0,p_tot)

#Prior specification
nu_simul <- 5
Psi_simul <- diag(p_tot)

#Simulate data
Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

K_N_simul <- 3
c_simul <- sample.int(K_N_simul, size = N, replace = TRUE, prob = c(1/3,1/6,1/2))
phi_simul <- rmvnorm(n = K_N_simul, mean = c(mu_simul), sigma = solve(Omega_simul))
lambda_simul <- exp(phi_simul)

#Number of visits per patient per process (for msm package)
n_times <- matrix(9 + rpois(N*p0, 1), N, p0)
#Times of visit (0 and Tmax included!)
T_max <- 30
times <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i,]
  times_i <- list()
  for(ip in 1:p0){
    times_i[[ip]] <- sort(c(0,runif(n_times_i[ip]-2, min = 0, max = T_max), T_max))
  }
  times[[i]] <- times_i
}

#Transition rates for Y are covariate- and time-dependent
lambda_Y_tilde_simul <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_tilde_simul_i <- list()
  for(h in 1:pY){
    lambda_Y_tilde_simul_i[[h]] <- lambda_simul[c_simul[i],(n_rates_cum[h]+1):n_rates_cum[h+1]]
  }
  lambda_Y_tilde_simul[[i]] <- lambda_Y_tilde_simul_i
}

#Constant covariates, but they can differ in each process i,...,pY
g <- c(2, 4, 3, 2, 1) #Simulate age and sex without intercept (it's the lambda_tilde!)
X <- list()

aux <- matrix(1,N,g[1])
aux[,1] <- rnorm(N, mean = 0, sd = 1) #standardized
aux[,2] <- (runif(N) <= 0.5)
X[[1]] <- aux

X[[2]] <- matrix(rgamma(N*g[2], 1, 1), N, g[2])

X[[3]] <- cbind(runif(N) <= 0.25, rnorm(N, mean = -1, sd = 1), runif(N) <= 0.1)

aux <- matrix(1,N,g[1])
aux[,1] <- rnorm(N, mean = 2, sd = 1) #standardized
aux[,2] <- (runif(N) > 0.3)
X[[4]] <- aux

X[[5]] <- matrix(runif(N) <= 0.15)


#Coefficient beta (different for each rate and process 1,...,pY)
beta_simul <- vector("list", length = pY)
for(h in 1:pY){
  beta_simul[[h]] <- matrix(rnorm(g[h]*dY[h]*(dY[h] - 1), sd = 0.1),g[h],dY[h]*(dY[h] - 1),byrow = TRUE)
}

#Time-varying covariates (not random for msm package, no need for averaging)
#Also depends on the processes
#Here we put one equal to zero
q <- c(2,3,2,2,2)
Z_msm <- vector("list", length = pY)
Z <- vector("list", length = pY)

#Careful with dimension q!
Z_msm_h <- vector("list", length = N)
Z_h <- vector("list", length = N)
h = 1
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm(times_i/T_max, sd = sqrt(0.5))
  Z2 <- rnorm(cos(times_i/T_max*2*pi), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h

h = 2
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rep(0,n_times_i)
  Z2 <- rep(0,n_times_i)
  Z3 <- rep(0,n_times_i)
  Z_msm_h[[i]] <- cbind(Z1,Z2,Z3)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z3 <- (Z3[1:(n_times_i-1)] + Z3[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2,Z3)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h

h = 3
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rnorm((times_i/T_max)^2, sd = sqrt(0.5))
  Z2 <- rnorm((times_i/T_max)^3, sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h


h = 4
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rnorm((times_i/T_max)^(.2), sd = sqrt(0.5))
  Z2 <- rnorm((times_i/T_max)^(.3), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h


h = 5
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #Zero-covariates to try
  Z1 <- rnorm(exp(times_i/T_max), sd = sqrt(0.5))
  Z2 <- rnorm(sin(times_i/T_max), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h



#Coefficient gamma (different for each rate and process 1,...,pY)
gamma_simul <- vector("list", length = pY)
for(h in 1:pY){
  gamma_simul[[h]] <- matrix(rnorm(q[h]*dY[h]*(dY[h] - 1), sd = 0.1),q[h],dY[h]*(dY[h] - 1),byrow = TRUE)
}

#Constant part of transition rates (needed in sim.msm, the time-varying covariates are included in the function used for the simulation)
lambda_Y_msm <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_msm_i <- list()
  for(h in 1:pY){
    lambda_Y_msm_i[[h]] <- lambda_Y_tilde_simul[[i]][[h]] * exp( X[[h]][i,] %*% beta_simul[[h]] )
  }
  lambda_Y_msm[[i]] <- lambda_Y_msm_i
}

#Generate the data (...a bit tricky with msm...
#...but we can do it independently and at different times for each process and patient...)

#Response process
Y <- vector("list", length = N)
lambda_Y_hat <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i,]
  
  lambda_Y_msm_i <- lambda_Y_msm[[i]]
  
  Y_i <- list()
  lambda_Y_hat_i <- list()
  for(h in 1:pY){
    
    lambda_Y_msm_i_h <- lambda_Y_msm_i[[h]]
    Q_Y <- matrix(0, dY[h], dY[h])
    for(i2 in 1:dY[h]){
      Q_Y[i2,c(1:dY[h])[-i2]] <- lambda_Y_msm_i_h[(i2-1)*(dY[h]-1) + c(1:(dY[h]-1))]
      Q_Y[i2,i2] <- -sum(Q_Y[i2,])
    }
    
    #Here we can include the time-varying covariates via sim.msm
    times_i_h <- times_i[[h]]
    n_times_i_h <- n_times_i[h]
    
    Z_msm_i_h <- Z_msm[[h]][[i]]
    gamma_simul_h <- gamma_simul[[h]]
    
    Y_i_h_msm <- sim.msm(qmatrix = Q_Y, maxtime = T_max, covs = Z_msm_i_h, beta = gamma_simul_h, obstimes = times_i_h, start = 2, mintime = 0)
    # Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
    
    out_aux <- outer(times_i_h[2:(n_times_i_h-1)], Y_i_h_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    Y_i[[h]] <- Y_i_h_msm$states[c(1,ind-1,length(Y_i_h_msm$times))] - 1 #(0,1)
    
    #For estimation of the initial transition rates we do not need any special q matrix
    #It only has to indicate id there are some constrained zeros in the matrix to be estimated (not for us I think)
    #Problems arise when there are no or only one transitions (zeros)
    data_Y_i_df <- data.frame(states = Y_i[[h]]+1, time = times_i_h)
    lambda_aux <- crudeinits.msm(states ~ time, data = data_Y_i_df, qmatrix = Q_Y)
    #If zeroes we impose arbitrary value
    lambda_aux[lambda_aux == 0] <- 0.5
    lambda_aux_vec <- lambda_aux[1,1]
    for(i2 in 1:dY[h]){
      lambda_aux_vec <- c(lambda_aux_vec, lambda_aux[i2,c(1:dY[h])[-i2]])
    }
    lambda_aux_vec <- lambda_aux_vec[-1]
    lambda_Y_hat_i[[h]] <- lambda_aux_vec
  }
  Y[[i]] <- Y_i
  lambda_Y_hat[[i]] <- lambda_Y_hat_i
}


#Put the estimated lambdas in a list
lambda_hat_all <- matrix(0,N,p_tot)
for(i in 1:N){
  lambda_hat_i <- 0
  for(h in 1:pY){
    lambda_hat_i <- c(lambda_hat_i, lambda_Y_hat[[i]][[h]])
  }
  lambda_hat_all[i,] <- lambda_hat_i[2:(p_tot+1)]
}


#Interval lengths
epsilon <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  epsilon_i <- vector("list", length = p0)
  for(ip in 1:p0){
    epsilon_i[[ip]] <- diff(times_i[[ip]])
  }
  epsilon[[i]] <- epsilon_i
}

data <- list(Y = Y, lambda_hat_all = lambda_hat_all, n_states = n_states, n_rates = n_rates, n_times = n_times, epsilon = epsilon, X = X, Z = Z)

## Now that we have the data, we need to initialize the parameters and priors
n_burn1 <- 10
n_burn2 <- 10
n_save <- 100
thin <- 1
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8

SM_alg <- TRUE

Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold, SM_alg = SM_alg, t_SM_split = 3, t_SM_merge = 3)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

#Prior on Lambda
update_Lambda <- TRUE
if(update_Lambda){# We provide hyperparameters for Lambda
  #These are based on mean and var for M (number of components)
  mm <- 1
  ss <- 1 #Variance quite small?
  a2 <- mm^2/ss
  b2 <- mm/ss
  Lambda <- mm
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda, a2 = a2, b2 = b2)
}else{# We fix the values of Lambda
  Lambda <- 1
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda)
}
Param_list <- modifyList(Lambda_list,Param_list)

update_gamma_S <- TRUE
if(update_gamma_S){# We provide hyperparameters for gamma_S
  mm <- 1
  ss <- 1
  a1 <- mm^2/ss
  b1 <- mm/ss
  gamma_S <- mm
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S, a1 = a1, b1 = b1)
}else{# We fix the values of gamma_S
  gamma_S <- 1
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S)
}
Param_list <- modifyList(gamma_S_list,Param_list)

#Prior on mu
update_mu <- TRUE
if(update_mu){# We provide hyperparameters for mu
  m_mu <- rep(0,p_tot)
  mu <- m_mu
  mu_list <- list(update_mu = update_mu, mu = mu, m_mu = m_mu)
}else{# We fix the values of mu
  mu <- mu_simul
  mu_list <- list(update_mu = update_mu, mu = mu)
}
Param_list <- modifyList(mu_list,Param_list)

#Prior on m_mu
update_m_mu <- FALSE
if(update_mu){# We provide hyperparameters for m_mu
  Sigma_m_mu <- diag(p_tot)
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu, Sigma_m_mu = Sigma_m_mu)
}else{# We fix the values of m_mu
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu)
}
Param_list <- modifyList(m_mu_list,Param_list)

#Prior on k0
update_k0 <- TRUE
if(update_k0){# We provide hyperparameters for k0
  mm <- 1
  ss <- 1
  a_k0 <- mm^2/ss
  b_k0 <- mm/ss
  k0 <- mm
  k0_list <- list(update_k0 = update_k0, k0 = k0, a_k0 = a_k0, b_k0 = b_k0)
}else{# We fix the values of k0
  k0 <- k0_simul
  k0_list <- list(update_k0 = update_k0, k0 = k0)
}
Param_list <- modifyList(k0_list,Param_list)

#Prior on graph
size_based_prior <- TRUE
if(size_based_prior){
  #The name of the parameters is the same as for d, because it's like marginalizeing wrt d
  a_d <- 1
  b_d <- 1
  G_list <- list(size_based_prior = size_based_prior, a_d = a_d, b_d = b_d)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_d <- TRUE
  if(update_d){
    #Different prior options for d
    d_beta_prior <- TRUE
    if(d_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_d <- mm*(mm*(1 - mm)/ss - 1)
      b_d <- (mm*(1 - mm)^2/ss - (1 - mm))
      d <- mm
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_d = a_d, b_d = b_d, d = d)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      d <- (1 + exp(-mm))^(-1)
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_lambda = a_lambda, d = d)
    }
  }else{# We fix the values of d
    d <- 0.25
    d_list <- list(update_d = update_d, d = d)
  }
  Param_list <- modifyList(d_list,Param_list)
}
Param_list <- modifyList(G_list,Param_list)


#Prior on beta
mu_beta <- vector("list", length = pY)
U_beta <- vector("list", length = pY)
V_beta <- vector("list", length = pY)
for(h in 1:pY){
  mu_beta[[h]] <- matrix(0,g[h],dY[h]*(dY[h] - 1))
  U_beta[[h]] <- diag(g[h])
  V_beta[[h]] <- diag(dY[h]*(dY[h] - 1))
}
beta_list <- list(mu_beta = mu_beta, U_beta = U_beta, V_beta = V_beta)
Param_list <- modifyList(beta_list, Param_list)

#Prior on gamma
mu_gamma <- vector("list", length = pY)
U_gamma <- vector("list", length = pY)
V_gamma <- vector("list", length = pY)
for(h in 1:pY){
  mu_gamma[[h]] <- matrix(0,q[h],dY[h]*(dY[h] - 1))
  U_gamma[[h]] <- diag(q[h])
  V_gamma[[h]] <- diag(dY[h]*(dY[h] - 1))
}
gamma_list <- list(mu_gamma = mu_gamma, U_gamma = U_gamma, V_gamma = V_gamma)
Param_list <- modifyList(gamma_list, Param_list)


is_DP <- FALSE
is_parametric <- FALSE
if(is_parametric){#We can fit a model for a given clustering
  c_init <- c(1:N)-1
  alpha_list <- list(update_alpha = FALSE, alpha = 0)
}else{
  if(is_DP){#We can fit a DP model
    #Prior on alpha
    # #Based on prior expected number of components
    alpha <- matrix(seq(0.001,2,length = 50000),50000,1)
    sum_alpha <- apply(alpha,1,function(alpha) sum(alpha/(alpha + c(1:N) - 1)))
    alpha_Mm <- alpha[which.min(abs(Mm - sum_alpha))]
    sum_alpha[which.min(abs(Mm - sum_alpha))]
    
    update_alpha <- TRUE
    if(update_alpha){# We provide hyperparameters for alpha
      mm <- alpha_Mm
      ss <- 1
      a_alpha <- mm^2/ss
      b_alpha <- mm/ss
      alpha <- mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha, a_alpha = a_alpha, b_alpha = b_alpha)
    }else{# We fix the values of alpha
      alpha <- alpha_Mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha)
    }
  }else{
    alpha_list <- list(update_alpha = FALSE, alpha = 0)
  }
  c_init <- rep(1,N)-1
}
Param_list <- modifyList(alpha_list,Param_list)

model_list <- list(is_parametric = is_parametric, is_DP = is_DP , c_init = c_init)

# Main Gibbs sampler
sample.bdmcmc <- DisProgr( data = data, model_list = model_list, Alg_list = Alg_list, Param_list = Param_list)




#Produce plots for report

Graph_list <- sample.bdmcmc$Graph_List

Omega_out <- Graph_list$Omega_out
G_out <- Graph_list$G_out
G0_out <- Graph_list$G0_out

Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

layout(matrix(c(1,2),nrow = 1, ncol = 2))

col_val <- c(Omega_mean, Omega_simul)
col_map <- colorRampPalette(c("white", "red"))(99)

image(Omega_mean, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_mean, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))
image(Omega_simul, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_simul, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))

image(G_mean)
image(G_simul)

image(G0_mean)
image(G0_simul)



sum_weights_out <- 1/Graph_list$sum_weights_out
Omega_meanW <- matrix(0,p_tot,p_tot)
G_meanW <- matrix(0,p_tot,p_tot)
G0_meanW <- matrix(0,p0,p0)
for(it in 1:n_save){
  Omega_meanW = Omega_meanW + Omega_out[,,it] * sum_weights_out[it]
  G_meanW = G_meanW + G_out[,,it] * sum_weights_out[it]
  G0_meanW = G0_meanW + G0_out[,,it] * sum_weights_out[it]
}
Omega_meanW <- Omega_meanW / sum(sum_weights_out)
G_meanW <- G_meanW / sum(sum_weights_out)
G0_meanW <- G0_meanW / sum(sum_weights_out)

image(Omega_meanW)
image(Omega_simul)

image(G_meanW)
image(G_simul)

image(G0_meanW)
image(G0_simul)

dev.off()

#d
d_out <- Graph_list$d_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20)
abline(v = d_simul, lwd = 3, col = "red")
plot(d_out, type = "l")

#mu
mu_out <- sample.bdmcmc$mu_out
matplot(t(mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~mu), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#m_mu
m_mu_out <- sample.bdmcmc$m_mu_out
matplot(t(m_mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~m[mu]), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(m_mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), m_mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(m_mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#k0
k0_out <- sample.bdmcmc$k0_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), breaks = 20)
abline(v = k0_simul, lwd = 3, col = "red")
plot(k0_out, type = "l")



#Phi
phi_star_out <- sample.bdmcmc$phi_star_out
c_out <- sample.bdmcmc$c_out
phi_mean <- matrix(0, N, p_tot)
for(it in 1:n_save){
  phi_mean <- phi_mean + phi_star_out[[it]][c_out[it,]+1,]
}
phi_mean <- phi_mean/n_save
matplot(c(1:p_tot), t(phi_mean), col = "grey", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~phi), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
matplot(c(1:p_tot), t(phi_simul), col = "red", lwd = 1.5, lty = 1, type = "l", add = TRUE)
legend("topright", legend = c("MCMC", "truth"), pch = 19, col = c("grey", "red"), bty = "n", cex = 2)


#Number of clusters
K_N_out <- rep(0,n_save)
for(it in 1:n_save){
  K_N_out[it] <- dim(phi_star_out[[it]])[1]
}
plot(table(K_N_out)/n_save, col = "blue", lwd = 2)


#beta
beta_out <- sample.bdmcmc$beta_out
for(h in 1:pY){
  beta_out_h <- array(0, dim = c(g[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    beta_out_h[,,it] <- beta_out[[it]][[h]]
  }
  
  beta_df <- data.frame(dim_g = rep(c(1:(g[h]*dY[h]*(dY[h] - 1))), n_save), beta_vec = c(beta_out_h))
  beta_true_df <- data.frame( x = c(1:(g[h]*dY[h]*(dY[h] - 1))), y = c(beta_simul[[h]]) ) 
  
  ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
    geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
    labs(x = "g x dY", y = "", title = bquote("Posterior of "~beta)) +
    theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
}


#gamma
gamma_out <- sample.bdmcmc$gamma_out
for(h in 1:pY){
  gamma_out_h <- array(0, dim = c(q[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    gamma_out_h[,,it] <- gamma_out[[it]][[h]]
  }
  
  gamma_df <- data.frame(dim_g = rep(c(1:(q[h]*dY[h]*(dY[h] - 1))), n_save), gamma_vec = c(gamma_out_h))
  gamma_true_df <- data.frame( x = c(1:(q[h]*dY[h]*(dY[h] - 1))), y = c(gamma_simul[[h]]) ) 
  
  ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
    geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
    labs(x = "g x dY", y = "", title = bquote("Posterior of "~gamma)) +
    theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
}



#Lambda
Lambda_out <- sample.bdmcmc$Lambda_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(Lambda_out, col = "lightblue", main = bquote("Posterior of "~Lambda), breaks = 20)
plot(Lambda_out, type = "l")


#gamma_S
gamma_S_out <- sample.bdmcmc$gamma_S_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(gamma_S_out, col = "lightblue", main = bquote("Posterior of "~gamma[S]), breaks = 20)
plot(gamma_S_out, type = "l")



#To see how the transition rates are estimated I want to show the evolution of trajectories

phi_star_out <- sample.bdmcmc$phi_star_out
c_out <- sample.bdmcmc$c_out
beta_out <- sample.bdmcmc$beta_out
gamma_out <- sample.bdmcmc$gamma_out

#Reconstruct transition rates for plots
#It's parametric so phi_star_out is a matrix
lambda_Y_out <- vector("list", length = n_save)
for(it in 1:n_save){
  phi_star_it_out <- phi_star_out[[it]]
  
  beta_out_it <- beta_out[[it]]
  gamma_out_it <- gamma_out[[it]]
  
  lambda_Yi_out <- vector("list", length = N)
  for(i in 1:N){
    n_times_i <- n_times[i,]
    lambda_Y_i <- vector("list", length = pY)
    for(h in 1:pY){
      ind_dY <- (n_rates_cum[h]+1):n_rates_cum[h+1]
      
      n_times_i_h <- n_times_i[h]
      lambda_Y_i_h <- matrix(0, n_times_i_h, dY[h]*(dY[h]-1))
      
      X_h <- X[[h]]
      Z_h <- Z[[h]]
      Z_h_i <- Z_h[[i]]
      beta_it_h_out <- beta_out_it[[h]]
      gamma_it_h_out <- gamma_out_it[[h]]
      
      lambda_Y_i[[h]] <- matrix(exp( phi_star_it_out[c_out[it,i]+1,ind_dY] + X_h[i,] %*% beta_it_h_out), nrow = n_times_i[h]-1, ncol = dY[h]*(dY[h]-1), byrow = TRUE) * exp( Z_h_i %*% gamma_it_h_out )
    }
    lambda_Yi_out[[i]] <- lambda_Y_i
  }
  lambda_Y_out[[it]] <- lambda_Yi_out
}

#Reconstruct simulated transition rates
lambda_Y_simul <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i,]
  lambda_Y_simul_i <- vector("list", length = pY)
  for(h in 1:pY){
    n_times_i_h <- n_times_i[h]
    X_h <- X[[h]]
    Z_h <- Z[[h]]
    Z_h_i <- Z_h[[i]]
    beta_simul_h <- beta_simul[[h]]
    gamma_simul_h <- gamma_simul[[h]]
    lambda_Y_simul_i[[h]] <- matrix(lambda_Y_tilde_simul[[c_simul[i]]][[h]] * exp( X_h[i,] %*% beta_simul_h ), n_times_i_h - 1, dY[h] * (dY[h] - 1), byrow = TRUE) * exp( Z_h_i %*% gamma_simul_h )
  }
  lambda_Y_simul[[i]] <- lambda_Y_simul_i
}


#Select patient
i_plot <- 28
Y_i <- Y[[i_plot]]
times_i_plot <- times[[i_plot]]
n_times_i_plot <- n_times[i_plot,]


#Simulated transition probabilities
prs_simul_log <- vector("list", length(p0))

lambda_Y_simul_i_plot <- lambda_Y_simul[[i_plot]]
epsilon_i_plot <- epsilon[[i_plot]]
for(h in 1:pY){
  n_times_i_h_plot <- n_times_i_plot[h]
  prs_simul_log_h <- rep(0,n_times_i_h_plot)
  Y_i_h <- Y_i[[h]]
  lambda_Y_simul_i_h_plot <- lambda_Y_simul_i_plot[[h]]
  epsilon_i_h_plot <- epsilon_i_plot[[h]]
  for(j in 2:n_times_i_h_plot){
    Yihj <- Y_i_h[j-1]
    prs_simul_log_h[j] <- log(lambda_Y_simul_i_h_plot[j-1,Yihj+1]) - log(sum(lambda_Y_simul_i_h_plot[j-1,])) + log( 1 - exp(- ((sum(lambda_Y_simul_i_h_plot[j-1,])) * epsilon_i_h_plot[j-1])) )
  }
  prs_simul_log[[h]] <- prs_simul_log_h
}


#Estimated transition probabilities
prs_out_log <- vector("list", length = p0)

for(h in 1:pY){
  n_times_i_h_plot <- n_times_i_plot[h]
  prs_out_log_h <- matrix(0,n_times_i_h_plot,n_save)
  Y_i_h <- Y_i[[h]]
  epsilon_i_h_plot <- epsilon_i_plot[[h]]
  for(j in 2:n_times_i_h_plot){
    Yihj <- Y_i_h[j-1]
    for(it in 1:n_save){
      lambda_Y_it_i_plot_h_out <- lambda_Y_out[[it]][[i_plot]][[h]]
      prs_out_log_h[j,it] <- log(lambda_Y_it_i_plot_h_out[j-1,Yihj+1]) - log(sum(lambda_Y_it_i_plot_h_out[j-1,])) + log( 1 - exp(- (sum(lambda_Y_it_i_plot_h_out[j-1,]) * epsilon_i_h_plot[j-1])) )
    }
  }
  prs_out_log[[h]] <- prs_out_log_h
}



#Colours of processes
col_proc <- rgb(255/255,51/255,51/255)

scale_red1 <- seq(0, 200, length = pY)
scale_red2 <- seq(0, 255, length = pY)
for(h in 1:pY){
  col_proc <- c(col_proc, rgb(255/255,scale_red1[h]/255,scale_red2[h]/255))
}

#Labels of processes
proc_lab <- paste0("Y", c(1:pY))

#Credibility intervals
prs_CI_low <- vector("list", length = p0)
prs_CI_up <- vector("list", length = p0)
for(ip in 1:p0){
  prs_CI_low[[ip]] <- apply(prs_out_log[[ip]],1, quantile, probs = 0.025)
  prs_CI_up[[ip]] <- apply(prs_out_log[[ip]],1, quantile, probs = 0.975)
}
prs_CI_df <- data.frame(prs_CI_low = unlist(prs_CI_low), prs_CI_up = unlist(prs_CI_up), proc = unlist(mapply(rep,proc_lab,n_times_i_plot), use.names = FALSE), times_vec = round(unlist(times_i_plot),3) )
prs_simul_df <- data.frame(x_times = unlist(times_i_plot), prs_simul_vec = unlist(prs_simul_log), proc = unlist(mapply(rep,proc_lab,n_times_i_plot), use.names = FALSE), times_vec = round(unlist(times_i_plot),3) )
prs_df <- data.frame(x_times = unlist(mapply(rep,times_i_plot,n_save)), prs_vec = unlist(prs_out_log), proc = unlist(mapply(rep,proc_lab,n_times_i_plot * n_save), use.names = FALSE), times_vec = round(unlist(mapply(rep,times_i_plot,n_save)),3) )


ggplot(prs_df, aes(x = x_times, y = prs_vec, col = proc)) + 
  geom_line(data = prs_simul_df, aes(x = x_times, y = prs_simul_vec, col = proc), size = 1.25, linetype = 2) + 
  stat_summary(fun.y = 'mean', geom = 'line', size = 1.25) +
  stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2) + 
  labs(x = "time", y = bquote(log(p[t]^{rs})), title = bquote("Posterior of "~log(p[t]^{rs})~" for subject "~.(i_plot))) +
  theme(text = element_text(size=35), axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 15)) +
  facet_grid(rows = vars(proc), scales="free") +
  # the labels must match what you specified above
  scale_fill_manual(name = proc_lab, values = col_proc) +
  scale_color_manual(name = proc_lab, values = col_proc) +
  theme_bw()





#####Binder loss function equal costs
K_Binder <- 1/2
#For selected combinations of parameters
c_out <- sample.bdmcmc$c_out

pij <- matrix(0,N,N)
cij <- vector("list", length = n_save)
for(g in 1:n_save){
  cij[[g]] <- outer(c_out[g,], c_out[g,], "==")
  pij <- pij + cij[[g]]
}
pij <- pij/n_save


Binder_f <- rep(0,n_save)
for(g in 1:n_save){
  aux <- (pij - K_Binder) * as.matrix(cij[[g]])
  aux <-  aux[upper.tri(aux)]
  Binder_f[g] <- sum(aux)
}
Binder_ind <- which.max(Binder_f)
Binder_out <- c_out[Binder_ind,] + 1

graphics.off()
library("gplots")
Binder_sort <- sort(Binder_out, index.return = TRUE)
Binder_sort <- Binder_sort$ix

heatmap.2(pij[Binder_sort,Binder_sort], ColSideColor = rainbow(K_N_simul)[c_simul[Binder_sort]], col=rev(heat.colors(100)), trace = "none", key  = TRUE, Rowv = FALSE, Colv = FALSE, dendrogram = "none", main = "Posterior co-clustering probabilities", cex.main = 3, cex.lab = 2, cex.axis = 1.5)




































#################################################################################################################
# EXAMPLE FROM MARIA (IT WORKS) #
#################################################################################################################

## Here we do a fuller simulation from a semi-parametric model where we also have the Markov processes
## We will try to handle directly different time points for different processes and patients
## CLUSTERING from mixture with random number of components

rm(list = ls())
set.seed(123)
library("BayesDiseaseProgr")
library("ggplot2")
library("fields")
library("msm")


######
#Simulate sensible data
######

#Number of patients
N = 150;
#Number of binary covariates observed
p = 2;
#Tot number of porcesses
p0 = p+1;

#Set graph structure(s)
# d_simul <- 0.25
G0_simul <- matrix(0, p0, p0)
G0_simul[1,2] <- 1
# G0_simul[lower.tri(G0_simul)] <- c(runif(p0*(p0-1)/2) <= d_simul)
G0_simul <- G0_simul + t(G0_simul)

#Assume 2-state responses
dY <- 2
n_rates <- c(dY, rep(2,p0-1))
n_rates_cum <- c(0, cumsum(n_rates))
p_tot <- sum(n_rates)
G_simul <- matrix(0, nrow = p_tot, ncol = p_tot)
diag(G0_simul) <- rep(1,p0)
for(i in 1:p0){
  for(j in i:p0){
    if(G0_simul[i,j]){
      G_simul[c((n_rates_cum[i]+1) : n_rates_cum[i+1]), c((n_rates_cum[j]+1) : n_rates_cum[j+1])] <- 1
      G_simul[c((n_rates_cum[j]+1) : n_rates_cum[j+1]), c((n_rates_cum[i]+1) : n_rates_cum[i+1])] <- 1 #For symmetry
    }
  }
}
diag(G0_simul) <- rep(0,p0)
diag(G_simul) <- rep(0,p_tot)

#Prior specification
nu_simul <- 5
Psi_simul <- diag(p_tot)

#Simulate data
Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

#We want 3 clusters
K_N_simul <- 2
# phi_star_simul <- rmvnorm(n = KN_simul, mean = c(mu_simul), sigma = solve(Omega_simul))
phi_star_simul <- matrix(0, K_N_simul, p_tot)
phi_star_simul[,c(5,6)] <- log(0.5)
phi_star_simul[1,c(1:4)] <- log(c(0.01,0.9,0.25,0.03))
phi_star_simul[2,c(1:4)] <- log(c(0.7,0.25,0.02,0.01))

lambda_star_simul <- exp(phi_star_simul)
s_simul <- c(rep(1,100), rep(2,50))

#Number of visits per patient (for msm package)
n_times <- 9 + rpois(N, 1)
#Times of visit (0 and Tmax included!)
T_max <- 30
times <- list()
for(i in 1:N){
  n_times_i <- n_times[i]
  times[[i]] <- sort(c(0,runif(n_times_i-2, min = 0, max = T_max), T_max))
}

#Transition rates for Y are covariate- and time-dependent
#Transition rates for processes H are constant (and H are binary processes)
lambda_Y_tilde_simul <- matrix(0,N,dY*(dY - 1))
lambda_H_simul <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_tilde_simul[i,] <- lambda_star_simul[s_simul[i],1:dY*(dY - 1)]
  lambda_H_simul[[i]] <- matrix(lambda_star_simul[s_simul[i],(dY*(dY - 1)+1):p_tot], nrow = p, ncol = 2, byrow = TRUE)
}

#Constant covariates
g <- 2 #Simulate age and sex without intercept (it's the lambda_tilde!)
X <- matrix(0,N,g)
X[,1] <- rnorm(N, mean = 0, sd = 1) #standardized
X[,2] <- (runif(N) <= 0.5)
# #Zero-covariates to try
# X <- matrix(0,N,g)


#Coefficient beta (different for each rate)
beta_simul <- matrix(rnorm(g*dY*(dY - 1), sd = 1),g,dY*(dY - 1),byrow = TRUE)

#Time-varying covariates (not random for msm package, no need for averaging)
q <- 2
Z_msm <- vector("list", length = N)
Z <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i]
  times_i <- times[[i]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm(times_i/T_max, sd = 1)
  Z2 <- rnorm(cos(times_i/T_max*2*pi), sd = 1)
  Z_msm[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z[[i]] <- cbind(Z1,Z2)
  
  
  # #Zero-covariates to try
  # Z1 <- rep(0,n_times_i)
  # Z2 <- rep(0,n_times_i)
  # Z_msm[[i]] <- cbind(Z1,Z2)
  # #For the model we need averages to include in transition probabilities
  # Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  # Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  # Z[[i]] <- cbind(Z1,Z2)
}

#Coefficient gamma (different for each rate)
gamma_simul <- matrix(rnorm(q*dY*(dY - 1), sd = 1),q,dY*(dY - 1),byrow = TRUE)

#Constant part of transition rates
lambda_Y_msm <- matrix(0,N,dY*(dY - 1))
for(i in 1:N){
  n_times_i <- n_times[i]
  lambda_Y_msm[i,] <- lambda_Y_tilde_simul[i,] * exp( X[i,] %*% beta_simul )
}

#Generate the data (...a bit tricky with msm...
#...but we can do it independently and at different times for each process and patient...)

#Response process
Y <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i]
  #No covariates, so use time 1
  Q_Y <- rbind(c(-lambda_Y_msm[i,1], lambda_Y_msm[i,1]), c(lambda_Y_msm[i,2], -lambda_Y_msm[i,2]))
  Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = Z_msm[[i]], beta = gamma_simul, obstimes = times_i, start = 2, mintime = 0)
  # Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
  
  out_aux <- outer(times_i[2:(n_times_i-1)], Yi_msm$times, "<=")
  ind <- apply(out_aux, 1, function(x){which(x)[1]})
  
  Yi <- Yi_msm$states[c(1,ind-1,length(Yi_msm$times))]  - 1 #(0,1)
  
  Y[[i]] <- Yi
}

#Covariate processes
H <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i]
  
  Hi <-vector("list", length = p)
  for(l in 1:p){
    #No covariates
    Q_H <- rbind(c(-lambda_H_simul[[i]][l,1], lambda_H_simul[[i]][l,1]), c(lambda_H_simul[[i]][l,2], -lambda_H_simul[[i]][l,2]))
    Hil_msm <- sim.msm(Q_H, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
    
    #Exclude first and last time points
    out_aux <- outer(times_i[2:(n_times_i-1)], Hil_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    Hil <- Hil_msm$states[c(1,ind-1,length(Hil_msm$times))]  - 1 #(0,1)
    
    Hi[[l]] <- Hil
    
  }
  H[[i]] <- Hi
}

#Interval lengths
epsilon <- lapply(times, diff)

data <- list(Y = Y, H = H, n_rates = n_rates, n_times = n_times, epsilon = epsilon, X = X, Z = Z)

## Now that we have the data, we need to initialize the parameters and priors
n_burn1 <- 1000
n_burn2 <- 5000
n_save <- 1000
thin <- 1
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8
SM_alg <- FALSE

Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold, SM_alg = SM_alg)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

#A-priori expected number of components
Mm <- 5

#Prior on Lambda
update_Lambda <- TRUE
if(update_Lambda){# We provide hyperparameters for Lambda
  #These are based on mean and var for M (number of components)
  mm <- Mm
  ss <- mm * mm/2 #Variance quite small?
  a2 <- mm^2/(ss - mm)
  b2 <- mm/(ss - mm)
  Lambda <- mm
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda, a2 = a2, b2 = b2)
}else{# We fix the values of Lambda
  Lambda <- 1
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda)
}
Param_list <- modifyList(Lambda_list,Param_list)

#Prior on gamma_S
# #Conditionally on M, we have the following relationship between the variance of w and gamma
# gamma_s_grid <- seq(0.001,1,length = 1000)
# plot(0,0,col = "white", xlim = c(0,1), ylim = c(0,0.5))
# for(m in 2:10){
#   lines(gamma_s_grid, (m-1)/(m^2*(gamma_s_grid*M + 1)), col = rainbow(10)[m], lwd = 1.5)
# }
update_gamma_S <- FALSE
if(update_gamma_S){# We provide hyperparameters for gamma_S
  mm <- 0.1
  ss <- 0.1
  a1 <- mm^2/ss
  b1 <- mm/ss
  gamma_S <- mm
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S, a1 = a1, b1 = b1)
}else{# We fix the values of gamma_S
  gamma_S <- 0.1
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S)
}
Param_list <- modifyList(gamma_S_list,Param_list)

#Prior on mu
update_mu <- TRUE
if(update_mu){# We provide hyperparameters for mu
  m_mu <- rep(0,p_tot)
  mu <- m_mu
  mu_list <- list(update_mu = update_mu, mu = mu, m_mu = m_mu)
}else{# We fix the values of mu
  mu <- mu_simul
  mu_list <- list(update_mu = update_mu, mu = mu)
}
Param_list <- modifyList(mu_list,Param_list)

#Prior on m_mu
update_m_mu <- FALSE
if(update_mu){# We provide hyperparameters for m_mu
  Sigma_m_mu <- diag(p_tot)
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu, Sigma_m_mu = Sigma_m_mu)
}else{# We fix the values of m_mu
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu)
}
Param_list <- modifyList(m_mu_list,Param_list)

#Prior on k0
update_k0 <- TRUE
if(update_k0){# We provide hyperparameters for k0
  mm <- 1
  ss <- 1
  a_k0 <- mm^2/ss
  b_k0 <- mm/ss
  k0 <- mm
  k0_list <- list(update_k0 = update_k0, k0 = k0, a_k0 = a_k0, b_k0 = b_k0)
}else{# We fix the values of k0
  k0 <- k0_simul
  k0_list <- list(update_k0 = update_k0, k0 = k0)
}
Param_list <- modifyList(k0_list,Param_list)

#Prior on graph
size_based_prior <- FALSE
if(size_based_prior){
  #The name of the parameters is the same as for d, because it's like marginalizeing wrt d
  a_d <- 1
  b_d <- 1
  G_list <- list(size_based_prior = size_based_prior, a_d = a_d, b_d = b_d)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_d <- TRUE
  if(update_d){
    #Different prior options for d
    d_beta_prior <- FALSE
    if(d_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_d <- mm*(mm*(1 - mm)/ss - 1)
      b_d <- (mm*(1 - mm)^2/ss - (1 - mm))
      d <- mm
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_d = a_d, b_d = b_d, d = d)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      d <- (1 + exp(-mm))^(-1)
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_lambda = a_lambda, d = d)
    }
  }else{# We fix the values of d
    d <- 0.25
    d_list <- list(update_d = update_d, d = d)
  }
  Param_list <- modifyList(d_list,Param_list)
}
Param_list <- modifyList(G_list,Param_list)

#Prior on beta
mu_beta <- matrix(0,g,dY*(dY - 1))
U_beta <- diag(g)
V_beta <- diag(dY)
beta_list <- list(mu_beta = mu_beta, U_beta = U_beta, V_beta = V_beta)
Param_list <- modifyList(beta_list, Param_list)

#Prior on gamma
mu_gamma <- matrix(0,q,dY*(dY - 1))
U_gamma <- diag(q)
V_gamma <- diag(dY)
gamma_list <- list(mu_gamma = mu_gamma, U_gamma = U_gamma, V_gamma = V_gamma)
Param_list <- modifyList(gamma_list, Param_list)

is_DP <- TRUE
is_parametric <- FALSE
if(is_parametric){#We can fit a model for a given clustering
  c_init <- c(1:N)-1
  alpha_list <- list(update_alpha = FALSE, alpha = 0)
}else{
  if(is_DP){#We can fit a DP model
    #Prior on alpha
    # #Based on prior expected number of components
    alpha <- matrix(seq(0.001,2,length = 50000),50000,1)
    sum_alpha <- apply(alpha,1,function(alpha) sum(alpha/(alpha + c(1:N) - 1)))
    alpha_Mm <- alpha[which.min(abs(Mm - sum_alpha))]
    sum_alpha[which.min(abs(Mm - sum_alpha))]
    
    update_alpha <- TRUE
    if(update_alpha){# We provide hyperparameters for alpha
      mm <- alpha_Mm
      ss <- 1
      a_alpha <- mm^2/ss
      b_alpha <- mm/ss
      alpha <- mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha, a_alpha = a_alpha, b_alpha = b_alpha)
    }else{# We fix the values of alpha
      alpha <- alpha_Mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha)
    }
  }else{
    alpha_list <- list(update_alpha = FALSE, alpha = 0)
  }
  c_init <- c(1:N)-1
}
Param_list <- modifyList(alpha_list,Param_list)

model_list <- list(is_parametric = is_parametric, is_DP = is_DP , c_init = c_init)

# Main Gibbs sampler
sample.bdmcmc <- DisProgr( data = data, model_list = model_list, Alg_list = Alg_list, Param_list = Param_list)


# #Save output
# save.image("OUTPUT_DiseaseProgr_Mix.RData")



#Produce plots for report

Graph_list <- sample.bdmcmc$Graph_List

Omega_out <- Graph_list$Omega_out
G_out <- Graph_list$G_out
G0_out <- Graph_list$G0_out

Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

layout(matrix(c(1,2),nrow = 1, ncol = 2))

col_val <- c(Omega_mean, Omega_simul)
col_map <- colorRampPalette(c("white", "red"))(99)

image(Omega_mean, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_mean, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))
image(Omega_simul, col = col_map, breaks = seq(min(col_val), max(col_val), length = 100))
image.plot(Omega_simul, col = col_map, horizontal = TRUE, legend.only = TRUE, breaks = seq(min(col_val), max(col_val), length = 100))

image(G_mean)
image(G_simul)

image(G0_mean)
image(G0_simul)



sum_weights_out <- 1/Graph_list$sum_weights_out
Omega_meanW <- matrix(0,p_tot,p_tot)
G_meanW <- matrix(0,p_tot,p_tot)
G0_meanW <- matrix(0,p0,p0)
for(it in 1:n_save){
  Omega_meanW = Omega_meanW + Omega_out[,,it] * sum_weights_out[it]
  G_meanW = G_meanW + G_out[,,it] * sum_weights_out[it]
  G0_meanW = G0_meanW + G0_out[,,it] * sum_weights_out[it]
}
Omega_meanW <- Omega_meanW / sum(sum_weights_out)
G_meanW <- G_meanW / sum(sum_weights_out)
G0_meanW <- G0_meanW / sum(sum_weights_out)

image(Omega_meanW)
image(Omega_simul)

image(G_meanW)
image(G_simul)

image(G0_meanW)
image(G0_simul)

dev.off()

#d
d_out <- Graph_list$d_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20)
abline(v = d_simul, lwd = 3, col = "red")
plot(d_out, type = "l")

#mu
mu_out <- sample.bdmcmc$mu_out
matplot(t(mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~mu), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#m_mu
m_mu_out <- sample.bdmcmc$m_mu_out
matplot(t(m_mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~m[mu]), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(m_mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), m_mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(m_mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#k0
k0_out <- sample.bdmcmc$k0_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), breaks = 20)
abline(v = k0_simul, lwd = 3, col = "red")
plot(k0_out, type = "l")



#Phi
phi_star_out <- sample.bdmcmc$phi_star_out
phi_all <- rep(0, p_tot)
for(it in 1:n_save){
  phi_all <- rbind(phi_all, phi_star_out[[it]])
}
phi_all <- phi_all[-1,]
matplot(c(1:p_tot), t(phi_all), col = "grey", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~phi), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
matplot(c(1:p_tot), t(phi_star_simul), col = "red", lwd = 1.5, lty = 1, type = "l", add = TRUE)
legend("topright", legend = c("MCMC", "truth"), pch = 19, col = c("grey", "red"), bty = "n", cex = 2)


#beta
beta_out <- sample.bdmcmc$beta_out
beta_df <- data.frame(dim_g = rep(c(1:(g*dY*(dY-1))), n_save), beta_vec = c(beta_out))
beta_true_df <- data.frame( x = c(1:(g*dY*(dY-1))), y = c(beta_simul) ) 

ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
  geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
  labs(x = "g x dY", y = "", title = bquote("Posterior of "~beta)) +
  theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 


beta_out <- sample.bdmcmc$beta_out
matplot(c(1:g), beta_simul, ylim = c(-5,2), col = "red", lwd = 3, lty = 1, type = "l", main = bquote("Posterior of "~beta), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(it in 1:n_save){
  matplot(c(1:g), beta_out[,,it], add = TRUE, col = "grey", lwd = 1, type= "l", lty = 1)
}
for(j1 in 1:g){
  for(j2 in 1:(dY*(dY-1))){
    # boxplot(j1, beta_out[j1,j2,], col = "lightblue", add = TRUE)
    points(j1, beta_simul[j1,j2], col = "red", pch = 19)
  }
}
matplot(c(1:g), beta_simul, col = "red", lwd = 3, lty = 1, type = "l", add = TRUE)
beta_mean <- apply(beta_out,c(1,2), mean)
matplot(c(1:g), beta_mean, col = "yellow", lwd = 1.5, lty = 2, type = "l", add = TRUE)
legend("bottomright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)



#gamma
gamma_out <- sample.bdmcmc$gamma_out
gamma_df <- data.frame(dim_g = rep(c(1:(g*dY*(dY-1))), n_save), gamma_vec = c(gamma_out))
gamma_true_df <- data.frame( x = c(1:(g*dY*(dY-1))), y = c(gamma_simul) ) 

ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
  geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
  labs(x = "g x dY", y = "", title = bquote("Posterior of "~gamma)) +
  theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 

gamma_out <- sample.bdmcmc$gamma_out
matplot(c(1:q), gamma_simul, ylim = c(-2,0.5), col = "red", lwd = 3, lty = 1, type = "l", main = bquote("Posterior of "~gamma), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(it in 1:n_save){
  matplot(c(1:q), gamma_out[,,it], add = TRUE, col = "grey", lwd = 1, type= "l", lty = 1)
}
for(j1 in 1:q){
  for(j2 in 1:(dY*(dY-1))){
    lines(c(j1,j1), quantile(gamma_out[j1,j2,], probs = c(0.025,0.975)), col = "blue", lwd = 2, lty = 1)  
  }
}
matplot(c(1:q), gamma_simul, col = "red", lwd = 3, lty = 1, type = "l", add = TRUE)
gamma_mean <- apply(gamma_out,c(1,2), mean)
matplot(c(1:q), gamma_mean, col = "yellow", lwd = 1.5, lty = 2, type = "l", add = TRUE)
legend("bottomright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#To see how the transition rates are estimated I want to show the evolution of trajectories

phi_star_out <- sample.bdmcmc$phi_star_out
beta_out <- sample.bdmcmc$beta_out
gamma_out <- sample.bdmcmc$gamma_out
c_out <- sample.bdmcmc$c_out

#Reconstruct transition rates for plots
lambda_Y_out <- vector("list", length = n_save)
lambda_H_out <- vector("list", length = n_save)
for(it in 1:n_save){
  lambda_Yi_out <- vector("list", length = N)
  lambda_Hi_out <- vector("list", length = N)
  for(i in 1:N){
    n_times_Y_i <- n_times[i]
    lambda_Yi_out[[i]] <- matrix(exp( phi_star_out[[it]][c_out[it,i],1:2] + X[i,] %*% beta_out[,,it]), nrow = n_times_Y_i-1, ncol = 2, byrow = TRUE) * exp( Z[[i]] %*% gamma_out[,,it] )
    lambda_Hi_out[[i]] <- matrix(exp( phi_star_out[[it]][c_out[it,i],3:p_tot]), nrow = p, ncol = 2, byrow = TRUE)
  }
  lambda_Y_out[[it]] <- lambda_Yi_out
  lambda_H_out[[it]] <- lambda_Hi_out
}

#Select patient
i_plot <- 68
Yi <- Y[[i_plot]]
Hi <- H[[i_plot]]
times_i_plot <- times[[i_plot]]
n_times_i_plot <- n_times[i_plot]

#Simulated transition probabilities
prs_simul_log <- matrix(0,p0,n_times_i_plot)
for(j in 2:n_times_i_plot){
  Yij <- Yi[j-1]
  prs_simul_log[1,j] <- log(lambda_Y_simul[[i_plot]][j-1,Yij+1]) - log(lambda_Y_simul[[i_plot]][j-1,1] + lambda_Y_simul[[i_plot]][j-1,2]) + log( 1 - exp(- ((lambda_Y_simul[[i_plot]][j-1,1] + lambda_Y_simul[[i_plot]][j-1,2]) * epsilon[[i_plot]][j-1])) )
}

for(l in 1:p){
  for(j in 2:n_times_i_plot){
    Hil <- Hi[[l]]
    Hilj <- Hil[j-1]
    prs_simul_log[l+1,j] <- log(lambda_H_simul[[i_plot]][l,Hilj+1]) - log(lambda_H_simul[[i_plot]][l,1] + lambda_H_simul[[i_plot]][l,2]) + log( 1 - exp(- ((lambda_H_simul[[i_plot]][l,1] + lambda_H_simul[[i_plot]][l,2]) * epsilon[[i_plot]][j-1])) )
  }
}


prs_out_log <- array(0, dim = c(p0,n_times_i_plot,n_save))
for(j in 2:n_times_i_plot){
  Yij <- Yi[j-1]
  for(it in 1:n_save){
    prs_out_log[1,j,it] <- log(lambda_Y_out[[it]][[i_plot]][j-1,Yij+1]) - log(lambda_Y_out[[it]][[i_plot]][j-1,1] + lambda_Y_out[[it]][[i_plot]][j-1,2]) + log( 1 - exp(- ((lambda_Y_out[[it]][[i_plot]][j-1,1] + lambda_Y_out[[it]][[i_plot]][j-1,2]) * epsilon[[i_plot]][j-1])) )
  }
}

for(l in 1:p){
  for(j in 2:n_times_i_plot){
    Hil <- Hi[[l]]
    Hilj <- Hil[j-1]
    for(it in 1:n_save){
      prs_out_log[l+1,j,it] <- log(lambda_H_out[[it]][[i_plot]][l,Hilj+1]) - log(lambda_H_out[[it]][[i_plot]][l,1] + lambda_H_out[[it]][[i_plot]][l,2]) + log( 1 - exp(- ((lambda_H_out[[it]][[i_plot]][l,1] + lambda_H_out[[it]][[i_plot]][l,2]) * epsilon[[i_plot]][j-1])) )
    }
  }
}

scale_blue1 <- seq(20, 200, length = p)
scale_blue2 <- seq(120, 255, length = p)
col_proc <- rgb(255/255,51/255,51/255)
for(ip in 1:p){
  col_proc <- c(col_proc, rgb(scale_blue1[ip]/255,scale_blue2[ip]/255,255/255))
}
proc_lab <- c("Y", paste0("H", c(1:p)))

prs_CI_low <- matrix(0,p0,n_times_i_plot)
prs_CI_up <- matrix(0,p0,n_times_i_plot)
for(ip in 1:p0){
  prs_CI_low[ip,] <- apply(prs_out_log[ip,,],1, quantile, probs = 0.025)
  prs_CI_up[ip,] <- apply(prs_out_log[ip,,],1, quantile, probs = 0.975)
}
prs_CI_df <- data.frame(prs_CI_low = c(prs_CI_low), prs_CI_up = c(prs_CI_up), proc = gl(p0, 1, length = p0*n_times_i_plot, labels = proc_lab), times_vec = gl(n_times_i_plot, p0, length = p0*n_times_i_plot, labels = round(times_i_plot,2)) )
prs_simul_df <- data.frame(x_times = c(t(matrix(c(times_i_plot), nrow = n_times_i_plot, ncol = p0))), prs_simul_vec = c(prs_simul_log), proc = gl(p0, 1, length = p0*n_times_i_plot, labels = proc_lab), times_vec = gl(n_times_i_plot, p0, length = p0*n_times_i_plot, labels = round(times_i_plot,2)) )
prs_df <- data.frame(x_times = rep(c(t(matrix(c(times_i_plot), nrow = n_times_i_plot, ncol = p0))), n_save), prs_vec = c(prs_out_log), proc = gl(p0, 1, length = p0*n_times_i_plot*n_save, labels = proc_lab), times_vec = gl(n_times_i_plot, p0, length = p0*n_times_i_plot*n_save, labels = round(times_i_plot,2)) )


ggplot(prs_df, aes(x = x_times, y = prs_vec, col = proc)) + 
  geom_line(data = prs_simul_df, aes(x = x_times, y = prs_simul_vec, col = proc), size = 1.25, linetype = 2) + 
  stat_summary(fun.y = 'mean', geom = 'line', size = 1.25) +
  stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2) + 
  labs(x = "time", y = bquote(log(p[t]^{rs})), title = bquote("Posterior of "~log(p[t]^{rs})~" for subject "~.(i_plot))) +
  theme(text = element_text(size=35), axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 15)) +
  facet_grid(rows = vars(proc), scales="free") +
  # the labels must match what you specified above
  scale_fill_manual(name = proc_lab, values = col_proc) +
  scale_color_manual(name = proc_lab, values = col_proc) +
  theme_bw()





#CLUSTERING

gamma_S_out <- sample.bdmcmc$gamma_S_out
plot(gamma_S_out, type = "l", lwd = 2)

Lambda_out <- sample.bdmcmc$Lambda_out
plot(Lambda_out, type = "l", lwd = 2)


K_N_out <- sample.bdmcmc$K_N_out
layout(matrix(c(1,2),1,2))
plot(K_N_out, type = "l", lwd = 2)
plot(table(K_N_out)/n_save, col = "blue", xlab = bquote(K[N]), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)

M_out <- sample.bdmcmc$M_out
layout(matrix(c(1,2),1,2))
plot(M_out, type = "l", lwd = 2)
plot(table(M_out)/n_save, col = "blue", xlab = bquote(M), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)

M_na_out <- M_out - K_N_out
layout(matrix(c(1,2),1,2))
plot(M_na_out, type = "l", lwd = 2)
plot(table(M_na_out)/n_save, col = "blue", xlab = bquote(M[na]), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)


#Binder loss function equal costs
Const_Binder <- 1/2

c_out <- sample.bdmcmc$c_out

pij <- matrix(0,N,N)
cij <- vector("list", length = n_save)
for(it in 1:n_save){
  cij[[it]] <- outer(c_out[it,], c_out[it,], "==")
  pij <- pij + cij[[it]]
}
pij <- pij/n_save

Binder_f <- rep(0,n_save)
for(it in 1:n_save){
  aux <- (pij - Const_Binder) * as.matrix(cij[[it]])
  aux <-  aux[upper.tri(aux)]
  Binder_f[it] <- sum(aux)
}
plot(Binder_f)
Binder_ind <- which.max(Binder_f)
Binder_out <- c_out[Binder_ind,]

sort_ind <- sort(Binder_out, index.return = TRUE)
sort_ind <- sort_ind$ix

#PPI's sorted by Binder clustering
image(pij[sort_ind,sort_ind], pty="s", col = rev(heat.colors(100)), main = "Posterior Inclusion Probabilities", cex.main = 3, cex.lab = 2, axes = FALSE)
image.plot(pij[sort_ind,sort_ind], pty="s", col = rev(heat.colors(100)), legend.only = TRUE, horizontal = TRUE, axis.args = list(cex.axis = 2), legend.width = 1)
dev.off()

#Plot of posterior distribution of parameters in the clusters
#At each iteration take the mean of rho in each cluster and overlap
phi_star_out <- sample.bdmcmc$phi_star_out
K_Binder <- length(unique(Binder_out))
phi_mean <- array(0,dim = c(K_Binder,p_tot, n_save))
for(it in 1:n_save){
  phi_now <- phi_star_out[[it]]
  for(j in 1:K_Binder){
    if(length(c_out[it,c(1:N)[Binder_out == j]]) > 1){
      phi_mean[j,,it] <- colMeans(phi_now[c_out[it,c(1:N)[Binder_out == j]],])
    }else{
      phi_mean[j,,it] <- phi_now[c_out[it,c(1:N)[Binder_out == j]],]
    }
  }
}

layout(matrix(c(1:K_Binder), 1, K_Binder, byrow = TRUE))
par(mar = c(5,5,5,5))
for(j in 1:K_Binder){
  matplot(c(1:p_tot), phi_mean[j,,], axes = FALSE, ylim = c(min(phi_mean[j,,]),max(phi_mean[j,,])), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote(hat(phi)[.(j)]), xlab = "",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
  lines(c(1:p_tot), rowMeans(phi_mean[j,,]), col = "yellow", lwd = 2, lty = 2)
  axis(2)
  axis(1, labels = c(1:p_tot), at = c(1:p_tot), las = 1)
  matplot(c(1:p_tot), t(phi_star_simul), col = "red", type = "l", lty = 1, lwd = 2, add = TRUE)
}
legend("bottomright", legend = c("truth", "MCMC", "post mean"), pch = 19, col = c("red", "grey", "yellow"), bty = "n", cex = 2)













#################################################################################################################
# TRICKIER EXAMPLE #
#################################################################################################################

## Here we do a fuller simulation from a semi-parametric model where we also have the Markov processes
## We will try to handle directly different time points for different processes and patients
## CLUSTERING from mixture with random number of components

rm(list = ls())
set.seed(123)
library("BayesDiseaseProgr")
library("ggplot2")
library("msm")
library("fields")

######
#Simulate sensible data
######

#Number of patients
N = 150;
#Number of binary covariates observed
p = 4;
#Tot number of porcesses
p0 = p+1;

#Set graph structure(s)
# d_simul <- 0.25
G0_simul <- matrix(0, p0, p0)
# G0_simul[lower.tri(G0_simul)] <- c(runif(p0*(p0-1)/2) <= d_simul)
G0_simul[1,2] <- 1
G0_simul[1,4] <- 1
G0_simul[2,3] <- 1
G0_simul[3,4] <- 1
G0_simul[4,5] <- 1
G0_simul[2,5] <- 1

G0_simul <- G0_simul + t(G0_simul)

#Assume 2-state responses
dY <- 2
n_rates <- c(dY, rep(2,p0-1))
n_rates_cum <- c(0, cumsum(n_rates))
p_tot <- sum(n_rates)
G_simul <- matrix(0, nrow = p_tot, ncol = p_tot)
diag(G0_simul) <- rep(1,p0)
for(i in 1:p0){
  for(j in i:p0){
    if(G0_simul[i,j]){
      G_simul[c((n_rates_cum[i]+1) : n_rates_cum[i+1]), c((n_rates_cum[j]+1) : n_rates_cum[j+1])] <- 1
      G_simul[c((n_rates_cum[j]+1) : n_rates_cum[j+1]), c((n_rates_cum[i]+1) : n_rates_cum[i+1])] <- 1 #For symmetry
    }
  }
}
diag(G0_simul) <- rep(0,p0)
diag(G_simul) <- rep(0,p_tot)

#Prior specification
nu_simul <- 5
Psi_simul <- diag(p_tot)

#Simulate data
Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

#We want 3 clusters
K_N_simul <- 3
# phi_star_simul <- rmvnorm(n = K_N_simul, mean = c(mu_simul), sigma = solve(Omega_simul))
phi_star_simul <- matrix(0, K_N_simul, p_tot)
phi_star_simul[,c(5,6,9,10)] <- log(0.5)
phi_star_simul[1,c(1,2)] <- log(c(0.01,0.9))
phi_star_simul[2,c(1,2)] <- log(c(0.75,0.02))
phi_star_simul[3,c(1,2)] <- log(c(0.3,0.3))
phi_star_simul[1,c(3,4,7,8)] <- log(0.01)
phi_star_simul[2,c(3,4,7,8)] <- log(0.8)
phi_star_simul[3,c(3,4,7,8)] <- log(0.25)

lambda_star_simul <- exp(phi_star_simul)

# matplot(c(1:p_tot), t(rmvnorm(n = 100, mean = c(mu_simul), sigma = solve(Omega_simul))),type = "l", col = "grey", lty = 1)
# matplot(c(1:p_tot), t(phi_star_simul),type = "l", col = "red", lty = 1, add = TRUE)
# matplot(c(1:p_tot), t(lambda_star_simul),type = "l", col = "blue", lty = 1, add = TRUE)


s_simul <- c(rep(1,50), rep(2,50), rep(3,50))
isTRUE(length(s_simul) == N)

#Number of visits per patient (for msm package)
n_times <- 14 + rpois(N, 1)
#Times of visit (0 and Tmax included!)
T_max <- 30
times <- list()
for(i in 1:N){
  n_times_i <- n_times[i]
  times[[i]] <- sort(c(0,runif(n_times_i-2, min = 0, max = T_max), T_max))
}

#Transition rates for Y are covariate- and time-dependent
#Transition rates for processes H are constant (and H are binary processes)
lambda_Y_tilde_simul <- matrix(0,N,dY*(dY - 1))
lambda_H_simul <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_tilde_simul[i,] <- lambda_star_simul[s_simul[i],1:dY*(dY - 1)]
  lambda_H_simul[[i]] <- matrix(lambda_star_simul[s_simul[i],(dY*(dY - 1)+1):p_tot], nrow = p, ncol = 2, byrow = TRUE)
}

#Constant covariates
g <- 2 #Simulate standardized age and sex without intercept (it's the lambda_tilde!)
X <- matrix(0,N,g)
X[,1] <- rnorm(N, mean = -0.1, sd = sqrt(0.1))
X[,2] <- (runif(N) <= 0.5)


#Coefficient beta (different for each rate)
beta_simul <- matrix(rnorm(g*dY*(dY - 1), sd = sqrt(0.1)),g,dY*(dY - 1),byrow = TRUE)

#Time-varying covariates (not random for msm package, no need for averaging)
q <- 2
Z_msm <- vector("list", length = N)
Z <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i]
  times_i <- times[[i]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm(times_i/T_max, sd = sqrt(0.1))
  Z2 <- rnorm(cos(times_i/T_max*2*pi), sd = sqrt(0.1))
  Z_msm[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z[[i]] <- cbind(Z1,Z2)
  
  # #Zero-covariates to try
  # Z1 <- rep(0,n_times_i)
  # Z2 <- rep(0,n_times_i)
  # Z_msm[[i]] <- cbind(Z1,Z2)
  # #For the model we need averages to include in transition probabilities
  # Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  # Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  # Z[[i]] <- cbind(Z1,Z2)
}

#Coefficient gamma (different for each rate)
gamma_simul <- matrix(rnorm(q*dY*(dY - 1), sd = sqrt(0.1)),q,dY*(dY - 1),byrow = TRUE)

#Constant part of transition rates
lambda_Y_msm <- matrix(0,N,dY*(dY - 1))
for(i in 1:N){
  n_times_i <- n_times[i]
  lambda_Y_msm[i,] <- lambda_Y_tilde_simul[i,] * exp( X[i,] %*% beta_simul )
}

#Generate the data (...a bit tricky with msm...
#...but we can do it independently and at different times for each process and patient...)

#Response process
Y <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i]
  
  Q_Y <- rbind(c(-lambda_Y_msm[i,1], lambda_Y_msm[i,1]), c(lambda_Y_msm[i,2], -lambda_Y_msm[i,2]))
  Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = Z_msm[[i]], beta = gamma_simul, obstimes = times_i, start = 2, mintime = 0)
  # Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
  
  out_aux <- outer(times_i[2:(n_times_i-1)], Yi_msm$times, "<=")
  ind <- apply(out_aux, 1, function(x){which(x)[1]})
  
  Yi <- Yi_msm$states[c(1,ind-1,length(Yi_msm$times))]  - 1 #(0,1)
  
  Y[[i]] <- Yi
}

#Covariate processes
H <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i]
  
  Hi <-vector("list", length = p)
  for(l in 1:p){
    #No covariates
    Q_H <- rbind(c(-lambda_H_simul[[i]][l,1], lambda_H_simul[[i]][l,1]), c(lambda_H_simul[[i]][l,2], -lambda_H_simul[[i]][l,2]))
    Hil_msm <- sim.msm(Q_H, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
    
    #Exclude first and last time points
    out_aux <- outer(times_i[2:(n_times_i-1)], Hil_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    Hil <- Hil_msm$states[c(1,ind-1,length(Hil_msm$times))]  - 1 #(0,1)
    
    Hi[[l]] <- Hil
    
  }
  H[[i]] <- Hi
}

#Interval lengths
epsilon <- lapply(times, diff)

data <- list(Y = Y, H = H, n_rates = n_rates, n_times = n_times, epsilon = epsilon, X = X, Z = Z)

## Now that we have the data, we need to initialize the parameters and priors
n_burn1 <- 10
n_burn2 <- 50
n_save <- 10
thin <- 1
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8
SM_alg <- FALSE

Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold, SM_alg = SM_alg)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

Mm <- 5

#Prior on Lambda
update_Lambda <- TRUE
if(update_Lambda){# We provide hyperparameters for Lambda
  #These are based on mean and var for M (number of components)
  mm <- Mm
  ss <- mm * mm/2 #Variance quite small?
  a2 <- mm^2/(ss - mm)
  b2 <- mm/(ss - mm)
  Lambda <- mm
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda, a2 = a2, b2 = b2)
}else{# We fix the values of Lambda
  Lambda <- 1
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda)
}
Param_list <- modifyList(Lambda_list,Param_list)

#Prior on gamma_S
# #Conditionally on M, we have the following relationship between the variance of w and gamma
# gamma_s_grid <- seq(0.001,1,length = 1000)
# plot(0,0,col = "white", xlim = c(0,1), ylim = c(0,0.5))
# for(m in 2:10){
#   lines(gamma_s_grid, (m-1)/(m^2*(gamma_s_grid*M + 1)), col = rainbow(10)[m], lwd = 1.5)
# }
update_gamma_S <- FALSE
if(update_gamma_S){# We provide hyperparameters for gamma_S
  mm <- 0.1
  ss <- 0.1
  a1 <- mm^2/ss
  b1 <- mm/ss
  gamma_S <- mm
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S, a1 = a1, b1 = b1)
}else{# We fix the values of gamma_S
  gamma_S <- 0.5
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S)
}
Param_list <- modifyList(gamma_S_list,Param_list)

#Prior on mu
update_mu <- TRUE
if(update_mu){# We provide hyperparameters for mu
  m_mu <- rep(0,p_tot)
  mu <- m_mu
  mu_list <- list(update_mu = update_mu, mu = mu, m_mu = m_mu)
}else{# We fix the values of mu
  mu <- mu_simul
  mu_list <- list(update_mu = update_mu, mu = mu)
}
Param_list <- modifyList(mu_list,Param_list)

#Prior on m_mu
update_m_mu <- FALSE
if(update_mu){# We provide hyperparameters for m_mu
  Sigma_m_mu <- 0.1*diag(p_tot)
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu, Sigma_m_mu = Sigma_m_mu)
}else{# We fix the values of m_mu
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu)
}
Param_list <- modifyList(m_mu_list,Param_list)

#Prior on k0
update_k0 <- TRUE
if(update_k0){# We provide hyperparameters for k0
  mm <- 1/(nu - 2)
  ss <- 1
  a_k0 <- mm^2/ss
  b_k0 <- mm/ss
  k0 <- mm
  k0_list <- list(update_k0 = update_k0, k0 = k0, a_k0 = a_k0, b_k0 = b_k0)
}else{# We fix the values of k0
  k0 <- k0_simul
  k0_list <- list(update_k0 = update_k0, k0 = k0)
}
Param_list <- modifyList(k0_list,Param_list)

#Prior on graph
size_based_prior <- FALSE
if(size_based_prior){
  #The name of the parameters is the same as for d, because it's like marginalizeing wrt d
  a_d <- 1
  b_d <- 1
  G_list <- list(size_based_prior = size_based_prior, a_d = a_d, b_d = b_d)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_d <- TRUE
  if(update_d){
    #Different prior options for d
    d_beta_prior <- FALSE
    if(d_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_d <- mm*(mm*(1 - mm)/ss - 1)
      b_d <- (mm*(1 - mm)^2/ss - (1 - mm))
      d <- mm
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_d = a_d, b_d = b_d, d = d)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      d <- (1 + exp(-mm))^(-1)
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_lambda = a_lambda, d = d)
    }
  }else{# We fix the values of d
    d <- 0.25
    d_list <- list(update_d = update_d, d = d)
  }
  Param_list <- modifyList(d_list,Param_list)
}
Param_list <- modifyList(G_list,Param_list)


#Prior on beta
mu_beta <- matrix(0,g,dY*(dY - 1))
U_beta <- diag(g)
V_beta <- diag(dY*(dY - 1))
beta_list <- list(mu_beta = mu_beta, U_beta = U_beta, V_beta = V_beta)
Param_list <- modifyList(beta_list, Param_list)

#Prior on gamma
mu_gamma <- matrix(0,q,dY*(dY - 1))
U_gamma <- diag(q)
V_gamma <- diag(dY*(dY - 1))
gamma_list <- list(mu_gamma = mu_gamma, U_gamma = U_gamma, V_gamma = V_gamma)
Param_list <- modifyList(gamma_list, Param_list)

is_DP <- FALSE
is_parametric <- FALSE
if(is_parametric){#We can fit a model for a given clustering
  c_init <- c(1:N)-1
  alpha_list <- list(update_alpha = FALSE, alpha = 0)
}else{
  if(is_DP){#We can fit a DP model
    #Prior on alpha
    # #Based on prior expected number of components
    alpha <- matrix(seq(0.001,2,length = 50000),50000,1)
    sum_alpha <- apply(alpha,1,function(alpha) sum(alpha/(alpha + c(1:N) - 1)))
    alpha_Mm <- alpha[which.min(abs(Mm - sum_alpha))]
    sum_alpha[which.min(abs(Mm - sum_alpha))]
    
    update_alpha <- TRUE
    if(update_alpha){# We provide hyperparameters for alpha
      mm <- alpha_Mm
      ss <- 1
      a_alpha <- mm^2/ss
      b_alpha <- mm/ss
      alpha <- mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha, a_alpha = a_alpha, b_alpha = b_alpha)
    }else{# We fix the values of alpha
      alpha <- alpha_Mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha)
    }
  }else{
    alpha_list <- list(update_alpha = FALSE, alpha = 0)
  }
  c_init <- rep(1,N)-1
}
Param_list <- modifyList(alpha_list,Param_list)

model_list <- list(is_parametric = is_parametric, is_DP = is_DP , c_init = c_init)

# Main Gibbs sampler
sample.bdmcmc <- DisProgr( data = data, model_list = model_list, Alg_list = Alg_list, Param_list = Param_list)



#Produce plots for report

Omega_out <- sample.bdmcmc$Omega_out
G_out <- sample.bdmcmc$G_out
G0_out <- sample.bdmcmc$G0_out

Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

layout(matrix(c(1,2),nrow = 1, ncol = 2))

image(Omega_mean)
image(Omega_simul)

image(G_mean)
image(G_simul)

image(G0_mean)
image(G0_simul)



sum_weights_out <- 1/MCMC_output$sum_weights_out
Omega_meanW <- matrix(0,p_tot,p_tot)
G_meanW <- matrix(0,p_tot,p_tot)
G0_meanW <- matrix(0,p0,p0)
for(it in 1:n_save){
  Omega_meanW = Omega_meanW + Omega_out[,,it] * sum_weights_out[it]
  G_meanW = G_meanW + G_out[,,it] * sum_weights_out[it]
  G0_meanW = G0_meanW + G0_out[,,it] * sum_weights_out[it]
}
Omega_meanW <- Omega_meanW / sum(sum_weights_out)
G_meanW <- G_meanW / sum(sum_weights_out)
G0_meanW <- G0_meanW / sum(sum_weights_out)

image(Omega_meanW)
image(Omega_simul)

image(G_meanW)
image(G_simul)

image(G0_meanW)
image(G0_simul)

dev.off()

#d
d_out <- sample.bdmcmc$d_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20, freq = FALSE)
abline(v = d_simul, lwd = 3, col = "red")
plot(d_out, type = "l")

#mu
mu_out <- sample.bdmcmc$mu_out
matplot(t(mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~mu), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#m_mu
m_mu_out <- sample.bdmcmc$m_mu_out
matplot(t(m_mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~m[mu]), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(m_mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), m_mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(m_mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#k0
k0_out <- sample.bdmcmc$k0_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), breaks = 20)
abline(v = k0_simul, lwd = 3, col = "red")
plot(k0_out, type = "l")



#Phi
phi_star_out <- sample.bdmcmc$phi_star_out
phi_all <- rep(0, p_tot)
for(it in 1:n_save){
  phi_all <- rbind(phi_all, phi_star_out[[it]])
}
phi_all <- phi_all[-1,]
matplot(c(1:p_tot), t(phi_all), col = "grey", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~phi), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
matplot(c(1:p_tot), t(phi_star_simul), col = "red", lwd = 1.5, lty = 1, type = "l", add = TRUE)
legend("topright", legend = c("MCMC", "truth"), pch = 19, col = c("grey", "red"), bty = "n", cex = 2)


#beta
beta_out <- sample.bdmcmc$beta_out
beta_df <- data.frame(dim_g = rep(c(1:(g*dY*(dY-1))), n_save), beta_vec = c(beta_out))
beta_true_df <- data.frame( x = c(1:(g*dY*(dY-1))), y = c(beta_simul) ) 

ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
  geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
  labs(x = "g x dY", y = "", title = bquote("Posterior of "~beta)) +
  theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 


beta_out <- sample.bdmcmc$beta_out
matplot(c(1:g), beta_simul, ylim = c(-5,2), col = "red", lwd = 3, lty = 1, type = "l", main = bquote("Posterior of "~beta), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(it in 1:n_save){
  matplot(c(1:g), beta_out[,,it], add = TRUE, col = "grey", lwd = 1, type= "l", lty = 1)
}
for(j1 in 1:g){
  for(j2 in 1:(dY*(dY-1))){
    # boxplot(j1, beta_out[j1,j2,], col = "lightblue", add = TRUE)
    points(j1, beta_simul[j1,j2], col = "red", pch = 19)
  }
}
matplot(c(1:g), beta_simul, col = "red", lwd = 3, lty = 1, type = "l", add = TRUE)
beta_mean <- apply(beta_out,c(1,2), mean)
matplot(c(1:g), beta_mean, col = "yellow", lwd = 1.5, lty = 2, type = "l", add = TRUE)
legend("bottomright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)



#gamma
gamma_out <- sample.bdmcmc$gamma_out
gamma_df <- data.frame(dim_g = rep(c(1:(g*dY*(dY-1))), n_save), gamma_vec = c(gamma_out))
gamma_true_df <- data.frame( x = c(1:(g*dY*(dY-1))), y = c(gamma_simul) ) 

ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
  geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
  labs(x = "g x dY", y = "", title = bquote("Posterior of "~gamma)) +
  theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 

gamma_out <- sample.bdmcmc$gamma_out
matplot(c(1:q), gamma_simul, ylim = c(-2,0.5), col = "red", lwd = 3, lty = 1, type = "l", main = bquote("Posterior of "~gamma), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(it in 1:n_save){
  matplot(c(1:q), gamma_out[,,it], add = TRUE, col = "grey", lwd = 1, type= "l", lty = 1)
}
for(j1 in 1:q){
  for(j2 in 1:(dY*(dY-1))){
    lines(c(j1,j1), quantile(gamma_out[j1,j2,], probs = c(0.025,0.975)), col = "blue", lwd = 2, lty = 1)  
  }
}
matplot(c(1:q), gamma_simul, col = "red", lwd = 3, lty = 1, type = "l", add = TRUE)
gamma_mean <- apply(gamma_out,c(1,2), mean)
matplot(c(1:q), gamma_mean, col = "yellow", lwd = 1.5, lty = 2, type = "l", add = TRUE)
legend("bottomright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#To see how the transition rates are estimated I want to show the evolution of trajectories

phi_star_out <- sample.bdmcmc$phi_star_out
beta_out <- sample.bdmcmc$beta_out
gamma_out <- sample.bdmcmc$gamma_out
c_out <- sample.bdmcmc$c_out

#Reconstruct transition rates for plots
lambda_Y_out <- vector("list", length = n_save)
lambda_H_out <- vector("list", length = n_save)
for(it in 1:n_save){
  lambda_Yi_out <- vector("list", length = N)
  lambda_Hi_out <- vector("list", length = N)
  for(i in 1:N){
    n_times_Y_i <- n_times[i]
    lambda_Yi_out[[i]] <- matrix(exp( phi_star_out[[it]][c_out[it,i],1:2] + X[i,] %*% beta_out[,,it]), nrow = n_times_Y_i-1, ncol = 2, byrow = TRUE) * exp( Z[[i]] %*% gamma_out[,,it] )
    lambda_Hi_out[[i]] <- matrix(exp( phi_star_out[[it]][c_out[it,i],3:p_tot]), nrow = p, ncol = 2, byrow = TRUE)
  }
  lambda_Y_out[[it]] <- lambda_Yi_out
  lambda_H_out[[it]] <- lambda_Hi_out
}

#Select patient
i_plot <- 68
Yi <- Y[[i_plot]]
Hi <- H[[i_plot]]
times_i_plot <- times[[i_plot]]
n_times_i_plot <- n_times[i_plot]

#Simulated transition probabilities
prs_simul_log <- matrix(0,p0,n_times_i_plot)
for(j in 2:n_times_i_plot){
  Yij <- Yi[j-1]
  prs_simul_log[1,j] <- log(lambda_Y_simul[[i_plot]][j-1,Yij+1]) - log(lambda_Y_simul[[i_plot]][j-1,1] + lambda_Y_simul[[i_plot]][j-1,2]) + log( 1 - exp(- ((lambda_Y_simul[[i_plot]][j-1,1] + lambda_Y_simul[[i_plot]][j-1,2]) * epsilon[[i_plot]][j-1])) )
}

for(l in 1:p){
  for(j in 2:n_times_i_plot){
    Hil <- Hi[[l]]
    Hilj <- Hil[j-1]
    prs_simul_log[l+1,j] <- log(lambda_H_simul[[i_plot]][l,Hilj+1]) - log(lambda_H_simul[[i_plot]][l,1] + lambda_H_simul[[i_plot]][l,2]) + log( 1 - exp(- ((lambda_H_simul[[i_plot]][l,1] + lambda_H_simul[[i_plot]][l,2]) * epsilon[[i_plot]][j-1])) )
  }
}


prs_out_log <- array(0, dim = c(p0,n_times_i_plot,n_save))
for(j in 2:n_times_i_plot){
  Yij <- Yi[j-1]
  for(it in 1:n_save){
    prs_out_log[1,j,it] <- log(lambda_Y_out[[it]][[i_plot]][j-1,Yij+1]) - log(lambda_Y_out[[it]][[i_plot]][j-1,1] + lambda_Y_out[[it]][[i_plot]][j-1,2]) + log( 1 - exp(- ((lambda_Y_out[[it]][[i_plot]][j-1,1] + lambda_Y_out[[it]][[i_plot]][j-1,2]) * epsilon[[i_plot]][j-1])) )
  }
}

for(l in 1:p){
  for(j in 2:n_times_i_plot){
    Hil <- Hi[[l]]
    Hilj <- Hil[j-1]
    for(it in 1:n_save){
      prs_out_log[l+1,j,it] <- log(lambda_H_out[[it]][[i_plot]][l,Hilj+1]) - log(lambda_H_out[[it]][[i_plot]][l,1] + lambda_H_out[[it]][[i_plot]][l,2]) + log( 1 - exp(- ((lambda_H_out[[it]][[i_plot]][l,1] + lambda_H_out[[it]][[i_plot]][l,2]) * epsilon[[i_plot]][j-1])) )
    }
  }
}

scale_blue1 <- seq(20, 200, length = p)
scale_blue2 <- seq(120, 255, length = p)
col_proc <- rgb(255/255,51/255,51/255)
for(ip in 1:p){
  col_proc <- c(col_proc, rgb(scale_blue1[ip]/255,scale_blue2[ip]/255,255/255))
}
proc_lab <- c("Y", paste0("H", c(1:p)))

prs_CI_low <- matrix(0,p0,n_times_i_plot)
prs_CI_up <- matrix(0,p0,n_times_i_plot)
for(ip in 1:p0){
  prs_CI_low[ip,] <- apply(prs_out_log[ip,,],1, quantile, probs = 0.025)
  prs_CI_up[ip,] <- apply(prs_out_log[ip,,],1, quantile, probs = 0.975)
}
prs_CI_df <- data.frame(prs_CI_low = c(prs_CI_low), prs_CI_up = c(prs_CI_up), proc = gl(p0, 1, length = p0*n_times_i_plot, labels = proc_lab), times_vec = gl(n_times_i_plot, p0, length = p0*n_times_i_plot, labels = round(times_i_plot,2)) )
prs_simul_df <- data.frame(x_times = c(t(matrix(c(times_i_plot), nrow = n_times_i_plot, ncol = p0))), prs_simul_vec = c(prs_simul_log), proc = gl(p0, 1, length = p0*n_times_i_plot, labels = proc_lab), times_vec = gl(n_times_i_plot, p0, length = p0*n_times_i_plot, labels = round(times_i_plot,2)) )
prs_df <- data.frame(x_times = rep(c(t(matrix(c(times_i_plot), nrow = n_times_i_plot, ncol = p0))), n_save), prs_vec = c(prs_out_log), proc = gl(p0, 1, length = p0*n_times_i_plot*n_save, labels = proc_lab), times_vec = gl(n_times_i_plot, p0, length = p0*n_times_i_plot*n_save, labels = round(times_i_plot,2)) )


ggplot(prs_df, aes(x = x_times, y = prs_vec, col = proc)) + 
  geom_line(data = prs_simul_df, aes(x = x_times, y = prs_simul_vec, col = proc), size = 1.25, linetype = 2) + 
  stat_summary(fun.y = 'mean', geom = 'line', size = 1.25) +
  stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2) + 
  labs(x = "time", y = bquote(log(p[t]^{rs})), title = bquote("Posterior of "~log(p[t]^{rs})~" for subject "~.(i_plot))) +
  theme(text = element_text(size=35), axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 15)) +
  facet_grid(rows = vars(proc), scales="free") +
  # the labels must match what you specified above
  scale_fill_manual(name = proc_lab, values = col_proc) +
  scale_color_manual(name = proc_lab, values = col_proc) +
  theme_bw()





#CLUSTERING

gamma_S_out <- sample.bdmcmc$gamma_S_out
plot(gamma_S_out, type = "l", lwd = 2)

Lambda_out <- sample.bdmcmc$Lambda_out
plot(Lambda_out, type = "l", lwd = 2)


K_N_out <- sample.bdmcmc$K_N_out
layout(matrix(c(1,2),1,2))
plot(K_N_out, type = "l", lwd = 2)
plot(table(K_N_out)/n_save, col = "blue", xlab = bquote(K[N]), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)

M_out <- sample.bdmcmc$M_out
layout(matrix(c(1,2),1,2))
plot(M_out, type = "l", lwd = 2)
plot(table(M_out)/n_save, col = "blue", xlab = bquote(M), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)

M_na_out <- M_out - K_N_out
layout(matrix(c(1,2),1,2))
plot(M_na_out, type = "l", lwd = 2)
plot(table(M_na_out)/n_save, col = "blue", xlab = bquote(M[na]), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)


#Binder loss function equal costs
Const_Binder <- 1/2

c_out <- sample.bdmcmc$c_out

pij <- matrix(0,N,N)
cij <- vector("list", length = n_save)
for(it in 1:n_save){
  cij[[it]] <- outer(c_out[it,], c_out[it,], "==")
  pij <- pij + cij[[it]]
}
pij <- pij/n_save

Binder_f <- rep(0,n_save)
for(it in 1:n_save){
  aux <- (pij - Const_Binder) * as.matrix(cij[[it]])
  aux <-  aux[upper.tri(aux)]
  Binder_f[it] <- sum(aux)
}
plot(Binder_f)
Binder_ind <- which.max(Binder_f)
Binder_out <- c_out[Binder_ind,]

sort_ind <- sort(Binder_out, index.return = TRUE)
sort_ind <- sort_ind$ix

#PPI's sorted by Binder clustering
image(pij[sort_ind,sort_ind], pty="s", col = rev(heat.colors(100)), main = "Posterior Inclusion Probabilities", cex.main = 3, cex.lab = 2, axes = FALSE)
image.plot(pij[sort_ind,sort_ind], pty="s", col = rev(heat.colors(100)), legend.only = TRUE, horizontal = TRUE, axis.args = list(cex.axis = 2), legend.width = 1)
dev.off()

#Plot of posterior distribution of parameters in the clusters
#At each iteration take the mean of rho in each cluster and overlap
phi_star_out <- sample.bdmcmc$phi_star_out
K_Binder <- length(unique(Binder_out))
phi_mean <- array(0,dim = c(K_Binder,p_tot, n_save))
for(it in 1:n_save){
  phi_now <- phi_star_out[[it]]
  for(j in 1:K_Binder){
    if(length(c_out[it,c(1:N)[Binder_out == j]]) > 1){
      phi_mean[j,,it] <- colMeans(phi_now[c_out[it,c(1:N)[Binder_out == j]],])
    }else{
      phi_mean[j,,it] <- phi_now[c_out[it,c(1:N)[Binder_out == j]],]
    }
  }
}

layout(matrix(c(1:K_Binder), 1, K_Binder, byrow = TRUE))
par(mar = c(5,5,5,5))
for(j in 1:K_Binder){
  matplot(c(1:p_tot), phi_mean[j,,], axes = FALSE, ylim = c(min(phi_mean[j,,]),max(phi_mean[j,,])), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote(hat(phi)[.(j)]), xlab = "",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
  lines(c(1:p_tot), rowMeans(phi_mean[j,,]), col = "yellow", lwd = 2, lty = 2)
  axis(2)
  axis(1, labels = c(1:p_tot), at = c(1:p_tot), las = 1)
  matplot(c(1:p_tot), t(phi_star_simul), col = "red", type = "l", lty = 1, lwd = 2, add = TRUE)
}
legend("bottomright", legend = c("truth", "MCMC", "post mean"), pch = 19, col = c("red", "grey", "yellow"), bty = "n", cex = 2)
























#################################################################################################################
# EXAMPLE FROM MARIA #
#################################################################################################################

## Here we do a fuller simulation from a semi-parametric model where we also have the Markov processes
## We will try to handle directly different time points for different processes and patients
## CLUSTERING from mixture with random number of components

rm(list = ls())
set.seed(123)
library("BayesDiseaseProgr")
library("ggplot2")
library("fields")
library("msm")


######
#Simulate sensible data
######

#Number of patients
N = 150;
#Number of binary covariates observed
p = 2;
#Tot number of porcesses
p0 = p+1;

#Set graph structure(s)
# d_simul <- 0.25
G0_simul <- matrix(0, p0, p0)
G0_simul[1,2] <- 1
# G0_simul[lower.tri(G0_simul)] <- c(runif(p0*(p0-1)/2) <= d_simul)
G0_simul <- G0_simul + t(G0_simul)

#Assume 2-state responses
dY <- 2
n_rates <- c(dY, rep(2,p0-1))
n_rates_cum <- c(0, cumsum(n_rates))
p_tot <- sum(n_rates)
G_simul <- matrix(0, nrow = p_tot, ncol = p_tot)
diag(G0_simul) <- rep(1,p0)
for(i in 1:p0){
  for(j in i:p0){
    if(G0_simul[i,j]){
      G_simul[c((n_rates_cum[i]+1) : n_rates_cum[i+1]), c((n_rates_cum[j]+1) : n_rates_cum[j+1])] <- 1
      G_simul[c((n_rates_cum[j]+1) : n_rates_cum[j+1]), c((n_rates_cum[i]+1) : n_rates_cum[i+1])] <- 1 #For symmetry
    }
  }
}
diag(G0_simul) <- rep(0,p0)
diag(G_simul) <- rep(0,p_tot)

#Prior specification
nu_simul <- 5
Psi_simul <- diag(p_tot)

#Simulate data
Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

#Simulated clusters
K_N_simul <- 2
s_simul <- c(rep(1,100), rep(2,50))

# phi_star_simul <- rmvnorm(n = KN_simul, mean = c(mu_simul), sigma = solve(Omega_simul))
phi_star_simul <- matrix(0, K_N_simul, p_tot)
phi_star_simul[,c(5,6)] <- log(0.5 * 10)
phi_star_simul[1,c(1:4)] <- log(c(0.1,0.9,0.25,0.02) * 10)
phi_star_simul[2,c(1:4)] <- log(c(0.7,0.2,0.02,0.25) * 10)

lambda_star_simul <- exp(phi_star_simul)

matplot(c(1:p_tot), t(rmvnorm(n = 5000, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))), type = "l", lty = 1, col = "grey")
matplot(c(1:p_tot), t(phi_star_simul), type = "l", lty = 1, col = "red", add = TRUE)
lines(c(1:p_tot), mu_simul, col = "blue", lty = 1)


#Number of visits per patient (for msm package)
n_times <- 9 + rpois(N, 1)
#Times of visit (0 and Tmax included!)
T_max <- 30
times <- list()
for(i in 1:N){
  n_times_i <- n_times[i]
  times[[i]] <- sort(c(0,runif(n_times_i-2, min = 0, max = T_max), T_max))
}

#Transition rates for Y are covariate- and time-dependent
#Transition rates for processes H are constant (and H are binary processes)
lambda_Y_tilde_simul <- matrix(0,N,dY*(dY - 1))
lambda_H_simul <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_tilde_simul[i,] <- lambda_star_simul[s_simul[i],1:dY*(dY - 1)]
  lambda_H_simul[[i]] <- matrix(lambda_star_simul[s_simul[i],(dY*(dY - 1)+1):p_tot], nrow = p, ncol = 2, byrow = TRUE)
}

#Constant covariates
g <- 2 #Simulate age and sex without intercept (it's the lambda_tilde!)
X <- matrix(0,N,g)
X[,1] <- rnorm(N, mean = 0, sd = 1) #standardized
X[,2] <- (runif(N) <= 0.5)
# #Zero-covariates to try
# X <- matrix(0,N,g)


#Coefficient beta (different for each rate)
beta_simul <- matrix(rnorm(g*dY*(dY - 1), sd = sqrt(0.5)),g,dY*(dY - 1),byrow = TRUE)

#Time-varying covariates (not random for msm package, no need for averaging)
q <- 2
Z_msm <- vector("list", length = N)
Z <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i]
  times_i <- times[[i]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm(times_i/T_max, sd = 1)
  Z2 <- rnorm(cos(times_i/T_max*2*pi), sd = 1)
  Z_msm[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z[[i]] <- cbind(Z1,Z2)
  
  
  # #Zero-covariates to try
  # Z1 <- rep(0,n_times_i)
  # Z2 <- rep(0,n_times_i)
  # Z_msm[[i]] <- cbind(Z1,Z2)
  # #For the model we need averages to include in transition probabilities
  # Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  # Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  # Z[[i]] <- cbind(Z1,Z2)
}

#Coefficient gamma (different for each rate)
gamma_simul <- matrix(rnorm(q*dY*(dY - 1), sd = sqrt(0.5)),q,dY*(dY - 1),byrow = TRUE)

#Constant part of transition rates
lambda_Y_msm <- matrix(0,N,dY*(dY - 1))
for(i in 1:N){
  n_times_i <- n_times[i]
  lambda_Y_msm[i,] <- lambda_Y_tilde_simul[i,] * exp( X[i,] %*% beta_simul )
}

#Generate the data (...a bit tricky with msm...
#...but we can do it independently and at different times for each process and patient...)

#Response process
Y <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i]
  #No covariates, so use time 1
  Q_Y <- rbind(c(-lambda_Y_msm[i,1], lambda_Y_msm[i,1]), c(lambda_Y_msm[i,2], -lambda_Y_msm[i,2]))
  Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = Z_msm[[i]], beta = gamma_simul, obstimes = times_i, start = 2, mintime = 0)
  # Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
  
  out_aux <- outer(times_i[2:(n_times_i-1)], Yi_msm$times, "<=")
  ind <- apply(out_aux, 1, function(x){which(x)[1]})
  
  Yi <- Yi_msm$states[c(1,ind-1,length(Yi_msm$times))]  - 1 #(0,1)
  
  Y[[i]] <- Yi
}

#Covariate processes
H <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i]
  
  Hi <-vector("list", length = p)
  for(l in 1:p){
    #No covariates
    Q_H <- rbind(c(-lambda_H_simul[[i]][l,1], lambda_H_simul[[i]][l,1]), c(lambda_H_simul[[i]][l,2], -lambda_H_simul[[i]][l,2]))
    Hil_msm <- sim.msm(Q_H, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
    
    #Exclude first and last time points
    out_aux <- outer(times_i[2:(n_times_i-1)], Hil_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    Hil <- Hil_msm$states[c(1,ind-1,length(Hil_msm$times))]  - 1 #(0,1)
    
    Hi[[l]] <- Hil
    
  }
  H[[i]] <- Hi
}

#Interval lengths
epsilon <- lapply(times, diff)

data <- list(Y = Y, H = H, n_rates = n_rates, n_times = n_times, epsilon = epsilon, X = X, Z = Z)

## Now that we have the data, we need to initialize the parameters and priors
#Times and sizes are the same of course
#Run BDMCMC algorithm for GGMs
n_burn1 <- 1000
n_burn2 <- 1500
n_save <- 1000
thin <- 2
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8
SM_alg <- TRUE

Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold, SM_alg = SM_alg)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

#A-priori expected number of components
Mm <- 3

#Prior on Lambda
update_Lambda <- TRUE
if(update_Lambda){# We provide hyperparameters for Lambda
  #These are based on mean and var for M (number of components)
  mm <- Mm
  ss <- mm * 1.01 #Variance quite small?
  a2 <- mm^2/(ss - mm)
  b2 <- mm/(ss - mm)
  Lambda <- mm
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda, a2 = a2, b2 = b2)
}else{# We fix the values of Lambda
  Lambda <- 1
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda)
}
Param_list <- modifyList(Lambda_list,Param_list)

#Prior on gamma_S
update_gamma_S <- FALSE
if(update_gamma_S){# We provide hyperparameters for gamma_S
  mm <- 0.1
  ss <- 0.1
  a1 <- mm^2/ss
  b1 <- mm/ss
  gamma_S <- mm
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S, a1 = a1, b1 = b1)
}else{# We fix the values of gamma_S
  #Find gamma given the number of components M = mm
  mm <- Mm
  #Fix lower bound for var(w_m) <= 1/M(1 - 1/M), close to upper bound because we want variability in weights
  #Divide by expected number of components to "distribute equally"
  low_var <- 1/mm * (1 - 1/mm) / Mm
  gamma_S <- (mm - 1 - mm^2*low_var)/(mm^3 * low_var)
  
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S)
}
Param_list <- modifyList(gamma_S_list,Param_list)

#Prior on mu
update_mu <- TRUE
if(update_mu){# We provide hyperparameters for mu
  m_mu <- rep(0,p_tot)
  mu <- m_mu
  mu_list <- list(update_mu = update_mu, mu = mu, m_mu = m_mu)
}else{# We fix the values of mu
  mu <- mu_simul
  mu_list <- list(update_mu = update_mu, mu = mu)
}
Param_list <- modifyList(mu_list,Param_list)

#Prior on m_mu
update_m_mu <- FALSE
if(update_mu){# We provide hyperparameters for m_mu
  Sigma_m_mu <- 0.1*diag(p_tot)
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu, Sigma_m_mu = Sigma_m_mu)
}else{# We fix the values of m_mu
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu)
}
Param_list <- modifyList(m_mu_list,Param_list)

#Prior on k0
update_k0 <- TRUE
if(update_k0){# We provide hyperparameters for k0
  mm <- 1/(nu - 2)
  ss <- 1
  a_k0 <- mm^2/ss
  b_k0 <- mm/ss
  k0 <- mm
  k0_list <- list(update_k0 = update_k0, k0 = k0, a_k0 = a_k0, b_k0 = b_k0)
}else{# We fix the values of k0
  k0 <- k0_simul
  k0_list <- list(update_k0 = update_k0, k0 = k0)
}
Param_list <- modifyList(k0_list,Param_list)

#Prior on graph
size_based_prior <- FALSE
if(size_based_prior){
  #The name of the parameters is the same as for d, because it's like marginalizeing wrt d
  a_d <- 1
  b_d <- 1
  G_list <- list(size_based_prior = size_based_prior, a_d = a_d, b_d = b_d)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_d <- TRUE
  if(update_d){
    #Different prior options for d
    d_beta_prior <- FALSE
    if(d_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_d <- mm*(mm*(1 - mm)/ss - 1)
      b_d <- (mm*(1 - mm)^2/ss - (1 - mm))
      d <- mm
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_d = a_d, b_d = b_d, d = d)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      d <- (1 + exp(-mm))^(-1)
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_lambda = a_lambda, d = d)
    }
  }else{# We fix the values of d
    d <- 0.25
    d_list <- list(update_d = update_d, d = d)
  }
  Param_list <- modifyList(d_list,Param_list)
}
Param_list <- modifyList(G_list,Param_list)


#Prior on beta
mu_beta <- matrix(0,g,dY*(dY - 1))
U_beta <- diag(g)
V_beta <- diag(dY*(dY - 1))
beta_list <- list(mu_beta = mu_beta, U_beta = U_beta, V_beta = V_beta)
Param_list <- modifyList(beta_list, Param_list)

#Prior on gamma
mu_gamma <- matrix(0,q,dY*(dY - 1))
U_gamma <- diag(q)
V_gamma <- diag(dY*(dY - 1))
gamma_list <- list(mu_gamma = mu_gamma, U_gamma = U_gamma, V_gamma = V_gamma)
Param_list <- modifyList(gamma_list, Param_list)

is_DP <- FALSE
is_parametric <- FALSE
if(is_parametric){#We can fit a model for a given clustering
  c_init <- c(1:N)-1
  alpha_list <- list(update_alpha = FALSE, alpha = 0)
}else{
  if(is_DP){#We can fit a DP model
    #Prior on alpha
    # #Based on prior expected number of components
    alpha <- matrix(seq(0.001,2,length = 50000),50000,1)
    sum_alpha <- apply(alpha,1,function(alpha) sum(alpha/(alpha + c(1:N) - 1)))
    alpha_Mm <- alpha[which.min(abs(Mm - sum_alpha))]
    sum_alpha[which.min(abs(Mm - sum_alpha))]
    
    update_alpha <- TRUE
    if(update_alpha){# We provide hyperparameters for alpha
      mm <- alpha_Mm
      ss <- 1
      a_alpha <- mm^2/ss
      b_alpha <- mm/ss
      alpha <- mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha, a_alpha = a_alpha, b_alpha = b_alpha)
    }else{# We fix the values of alpha
      alpha <- alpha_Mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha)
    }
  }else{
    alpha_list <- list(update_alpha = FALSE, alpha = 0)
  }
  c_init <- c(1:N)-1
}
Param_list <- modifyList(alpha_list,Param_list)

model_list <- list(is_parametric = is_parametric, is_DP = is_DP , c_init = c_init)

# Main Gibbs sampler
sample.bdmcmc <- DisProgr(data = data, model_list = model_list, Alg_list = Alg_list, Param_list = Param_list)


# #Save output
# save.image("OUTPUT_DiseaseProgr_Mix.RData")



#Produce plots for report

Omega_out <- sample.bdmcmc$Omega_out
G_out <- sample.bdmcmc$G_out
G0_out <- sample.bdmcmc$G0_out

Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

layout(matrix(c(1,2),nrow = 1, ncol = 2))

image(Omega_mean)
image(Omega_simul)

image(G_mean)
image(G_simul)

image(G0_mean)
image(G0_simul)



sum_weights_out <- 1/sample.bdmcmc$sum_weights_out
Omega_meanW <- matrix(0,p_tot,p_tot)
G_meanW <- matrix(0,p_tot,p_tot)
G0_meanW <- matrix(0,p0,p0)
for(it in 1:n_save){
  Omega_meanW = Omega_meanW + Omega_out[,,it] * sum_weights_out[it]
  G_meanW = G_meanW + G_out[,,it] * sum_weights_out[it]
  G0_meanW = G0_meanW + G0_out[,,it] * sum_weights_out[it]
}
Omega_meanW <- Omega_meanW / sum(sum_weights_out)
G_meanW <- G_meanW / sum(sum_weights_out)
G0_meanW <- G0_meanW / sum(sum_weights_out)

image(Omega_meanW)
image(Omega_simul)

image(G_meanW)
image(G_simul)

image(G0_meanW)
image(G0_simul)





k0_out <- 1/sample.bdmcmc$k0_out
Omega_mean_k0 <- matrix(0,p_tot,p_tot)
for(it in 1:n_save){
  Omega_mean_k0 = Omega_mean_k0 + Omega_out[,,it] * k0_out[it]
}
Omega_mean_k0 <- Omega_mean_k0 / n_save

image(Omega_mean_k0)
image(k0_simul * Omega_simul)

dev.off()

#d
d_out <- sample.bdmcmc$d_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20)
abline(v = d_simul, lwd = 3, col = "red")
plot(d_out, type = "l")

#mu
mu_out <- sample.bdmcmc$mu_out
matplot(t(mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~mu), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#m_mu
m_mu_out <- sample.bdmcmc$m_mu_out
matplot(t(m_mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~m[mu]), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(m_mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), m_mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(m_mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#k0
k0_out <- sample.bdmcmc$k0_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), breaks = 20)
abline(v = k0_simul, lwd = 3, col = "red")
plot(k0_out, type = "l")
quantile(k0_out)


#Phi
phi_star_out <- sample.bdmcmc$phi_star_out
phi_all <- rep(0, p_tot)
for(it in 1:n_save){
  phi_all <- rbind(phi_all, phi_star_out[[it]])
}
phi_all <- phi_all[-1,]
matplot(c(1:p_tot), t(phi_all), col = "grey", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~phi), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
matplot(c(1:p_tot), t(phi_star_simul), col = "red", lwd = 1.5, lty = 1, type = "l", add = TRUE)
legend("topright", legend = c("MCMC", "truth"), pch = 19, col = c("grey", "red"), bty = "n", cex = 2)


#beta
beta_out <- sample.bdmcmc$beta_out
beta_df <- data.frame(dim_g = rep(c(1:(g*dY*(dY-1))), n_save), beta_vec = c(beta_out))
beta_true_df <- data.frame( x = c(1:(g*dY*(dY-1))), y = c(beta_simul) ) 

ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
  geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
  labs(x = "g x dY", y = "", title = bquote("Posterior of "~beta)) +
  theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 


beta_out <- sample.bdmcmc$beta_out
matplot(c(1:g), beta_simul, ylim = c(-5,2), col = "red", lwd = 3, lty = 1, type = "l", main = bquote("Posterior of "~beta), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(it in 1:n_save){
  matplot(c(1:g), beta_out[,,it], add = TRUE, col = "grey", lwd = 1, type= "l", lty = 1)
}
for(j1 in 1:g){
  for(j2 in 1:(dY*(dY-1))){
    # boxplot(j1, beta_out[j1,j2,], col = "lightblue", add = TRUE)
    points(j1, beta_simul[j1,j2], col = "red", pch = 19)
  }
}
matplot(c(1:g), beta_simul, col = "red", lwd = 3, lty = 1, type = "l", add = TRUE)
beta_mean <- apply(beta_out,c(1,2), mean)
matplot(c(1:g), beta_mean, col = "yellow", lwd = 1.5, lty = 2, type = "l", add = TRUE)
legend("bottomright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)



#gamma
gamma_out <- sample.bdmcmc$gamma_out
gamma_df <- data.frame(dim_g = rep(c(1:(g*dY*(dY-1))), n_save), gamma_vec = c(gamma_out))
gamma_true_df <- data.frame( x = c(1:(g*dY*(dY-1))), y = c(gamma_simul) ) 

ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
  geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
  labs(x = "g x dY", y = "", title = bquote("Posterior of "~gamma)) +
  theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 

gamma_out <- sample.bdmcmc$gamma_out
matplot(c(1:q), gamma_simul, ylim = c(-2,0.5), col = "red", lwd = 3, lty = 1, type = "l", main = bquote("Posterior of "~gamma), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(it in 1:n_save){
  matplot(c(1:q), gamma_out[,,it], add = TRUE, col = "grey", lwd = 1, type= "l", lty = 1)
}
for(j1 in 1:q){
  for(j2 in 1:(dY*(dY-1))){
    lines(c(j1,j1), quantile(gamma_out[j1,j2,], probs = c(0.025,0.975)), col = "blue", lwd = 2, lty = 1)  
  }
}
matplot(c(1:q), gamma_simul, col = "red", lwd = 3, lty = 1, type = "l", add = TRUE)
gamma_mean <- apply(gamma_out,c(1,2), mean)
matplot(c(1:q), gamma_mean, col = "yellow", lwd = 1.5, lty = 2, type = "l", add = TRUE)
legend("bottomright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#To see how the transition rates are estimated I want to show the evolution of trajectories

phi_star_out <- sample.bdmcmc$phi_star_out
beta_out <- sample.bdmcmc$beta_out
gamma_out <- sample.bdmcmc$gamma_out
c_out <- sample.bdmcmc$c_out

#Reconstruct transition rates for plots
lambda_Y_out <- vector("list", length = n_save)
lambda_H_out <- vector("list", length = n_save)
for(it in 1:n_save){
  lambda_Yi_out <- vector("list", length = N)
  lambda_Hi_out <- vector("list", length = N)
  for(i in 1:N){
    n_times_Y_i <- n_times[i]
    lambda_Yi_out[[i]] <- matrix(exp( phi_star_out[[it]][c_out[it,i],1:2] + X[i,] %*% beta_out[,,it]), nrow = n_times_Y_i-1, ncol = 2, byrow = TRUE) * exp( Z[[i]] %*% gamma_out[,,it] )
    lambda_Hi_out[[i]] <- matrix(exp( phi_star_out[[it]][c_out[it,i],3:p_tot]), nrow = p, ncol = 2, byrow = TRUE)
  }
  lambda_Y_out[[it]] <- lambda_Yi_out
  lambda_H_out[[it]] <- lambda_Hi_out
}

#Select patient
i_plot <- 68
Yi <- Y[[i_plot]]
Hi <- H[[i_plot]]
times_i_plot <- times[[i_plot]]
n_times_i_plot <- n_times[i_plot]

#Simulated transition probabilities
prs_simul_log <- matrix(0,p0,n_times_i_plot)
for(j in 2:n_times_i_plot){
  Yij <- Yi[j-1]
  prs_simul_log[1,j] <- log(lambda_Y_simul[[i_plot]][j-1,Yij+1]) - log(lambda_Y_simul[[i_plot]][j-1,1] + lambda_Y_simul[[i_plot]][j-1,2]) + log( 1 - exp(- ((lambda_Y_simul[[i_plot]][j-1,1] + lambda_Y_simul[[i_plot]][j-1,2]) * epsilon[[i_plot]][j-1])) )
}

for(l in 1:p){
  for(j in 2:n_times_i_plot){
    Hil <- Hi[[l]]
    Hilj <- Hil[j-1]
    prs_simul_log[l+1,j] <- log(lambda_H_simul[[i_plot]][l,Hilj+1]) - log(lambda_H_simul[[i_plot]][l,1] + lambda_H_simul[[i_plot]][l,2]) + log( 1 - exp(- ((lambda_H_simul[[i_plot]][l,1] + lambda_H_simul[[i_plot]][l,2]) * epsilon[[i_plot]][j-1])) )
  }
}


prs_out_log <- array(0, dim = c(p0,n_times_i_plot,n_save))
for(j in 2:n_times_i_plot){
  Yij <- Yi[j-1]
  for(it in 1:n_save){
    prs_out_log[1,j,it] <- log(lambda_Y_out[[it]][[i_plot]][j-1,Yij+1]) - log(lambda_Y_out[[it]][[i_plot]][j-1,1] + lambda_Y_out[[it]][[i_plot]][j-1,2]) + log( 1 - exp(- ((lambda_Y_out[[it]][[i_plot]][j-1,1] + lambda_Y_out[[it]][[i_plot]][j-1,2]) * epsilon[[i_plot]][j-1])) )
  }
}

for(l in 1:p){
  for(j in 2:n_times_i_plot){
    Hil <- Hi[[l]]
    Hilj <- Hil[j-1]
    for(it in 1:n_save){
      prs_out_log[l+1,j,it] <- log(lambda_H_out[[it]][[i_plot]][l,Hilj+1]) - log(lambda_H_out[[it]][[i_plot]][l,1] + lambda_H_out[[it]][[i_plot]][l,2]) + log( 1 - exp(- ((lambda_H_out[[it]][[i_plot]][l,1] + lambda_H_out[[it]][[i_plot]][l,2]) * epsilon[[i_plot]][j-1])) )
    }
  }
}

scale_blue1 <- seq(20, 200, length = p)
scale_blue2 <- seq(120, 255, length = p)
col_proc <- rgb(255/255,51/255,51/255)
for(ip in 1:p){
  col_proc <- c(col_proc, rgb(scale_blue1[ip]/255,scale_blue2[ip]/255,255/255))
}
proc_lab <- c("Y", paste0("H", c(1:p)))

prs_CI_low <- matrix(0,p0,n_times_i_plot)
prs_CI_up <- matrix(0,p0,n_times_i_plot)
for(ip in 1:p0){
  prs_CI_low[ip,] <- apply(prs_out_log[ip,,],1, quantile, probs = 0.025)
  prs_CI_up[ip,] <- apply(prs_out_log[ip,,],1, quantile, probs = 0.975)
}
prs_CI_df <- data.frame(prs_CI_low = c(prs_CI_low), prs_CI_up = c(prs_CI_up), proc = gl(p0, 1, length = p0*n_times_i_plot, labels = proc_lab), times_vec = gl(n_times_i_plot, p0, length = p0*n_times_i_plot, labels = round(times_i_plot,2)) )
prs_simul_df <- data.frame(x_times = c(t(matrix(c(times_i_plot), nrow = n_times_i_plot, ncol = p0))), prs_simul_vec = c(prs_simul_log), proc = gl(p0, 1, length = p0*n_times_i_plot, labels = proc_lab), times_vec = gl(n_times_i_plot, p0, length = p0*n_times_i_plot, labels = round(times_i_plot,2)) )
prs_df <- data.frame(x_times = rep(c(t(matrix(c(times_i_plot), nrow = n_times_i_plot, ncol = p0))), n_save), prs_vec = c(prs_out_log), proc = gl(p0, 1, length = p0*n_times_i_plot*n_save, labels = proc_lab), times_vec = gl(n_times_i_plot, p0, length = p0*n_times_i_plot*n_save, labels = round(times_i_plot,2)) )


ggplot(prs_df, aes(x = x_times, y = prs_vec, col = proc)) + 
  geom_line(data = prs_simul_df, aes(x = x_times, y = prs_simul_vec, col = proc), size = 1.25, linetype = 2) + 
  stat_summary(fun.y = 'mean', geom = 'line', size = 1.25) +
  stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2) + 
  labs(x = "time", y = bquote(log(p[t]^{rs})), title = bquote("Posterior of "~log(p[t]^{rs})~" for subject "~.(i_plot))) +
  theme(text = element_text(size=35), axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 15)) +
  facet_grid(rows = vars(proc), scales="free") +
  # the labels must match what you specified above
  scale_fill_manual(name = proc_lab, values = col_proc) +
  scale_color_manual(name = proc_lab, values = col_proc) +
  theme_bw()





#CLUSTERING
c_out <- sample.bdmcmc$c_out
colMeans(c_out)

alpha_out <- sample.bdmcmc$alpha_out
plot(alpha_out, type = "l", lwd = 2)

K_N_out <- sample.bdmcmc$K_N_out
layout(matrix(c(1,2),1,2))
plot(K_N_out, type = "l", lwd = 2)
plot(table(K_N_out)/n_save, col = "blue", xlab = bquote(K[N]), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)

#Binder loss function equal costs
Const_Binder <- 1/2

c_out <- sample.bdmcmc$c_out

pij <- matrix(0,N,N)
cij <- vector("list", length = n_save)
for(it in 1:n_save){
  cij[[it]] <- outer(c_out[it,], c_out[it,], "==")
  pij <- pij + cij[[it]]
}
pij <- pij/n_save

Binder_f <- rep(0,n_save)
for(it in 1:n_save){
  aux <- (pij - Const_Binder) * as.matrix(cij[[it]])
  aux <-  aux[upper.tri(aux)]
  Binder_f[it] <- sum(aux)
}
plot(Binder_f)
Binder_ind <- which.max(Binder_f)
Binder_out <- c_out[Binder_ind,]

sort_ind <- sort(Binder_out, index.return = TRUE)
sort_ind <- sort_ind$ix

#PPI's sorted by Binder clustering
image(pij[sort_ind,sort_ind], pty="s", col = rev(heat.colors(100)), main = "Posterior Inclusion Probabilities", cex.main = 3, cex.lab = 2, axes = FALSE)
image.plot(pij[sort_ind,sort_ind], pty="s", col = rev(heat.colors(100)), legend.only = TRUE, horizontal = TRUE, axis.args = list(cex.axis = 2), legend.width = 1)
dev.off()

#Plot of posterior distribution of parameters in the clusters
#At each iteration take the mean of rho in each cluster and overlap
phi_star_out <- sample.bdmcmc$phi_star_out
K_Binder <- length(unique(Binder_out))
phi_mean <- array(0,dim = c(K_Binder,p_tot, n_save))
for(it in 1:n_save){
  phi_now <- phi_star_out[[it]]
  for(j in 1:K_Binder){
    if(length(c_out[it,c(1:N)[Binder_out == j]]) > 1){
      phi_mean[j,,it] <- colMeans(phi_now[c_out[it,c(1:N)[Binder_out == j]],])
    }else{
      phi_mean[j,,it] <- phi_now[c_out[it,c(1:N)[Binder_out == j]],]
    }
  }
}

layout(matrix(c(1:K_Binder), 1, K_Binder, byrow = TRUE))
par(mar = c(5,5,5,5))
for(j in 1:K_Binder){
  matplot(c(1:p_tot), phi_mean[j,,], axes = FALSE, ylim = c(min(phi_mean[j,,]),max(phi_mean[j,,])), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote(hat(phi)[.(j)]), xlab = "",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
  lines(c(1:p_tot), rowMeans(phi_mean[j,,]), col = "yellow", lwd = 2, lty = 2)
  axis(2)
  axis(1, labels = c(1:p_tot), at = c(1:p_tot), las = 1)
  matplot(c(1:p_tot), t(phi_star_simul), col = "red", type = "l", lty = 1, lwd = 2, add = TRUE)
}
legend("bottomright", legend = c("truth", "MCMC", "post mean"), pch = 19, col = c("red", "grey", "yellow"), bty = "n", cex = 2)































#################################################################################################################
# EXAMPLE with ptot high to test SM algorithm
# also dH0 = 0
#################################################################################################################

rm(list = ls())
set.seed(123)
library("BayesDiseaseProgr")
library("ggplot2")
library("fields")
library("msm")


######
#Simulate sensible data
######

#Number of patients
N = 150;
#Number of binary covariates observed
p = 5;
#Tot number of porcesses
p0 = p+1;

#Set graph structure(s)
d_simul <- 0.05
G0_simul <- matrix(0, p0, p0)
G0_simul[lower.tri(G0_simul)] <- c(runif(p0*(p0-1)/2) <= d_simul)
G0_simul <- G0_simul + t(G0_simul)

#Assume 2-state responses
dY <- 2
n_rates <- c(dY, rep(2,p0-1))
n_rates_cum <- c(0, cumsum(n_rates))
p_tot <- sum(n_rates)
G_simul <- matrix(0, nrow = p_tot, ncol = p_tot)
diag(G0_simul) <- rep(1,p0)
for(i in 1:p0){
  for(j in i:p0){
    if(G0_simul[i,j]){
      G_simul[c((n_rates_cum[i]+1) : n_rates_cum[i+1]), c((n_rates_cum[j]+1) : n_rates_cum[j+1])] <- 1
      G_simul[c((n_rates_cum[j]+1) : n_rates_cum[j+1]), c((n_rates_cum[i]+1) : n_rates_cum[i+1])] <- 1 #For symmetry
    }
  }
}
diag(G0_simul) <- rep(0,p0)
diag(G_simul) <- rep(0,p_tot)

#Prior specification
nu_simul <- 5
Psi_simul <- diag(p_tot)

#Simulate data
Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

#Simulated clusters
K_N_simul <- 2
s_simul <- c(rep(1,100), rep(2,50))

# phi_star_simul <- rmvnorm(n = KN_simul, mean = c(mu_simul), sigma = solve(Omega_simul))
phi_star_simul <- matrix(0, K_N_simul, p_tot)
phi_star_simul[,c(5,6)] <- log(0.5 * 10)
phi_star_simul[1,c(1:4)] <- log(c(0.1,0.9,0.25,0.02) * 10)
phi_star_simul[2,c(1:4)] <- log(c(0.7,0.2,0.02,0.25) * 10)

lambda_star_simul <- exp(phi_star_simul)

matplot(c(1:p_tot), t(rmvnorm(n = 5000, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))), type = "l", lty = 1, col = "grey")
matplot(c(1:p_tot), t(phi_star_simul), type = "l", lty = 1, col = "red", add = TRUE)
lines(c(1:p_tot), mu_simul, col = "blue", lty = 1)


#Number of visits per patient (for msm package)
n_times <- 9 + rpois(N, 1)
#Times of visit (0 and Tmax included!)
T_max <- 30
times <- list()
for(i in 1:N){
  n_times_i <- n_times[i]
  times[[i]] <- sort(c(0,runif(n_times_i-2, min = 0, max = T_max), T_max))
}

#Transition rates for Y are covariate- and time-dependent
#Transition rates for processes H are constant (and H are binary processes)
lambda_Y_tilde_simul <- matrix(0,N,dY*(dY - 1))
lambda_H_simul <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_tilde_simul[i,] <- lambda_star_simul[s_simul[i],1:dY*(dY - 1)]
  lambda_H_simul[[i]] <- matrix(lambda_star_simul[s_simul[i],(dY*(dY - 1)+1):p_tot], nrow = p, ncol = 2, byrow = TRUE)
}

#Constant covariates
g <- 2 #Simulate age and sex without intercept (it's the lambda_tilde!)
X <- matrix(0,N,g)
X[,1] <- rnorm(N, mean = 0, sd = 1) #standardized
X[,2] <- (runif(N) <= 0.5)
# #Zero-covariates to try
# X <- matrix(0,N,g)


#Coefficient beta (different for each rate)
beta_simul <- matrix(rnorm(g*dY*(dY - 1), sd = sqrt(0.5)),g,dY*(dY - 1),byrow = TRUE)

#Time-varying covariates (not random for msm package, no need for averaging)
q <- 2
Z_msm <- vector("list", length = N)
Z <- vector("list", length = N)
for(i in 1:N){
  n_times_i <- n_times[i]
  times_i <- times[[i]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm(times_i/T_max, sd = 1)
  Z2 <- rnorm(cos(times_i/T_max*2*pi), sd = 1)
  Z_msm[[i]] <- cbind(Z1,Z2)
  #For the model we need averages to include in transition probabilities
  Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  Z[[i]] <- cbind(Z1,Z2)
  
  
  # #Zero-covariates to try
  # Z1 <- rep(0,n_times_i)
  # Z2 <- rep(0,n_times_i)
  # Z_msm[[i]] <- cbind(Z1,Z2)
  # #For the model we need averages to include in transition probabilities
  # Z1 <- (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2
  # Z2 <- (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2
  # Z[[i]] <- cbind(Z1,Z2)
}

#Coefficient gamma (different for each rate)
gamma_simul <- matrix(rnorm(q*dY*(dY - 1), sd = sqrt(0.5)),q,dY*(dY - 1),byrow = TRUE)

#Constant part of transition rates
lambda_Y_msm <- matrix(0,N,dY*(dY - 1))
for(i in 1:N){
  n_times_i <- n_times[i]
  lambda_Y_msm[i,] <- lambda_Y_tilde_simul[i,] * exp( X[i,] %*% beta_simul )
}

#Generate the data (...a bit tricky with msm...
#...but we can do it independently and at different times for each process and patient...)

#Response process
Y <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i]
  #No covariates, so use time 1
  Q_Y <- rbind(c(-lambda_Y_msm[i,1], lambda_Y_msm[i,1]), c(lambda_Y_msm[i,2], -lambda_Y_msm[i,2]))
  Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = Z_msm[[i]], beta = gamma_simul, obstimes = times_i, start = 2, mintime = 0)
  # Yi_msm <- sim.msm(Q_Y, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
  
  out_aux <- outer(times_i[2:(n_times_i-1)], Yi_msm$times, "<=")
  ind <- apply(out_aux, 1, function(x){which(x)[1]})
  
  Yi <- Yi_msm$states[c(1,ind-1,length(Yi_msm$times))]  - 1 #(0,1)
  
  Y[[i]] <- Yi
}

#Covariate processes
H <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i]
  
  Hi <-vector("list", length = p)
  for(l in 1:p){
    #No covariates
    Q_H <- rbind(c(-lambda_H_simul[[i]][l,1], lambda_H_simul[[i]][l,1]), c(lambda_H_simul[[i]][l,2], -lambda_H_simul[[i]][l,2]))
    Hil_msm <- sim.msm(Q_H, maxtime = T_max, covs = NULL, beta = NULL, obstimes = times_i, start = 2, mintime = 0)
    
    #Exclude first and last time points
    out_aux <- outer(times_i[2:(n_times_i-1)], Hil_msm$times, "<=")
    ind <- apply(out_aux, 1, function(x){which(x)[1]})
    
    Hil <- Hil_msm$states[c(1,ind-1,length(Hil_msm$times))]  - 1 #(0,1)
    
    Hi[[l]] <- Hil
    
  }
  H[[i]] <- Hi
}


# #Add some NA
# Y[[1]][2]
# Y[[1]][2] <- NA
# Y[[10]][4]
# Y[[10]][4] <- NA
# Y[[100]][6]
# Y[[100]][6] <- NA
# H[[1]][[5]][2]
# H[[1]][[5]][2] <- NA
# H[[2]][[4]][4]
# H[[2]][[4]][4] <- NA
# H[[3]][[3]][6]
# H[[3]][[3]][6] <- NA


#Interval lengths
epsilon <- lapply(times, diff)

data <- list(Y = Y, H = H, n_rates = n_rates, n_times = n_times, epsilon = epsilon, X = X, Z = Z)

## Now that we have the data, we need to initialize the parameters and priors
n_burn1 <- 100
n_burn2 <- 150
n_save <- 100
thin <- 1
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8
SM_alg <- TRUE

Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold, SM_alg = SM_alg)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

#A-priori expected number of components
Mm <- 3

#Prior on Lambda
update_Lambda <- TRUE
if(update_Lambda){# We provide hyperparameters for Lambda
  #These are based on mean and var for M (number of components)
  mm <- Mm
  ss <- mm * 1.01 #Variance quite small?
  a2 <- mm^2/(ss - mm)
  b2 <- mm/(ss - mm)
  Lambda <- mm
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda, a2 = a2, b2 = b2)
}else{# We fix the values of Lambda
  Lambda <- 1
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda)
}
Param_list <- modifyList(Lambda_list,Param_list)

#Prior on gamma_S
update_gamma_S <- FALSE
if(update_gamma_S){# We provide hyperparameters for gamma_S
  mm <- 0.1
  ss <- 0.1
  a1 <- mm^2/ss
  b1 <- mm/ss
  gamma_S <- mm
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S, a1 = a1, b1 = b1)
}else{# We fix the values of gamma_S
  #Find gamma given the number of components M = mm
  mm <- Mm
  #Fix lower bound for var(w_m) <= 1/M(1 - 1/M), close to upper bound because we want variability in weights
  #Divide by expected number of components to "distribute equally"
  low_var <- 1/mm * (1 - 1/mm) / Mm
  gamma_S <- (mm - 1 - mm^2*low_var)/(mm^3 * low_var)
  
  gamma_S_list <- list(update_gamma_S = update_gamma_S, gamma_S = gamma_S)
}
Param_list <- modifyList(gamma_S_list,Param_list)

#Prior on mu
update_mu <- TRUE
if(update_mu){# We provide hyperparameters for mu
  m_mu <- rep(0,p_tot)
  mu <- m_mu
  mu_list <- list(update_mu = update_mu, mu = mu, m_mu = m_mu)
}else{# We fix the values of mu
  mu <- mu_simul
  mu_list <- list(update_mu = update_mu, mu = mu)
}
Param_list <- modifyList(mu_list,Param_list)

#Prior on m_mu
update_m_mu <- FALSE
if(update_mu){# We provide hyperparameters for m_mu
  Sigma_m_mu <- 0.1*diag(p_tot)
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu, Sigma_m_mu = Sigma_m_mu)
}else{# We fix the values of m_mu
  m_mu <- rep(0,p_tot)
  m_mu_list <- list(update_m_mu = update_m_mu, m_mu = m_mu)
}
Param_list <- modifyList(m_mu_list,Param_list)

#Prior on k0
update_k0 <- TRUE
if(update_k0){# We provide hyperparameters for k0
  mm <- 1/(nu - 2)
  ss <- 1
  a_k0 <- mm^2/ss
  b_k0 <- mm/ss
  k0 <- mm
  k0_list <- list(update_k0 = update_k0, k0 = k0, a_k0 = a_k0, b_k0 = b_k0)
}else{# We fix the values of k0
  k0 <- k0_simul
  k0_list <- list(update_k0 = update_k0, k0 = k0)
}
Param_list <- modifyList(k0_list,Param_list)


#Prior on graph
size_based_prior <- FALSE
if(size_based_prior){
  #The name of the parameters is the same as for d, because it's like marginalizeing wrt d
  a_d <- 1
  b_d <- 1
  G_list <- list(size_based_prior = size_based_prior, a_d = a_d, b_d = b_d)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_d <- TRUE
  if(update_d){
    #Different prior options for d
    d_beta_prior <- TRUE
    if(d_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_d <- mm*(mm*(1 - mm)/ss - 1)
      b_d <- (mm*(1 - mm)^2/ss - (1 - mm))
      d <- mm
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_d = a_d, b_d = b_d, d = d)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      d <- (1 + exp(-mm))^(-1)
      d_list <- list(update_d = update_d, d_beta_prior = d_beta_prior, a_lambda = a_lambda, d = d)
    }
  }else{# We fix the values of d
    d <- 0.25
    d_list <- list(update_d = update_d, d = d)
  }
  Param_list <- modifyList(d_list,Param_list)
}
Param_list <- modifyList(G_list,Param_list)


#Prior on beta
mu_beta <- matrix(0,g,dY*(dY - 1))
U_beta <- diag(g)
V_beta <- diag(dY*(dY - 1))
beta_list <- list(mu_beta = mu_beta, U_beta = U_beta, V_beta = V_beta)
Param_list <- modifyList(beta_list, Param_list)

#Prior on gamma
mu_gamma <- matrix(0,q,dY*(dY - 1))
U_gamma <- diag(q)
V_gamma <- diag(dY*(dY - 1))
gamma_list <- list(mu_gamma = mu_gamma, U_gamma = U_gamma, V_gamma = V_gamma)
Param_list <- modifyList(gamma_list, Param_list)

is_DP <- FALSE
is_parametric <- FALSE
if(is_parametric){#We can fit a model for a given clustering
  c_init <- c(1:N)-1
  alpha_list <- list(update_alpha = FALSE, alpha = 0)
}else{
  if(is_DP){#We can fit a DP model
    #Prior on alpha
    # #Based on prior expected number of components
    alpha <- matrix(seq(0.001,2,length = 50000),50000,1)
    sum_alpha <- apply(alpha,1,function(alpha) sum(alpha/(alpha + c(1:N) - 1)))
    alpha_Mm <- alpha[which.min(abs(Mm - sum_alpha))]
    sum_alpha[which.min(abs(Mm - sum_alpha))]
    
    update_alpha <- TRUE
    if(update_alpha){# We provide hyperparameters for alpha
      mm <- alpha_Mm
      ss <- 1
      a_alpha <- mm^2/ss
      b_alpha <- mm/ss
      alpha <- mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha, a_alpha = a_alpha, b_alpha = b_alpha)
    }else{# We fix the values of alpha
      alpha <- alpha_Mm
      alpha_list <- list(update_alpha = update_alpha, alpha = alpha)
    }
  }else{
    alpha_list <- list(update_alpha = FALSE, alpha = 0)
  }
  c_init <- c(1:N)-1
}
Param_list <- modifyList(alpha_list,Param_list)

model_list <- list(is_parametric = is_parametric, is_DP = is_DP , c_init = c_init)

# Main Gibbs sampler
sample.bdmcmc <- DisProgr(data = data, model_list = model_list, Alg_list = Alg_list, Param_list = Param_list)



# #Add some NA
# Y[[1]][2] was 0
# Y[[1]][2] <- NA
# Y[[10]][4] was 0
# Y[[10]][4] <- NA
# Y[[100]][6] was 0
# Y[[100]][6] <- NA
# H[[1]][[5]][2] was 1
# H[[1]][[5]][2] <- NA
# H[[2]][[4]][4] was 1
# H[[2]][[4]][4] <- NA
# H[[3]][[3]][6] was 0
# H[[3]][[3]][6] <- NA

Y_out <- sample.bdmcmc$Y_out
Y_miss <- rep(0,n_save)
for(it in 1:n_save){
  Y_miss[g] <- Y_out[[it]][[1]][2]
}
table(Y_miss)/n_save

H_out <- sample.bdmcmc$H_out
H_miss <- rep(0,n_save)
for(it in 1:n_save){
  H_miss[g] <- H_out[[it]][[3]][[3]][6]
}
table(H_miss)/n_save

# #Save output
# save.image("OUTPUT_DiseaseProgr_Mix.RData")



#Produce plots for report

Omega_out <- sample.bdmcmc$Omega_out
G_out <- sample.bdmcmc$G_out
G0_out <- sample.bdmcmc$G0_out

Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

layout(matrix(c(1,2),nrow = 1, ncol = 2))

image(Omega_mean)
image(Omega_simul)

image(G_mean)
image(G_simul)

image(G0_mean)
image(G0_simul)





k0_out <- 1/sample.bdmcmc$k0_out
Omega_mean_k0 <- matrix(0,p_tot,p_tot)
for(it in 1:n_save){
  Omega_mean_k0 = Omega_mean_k0 + Omega_out[,,it] * k0_out[it]
}
Omega_mean_k0 <- Omega_mean_k0 / n_save

image(Omega_mean_k0)
image(k0_simul * Omega_simul)

dev.off()

#d
d_out <- sample.bdmcmc$d_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20)
abline(v = d_simul, lwd = 3, col = "red")
plot(d_out, type = "l")

#mu
mu_out <- sample.bdmcmc$mu_out
matplot(t(mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~mu), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#m_mu
m_mu_out <- sample.bdmcmc$m_mu_out
matplot(t(m_mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~m[mu]), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(j in 1:p_tot){
  lines(c(j,j), quantile(m_mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
}
lines(c(1:p_tot), m_mu_simul, col = "red", lwd = 3)
lines(c(1:p_tot), colMeans(m_mu_out), col = "yellow", lwd = 2, lty = 2)
par(xpd=TRUE)
legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)


#k0
k0_out <- sample.bdmcmc$k0_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), breaks = 20)
abline(v = k0_simul, lwd = 3, col = "red")
plot(k0_out, type = "l")
quantile(k0_out)


#Phi
phi_star_out <- sample.bdmcmc$phi_star_out
phi_all <- rep(0, p_tot)
for(it in 1:n_save){
  phi_all <- rbind(phi_all, phi_star_out[[it]])
}
phi_all <- phi_all[-1,]
matplot(c(1:p_tot), t(phi_all), col = "grey", lwd = 2, lty = 1, type = "l", main = bquote("Posterior of "~phi), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
matplot(c(1:p_tot), t(phi_star_simul), col = "red", lwd = 1.5, lty = 1, type = "l", add = TRUE)
legend("topright", legend = c("MCMC", "truth"), pch = 19, col = c("grey", "red"), bty = "n", cex = 2)


#beta
beta_out <- sample.bdmcmc$beta_out
beta_df <- data.frame(dim_g = rep(c(1:(g*dY*(dY-1))), n_save), beta_vec = c(beta_out))
beta_true_df <- data.frame( x = c(1:(g*dY*(dY-1))), y = c(beta_simul) ) 

ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
  geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
  labs(x = "g x dY", y = "", title = bquote("Posterior of "~beta)) +
  theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 


beta_out <- sample.bdmcmc$beta_out
matplot(c(1:g), beta_simul, ylim = c(-5,2), col = "red", lwd = 3, lty = 1, type = "l", main = bquote("Posterior of "~beta), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(it in 1:n_save){
  matplot(c(1:g), beta_out[,,it], add = TRUE, col = "grey", lwd = 1, type= "l", lty = 1)
}
for(j1 in 1:g){
  for(j2 in 1:(dY*(dY-1))){
    # boxplot(j1, beta_out[j1,j2,], col = "lightblue", add = TRUE)
    points(j1, beta_simul[j1,j2], col = "red", pch = 19)
  }
}
matplot(c(1:g), beta_simul, col = "red", lwd = 3, lty = 1, type = "l", add = TRUE)
beta_mean <- apply(beta_out,c(1,2), mean)
matplot(c(1:g), beta_mean, col = "yellow", lwd = 1.5, lty = 2, type = "l", add = TRUE)
legend("bottomright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)



#gamma
gamma_out <- sample.bdmcmc$gamma_out
gamma_df <- data.frame(dim_g = rep(c(1:(g*dY*(dY-1))), n_save), gamma_vec = c(gamma_out))
gamma_true_df <- data.frame( x = c(1:(g*dY*(dY-1))), y = c(gamma_simul) ) 

ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
  geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
  labs(x = "g x dY", y = "", title = bquote("Posterior of "~gamma)) +
  theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 

gamma_out <- sample.bdmcmc$gamma_out
matplot(c(1:q), gamma_simul, ylim = c(-2,0.5), col = "red", lwd = 3, lty = 1, type = "l", main = bquote("Posterior of "~gamma), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
for(it in 1:n_save){
  matplot(c(1:q), gamma_out[,,it], add = TRUE, col = "grey", lwd = 1, type= "l", lty = 1)
}
for(j1 in 1:q){
  for(j2 in 1:(dY*(dY-1))){
    lines(c(j1,j1), quantile(gamma_out[j1,j2,], probs = c(0.025,0.975)), col = "blue", lwd = 2, lty = 1)  
  }
}
matplot(c(1:q), gamma_simul, col = "red", lwd = 3, lty = 1, type = "l", add = TRUE)
gamma_mean <- apply(gamma_out,c(1,2), mean)
matplot(c(1:q), gamma_mean, col = "yellow", lwd = 1.5, lty = 2, type = "l", add = TRUE)
legend("bottomright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)



#CLUSTERING
c_out <- sample.bdmcmc$c_out
colMeans(c_out)

alpha_out <- sample.bdmcmc$alpha_out
plot(alpha_out, type = "l", lwd = 2)

K_N_out <- sample.bdmcmc$K_N_out
layout(matrix(c(1,2),1,2))
plot(K_N_out, type = "l", lwd = 2)
plot(table(K_N_out)/n_save, col = "blue", xlab = bquote(K[N]), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)

#Binder loss function equal costs
Const_Binder <- 1/2

c_out <- sample.bdmcmc$c_out

pij <- matrix(0,N,N)
cij <- vector("list", length = n_save)
for(it in 1:n_save){
  cij[[it]] <- outer(c_out[it,], c_out[it,], "==")
  pij <- pij + cij[[it]]
}
pij <- pij/n_save

Binder_f <- rep(0,n_save)
for(it in 1:n_save){
  aux <- (pij - Const_Binder) * as.matrix(cij[[it]])
  aux <-  aux[upper.tri(aux)]
  Binder_f[it] <- sum(aux)
}
plot(Binder_f)
Binder_ind <- which.max(Binder_f)
Binder_out <- c_out[Binder_ind,]

sort_ind <- sort(Binder_out, index.return = TRUE)
sort_ind <- sort_ind$ix

#PPI's sorted by Binder clustering
image(pij[sort_ind,sort_ind], pty="s", col = rev(heat.colors(100)), main = "Posterior Inclusion Probabilities", cex.main = 3, cex.lab = 2, axes = FALSE)
image.plot(pij[sort_ind,sort_ind], pty="s", col = rev(heat.colors(100)), legend.only = TRUE, horizontal = TRUE, axis.args = list(cex.axis = 2), legend.width = 1)
dev.off()

#Plot of posterior distribution of parameters in the clusters
#At each iteration take the mean of rho in each cluster and overlap
phi_star_out <- sample.bdmcmc$phi_star_out
K_Binder <- length(unique(Binder_out))
phi_mean <- array(0,dim = c(K_Binder,p_tot, n_save))
for(it in 1:n_save){
  phi_now <- phi_star_out[[it]]
  for(j in 1:K_Binder){
    if(length(c_out[it,c(1:N)[Binder_out == j]]) > 1){
      phi_mean[j,,it] <- colMeans(phi_now[c_out[it,c(1:N)[Binder_out == j]],])
    }else{
      phi_mean[j,,it] <- phi_now[c_out[it,c(1:N)[Binder_out == j]],]
    }
  }
}

layout(matrix(c(1:K_Binder), 1, K_Binder, byrow = TRUE))
par(mar = c(5,5,5,5))
for(j in 1:K_Binder){
  matplot(c(1:p_tot), phi_mean[j,,], axes = FALSE, ylim = c(min(phi_mean[j,,]),max(phi_mean[j,,])), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote(hat(phi)[.(j)]), xlab = "",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
  lines(c(1:p_tot), rowMeans(phi_mean[j,,]), col = "yellow", lwd = 2, lty = 2)
  axis(2)
  axis(1, labels = c(1:p_tot), at = c(1:p_tot), las = 1)
  matplot(c(1:p_tot), t(phi_star_simul), col = "red", type = "l", lty = 1, lwd = 2, add = TRUE)
}
legend("bottomright", legend = c("truth", "MCMC", "post mean"), pch = 19, col = c("red", "grey", "yellow"), bty = "n", cex = 2)

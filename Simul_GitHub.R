#################################
# SUMS model: simulated example #
#################################

## In this simulated example, we present a simulated example to demonstrate the SUMS model performance and applicability

# Data are sampled using the msm package (Jackson C.H. 2011) from Multi-State processes whose transition rates are covariate-dependent
# We employ simulated covariates of both time-homogeneous and time-varying (continuous) types
# The data are then fitted using the proposed SUMS package
# R commands to produce some summary plots aimed at posterior inference are available

rm(list = ls())
set.seed(123)
library("SUMS")
library("msm")
library("TraMineR")


######
#Simulate sensible data
######

#Number of subjects
N = 150;
#Number of response processes
pY = 2;
#Number of covariate processes
pH = 3;
#Tot number of processes
p0 = pY+pH;

#Set graph structure(s)
G0_simul <- matrix(0, p0, p0)
G0_simul[1,2] <- 1
G0_simul[1,3] <- 1
G0_simul[2,4] <- 1
G0_simul[2,5] <- 1
G0_simul[4,5] <- 1

G0_simul <- G0_simul + t(G0_simul)

#Number of states for each process
dY <- c(2, 3)
dH <- c(2, 3, 2)
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

Omega_simul <- rgwish(nu = nu_simul, Psi = Psi_simul, G = G_simul)
m_mu_simul <- rep(0,p_tot)
k0_simul <- 0.1
mu_simul <- rmvnorm(n = 1, mean = c(m_mu_simul), sigma = solve(k0_simul * Omega_simul))

#Sepcify transition rate for 3 clusters
K_N_simul <- 3
c_simul <- sample.int(K_N_simul, N, replace = TRUE, prob = c(1/3, 1/6, 1/2))
# phi_simul <- rmvnorm(n = K_N_simul, mean = c(mu_simul), sigma = solve(Omega_simul))
phi_simul <- matrix(0, K_N_simul, p_tot)
phi_simul[1,c(1:n_rates_cum[2])] <- log(c(0.1,1.9))
phi_simul[2,c(1:n_rates_cum[2])] <- log(c(1.25,0.2))
phi_simul[3,c(1:n_rates_cum[2])] <- log(c(0.3,1.3))
phi_simul[1,c((n_rates_cum[2]+1):n_rates_cum[3])] <- log(c(2,1.8,1.5,1,0.75,0.5))
phi_simul[2,c((n_rates_cum[2]+1):n_rates_cum[3])] <- log(c(1,0.85,1.15,0.9,1.1,1))
phi_simul[3,c((n_rates_cum[2]+1):n_rates_cum[3])] <- log(c(0.2,0.4,0.75,1.2,1.6,2))
phi_simul[1,c((n_rates_cum[3]+1):n_rates_cum[4])] <- log(0.5)
phi_simul[2,c((n_rates_cum[3]+1):n_rates_cum[4])] <- log(1)
phi_simul[3,c((n_rates_cum[3]+1):n_rates_cum[4])] <- log(2)
phi_simul[1,c((n_rates_cum[4]+1):n_rates_cum[5])] <- log(c(0.5,0.5,0.75,0.75,1.5,1.5))
phi_simul[2,c((n_rates_cum[4]+1):n_rates_cum[5])] <- log(c(0.75,1.25,0.8,2,0.5,1))
phi_simul[3,c((n_rates_cum[4]+1):n_rates_cum[5])] <- log(c(1.75,1.5,1,0.75,0.5,0.45))
phi_simul[1,c((n_rates_cum[5]+1):n_rates_cum[6])] <- log(2.5)
phi_simul[2,c((n_rates_cum[5]+1):n_rates_cum[6])] <- log(1.5)
phi_simul[3,c((n_rates_cum[5]+1):n_rates_cum[6])] <- log(0.5)

matplot(t(phi_simul), type = "l", lty = 1, lwd = 2)
lambda_simul <- exp(phi_simul)


# Select times of observations
# In the following, select only some times of observations (Unif <= 0.8), resembling a "missed visit" pattern that could be observed in real data
T_max <- 30
n_times_simul <- 8
times_simul <- round(seq(0, T_max, length = n_times_simul), 2)

#Times of visit (0 and Tmax included!) sampled among the possible ones (maybe all)
times <- vector("list", length = N)
n_times <- matrix(NA, N, p0)
for(i in 1:N){
  times_i <- list()
  for(ip in 1:p0){
    times_i[[ip]] <- c(0, times_simul[2:n_times_simul][runif(n_times_simul-1) <= 0.8]) 
  }
  n_times[i,] <- sapply(times_i, length)
  times[[i]] <- times_i
}

#Transition rates for Y are covariate- and time-dependent
#Transition rates for processes H are constant
lambda_Y_tilde_simul <- vector("list", length = N)
lambda_H_simul <- vector("list", length = N)
for(i in 1:N){
  lambda_Y_tilde_simul_i <- list()
  for(h in 1:pY){
    lambda_Y_tilde_simul_i[[h]] <- lambda_simul[c_simul[i],(n_rates_cum[h]+1):n_rates_cum[h+1]]
  }
  lambda_Y_tilde_simul[[i]] <- lambda_Y_tilde_simul_i
  
  lambda_H_simul_i <- list()
  for(l in 1:pH){
    lambda_H_simul_i[[l]] <- lambda_simul[c_simul[i],(n_rates_cum[pY+l]+1):n_rates_cum[pY+l+1]]
  }
  lambda_H_simul[[i]] <- lambda_H_simul_i
}

#Constant covariates, but they can differ in each process 1,...,pY
g <- c(2,4)
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
  Z1 <- c(Z1[1], (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2)
  Z2 <- c(Z2[1], (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2)
  Z_h[[i]] <- cbind(Z1,Z2)
}
Z_msm[[h]] <- Z_msm_h
Z[[h]] <- Z_h

h = 2
for(i in 1:N){
  
  n_times_i <- n_times[i,h]
  times_i <- times[[i]][[h]]
  #For sampling we need observed covariates at observed times
  Z1 <- rnorm((times_i/T_max)^2, sd = sqrt(0.5))
  Z2 <- rnorm((times_i/T_max)^3, sd = sqrt(0.5))
  Z3 <- rnorm(exp(times_i/T_max), sd = sqrt(0.5))
  Z_msm_h[[i]] <- cbind(Z1,Z2,Z3)
  #For the model we need averages to include in transition probabilities
  Z1 <- c(Z1[1], (Z1[1:(n_times_i-1)] + Z1[2:n_times_i])/2)
  Z2 <- c(Z2[1], (Z2[1:(n_times_i-1)] + Z2[2:n_times_i])/2)
  Z3 <- c(Z3[1], (Z3[1:(n_times_i-1)] + Z3[2:n_times_i])/2)
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
    lambda_Y_msm_i[[h]] <- lambda_Y_tilde_simul[[i]][[h]] * exp( X[[h]][i,] %*% beta_simul[[h]] )
  }
  lambda_Y_msm[[i]] <- lambda_Y_msm_i
}

#Generate the data (...a bit tricky with msm...but we can do it independently and at different times for each process and subject)

#Simulate response process Y
Y <- vector("list", length = N)
Y_mat <- array(NA, dim = c(pY, N, n_times_simul))
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i,]
  
  lambda_Y_msm_i <- lambda_Y_msm[[i]]
  
  Y_i <- list()
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
    
    count_t <- 0
    for(t in 1:n_times_simul){
      if(times_simul[t] %in% times_i_h ){
        count_t <- count_t + 1
        Y_mat[h,i,t] <- Y_i[[h]][count_t]
      }
    }
  }
  Y[[i]] <- Y_i
}


#Simulate covariate processes H
H <- vector("list", length = N)
H_mat <- array(NA, dim = c(pH, N, n_times_simul))
for(i in 1:N){
  times_i <- times[[i]]
  n_times_i <- n_times[i,]
  
  lambda_H_simul_i <- lambda_H_simul[[i]]
  
  H_i <-vector("list", length = pH)
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
    
    count_t <- 0
    for(t in 1:n_times_simul){
      if(times_simul[t] %in% times_i_l ){
        count_t <- count_t + 1
        H_mat[l,i,t] <- H_i[[l]][count_t]
      }
    }
  } 
  H[[i]] <- H_i
}


#Use optimal matching analysis to compute distance between trajectories
YH_om_sum <- matrix(0, N, N)
for(h in 1:pY){
  n_states_h <- n_states[[1]][h]
  Y_h_alphab <- c(0:(n_states_h-1))
  
  Y_h_seq <- seqdef(Y_mat[h,,], alphabet = Y_h_alphab)
  Y_h_om <- seqdist(Y_h_seq, method = "OM", indel = 1, sm = "TRATE", with.missing = TRUE)
  
  YH_om_sum <- YH_om_sum + Y_h_om
}
for(l in 1:pH){
  n_states_l <- n_states[[2]][l]
  H_l_alphab <- c(0:(n_states_l-1))
  
  H_l_seq <- seqdef(H_mat[l,,], alphabet = H_l_alphab)
  H_l_om <- seqdist(H_l_seq, method = "OM", indel = 1, sm = "TRATE", with.missing = TRUE)
  
  YH_om_sum <- YH_om_sum + H_l_om
}

#Construct matrix of dissimilarities as average over the 5 processes
lambda_dist_mat <- YH_om_sum / p0


#Lengths of the time intervals
epsilon <- vector("list", length = N)
for(i in 1:N){
  times_i <- times[[i]]
  epsilon_i <- vector("list", length = p0)
  for(ip in 1:p0){
    epsilon_i[[ip]] <- diff(times_i[[ip]])
  }
  epsilon[[i]] <- epsilon_i
}

data <- list(Y = Y, H = H, lambda_dist_mat = lambda_dist_mat, n_states = n_states, n_rates = n_rates, n_times = n_times, epsilon = epsilon, X = X, Z = Z)



## Now that we have the data, we need to initialize the algorithm specs, parameters and priors
n_burn1 <- 100
n_burn2 <- 5000
n_save <- 2500
thin <- 2
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8

# For update of allocation variables (Split and Merge algorithm)
update_c <- TRUE
SM_alg <- TRUE
Gibbs_its <- 1
use_Sigma <- TRUE
update_phij <- TRUE
s_SM <- 1
Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold, update_c = update_c, SM_alg = SM_alg, Gibbs_its = Gibbs_its, use_Sigma = use_Sigma, update_phij = update_phij, s_SM = s_SM)

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
  ss <- mm * 1.01 #Variance quite small?
  a2 <- mm^2/(ss - mm)
  b2 <- mm/(ss - mm)
  Lambda <- mm
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda, a2 = a2, b2 = b2)
}else{# We fix the values of Lambda
  Lambda <- 0.1
  Lambda_list <- list(update_Lambda = update_Lambda, Lambda = Lambda)
}
Param_list <- modifyList(Lambda_list,Param_list)

#Prior on gamma_S
update_gamma_S <- TRUE
#Find gamma given the number of components M = mm
mm <- Mm
if(update_gamma_S){# We provide hyperparameters for gamma_S
  mm <- 1 - 1/mm
  ss <- 1
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
  
  update_d <- FALSE
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
    d <- 0.1
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
set.seed(8)
MCMC_output <- DisProgr( data = data, model_list = model_list, Alg_list = Alg_list, Param_list = Param_list)


# #Save output
# save(MCMC_output, file = "OUTPUT_SUMS.RData")






##############
# Some plots #
#''''''''''''#

library("fields")

Graph_list <- MCMC_output$Graph_List

Omega_out <- Graph_list$Omega_out
G_out <- Graph_list$G_out
G0_out <- Graph_list$G0_out

# Precision matrix
Omega_mean <- apply(Omega_out, c(1,2), mean)
G_mean <- apply(G_out, c(1,2), mean)
G0_mean <- apply(G0_out, c(1,2), mean)

col_val <- c(Omega_mean, Omega_simul)
col_map <- colorRampPalette(c("white", "red"))(99)
breaks_map <- seq(min(col_val), max(col_val), length = 100)

layout(matrix(c(1:2),nrow = 1, ncol = 2))
par(mar = c(10.5,0.5,10.5,0.5))
image(Omega_mean, col = col_map, breaks = breaks_map, axes = FALSE, main = "Estimated", cex.main = 3)
image.plot(Omega_mean, col = col_map, breaks = breaks_map, horizontal = TRUE, lengend.only = TRUE, add = TRUE)
image(Omega_simul, col = col_map, breaks = breaks_map, axes = FALSE, main = "Simulated", cex.main = 3)
image.plot(Omega_simul, col = col_map, breaks = breaks_map, horizontal = TRUE, lengend.only = TRUE, add = TRUE)

# Graph G
col_val <- c(G_mean, G_simul)
col_map <- colorRampPalette(c("white", "red"))(99)
breaks_map <- seq(min(col_val), max(col_val), length = 100)

layout(matrix(c(1:2),nrow = 1, ncol = 2))
par(mar = c(10.5,0.5,10.5,0.5))
image(G_mean, col = col_map, breaks = breaks_map, axes = FALSE, main = "Estimated", cex.main = 3)
image.plot(G_mean, col = col_map, breaks = breaks_map, horizontal = TRUE, lengend.only = TRUE, add = TRUE)
image(G_simul, col = col_map, breaks = breaks_map, axes = FALSE, main = "Simulated", cex.main = 3)
image.plot(G_simul, col = col_map, breaks = breaks_map, horizontal = TRUE, lengend.only = TRUE, add = TRUE)

# Graph G0
col_val <- c(G0_mean, G0_simul)
col_map <- colorRampPalette(c("white", "red"))(99)
breaks_map <- seq(min(col_val), max(col_val), length = 100)

layout(matrix(c(1:2),nrow = 1, ncol = 2))
par(mar = c(10.5,0.5,10.5,0.5))
image(G0_mean, col = col_map, breaks = breaks_map, axes = FALSE, main = "Estimated", cex.main = 3)
image.plot(G0_mean, col = col_map, breaks = breaks_map, horizontal = TRUE, lengend.only = TRUE, add = TRUE)
image(G0_simul, col = col_map, breaks = breaks_map, axes = FALSE, main = "Simulated", cex.main = 3)
image.plot(G0_simul, col = col_map, breaks = breaks_map, horizontal = TRUE, lengend.only = TRUE, add = TRUE)







dev.off()
# Posterior distribution of d (if random)
if(update_d){
  d_out <- Graph_list$d_out
  layout(matrix(c(1,2),nrow = 1, ncol = 2))
  hist(d_out, col = "lightblue", main = bquote("Posterior of "~d), breaks = 20)
  plot(d_out, type = "l")
  abline(h = sum(G0_simul)/(2*p0*(p0-1)), lwd = 3, col = "red")
  dev.off()
}

# Posterior distribution of mu (if random)
if(update_mu){
  mu_out <- MCMC_output$mu_out
  matplot(t(mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~mu), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
  for(j in 1:p_tot){
    lines(c(j,j), quantile(mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
  }
  # lines(c(1:p_tot), mu_simul, col = "red", lwd = 3)
  lines(c(1:p_tot), colMeans(mu_out), col = "yellow", lwd = 2, lty = 2)
  par(xpd=TRUE)
  legend("topright", legend = c("MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)
}


# Posterior distribution of m_mu (if random)
if(update_m_mu){
  m_mu_out <- MCMC_output$m_mu_out
  matplot(t(m_mu_out), type = "l", lty = 1, col = "grey", lwd = 0.75, main = bquote("Posterior samples of "~m[mu]), xlab = "dim",  ylab = "", cex.main = 3, cex.lab = 2, cex.axis = 2)
  for(j in 1:p_tot){
    lines(c(j,j), quantile(m_mu_out[,j], probs = c(0.025, 0.975)), col = "blue", lwd = 2, lty = 1)
  }
  lines(c(1:p_tot), m_mu_simul, col = "red", lwd = 3)
  lines(c(1:p_tot), colMeans(m_mu_out), col = "yellow", lwd = 2, lty = 2)
  par(xpd=TRUE)
  legend("topright", legend = c("truth", "MCMC", "post mean", "95% CI"), pch = 19, col = c("red", "grey", "yellow", "blue"), bty = "n", cex = 2)
}


# Posterior distribution of k0 (if random)
if(update_k0){
  k0_out <- MCMC_output$k0_out
  layout(matrix(c(1,2),nrow = 1, ncol = 2))
  hist(k0_out, col = "lightblue", main = bquote("Posterior of "~k[0]), xlab = "", breaks = 20, freq = FALSE, cex.main = 3, cex.lab = 2, cex.axis = 2)
  plot(k0_out, type = "l")
  abline(h = k0_simul, lwd = 3, col = "red")
}



## Posterior distribution of regression coefficients beta
library("ggplot2")

beta_out <- MCMC_output$XZ_List$beta_out
for(h in 1:pY){
  beta_out_h <- array(0, dim = c(g[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    beta_out_h[,,it] <- beta_out[[it]][[h]]
  }
  
  beta_df <- data.frame(dim_g = rep(c(1:(g[h]*dY[h]*(dY[h] - 1))), n_save), beta_vec = c(beta_out_h))
  beta_true_df <- data.frame( x = c(1:(g[h]*dY[h]*(dY[h] - 1))), y = c(beta_simul[[h]]) ) 
  
  print(
    ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
      geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
      labs(x = "g^h x dY^h * (dY^h - 1)", y = "", title = bquote("Posterior of "~beta)) +
      theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
  )
}


## Posterior distribution of regression coefficients gamma
gamma_out <- MCMC_output$XZ_List$gamma_out
for(h in 1:pY){
  gamma_out_h <- array(0, dim = c(q[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    gamma_out_h[,,it] <- gamma_out[[it]][[h]]
  }
  
  gamma_df <- data.frame(dim_g = rep(c(1:(q[h]*dY[h]*(dY[h] - 1))), n_save), gamma_vec = c(gamma_out_h))
  gamma_true_df <- data.frame( x = c(1:(q[h]*dY[h]*(dY[h] - 1))), y = c(gamma_simul[[h]]) ) 
  
  print(
    ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
      geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
      labs(x = "q^h x dY^h * (dY^h - 1)", y = "", title = bquote("Posterior of "~gamma)) +
      theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
  )
}



#CLUSTERING

# Posterior distribution of gamma_S (if random)
if(update_gamma_S){
  gamma_S_out <- MCMC_output$gamma_S_out
  layout(matrix(c(1,2),nrow = 1, ncol = 2))
  hist(gamma_S_out, col = "lightblue", main = bquote("Posterior of "~gamma[S]), breaks = 20)
  plot(gamma_S_out, type = "l")
}

# Posterior distribution of Lambda (if random)
if(update_Lambda){
  Lambda_out <- MCMC_output$Lambda_out
  layout(matrix(c(1,2),nrow = 1, ncol = 2))
  hist(Lambda_out, col = "lightblue", main = bquote("Posterior of "~Lambda), breaks = 20)
  plot(Lambda_out, type = "l")
}

# Posterior distribution of K_N
K_N_out <- MCMC_output$K_N_out
plot(K_N_out, type = "l", lwd = 2, xlab = bquote(K[N]), ylab = "", main = "")
plot(table(K_N_out)/n_save, lwd = 2, col = "blue", xlab = bquote(K[N]), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)

# Posterior distribution of M
M_out <- MCMC_output$M_out
plot(M_out, type = "l", lwd = 2, xlab = bquote(M), ylab = "", main = "")
plot(table(M_out)/n_save, lwd = 2, col = "blue", xlab = bquote(M), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)

# Posterior distribution of M_na
M_na_out <- M_out - K_N_out
plot(M_na_out, type = "l", lwd = 2, xlab = bquote(M[na]), ylab = "", main = "")
plot(table(M_na_out)/n_save, lwd = 2, col = "blue", xlab = bquote(M[na]), ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main = 3)

#Plot of entropy
Entropy <- MCMC_output$Binder_List$Entropy_out
layout(matrix(c(1,2),nrow = 1, ncol = 2))
hist(Entropy, col = "lightblue", main = bquote("Entropy"), breaks = 20)
plot(Entropy, type = "l")


#Binder loss function equal costs
pij <- MCMC_output$Binder_List$pij
pij <- (pij + t(pij))
diag(pij) <- 1

Binder_f <- MCMC_output$Binder_List$Binder_f
Binder_out <- c(MCMC_output$Binder_List$Binder_est)

sort_ind <- sort(Binder_out, index.return = TRUE)
sort_ind <- sort_ind$ix


#PPI's sorted by Binder clustering
image(pij[sort_ind,sort_ind], pty="s", col = rev(heat.colors(100)), main = "", cex.main = 3, cex.lab = 2, axes = FALSE)
image.plot(pij[sort_ind,sort_ind], pty="s", col = rev(heat.colors(100)), legend.only = TRUE, horizontal = TRUE, axis.args = list(cex.axis = 2), legend.width = 1)



# To check for misclassification
dev.off()
library('gplots')

c_simul_col <- rep(0, N)
c_simul_col[c_simul == 1] <- rgb(0/255, 204/255, 102/255)
c_simul_col[c_simul == 2] <- rgb(127/255, 0/255, 255/255)
c_simul_col[c_simul == 3] <- rgb(0/255, 128/255, 255/255)

c_simul_col <- c_simul_col[sort_ind]
pij_sorted <- pij[sort_ind,sort_ind]
heatmap.2(pij_sorted, dendrogram = "none", labRow = FALSE, labCol = FALSE, Rowv = FALSE, Colv = FALSE, ColSideColors = c_simul_col, RowSideColors = c_simul_col, na.rm = TRUE, key = FALSE, main = "", trace="none", lwid=c(1,5), lhei=c(1,5), col = rev(heat.colors(100)))
dev.off()






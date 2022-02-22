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
library("truncnorm")

######
#Simulate sensible data
######

#Number of subjects
N = 200;
#Number of response processes
pY = 1;
#Number of covariate processes
pH = 2;
#Tot number of processes
p0 = pY+pH;

#Set graph structure(s)
G0_simul <- matrix(0, p0, p0)
G0_simul[1,2] <- 1
G0_simul[1,3] <- 1

G0_simul <- G0_simul + t(G0_simul)

#Number of states for each process
dY <- 2
dH <- c(2, 2)
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

#Simulate phi's from the prior (2 clusters)
Omega_simul <- diag(p_tot) * 10
Omega_simul[G_simul == 1] <- 2
chol(Omega_simul)

aux <- c(0.12,0.37,0.11,0.26,0.21,0.34)
mu_simul <- log(aux)

#Sepcify transition rate for 3 clusters
K_N_simul <- 2
c_simul <- sample.int(K_N_simul, N, replace = TRUE, prob = c(1/3, 2/3))
phi_simul <- matrix(0, K_N_simul, p_tot)
phi_simul[1,] <- rmvnorm(n = 1, mean = c(mu_simul), sigma = solve(Omega_simul))
phi_simul[2,] <- rmvnorm(n = 1, mean = c(mu_simul) + 1, sigma = solve(Omega_simul))

matplot(t(phi_simul), type = "l", lty = 1, lwd = 2)
lambda_simul <- exp(phi_simul)
matplot(t(lambda_simul), type = "l", lty = 1, lwd = 2)


# Select times of observations (0 included!)
T_max <- 10
t_init <- 0

#SAME FOR ALL PROCESSES (for fitting of De Iorio 2018)
#This can be generalised
times <- vector("list", length = N)
n_times <- matrix(NA, N, p0)
for(i in 1:N){
  times_i <- list()
  times_i_1 <- t_init
  time_now <- t_init
  while(time_now <= T_max){
    time_now <- time_now + rtruncnorm(1, a=0.50, b=Inf, mean = 1, sd = 1)
    times_i_1 <- c(times_i_1, time_now)
  }
  times_i_1 <- times_i_1[-length(times_i_1)]
  
  times_i[[1]] <- times_i_1
  
  # Same time points for the time-varying categorical covariates (although we sample them from msm)
  # One binary covariate for each process in this simple example
  times_i[[2]] <- times_i[[1]]
  times_i[[3]] <- times_i[[1]]
  
  n_times[i,] <- sapply(times_i, length)
  times[[i]] <- times_i
}
mean(n_times)


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

#Constant covariates
# They cannot differ between process in De Iorio 2018...
g <- 2
X <- list()

aux <- matrix(1,N,g[1])
aux[,1] <- rnorm(N, mean = 0, sd = 1) #standardized
aux[,2] <- (runif(N) <= 0.5)

X[[1]] <- aux


#Coefficient beta
beta_simul <- vector("list", length = pY)
beta_simul[[1]] <- matrix(c(1,1,-1,-1), g[1], dY[1]*(dY[1] - 1),byrow = TRUE)

#Time-varying covariates (not random for msm package, no need for averaging)
# They cannot differ between process in De Iorio 2018...
q <- 2
Z_msm <- vector("list", length = pY)
Z <- vector("list", length = pY)

#Careful with dimension q!
for(h in 1:pY){
  Z_msm_h <- vector("list", length = N)
  Z_h <- vector("list", length = N)
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
}


#Coefficient gamma
gamma_simul <- vector("list", length = pY)
gamma_simul[[1]] <- matrix(c(0.75,-1.25,1.25,-0.75), q[1], dY[1]*(dY[1] - 1),byrow = TRUE)

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
    
  }
  Y[[i]] <- Y_i
}


#Simulate covariate processes H
H <- vector("list", length = N)
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
    
  } 
  H[[i]] <- H_i
}


#Construct matrix of dissimilarities as average over the 5 processes
lambda_dist_mat <- matrix(0, N, N)


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
n_burn1 <- 10
n_burn2 <- 10
n_save <- 10
thin <- 1
n_edges <- 1 #Lets' leave it like this!!!
threshold <- 1e-8

# For update of allocation variables
update_c <- TRUE
Alg_list <- list(n_save = n_save, n_burn1 = n_burn1, n_burn2 = n_burn2, thin = thin, n_edges = n_edges, threshold = threshold, update_c = update_c)

#Parameters and hyperparameters
nu = nu_simul
Psi = Psi_simul
Param_list <- list(nu = nu, Psi = Psi)

#A-priori expected number of components
Mm <- 5

#Prior on Lambda
update_Lambda <- FALSE
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
update_gamma_S <- FALSE
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
  #The name of the parameters is the same as for eta, because it's like marginalizing wrt d
  a_eta <- 1
  b_eta <- 1
  G_list <- list(size_based_prior = size_based_prior, a_eta = a_eta, b_eta = b_eta)
}else{#Prior on d
  G_list <- list(size_based_prior = size_based_prior)
  
  update_eta <- FALSE
  if(update_eta){
    #Different prior options for eta
    eta_beta_prior <- TRUE
    if(eta_beta_prior){
      mm <- 0.5
      ss <- mm * (1 - mm) / 2 #non-informative enough?
      a_eta <- mm*(mm*(1 - mm)/ss - 1)
      b_eta <- (mm*(1 - mm)^2/ss - (1 - mm))
      eta <- mm
      eta_list <- list(update_eta = update_eta, eta_beta_prior = eta_beta_prior, a_eta = a_eta, b_eta = b_eta, eta = eta)
    }else{
      a_lambda <- 1
      mm <- 0 #Mean of normal
      eta <- (1 + exp(-mm))^(-1)
      eta_list <- list(update_eta = update_eta, eta_beta_prior = eta_beta_prior, a_lambda = a_lambda, eta = eta)
    }
  }else{# We fix the values of eta
    eta <- 0.1
    eta_list <- list(update_eta = update_eta, eta = eta)
  }
  Param_list <- modifyList(eta_list,Param_list)
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


#Save output
save(MCMC_output, file = "OUTPUT_SUMS.RData")






##############
# Some plots #
#''''''''''''#

proc_names_Y <- paste0("Y", c(1:pY))
proc_names_H <- paste0("H", c(1:pH))



#Produce plots for comaprison section of the paper
library("fields")
library("igraph")
library("gplots")
library("ggplot2")



#####################
#CLUSTERING
Binder_List <- MCMC_output$Binder_List


K_N_out <- MCMC_output$K_N_out
M_out <- MCMC_output$M_out
postscript("K_N_M_hist.eps")
par(mar = c(5,7,10,7), mgp=c(5.5,2,0))
plot(0, 0, lwd = 5, col = "white", xlim = c(0.5,5.5), ylim = c(0,1), xlab = "", ylab = "", main = "", cex.axis = 3, cex.lab = 4, cex.main = 3, axes = FALSE)
axis(2, labels = c(0,0.4,0.8,1), at = c(0,0.4,0.8,1), cex.axis = 1.5)
axis(1, labels = c(1:5), at = c(1:5), cex.axis = 2)

t1 <- table(K_N_out)/n_save
t2 <- table(M_out)/n_save
segments(x0 = sort(unique(K_N_out) - 0.05), x1 = sort(unique(K_N_out) - 0.05), y0 = rep(0,length(unique(K_N_out))), y1 = as.vector(t1), lwd = 5, col = "blue")
segments(x0 = sort(unique(M_out) + 0.05), x1 = sort(unique(M_out) + 0.05), y0 = rep(0,length(unique(M_out))), y1 = as.vector(t2), lwd = 5, col = "red")
legend("topright", legend = c(expression(K[N]), "M"), cex = 2, bty = "n", col = c("blue", "red"), pch = 19)
dev.off()




Graph_list <- MCMC_output$Graph_List
sum_weights_out <- 1/Graph_list$sum_weights_out
norm_weights_out <- sum_weights_out / sum(sum_weights_out)


#Binder loss function equal costs
Const_Binder <- 1/2

c_out <- MCMC_output$c_out + 1
pij <- matrix(0,N,N)
for(it in 1:n_save){
  pij <- pij + outer(c_out[it,], c_out[it,], "==") * norm_weights_out[it]
}

Binder_f <- rep(0,n_save)
for(it in 1:n_save){
  cij <- outer(c_out[it,], c_out[it,], "==")
  aux <- (pij - Const_Binder) * as.matrix(cij)
  aux <-  aux[upper.tri(aux)]
  Binder_f[it] <- sum(aux)
}
Binder_ind <- which.max(Binder_f)
Binder_out <- c_out[Binder_ind,]
sort_ind <- sort(Binder_out, index.return = TRUE)
sort_ind <- sort_ind$ix


pij_sorted <- pij[sort_ind,sort_ind]

K_Binder <- length(unique(Binder_out))

Cluster_col <- c(rgb(0/255, 204/255, 102/255), rgb(127/255, 0/255, 255/255), rgb(0/255, 128/255, 255/255))

c_simul_col <- rep(0, N)
for(j in 1:K_N_simul){
  c_simul_col[c_simul == j] <- Cluster_col[j]
}
c_simul_col_sorted <- c_simul_col[sort_ind]


postscript("CoClusteringProbs.eps")
par(mar = c(5,7.5,5,7.5))
image(pij_sorted, pty="s", col = rev(heat.colors(100)), main = "", cex.main = 3, cex.lab = 2, axes = FALSE)
image.plot(pij_sorted, pty="s", col = rev(heat.colors(100)), legend.only = TRUE, horizontal = TRUE, axis.args = list(cex.axis = 2), legend.width = 1)
dev.off()
pdf("CoClusteringProbs.pdf")
image(pij_sorted, pty="s", col = rev(heat.colors(100)), main = "", cex.main = 3, cex.lab = 2, axes = FALSE)
image.plot(pij_sorted, pty="s", col = rev(heat.colors(100)), legend.only = TRUE, horizontal = TRUE, axis.args = list(cex.axis = 2), legend.width = 1)
dev.off()

postscript("CoClusteringProbs_withtruth.eps")
heatmap.2(pij_sorted, dendrogram = "none", labRow = FALSE, labCol = FALSE, Rowv = FALSE, Colv = FALSE, ColSideColors = c_simul_col_sorted, RowSideColors = c_simul_col_sorted, na.rm = TRUE, key = FALSE, main = "", trace="none", lwid=c(1,5), lhei=c(1,5), col = rev(heat.colors(100)))
dev.off()
pdf("CoClusteringProbs_withtruth.pdf")
heatmap.2(pij_sorted, dendrogram = "none", labRow = FALSE, labCol = FALSE, Rowv = FALSE, Colv = FALSE, ColSideColors = c_simul_col_sorted, RowSideColors = c_simul_col_sorted, na.rm = TRUE, key = FALSE, main = "", trace="none", lwid=c(1,5), lhei=c(1,5), col = rev(heat.colors(100)))
dev.off()





## Posterior distribution of regression coefficients beta
beta_out <- MCMC_output$XZ_List$beta_out

for(h in 1:pY){
  beta_out_h <- array(0, dim = c(g[h],dY[h]*(dY[h] - 1),n_save))
  for(it in 1:n_save){
    beta_out_h[,,it] <- beta_out[[it]][[h]]
  }
  
  beta_df <- data.frame(dim_g = rep(c(1:(g[h]*dY[h]*(dY[h] - 1))), n_save), beta_vec = c(beta_out_h))
  beta_true_df <- data.frame( x = c(1:(g[h]*dY[h]*(dY[h] - 1))), y = c(beta_simul[[h]]) ) 
  
  postscript(paste("beta_", h, "_violin.eps", sep = ""))
  
  print(
    ggplot(beta_df, aes(y = beta_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
      geom_point(data = beta_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
      labs(x = "Coefficient", y = "", title = bquote("Posterior of "~beta^{(.(h))})) +
      theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) +
      theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"))
  )
  dev.off()
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
  
  postscript(paste("gamma_", h, "_violin.eps", sep = ""))
  
  print(
    ggplot(gamma_df, aes(y = gamma_vec, x = factor(dim_g))) + geom_violin(trim = FALSE) +
      geom_point(data = gamma_true_df, mapping = aes(x = x, y = y, colour = "red", size = 5), show.legend = FALSE) +
      labs(x = "Coefficient", y = "", title = bquote("Posterior of "~gamma^{(.(h))})) +
      theme(text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 12.5)) 
  )
  dev.off()
}


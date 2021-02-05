# SUMS
Seemingly-Unrelated Multi-State processes: a Bayesian semiparametric approach

# Authors

>  **Cremaschi Andrea**,  Singapore Institute for Clinical Sciences, A*STAR, Singapore

>  De Iorio Maria, Yale-NUS College, Singapore

>  Argiento Raffaele, Universita' Cattolica del Sacro Cuore, Milano

>  Meaney Michael,  Singapore Institute for Clinical Sciences, A*STAR, Singapore

>  Kee Michelle,  Singapore Institute for Clinical Sciences, A*STAR, Singapore

# Description
We present a joint semi-parametric model for several possibly related multi-state processes (Seemingly Unrelated Multi-State Models, SUMS), assuming a Markov structure for the transition model. The dependence between different processes is captured by specifying a joint random effect distribution on the transition rates of each process. We assume a nonparametric prior for the random effect distribution, to allow for clustering of the individuals, overdispersions and outliers. Moreover, the precision matrix of the random effect distribution is modelled conditionally to a graph which describes the dependence structure between processes, exploiting tools from the Gaussian Graphical model literature. In this framework, it is also possible to include covariates effect. In this repository, an R/C++ package for the implementation of SUMS is provided, with a simulated data example.

# Contents
1) R package **SUMS**
The current package allows for additional modelling framework than the one presented in the above-mentioned manuscript. In particular, it allows the implementation of

2) R file SUMS_Simul_GitHub.R with simulated example. In this file, different modelling choices, prior elicitation and hyperparameter selection are available.

/#################################
/# **SUMS model: simulated example** #
/#################################

In this simulated example, we demonstrate the SUMS model performance and applicability.

Data are sampled using the msm package (Jackson C.H. 2011) from Multi-State processes whose transition rates are (if needed) covariate-dependent.
We employ simulated covariates of both time-homogeneous (categorical and continuous) and time-varying (continuous) types.
The data are then fitted using the proposed SUMS package.

Along the setting of the package, the user can specify three different modelling options for the random effect distirbution:
(i)   Mixture with random number of components (Argiento and De Iorio 2019)__
(ii)  DP-mixture__
(iii) Parametric Mixture with one component

Furthermore, prior elicitation can be selected for the graph inclusion probability eta between:
(i)   Beta prior
(ii)  Size-based prior (Armstrong et al., 2009)
(iii) Fixed value of eta

R commands to produce some summary plots aimed at posterior inference are available.



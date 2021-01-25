# SUMS
Seemingly-Unrelated Multi-State processes: a Bayesian semiparametric approach

# Authors

>  Cremaschi Andrea, AStar, SICS, Singapore\\
>  De Iorio M., Yale-NUS College, Singapore\\
>  Argiento R., Universit\'{a} Cattolica del Sacro Cuore, Milano\\
>  Meaney M., AStar, SICS, Singapore \\
>  Kee M., AStar, SICS, Singapore

# Description
We present a joint semi-parametric model for several possibly related multi-state processes (Seemingly Unrelated Multi-State Models, SUMS), assuming a Markov structure for the transition model. The dependence between different processes is captured by specifying a joint random effect distribution on the transition rates of each process. We assume a nonparametric prior for the random effect distribution, to allow for clustering of the individuals, overdispersions and outliers. Moreover, the precision matrix of the random effect distribution is modelled conditionally to a graph which describes the dependence structure between processes, exploiting tools from the Gaussian Graphical model literature. In this framework, it is also possible to include covariates effect. In this repository, an R/C++ package for th eimplementation of SUMS is provided.

# Contents
1) R package SUMS
2) R file SUMS_Simul_GitHub.R with simulated example

/#################################
/# SUMS model: simulated example #
/#################################

In this simulated example, we present a simulated example to demonstrate the SUMS model performance and applicability

Data are sampled using the msm package (Jackson C.H. 2011) from Multi-State processes whose transition rates are covariate-dependent
We employ simulated covariates of both time-homogeneous and time-varying (continuous) types
The data are then fitted using the proposed SUMS package
R commands to produce some summary plots aimed at posterior inference are available



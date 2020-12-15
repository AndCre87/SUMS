#ifndef RGWISH_H
#define RGWISH_H

#include <RcppArmadillo.h>

arma::mat rgwish_c(double nu, arma::mat Psi, arma::mat G, double threshold);

#endif

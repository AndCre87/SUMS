// Code from Tesch1, slightly changed to accommodate types that I know how to use.
// When compared to expm package performs very very very similarly
// Much better than arma function for high values of lambda (rates in my specific problem)

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//  Reference:
//
//    Cleve Moler, Charles VanLoan,
//    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
//    Twenty-Five Years Later,
//    SIAM Review,
//    Volume 45, Number 1, March 2003, pages 3-49.
//    This code is based on John Burkhard's matrix_exponential.cpp
//

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "expmat_Tesch1.h"

// [[Rcpp::export]]
arma::mat expmat_Tesch1(const arma::mat & a)
{
  int n = a.n_rows;
  double c;
  arma::mat a2, d, e, x;
  int ee;
  int k;
  int p;
  const int q = 6;
  int s;
  
  ee = (int) log2(norm(a, "inf")) + 1;
  s = std::max(0, ee + 1);
  a2 = a * (1.0 / pow(2.0, s));
  x = a2;
  c = 0.5;
  e = arma::mat(n,n,arma::fill::eye) + c * a2;
  d = arma::mat(n,n,arma::fill::eye) + -c * a2;
  p = 1;
  
  for (k = 2; k <= q; k++) {
    c = c * (q - k + 1.0) / (k * (2.0 * q - k + 1.0));
    x = a2 * x;
    e += c * x;
    if (p) {
      d += c * x;
    }
    else {
      d += -c * x;
    }
    p = !p;
  }
  e = solve(d, e);
  for (k = 1; k <= s; k++) {
    e = e * e;
  }
  
  return e;
}

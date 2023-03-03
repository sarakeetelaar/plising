#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double sumSigma(NumericMatrix sigma,
                NumericMatrix x,
                int index,
                int v) {
  int p = x.ncol();
  double sumsig = 0.0;
  for (int q = 0; q < p; q++) {
    if (q != index) {
      sumsig += sigma(index, q) * x(v, q);
    }
  }
  return sumsig;
}
// [[Rcpp::export]]
NumericMatrix derivativeHelp(NumericMatrix x, 
                             NumericVector mu,
                             NumericMatrix sigma
) {
  int p = x.ncol();
  int n = x.nrow();
  int nParam = p * (p + 1) / 2;
  NumericMatrix derivOuter(nParam, nParam);
  for (int v = 0; v < n; v++) {
    NumericVector deriv( nParam );
    for (int ind = 0; ind < p; ind++) {
      double sigmaSum = sumSigma(sigma, x, ind, v);
      deriv[ind] = x(v, ind) - std::exp(mu[ind] + sigmaSum)/(1 + std::exp(mu[ind] + sigmaSum));
    }
    int count = x.ncol();
    for (int ind = 0; ind < (p - 1); ind++) {
      double iSum = sumSigma(sigma, x, ind, v);
      for (int jInd = (ind + 1); jInd < p; jInd++) {
        double jSum = sumSigma(sigma, x, jInd, v);
        deriv[count] = 2 * x(v, jInd) * x(v, ind) - x(v, ind) * std::exp(mu[jInd] + jSum) / (1 + std::exp(mu[jInd] + jSum)) 
          - x(v, jInd) * std::exp(mu[ind] + iSum)/(1 + std::exp(mu[ind] + iSum));
        count++;
      }
    }
    NumericMatrix outer(nParam, nParam);
    for (int row = 0; row < nParam; row++) {
      for (int col = 0; col < nParam; col++) {
        outer(row, col) = deriv[row] * deriv[col];
      }
    }
    derivOuter += outer;
  }
  return derivOuter;
}

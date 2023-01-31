#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp: export]]
void normalizingConstant(NumericMatrix theta, int p) {
  double Z = 0;
  double Esuf[p*(p+1)/2];
  double E_ss[p*(p+1)/2][p*(p+1)/2];
  arma::vec y(p, arma::fill::zeros);
  int sum_y = 0;
  for (int i =0; i < p; i++) {
    sum_y = sum_y + y(i);
  }
  arma::vec suf;
  while (sum_y < p) {
    if (y[0] == 0) {
      y[0] = 1;
    }
    else {
      for (int j = 1; j < p; j++) {
        if (y[j] == 0) {
          y[j] = 1;
          for (int k = 0; k <= (j-1); k++){
            y[k] = 0;
          }
          break;
        }
      }
    }
    int count = 0;
    for (int i =0; i < p; i++) {
      for (int j = i; j < p; j++) {
        if (j == i) {
          suf(count) = y(i)*y(j)
        }
        else {
        suf(count) = 2*y(i)*y(j);
        }
        count++;
      }
    }
    double weight = exp(suf.t()*theta);
    Z = Z + weight;
    E_suf = E_suf + suf*weight;
    E_ss = E_ss + weight* (suf.t() * suf);
  }
  return 0;
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/

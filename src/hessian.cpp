#include <Rcpp.h>
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
double expprob(NumericMatrix x, 
               NumericMatrix sigma,
               NumericVector mu,
               int v_ind,
               int q_ind) {
  double sumsigma = 0;
  int p = x.ncol();
  for (int q = 0; q < p; q++) {
    if (q != q_ind) {
      sumsigma += sigma(q_ind, q) * x(v_ind, q);
    }
  }
  double pr = std::exp(mu[q_ind] + sumsigma) / (1 + std::exp(mu[q_ind] + sumsigma));
  return pr;
}

// [[Rcpp::export]]
NumericMatrix muHessian(NumericMatrix x,
                             NumericMatrix sigma,
                             NumericVector mu) {
  int p = x.ncol();
  int N = x.nrow();
  NumericMatrix muhes(p, p);
  for (int index = 0; index < p; index ++) {
    double diag = 0;
    for (int vind = 0; vind < N; vind ++) {
      double sigsum = 0;
      for (int q = 0; q < p; q ++) {
        if (q != index) {
          sigsum += sigma(index, q) * x(vind, q);
        }
      }
      double thispr = std::exp(mu[index] + sigsum) / (1 + std::exp(mu[index] + sigsum));
      diag += thispr * (1 - thispr);
    }
    muhes(index, index) = diag;
  }
  return muhes;
}

// [[Rcpp::export]]
NumericVector sigmaHessian(NumericMatrix x,
                           NumericMatrix sigma, 
                           NumericVector mu) {
  int p = x.ncol();
  int N = x.nrow();
  int nTheta = p * (p - 1) / 2;
  NumericMatrix sighes(nTheta, nTheta);
  
  int rowNum = -1;
  int colNum = -1;
  int colShift = -1;
  for (int i1 = 0; i1 < (p-1); i1++) {
    for (int j1 = (i1 + 1); j1 < p; j1++) {
      for (int r = i1; r < (p-1); r++) {
        for (int k = (r + 1); k < p; k ++) {
          double el = 0;
          //case I
          if ((i1 == r) & (j1 == k)) {
            for (int v = 0; v < N; v++) {
              double pi = expprob(x, sigma, mu, v, i1);
              double pj = expprob(x, sigma, mu, v, j1);
              double term = - x(v, j1) * pi * (1 - pi) - x(v, i1) * pj * (1 - pj);

              el += term;
            }
            rowNum ++;
            colShift++;
            colNum = colShift;
            
            sighes(rowNum, colNum) = el;
          }
          else if ((i1 == r) & (j1 != k)) {
            for (int v = 0; v < N; v++) {
              double pi = expprob(x, sigma, mu, v, i1);
              double term = -x(v, j1) * x(v, k) * pi * (1 - pi);

              el += term;
            }
            colNum++;
            
            sighes(rowNum, colNum) = el;
          }
          else if ((i1 != k) & (j1 == r)) {
            for (int v = 0; v < N; v ++) {
              double pj = expprob(x, sigma, mu, v, j1);
              double term = -x(v, i1) * x(v, k) * pj * (1 - pj);

              el += term;
            }
            colNum++;
            sighes(rowNum, colNum) = el;
          }
          else if ((i1 != r) & (j1 == k)) {
            for (int v = 0; v < N; v ++) {
              double pj = expprob(x, sigma, mu, v, j1);
              double term = -x(v, i1) * x(v, r) * pj * (1 - pj);

              el += term;
            }
            colNum++;
            sighes(rowNum, colNum) = el;
          }
          else {
            colNum++;
            sighes(rowNum, colNum) = 0;
          }
        }
      }
    }
  }
  // symmetrize matrix
  for (int row = 0; row < (nTheta - 1);  row ++) {
    for (int col = (row + 1); col < nTheta; col++) {
      sighes(col, row) = sighes(row, col);
    }
  }
  return sighes;
}

// [[Rcpp::export]]
NumericMatrix crossHessian(NumericMatrix x,
                           NumericMatrix sigma,
                           NumericVector mu) {
  int p = x.ncol();
  int N = x.nrow();
  int nTheta = p * (p - 1)/ 2;
  int count = 0;
  NumericMatrix crossHes(p, nTheta);
  for (int mu_ind = 0; mu_ind < p; mu_ind ++) {
    for (int sig_1 = 0; sig_1 < (p - 1); sig_1 ++) {
     for (int sig_2 = (sig_1 + 1); sig_2 < p; sig_2 ++) {
       double el = 0;
       if (sig_1 == mu_ind) {
         for (int v = 0; v < N; v ++) {
           double pi = expprob(x, sigma, mu, v, mu_ind);
           double term = - x(v, sig_2) * pi * (1 - pi);
           el += term;
         }
        
       }
       else if (sig_2 == mu_ind) {
         for (int v = 0; v < N; v ++) {
           double pi = expprob(x, sigma, mu, v, mu_ind);
           double term = -x(v, sig_1) * pi * (1 - pi);
           el += term;
         }
       }
       int rowNum = floor(count / nTheta);
       int colNum = count % nTheta;
       count ++; 
       crossHes(rowNum, colNum) = el;
       
     } 
    }
  }
  
  return crossHes;
}

//Gradient functions 

// [[Rcpp::export]]
NumericVector singleGradientPL(NumericMatrix x,
                         NumericMatrix sigma,
                         NumericVector mu,
                         int obs) {
  int p = x.ncol();
  int nparam = p * (p + 1) / 2;
  NumericVector grad( nparam );
  for (int ind = 0; ind < p; ind ++) {
    double pi = expprob(x, sigma, mu, obs, ind);
    grad[ind] = x(obs, ind) - pi; 
  }
  int index = p;
  for (int ind1 = 0; ind1 < (p - 1); ind1 ++) {
    for (int ind2 = (ind1 + 1); ind2 < p; ind2 ++) {
      double pi = expprob(x, sigma, mu, obs, ind1);
      double pj = expprob(x, sigma, mu, obs, ind2);
      grad[index] = 2 * x(obs, ind1) * x(obs, ind2) - x(obs, ind2) * pi - x(obs, ind1) * pj;
      index++;
    }
  }
  return grad;
}

// [[Rcpp::export]]
NumericMatrix outerGradient(NumericMatrix x,
                            NumericMatrix sigma,
                            NumericVector mu) {
  int p = x.ncol();
  int N = x.nrow();
  int nparam = p * (p + 1) /2;
  NumericMatrix outergrad (nparam, nparam);
  for (int obs = 0; obs < N; obs ++) {
    NumericVector gradient = singleGradientPL(x, sigma, mu, obs);
    NumericMatrix outer(nparam, nparam);
    for (int i_row = 0; i_row < nparam; i_row ++) {
      for (int i_col = 0; i_col < nparam; i_col ++) {
        outer(i_row, i_col) = gradient[i_row] * gradient[i_col];
      }
    }
    outergrad += outer;
  }
  return outergrad;
}



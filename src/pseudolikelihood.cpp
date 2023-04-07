#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double expprobPL(NumericMatrix x, 
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

// [[Rcpp::export]]
int mapSigmaIndex(int index1, int index2, int p) {
  int counter = 0;
  for (int c1 = 0; c1 < (p - 1); c1 ++) {
    for (int c2 = (c1 + 1); c2 < p; c2++) {
      if ((c1 == index1 && c2 == index2) || (c1 == index2 && c2 == index1)) {
        return counter;
      }
      counter ++;
    }
  }
  return counter;
}

// [[Rcpp::export]]
NumericMatrix createSigmaHessian(NumericMatrix sig,
                                 NumericVector mu,
                                 NumericMatrix x) {
  
  int p = x.ncol();
  int N = x.nrow();
  int nSigma = p * (p - 1) / 2;
  NumericMatrix sigHessian( nSigma, nSigma );
  for (int i = 0 ; i < (p - 1) ;i ++) {
    for (int j = (i + 1); j < p; j ++ ) {
      int rowNum = mapSigmaIndex(i, j, p);
      
      for (int r = 0; r < (p - 1); r ++) {
        for (int q = (r + 1); q < p; q ++) {
          int colNum = mapSigmaIndex(r, q, p);
          
          if ((i == r) & (j != q)) {
            double thisterm = 0;
            for (int v = 0; v < N; v++) {
              double pi = expprobPL(x, sig, mu, v, i);
              thisterm -= x(v, j) * x(v, q) * pi * (1 - pi);
            }
            sigHessian(rowNum, colNum) = thisterm;
          }
          else if ((i == q) & (j != r)) {
            double thisterm = 0;
            for (int v = 0; v < N ; v ++) {
              double pi = expprobPL(x, sig, mu, v, i);
              thisterm -= x(v, j) * x(v, q) * pi * (1 - pi);
            }
            sigHessian(rowNum,colNum) = thisterm;
          }
          
          else if ((j == r) & (i != q)) {
            double thisterm = 0;
            for (int v = 0; v < N; v++) {
              double pj = expprobPL(x, sig, mu, v, j);
              thisterm -= x(v, i) * x(v, q) * pj * (1 - pj);
            }
            sigHessian(rowNum, colNum) = thisterm;
            
          }
          else if ((j == q) & (i != r)) {
            double thisterm = 0;
            for (int v = 0; v < N ; v ++) {
              double pj = expprobPL(x, sig, mu, v, j);
              thisterm -= x(v, i) * x(v, r) * pj * (1 - pj);
            }
            sigHessian(rowNum, colNum) = thisterm;
          }
          
          else if ((i == r) & (j == q)) {
            double thisterm = 0;
            for (int v = 0; v < N; v ++) {
              double pi = expprobPL(x, sig, mu, v, j);
              double pj = expprobPL(x, sig, mu, v, j);
              thisterm -= x(v, i) * pj * (1 - pj) + x(v, j) * pi * (1 - pi);
            }
            sigHessian(rowNum, colNum) = thisterm;
            
          }
        }
      }
    }
  }
  return sigHessian;
}


// [[Rcpp::export]]
NumericMatrix createMuHessian(NumericMatrix sig,
                              NumericVector mu,
                              NumericMatrix x) {
  int p = x.ncol();
  int N = x.nrow();
  
  NumericMatrix muhes(p, p);
  
  for (int ind = 0; ind < p; ind++) {
    double thisterm = 0;
    for (int v = 0; v < N; v ++) {
      double pi = expprobPL(x, sig, mu, v, ind);
      thisterm -= pi * (1 - pi);
    }
    muhes(ind, ind) = thisterm;
  }
  
  return muhes;
  
}


// [[Rcpp::export]]
NumericMatrix createCrossHessian(NumericMatrix sig,
                                 NumericVector mu,
                                 NumericMatrix x) {
  int p = x.ncol();
  int N = x.nrow();
  int nSigma = p * (p - 1) / 2;
  
  NumericMatrix crosshes(p, nSigma);
  
  for (int i = 0; i < p; i++) {
    int rowNum = i;
    for (int j = 0; j < p; j++) {
      if (j == i)
        continue;
      int colNum = mapSigmaIndex(i, j, p);
      double thisterm = 0;
      for (int v = 0; v < N; v++) {
        double pi = expprobPL(x, sig, mu, v, i);
        thisterm -= x(v, j) * pi * (1 - pi);
      }
      
      crosshes(rowNum, colNum) = thisterm;
    }
  }
  return crosshes;
}
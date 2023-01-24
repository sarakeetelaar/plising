# include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double normalizingConstant(
  NumericVector theta,
  NumericVector Esuf,
  NumericMatrix Ess,
  double Z,
  NumericVector y
) {
  int p = y.length();
  int nParam = p*(p+1)/2;
  int sumY = 0;
  int count = 1;
  do {
    if (y[0] == 0) {
      y[0] = 1;
    } else {
      for (int q = 1; q < p; q++) {
        if (y[q] == 0) {
          y[q] = 1;
          for (int w = 0; w < q; w++) {
            y[w] = 0;
          }
          break;
        }
      }
    }
    NumericVector suf( nParam );
    int counter = 0;
    for (int i = 0; i < p; i++) {
      for (int j = i; j < p; j++) {
        if (i == j) {
          suf[counter] = y[i]*y[j];
        }
        else {
          suf[counter] = 2*y[i]*y[j];
        }
        counter++;
      }
    }
    double thetaY = 0;
    for (int i =0; i < nParam; i++) {
      double term = suf[i] * theta[i];
      thetaY += term;
    }
    double weight = std::exp(thetaY);
    Z += weight;
    Esuf += weight*suf;
    
    NumericMatrix suftsuf(nParam, nParam);
    for (int i = 0; i < nParam; i++) {
      for (int j = 0; j < nParam; j++) {
        suftsuf(i,j) = suf[i] * suf[j];
      }
    }
    Ess += weight * suftsuf;
    sumY = 0;
    for (int i = 0; i < p; i++) {
      sumY += y[i];
    }
    count++;
  } while (
    !(sumY == p)
  );
  return Z;
}

// [[Rcpp::export]]
double hessenNorm(
  NumericVector thetaH,
  NumericVector EsufH,
  NumericMatrix EssH,
  double ZH,
  IntegerMatrix yH
) {
  int p = yH.ncol();
  int nParam = p*(p+1)/2;
  int nRow = yH.nrow();
  
  for (int row = 0; row < nRow; row ++) {
    IntegerVector yHr = yH(row, _ );
    NumericVector suf( nParam );
    int counter = 0;
    for (int i = 0; i < p; i++) {
      for (int j = i; j < p; j++) {
        if (i == j) {
          suf[counter] = yHr[i]*yHr[j];
        }
        else {
          suf[counter] = 2*yHr[i]*yHr[j];
        }
        counter++;
      }
    }
    double thetaY = 0;
    for (int i =0; i < nParam; i++) {
      double term = suf[i] * thetaH[i];
      thetaY += term;
    }
    double weight = std::exp(thetaY);
    ZH += weight;
    EsufH += weight*suf;
    
    NumericMatrix suftsuf(nParam, nParam);
    for (int i = 0; i < nParam; i++) {
      for (int j = 0; j < nParam; j++) {
        suftsuf(i,j) = suf[i] * suf[j];
      }
    }
    EssH += weight * suftsuf;
    
  }
  return ZH;
}


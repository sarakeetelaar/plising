# include <Rcpp.h>
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



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// [[Rcpp::export]]
void normalizingConstant(
  NumericVector theta,
  NumericVector Esuf,
  NumericMatrix Ess,
  double Z,
  IntegerVector y
) {
  int p = y.length();
  int nParam = p*(p+1)/2;
  int sumY = 0;
  while (sumY != p) {
    
    if (y[0] == 0) {
      y[0] = 1;
    } else{
      for (int i = 1; i < p; i++) {
        if (y[i] == 0) {
          y[i] = 1;
          
          for(int j = 0; j < (i-1); j++) {
            y[j] = 0;
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
    std::cout << Z;
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
  }
}


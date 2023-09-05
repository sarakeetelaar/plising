#include <Rcpp.h>
using namespace Rcpp;
// File contains general functions for the analysis


// [[Rcpp::export]]
NumericMatrix symmetrizeLR(NumericMatrix theta) {
  // transforms output of the logistic regressions 
  // to a symmetric matrix
  int p = theta.ncol();
  for (int i = 0; i < (p - 1); i ++) {
    for (int j = (i + 1); j < p; j ++) {
      theta(i, j) = theta(i, j) + theta(j, i) / 2;
      theta(j, i) = theta(i, j);
    }
  }
  return theta;
}

// [[Rcpp::export]]
NumericVector bootstrapVariances(IntegerMatrix data,
                                 Function func,
                                 Function sampler,
                                 Function removemissing) {
  // calculates bootstrap variances for given measure func
  // sampler gives the sampling method for bootstrap
  // removemissing makes sure that the errors are not included
  int p = data.ncol();
  int n = data.nrow();
  int numParam = p * (p + 1) / 2;
  NumericMatrix allresults(numParam, 1000);
  IntegerMatrix bootSample(n, p);
  for (int rep = 0; rep < 1000; rep++) {
    IntegerVector bootIndices = as<IntegerVector>(sampler(n));
    for (int row = 0; row < n; row ++ ) {
      int thisIndex = bootIndices[row] - 1;
      bootSample(row, _) = data(thisIndex, _);
    }
    NumericVector param = as<NumericVector>(func(bootSample));
    allresults(_, rep) = param;
  }
  allresults = as<NumericMatrix>(removemissing(allresults));
  int numRep = allresults.ncol();
  NumericVector means( numParam );
  for (int rep = 0; rep < numRep; rep++) {
    for (int i = 0; i < numParam; i ++) {
      means[i] += allresults(i, rep) / numRep;
    }
  }

  NumericVector variances( numParam );
  for (int i = 0; i < numParam; i++) {
    for (int rep = 0; rep < numRep; rep++) {
      double term = allresults(i, rep) - means[i];
      variances[i] += term * term;
      
    }
  }
  for (int i = 0; i < numParam; i++) {
    variances[i] /= numRep;
  }
  
  return variances;
}

// [[Rcpp::export]]
NumericMatrix symmetrizeMatrix(NumericMatrix theta) {
  int p = theta.ncol();
  for (int ind1 = 0; ind1 < (p-1); ind1++) {
    for (int ind2 = (ind1+1); ind2 < p; ind2++) {
      theta(ind1, ind2) = (theta(ind1, ind2) + theta(ind2, ind1))/2;
      theta(ind2, ind1) = theta(ind1, ind2);
    }
  }
  return theta;
}


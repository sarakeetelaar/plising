###############################################################################
# This file contains all functions related to the method of independent 
# logistic regressions
###############################################################################

# logistic regression for a single variable, regressed on the others
log_reg <- function(data, index) {
  
  logreg <- glm(formula = data[, index] ~.,
               data = data[, -index], 
               family = binomial(link = "logit"))
  coefs <- logreg$coefficients
  coefs[-1] <- coefs[-1] / 2
  return(coefs)
}

# logistic regression function used for the bootstrap variances
# performs logistic regression given input data frame
lr <- function(data) {
  tryCatch({
    p <- ncol(data)
    coefs <- matrix(0, p, p)
    data <- as.data.frame(data)
    for (i in 1:p) {
      logreg <- glm(formula = data[, i] ~., data = data[, -i], 
                   family = binomial(link = "logit"))
      params <- logreg$coefficients
      params[-1] <- params[-1] / 2
      coefs[i, i] <- params[i]
      coefs[i, -i] <- params[-i]
    }
    coefs <- symmetrizeLR(coefs)
    return(coefs[lower.tri(coefs, diag = TRUE)])
  }, error = function(cond) {
    
    # if it fails return a vector containing -1 
    # these columns will be removed from the bootstrap results
    message(cond)
    p <- ncol(data)
    l <- p * (p + 1) / 2
    return(rep(-1, l))
  })
}

# computes bootstrapped variances for the logistic regression method
bootstrap_variances_lr <- function(x) {
  
  bootvars <- bootstrapVariances(x, lr, sampler, remove_missing_rows)
  muind <- re_index(ncol(x))
  return(list(mu = bootvars[muind], sigma = bootvars[-muind]))
}

# Main function, that performs logistic regression on the data 
# and calculates the variances
logistic_regression <- function(data) {
  
  out = tryCatch({
    p <- ncol(data)
    N <- nrow(data)
    sig <- matrix(0, p, p)
    m <- matrix(0, p, 1)
    df <- data.frame(data)
    
    for (i in 1:p) {
      params <- log_reg(df, i) 
      sig[i, -i] <- params[-1]
      m[i] <- params[1]
    }
    
    variances <- bootstrap_variances_lr(data)
    
    varmu <- variances$mu
    varsig <- variances$sig
    sig <- symmetrizeMatrix(sig)
    
    mu <- list(est = m, var = varmu)
    sigma <- list(est=sig, var = varsig)
    return(list(mu = mu, sigma = sigma))
  },
  error = function(cond){
    # if it fails show that this failed and also 
    # show the reason of failure
    message("logistic regression failed")
    message(cond)
    return(NA)
  })
}


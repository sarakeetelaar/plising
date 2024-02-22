##############################################################################
# File contains all functions related to the pseudolikelihood estimation
##############################################################################

# helper function calculates the log likelihood given data x 
# and parameter matrix sigma with thresholds on diagonal
log_pseudolikelihood <- function(x, sigma) {
  
  p <- ncol(x)
  log_ps <- 0
  for (i in 1:p) {
    f <- sigma[i, i] + x[, -i] %*% sigma[-i, i]
    
    log_ps <- log_ps + x[, i] %*% f
    log_ps <- log_ps - sum(log(1 + exp(f)))
  }
  return(log_ps)
}

# updates all pseudolikelihood values
# one at a time, used for starting values
onedimensional_update = function(x, sigma, suff_stat) {
  p = ncol(x)
  
  # update main effects 
  for(s in 1:p) {
    tmp <- sigma[s, s] + x[, -s] %*% sigma[s, -s]
    prob <- expit(tmp)
    first_derivative <- suff_stat[s, s] - sum(prob)
    second_derivative <- -sum(prob * (1 - prob)) 
    sigma[s, s] <- sigma[s, s] - first_derivative / second_derivative
  }
  
  # update interaction effects 
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      tmp <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
      prob_s <- expit(tmp)
      
      tmp <- sigma[t, t] + x[, -t] %*% sigma[-t, t]
      prob_t <- expit(tmp)
      
      first_derivative <- 2 * suff_stat[s, t] -  x[, t] %*% prob_s
      first_derivative <- first_derivative - x[, s] %*% prob_t
      
      tmp <- prob_s * (1 - prob_s)
      tmp2 <- prob_t * (1 - prob_t)
      second_dervative <- -x[, t] %*% tmp - x[, s] %*% tmp2
      
      sigma[s, t] <- sigma[s, t] - first_derivative / second_dervative
      sigma[t, s] <- sigma[s, t]
    }
  }
  return(sigma)
}

# Performs the newton raphson optimiation for given data set 
optimize_pseudolikelihood <- function(x, iteration_max = 100) {
  p <- ncol(x)
  n <- nrow(x)
  # parameter matrix and sufficient statistics 
  sigma <- matrix(0, nrow = p, ncol = p)
  suff_stat <- t(x) %*% x
  
  # compute starting values 
  for(iteration in 1:5) {
    sigma <- onedimensional_update(sigma = sigma,
                                   x = x,
                                   suff_stat = suff_stat)
  }
  
  log_pseudolikelihood <- log_pseudolikelihood(sigma = sigma,
                                               x = x)
  # index allows us to convert sigma (matrix) to eta (vector) and revert back -
  index <- indexing(p)
  
  log_pseudolikelihood_storage <- c()
  
  for(iteration in 1:iteration_max) {
    # update pseudolikelihood parameters 
    updated <- multidimensional_update(sigma = sigma,
                                       index = index,
                                       x = x,
                                       suff_stat = suff_stat)
    sigma <- updated$sigma

    # update log of pseudolikelihood
    log_pseudolikelihood_old <- log_pseudolikelihood
    log_pseudolikelihood <- log_pseudolikelihood(sigma = sigma,
                                                 x = x)
    # store values
    log_pseudolikelihood_storage[iteration] <- log_pseudolikelihood 
    
    difference <- abs(log_pseudolikelihood - log_pseudolikelihood_old)
    
    if(difference < sqrt(.Machine$double.eps)) break
    
    if(iteration == iteration_max) {
      warning(paste("The optimization procedure did not convergence in", iteration_max, "iterations.",
                    sep = " "), call. = FALSE)
      stop("No convergence")
    }
  }
  
  mu <- diag(sigma)
  diag(sigma) <- 0
  sigma <- sigma/2
  
  var_pl <- get_all_variances(x, sigma, mu, index)


  raw_var <- var_pl$raw
  sandwich_var <- var_pl$sandwich

  if ((length(raw_var[raw_var < 0]) > 0) | (length(sandwich_var[sandwich_var < 0]) > 0)) {
    stop("negative variances")
  }
  return(list(mu = mu, sigma = sigma, var = raw_var, sandwich = sandwich_var))
}

# performs an iteration of newton raphson and updates all parameters together
multidimensional_update <- function(x, sigma, index, suff_stat) {
  n <- nrow(x)
  p <- ncol(x)
  
  # compute vector of first-order derivatives (main effects) 
  prob <- matrix(0, nrow = n, ncol = p)
  derivatives <- vector(length = p * (p + 1) / 2)
  for(s in 1:p) {
    phi <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
    prob[, s] <- expit(phi)
    derivatives[s] <- suff_stat[s, s] - sum(prob[, s])
    derivatives[s] <- derivatives[s] 
  }
  
  # compute vector of first-order derivatives (interaction effects) 
  for(s in 1:(p - 1)) {
    for(t in (s + 1): p) {
      row <- p + index[index[, 2] == s & index[, 3] == t, 1]
      derivatives[row] <- 2 * suff_stat[s, t] - x[, t] %*% prob[, s]
      derivatives[row] <- derivatives[row] - x[, s] %*% prob[, t]
      derivatives[row] <- derivatives[row] 
    }
  }
  # compute the inverse Hessian
  inv_hessian <- invert_hessian(sigma = sigma,
                               index = index,
                               x = x)
  hessian <- inv_hessian$hessian
  inv_hessian <- inv_hessian$inv_hessian


  # convert to eta values 
  eta <- vector(length = p * (p + 1) / 2)
  eta[1:p] <- diag(sigma)
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row <- index[index[, 2] == s & index[, 3] == t, 1] + p
      eta[row] <- sigma[s, t]
    }
  }
  
  # Newton - Raphson update 
  eta <- eta - inv_hessian %*% derivatives
  
  # revert to sigma values 
  diag(sigma) <- eta[1:p]
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row <- p + index[index[,2] == s & index[, 3] == t, 1]
      sigma[s, t] <- eta[row]
      sigma[t, s] <- eta[row]
    }
  }
  return(list(sigma = sigma))
}

# Calculates hessian and inverse hessian
# Hessian and inverse hessian are built up in blocks
# For inverting the schur complement is used
invert_hessian = function(x, sigma, index) {
  # the Hessian matrix is build up as
  #           (main_hessian    , cross_hessian)
  #           (cross_hessian^T , int_hessian)
  # main_hessian contains partial derivatives w.r.t the main effects.
  # int_hessian contains partial derivatives w.r.t the interaction effects.
  tryCatch({
    p <- ncol(x)
    n <- nrow(x)
    # precomputation 
    pq <- matrix(0, nrow = n, ncol = p)
    for(s in 1:p) {
      phi <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
      pq[, s] <- nexpit(phi)
    }
    
    # compute int_hessian
    main_hessian <- matrix(0, nrow = p, ncol = p)
    for (s in 1:p) {
      main_hessian[s,s] <- -sum(pq[, s] * (1 - pq[, s]))
    }
    
    int_hessian <- matrix(0, nrow = p * (p - 1) / 2, ncol = p * (p - 1) / 2)
    for(row in 1:(p * (p - 1) / 2)) {
      s <- index[row, 2]
      t <- index[row, 3]
      
      # partial derivative w.r.t. sigma[s, t]; s < t 
      int_hessian[row, row] <- -x[, t] %*% pq[, s] - x[, s] %*% pq[, t]
      int_hessian[row, row] <- int_hessian[row, row]
      
      # partial derivative w.r.t. sigma[s, t] & sigma [s, r]; s < t 
      I <- which(index[, 2] == s & index[, 3] != t)
      for(i in I) {
        r <- index[i, 3]
        column <- index[i, 1]
        int_hessian[row, column] <- -(x[, t] * x[, r]) %*% pq[, s]
      }
      
      # partial derivative w.r.t. sigma[s, t] & sigma [r, s]; s < t 
      I <- which(index[, 3] == s)
      for(i in I) {
        r <- index[i, 2]
        column <- index[i, 1]
        int_hessian[row, column] <- -(x[, t] * x[, r]) %*% pq[, s]
      }
      
      # partial derivative w.r.t. sigma[s, t] & sigma [k, t]; s < t 
      I <- which(index[, 3] == t & index[, 2] != s)
      for(i in I) {
        k <- index[i, 2]
        column <- index[i, 1]
        int_hessian[row, column] <- -(x[, s] * x[, k]) %*% pq[, t]
      }
      
      # partial derivative w.r.t. sigma[s, t] & sigma [t, k]; s < t
      I <- which(index[, 2] == t)
      for(i in I) {
        k <- index[i, 3]
        column <- index[i, 1]
        int_hessian[row, column] <- -(x[, s] * x[, k]) %*% pq[, t]
      }
    }
    
    # compute cross_hessian 
    cross_hessian <- matrix(0, nrow = p, ncol = p * (p - 1) / 2)
    for(s in 1:p) {
      I <- which(index[, 2] == s)
      for(i in I) {
        t <- index[i, 3]
        column <- index[i, 1]
        cross_hessian[s, column] <- - x[, t] %*% pq[, s]
      }
      
      I <- which(index[, 3] == s)
      for(i in I) {
        t <- index[i, 2]
        column <- index[i, 1]
        cross_hessian[s, column] <- -x[, t] %*% pq[, s]
      }
    }
    
    
    # compute inv_hessian 
    inv_hessian <- matrix(data = 0,
                         nrow = p * (p + 1) / 2,
                         ncol = p * (p + 1) / 2)
    
    # indices for main effects and interactions
    index_main <- 1:p
    index_int <- (1:(p * (p + 1) / 2))[-index_main]
    
    hes <- matrix(nrow = p * (p + 1) / 2, ncol = p * (p + 1) / 2)
    hes[index_main, index_main] <- main_hessian
    hes[index_main, index_int] <- cross_hessian
    hes[index_int, index_main] <- t(cross_hessian)
    hes[index_int, index_int] <- int_hessian
    # inverse of main for Schur complement 
    inv_main <- matrix(0, nrow = p, ncol = p)
    for(s in 1:p) {
      inv_main[s, s] <- 1 / (-sum(pq[, s]))
    }
  
    # use Schur complement for inverting 
    inv_hessian[index_int, index_int] <- int_hessian -
      t(cross_hessian) %*% inv_main %*% cross_hessian
    
    inv_hessian[index_int, index_int] <- solve(inv_hessian[index_int, index_int])
    
    inv_hessian[index_main, index_int] <-
      -inv_main %*% cross_hessian %*% inv_hessian[index_int, index_int]
    
    inv_hessian[index_int, index_main] <- t(inv_hessian[index_main, index_int])
    inv_hessian[index_main, index_main] <- inv_main -
      inv_hessian[index_main, index_int] %*% t(cross_hessian) %*% inv_main
    
    return(list(inv_hessian = inv_hessian, hessian = hes))
  }, 
  error = function(cond) {
    message(cond)
    return(list(inv_hessian = NA, hessian = NA))
  })

}

# Calculates variances from hessian and also the sandwich corrected variances
get_all_variances <- function(x, sigma, mu, index) {
  
  hes <- invert_hessian(x, sigma, index)
  inv_hes <- hes$inv_hessian
  min_inv_hes <- -inv_hes
  raw_variance <- diag(min_inv_hes)
  sandwich_variance <- diag(min_inv_hes %*% outerGradient(x, sigma, mu) %*% min_inv_hes)
  return(list(raw = raw_variance, sandwich = sandwich_variance))
  
}

#helper for bootstrap calculates pseudolikelihood for a bootstrap iteration
pl_estimation  <- function(x) {
  tryCatch({
    res <- optimize_pseudolikelihood(x)
    theta <- res$sigma
    diag(theta) <- res$mu
    return(theta[lower.tri(theta, diag=TRUE)])
  }, error = function(cond) {
    # if it fails return vector of -1
    # this will be removed from the bootstrap estimates
    p = ncol(x)
    l = p * (p + 1) / 2
    return(rep(-1, l))
  })
}

# Main function of the pseudolikelihood returns estimators and 
# three different variance measures
pseudolikelihood <- function(data) {
  out = tryCatch({
    p <- ncol(data)
    # estimation 
    estimators <- optimize_pseudolikelihood(data)
    mu_ind <- re_index(p)
    
    # extract variance measures
    var_raw <- estimators$var
    sandwich_var <- estimators$sandwich
    
    var_mu_r <- var_raw[1:p]
    var_mu_s <- sandwich_var[1:p]
    var_sig_r <- var_raw[-(1:p)]
    var_sig_s <- sandwich_var[-(1:p)]

    var_sig_sand <- matrix(0, p, p)
    var_sig_sand[lower.tri(var_sig_sand)] <- var_sig_s
    var_sig_sand <- var_sig_sand + t(var_sig_sand)
    
    var_sig_raw <- matrix(0, p, p)
    var_sig_raw[lower.tri(var_sig_raw)] <- var_sig_r
    var_sig_raw <- var_sig_raw + t(var_sig_raw)
    
    varmu <- list(raw = var_mu_r, sandwich = var_mu_s)
    varsig <- list(raw = var_sig_raw, sandwich = var_sig_sand)
    
    return(list(mu = list(est = estimators$mu, var = varmu),
                sigma = list(est = estimators$sigma, var = varsig)))
  },
  error = function(cond) {
    # show error message when fails and catch error
    message("pseudolikelihood failed")
    message(cond)
    return(NA)
  })
}

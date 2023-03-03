# inputs of the function x: data set sigma parameters
log_pseudolikelihood = function(x, sigma) {
  p = ncol(x)
  
  log_ps = 0
  
  for (i in 1:p) {
    f = sigma[i,i] + x[, -i] %*% sigma[-i, i]
    
    log_ps = log_ps + x[, i] %*% f
    log_ps = log_ps - sum(log(1 + exp(f)))
  }
  return(log_ps)
}

# updates all pseudolikelihood values
# one at a time to provide reasonable starting values for the function number 6
onedimensional_update = function(x, sigma, suff_stat, prior_var = Inf) {
  p = ncol(x)
  
  # update main effects 
  for(s in 1:p) {
    tmp = sigma[s, s] + x[, -s] %*% sigma[s, -s]
    prob = expit(tmp)
    
    first_derivative = suff_stat[s, s] - sum(prob)
    first_derivative = first_derivative - sigma[s, s] / prior_var
    
    second_derivative = -sum(prob * (1 - prob)) - 1 / prior_var
    
    sigma[s, s] = sigma[s, s] - first_derivative / second_derivative
  }
  
  # update interaction effects 
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      tmp = sigma[s, s] + x[, -s] %*% sigma[-s, s]
      prob_s = expit(tmp)
      
      tmp = sigma[t, t] + x[, -t] %*% sigma[-t, t]
      prob_t = expit(tmp)
      
      first_derivative = 2 * suff_stat[s, t] -  x[, t] %*% prob_s
      first_derivative = first_derivative - x[, s] %*% prob_t
      first_derivative = first_derivative - sigma[s, t] / prior_var
      
      tmp = prob_s * (1 - prob_s)
      tmp2 = prob_t * (1 - prob_t)
      second_dervative = -x[, t] %*% tmp - x[, s] %*% tmp2
      second_dervative = second_derivative - 1 / prior_var
      
      sigma[s, t] = sigma[s, t] - first_derivative / second_dervative
      sigma[t, s] = sigma[s, t]
    }
  }
  
  return(sigma)
  
}


# 5 ----------------------------------------------------------------------------

# computes the inverse Hessian of the log
# pseudolikelihood distribution which is used by function 6

invert_hessian = function(x, sigma, index, prior_var = Inf) {
  # the Hessian matrix is build up as
  #           (main_hessian    , cross_hessian)
  #           (cross_hessian^T , int_hessian)
  # main_hessian contains partial derivatives w.r.t the main effects.
  # int_hessian contains partial derivatives w.r.t the interaction effects.
  
  p = ncol(x)
  n = nrow(x)
  # precomputation 
  pq = matrix(0, nrow = n, ncol = p)
  for(s in 1:p) {
    phi = sigma[s, s] + x[, -s] %*% sigma[-s, s]
    pq[, s] = nexpit(phi)
  }
  
  # compute int_hessian
  main_hessian = matrix(0, nrow = p, ncol = p)
  for (s in 1:p) {
    main_hessian[s,s] = -sum(pq[, s]*(1-pq[, s])) - 1
  }
  
  int_hessian = matrix(0, nrow = p * (p - 1) / 2, ncol = p * (p-1) / 2)
  for(row in 1:(p * (p - 1) / 2)) {
    s = index[row, 2]
    t = index[row, 3]
    
    # partial derivative w.r.t. sigma[s, t]; s < t 
    int_hessian[row, row] = -x[, t] %*% pq[, s] - x[, s] %*% pq[, t]
    int_hessian[row, row] = int_hessian[row, row] - 1 / prior_var
    
    # partial derivative w.r.t. sigma[s, t] & sigma [s, r]; s < t 
    I = which(index[, 2] == s & index[, 3] != t)
    for(i in I) {
      r = index[i, 3]
      column = index[i, 1]
      int_hessian[row, column] = -(x[, t] * x[, r]) %*% pq[, s]
    }
    
    # partial derivative w.r.t. sigma[s, t] & sigma [r, s]; s < t 
    I = which(index[, 3] == s)
    for(i in I) {
      r = index[i, 2]
      column = index[i, 1]
      int_hessian[row, column] = -(x[, t] * x[, r]) %*% pq[, s]
    }
    
    # partial derivative w.r.t. sigma[s, t] & sigma [k, t]; s < t 
    I = which(index[, 3] == t & index[, 2] != s)
    for(i in I) {
      k = index[i, 2]
      column = index[i, 1]
      int_hessian[row, column] = -(x[, s] * x[, k]) %*% pq[, t]
    }
    
    # partial derivative w.r.t. sigma[s, t] & sigma [t, k]; s < t
    I = which(index[, 2] == t)
    for(i in I) {
      k = index[i, 3]
      column = index[i, 1]
      int_hessian[row, column] = -(x[, s] * x[, k]) %*% pq[, t]
    }
  }
  
  # compute cross_hessian 
  cross_hessian = matrix(0, nrow = p, ncol = p * (p - 1) / 2)
  for(s in 1:p) {
    I = which(index[, 2] == s)
    for(i in I) {
      t = index[i, 3]
      column = index[i, 1]
      cross_hessian[s, column] = - x[, t] %*% pq[, s]
    }
    
    I = which(index[, 3] == s)
    for(i in I) {
      t = index[i, 2]
      column = index[i, 1]
      cross_hessian[s, column] = - x[, t] %*% pq[, s]
    }
  }
  
  
  # compute inv_hessian 
  inv_hessian = matrix(data = 0,
                        nrow = p * (p + 1) / 2,
                        ncol = p * (p + 1) / 2)
  
  # indices for main effects and interactions
  index_main = 1:p
  index_int = (1:(p * (p + 1) / 2))[-index_main]
  
  hes = matrix(nrow = p * (p + 1) / 2, ncol = p * (p + 1) / 2)
  hes[index_main, index_main] = main_hessian
  hes[index_main, index_int] = cross_hessian
  hes[index_int, index_main] = t(cross_hessian)
  hes[index_int, index_int] = int_hessian
  # inverse of main for Schur complement 
  inv_main = matrix(0, nrow = p, ncol = p)
  for(s in 1:p)
    inv_main[s, s] = 1 / (-sum(pq[, s]) - 1 / prior_var)
  
  # use Schur complement for inverting 
  inv_hessian[index_int, index_int] = int_hessian -
    t(cross_hessian) %*% inv_main %*% cross_hessian
  inv_hessian[index_int, index_int] = solve(inv_hessian[index_int, index_int])
  inv_hessian[index_main, index_int] =
    -inv_main %*% cross_hessian %*% inv_hessian[index_int, index_int]
  inv_hessian[index_int, index_main] = t(inv_hessian[index_main, index_int])
  inv_hessian[index_main, index_main] = inv_main -
    inv_hessian[index_main, index_int] %*% t(cross_hessian) %*% inv_main
  
  
  return(list(inv_hessian=inv_hessian, hessian=hes))
}

# 6 -----------------------------------------------------------------------------

# Multidimensional Newton Raphson algorithm pseudolikelihood updates used within the main
# optimization function to optimize the pseudologlikelihood parameters
multidimensional_update = function(x, sigma, index, suff_stat, prior_var = Inf) {
  n = nrow(x)
  p = ncol(x)
  
  # compute vector of first-order derivatives (main effects) 
  prob = matrix(0, nrow = n, ncol = p)
  derivatives = vector(length = p * (p + 1) / 2)
  for(s in 1:p) {
    phi = sigma[s, s] + x[, -s] %*% sigma[-s, s]
    prob[, s] = expit(phi)
    derivatives[s] = suff_stat[s, s] - sum(prob[, s])
    derivatives[s] = derivatives[s] - sigma[s, s] / prior_var
  }
  
  # compute vector of first-order derivatives (interaction effects) 
  for(s in 1:(p - 1)) {
    for(t in (s + 1): p) {
      row = p + index[index[, 2] == s & index[, 3] == t, 1]
      derivatives[row] = 2 * suff_stat[s, t] - x[, t] %*% prob[, s]
      derivatives[row] = derivatives[row] - x[, s] %*% prob[, t]
      derivatives[row] = derivatives[row] - sigma[s, t] / prior_var
    }
  }
  
  # compute the inverse Hessian
  inv_hessian = invert_hessian(sigma = sigma,
                                index = index,
                                x = x,
                                prior_var = prior_var)
  hessian = inv_hessian$hessian
  inv_hessian = inv_hessian$inv_hessian
  sig = sigma
  diag(sig) = 0
  mu = diag(sigma)
  derivs = derivativeHelp(x, mu, sig)
  sandwich = (solve(-hessian) %*% derivs %*% solve(-inv_hessian)) / n
  SE = diag(sandwich)
  # convert to eta values 
  eta = vector(length = p * (p + 1) / 2)
  eta[1:p] = diag(sigma)
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row = index[index[,2] == s & index[, 3] == t, 1] + p
      eta[row] = sigma[s, t]
    }
  }
  
  # Newton - Raphson update 
  eta = eta - inv_hessian %*% derivatives
  
  # revert to sigma values 
  diag(sigma) = eta[1:p]
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row = p + index[index[,2] == s & index[, 3] == t, 1]
      sigma[s, t] = eta[row]
      sigma[t, s] = eta[row]
    }
  }
  return(list(sigma=sigma, SE=SE))
}

optimize_pseudolikelihood = function(x, iteration_max = 1e2, prior_var = Inf) {
  out = tryCatch(
    {
      p = ncol(x)
      n = nrow(x)
      # parameter matrix and sufficient statistics 
      sigma = matrix(0, nrow = p, ncol = p)
      suff_stat = t(x) %*% x
      
      # compute starting values 
      for(iteration in 1:5) {
        sigma = onedimensional_update(sigma = sigma,
                                       x = x,
                                       suff_stat = suff_stat,
                                       prior_var = prior_var)
      }
      
      log_pseudolikelihood = log_pseudolikelihood(sigma = sigma,
                                                   x = x)
      
      # index allows us to convert sigma (matrix) to eta (vector) and revert back -
      index = indexing (p)
      
      log_pseudolikelihood_storage = c()
      
      for(iteration in 1:iteration_max) {
        # update pseudolikelihood parameters 
        updated = multidimensional_update(sigma = sigma,
                                           index = index,
                                           x = x,
                                           suff_stat = suff_stat,
                                           prior_var = prior_var)
        sigma = updated$sigma
        
        # update log of pseudolikelihood
        log_pseudolikelihood_old = log_pseudolikelihood
        log_pseudolikelihood = log_pseudolikelihood(sigma = sigma,
                                                     x = x)
        log_pseudolikelihood_storage[iteration] = log_pseudolikelihood #store the pseudolik values for later
        
        difference = abs(log_pseudolikelihood - log_pseudolikelihood_old)
        if(difference < sqrt(.Machine$double.eps)) break
        
        if(iteration == iteration_max)
          warning(paste("The optimization procedure did not convergence in", iteration_max, "iterations.",
                        sep = " "), call. = FALSE)
      }
      se = updated$SE
      mu = diag(sigma)
      diag(sigma) = 0
      sigma = sigma/2
      return(list(se=se, mu = mu, sigma = sigma))
    },
    error = function(cond) {
      message(cond)
      return(NA)
    })
  
}


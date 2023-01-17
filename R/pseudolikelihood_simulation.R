# data generation
library("HelpersMG")
data_generation = function(N=1000,p=10, no.reps=1) {
  ##
  set.seed(2)
  Mu = matrix(runif(p, -1, 1), nrow = p, ncol = 1) 
  Sigma = runif(p, -1, 1) / 2
  Sigma = Sigma %*% t(Sigma)
  diag(Sigma) = 0
  Theta = Sigma
  diag(Theta) = Mu
  Theta = Theta[lower.tri(Theta, diag=TRUE)]
  X_list = list()
  
  for (r in 1: no.reps) {
    X = matrix(0, nrow=N, ncol=p)
    for(iter in 1:1e3) {
      for(i in 1:p) {
        X[, i] = 1 * (rlogis(N) <= Mu[i] + 2 * X %*% Sigma[, i])
      }
    }
    X_list[[r]] = X
  }
  return(list(X=X_list, mu = Mu, sigma=Sigma))
}


#helper function to get indices of mu
re_index = function(p=10) {
  start = 1
  indices = rep(0, p)
  count = 0
  while(count<(p-1)) {
    indices[count+2] = start + (p - count)
    start = start + (p-count)
    count = count+1
    
  }
  indices[1] = 1 
  return(indices)
}
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

# optimize pseudolikelihood
# logistic function 
expit <- function(phi) {
  return(exp(phi) / (1 + exp(phi)))
}


# 2 ----------------------------------------------------------------------------
# negative logistic function

nexpit <- function(phi) {
  return(exp(phi) / (1 + exp(phi)) ^ 2)
}

# 3 ----------------------------------------------------------------------------
# indexing function used to convert the pxp
# parameter matrix sigma to a {p * (p + 1) / 2} vector

indexing <- function(p) {
  index <- matrix(0, nrow = p * (p - 1) / 2, ncol = 3)
  counter <- 0
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      counter <- counter + 1
      index[counter, 1] <- counter
      index[counter, 2] <- s
      index[counter, 3] <- t
    }
  }
  return(index)
}

# 4 ----------------------------------------------------------------------------

# updates all pseudolikelihood values
# one at a time to provide reasonable starting values for the function number 6

onedimensional_update <- function(x, sigma, suff_stat, prior_var = Inf) {
  p <- ncol(x)
  
  # update main effects 
  for(s in 1:p) {
    tmp <- sigma[s, s] + x[, -s] %*% sigma[s, -s]
    prob <- expit(tmp)
    
    first_derivative <- suff_stat[s, s] - sum(prob)
    first_derivative <- first_derivative - sigma[s, s] / prior_var
    
    second_derivative <- -sum(prob * (1 - prob)) - 1 / prior_var
    
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
      first_derivative <- first_derivative - sigma[s, t] / prior_var
      
      tmp <- prob_s * (1 - prob_s)
      tmp2 <- prob_t * (1 - prob_t)
      second_dervative <- -x[, t] %*% tmp - x[, s] %*% tmp2
      second_dervative <- second_derivative - 1 / prior_var
      
      sigma[s, t] <- sigma[s, t] - first_derivative / second_dervative
      sigma[t, s] <- sigma[s, t]
    }
  }
  
  return(sigma)
  
}


# 5 ----------------------------------------------------------------------------

# computes the inverse Hessian of the log
# pseudolikelihood distribution which is used by function 6

invert_hessian <- function(x, sigma, index, prior_var = Inf) {
  # the Hessian matrix is build up as
  #           (main_hessian    , cross_hessian)
  #           (cross_hessian^T , int_hessian)
  # main_hessian contains partial derivatives w.r.t the main effects.
  # int_hessian contains partial derivatives w.r.t the interaction effects.
  
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
    main_hessian[s,s] <- -sum(pq[, s]*(1-pq[,s])) - 1
  }
  
  int_hessian <- matrix(0, nrow = p * (p - 1) / 2, ncol = p * (p-1) / 2)
  for(row in 1:(p * (p - 1) / 2)) {
    s <- index[row, 2]
    t <- index[row, 3]
    
    # partial derivative w.r.t. sigma[s, t]; s < t 
    int_hessian[row, row] <- -x[, t] %*% pq[, s] - x[, s] %*% pq[, t]
    int_hessian[row, row] <- int_hessian[row, row] - 1 / prior_var
    
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
      cross_hessian[s, column] <- - x[, t] %*% pq[, s]
    }
  }
  
  
  # compute inv_hessian 
  inv_hessian <- matrix(data = 0,
                        nrow = p * (p + 1) / 2,
                        ncol = p * (p + 1) / 2)
  
  # indices for main effects and interactions
  index_main <- 1:p
  index_int <- (1:(p * (p + 1) / 2))[-index_main]
  
  hes = matrix(nrow=p*(p+1)/2, ncol=p*(p+1)/2)
  hes[index_main, index_main] = main_hessian
  hes[index_main, index_int] = cross_hessian
  hes[index_int, index_main] = t(cross_hessian)
  hes[index_int, index_int] = int_hessian
  # inverse of main for Schur complement 
  inv_main <- matrix(0, nrow = p, ncol = p)
  for(s in 1:p)
    inv_main[s, s] <- 1 / (-sum(pq[, s]) - 1 / prior_var)
  
  # use Schur complement for inverting 
  inv_hessian[index_int, index_int] <- int_hessian -
    t(cross_hessian) %*% inv_main %*% cross_hessian
  inv_hessian[index_int, index_int] <- solve(inv_hessian[index_int, index_int])
  inv_hessian[index_main, index_int] <-
    -inv_main %*% cross_hessian %*% inv_hessian[index_int, index_int]
  inv_hessian[index_int, index_main] <- t(inv_hessian[index_main, index_int])
  inv_hessian[index_main, index_main] <- inv_main -
    inv_hessian[index_main, index_int] %*% t(cross_hessian) %*% inv_main
  
  
  return(list(inv_hessian=inv_hessian, hessian=hes))
}

# 6 -----------------------------------------------------------------------------

# Multidimensional Newton Raphson algorithm pseudolikelihood updates used within the main
# optimization function to optimize the pseudologlikelihood parameters
multidimensional_update <- function(x, sigma, index, suff_stat, prior_var = Inf) {
  n <- nrow(x)
  p <- ncol(x)
  
  # compute vector of first-order derivatives (main effects) 
  prob <- matrix(0, nrow = n, ncol = p)
  derivatives <- vector(length = p * (p + 1) / 2)
  for(s in 1:p) {
    phi <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
    prob[, s] <- expit(phi)
    derivatives[s] <- suff_stat[s, s] - sum(prob[, s])
    derivatives[s] <- derivatives[s] - sigma[s, s] / prior_var
  }
  
  # compute vector of first-order derivatives (interaction effects) 
  for(s in 1:(p - 1)) {
    for(t in (s + 1): p) {
      row <- p + index[index[, 2] == s & index[, 3] == t, 1]
      derivatives[row] <- 2 * suff_stat[s, t] - x[, t] %*% prob[, s]
      derivatives[row] <- derivatives[row] - x[, s] %*% prob[, t]
      derivatives[row] <- derivatives[row] - sigma[s, t] / prior_var
    }
  }
  
  # compute the inverse Hessian
  inv_hessian <- invert_hessian(sigma = sigma,
                                index = index,
                                x = x,
                                prior_var = prior_var)
  hessian = inv_hessian$hessian
  inv_hessian = inv_hessian$inv_hessian
  
  SE = SEfromHessian(hessian)
  # convert to eta values 
  eta <- vector(length = p * (p + 1) / 2)
  eta[1:p] <- diag(sigma)
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row <- index[index[,2] == s & index[, 3] == t, 1] + p
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
  return(list(sigma=sigma, SE=SE))
}

optimize_pseudolikelihood <- function(x, iteration_max = 1e2, prior_var = Inf) {
  p <- ncol(x)
  n <- nrow(x)
  
  # parameter matrix and sufficient statistics 
  sigma <- matrix(0, nrow = p, ncol = p)
  suff_stat <- t(x) %*% x
  
  # compute starting values 
  for(iteration in 1:5) {
    sigma <- onedimensional_update(sigma = sigma,
                                   x = x,
                                   suff_stat = suff_stat,
                                   prior_var = prior_var)
  }
  
  log_pseudolikelihood <- log_pseudolikelihood(sigma = sigma,
                                               x = x)
  
  # index allows us to convert sigma (matrix) to eta (vector) and revert back -
  index <- indexing (p)
  
  log_pseudolikelihood_storage <- c()
  
  for(iteration in 1:iteration_max) {
    # update pseudolikelihood parameters 
    updated <- multidimensional_update(sigma = sigma,
                                       index = index,
                                       x = x,
                                       suff_stat = suff_stat,
                                       prior_var = prior_var)
    sigma = updated$sigma
    
    # update log of pseudolikelihood
    log_pseudolikelihood_old <- log_pseudolikelihood
    log_pseudolikelihood <- log_pseudolikelihood(sigma = sigma,
                                                 x = x)
    log_pseudolikelihood_storage[iteration] <- log_pseudolikelihood #store the pseudolik values for later
    
    difference <- abs(log_pseudolikelihood - log_pseudolikelihood_old)
    if(difference < sqrt(.Machine$double.eps)) break
    
    if(iteration == iteration_max)
      warning(paste("The optimization procedure did not convergence in", iteration_max, "iterations.",
                    sep = " "), call. = FALSE)
  }
  ll = log_pseudolikelihood
  se <- updated$SE
  mu <- diag(sigma)
  diag(sigma) <- 0
  sigma = sigma/2
  mse_mu = mean((mu - Mu) ^2)
  mse_sigma = mean((sigma[lower.tri(sigma)] - Sigma[lower.tri(Sigma)]) ^2)
  return(list(likelihood=ll,se =se, mu = mu, sigma = sigma, mse_mu = mse_mu, mse_sigma = mse_sigma, log_pseudolikelihood = log_pseudolikelihood_storage))
  
}


#estimators_pl = optimize_pseudolikelihood(X)

pseudolikelihood = function(data) {
  estimators = optimize_pseudolikelihood(data)
  
  sigma = estimators$sigma
  mu = estimators$mu
  mse_m = estimators$mse_mu
  mse_s = estimators$mse_sigma
  se = estimators$se
  p = length(mu)
  se_mu = se[1:p]
  se_sigma = se[-(1:p)]
  ll = estimators$likelihood
  return(list(likelihood=ll, sigma_pl = sigma, mu_pl = mu, mse_sigma = mse_s, mse_mu = mse_m, se_mu = se_mu, se_sigma = se_sigma))
}





exact_likelihood = function(data, theta_0) {
  N = nrow(data)
  O = t(data) %*% data * 2
  diag(O) = diag(O)/2
  O = O[lower.tri(O, diag=TRUE)]
  
  p = ncol(data)
  theta = theta_0[lower.tri(theta_0, diag=TRUE)]
  for(iter in 1:1e3) {
    theta_old = theta
    
    E_suf = matrix(0, nrow = p*(p+1)/2, ncol=1)
    E_ss = matrix(0, nrow = p*(p+1)/2, ncol = p*(p+1)/2)
    
    
    Z = 0
    y = matrix(0, nrow = p, ncol = 1)
    
    repeat {
      if(y[1] == 0) {
        y[1] = 1
      } else {
        for(i in 2:p) {
          if(y[i] == 0) {
            y[i] = 1
            y[1:(i-1)] = 0
            break
          }
        }
      }
      suf = 2 * y %*% t(y)
      diag(suf) = diag(suf)/2
      suf = suf[lower.tri(suf, diag=TRUE)]
      
      weight = as.numeric(exp(t(suf)%*%theta))
      
      Z = Z + weight
      
      E_suf = E_suf + suf * weight
      
      E_ss = E_ss + suf %*% t(suf) * weight
      
      
      if(sum(y) == p)
        break
    }
    
    E_suf = E_suf/Z
    E_ss = E_ss/Z
    
    H = -N * (E_ss - E_suf %*% t(E_suf))
    
    grad = O - N*E_suf
    
    
    theta = theta - solve(H)%*%grad
    
    ll = t(O) %*% theta - N * log(Z)
    D = sum(abs(theta-theta_old))
    if(D < sqrt(.Machine$double.eps)) {
      break
    }
  }
  #SE = sqrt(diag(solve(-H)))
  ll_value = ll
  SE = SEfromHessian(-H)
  indices = re_index(p)
  
  SE_mu = SE[indices]
  SE_sigma = SE[-indices]
  
  theta_l = theta
  
  theta_l = matrix(0, nrow=p, ncol=p)
  theta_l[lower.tri(theta_l,diag=TRUE)] = theta
  
  for (i in 1:(p-1)) {
    for (j in (i+1): p) {
      theta_l[i,j] = theta_l[j,i]
    }
  }
  mu_ml = diag(theta_l)
  diag(theta_l) = 0
  
  
  mse_muML = mean((mu_ml - Mu)^2)
  mse_sigML = mean((theta_l[lower.tri(theta_l)] - Sigma[lower.tri(Sigma)])^2)
  
  return(list(likelihood=ll, mu_ML = mu_ml, sigma_ML = theta_l, mse_mu = mse_muML, mse_sigma = mse_sigML, SE_mu = SE_mu, SE_sigma = SE_sigma))
  
}
#------ Actual performance of both the exact likelihood and the pseudolikelihood

# # 
# set.seed(1)
# N_list = seq(100, 1000, 100)
# l = length(N_list)
# mse_muPL = rep(0, l)
# mse_muML = rep(0, l)
# mse_sigPL = relative_bias_mu = rep(0, l)
# mse_sigML = relative_bias_sig = rep(0, l)
# likelihood_values = pseudolikelihood_values = rep(0,l)
# 
# se_muPL = rep(0,l)
# se_sigPL = rep(0,l)
# se_muML = rep(0, l)
# se_sigML = rep(0,l)
# time_PL = rep(0, l)
# time_ML = rep(0, l)
# no.reps = 100
# 
# for (i in 1:l) {
#   n = N_list[i]
#   print(paste("n = ", n))
#   data_gen = data_generation(N=n, no.reps=no.reps)
#   data = data_gen$X
#   Mu = data_gen$mu
#   Sigma = data_gen$sigma
# 
# 
#   for (rep in 1:no.reps) {
#     this_data = data[[rep]]
#     time1 = Sys.time()
#     est_pl = pseudolikelihood(this_data)
#     time2 = Sys.time()
# 
#     mu_pl = est_pl$mu_pl
#     sig_pl = est_pl$sigma_pl
# 
#     diag(sig_pl) = mu_pl
# 
#     time22 = Sys.time()
#     est_ml = exact_likelihood(this_data, sig_pl)
#     time3 = Sys.time()
# 
#     time_PL[i] = time_PL[i] + (time2 - time1)/no.reps
#     time_ML[i] = time_ML[i] + (time3 - time22)/no.reps
# 
#     mse_muPL[i] = mse_muPL[i] + est_pl$mse_mu/no.reps
#     mse_sigPL[i] = mse_sigPL[i] + est_pl$mse_sigma/no.reps
# 
#     mse_muML[i] = mse_muML[i] + est_ml$mse_mu/no.reps
#     mse_sigML[i] = mse_sigML[i] + est_ml$mse_sigma/no.reps
# 
#     se_muPL[i] = se_muPL[i] + mean(est_pl$se_mu)/no.reps
#     se_sigPL[i] = se_sigPL[i] + mean(est_pl$se_sigma)/no.reps
# 
#     se_muML[i] = se_muML[i] + mean(est_ml$SE_mu)/no.reps
#     se_sigML[i] = se_sigML[i] + mean(est_ml$SE_sigma)/no.reps
#     
#     relative_bias_mu[i] = relative_bias_mu[i] + (est_pl$mse_mu/est_ml$mse_mu)/no.reps
#     relative_bias_sig[i] = relative_bias_sig[i] + (est_pl$mse_sigma/est_ml$mse_sigma)/no.reps
#     
#   }
# 
# }
# 
# 
# 
# parameter_plots_n = function() {
#   Ns = c(100, 1000, 2500, 5000)
#   par(mfrow=c(2,2))
#   no.reps=100
#   for (n in 1:length(Ns)) {
#     N = Ns[n]
#     print(paste("N =", N))
#     set.seed(n)
#     data = data_generation(N=N, no.reps=no.reps)
#     X = data$X
#     Mu = data$mu
#     Sigma = data$sigma
#     p = length(Mu)
#     sigma_pl = sigma_ml = matrix(0, p, p)
#     mu_pl = mu_ml = se_mu_ml = se_mu_pl =rep(0,p)
#     se_sig_ml = se_sig_pl = rep(0, p*(p-1)/2)
#     for (rep in 1:no.reps) {
#       print(paste("rep =", rep))
#       x = X[[rep]]
#       pseu = pseudolikelihood(x)
#       thetapl = pseu$sigma_pl
#       diag(thetapl) = pseu$mu_pl
#       sigma_pl = sigma_pl + pseu$sigma_pl/no.reps
#       mu_pl = mu_pl + pseu$mu_pl/no.reps
#       
#       exact = exact_likelihood(x, thetapl)
#       sigma_ml = sigma_ml + exact$sigma/no.reps
#       mu_ml = mu_ml + exact$mu/no.reps
#       
#       se_mu_pl = se_mu_pl + pseu$se_mu/no.reps
#       se_sig_pl = se_sig_pl + pseu$se_sigma/no.reps
#       
#       se_mu_ml = se_mu_ml + exact$SE_mu/no.reps
#       se_sig_ml = se_sig_ml + exact$SE_sigma/no.reps
#     }
#     plot(mu_ml, mu_pl, main=paste("N =", N), xlab="Mu ML", ylab="Mu PL");abline(0,1)
#     
#   }
# }


# # -------------------- varying number of variables -------------------------

run_comparison_p = function() {
  p_list = 5:14
  no.reps = 100
  len = length(p_list)
  mse_mu_ml = mse_sig_ml = mse_mu_pl = mse_sig_pl = mse_mu_H = mse_sig_H = rep(0, len)
  mse_pl_p = mse_ml_p = times_ml = times_pl = times_hes = rep(0, len)
  
  se_pl_mu = se_pl_sig = se_ml_mu = se_ml_sig = rep(0, len)
  ll_values = pl_values = rep(0, len)
  for (i in 1:len) {
    p = p_list[i]
    data = data_generation(p=p, no.reps=no.reps)
    X = data$X
    Mu = data$mu
    Sigma = data$sigma
    
    for (rep in 1:no.reps) {
      x = X[[rep]]
      time_1 = Sys.time()
      pl = pseudolikelihood(x)
      
      mu_PL = pl$mu_pl
      sigma_PL = pl$sigma_pl
      
      diag(sigma_PL) = mu_PL
      time_2 = Sys.time()
      ml = exact_likelihood(x, sigma_PL)
      time_3 = Sys.time()
      # 
      # mse_mu_ml[i] = mse_mu_ml[i] + ml$mse_mu/no.reps
      # mse_sig_ml[i] = mse_sig_ml[i] + ml$mse_sigma/no.reps
      # mse_mu_pl[i] = mse_mu_pl[i] + pl$mse_mu/no.reps
      # mse_sig_pl[i] = mse_sig_pl[i] + pl$mse_sigma/no.reps
      ll_values[i] = ll_values[i] + ml$likelihood/no.reps
      pl_values[i] = pl_values[i] + pl$likelihood/no.reps
      
      # se_pl_mu[i] = se_pl_mu[i] + mean(pl$se_mu)/no.reps
      # se_pl_sig[i] = se_pl_sig[i] + mean(pl$se_sigma)/no.reps
      # se_ml_mu[i] = se_ml_mu[i] + mean(ml$SE_mu)/no.reps
      # se_ml_sig[i] = se_ml_sig[i] + mean(ml$SE_sigma)/no.reps
      # times_pl[i] = times_pl[i] + (time_2 - time_1)/no.reps
      # times_ml[i] = times_ml[i] + (time_3 - time_2)/no.reps
    }
  }
  par(mfrow=c(2,1))
  fractions = pl_values/ll_values
  ymin = min(c(ll_values, pl_values))
  ymax = max(c(ll_values, pl_values))
  
  plot(p_list, ll_values, ylim=c(ymin, ymax), xlab="p", ylab="likelihood", type="l", col="red")
  lines(p_list, pl_values, type="l", col="blue")
  
  plot(p_list, fractions, type="l", xlab="p", ylab="fraction likelihood")
}

parameter_plots_p = function() {
  Ps = c(5, 7, 10, 15)
  par(mfrow=c(2,2))
  no.reps = 100
  for (p in Ps) {
    print(paste("p =", p))
    data_gen = data_generation(p=p, no.reps=no.reps)
    muex = mupl =  rep(0, p)
    sigex = sigpl = matrix(0, p, p)
    se_mupl = se_muml = rep(0, p)
    se_sigpl = se_sigml = rep(0, p*(p-1)/2)
    X = data_gen$X
    Mu = data_gen$mu
    Sigma = data_gen$sigma
    for (rep in 1:no.reps) {
      if (rep%%10 == 0) {
        print(paste("rep", rep))
      }
      x = X[[rep]]
      pseu = pseudolikelihood(x)
      print("pseudolikelihood finished")
      mupl = mupl + pseu$mu_pl/no.reps
      thistheta = pseu$sigma_pl
      diag(thistheta) = pseu$mu_pl
      sigpl = sigpl + pseu$sigma_pl/no.reps
      diag(sigpl) = mupl
      
      exact = exact_likelihood(x, thistheta)
      muex = muex + exact$mu/no.reps
      sigex = sigex + exact$sigma/no.reps
      
      se_mupl = se_mupl + pseu$se_mu/no.reps
      se_sigpl = se_sigpl + pseu$se_sigma/no.reps
      
      se_muml = se_muml + exact$SE_mu/no.reps
      se_sigml = se_sigml + exact$SE_sigma/no.reps
    }
    plot(se_sigml, se_sigpl, main=paste("p =", p), xlab= "ML Sigma", ylab="PL Sigma");abline(0,1)
    
    
  }
}
p = 10
y = rep(0,p)
theta = rep(0, p*(p+1)/2)
Esuf = matrix(0, nrow=p*(p+1)/2, 1)
Ess = matrix(0, nrow=p*(p+1)/2, ncol=p*(p+1)/2)
Z = 0.0



library(igraph)
random_network = function(p, m) {
  gr = random.graph.game(p, m, type="gnm", loops=F)
  A = as_adjacency_matrix(gr)
  return(as.matrix(A))
}

create_small_world = function(p, nei=1, pr=.1) {
  world = sample_smallworld(1, p, nei, pr)
  A = as_adjacency_matrix(world)
  return(as.matrix(A))
}

scale_free_network = function(p=10, m=NULL) {
  graph = barabasi.game(p, 1, m, directed=FALSE)
  A = as_adjacency_matrix(graph)
  return(as.matrix(A))
}

data_generation = function(N=1000,p=10, no.reps=1, A=NULL, Mu, sigma_values) {
  ##
  set.seed(2)
  Sigma = matrix(0, p,p)
  Sigma[lower.tri(Sigma)] = sigma_values
  
  #apply graph structure on parameters
  if (!is.null(A))
    Sigma[A==0] = 0
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      Sigma[i,j] = Sigma[j,i]
    }
  } 
  X_list = list()
  
  for (r in 1: no.reps) {
    X = matrix(0, nrow=N, ncol=p)
    for(iter in 1:1e2) {
      for(i in 1:p) {
        X[, i] = 1 * (rlogis(N) <= Mu[i] + 2 * X %*% Sigma[, i])
      }
    }
    X_list[[r]] = X
  }
  return(list(X=X_list, mu=Mu, sigma=Sigma))
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
  SE = diag(inv_hessian)
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
  return(list(likelihood=ll,se =se, mu = mu, sigma = sigma))
  
}


#estimators_pl = optimize_pseudolikelihood(X)

pseudolikelihood = function(data) {
  out = tryCatch(
    {
    estimators = optimize_pseudolikelihood(data)
    sigma = estimators$sigma
    mu = estimators$mu
    
    se = estimators$se
    p = length(mu)
    se_mu = se[1:p]
    se_sigma = se[-(1:p)]
    ll = estimators$likelihood
    return(list(sigma_pl = sigma, mu_pl = mu, se_mu = se_mu, se_sigma = se_sigma))
  },
  error = function(cond) {
    message(cond)
    return(NA)
  }
  )
}

hessen = function(data, theta_0=matrix(0, p, p)) {
  out = tryCatch( 
    {
      N = nrow(data)
      O = t(data) %*% data * 2
      diag(O) = diag(O)/2
      O = O[lower.tri(O, diag=TRUE)]
      
      ys = unique(data)
      p = ncol(data)
      theta = theta_0[lower.tri(theta_0, diag=TRUE)]
      for(iter in 1:1e3) {
        theta_old = theta
        
        E_suf = matrix(0, nrow = p*(p+1)/2, ncol=1)
        E_ss = matrix(0, nrow = p*(p+1)/2, ncol = p*(p+1)/2)
        
        
        Z = 0.0
        Z = hessenNorm(theta=theta, EsufH=E_suf, EssH=E_ss, ZH=Z, yH=ys)
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
      SE = diag(solve(-H))
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
      
      return(list(likelihood=ll, mu_ML = mu_ml, sigma_ML = theta_l, SE_mu = SE_mu, SE_sigma = SE_sigma))
    },
    error = function(cond) {
      message(cond)
      return(NA)
    }
  )
  
}




exact_likelihood = function(data, theta_0=matrix(0, p, p)) {
  out = tryCatch({
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
      
      
      Z = 0.0
      y = matrix(0, nrow = p, ncol = 1)
      Z = normalizingConstant(theta, E_suf, E_ss, Z, y)
      
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
    SE = diag(solve(-H))
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
    
    return(list(mu_ML = mu_ml, sigma_ML = theta_l, SE_mu = SE_mu, SE_sigma = SE_sigma))
  },
  error = function(cond) {
    message(cond)
    return(NA)
  }
  )
}

logistic_regression = function(data) {
  p = ncol(data)
  N = nrow(data)
  sigma = matrix(0, p, p)
  mu = matrix(0, p, 1)
  df = data.frame(data)
  for (i in 1:p) {
    logreg = glm(formula=df[,i]~., data=df[,-i], family=binomial(link="logit"))
    params = logreg$coefficients
    sigma[i, -i] = params[-1]
    mu[i] = params[1]
  }
  sigma = symmetrizeMatrix(sigma)

  return(list(mu=mu, sigma=sigma))
  
}

#function performance of the pseudolikelihood
# now only implemented for the scale free graph
perform_comparison = function(p_options, N_options, no.reps, graph_type='full', frac=NULL) {
  sz_p = length(p_options)
  sz_N = length(N_options)
  
  # to be adjusted
  graphs = list()
  
  real_mu= list()
  real_sigma = list()
  
  mu_pl = mu_ml = sigma_pl = sigma_ml = list()
  var_mu_pl = var_sig_pl = var_mu_ml = var_sig_ml = list()
  
  for (p_ind in 1:sz_p) {
    p = p_options[p_ind]
    if (graph_type == 'sf') {
      graphs[[p_ind]] = scale_free_network(p=p)
    }
    else if (graph_type == 'sm') {
      graphs[[p_ind]] = create_small_world(p=p)
    }
    else if (graph_type=='empty') {
      graphs[[p_ind]] = matrix(0, p, p)
    }
    else if (!is.null(frac)) {
      gr = rbinom(p*(p-1)/2, 1, frac)
      g = matrix(0, p, p)
      g[lower.tri(g)] = gr
      for (i in 1:(p-1)) {
        for (j in (i+1):p) {
          g[i,j] = g[j,i]
        }
      }
      graphs[[p_ind]] = g
    }
    else {
      g = matrix(1, p, p)
      diag(g) = 0
      graphs[[p_ind]] = g
    }
    
    mu = runif(p, -1,1)
    sigma_val = matrix(runif(p*(p-1)/2, -1, 1)/2, nrow=p*(p-1)/2, ncol=1)
    # sw = create_small_world(p=p)
    # small_worlds[[p_ind]] = sw
    
    mu_pl[[p_ind]] = mu_ml[[p_ind]] = sigma_pl[[p_ind]] = sigma_ml[[p_ind]] = list()
    var_mu_pl[[p_ind]] = var_sig_pl[[p_ind]] = var_mu_ml[[p_ind]] = var_sig_ml[[p_ind]] = list()
    
    for (n_ind in 1:sz_N) {

      n = N_options[n_ind]
      print(paste("p =", p, "n =", n))
      data = data_generation(N=n, p=p, no.reps=no.reps, Mu=mu, sigma_values=sigma_val, A=graphs[[p_ind]])
      X_list = data$X
      mu = data$mu
      sigma = data$sigma
      
      if (n_ind == 1) {
        real_mu[[p_ind]] = mu
        real_sigma[[p_ind]] = sigma
      }
      
      mu_PL = mu_ML = matrix(0, p, 1)
      sig_PL = sig_ML = matrix(0, p*(p-1)/2, 1)
      se_mu_PL = se_mu_ML = matrix(0, p, 1)
      se_sig_PL = se_sig_ML = matrix(0, p*(p-1)/2, 1)
      
      for (rep in 1:no.reps) {
        x = X_list[[rep]]
        pseudo = pseudolikelihood(x)
        print("pseudolikelihood finished")
        sigpl = pseudo$sigma_pl
        mupl = pseudo$mu_pl
        mu_PL = mu_PL + mupl/no.reps
        sig_PL = sig_PL + sigpl[lower.tri(sigpl)]/no.reps
        se_mu_PL = se_mu_PL + pseudo$se_mu/no.reps
        se_sig_PL = se_sig_PL + pseudo$se_sigma/no.reps
        
        diag(sigpl) = mupl
        likeli = exact_likelihood(x, sigpl)
        print("likelihood finished")
        sigml = likeli$sigma_ML
        muml = likeli$mu_ML
        
        sig_ML = sig_ML + sigml[lower.tri(sigml)]/no.reps
        mu_ML = mu_ML + muml/no.reps
        se_mu_ML = se_mu_ML + likeli$SE_mu/no.reps
        se_sig_ML = se_sig_ML + likeli$SE_sigma/no.reps
        
        
        
      }
      mu_pl[[p_ind]][[n_ind]] = mu_PL
      sigma_pl[[p_ind]][[n_ind]] = sig_PL
      mu_ml[[p_ind]][[n_ind]] = mu_ML
      sigma_ml[[p_ind]][[n_ind]] = sig_ML
      
      var_mu_pl[[p_ind]][[n_ind]] = se_mu_PL
      var_sig_pl[[p_ind]][[n_ind]] = se_sig_PL
      var_mu_ml[[p_ind]][[n_ind]] = se_mu_ML
      var_sig_ml[[p_ind]][[n_ind]] = se_sig_ML
    }
    
  }
  return(list(sigma_true=real_sigma, mu_true=real_mu, mu_pl=mu_pl, sigma_pl=sigma_pl,
              mu_ml=mu_ml, sigma_ml=sigma_ml, var_mu_pl=var_mu_pl, var_sigma_pl=var_sig_pl,
              var_mu_ml=var_mu_ml, var_sigma_ml=var_sig_ml, graphs=graphs))
}


pseudolikelihood = function(data) {
  out = tryCatch({
      
      estimators = optimize_pseudolikelihood(data)
      var_raw = estimators$var
      sandwich_var = estimators$sandwich
      #var_boot = bootstrap_var_pl(data, B=1e2)
      
      p = ncol(data)
      

      # mu = list(est=estimators$mu, var=var_mu)
      # sigma = list(est=estimators$sigma, var=var_sigma)
      # 
      return(list(mu=estimators$mu, sigma=estimators$sigma, var_raw=var_raw, var_sandwich=sandwich_var))
    },
  error = function(cond) {
    message("pseudolikelihood failed")
    message(cond)
    return(NA)
  })
}

reduced_population = function(data) {
  out = tryCatch( 
    {
      mu = list()
      sigma = list()
      N = nrow(data)
      O = t(data) %*% data
      O = O[lower.tri(O, diag=TRUE)]
      
      ys = unique(data)
      p = ncol(data)
      theta = matrix(0, p * (p + 1) / 2, 1)
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
        print(paste("iteration", iter))
        theta = theta - solve(H)%*%grad
        
        ll = t(O) %*% theta - N * log(Z)
        D = sum(abs(theta-theta_old))
        if(D < sqrt(.Machine$double.eps)) {
          break
        }
      }

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
      
      mu$est = mu_ml
      sigma$est = theta_l/2
      
      mu$var = SE_mu
      sigma$var = fill_matrix(SE_sigma)
      
      return(list(mu=mu, sigma=sigma))
    },
    error = function(cond) {
      message("observed population method failed")
      message(cond)
      return(NA)
    }
  )
  
}

exact_likelihood = function(data, theta_0=matrix(0, p, p), max.nodes=20) {
  out = tryCatch({
    if (ncol(data) > max.nodes) {
      message("max number of variables for this method exceeded")
      return(NULL)
    }
    mu = sigma = list()
    N = nrow(data)
    O = t(data) %*% data
    diag(O) = diag(O)
    O = O[lower.tri(O, diag=TRUE)]
    
    p = ncol(data)
    theta = theta_0[lower.tri(theta_0, diag=TRUE)]
    ll = 1e10
    iteration = 0
    for(iter in 1:1e3) {
      ll_old = ll
      theta_old = theta
      
      E_suf = matrix(0, nrow = p*(p+1)/2, ncol=1)
      E_ss = matrix(0, nrow = p*(p+1)/2, ncol = p*(p+1)/2)
      
      
      Z = 0.0
      y = matrix(0, nrow = p, ncol = 1)
   
      Z = normalizingConstant(theta, E_suf, E_ss, Z, y)
      E_suf = E_suf
      E_ss = E_ss
      
      H = -N/Z * (E_ss - 1/Z*E_suf %*% t(E_suf))
      grad = O - N/Z *E_suf
      
      theta = theta - solve(H)%*%grad / 5
      ll = t(O) %*% theta - N * log(Z)
      D = sum(abs(ll - ll_old))
      iteration = iteration+1
      if(D < sqrt(.Machine$double.eps)) {
        break
      }
    }
    print(paste("solved in", iteration, "iterations"))
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
    
    mu$est = mu_ml
    mu$var = SE_mu
    sigma$est = theta_l/2
    sigma$var = SE_sigma
    
    return(list(mu=mu, sigma=sigma))
  },
  error = function(cond) {
    # estimators_sacha = psychonetrics::Ising(data, vars=1:p)
    # estimators = estimators_sacha %>% runmodel
    # ests = estimators@parameters[["est"]]
    # mu = ests[1:p]
    # sigvals = ests[-(1:p)]
    # sigvals = ests[-length(sigvals)]
    # sigma = matrix(0, p, p)
    # sigma[lower.tri(sigma)] = sigvals
    # sigma = (sigma + t(sigma)) / 2
    # ses = estimators@parameters[["se"]]
    # se_mu = ses[1:p]^2
    # se_sig = ses[-(1:p)]
    # se_sig = se_sig[-length(se_sig)]^2
    # 
    # mu = list(est=mu, var=se_mu)
    # sigma = list(est=sigma, var=se_sig)
    # return(list(mu=mu, sigma=sigma))
    
    message("exact likelihood failed")
    message(cond)
    return(NA)
  }
  )
}

maximum_likelihood = function(data) {
  p = ncol(data)
  model = psychonetrics::Ising(data)
  ml = model %>% runmodel
  
  results = transform_psych_results(ml, p)
  return(results)
}

logistic_regression = function(data) {
  out = tryCatch({
    p = ncol(data)
    N = nrow(data)
    sig = matrix(0, p, p)
    m = matrix(0, p, 1)
    varmu = matrix(0, p, 1)
    varsig = matrix(0, p, p)
    df = data.frame(data)
    for (i in 1:p) {
      params = log_reg(df, i) 
      vars = bootstrap_var(df, i)
      varsig[i, -i] = vars[-1]
      varmu[i] = vars[1]
      sig[i, -i] = params[-1]
      m[i] = params[1]
    }
    
    sig = symmetrizeMatrix(sig)
    varsig = symmetrizeMatrix(varsig)
    
    mu = list(est = m, var = varmu)
    sigma = list(est=sig, var = varsig)
    
    return(list(mu=mu, sigma=sigma))
    },
    error = function(cond){
      message("logistic regression failed")
      message(cond)
      return(NA)
    }
  )
}

estimate_ising = function(data, method, start_val=NULL) {
  if (method == "pseudo") 
    return(pseudolikelihood(data))
  else if (method == "exact"){
    if (is.null(start_val)) 
      return(exact_likelihood(data))
    else 
      return(exact_likelihood(data, start_val))
  }
  else if (method == "independent")
    return (logistic_regression(data))
  else if (method == "reduced") {
    if (is.null(start_val))
      return (reduced_population(data))
    else
      return(reduced_population(data, theta_0=start_val))
  }
  else
    message("This method is not included in the package")
} 

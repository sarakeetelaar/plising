pseudolikelihood = function(data) {
  # out = tryCatch({
    
    estimators = optimize_pseudolikelihood(data)
    var_pl = estimators$var
    
    p = ncol(data)
    
    var_sigma = matrix(0, nrow = p, ncol = p)
    
    var_mu = var_pl[1:p]
    var_sig = var_pl[-(1:p)]
    var_sigma[lower.tri(var_sigma)] = var_sig
    for (i in 1:p) {
      for (j in 1:p) {
        var_sigma[i, j] = var_sigma[j, i]
      }
    }
    
    mu = list(est=estimators$mu, var=var_mu)
    sigma = list(est=estimators$sigma, var=var_sigma)
    
    return(list(mu=mu, sigma=sigma))
  }
  # error = function(cond) {
  #   message("pseudolikelihood failed")
  #   message(cond)
  #   return(NA)
  # })
#}

reduced_population = function(data, theta_0=matrix(0, p, p)) {
  out = tryCatch( 
    {
      mu = list()
      sigma = list()
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
      sigma$est = theta_l
      
      mu$var = SE_mu
      sigma$var = SE_sigma
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
    O = t(data) %*% data * 2
    diag(O) = diag(O)/2
    O = O[lower.tri(O, diag=TRUE)]
    
    p = ncol(data)
    theta = theta_0[lower.tri(theta_0, diag=TRUE)]
    for(iter in 1:1e3) {
      print(paste("iteration", iter))
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
      print(paste("gradient", grad))
      print(paste("hessian", H))
      
      theta = theta - (1/2) * solve(H)%*%grad
      print(paste("theta = ", theta))
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
    mu$var = SE_mu
    sigma$est = theta_l
    sigma$var = SE_sigma
    
    return(list(mu=mu, sigma=sigma))
  },
  error = function(cond) {
    message("exact likelihood failed")
    message(cond)
    return(NA)
  }
  )
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
      thisy = df[,i]
      thisx = df[,-i]
      thisdata = data.frame(cbind(thisy, thisx))
      
      bootstr = boot::boot(data=thisdata, statistic=bootstrap_lr, R=1e3)
      vars = diag(cov(bootstr$t))
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

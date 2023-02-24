

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
    return(list(sigma_pl = sigma, mu_pl = mu, var_mu = se_mu, var_sigma = se_sigma))
  },
  error = function(cond) {
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
      
      return(list(likelihood=ll, mu_ML = mu_ml, sigma_ML = theta_l, var_mu = SE_mu, var_sigma = SE_sigma))
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
    
    return(list(mu_ML = mu_ml, sigma_ML = theta_l, var_mu = SE_mu, var_sigma = SE_sigma))
  },
  error = function(cond) {
    return(NA)
  }
  )
}



logistic_regression = function(data) {
  out = tryCatch(
    {
    p = ncol(data)
    N = nrow(data)
    sigma = matrix(0, p, p)
    mu = matrix(0, p, 1)
    varmu = matrix(0, p, 1)
    varsig = matrix(0, p, p)
    df = data.frame(data)
    for (i in 1:p) {
      params = log_reg(df, i)
      thisy = df[,i]
      thisx = df[,-i]
      thisdata = data.frame(cbind(thisy, thisx))
      
      bootstr = boot::boot(data=thisdata, statistic=bootstrap_lr, R=1e3)
      vars = diag(var(bootstr$t))
      
      varsig[i, -i] = vars[-1]
      varmu[i] = vars[1]
      sigma[i, -i] = params[-1]
      mu[i] = params[1]
    }
    sigma = symmetrizeMatrix(sigma)
    varsig = symmetrizeMatrix(varsig)
    theta = sigma
    diag_sigma = mu
    theta = theta[lower.tri(theta, diag=T)]
    return(list(mu=mu, sigma=sigma, theta=theta, var_mu=varmu, var_sigma=varsig))
    },
    error = function(cond){
      return(NA)
    }
  )
}

estimate_ising = function(data, method) {
  if (method == "pl") 
    return(pseudolikelihood(data))
  else if (method == "ml")
    return(exact_likelihood(data))
  else if (method == "lr")
    return (logistic_regression(data))
  else if (method =="h")
    return (hessen(data))
  else
    message("This method is not (yet) included in the package")
} 

compare_all = function(N, p, no.reps, mu_val, sigma_val, graph=NULL) {
  data = data_generation(N, p, no.reps, mu_val, sigma_val, A=graph) 
  all_X = data$X
  mu = data$mu
  sigma = data$sigma
  
  theta = sigma
  diag(theta) = mu
  theta = theta[lower.tri(theta, diag=T)]
  theta_PL = theta_ML = theta_H = theta_LR = matrix(0, p, p)
  dif_pl = dif_ml = dif_h = dif_lr = matrix(0, p, p)
  pl_old = ml_old = h_old = lr_old = NULL
  countPL = countML = countLR = countH = 0
  countPL2 = countML2 = countLR2 = countH2 = 0
  for (rep in 1:no.reps) {
    x = all_X[[rep]]
    pseudo = pseudolikelihood(x)
    
    if (is.list(pseudo)) {
      countPL = countPL+1
      theta_pl = pseudo$sigma_pl
      diag(theta_pl) = pseudo$mu_pl
      theta_PL = theta_PL + theta_pl
      if (!is.null(pl_old)) {
        countPL2 = countPL2 + 1
        
        dif_pl = dif_pl + (theta_pl - pl_old)^2
      }
      pl_old = theta_pl
    }
    
    
    likelihood = exact_likelihood(x, theta_pl)
    if (is.list(likelihood)){
      
      countML = countML+1
      theta_ml = likelihood$sigma_ML
      diag(theta_ml) = likelihood$mu_ML
      theta_ML = theta_ML + theta_ml
      if (!is.null(ml_old)) {
        dif_pl = dif_pl + (theta_ml - ml_old)^2
        countML2 = countML2 + 1
      }
      ml_old = theta_ml
    }
    hes = hessen(x, theta_pl)
    if(is.list(hes)){
      countH = countH + 1
      theta_h = hes$sigahes$mu_ML
      theta_H = theta_H + theta_h
      if (!is.null(h_old)) {
        dif_h = dif_h + (theta_h - h_old)^2
        countH2 = countH2 + 1
      }
      h_old = theta_h
    }
    
    logreg = logistic_regression(x)
    
    if (is.list(logreg)) {
      countLR = countLR + 1  
      theta_l = logreg$sigma/2
      diag(theta_l) = logreg$mu
      theta_LR = theta_LR + theta_l
      if (!is.null(lr_old)) {
        dif_lr = dif_lr + (theta_l - lr_old)^2
        countLR2 = countLR2 + 1
      }
      lr_old = theta_l
    }
  }
  dif_pl = dif_pl/(countPL2)
  dif_ml = dif_ml/(countML2)
  dif_h = dif_h/(countH2)
  dif_lr = dif_lr/(countLR2)
  
  theta_H = theta_H/countH
  theta_ML = theta_ML/countML
  theta_PL = theta_PL/countPL
  theta_LR = theta_LR/countLR
  
  theta_PL = theta_PL[lower.tri(theta_PL, diag=T)]
  theta_ML = theta_ML[lower.tri(theta_ML, diag=T)]
  theta_H = theta_H[lower.tri(theta_H, diag=T)]
  theta_LR = theta_LR[lower.tri(theta_LR, diag=T)]
  
  print(paste("For N =", N, "and p =", p, "the number of succeeded replications for the 
              pseudolikelihood is", countPL, ", for the exact likelihood", countML,
              "for Hessen", countH, "and for logistic regression", countLR))
  return(list(pl=theta_PL, ml=theta_ML, hes=theta_H, lr=theta_LR, theta=theta,
              pl_var=dif_pl, ml_var=dif_ml, hes_var=dif_h, lr_var=dif_lr))
}

#for full graph
comparison_complete_graph = function(n_list, p_list) {
  thetapl = thetaml = thetah = thetalr = list()
  diffpl = diffml = diffhes = difflr = list()
  theta_true = list()
  
  for (ip in 1:length(p_list)) {
    thetaPL = thetaML = thetaH = thetaLR = list()
    diffPL = diffML = diffH = diffLR = list()
    
    p = p_list[ip]
    
    mu = runif(p, -1,1)
    sig_val = runif(p*(p-1)/2, -1, 1)/2
    theta = matrix(0, p,p)
    theta[lower.tri(theta)] = sig_val
    diag(theta) = mu
    theta = theta[lower.tri(theta, diag=T)]
    theta_true[[ip]] = theta
    for (i_n in 1:length(N_list)) {
      
      n = N_list[i_n]
      results = compare_all(N=n, p=p, no.reps=100, mu=mu, sigma_val=sig_val)
      thetaPL[[i_n]] = results$pl
      thetaML[[i_n]] = results$ml
      thetaH[[i_n]] = results$hes
      thetaLR[[i_n]] = results$lr
      
      diffPL[[i_n]] = results$pl_var
      diffML[[i_n]] = results$ml_var
      diffH[[i_n]] = results$h_var
      diffLR[[i_n]] = results$lr_var
    }
    
    diffpl[[ip]] = diffPL
    diffml[[ip]] = diffML
    diffhes[[ip]] = diffH
    difflr[[ip]] = diffLR
    
    thetapl[[ip]] = thetaPL
    thetaml[[ip]] = thetaML
    thetah[[ip]] = thetaH
    thetalr[[ip]] = thetaLR
  }
  return(theta=theta_true, est_pl=thetapl, est_ml=thetaml, est_hes=thetah, est_lr=thetalr,
         var_pl=diffpl, var_ml=diffml, var_hes=diffhes, var_lr=difflr)
}

comparison_random_graph = function(n_list, p_list, degree) {
  thetapl = thetaml = thetah = thetalr = list()
  diffpl = diffml = diffhes = difflr = list()
  theta_true = list()
  
  for (ip in 1:length(p_list)) {
    thetaPL = thetaML = thetaH = thetaLR = list()
    diffPL = diffML = diffH = diffLR = list()
    
    p = p_list[ip]
    
    mu = runif(p, -1,1)
    sig_val = runif(p*(p-1)/2, -1, 1)/2
    theta = matrix(0, p,p)
    theta[lower.tri(theta)] = sig_val
    diag(theta) = mu
    theta = theta[lower.tri(theta, diag=T)]
    theta_true[[ip]] = theta
    A = random_network(p=p, m=p*degree)
    
    for (i_n in 1:length(N_list)) {
      
      n = N_list[i_n]
      results = compare_all(N=n, p=p, no.reps=100, mu=mu, sigma_val=sig_val, graph=A)
      thetaPL[[i_n]] = results$pl
      thetaML[[i_n]] = results$ml
      thetaH[[i_n]] = results$hes
      thetaLR[[i_n]] = results$lr
      
      diffPL[[i_n]] = results$pl_var
      diffML[[i_n]] = results$ml_var
      diffH[[i_n]] = results$h_var
      diffLR[[i_n]] = results$lr_var
    }
    
    diffpl[[ip]] = diffPL
    diffml[[ip]] = diffML
    diffhes[[ip]] = diffH
    difflr[[ip]] = diffLR
    
    thetapl[[ip]] = thetaPL
    thetaml[[ip]] = thetaML
    thetah[[ip]] = thetaH
    thetalr[[ip]] = thetaLR
  }
  return(theta=theta_true, est_pl=thetapl, est_ml=thetaml, est_hes=thetah, est_lr=thetalr,
         var_pl=diffpl, var_ml=diffml, var_hes=diffhes, var_lr=difflr)
}

comparison_small_world = function(n_list, p_list, nei, pr) {
  thetapl = thetaml = thetah = thetalr = list()
  diffpl = diffml = diffhes = difflr = list()
  theta_true = list()
  
  for (ip in 1:length(p_list)) {
    thetaPL = thetaML = thetaH = thetaLR = list()
    diffPL = diffML = diffH = diffLR = list()
    
    p = p_list[ip]
    
    mu = runif(p, -1,1)
    sig_val = runif(p*(p-1)/2, -1, 1)/2
    theta = matrix(0, p,p)
    theta[lower.tri(theta)] = sig_val
    diag(theta) = mu
    theta = theta[lower.tri(theta, diag=T)]
    theta_true[[ip]] = theta
    
    A = create_small_world(p, nei, pr)
    for (i_n in 1:length(N_list)) {
      
      n = N_list[i_n]
      results = compare_all(N=n, p=p, no.reps=100, mu=mu, sigma_val=sig_val, graph=A)
      thetaPL[[i_n]] = results$pl
      thetaML[[i_n]] = results$ml
      thetaH[[i_n]] = results$hes
      thetaLR[[i_n]] = results$lr
      
      diffPL[[i_n]] = results$pl_var
      diffML[[i_n]] = results$ml_var
      diffH[[i_n]] = results$h_var
      diffLR[[i_n]] = results$lr_var
    }
    
    diffpl[[ip]] = diffPL
    diffml[[ip]] = diffML
    diffhes[[ip]] = diffH
    difflr[[ip]] = diffLR
    
    thetapl[[ip]] = thetaPL
    thetaml[[ip]] = thetaML
    thetah[[ip]] = thetaH
    thetalr[[ip]] = thetaLR
  }
  return(theta=theta_true, est_pl=thetapl, est_ml=thetaml, est_hes=thetah, est_lr=thetalr,
         var_pl=diffpl, var_ml=diffml, var_hes=diffhes, var_lr=difflr)
}

comparison_empty_graph = function(n_list, p_list) {
  thetapl = thetaml = thetah = thetalr = list()
  diffpl = diffml = diffhes = difflr = list()
  theta_true = list()
  
  for (ip in 1:length(p_list)) {
    thetaPL = thetaML = thetaH = thetaLR = list()
    diffPL = diffML = diffH = diffLR = list()
    
    p = p_list[ip]
    
    A = matrix(0, p, p)
    mu = runif(p, -1,1)
    sig_val = runif(p*(p-1)/2, -1, 1)/2
    theta = matrix(0, p,p)
    theta[lower.tri(theta)] = sig_val
    diag(theta) = mu
    theta = theta[lower.tri(theta, diag=T)]
    theta_true[[ip]] = theta
    for (i_n in 1:length(N_list)) {
      
      n = N_list[i_n]
      results = compare_all(N=n, p=p, no.reps=100, mu=mu, sigma_val=sig_val, graph=A)
      thetaPL[[i_n]] = results$pl
      thetaML[[i_n]] = results$ml
      thetaH[[i_n]] = results$hes
      thetaLR[[i_n]] = results$lr
      
      diffPL[[i_n]] = results$pl_var
      diffML[[i_n]] = results$ml_var
      diffH[[i_n]] = results$h_var
      diffLR[[i_n]] = results$lr_var
    }
    
    diffpl[[ip]] = diffPL
    diffml[[ip]] = diffML
    diffhes[[ip]] = diffH
    difflr[[ip]] = diffLR
    
    thetapl[[ip]] = thetaPL
    thetaml[[ip]] = thetaML
    thetah[[ip]] = thetaH
    thetalr[[ip]] = thetaLR
  }
  return(theta=theta_true, est_pl=thetapl, est_ml=thetaml, est_hes=thetah, est_lr=thetalr,
         var_pl=diffpl, var_ml=diffml, var_hes=diffhes, var_lr=difflr)
}
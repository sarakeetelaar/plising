compare_all = function(N, p, no.reps, mu_val, sigma_val, graph=NULL) {
  data = data_generation(N, p, no.reps, mu_val, sigma_val, A=graph) 
  all_X = data$X
  
  mu = data$mu
  sigma = data$sigma
  
  theta = sigma
  diag(theta) = mu
  theta = theta[lower.tri(theta, diag=T)]
  for (rep in 1:no.reps) {
    x = X[[rep]]
    pseudo = pseudolikelihood(x)
    
    if (is.list(pseudo)) {
      countPL = countPL+1
      theta_pl = pseudo$sigma_pl
      diag(theta_pl) = pseudo$mu_pl
      print("pseudolikelihood finished")
      theta_PL = theta_PL + theta_pl
      if (exists("pl_old")) {
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
      print("Exact likelihood finished")
      theta_ML = theta_ML + theta_ml
      if (exists("ml_old")) {
        dif_pl = dif_pl + (theta_ml - ml_old)^2
        countML2 = countML2 + 1
      }
      ml_old = theta_ml
    }
    hes = hessen(x, theta_pl)
    if(is.list(hes)){
      countH = countH + 1
      theta_h = hes$sigma_ML
      diag(theta_h) = hes$mu_ML
      print("hessen method finished")
      theta_H = theta_H + theta_h
      if (exists("h_old")) {
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
      print("logistic regression finished")
      theta_LR = theta_LR + theta_l
      if (exists("lr_old")) {
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
  
  return(list(pl=theta_PL, ml=theta_ML, hes=theta_H, lr=theta_LR, theta=theta,
              pl_var=dif_pl, ml_var=dif_ml, hes_var=dif_h, lr_var=dif_lr))
}




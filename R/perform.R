compare_all = function(N, p, no.reps, mu_val, sigma_val, graph=NULL) {
  mu_index = re_index(p)
  data = data_generation(N, p, no.reps, A=graph) 
  mu = data$mu
  sigma = data$sigma
  all_X = data$X
  theta = sigma
  diag(theta) = mu
  theta = theta[lower.tri(theta, diag=T)]
  theta_PL = theta_ML = theta_H = theta_LR = var_LR = matrix(0, p, p)
  var_PL = var_ML = var_H= matrix(0, p*(p+1)/2, 1)
  dif_pl = dif_ml = dif_h = dif_lr = matrix(0, p, p)
  pl_old = ml_old = h_old = lr_old = NULL
  countPL = countML = countLR = countH = 0
  countPL2 = countML2 = countLR2 = countH2 = 0
  for (rep in 1:no.reps) {
    x = all_X[[rep]]
    pseudo = pseudolikelihood(x)
    theta_pl = matrix(0, p, p)
    if (is.list(pseudo)) {
      countPL = countPL+1
      theta_pl = pseudo$sigma_pl
      var_pl = rep(0, p*(p+1)/2)
      var_pl[-mu_index] = pseudo$var_sigma
      var_pl[mu_index] = pseudo$var_mu
      diag(theta_pl) = pseudo$mu_pl
      theta_PL = theta_PL + theta_pl
      var_PL = var_PL + var_pl
    }
    
    
    likelihood = exact_likelihood(x, theta_pl)
    if (is.list(likelihood)){
      
      countML = countML+1
      theta_ml = likelihood$sigma_ML
      diag(theta_ml) = likelihood$mu_ML
      
      var_ml = matrix(0, p*(p+1)/2, 1)
      var_ml[-mu_index] = likelihood$var_sigma
      var_ml[mu_index] = likelihood$var_mu
      theta_ML = theta_ML + theta_ml
      var_ML = var_ML + var_ml
      if (!is.null(ml_old)) {
        dif_pl = dif_pl + (theta_ml - ml_old)^2
        countML2 = countML2 + 1
      }
      ml_old = theta_ml
    }
    hes = hessen(x, theta_pl)
    if(is.list(hes)){
      countH = countH + 1
      theta_h = hes$sigma_ML
      diag(theta_H) = hes$mu_ML
      
      var_h = matrix(0, p*(p+1)/2, 1)
      var_h[-mu_index] = hes$var_sigma
      var_h[mu_index] = hes$var_mu
      
      theta_H = theta_H + theta_h
      var_H = var_H + var_h
      h_old = theta_h
    }
    
    logreg = logistic_regression(x)
    
    if (is.list(logreg)) {
      countLR = countLR + 1  
      theta_l = logreg$sigma/2
      diag(theta_l) = logreg$mu
      
      var_lr = logreg$var_sigma
      diag(var_lr) = logreg$var_mu
      
      theta_LR = theta_LR + theta_l
      var_LR = var_LR + var_lr
      
      lr_old = theta_l
    }
  }
  if (countH != 0){
    theta_H = theta_H/countH
    var_H = var_H/countH
  }
  else
    theta_H = NA
  theta_ML = theta_ML/countML
  var_ML = var_ML/countML
  
  theta_PL = theta_PL/countPL
  var_PL = var_PL/countPL
  
  theta_LR = theta_LR/countLR
  var_LR = var_LR/countLR2
  
  theta_PL = theta_PL[lower.tri(theta_PL, diag=T)]
  
  theta_ML = theta_ML[lower.tri(theta_ML, diag=T)]
  
  theta_H = theta_H[lower.tri(theta_H, diag=T)]
  
  theta_LR = theta_LR[lower.tri(theta_LR, diag=T)]
  var_LR = var_LR[lower.tri(var_LR, diag=T)]
  
  print(paste("For N =", N, "and p =", p, "the number of succeeded replications for the 
              pseudolikelihood is", countPL, ", for the exact likelihood", countML,
              "for Hessen", countH, "and for logistic regression", countLR))
  return(list(pl=theta_PL, ml=theta_ML, hes=theta_H, lr=theta_LR, theta=theta,
              var_pl=var_PL, var_ml=var_ML, var_hes=var_H, var_lr=var_LR))
}

#for full graph
comparison_complete_graph = function(n_list, p_list, no.reps) {
  thetapl = thetaml = thetah = thetalr = list()
  varpl = varml = varhes = varlr = list()
  theta_true = list()
  
  for (ip in 1:length(p_list)) {
    thetaPL = thetaML = thetaH = thetaLR = list()
    varPL = varML = varH = varLR = list()
    
    p = p_list[ip]
    
    
    sig_val = runif(p*(p-1)/2, -1, 1)/2
    theta = matrix(0, p,p)
    theta[lower.tri(theta)] = sig_val
    diag(theta) = mu
    theta = theta[lower.tri(theta, diag=T)]
    theta_true[[ip]] = theta
    for (i_n in 1:length(n_list)) {
      
      n = n_list[i_n]
      results = compare_all(N=n, p=p, no.reps=no.reps, mu=mu, sigma_val=sig_val)
      thetaPL[[i_n]] = results$pl
      thetaML[[i_n]] = results$ml
      thetaH[[i_n]] = results$hes
      thetaLR[[i_n]] = results$lr
      
      varPL[[i_n]] = results$var_pl
      varML[[i_n]] = results$var_ml
      varH[[i_n]] = results$var_hes
      varLR[[i_n]] = results$var_lr
    }
    
    varpl[[ip]] = varPL
    varml[[ip]] = varML
    varhes[[ip]] = varH
    varlr[[ip]] = varLR
    
    thetapl[[ip]] = thetaPL
    thetaml[[ip]] = thetaML
    thetah[[ip]] = thetaH
    thetalr[[ip]] = thetaLR
  }
  return(list(theta=theta_true, est_pl=thetapl, est_ml=thetaml, est_hes=thetah, est_lr=thetalr,
         var_pl=varpl, var_ml=varml, var_hes=varhes, var_lr=varlr))
}

comparison_random_graph = function(n_list, p_list, degree, no.reps) {
  thetapl = thetaml = thetah = thetalr = list()
  varpl = varml = varhes = varlr = list()
  theta_true = list()
  graphs = list()
  for (ip in 1:length(p_list)) {
    thetaPL = thetaML = thetaH = thetaLR = list()
    varPL = varML = varH = varLR = list()
    
    p = p_list[ip]
    
    mu = runif(p, -1,1)
    sig_val = runif(p*(p-1)/2, -1, 1)/2
    A = random_network(p=p, m=p*degree)
    graphs[[ip]] = graph_from_adjacency_matrix(A, mode="undirected")
    for (i_n in 1:length(n_list)) {
      
      n = n_list[i_n]
      results = compare_all(N=n, p=p, no.reps=no.reps, mu=mu, sigma_val=sig_val, graph=A)
      thetaPL[[i_n]] = results$pl
      thetaML[[i_n]] = results$ml
      thetaH[[i_n]] = results$hes
      thetaLR[[i_n]] = results$lr
      
      varPL[[i_n]] = results$var_pl
      varML[[i_n]] = results$var_ml
      varH[[i_n]] = results$var_hes
      varLR[[i_n]] = results$var_lr
      
      if (i_n == 1) {
        theta_true[[ip]] = results$theta
      }
    }
    
    varpl[[ip]] = varPL
    varml[[ip]] = varML
    varhes[[ip]] = varH
    varlr[[ip]] = varLR
    
    thetapl[[ip]] = thetaPL
    thetaml[[ip]] = thetaML
    thetah[[ip]] = thetaH
    thetalr[[ip]] = thetaLR
  }
  return(list(theta=theta_true, est_pl=thetapl, est_ml=thetaml, est_hes=thetah, est_lr=thetalr,
         var_pl=varpl, var_ml=varml, var_hes=varhes, var_lr=varlr, graphs=graphs))
}

comparison_small_world = function(n_list, p_list, nei, pr, no.reps) {
  thetapl = thetaml = thetah = thetalr = list()
  varpl = varml = varhes = varlr = list()
  theta_true = graphs = list()
  
  for (ip in 1:length(p_list)) {
    thetaPL = thetaML = thetaH = thetaLR = list()
    varPL = varML = varH = varLR = list()
    
    p = p_list[ip]
    
    mu = runif(p, -1,1)
    sig_val = runif(p*(p-1)/2, -1, 1)/2
    
    A = create_small_world(p, nei, pr)
    graphs[[ip]] = graph_from_adjacency_matrix(A, mode="undirected")
    for (i_n in 1:length(n_list)) {
      
      n = n_list[i_n]
      results = compare_all(N=n, p=p, no.reps=no.reps, mu=mu, sigma_val=sig_val, graph=A)
      thetaPL[[i_n]] = results$pl
      thetaML[[i_n]] = results$ml
      thetaH[[i_n]] = results$hes
      thetaLR[[i_n]] = results$lr
      
      varPL[[i_n]] = results$var_pl
      varML[[i_n]] = results$var_ml
      varH[[i_n]] = results$var_hes
      varLR[[i_n]] = results$var_lr
      
      if (i_n == 1) {
        theta_true[[ip]] = results$theta
      }
    }
    
    varpl[[ip]] = varPL
    varml[[ip]] = varML
    varhes[[ip]] = varH
    varlr[[ip]] = varLR
    
    thetapl[[ip]] = thetaPL
    thetaml[[ip]] = thetaML
    thetah[[ip]] = thetaH
    thetalr[[ip]] = thetaLR
  }
  return(list(theta=theta_true, est_pl=thetapl, est_ml=thetaml, est_hes=thetah, est_lr=thetalr,
         var_pl=varpl, var_ml=varml, var_hes=varhes, var_lr=varlr, graphs=graphs))
}

comparison_empty_graph = function(n_list, p_list, no.reps) {
  thetapl = thetaml = thetah = thetalr = list()
  varpl = varml = varhes = varlr = list()
  theta_true = list()
  
  for (ip in 1:length(p_list)) {
    thetaPL = thetaML = thetaH = thetaLR = list()
    varPL = varML = varH = varLR = list()
    
    p = p_list[ip]
    
    A = matrix(0, p, p)
    mu = runif(p, -1,1)
    sig_val = runif(p*(p-1)/2, -1, 1)/2
    for (i_n in 1:length(n_list)) {
      
      n = n_list[i_n]
      results = compare_all(N=n, p=p, no.reps=no.reps, mu=mu, sigma_val=sig_val, graph=A)
      thetaPL[[i_n]] = results$pl
      thetaML[[i_n]] = results$ml
      thetaH[[i_n]] = results$hes
      thetaLR[[i_n]] = results$lr
      
      varPL[[i_n]] = results$var_pl
      varML[[i_n]] = results$var_ml
      varH[[i_n]] = results$var_hes
      varLR[[i_n]] = results$var_lr
      
      if (i_n == 1) {
        theta_true[[ip]] = results$theta
      }
    }
    
    varpl[[ip]] = varPL
    varml[[ip]] = varML
    varhes[[ip]] = varH
    varlr[[ip]] = varLR
    
    thetapl[[ip]] = thetaPL
    thetaml[[ip]] = thetaML
    thetah[[ip]] = thetaH
    thetalr[[ip]] = thetaLR
  }
  return(list(theta=theta_true, est_pl=thetapl, est_ml=thetaml, est_hes=thetah, est_lr=thetalr,
         var_pl=varpl, var_ml=varml, var_hes=varhes, var_lr=varlr))
}

comparison_plots = function(results, pind, nind, poptions, noptions) {
  p = poptions[pind]
  n = noptions[nind]
  
  pl = results$est_pl[[pind]][[nind]]
  ml = results$est_ml[[pind]][[nind]]
  hes = results$est_hes[[pind]][[nind]]
  lr = results$est_lr[[pind]][[nind]]
  real = results$theta[[pind]]
  
  muind = re_index(p)
  cols = rep('black', length(pl))
  cols[muind] = 'red'
  par(mfrow=c(2,2), xpd=F)
  
  plot(pl, real, xlab="PL", ylab="theta", main="true value", col=cols);abline(0,1)
  mtext(paste("p =", p, ", N =", N), side=3, line=-2, outer=T)
  plot(pl, lr, xlab="PL", ylab="LR", main = "logistic regression", col=cols);abline(0,1)
  tryCatch({
  plot(pl, ml, xlab="PL", ylab="ML", main="exact likelihood", col=cols);abline(0,1)},
  error = function(cond){
    message(cond)
  }
  )

  tryCatch( {
    plot(pl, hes, xlab="PL", ylab="hessen", main="hessen", col=cols);abline(0,1)},
    error=function(cond) {
      message(cond)
    }
  )

}

print_bias = function(results, p_list, n_list) {
  for (indp in 1:length(p_list)) {
    thisp = p_list[indp]
    true_theta = results$theta[[indp]]
    mu_ind = re_index(thisp)
    true_mu = true_theta[mu_ind]
    true_sig = true_theta[-mu_ind]
    
    for (indn in 1:length(n_list)) {
      pl = results$est_pl[[indp]][[indn]]
      ml = results$est_ml[[indp]][[indn]]
      h = results$est_hes[[indp]][[indn]]
      lr = results$est_lr[[indp]][[indn]]
      
      mupl = pl[mu_ind]
      sigpl = pl[-mu_ind]
      muml = ml[mu_ind]
      sigml = ml[-mu_ind]
      muh = h[mu_ind]
      sigh = h[-mu_ind]
      mulr = lr[mu_ind]
      siglr = lr[-mu_ind]
      
      biasmupl = mean((mupl-true_mu)^2)
      biassigpl = mean((sigpl-true_sig)^2)
      biasmuml = mean((muml-true_mu)^2)
      biassigml = mean((sigml-true_sig)^2)
      biasmuh = mean((muh-true_mu)^2)
      biassigh = mean((sigh-true_sig)^2)
      biasmulr = mean((mulr-true_mu)^2)
      biassiglr = mean((siglr-true_sig)^2)
      print(paste("N=", n_list[indn], "p=", p_list[indp]))
      print(paste("the squared bias for PL: ", biasmupl, "sig: ", biassigpl))
      print(paste("the squared bias for ML: ", biasmuml, "sig: ", biassigml))
      print(paste("the squared bias for H: ", biasmuh, "sig: ", biassigh))
      print(paste("the squared bias for LR: ", biasmulr, "sig: ", biassiglr))
    }
  }
  
}


get_bias = function(results, p_ind, p_options, n_options) {
  len_n = length(n_options)
  p = p_options[p_ind]
  true_theta = results$theta[[p_ind]]
  mu_index = re_index(p)
  mu = true_theta[mu_index]
  sig = true_theta[-mu_index]
  
  pl = results$est_pl[[p_ind]]
  ml = results$est_ml[[p_ind]]
  h = results$est_hes[[p_ind]]
  lr = results$est_lr[[p_ind]]
  
  bias_sigpl = bias_sigml = bias_sigh = bias_siglr = rep(0, len_n)
  bias_mupl = bias_muml = bias_muh = bias_mulr = rep(0, len_n)
  for (i in 1:len_n) {
    thispl = pl[[i]]
    bias_sigpl[i] = mean((thispl[-mu_index]-sig)^2)
    bias_mupl[i] = mean((thispl[mu_index]-mu)^2)
    
    thisml = ml[[i]]
    bias_sigml[i] = mean((thisml[-mu_index] - sig)^2)
    bias_muml[i] = mean((thisml[mu_index]-mu)^2)
    
    thish = h[[i]]
    bias_sigh[i] = mean((thish[-mu_index] - sig)^2)
    bias_muh[i] = mean((thish[mu_index] - mu)^2)
    
    thislr = lr[[i]]
    bias_siglr[i] = mean((thislr[-mu_index] -sig)^2)
    bias_mulr[i] = mean((thislr[mu_index] - mu)^2)
  }
  mulist = c(bias_mupl, bias_muml, bias_muh, bias_mulr)
  siglist = c(bias_sigpl, bias_sigml, bias_sigh, bias_siglr)
  
  range_mu=c(min(mulist, na.rm=T), max(mulist, na.rm=T))
  range_sig = c(min(siglist, na.rm=T), max(siglist, na.rm=T))
  return(list(bias_mupl=bias_mupl, bias_sigpl=bias_sigpl, bias_sigml=bias_sigml, 
              bias_muml=bias_muml, bias_sigh=bias_sigh,
              bias_muh=bias_muh, bias_mulr=bias_mulr, bias_siglr=bias_siglr,
              range_mu=range_mu, range_sigma=range_sig))
}

plot_bias = function(results, p_options, p_ind, n_options) {
  p = p_options[p_ind]
  par(mfrow=c(2, 1))
  bias = get_bias(results, p_ind, p_options, n_options)
  plot(n_options, bias$bias_mulr, main=paste("Mu, p =", p), col="orange",type="l", xlab="N", ylab= "squared bias", ylim=bias$range_mu)
  lines(n_options, bias$bias_muml, col="red")
  lines(n_options, bias$bias_muh, col="darkgreen")
  lines(n_options, bias$bias_mupl, col="blue")
  legend(900, bias$range_mu[2], legend=c("PL", "ML", "Hessen", "Log Reg"), col=c("blue", "red", "darkgreen", "orange"), lty=1, cex=.8)
  
  plot(n_options, bias$bias_siglr, main=paste("Sigma, p =", p), col="orange", type="l", xlab="N", ylab="squared bias", ylim=bias$range_sigma)
  lines(n_options, bias$bias_sigml, col="red")
  lines(n_options, bias$bias_sigh, col='darkgreen')
  lines(n_options, bias$bias_sigpl, col='blue')
  legend(900, bias$range_sigma[2], legend=c("PL", "ML", "Hessen", "Log Reg"), col=c("blue", "red", "darkgreen", "orange"), lty=1, cex=.8)
}



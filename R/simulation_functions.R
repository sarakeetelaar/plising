complete_results = function(p_list, n_options, no.reps) {
  results =  list()
  for (p in p_list) {
    params = parameter_generation(p=p)
    
    sim_results = perform_simulation(params, n_options, no.reps)
    
    results[[p - p_list[1] + 1]] = list(p=p, parameters=params, results=sim_results)
  }
  
  return(results)
}

random_results = function(p_list, n_options, no.reps) {
  results = list()
  
  for (p in p_list) {
    gr = random_graph(p=p)
    
    params = parameter_generation(p=p, graph=gr)
    
    sim_results = perform_simulation(params, n_options, no.reps)
    
    results[[p - p_list[1] + 1]] = list(p=p, parameters=params, results=sim_results)
  }
  
  return(results)
}

small_world_results = function(p_list, n_options, no.reps) {
  results = list()
  
  for (p in p_list) {
    gr = small_world(p=p)
    
    params = parameter_generation(p=p, graph=gr)
    
    sim_results = perform_simulation(params, n_options, no.reps)
    
    results[[p - p_list[1] + 1]] = list(p=p, parameters=params, results=sim_results)
  }
  
  return(results)
}
perform_simulation = function(parameters, n_options, no.reps) {
  mu = parameters$mu
  sigma = parameters$sigma
  p = length(mu)
  
  results = list()
  
  for (ind_n in 1:length(n_options)) {
    n = n_options[ind_n]

    pl_sigma_est = ml_sigma_est = lr_sigma_est = red_sigma_est = matrix(0, p, p)
    pl_sigma_var = ml_sigma_var = lr_sigma_var = red_sigma_var = matrix(0, p, p)
    pl_mu_est = ml_mu_est = lr_mu_est = red_mu_est = matrix(0, p, 1)
    pl_mu_var = ml_mu_var = lr_mu_var = red_mu_var = matrix(0, p, 1)
    countML = countPL = 0
    countR = 0
    for (rep in 1:no.reps) {
      x = full_data_generation(parameters, n)
      
      pl = pseudolikelihood(x)
      if (is.list(pl)) {
        pl_sigma_est = pl_sigma_est + pl$sigma$est
        pl_mu_est = pl_mu_est + pl$mu$est
        pl_sigma_var = pl_sigma_var = pl$sigma$var
        pl_mu_var = pl_mu_var + pl$mu$var
        
        countPL = countPL + 1
      }
      
      ml = maximum_likelihood(x)
      if (is.list(ml)) {
        ml_sigma_est = ml_sigma_est + ml$sigma$est/2
        ml_mu_est = ml_mu_est + ml$mu$est
        ml_sigma_var = ml_sigma_var + ml$sigma$var
        ml_mu_var = ml_mu_var + ml$mu$var
        countML = countML + 1
      }
      
      lr = logistic_regression(x)
      lr_sigma_est = lr_sigma_est + lr$sigma$est
      lr_mu_est = lr_mu_est + lr$mu$est
      lr_sigma_var = lr_sigma_var + lr$sigma$var
      lr_mu_var = lr_mu_var + lr$mu$var
      # 
      # reduced = reduced_population(x)
      # if (is.list(reduced)) {
      #   red_sigma_est = red_sigma_est + reduced$sigma$est
      #   red_mu_est = red_mu_est + reduced$mu$est
      #   red_sigma_var = red_sigma_var + reduced$sigma$var
      #   red_mu_var = red_mu_var + reduced$mu$var
      #   countR = countR + 1
      # }
      
      
    }
    pl = list(sigma=list(est=pl_sigma_est/countPL, var=pl_sigma_var/countPL),
              mu = list(est=pl_mu_est/countPL, var=pl_mu_var/countPL))
    ml = list(sigma=list(est=ml_sigma_est/countML, var=ml_sigma_var/countML),
              mu = list(est=ml_mu_est/countML, var=ml_mu_var/countML))
    lr = list(sigma=list(est=lr_sigma_est/no.reps, var=lr_sigma_var/no.reps),
              mu = list(est=lr_mu_est/no.reps, var=lr_mu_var/no.reps))
    # if (countR != 0)
    #   reduced = list(sigma=list(est=red_sigma_est/countR, var=red_sigma_var/countR),
    #                  mu=list(est=red_mu_est/countR, var=red_mu_var/countR))
    # else
    #   reduced=NA
    results[[ind_n]] = list(N=n, pl=pl, ml=ml, lr=lr)
  }
  results$parameters = parameters
  return(results)
}
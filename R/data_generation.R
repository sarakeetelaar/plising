data_generation = function(N=1000, graph, mu, sigma=1) {
  p = ncol(graph)
  X = matrix(0, nrow=N, ncol=p)
  for(iter in 1:1e3) {
    for(i in 1:p) {
      X[, i] = 1 * (rlogis(N) <= mu[i] + 2*sigma*X %*% graph[, i])
    }
  }
  return(X)
}

estimate_full_parameters = function(data) {
  data = data.matrix(data)
  est = IsingSampler::EstimateIsing(data, responses=c(0L, 1L), method="pl")
  
  return(est)
}


subset_parameters = function(est, p) {
  full_sigma =  est$graph
  full_mu = est$thresholds
  
  selected_vars = sample(1:length(full_mu), p)
  sigma = full_sigma[selected_vars, selected_vars]
  mu = full_mu[selected_vars]
  
  return(list(mu=mu, sigma=sigma/2))
}
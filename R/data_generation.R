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

estimate_parameters = function(data, N, p) {
  options = 6:41
  rownumbers = sample(1:26715, N, replace=F)
  colnumbers = sample(6:41, p, replace=F)
  sub_data = data[rownumbers, colnumbers]
  df = apply(sub_data, 2, strtoi)
  print(df)
  estimate = estimate_ising(df, "pseudo")
  sigma = estimate$sigma$est
  mu = estimate$mu$est
  return(list(sigma=sigma, mu=mu))
}

process_data = function(data, N, p) {
  data = data.matrix(data)
  rows = sample(1:nrow(data), N)
  cols = sample(1:ncol(data), p)
  subset = data[rows, cols]
  
  df = apply(subset, 2, dichotomize)
}

dichotomize = function(el) {
  if (el < 2)
    return(0)
  return(1)
}

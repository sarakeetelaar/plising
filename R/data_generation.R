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

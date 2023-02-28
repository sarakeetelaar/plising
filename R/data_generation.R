data_generation = function(N=1000, graph, mu, sigma=1, no.reps=1e3) {
  data = list()
  p = ncol(graph)
  for (r in 1: no.reps) {
    X = matrix(0, nrow=N, ncol=p)
    for(iter in 1:1e2) {
      for(i in 1:p) {
        X[, i] = 1 * (rlogis(N) <= mu[i] + 2*sigma*X %*% graph[, i])
      }
    }
    data[[r]] = X
  }
  return(data)
}

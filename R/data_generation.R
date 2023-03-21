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

data_preparation = function(data=Wenchuan) {
  df = na.omit(data)
  df = ifelse(df < 3, 0, 1)
  
  return(df)
}

parameter_generation = function(p, df_prep=wenchuan[, 1:p],  graph=matrix(1, p, p)) {
  res = IsingSampler::EstimateIsing(df_prep, responses=c(0L, 1L), method="uni", adj=graph)
  
  sigma = res$graph/2
  mu = res$thresholds
  
  return(list(mu=mu, sigma=sigma))
}

small_world = function(p) {
  sw = igraph::sample_smallworld(1, p, 2, .2)
  adj = as.matrix(igraph::as_adjacency_matrix(sw))
  
  return(adj)
}

random_graph = function(p, pr=.3) {
  rg = igraph::erdos.renyi.game(p, pr)
  
  adj = as.matrix(igraph::as_adjacency_matrix(rg))
  
}
#combines the methods above to generate data given p and N and graph structure

full_data_generation = function(parameters, N) {
  gr = parameters$sigma
  mu = parameters$mu
  
  x = data_generation(N=N, graph=gr, mu=mu)
  return(x)
}

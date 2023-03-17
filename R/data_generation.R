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

parameter_generation = function(df_prep, p) {
  df = df_prep[, 1:p]
  res = pseudolikelihood(df)
  
  sigma = res$sigma$est
  mu = res$mu$est
  
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
full_data_generation = function(data=Wenchuan, p, N, graph=NULL) {
  
  if (is.null(graph)) {
    graph = matrix(1, p, p)
  }
  prep_df = data_preparation(data)
  params = parameter_generation(prep_df, p)
  
  gr = params$sigma * graph
  mu = params$mu
  x = data_generation(N = N, graph = gr, mu = mu)
  return(list(x=x, mu=mu, sigma=gr))
}

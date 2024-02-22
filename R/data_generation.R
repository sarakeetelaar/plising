###############################################################################
# File contains all functions related to the data generation process of the 
# simulation
##############################################################################

# generates Ising data using input parameters and population size N
full_data_generation = function(parameters, N) {
  
  mu <- parameters$mu
  sigma <- parameters$sigma
  p <- ncol(sigma)
  X <- matrix(0, nrow = N, ncol = p)
  # generate Ising data
  for(iter in 1:1e3) {
    for(i in 1:p) {
      X[, i] <- 1 * (rlogis(N) <= mu[i] + 2 * X %*% sigma[, i])
    }
  }
  return(X)
}

# function dichotomizes and removes na
data_preparation = function(data) {
  df <- na.omit(data)
  df <- ifelse(df < 3, 0, 1)
  return(df)
}

# Generates adjacency graph of a small world graph given size p
# Generates adjacency graph of a small world graph given size p
small_world = function(p) {
  
  sw <- igraph::sample_smallworld(1, p, 2, .2)
  check = qgraph::smallworldness(sw)
  check = check[1]
  pr = .3
  while ((check < 1) & (pr < 1)) {
    sw = igraph::sample_smallworld(1, p, 2, pr)
    check = qgraph::smallworldness(sw)
    check = check[1]
    pr = pr + .1
  }
  adj <- as.matrix(igraph::as_adjacency_matrix(sw))
  return(adj)
}
# Generate adjacency graph of a random graph given size p
random_graph = function(p) {
  
  rg <- igraph::erdos.renyi.game(p, .3)
  adj <- as.matrix(igraph::as_adjacency_matrix(rg))
  return(adj)
  
}

#Generate all parameters for all options of p for a complete graph
# using a given data set to estimate from
generate_all_parameters_complete = function(p_list, df) {
  
  params <- list()
  idx <- 1
  for (p in p_list) {
    df_prep <- df[, 1:p]
    res <- psychonetrics::Ising(df_prep)
    res <- res %>% runmodel
    res <- transform_psych_results(res, p)
    sigma <- res$sigma$est
    mu <- res$mu$est
    
    param_tmp <- list(mu = mu, sigma = sigma / 2)
    params[[idx]] <- param_tmp
    idx <- idx + 1
  }
  return(params)
}

#Generate all parameters for all options of p for a small world graph
# using a given data set to estimate from
generate_all_parameters_smallworld = function(p_list, df) {
  
  params <- list()
  idx <- 1
  for (p in p_list) {
    df_prep <- df[, 1:p]
    graph <- small_world(p)
    res <- psychonetrics::Ising(df_prep, omega = graph)
    res <- res %>% runmodel
    res <- transform_psych_results(res, p)
    
    params[[idx]] <- list(
      mu = res$mu$est, sigma = res$sigma$est / 2
    )
    idx <- idx + 1
  }
  return(params)
}

#Generate all parameters for all options of p for a random graph
# using a given data set to estimate from
generate_all_parameters_random = function(p_list, df) {
  
  params <- list()
  idx <- 1
  for (p in p_list) {
    df_prep <- df[, 1:p]
    graph <- random_graph(p)
    res <- psychonetrics::Ising(df_prep, omega = graph)
    res <- res %>% runmodel
    res <- transform_psych_results(res, p)
    
    params[[idx]] <- list(
      mu = res$mu$est, sigma = res$sigma$est / 2
    )
    idx <- idx + 1
  }
  return(params)
}

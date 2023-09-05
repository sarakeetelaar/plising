#####################################################################
# This file contains the functions used for the simulations, wrapping
# around the different methods
#####################################################################

# Function performs a single simulation run for one replications for given 
# parameters (hwich includes graph size and type), sample size and a seed
single_simulation <- function(params, N, seed) {
  set.seed(seed)
  x <- full_data_generation(params, N)
  mu <- params$mu
  sigma <- params$sigma
  p <- length(mu)
  
  results <- list()
  
  #pseudolikelihood
  pl <- pseudolikelihood(x)

  #exact likelihood using package psychonetrics
  ml <- maximum_likelihood(x)
  
  # disjoint pl (called lr)
  lr <- logistic_regression(x)
  return(list(pl = pl, ml = ml, lr = lr, sigma = sigma, mu = mu, N = N, p = p))
}

# Function performs exact maximum likelihood estimation using psychonetrics 
maximum_likelihood <- function(data) {
  
  out <- tryCatch({
    p <- ncol(data)
    model <- psychonetrics::Ising(data)
    mod <- model %>% runmodel
    ml <- transform_psych_results(mod, p)
    return(ml)
  }, error = function(cond) {
    message(cond)
    return(NA)
  })
}


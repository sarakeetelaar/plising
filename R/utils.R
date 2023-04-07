#helper function to get indices of mu
re_index = function(p=10) {
  start = 1
  indices = rep(0, p)
  count = 0
  while(count<(p-1)) {
    indices[count+2] = start + (p - count)
    start = start + (p-count)
    count = count+1
    
  }
  indices[1] = 1 
  return(indices)
}

# logistic function 
expit <- function(phi) {
  return(exp(phi) / (1 + exp(phi)))
}

nexpit <- function(phi) {
  return(exp(phi) / (1 + exp(phi)) ^ 2)
}

# indexing function used to convert the pxp
# parameter matrix sigma to a {p * (p + 1) / 2} vector
indexing <- function(p) {
  index <- matrix(0, nrow = p * (p - 1) / 2, ncol = 3)
  counter <- 0
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      counter <- counter + 1
      index[counter, 1] <- counter
      index[counter, 2] <- s
      index[counter, 3] <- t
    }
  }
  return(index)
}
transform_psych_results = function(res, p) {
  mlres = res@parameters[["est"]]
  mlses = res@parameters[["se"]]
  l_param = length(mlres)
  sig_index = (1:(l_param - 1))[-(1: p)]
  mu_index = 1: p
  
  muest = mlres[mu_index]
  muvar = mlses[mu_index] ^ 2
  
  sigest = mlres[sig_index]
  sigvar = mlses[sig_index] ^ 2
  
  
  sigest = fill_matrix(sigest, p)
  sigvar = fill_matrix(sigvar, p)
  
  return(list(mu=list(est=muest, var=muvar),
              sigma = list(est=sigest, var=sigvar)))
}

fill_matrix = function(sigvals, p) {
  sigma = matrix(0, p, p)
  sigma[lower.tri(sigma)] = sigvals
  sigma = sigma + t(sigma)
  
  return(sigma)
}
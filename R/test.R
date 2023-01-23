data = function(p=10, N=1) {
  n_param = p*(p+1)/2
  theta = matrix(0, p, p)
  params = runif(n_param, -.5, .5)
  theta[lower.tri(theta, diag=TRUE)] = params
  
  mu = diag(theta)
  sigma = theta
  diag(sigma) = 0
  
  for (i in 1:(p-1)) {
    for (j in (i+1): p) {
      theta[i,j] = theta[j,i]
    }
  }
  
  X = matrix(0, nrow=N, ncol=p)
  for (iter in 1:1e2) {
    for (i in 1:p) {
      X[,i] = 1 * (rlogis(N) <= mu[i] + 2 * X %*% sigma[, i])
    }
  }
  return(list(data=X, theta=params))
}

normalizing_constant = function(theta, x) {
  p = ncol(x)
  N = nrow(x)
  Z = 0
  full_y = expand.grid(rep(list(0:1), p))
  for (ind in 1:nrow(full_y)) {
    vals = unlist(full_y[ind,])
    y = matrix(vals, nrow=p, ncol=1)
    suf = 2 * y %*% t(y)
    diag(suf) = diag(suf)/2
    suf = suf[lower.tri(suf, diag=T)]
    
    weight = as.numeric(exp(t(suf) %*% theta))
    Z = Z + weight
  }  
  return(N * log(Z))    
}

normalizing_pseudo = function(theta, x) {
  p = ncol(x)
  N = nrow(x)
  theta_matrix = matrix(0, nrow=p,ncol=p)
  print(length(theta_matrix[lower.tri(theta_matrix, diag=TRUE)]))
  theta_matrix[lower.tri(theta_matrix, diag=T)] = theta
  mu = diag(theta_matrix)
  sigma = theta_matrix
  diag(sigma) = 0
  N = nrow(x)
  p = ncol(x)
  
  Z = 0 
  for (v in 1:N) {
    for (i in 1:p) {
      helpsum = x[v,-i] %*% sigma[-i, i]
      Z = Z + log(1 + exp(mu[i] + helpsum))
      
    }
  }
  return(Z)
}

set.seed(1)
# p_list = 2:40
# l = length(p_list)
# fractions = rep(0, l)
# counter = 1
# for (p in p_list) {
#   print(paste("p =", p))
#   data_generation = data(p=p)
#   x = data_generation$data
#   theta = data_generation$theta
#   Z = normalizing_constant(theta, x)
#   Z_pl = normalizing_pseudo(theta, x)
#   
#   
#   fractions[counter] = Z_pl/Z
#   counter = counter + 1
# }
# 
# plot(p_list, fractions, type="l")

n_list = seq(100, 10000, 100)
len = length(n_list)
fractions = rep(0, len)
counter = 1
for (i in 1:len) {
  print(counter)
  n=n_list[i]
  data_generation = data(p=3, N=n)
  x = data_generation$data
  theta = data_generation$theta
  Z = normalizing_constant(theta, x)
  Z_pl = normalizing_pseudo(theta, x)
  
  fractions[i] = Z_pl/Z
  counter = counter + 1
  
}
plot(n_list, fractions, type="l", xlab="N", ylab="fraction", main="Fraction of normalizing constants over N")



expprob2 = function(x, sigma, mu, v_ind, q_ind) {
  sumsigma = 0
  p = ncol(x)
  for (q in 1:p) {
    if (q != q_ind) {
      sumsigma = sumsigma + sigma[q_ind, q] * x[v_ind, q]
    }
  }
  pr = exp(mu[q_ind] + sumsigma) / (1 + exp(mu[q_ind] + sumsigma))
  return(pr)
}

hessian_pl = function(x, sigma, mu) {
  p = length(mu)
  n = nrow(x)
  nsigma = p * (p - 1) / 2
  nparam = p* (p + 1) / 2
  
  hess = matrix(0, nparam, nparam)
  mu_hessian= muhessian(x, sigma, mu)
  sig_hessian = sigmahessian(x, sigma, mu)
  cross_hessian = crosshessian(x, sigma, mu)
  
  hess[1:p, 1:p] = mu_hessian
  hess[1:p, (p+1): nparam] = cross_hessian
  hess[(p + 1): nparam, 1:p] = t(cross_hessian)
  hess[(p + 1): nparam, (p + 1): nparam] = sig_hessian
  

  pq = matrix(0, nrow = n, ncol = p)
  for(s in 1:p) {
    phi = sigma[s, s] + x[, -s] %*% sigma[-s, s]
    pq[, s] = nexpit(phi)
  }
  
  inv_hes = matrix(0, nrow = p * (p + 1) / 2, ncol = p *(p + 1) / 2)
  index_mu = 1:p
  index_sig = (1:(p * (p + 1) / 2))[-index_mu]
  
  inv_mu = matrix(0, p, p)
  for (s in 1:p) {
    inv_mu[s, s] = 1 / (-sum(pq[, s]))
  }
  
  inv_hes[index_sig, index_sig] = sig_hessian - t(cross_hessian) %*% inv_mu %*% cross_hessian
  inv_hes[index_sig, index_sig] = solve(inv_hes[index_sig, index_sig])
  inv_hes[index_mu, index_sig] = - inv_mu %*% cross_hessian %*% inv_hes[index_sig, index_sig]
  inv_hes[index_sig, index_mu] = t(inv_hes[index_mu, index_sig])
  inv_hes[index_mu, index_mu] = inv_mu - inv_hes[index_mu, index_sig] %*% t(cross_hessian) %*% inv_mu
  
  return(list(hessian=hess, inverse=inv_hes))
  
  
}

variance_pl = function(x, sigma, mu) {
  p = length(mu)
  hessian = hessian_pl(x, sigma, mu)
  inv_hessian  = hessian$inverse
  
  outergrad = outerGradient(x, sigma, mu)
  print(outergrad)
  print(inv_hessian)
  inv2 = solve(-hessian$hessian)
  vars = (inv2 %*% outergrad) %*% inv2
  
  variances = diag(vars)
  var_mu = variances[1:p]
  var_sig = variances[-(1:p)]
  
  return(list(mu=var_mu, sigma=var_sig))
  
}

muhessian = function(x, sigma, mu) {
  p = ncol(x)
  N = nrow(x)
  
  muhes = matrix(0, p, p)
  for (i in 1:p) {
    diag = 0
    for (v in 1:N) {
      sigsum = 0
      for (q in 1:p) {
        if (q != i) {
          sigsum = sigsum + sigma[i, q] * x[v, q]
        }
      }
      thispr = exp(mu[i] + sigsum)/ (1 + exp(mu[i] + sigsum))
      diag = diag + thispr * (1 - thispr)
    }
    muhes[i, i] = diag
  }
  return(muhes)
}

sigmahessian = function(x, sigma, mu) {
  p = ncol(x)
  N = nrow(x)
  nTheta = p * (p - 1) / 2
  sighes = matrix(0, nTheta, nTheta)
  uppertri = matrix(0, nTheta * (nTheta - 1) / 2, 1)
  index = 1
  
  for (i1 in 1:(p-1)) {
    for (j1 in (i1 + 1): p) {
      k = j1
      change = T
      for (r in i1:(p-1)) {
        if (change == F) {
          k = (r + 1)
        }
        while (k <= p) {

          el = 0
          if ((i1 == r) & (j1 == k)) {
            for (v in 1:N) {
              pi = expprob2(x, sigma, mu, v, i1)
              pj = expprob2(x, sigma, mu, v, j1)
              term = - x[v, j1] * pi * (1 - pi) - x[v, i1] * pj * (1 - pj)
              
              el = el + term
            }
          }
          else if ((i1 == r) & (j1 != k)) {
            for (v in 1:N) {
              pi = expprob2(x, sigma, mu, v, i1)
              term = -x[v, j1] * x[v, k] * pi * (1 - pi)
              
              el = el + term
            }
          }
          else if ((i1 != k) & (j1 == r)) {
            for (v in 1:N) {
              pj = expprob2(x, sigma, mu, v, j1)
              term = -x[v, i1] * x[v, k] * pj * (1 - pj)
              
              el = el + term
            }
          }
          uppertri[index] = el
          index = index + 1
          k = k + 1
          change = F
        }
      }
    }
  }
  print(length(uppertri))
  sighes[upper.tri(sighes, diag=T)] = uppertri
  
  for (i in 1:(nTheta-1)) {
    for (j in (i+1):nTheta) {
      sighes[j, i] = sighes[i, j]
    }
  }
  return(sighes)
}

crosshessian = function(x, sigma, mu) {
  p = ncol(x)
  N = nrow(x)
  nTheta = p * (p - 1) / 2
  count = 1
  crosshes = matrix(0, nrow=p, ncol=nTheta)
  
  for (mu_ind in 1:p) {
    for (sig_1 in 1:(p-1)) {
      for (sig_2 in (sig_1 + 1): p) {
        el = 0
        if (sig_1 == mu_ind) {
          for (v in 1:N) {
            pi = expprob2(x, sigma, mu, v, mu_ind)
            term = -x[v, sig_2] * pi * (1 - pi)
            el = el + term
          }
        }
        else if (sig_2 == mu_ind) {
          for (v in 1:N) {
            pi = expprob2(x, sigma, mu, v, mu_ind)
            term = -x[v, sig_1] * pi * (1 - pi)
            
            el = el + term
          }
        }
        rowNum = floor(count / nTheta)
        colNum = count %% nTheta
        count = count + 1
        
        crosshes[rowNum, colNum] = el
      }
    }
  }
  return(crosshes)
}















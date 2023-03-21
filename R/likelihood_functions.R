invert_hessian_ml = function(hessian, p) {
  inv_hessian = matrix(0, nrow = p * (p + 1) / 2, ncol = p * (p + 1) /2)
  index_mu = re_index(p)
  index_sig = (1:(p* (p + 1)/2))[-index_mu]
  
  muhes = hessian[index_mu, index:mu]
  sighes = hessian[index_sig, index_sig]
  crosshes = hessian[index_mu, index_sig]
}
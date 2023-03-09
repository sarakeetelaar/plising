hessian_pl = function(x, sigma, mu) {
  p = length(mu)
  nsigma = p * (p - 1) / 2
  nparam = p* (p + 1) / 2
  
  hess = matrix(0, nparam, nparam)
  mu_hessian = muHessian(x, sigma, mu)
  sig_hessian = sigmaHessian(x, sigma, mu)
  cross_hessian = crossHessian(x, sigma, mu)
  
  hess[1:p, 1:p] = mu_hessian
  hess[1:p, (p+1): nparam] = cross_hessian
  hess[(p + 1): nparam, 1:p] = t(cross_hessian)
  hess[(p + 1): nparam, (p + 1): nparam] = sig_hessian
  
  inv_hes = solve(hess)
  
  return(list(hessian=hess, inverse_hessian=inv_hes))
  
  
}
# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

expprob <- function(x, sigma, mu, v_ind, q_ind) {
    .Call(`_plising_expprob`, x, sigma, mu, v_ind, q_ind)
}

muHessian <- function(x, sigma, mu) {
    .Call(`_plising_muHessian`, x, sigma, mu)
}

sigmaHessian <- function(x, sigma, mu) {
    .Call(`_plising_sigmaHessian`, x, sigma, mu)
}

crossHessian <- function(x, sigma, mu) {
    .Call(`_plising_crossHessian`, x, sigma, mu)
}

singleGradientPL <- function(x, sigma, mu, obs) {
    .Call(`_plising_singleGradientPL`, x, sigma, mu, obs)
}

outerGradient <- function(x, sigma, mu) {
    .Call(`_plising_outerGradient`, x, sigma, mu)
}

normalizingConstant <- function(theta, Esuf, Ess, Z, y) {
    .Call(`_plising_normalizingConstant`, theta, Esuf, Ess, Z, y)
}

symmetrizeMatrix <- function(theta) {
    .Call(`_plising_symmetrizeMatrix`, theta)
}

hessenNorm <- function(thetaH, EsufH, EssH, ZH, yH) {
    .Call(`_plising_hessenNorm`, thetaH, EsufH, EssH, ZH, yH)
}

sumSigma <- function(sigma, x, index, v) {
    .Call(`_plising_sumSigma`, sigma, x, index, v)
}

derivativeHelp <- function(x, mu, sigma) {
    .Call(`_plising_derivativeHelp`, x, mu, sigma)
}


#' Optimize sGen
#'
#' This is a component function within sGenSolver that estimates the log-
#' likelihood of the relevant Ricker function.
#'
#' @param S A numeric of spawner abundances.
#' @param theta A numeric vector of Ricker stock recruit parameters: alpha,
#' beta, and sigma.
#' @param sMSY A numeric of spawner abundance estimated to result in maximum
#' sustainable yield.
#' @return Returns a list of prt (estimated productivity), epsilon (the
#' difference between sMSY and productivity), the negative log likelihood
#' used to estimate sGen, and S.
#'
#' @export

sGenOptimum <- function(S, theta, sMSY) {
  a = theta[1]
  b = theta[2]
  sig = exp(theta[3])
  prt <- S * exp(a - b * S)
  epsilon <- log(sMSY) - log(prt)
  nLogLike <- sum(dnorm(epsilon, 0, sig, log = T))

  return(list(prt = prt, epsilon = epsilon, nLogLike = nLogLike, S = S))
}

#______________________________________________________________________________

#' Solve for sGen
#'
#' This function solves for sGen based on sMSY and the log-likelihood estimated
#' in sGenOptimum.
#'
#' @param theta A numeric vector of Ricker stock recruit parameters: alpha,
#' beta, and sigma.
#' @param sMSY A numeric of spawner abundance estimated to result in maximum
#' sustainable yield.
#' @return Returns a numeric that is the spawner abundance that minimizes the
#' the log likelihood.
#'
#' @examples
#' Stock-recruit parameters approximate those of Fraser River sockeye salmon
#' Chilko CU.
#' alpha = 1.2
#' beta = 1.5
#' sigma = 0.8
#' theta = c(alpha, beta, sigma)
#' sMSY = 0.3
#' sGenSolver(theta, sMSY)
#' @export

sGenSolver <- function(theta, sMSY) {
  #gives the min Ricker log-likelihood
  fnSGen <- function(S, theta, sMSY) -1.0 * sGenOptimum(S, theta, sMSY)$nLogLike

  fit <- optimize(f = fnSGen, interval = c(0, ((theta[1] / theta[2]) * (0.5 - 0.07 * theta[1]))),
                 theta = theta, sMSY = sMSY)
  return(list(fit = fit$minimum))
}

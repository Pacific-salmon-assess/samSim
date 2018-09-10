#' Optimize Sgen
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

  return(list(prt = prt, epsilon = epsilon, nloglike = nloglike, S = S))
}

#______________________________________________________________________________

#' Function Sgen
#'
#' This function solves for sGen based on sMSY and the log-likelihood estimated
#' in sGenOptimum.
#'
#' @param S A numeric of spawner abundances.
#' @param theta A numeric vector of Ricker stock recruit parameters: alpha,
#' beta, and sigma.
#' @param sMSY A numeric of spawner abundance estimated to result in maximum
#' sustainable yield.
#' @return Returns a numeric that is the value minimizing the log likelihood.
#'
#' @export

sGenSolver <- function(S, theta, sMSY) {
  #gives the min Ricker log-likelihood
  fnSGen <- function(S, theta, sMSY) -1.0 * Sgen.optim(S, theta, sMSY)$nloglike

  fit <- optimize(f = fnSGen, interval = c(0, ((theta[1] / theta[2]) * (0.5 - 0.07 * theta[1]))),
                 theta = theta, sMSY = sMSY)
  return(list(fit = fit$minimum))
}



## NOTE SMSY estimate needs to be updated to reflect change to Scheuerell 2016 calc
fn.sgen <- function(s, theta, s.msy) -1.0*Sgen.optim(s, theta, s.msy)$nloglike	#gives the min Ricker LL, dd is a dummy
solver.sgen <- function(s, theta, s.msy)
{
  #fit=optimize(f=fn.sgen,interval=c(0,S.msy[m]), theta=theta, s.msy=s.msy)
  fit = optimize(f = fn.sgen,interval = c(0, ((theta[1] / theta[2]) * (0.5 - 0.07 * theta[1]))),
                 theta = theta, s.msy = s.msy)
  #if(fit$convergence!=0)print("Non-convergence for benchmark S.lower.3")
  return(list(fit = fit$minimum))
}


solver.sgen(theta = c(1.2, 0.5,
                          0.7),
            s.msy = 0.3)

Sgen.optim(theta = c(1.2, 0.5,
                     0.7),
           s.msy = 0.3)

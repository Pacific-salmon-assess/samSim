#' Generate log-normal error associated with proportions data
#'
#' This function generates proporations at age with multivariate logistic error (Schnute and
#' Richards 1995, eqns.S.9 and S.10). Includes an adjustment to skip calculations when dealing
#' with populations that have fixed age structures (e.g. pink salmon).
#'
#' @param ppnAgeVec A numeric vector of proportions, generally individuals returning at a given
#'  age.
#' @param nAges An integer specifying number of potential ages.
#' @param tau A numeric specifying the parameter representing average variation in proportions (e.g.
#'  observed interannual variation in age stracuture).
#' @param error A numeric vector of length \code{nAges} that specifies random values to deviate proportions
#'  by (generally generated from a uniform distribution bounded by 0 and 1).
#' @return Returns a numeric vector of length \code{nAges} representing proportions for each class.
#'
#' @examples
#' ppnAgeVec <- c(0.2, 0.4, 0.3, 0.1)
#' nAges <- length(ppnAgeVec)
#' tau <- 0.7
#' error <- runif(nAges, 0.0001, 0.9999)
#' ppnAgeErr(ppnAgeVec, nAges, tau, error)
#'
#' @export

ppnAgeErr <- function(ppnAgeVec, nAges, tau, error) {
  ppnAgeVec[is.na(ppnAgeVec)] <- 0 # NAs produced when dividing 0 by 0 for proportion of recruits at age; replace w/ 0s
  p <- rep(0, length.out = nAges)
  if (max(ppnAgeVec) == 1) { #if any of the average proportions equals one, assume the structure is fixed and adjust output
    pos <- which(ppnAgeVec == 1)
    p[pos] <- 1
  } else {
    ppnAgeVec[ppnAgeVec == 0] <- 0.0000001 #change 0s to very small values to calc ppn
    dum <- 0.
    x <- 0.
    dum <- log(ppnAgeVec) + tau * qnorm(error, 0, 1) #vectorized; no loops needed
    x <- log(ppnAgeVec)+ tau * qnorm(error, 0, 1) - (1 / nAges) * sum(dum)
    p <- exp(x) / sum(exp(x))
  }
  return(p)
}

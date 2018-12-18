#' Calculate multivariate logistric proportions
#'
#' This is a component function within getTau used to calculate log proportions.
#' Originally written by B. Dorner and C. Holt so C. Freshwater unable to
#' fully verify documentation.
#'
#' @param tau Value for tau to test
#' @param logProps Log-transformed vector of mean proportion data
#'
#' @return Vector of simulated proportion data based on specified tau and means.
#' @export
#'
#' @examples
#'
mvLogisticLogProp <- function(tau, logProps){
  p <- logProps[!is.na(logProps)]
  n <- length(p)
  rNum <- rnorm(n, sd=1)
  eps <- tau*rNum
  x <- p + eps - (1/n) * sum(p + eps)
  result <- rep(NA, n)
  result[!is.na(logProps)] <- exp(x)/sum(exp(x))
  result
}

#______________________________________________________________________________

#' Estimate tau
#'
#' Main function used to estimate variation in multivariate proportion data.
#' Uses an objective function to minimize difference between sum of squares of
#' observed vs. simulated standard deviations for proportions data. Note that
#' objective function is constrained to sample between 0 and 3 at intervals of
#' 0.1; however that should be sufficient for most proportion data.
#'
#' @param ppnMat Matrix of proportion data to apply grid search for tau to.
#' @param plotTaus If true (default), tau-specific values from objective
#' function are plotted allowing users to visually identify minimum.
#'
#' @return A list containing a vector of tau-specific objective function values
#' (*objFun*) and a numeric representing the tau value that minimizes the
#' objective function (*bestTau*).
#' @export
#'
#' @examples
#'
getTau <- function(ppnMat, plotTaus = TRUE){
  targetFun <- function(tau, expLogProp, targetSD, n=10000) {
    cat(". ")
    Tau <- rep(tau, n)
    simProp <- t(sapply(Tau, mvLogisticLogProp, logProp=expLogProp))
    simSD <- apply(simProp, 2, sd)
    return(sum((simSD - targetSD)^2, na.rm=TRUE))
  }#end targetFun

  if (any(rowSums(ppnMat) > 1.01)) {
    warning("Some proportions do not sum to one; check calculations")
  }

  # Replace 0 values that crash grid search
  ppnMat[ppnMat == 0] <- 1e-10

  ## Grid search for best tau:
  tau <- seq(from = 0, to = 3, by = 0.1)
  result <- data.frame(tau=tau, objective=rep(NA, length(tau)))
  cat("finding scale parameter for variability in relative recruit
      proportions:\n")
  for(t in tau){
    cat(". ")
    result[t*10+1, "objective"] <- targetFun(t, apply(log(ppnMat), 2, mean),
                                             apply(ppnMat, 2, sd))
  }
  cat("\n")
  bestTau <- tau[result[, "objective"] == min(result[, "objective"])]
  if (plotTaus == TRUE) {
    plot(tau, result$objective, xlab="tau", ylab="objective function")
  }
  return(list(objFun=result, bestTau=bestTau))
}#end getTau

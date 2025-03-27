#' Calculate multivariate logistic proportions
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
#' 
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
#'
#' ppnMat<-matrix(0,nrow=40,ncol=4)
#' for(i in 1:40){
#'   ppnMat[i,]<-ppnAgeErr(c(0.1,0.3,0.4,0.2), 0.25,
#'        error = runif(4, 0.0001, 0.9999)) 
#' }
#'
#' #estimate and recover tau
#' getTau(ppnMat, plotTaus = TRUE)
#'
#'
#'
#'
getTau <- function(ppnMat, gridstep = 0.01, plotTaus = TRUE){
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
  tau <- seq(from = 0, to = 3, by = gridstep)
  result <- data.frame(tau=tau, objective=rep(NA, length(tau)))
  cat("finding scale parameter for variability in relative recruit
      proportions:\n")
  for(ti in seq_along(tau)){
    cat(". ")
    result[ti, "objective"] <- targetFun(tau[ti], apply(log(ppnMat), 2, mean),
                                             apply(ppnMat, 2, sd))
  }
  cat("\n")
  bestTau <- tau[result[, "objective"] == min(result[, "objective"])]
  if (plotTaus == TRUE) {
    plot(tau, result$objective, xlab="tau", ylab="objective function")
  }
  return(list(objFun=result, bestTau=bestTau))
}#end getTau

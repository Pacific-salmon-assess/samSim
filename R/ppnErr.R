#' Generate log-normal error associated with proportions data
#'
#' This function generates proporations at age with multivariate logistic error
#' (Schnute and Richards 1995, eqns.S.9 and S.10). Includes an adjustment to
#' skip calculations when dealing with populations that have fixed age struc-
#' tures (e.g. pink salmon).
#'
#' @param ppnAgeVec A numeric vector of proportions of individuals returning
#' at a given age. \code{nAges = length(ppnAgeVec)}.
#' @param tau A numeric specifying the parameter representing average variation
#' in proportions (e.g. observed interannual variation in age stracuture).
#' @param error A numeric vector of length \code{nAges} that specifies random
#' values to deviate proportions by (generally generated from a uniform dis-
#' tribution bounded by 0 and 1).
#' @return Returns a numeric vector of length \code{nAges} representing pro-
#' portions for each class.
#'
#' @examples
#' ppnAgeVec <- c(0.2, 0.4, 0.3, 0.1)
#' nAges <- length(ppnAgeVec)
#' tau <- 0.7
#' error <- runif(nAges, 0.0001, 0.9999)
#' ppnAgeErr(ppnAgeVec, tau, error)
#'
#' @export

ppnAgeErr <- function(ppnAgeVec, tau, error) {
  nAges <- length(ppnAgeVec)
  #NAs produced when dividing 0 by 0 for ppn of recruits at age; replace w/ 0s
  ppnAgeVec[is.na(ppnAgeVec)] <- 0
  p <- rep(0, length.out = nAges)
  #if any of the average proportions equals one, assume the structure is fixed
  #and adjust output accordingly
  if (max(ppnAgeVec) == 1) {
    pos <- which(ppnAgeVec == 1)
    p[pos] <- 1
  } else {
    ppnAgeVec[ppnAgeVec == 0] <- 0.0000001
    dum <- 0.
    x <- 0.
    dum <- log(ppnAgeVec) + tau * qnorm(error, 0, 1) #vectorized
    x <- log(ppnAgeVec)+ tau * qnorm(error, 0, 1) - (1 / nAges) * sum(dum)
    p <- exp(x) / sum(exp(x))
  }
  return(p)
}

#______________________________________________________________________________

#' Generate log-normal error associated with catch data
#'
#' Like the \code{ppnAgeErr} function, this function generates assigns multi-
#' variate logistic error (Schnute and Richards 1995, eqns.S.9 and S.10) to
#' proportions data. Here it is intended to represent error in assigning mixed
#' stock catches to specific stocks.
#'
#' @param ppnCatchVec A numeric vector of proportions representing the true
#' relative contribution of each stock to total catch. \code{nStocks =
#' length(ppnAgeVec)}.
#' @param tau A numeric specifying the parameter that represents average var-
#' iation in proportions.
#' @return Returns a numeric vector representing proportions for each stock.
#'
#' @examples
#' ppnCatchVec <- c(0.2, 0.4, 0.3, 0.1)
#' tau <- 0.7
#' ppnCatchErr(ppnCatchVec, tau)
#'
#' @export

ppnCatchErr <- function(ppnCatchVec, tau){
  ppnCatchVec[is.na(ppnCatchVec)] <- 0
  nStocks <- length(ppnCatchVec)
  #unlike ppnAgeErr, error can be estimated internally because don't need to
  #align multiple years
  error <- runif(nStocks, 0.0001, 0.9999)
  p<-0
  ppnCatchVec[ppnCatchVec == 0] <- 0.0000001
  dum<-0.
  x<-0.
  dum <- log(ppnCatchVec) + tau * qnorm(error, 0, 1)
  x <- log(ppnCatchVec) + tau * qnorm(error, 0, 1) - (1/nStocks) * sum(dum)
  p <- exp(x) / sum(exp(x))
  return(p)
}

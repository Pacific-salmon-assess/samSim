#' Calculate realized catches
#'
#' This function calculates realized catch based on beta distributed outcome
#' uncertainty.
#'
#' Outcome uncertainty is often easier to parameterize using observed deviations
#' in target and realized exploitation rates rather than catches. Here target
#' total allowable catches are back-converted to exploitation rates and realized
#' catches are generated using  beta distributed outcome uncertainty (as in
#' Anderson et al. 2015 J. App. Ecol.).
#'
#' For simplicity's sake when TAC is zero, realized HR is zero. Therefore to
#' represent "significant" catch even when TAC is zero, pass large sigma values
#' and a small, non-zero TAC value.
#'
#' **Note** the sigma values representing the standard deviation of the shape
#' parameter for the beta distribution are **not** equivalent to the SD of a
#' normal distribution. Values larger than ~0.2 can result in U-shaped realized
#' exploitation rates. Distribution should be parameterized based on data or use
#' default value (0.1) from Pestes et al. 2008 for Cultus Lake sockeye salmon.
#'
#' **Note** small target harvest rates and/or high sigmas can produce negative
#' location parameters which result in NAs; replace with small values
#' (see betaDistributionBounds.Rmd for details) .
#'
#' @param rec A numeric vector of length nCU representing CU-specific
#' recruitment.
#' @param tac A numeric vector of length nCU representing CU-specific target
#' total allowable catch rates (generally passed from \code{calcTAC} function).
#' @param sigma A numeric representing the standard deviation of the shape
#' parameter.
#' @return Returns a vector of length nCU representing realized catches
#'
#' @examples
#' head(exampleHCRList)
#' exRec <- exampleHCRList$recRY
#' exTAC <- 0.2 * exRec
#' calcRealCatch(exRec, exTAC, sigma = 0.1)
#'
#' @export
calcRealCatch <- function(rec, tac, sigma = 0.1) {
  #hack to replace 0s and let function run without looping
  tempRec <- rec
  tempRec[tempRec == 0] <- 0.00001
  if (sigma == 0) {
    mu <- tac / tempRec
    realCatch <- mu * rec
  }

  if (sigma > 0) {
    ## Run as loop if 0s are present, otherwise vectorize
    if (any(tac == 0)) {
      nCU <- length(rec)
      realCatch <- rep(NA, length.out = nCU)
      for (k in seq_along(rec)) {
        #calc target ER (capped at 99% which may occur if other fisheries
        #overfish due to OU)
        mu <- min(0.99, tac[k] / tempRec[k])
        if (mu != 0) {
          location <- pmax(0.0001,
                           mu^2 * (((1 - mu) / sigma^2) - (1 / mu)))
          shape <- location * (1 / mu - 1)
          realER <- rbeta(1, location, shape, ncp = 0)
          realCatch[k] <- realER * rec[k]
        } else {
          #if TAC is 0 then realCatch is 0
          realCatch[k] <- 0
        }
      } #end CU loop
    } else {
      #calc target ER (capped at 99% which may occur if other fisheries
      #overfish due to OU)
      mu <- pmin(0.99, tac / tempRec)
      location <- pmax(0.0001,
                       mu^2 * (((1 - mu) / sigma^2) - (1 / mu)))
      shape <- location * (1 / mu - 1)
      realER <- rbeta(length(mu), location, shape, ncp = 0)
      realCatch <- realER * rec
    } #end vectorized version
  }

  #warning message for nonsensical values
  if (any(is.na(realCatch))) {
    stop("Realized ER can't be calculated; check that target HR and sigma
               are reasonable")
  }

  realCatch
}

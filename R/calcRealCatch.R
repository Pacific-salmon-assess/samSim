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
#' @param random A logical (default `FALSE`) used to restandardize random number
#' generator because `rbeta()` seems to produce variable draws (should be
#' corrected in future).
#' @param setSeedInput A numeric representing the combination MC trial and
#' simulation year to ensure that scenarios are sampling the seed, but allowing
#' draws to otherwise vary.
#' @return Returns a vector of length nCU representing realized catches
#'
#' @examples
#' head(exampleHCRList)
#' exRec <- exampleHCRList$recRY
#' exTAC <- 0.2 * exRec
#' calcRealCatch(exRec, exTAC, sigma = 0.1)
#'
#' @export
calcRealCatch <- function(rec, tac, sigma = 0.1, random =  FALSE,
                          setSeedInput = NULL) {
  #hack to replace 0s and let function run without looping
  tempRec <- rec
  tempRec[tempRec == 0] <- 0.00001
  if (sigma == 0) {
    mu <- tac / tempRec
    realCatch <- mu * rec
  }

  erMat <- matrix(NA, nrow = length(rec), ncol = 4)

  if (sigma > 0) {
    mu <- pmin(0.99, tac / tempRec)
    #Constraint location and shape parameters to be non-negative, otherwise NaN
    #produced
    location <- pmax(0.00001,
                     mu^2 * (((1 - mu) / sigma^2) - (1 / mu)))
    shape <- pmax(0.00001, location * (1 / mu - 1))
    ## Run as loop if 0s are present, otherwise vectorize
    if (any(tac == 0)) {
      nCU <- length(rec)
      realCatch <- rep(NA, length.out = nCU)
      for (k in seq_along(rec)) {
        #calc target ER (capped at 99% which may occur if other fisheries
        #overfish due to OU)
        if (mu[k] != 0) {
          realER <- rbeta(n = 1, shape1 = location[k], shape2 = shape[k])
          realCatch[k] <- realER * rec[k]
        } else {
          #if TAC is 0 then realCatch is 0
          realCatch[k] <- 0
          #make blank draw with dummy pars to balance random number generator
          #redundant giving final number reset but keep for now
          blank <- rbeta(1, 0.5, 0.5, ncp = 0)
        }
      } #end CU loop
    } else {
      #calc target ER (capped at 99% which may occur if other fisheries
      #overfish due to OU)
      realER <- rbeta(length(mu), location, shape)
      realCatch <- realER * rec
    } #end vectorized version
  }
  #get rid of absurdly small catches
  realCatch <- ifelse(realCatch < 1e-6, 0, realCatch)

  #warning message for nonsensical values
  if (any(is.na(realCatch))) {
    stop("Realized ER can't be calculated; check that target HR and sigma
               are reasonable")
  }

  #for whatever reason rbeta appears to draw a different number of randoms
  #depending on the shape number; until I can figure out why reset unless
  #running random chains intentionally
  if (random != TRUE) {
    set.seed(setSeedInput)
  }

  realCatch
}

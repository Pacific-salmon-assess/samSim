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
#' **Note** the sigma values representing the standard deviation of the shape
#' parameter for the beta distribution are **not** equivalent to the SD of a
#' normal distribution. Values larger than ~0.2 can result in U-shaped realized
#' exploitation rates. Distribution should be parameterized based on data or use
#' default value (0.1) from Pestes et al. 2008 for Cultus Lake sockeye salmon.
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
#' dum <- calcRealCatch(exRec, exTAC, sigma = 0.1)
#'
#' hist(dum$realER)
#' hist(dum$realCatch)
#'
#' @export
calcRealCatch <- function(rec, tac, sigma = 0.1) {
  #hack to replace 0s and let function run without looping
  tempRec <- rec
  tempRec[tempRec == 0] <- 0.00001
  #calc target ER
  mu <- pmax(0.00001, tac / tempRec)

  #calc realized harvest rates and catches
  if (sigma != 0) {
    #pmax necessary to prevent crashes with very small target TAC
    location <- pmax(0.1,
                     mu^2 * (((1 - mu) / sigma^2) - (1 / mu)))
    shape <- pmax(0.00001,
                  location * (1 / mu - 1))
    realER <- rbeta(length(mu), location, shape, ncp = 0)
    realCatch <- realER * rec
  }
  if (sigma == 0) {
    realCatch <- mu * rec
  }
  realCatch
}

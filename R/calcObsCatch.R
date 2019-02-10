#' Calculate observed catch
#'
#' This function adjusts realized catches by normally distributed observation
#' error to simulate observed catches. For mixed stock fisheries a tau
#' parameter representing multivariate logistic error is also incorporated to
#' account for uncertainty in the stock assignment process.
#'
#' @importFrom dplyr select group_by summarise_all summarise
#' @importFrom magrittr %>%
#'
#' @param catchVec A numeric vector representing realized CU-specific catches
#' (i.e. outcome uncertainty incorporated).
#' @param recVec A numeric vector representing CU-specific recruits (RY).
#' @param manUnit A character vector identifying the MU that each CU belongs to.
#' @param tauCatch A numeric representing multivariate logistic error associated
#' with assigning catch in mixed stock fisheries to the correct CU.
#' @param stkID A character vector of CU names.
#' @param catchObsErr A numeric representing log-nomrally distributed error in
#' catch observations.
#' @param extinctThresh A numeric representing the extinction threshold for the
#' aggregate.
#' @return Returns a numeric vector of observed catch.
#'
#' @examples
#' #Note that the function is intended to receive vectors rather than the DFs
#' #used in this example to increase efficiency.
#' head(exampleHCRList)
#'
#' catch <- exampleHCRList$mixCatch
#' rec <- exampleHCRList$recRY
#' manUnit <- exampleHCRList$mu
#' stock <- exampleHCRList$stock
#' calcObsCatch(catch, rec, manUnit, tauCatch = 0.1, stock, catchObsErr = 0.2,
#' extinctThresh = 0.0001)
#' @export
#'
calcObsCatch <- function(catchVec, recVec, manUnit, tauCatch, stkID, catchObsErr,
                         extinctThresh) {
  d <- data.frame(stkID = stkID,
                  mu = manUnit,
                  catch = catchVec,
                  rec = recVec,
                  obsErr = catchObsErr)
  d2 <- d %>%
    group_by(mu) %>%
    summarise(catchMU = sum(catch),
              recMU = sum(rec),
              n = length(rec))
  d3 <- merge(d, d2, by = "mu") %>%
    mutate(ppn = catch/catchMU,
           obsCatchOut = NA)
  muIndex <- unique(d2$mu)
  for (j in seq_along(muIndex)) {
    d4 <- subset(d3, mu == muIndex[j])
    ppnErr <- ppnCatchErr(d4$ppn, tauCatch)
    d3$obsCatchOut[which(d3$mu %in% d4$mu)] <- ppnErr * d4$catchMU * d4$obsErr
  }
  d3$obsCatchOut[which(d3$obsCatchOut < 1e-6)] <- 1e-6
  d3 <- with(d3, d3[order(stkID), ]) #reorder so output same as input
  return(d3$obsCatchOut)
}

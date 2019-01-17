#' Constrain mixed-stock harvest
#'
#' This function asseses whether mixed stock fisheries should be constrained
#' based on the abundance of co-migrating stocks. This process is modeled on
#' the current total allowable mortality rule used to manage Fraser River
#' sockeye salmon fisheries and is not a universal function.
#'
#' Total allowable catch (TAC) is set as a function of estimated abundance
#' at the management unit level. A given MU's TAC is constrained unless all
#' MUs with adjacent migration phenologies are above their upper fishery
#' reference point. Note that estimates of TAC and the application of
#' any constraints occurs in \code{calcTAC} function.
#'
#' **Note** Although all input vectors have length equal to \code{nCU}, the
#' number of unique values will be equal to \code{nMU} because this is the
#' scale at which in-season forecasts are available.
#'
#' **Note** Originally formatted with forecasted abundance but replaced with
#' true to avoid inflating outcome uncertainty.
#'
#' @param rec A numeric vector representing MU-specific return abundance.
#' @param highFRP A numeric vector representing MU-specific upper fishery
#' reference points. When return abundance is above this value
#' _adjacent_ MUs do not need to be constrained.
#' @param manAdjustment A numeric vector representing MU-specific management
#' adjustments. These values are used to adjust spawner abundance
#' to account for en route mortality (i.e. they increase the target escapement
#' goal).
#' @param manUnit A character vector identifying the MU that each CU belongs to.
#' @return Returns a two-element list of binary vectors. In \code{muAboveFRP}
#' ones represent MUs with return abundance above their FRP after incorporating
#' management adjustments. In \code{muConstrained} ones represent MUs that
#' should be constrained based on the abundance of _adjacent_ MUs.
#'
#' @examples
#' #Note that the function is intended to receive vectors rather than the DFs
#' #used in this example to increase efficiency.
#' head(exampleHCRList)
#' names(exampleHCRList)[4] <- "recRYMU"
#'
#' rec <- exampleHCRList$recRYMU
#' highFRP <- exampleHCRList$highFRP
#' manAdjustment <- exampleHCRList$adjustment
#' manUnit <- exampleHCRList$mu
#' constrain(rec, highFRP, manAdjustment, manUnit)
#' @export
constrain <- function(rec, highFRP, manAdjustment, manUnit) {
  muName <- unique(manUnit)
  nCU <- length(rec)
  nMU <- length(muName)
  muAboveFRP <- rep(0, nCU)
  conFinal <- rep(0, nCU)
  # Check 1: what is recruitment relative to reference point after adjusting
  # downwards w/ pMA
  for (k in 1:nCU) {
    if (rec[k] > highFRP[k] * (1 + manAdjustment[k])) {
      muAboveFRP[k] <- 1
    }
  }
  # Split into single MU value
  eStuAbove <- unique(muAboveFRP[which(manUnit %in% "EStu")])
  eSummAbove <- unique(muAboveFRP[which(manUnit %in% "ESumm")])
  summAbove <- unique(muAboveFRP[which(manUnit %in% "Summ")])
  lateAbove <- unique(muAboveFRP[which(manUnit %in% "Lat")])

  # Check 2: should each MU be constrained based on neighboring MUs status
  for (m in 1:nMU) {
    if (muName[m] == "EStu") {
      conFinal[which(manUnit %in% muName[m])] <- ifelse(eSummAbove == 1, 0, 1)
    }
    if (muName[m] == "ESumm") {
      conFinal[which(manUnit %in% muName[m])] <- ifelse(eStuAbove == 1 &
                                                          summAbove == 1, 0, 1)
    }
    if (muName[m] == "Summ") {
      conFinal[which(manUnit %in% muName[m])] <- ifelse(eSummAbove == 1 &
                                                          lateAbove == 1, 0, 1)
    }
    if (muName[m] == "Lat") {
      conFinal[which(manUnit %in% muName[m])] <- ifelse(summAbove == 1, 0, 1)
    }
  }
  constraintList <- list(muAboveFRP, conFinal)
  names(constraintList) <- c("muAboveFRP", "muConstrained")
  return(constraintList)
}

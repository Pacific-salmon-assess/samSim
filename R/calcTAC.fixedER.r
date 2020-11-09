#' Calculate total allowable catch for a fixed ER harvest policy
#'
#' This function calculates total allowable catch (TAC) for a given Conservation Unit (CU) using
#' a fixed exploitation rate. TAC is divided between one American and two Canadian
#' fisheries (mixed stock and single stock) based on \code{ppnMix} variable.
#'
#'
#' @param rec A numeric representing MU-specific return abundance.
#' @param canER A numeric representing the target Canadian exploitation rate
#' @param amER A numeric representing the target American exploitation rate
#' @param ppnMix A numeric representing the proportion of the Canadian TAC
#' allocated to mixed stock fisheries.
#' @return Returns a four element list of numeric vectors with length equal to
#' forecast:Total Canadian TAC, Canadian mixed fishery TAC, Canadian single fishery TAC, American TAC
#'
#' @examples
#' #Note that the function is intended to receive vectors rather than the DF
#' #used in this example to increase efficiency within the full closed-loop
#' simulation.
#' head(exampleHCRList)
#' names(exampleHCRList)[4] <- "recRYMU"
#' rec <- exampleHCRList$recRYMU

#' ## Fixed ER version
#' calcTAC(rec, canER = 0.4, canER=0.3, amER = 0.1, ppnMix = 1)
#'
#'
#' @export
calcTAC.fixedER <- function(rec, canER, amER, ppnMixVec) {

  canTAC <- canER * rec
  canMixTAC <- canTAC * ppnMixVec
  canSingTAC <- (canTAC * (1 -  ppnMixVec))
  amTAC<- amER * rec


  tacList <- list(canTAC, canMixTAC, canSingTAC, amTAC)
  tacList <- lapply(tacList, function (x){ #replace NAs with 0s
    x[is.na(x)] <- 0
    return(x)
  })
  names(tacList) <- c("canTAC", "canMixTAC", "canSingTAC", "amTAC")
  return(tacList)

}

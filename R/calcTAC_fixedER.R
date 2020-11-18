#' Title
#'
#' @param rec A numeric representing MU-specific return abundance
#' @param canER A numeric representing the target Canadian exploitation rate
#' @param amER A numeric representing the target American exploitation rate
#' @param ppnMixVec A numeric representing the proportion of the Canadian TAC
#' allocated to mixed stock fisheries
#'
#' @return Returns a four element list of numeric vectors with length equal to
#' forecast:Total Canadian TAC, Canadian mixed fishery TAC, Canadian single fishery TAC, American TAC
#' @export
#' @examples
#' calcTAC_fixedER(rec=1000, canER=0.3, amER = 0.1, ppnMix = 0.5)
#'
#'
calcTAC_fixedER <- function(rec, canER, amER, ppnMixVec) {

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

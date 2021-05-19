#' Title
#'
#' @param rec A numeric representing MU-specific return abundance
#' @param canER A numeric representing the target Canadian exploitation rate
#' @param amER A numeric representing the target American exploitation rate
#' @param ppnMixVec A numeric representing the proportion of the Canadian TAC
#' allocated to mixed stock fisheries
#' @param cvER A numeric representing annual variability in ER
#' @param randomVar A TRUE/FALSE variable indicating whether the Canadian ER should have
#' annual implementation error around the target
#' @param runif A vector of random numbers of length equal to the number of CUs
#' which provides constant CU-specific deviations in ER from annual values
#'
#' @return Returns a four element list of numeric vectors with length equal to
#' forecast:Total Canadian TAC, Canadian mixed fishery TAC, Canadian single fishery TAC, American TAC
#' @export
#' @examples
#' calcTAC_fixedER(rec=1000, canER=0.3, amER = 0.1, ppnMix = 0.5)
#'
#'
calcTAC_fixedER <- function(rec, canER, amER, ppnMixVec, cvER, randomVar=T, runif=NULL) {

  if (randomVar == F) {
    canTAC <- canER * rec
    canMixTAC <- canTAC * ppnMixVec
    canSingTAC <- (canTAC * (1 -  ppnMixVec))
    amTAC<- amER * rec
  }
  # At present, variable ERs are only applied to Canadian ER
  if (randomVar == T) {
  # calculate beta shape pars for can ER distribution

    sigCanER<-cvER*canER

    shape1<- canER^2 * (((1-canER)/sigCanER^2)-(1/canER))
    shape2<-shape1 * (1/canER-1)

    sampBeta<-function(stk) {
      x<-rbeta(1,shape1[stk],shape2[stk])
    }
    sampBetaRunif<-function(stk) {
      x<-qbeta(runif[stk],shape1[stk],shape2[stk])
    }

    # get realized ER
    if(is.null(runif)) canER.real<-sapply(1:length(sigCanER),sampBeta)

    if(!is.null(runif)) {
      canER.real <-  sapply(1:length(sigCanER),sampBetaRunif)
      # Align the random number seeds with the scenario is.null(runif)==TRUE
      # where rbeta is called. runif only uses 1 random seed/call, but rbeta
      # uses 2, so need to call another set of random numbers = nCU
      runif(length(cvER))
    }

    # if any CUs have a CV of 0, set to mean canER
    canER.real[sigCanER ==0]<-canER

    # calculate TACs
    canTAC <- canER.real * rec
    canMixTAC <- canTAC * ppnMixVec
    canSingTAC <- (canTAC * (1 -  ppnMixVec))
    amTAC<- amER * rec

  }

  #browser()

  tacList <- list(canTAC, canMixTAC, canSingTAC, amTAC, canER.real)
  tacList <- lapply(tacList, function (x){ #replace NAs with 0s
    x[is.na(x)] <- 0
    return(x)
  })
  names(tacList) <- c("canTAC", "canMixTAC", "canSingTAC", "amTAC", "canER.real")
  return(tacList)

}

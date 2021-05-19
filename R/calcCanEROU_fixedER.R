#' Title
#'
#' @param canER A numeric representing the target Canadian exploitation rate
#' @param cvERSMU A numeric representing annual variability in ER
#' @param randomVar A TRUE/FALSE variable indicating whether the Canadian ER should have
#' annual implementation error around the target
#'
#' @return Returns a four element list of numeric vectors with length equal to
#' forecast:Total Canadian TAC, Canadian mixed fishery TAC, Canadian single fishery TAC, American TAC
#' @export
#' @examples
#' calcCanEROU_fixedER(canER=0.3, cvERSMU=0.1)
#'
#'
calcCanEROU_fixedER <- function(canER,  cvERSMU, randomVar=T) {
  # At present, OU is only applied to Canadian ER

  if (randomVar == F) {
    canEROU <- canER
  }
  if (randomVar == T) {
    # calculate beta shape pars for can ER distribution

    sigCanER<-cvERSMU*canER

    shape1<- canER^2 * (((1-canER)/sigCanER^2)-(1/canER))
    shape2<-shape1 * (1/canER-1)

    sampBeta<-function(stk) {
      x<-rbeta(1,shape1[stk],shape2[stk])
    }
    # get realized ER
    canEROUR<-sapply(1:length(sigCanER),sampBeta)
    # if any CUs have a CV of 0, set to mean canER
    canEROU[sigCanER ==0]<-canER
  }
  return(canEROU)

}

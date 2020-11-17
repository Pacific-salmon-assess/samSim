#' Generate recruit abundance with Ricker model
#'
#' This function calculates recruitment from Ricker curve with AR(1) process
#' (according to Peterman et al. 2003; modified to take more recent parameter-
#' ization). Uses parameters from arima.mle (a, -b, sig, rho in log space) with
#' multivariate normally distributed errors. Note that internal \code{if}
#' statements prevent it from being vectorized so must be passed single values,
#' i.e. all vectors for inputs and outputs are length 1. Note that by default
#' prevErr and rho are NULL, resulting in a standard Ricker model.
#'
#' @param S A numeric vector of spawner abundances.
#' @param a A numeric vector of alpha values, i.e. productivity at low spawner
#' abundance.
#' @param b A numeric vector of beta values, i.e. density dependence para-
#' meter.
#' @param error A numeric vector of recruitment deviations, typically generated
#' using \code{rmvnorm()} and relevant process variance estimates (sigma).
#' @param rho A numeric vector of rho values, i.e. AR1 coefficient.
#' outside of model using multivariate normal (or equivalent) distribution.
#' @param prevErr A numeric vector representing recruitment deviations from
#' previous brood year.
#' @return Returns a list of R, a numeric representing recruit abundance, and
#' \code{errNext} which is used to generate subsequent process error (i.e. next
#' year's prevErr.
#'
#' @examples
#' #Spawner and recruit values represent millions of fish, stock-recruit
#' parameters approximate those of Fraser River sockeye salmon Chilko CU.
#'
#' #without autoregressive error
#' rickerModel(S = 1.1, a = 1.8, b = 1.2, error = 0.3)
#'
#' #with autoregressive error
#' rickerModel(S = 1.1, a = 1.8, b = 1.2, error = 0.3, rho = 0.2,
#' prevErr = 0.7)
#'
#' @export

rickerModel <- function(S, a, b, error, rho = NULL, prevErr = NULL) {
  if (is.null(rho)) {
    rho <- 0
  }
  if (is.null(prevErr)) {
    prevErr <- 0
  }
  err <- prevErr * rho + error
  # if (a >= 0) {
  #   if (b != 0 & S > 0) {
      R <- S * exp(a - b * S) * exp(err)
      errNext <- log(R / S) - (a - b * S)
  #   }
  #   if (b == 0 & S > 0) {
  #     R <- S * exp(err)
  #     errNext <- log(R / S)
  #   }
  # }
  # if (a < 0 & S > 0) {
  #   R <- S * exp(a) * exp(error)
  #   errNext <- log(R / S)
  # }
  # if (S == 0) {
  #   R <- 0
  #   errNext <- err
  # }
  return(list(R, errNext))
}

#------------------------------------------------------------------------------

#' Generate recruit abundance with Larkin model
#'
#' This function calculates recruitment from Larkin model (according to
#' Peterman et al. 2003, modified to take more recent parameter-
#' ization). Uses parameters in log space, like rickerMod, with
#' multivariate normally distributed errors, but cannot incorporate AR1 process
#' error because such models have not been validated.
#' @section Note: the log-normal bias correction has not been fixed for the
#' Larkin model.

#' @param S A numeric vector of spawner abundances.
#' @param Sm1,Sm2,Sm3 A numeric vector of spawner abundances at 1, 2 and 3 year
#' lags, respectively.
#' @param a A numeric vector of alpha values, i.e. productivity at low spawner
#' abundance.
#' @param b A numeric vector of beta values, i.e. density dependence para-
#' meter.
#' @param b1,b2,b3 A numeric vector of delayed density dependent effects at 1,
#' 2, and 3 yera lags, respectively.
#' @param error A numeric vector of recruitment deviations, typically generated
#' using \code{rmvnorm()} and relevant process variance estimates (sigma).
#' @return Returns a numeric representing recruit abundance.
#'
#' @examples
#' #Spawner and recruit values represent millions of fish, parameters
#' approximate those of Shuswap CU
#' larkinModel(S = 1.1, Sm1 = 0.4, Sm2 = 0.2, Sm3 = 0.15, a = 2.2, b = 0.29,
#' b1 = 0.42, b2 = 0.31, b3 = 0.21, error = 0.3)
#'
#' @export

larkinModel <- function(S, Sm1, Sm2, Sm3, a, b, b1, b2, b3, error) {
  R <- (S * exp(a - b * S - b1 * Sm1 - b2 * Sm2 - b3 * Sm3)) * exp(error)
  return(R)
}

#------------------------------------------------------------------------------

#' Generate recruit abundance with the Ricker Model with a marine survival co-variate
#'
#' This function calculates recruitment from Ricker curve with a marine survival covariate
#' that is specific to brood year added in.  Note that internal \code{if}
#' statements prevent it from being vectorized so must be passed single values,
#' i.e. all vectors for inputs and outputs are length 1.
#'
#' @param S A numeric of spawner abundance.
#' @param a A numeric of alpha value, i.e. productivity at low spawner
#' abundance.
#' @param b A numeric vector of beta value, i.e. density dependence para-
#' meter.
#' @param error A numeric recruitment deviation, typically generated
#' using \code{rmvnorm()} and relevant process variance estimates (sigma).
#' @param ppnAges A numeric vector of proportion of spawner abundance at age
#' @param gamma A numeric value represting the marine survival coefficient
#' @param mSurvAtAge A numeric vector representing marine survival covariates at age
#'
#' @return A list of recruitment at each age for the modelled brood year
#' (e.g., R2 = recruitment to age 2 from broodyear, R3 = recruitment to age 3, etc),
#' as well as recruitment for all return ages combined from the modelled brood year (RecBY)
#' @export
#'
#' @examples
#' rickerSurvModel(S = 1000, a = 2.1, b = 0.00001, error = 0.8, ppnAges = c(0,0.83,0.17,0,0), gamma=0.4, mSurvAtAge=c(0,0.013, 0.015,0,0))
#'
rickerSurvModel <- function(S, a, b, error, ppnAges, gamma, mSurvAtAge) {

  R2 <- ppnAges[1] * S * exp(a-b*S + gamma * mSurvAtAge[1])  * exp(error)
  R3 <- ppnAges[2] * S * exp(a-b*S + gamma * mSurvAtAge[2]) * exp(error)
  R4 <- ppnAges[3] * S * exp(a-b*S + gamma * mSurvAtAge[3]) * exp(error)
  R5 <- ppnAges[4] * S * exp(a-b*S + gamma * mSurvAtAge[4]) * exp(error)
  R6 <- ppnAges[5] * S * exp(a-b*S + gamma * mSurvAtAge[5]) * exp(error)

  RecBY<-sum(R2,R3,R4,R5,R6)

  return(list(R2,R3,R4,R5,R6,RecBY))
}



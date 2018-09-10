#' Generate recruit abundance with Ricker model
#'
#' This function calculates recruitment from Ricker curve with AR(1) process
#' (according to Peterman et al. 2003; modified to take more recent parameter-
#' ization). Uses parameters from arima.mle (a, -b, sig, rho in log space) with
#' multivariate normally distributed errors. Note that internal \code{if}
#' statements prevent it from being vectorized so must be passed single values,
#' i.e. all vectors for inputs and outputs are length 1. Note that by default
#' utminus1 and rho are NULL, resulting in a standard Ricker model.
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
#' @param utminus1 A numeric vector representing recruitment deviations from
#' previous brood year.
#' @return Returns a list of R, a numeric representing recruit abundance, and
#' \code{err.next} which is used to generate subsequent process error.
#'
#' @examples
#' #Spawner and recruit values represent millions of fish, parameters loosely
#' #based on Chilko CU
#' #without autoregressive error
#' rickerModel(S = 1.1, a = 1.8, b = 1.2, error = 0.3)

#' #with autoregressive error
#' rickerModel(S = 1.1, a = 1.8, b = 1.2, error = 0.3, rho = 0.2,
#' utminus1 = 0.7)
#'
#' @export

rickerModel <- function(S, a, b, error, rho = NULL, utminus1 = NULL) {
  err <- utminus1 * rho + error
  if (a >= 0) {
    if (b != 0 & S > 0) {
      R <- S * exp(a - b * S) * exp(err)
      err.next <- log(R / S) - (a - b * S) + error
    }
    if (b == 0 & S > 0) {
      R <- S * exp(err)
      err.next <- log(R / S) - 0
    }
  }
  if (a < 0 & S > 0) {
    R <- S * exp(a) * exp(error)
    err.next <- log(R / S) - 0 + error
  }
  if (S == 0) {
    R <- 0
    err.next <- err
  }
  return(list(R, err.next))
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

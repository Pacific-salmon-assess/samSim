#' Posterior estimates of Ricker stock recruit parameters
#'
#' A dataset containing posterior estimates of stock recruit parameters for CU-
#' specific Ricker models fit to Fraser River sockeye salmon conservation units
#' (CUs).
#'
#' @format A data frame with 19000 rows and 5 variables:
#' \describe{
#'   \item{stk}{stock id, approximately equivalent to a conservation unit}
#'   \item{alpha}{productivity parameter (i.e. recruits per spawner at low
#'   spawner abundances)}
#'   \item{beta0}{density dependence parameter}
#'   \item{sigma}{standard deviation of normally distributed process error with
#'   mean = 0}
#'   \item{deviance}{estimate of model fit equivalent to residual sum of
#'   squares for linear models}
#' }
#' @source Fraser River Salmon Spawning Initiative (Ann-Marie Huang (DFO) personal
#' communication; January 2018)
"rickerParameters"

#______________________________________________________________________________

#' Posterior estimates of Ricker stock recruit parameters
#'
#' A dataset containing posterior estimates of stock recruit parameters for CU-
#' specific Larkin models fit to Fraser River sockeye salmon conservation units
#' (CUs). Larkin models are extensions of Ricker models that account for de-
#' layed density dependent effects by incorporating additional beta parameters.
#'
#' @format A data frame with 19000 rows and 8 variables:
#' \describe{
#'   \item{stk}{stock id, approximately equivalent to a conservation unit}
#'   \item{alpha}{productivity parameter (i.e. recruits per spawner at low
#'   spawner abundances)}
#'   \item{beta0}{density dependence parameter for spawners in current year}
#'   \item{beta1}{density dependence parameter for spawners one year prior}
#'   \item{beta2}{density dependence parameter for spawners two years priod}
#'   \item{beta3}{density dependence parameter for spawners three years prior}
#'   \item{sigma}{standard deviation of normally distributed process error with
#'   mean = 0}
#'   \item{deviance}{estimate of model fit equivalent to residual sum of
#'   squares for linear models}
#' }
#' @source Fraser River Salmon Spawning Initiative (Ann-Marie Huang (DFO) personal
#' communication; January 2018)
"larkinParameters"

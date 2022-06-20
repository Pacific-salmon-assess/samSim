

#' compute length of regimes
#'
#' This funncion computes the sequence of regime scalars based on used defines, 
#' low, high scallars as well as  regime length
#'
#'
#' @param a1 Parameter value for the first regime. Value to be multiplied by the average parameter
#' @param a2  Parameter value for the second regime.
#' @param y length of regime trend
#' @param reglen length of each regime, default is 10
#' @return Returns a vector of length y with alternating sequence of regime scalars
#'
#' @examples
#' runRegime(.75,1.25,23,5)
#'
#' 
#' @export
runRegime <- function(a1,a2,y,reglen=10)( rep( c(rep(a1,reglen),rep(a2,reglen)),
                                       ceiling(y/(reglen*2)) )[1:y] )

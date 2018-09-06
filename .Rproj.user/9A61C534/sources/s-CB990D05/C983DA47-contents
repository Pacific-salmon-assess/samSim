#' Calculate quantiles
#'
#' Functions to calculate lower and upper quantiles.
#'
#' This is a generic wrapper function for quantile that defaults to specific values and removes
#' NAs.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector.
#'
#' @examples
#' qLow(c(3, 12, NA, 1, 7, 25))
#' qHigh(c(3, 12, NA, 1, 7, 25))
#'
#' @export
qLow <- function(x) {
  q <- quantile(x, probs = 0.10, na.rm = TRUE)
  return(q)
}

qHigh <- function(x) {
  q <- quantile(x, probs = 0.90, na.rm = TRUE)
  return(q)
}

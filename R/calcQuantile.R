#' Calculate lower quantile
#'
#' Two simple wrapper functions for \code{quantile()} that are used in plotting
#' functions. qLow defaults to the 5th percentile and qHigh to the 95th.
#'
#' @param x A numeric vector.
#' @return Returns a numeric representing either the 5th or 95th percentile.
#'
#' @export

qLow <- function(x) {
  q <- quantile(x, probs = 0.05, na.rm = TRUE)
  return(q)
}

#______________________________________________________________________________

#' Calculate upper quantile
#'
#' Two simple wrapper functions for \code{quantile()} that are used in plotting
#' functions. qLow defaults to the 5th percentile and qHigh to the 95th.
#'
#' @param x A numeric vector.
#' @return Returns a numeric representing either the 5th or 95th percentile.
#'
#' @export

qHigh <- function(x) {
  q <- quantile(x, probs = 0.95, na.rm = TRUE)
  return(q)
}

#' Calculate lower quantile
#'
#' Two simple wrapper functions for \code{quantile()} that are used in plotting
#' functions. qLow defaults to the 10th percentile and qHigh to the 90th.
#'
#' @param x A numeric vector.
#' @return Returns a numeric representing either the 10th or 90th percentile.
#'
#' @export

qLow <- function(x) {
q <- quantile(x, probs = 0.10, na.rm = TRUE)
return(q)
}

#______________________________________________________________________________

#' Calculate upper quantile
#'
#' Two simple wrapper functions for \code{quantile()} that are used in plotting
#' functions. qLow defaults to the 10th percentile and qHigh to the 90th.
#'
#' @param x A numeric vector.
#' @return Returns a numeric representing either the 10th or 90th percentile.
#'
#' @export

qHigh <- function(x) {
  q <- quantile(x, probs = 0.90, na.rm = TRUE)
  return(q)
}

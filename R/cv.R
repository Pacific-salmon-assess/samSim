#' cv
#'
#' Simple helper function calculate coefficient of variation.
#'
#' @param x A numeric vector.
#' @return Returns a numeric
#'
#' @examples
#' x <- rnorm(50, mean = 0, sd = 1)
#' cv(x)
#' @export
cv <- function(x){
  cv <- sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE)
  return(cv)
}

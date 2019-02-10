#' Calculate observed error
#'
#' This helper function calculates relative error rates with a few corrections
#' to account for missing or wonky data.
#'
#' @param obs A numeric vector representing observed data.
#' @param true A numeric vector representing true data.
#' @return Returns a numeric vector of absolute relative error.
#'
#' @export
#'
calcErr <- function(obs, true) {
  err <- ifelse((obs - true) / true == Inf |
                  (obs - true) / true == -Inf,
                NA, (obs - true) / true)
  err <- ifelse(abs(err) < 0.00001, 0, err)
  return(err)
}

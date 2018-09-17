#' Generate cycle lines
#'
#' This function generates cycle lines necessary to model Fraser River stocks
#' exhibiting delayed density dependence (e.g. Lower Shuswap). Note that it is
#' not a universal function - it can only generate four year cycles and assumes
#' that the observed data begins after 1900.
#'
#' Depends on the helper function \code{checkInteger}.
#'
#' @param firstObs A numeric representing the first year of observed data.
#' @param simLength A numeric representing the length of the simulation period.
#' @return Returns a numeric vector representing cycle line IDs for output data
#' with \code{length(simLength)}.
#'
#' @examples
#' genCycle(firstObs = 2009, simLength = 50)
#'
#' @export

genCycle <- function(firstObs, simLength){
  if (firstObs < 1900) {
    stop("First observation too early to adjust cycle lines. Start date must be
         after 1900.")
  }
  if (checkInteger((firstObs - 1900) / 4)) {
    cycle <- rep(c(4, 1, 2, 3), length.out = simLength)
  }
  if (checkInteger((firstObs + 1 - 1900) / 4)){
    cycle <- rep(c(3, 4, 1, 2), length.out = simLength)
  }
  if (checkInteger((firstObs + 2 - 1900) / 4)){
    cycle <- rep(c(2, 3, 4, 1), length.out = simLength)
  }
  if (checkInteger((firstObs + 3 - 1900) / 4)){
    cycle <- rep(c(1, 2, 3, 4), length.out = simLength)
  }
  return(cycle)
}

#______________________________________________________________________________

#' Check integer
#'
#' Helper function to calculate which cycle line an observed year belongs to. Used in
#' \code{genCygle}.
#'
#' @param x A numeric.
#' @return Returns a logical confirming whether x is or is not an integer.
#'
#' @examples
#' checkInteger(3.4)
#' checkInteger(3)
#'
#' @export
checkInteger <- function(x) {
  x == round(x)
}

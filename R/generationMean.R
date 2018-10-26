#' Generational mean
#'
#' Calculate running mean of a vector over generation length.
#'
#' @param vec A numeric vector, typically of spawner abundances.
#' @param gen A numeric representing generation length (years).
#' @return Returns a numeric vector of the running mean
#'
#' @examples
#' rec <- recMatrix[ , 1]
#' genMean(rec, gen = 4)
#'
#' @export

genMean <- function(vec, gen) {
  runningMean <- NA
  for (i in gen:length(vec)) {
    runningMean[i]<-(prod(vec[(i - gen + 1):i])) ^ (1 / gen)
  }
  return(runningMean)
}

#' Infill spawner/recruit abundance
#'
#' This function infills matrices based on observed geometric mean abundances
#' and relative size of other CUs. In other words abundance is a function of
#' average stock size and the observed return to observed stocks in a given
#' year.
#'
#' Generally used with stocks that have gappy data, particularly in recent
#' years because the closed-loop simulation model requires a seamless tran-
#' sition between retrospective and forward simulated data.
#'
#' Depends on helper function \code{geoMean} used to calculate geometric means.
#'
#' @param mat A matrix of observed abundances with \code{nrow} equal to number
#' of years and \code{ncol} equal to the number of stocks.
#' @return Returns a matrix of the same dimensions with NAs replaced by infilled
#' values
#'
#' @examples
#' set.seed(123)
#' #matrix of random time series representing stock abundance
#' popDat <- cbind(popA = round(exp(rnorm(20, 7, 1))),
#'                 popB = round(exp(rnorm(20, 8, 1))),
#'                 popC = round(exp(rnorm(20, 9, 1)))
#' )
#' #remove 20% of values
#' gaps <- floor(runif(n = 0.2 * length(popDat), min = 1, max = length(popDat)))
#' popDat[gaps] <- NA
#'
#' @export
infill <- function(mat) {
  tsLength <- min(25, nrow(mat))
  meanAbund <- apply(mat[(nrow(mat) - tsLength):nrow(mat), ], 2, geoMean)
  ppnAbund <- matrix(meanAbund/sum(meanAbund), nrow = nrow(mat),
                     ncol = ncol(mat), byrow = TRUE)
  present <- ifelse(is.na(mat), 0, 1)
  ppnPresent <- ppnAbund * present
  expansion <- apply(ppnPresent, 1, function(x) 1 / sum(x))
  expandedTotal <- apply(mat, 1, function(x) sum(x, na.rm = TRUE)) * expansion
  infillMat <- ppnAbund * expandedTotal
  return(infillMat)
}

#______________________________________________________________________________

#' Calculate geometric mean
#'
#' Helper function to calculate geometric means. Used in \code{infill}.
#'
#' @param x A numeric vector.
#' @return Returns a numeric equal to geometric mean
#'
#' @examples
#' x <- c(2, 5, 12, 9)
#' geoMean(x)

#' @export
geoMean <- function(x) {
  xx <- x[which(is.na(x) == F)]
  exp(mean(log(xx)))
}

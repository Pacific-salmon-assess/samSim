#' Truncated normal distribution
#'
#' This function calculates means across an array's third dimension, then its
#' second dimension.
#'
#' @param y A numeric array.
#' @return A numeric vector representing the mean of means with length equal to
#' \code{ncol(y)}.
#'
#' @examples
#' a <- matrix(rnorm(5 * 4, mean = 0, sd = 1), 5, 4)
#' b <- matrix(rnorm(5 * 4, mean = 1, sd = 1), 5, 4)
#' c <- matrix(rnorm(5 * 4, mean = 2, sd = 1), 5, 4)
#' y <- sapply(list(a, b, c), identity, simplify="array")
#' arrayMean(y)
#' @export
arrayMean <- function(y) {
  nCU <- dim(y)[2]
  datOut <- rep(NA, length.out = nCU)
  for (i in 1:nCU) {
    trialMeans <- apply(y[ , i, ], 1, mean)
    datOut[i] <- mean(trialMeans, na.rm = TRUE)
  }
  return(datOut)
}

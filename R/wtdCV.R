#' Calculate weighted component CV
#'
#' This function calculates an estimate of the mean coefficient of variation
#' among an aggregate's components. The estimate is weighted by each components
#' mean abundance. Based on equation 4 in Thibaut and Connolly 2013.
#'
#' @param datMat A numeric matrix, generally representing a time series of
#' component specific abundances.
#' @param weight A logical argument that specifies whether CV should be
#' weighted or not.
#' @return Returns a numeric representing mean CV.
#'
#' @examples
#' r <- recMatrix[1:10, ]
#' wtdCV(r, weight = TRUE)
#'
#' @export

wtdCV <- function(datMat, weight = TRUE) {
  # need to replace any all 0 columns with v. small numbers for following calc
  # to work, but want to make sure NAs aren't changed
  for (i in 1:ncol(datMat)) {
    dum <- datMat[ , i]
    if (sum(dum, na.rm = TRUE) == 0) {
      dum[dum == 0] <- 1e-9
    }
    datMat[ , i] <- dum
  }

  #temporal mean of aggregate abundance
  aggAbund <- sum(apply(datMat, 2, function (x) mean(x)))
  #wtd mean of aggregate abundance
  wtdAbund <- apply(datMat, 2, function (x) mean(x) / aggAbund)
  wtdCV <- sum(wtdAbund * apply(datMat, 2,
                                function(x) sqrt(var(x)) /
                                  mean(x)))
  unWtdCV <- mean(apply(datMat, 2,
                        function(x) sqrt(var(x, na.rm = TRUE)) /
                          mean(x, na.rm = TRUE)), na.rm = TRUE)
  if (weight == TRUE) {
    return(wtdCV)
  } else {
    return(unWtdCV)
  }
}

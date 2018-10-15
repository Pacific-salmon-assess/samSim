#' Calculate weighted component CV
#'
#' This function calculates an estimate of the mean coefficient of variation
#' among an aggregate's components. The estimate is weighted by each components
#' mean abundance. Based on equation 4 in Thibaut and Connolly 2013.
#'
#' @param datMat A numeric matrix, generally representing a time series of
#' component specific abundances.
#' @param weightMat An optional numeric matrix that is passed if mean CV is
#' being weighted by a different variable than that in datMat. **Note that**
#' **if weightMat is passed, aggregate variability will no longer equal**
#' **sqrt(phi) * CVc.**
#' @param weight A logical argument that specifies whether CV should be
#' weighted or not. **Note that if weightMat is passed, aggregate variability**
#' **will no longer equal sqrt(phi) * CVc.**
#' @return Returns a numeric representing mean CV.
#'
#' @examples
#' r <- recMatrix[1:10, ]
#' wtdCV(r, weightMat = NULL, weight = TRUE)
#' @export
wtdCV <- function(datMat, weightMat = NULL, weight = TRUE) {
  if (is.null(weightMat)) { #if weightMat is NULL assume datMat is a matrix of abundance
    weightMat <- datMat
  }
  if (ncol(datMat) != ncol(weightMat)){
    stop("Input matrices have unequal number of components")
  }
  if (any(is.na(datMat))) {
    warning("NAs present. This will affect estimates of weighted CV.")
  }
  #temporal mean of aggregate abundance
  aggAbund <- sum(apply(weightMat, 2, function (x) mean(x, na.rm = TRUE)))
  #wtd mean of aggregate abundance
  wtdAbund <- apply(weightMat, 2, function (x) mean(x, na.rm = TRUE) / aggAbund)
  wtdCV <- sum(wtdAbund * apply(datMat, 2,
                                function(x) sqrt(var(x, na.rm = TRUE)) /
                                  mean(x, na.rm = TRUE)), na.rm = TRUE)
  unWtdCV <- mean(apply(datMat, 2,
                        function(x) sqrt(var(x, na.rm = TRUE)) /
                          mean(x, na.rm = TRUE)), na.rm = TRUE)
  if (weight == TRUE) {
    return(wtdCV)
  } else {
    return(unWtdCV)
  }
}

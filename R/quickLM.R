#' Linear model
#'
#' This function uses C++ to quickly fit linear models. Depends on fastLm.cpp.
#'
#' @useDynLib samSim
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#'
#' @param xVec,yVec Numeric vectors of equal length representing predictor and
#' response variables in linear model.
#' @return Returns a list of model coefficients.
#'
#' @examples
#' x <- rnorm(10, 0, 1)
#' y <- 2*x + 1
#' quickLm(x, y)
#'
#' @export
quickLm <- function(xVec, yVec){
  #C++ won't accept mix of NAs so trim and convert to matrix w/ 1s for intercept
  xVec2 <- xVec[!is.na(xVec + yVec)]
  # xMat <- cbind(1, xVec2)
  xMat <- xVec2
  yMat <- as.numeric(yVec[!is.na(xVec + yVec)])
  mod <- RcppArmadillo::fastLm(yMat ~ xMat)

  return(mod$coefficients)
}

#______________________________________________________________________________

#' Unload DLL 
#'
#' This function unloads the dynamic library whenever the package is unloaded.

.onUnload <- function (libpath) {
  library.dynam.unload("samSim", libpath)
}

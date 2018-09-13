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
#' y <- 2x + 1
#' fastLm(x, y)

#' @export
#' #Function that cleans input data then uses Rcpp to provide faster lm estimates
quickLm <- function(xVec, yVec){
  #C++ won't accept mix of NAs so trim and convert to matrix w/ 1s for intercept
  xVec2 <- xVec[!is.na(xVec + yVec)]
  xMat <- matrix(c(rep(1, length(xVec2)), xVec2), ncol = 2)
  yMat <- yVec[!is.na(xVec + yVec)]
  mod <- fastLm(yMat, xMat)

  return(mod$coefficients)
}

#' Extinction check
#'
#' This function checks whether a population has gone quasi-extinct based on
#' spawner abundances, generation length, and the extinction threshold.
#'
#' @param y A numeric representing the current year within the simulation
#' (i.e. has a quasi-extinction occurred based on abundance in y-gen:y).
#' @param gen A numeric representing generation length for the species being
#' modeled.
#' @param extinctThresh A numeric representing the assumed quasi-extinction
#' threshold.
#' @param spwnMat A matrix of spawner abundances with ncol = nCU and minimum
#' nrow = gen.
#' @return Returns a binary vector of length = ncol(spwnMat) with ones
#' representing extinct stocks/CUs.
#'
#' @examples
#' nCUs <- 5
#' gen <- 4
#' yr <- 6
#' sMat <- matrix(0, nrow = yr, ncol = nCUs)
#' for (i in 1:nCUs) {
#'   sMat[ , i] <- exp(rnorm(6, mean = 0.5, sd = 1))
#' }
#' sMat[ , 5] <- c(0.1, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001)
#' extinctionCheck(y = yr, gen = gen, extinctThresh = 0.001, spwnMat = sMat)
#' @export
extinctionCheck <- function(y, gen, extinctThresh, spwnMat) {
  nPops <- ncol(spwnMat)
  extVec <- rep(0, length.out = nPops)
  for (k in 1:nPops) {
    if (gen == 2) {
      if (spwnMat[y, k] < extinctThresh & spwnMat[y - 1, k] < extinctThresh) {
        extVec[k] <- 1
      }
    }
    if (gen == 3) {
      if (spwnMat[y, k] < extinctThresh & spwnMat[y - 1, k] < extinctThresh
          & spwnMat[y-2, k] < extinctThresh) {
        extVec[k] <- 1
      }
    }
    if (gen == 4) {
      if (spwnMat[y, k] < extinctThresh & spwnMat[y - 1, k] < extinctThresh
          & spwnMat[y-2, k] < extinctThresh & spwnMat[y - 3, k] < extinctThresh) {
        extVec[k] <- 1
      }
    }
    if (gen == 5) {
      if (spwnMat[y, k] < extinctThresh & spwnMat[y - 1, k] < extinctThresh
          & spwnMat[y-2, k] < extinctThresh & spwnMat[y - 3, k] < extinctThresh
          & spwnMat[y - 4, k] < extinctThresh) {
        extVec[k] <- 1
      }
    }
    if (gen == 6) {
      if (spwnMat[y, k] < extinctThresh & spwnMat[y - 1, k] < extinctThresh
          & spwnMat[y-2, k] < extinctThresh & spwnMat[y - 3, k] < extinctThresh
          & spwnMat[y - 4, k] < extinctThresh & spwnMat[y - 5, k] < extinctThresh) {
        extVec[k] <- 1
      }
    }
  }
  return(extVec)
}

#' Summarize stock-recruit parameters
#'
#' This function samples relevant stock-recruit parameters from a dataframe
#' containing posterior samples from MCMC sampling. Sampling is based on per-
#' centile arguments.
#'
#' @param pars A dataframe containing stock-recruit parameter. Should be
#' formatted so that stock/population is in first column and each subsequent
#' column contains one parameter; each row represents one sample.
#' @param alphaOnly A logical argument. If TRUE then percentiles are only used
#' for alpha parameter, all others use medians.
#' @param highP,lowP Numerics that represent high and low percentiles that are
#' being estimated.
#' @param stks A vector that is used to identify a subset of stocks for which
#' are needed. Can be either a numeric or character, but should match the
#' format of the first column of the input dataframe.
#' @return Returns a numeric that is the spawner abundance that minimizes the
#' the log likelihood.
#'
#' @examples
#' Stock-recruit parameters approximate those of Fraser River sockeye salmon
#' Chilko CU.
#' alpha = 1.2
#' beta = 1.5
#' sigma = 0.8
#' theta = c(alpha, beta, sigma)
#' sMSY = 0.3
#' sGenSolver(theta, sMSY)
#' @export
#'
#' #This function samples appropriate SR parameters based on arguments passed from simPar.csv and
#a datafile containing Ricker and/or Larkin parameters; if alphaOnly = TRUE, betas and sigmas are
#set at median values
getSRPars <- function(pars, alphaOnly = TRUE, highP = 0.9, lowP = 0.1, stks = NULL) {
  srLow <- NULL
  srMed <- NULL
  srHigh <- NULL
  if (is.null(stks) == FALSE) {
    pars <- pars[pars$stk %in% stks, ]
  }
  stkKey <- unique(pars$stk)
  perc <- pars %>%
    group_by(stk) %>%
    summarise(alphaHigh = quantile(alpha, probs = highP, na.rm = TRUE),
              alphaMed = quantile(alpha, probs = 0.5, na.rm = TRUE),
              alphaLow = quantile(alpha, probs = lowP, na.rm = TRUE)
    )

  if (alphaOnly == TRUE) { #calculate alpha percentiles and other pars medians
    meds <- pars %>%
      dplyr::select(-c(alpha)) %>% #need dplyr because select sometimes masked by other packages
      group_by(stk) %>%
      summarise_all(median)
    srLow <- merge(perc[, c("stk", "alphaLow")], meds, by = "stk")
    srMed <- merge(perc[, c("stk", "alphaMed")], meds, by = "stk")
    srHigh <- merge(perc[, c("stk", "alphaHigh")], meds, by = "stk")
  }

  if (alphaOnly == FALSE) { #calculate alpha percentiles and sample other pars accordingly
    for (k in seq_along(stkKey)) {
      d <- subset(pars, stk == stkKey[k])
      dLow <- d[which.min(abs(d$alpha - perc$alphaLow[k])), ] #ID and pull row for each stock that most closely matches percentile
      dMed <- d[which.min(abs(d$alpha - perc$alphaMed[k])), ]
      dHigh <- d[which.min(abs(d$alpha - perc$alphaHigh[k])), ]
      srLow <- rbind(srLow, dLow)
      srMed <- rbind(srMed, dMed)
      srHigh <- rbind(srHigh, dHigh)
    }
  }
  srList <- list(srLow, srMed, srHigh)
  names(srList) <- c("pLow", "pMed", "pHigh")
  srList <- lapply(srList, function (x) {
    names(x)[2] <- "alpha"
    x
  })
  return(srList)
}

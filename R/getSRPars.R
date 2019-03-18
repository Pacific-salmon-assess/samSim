#' Summarize stock-recruit parameters
#'
#' This function samples relevant stock-recruit parameters from a dataframe
#' containing posterior samples from MCMC sampling. Sampling is based on
#' percentile arguments.
#'
#' @importFrom dplyr filter group_by select summarise summarise_all
#'
#' @param pars A dataframe containing stock-recruit parameters. Should be
#' formatted so that stock/population (\code{stk}) is in first column and each
#' subsequent column contains one parameter; each row represents one posterior
#' sample.
#' @param alphaOnly A logical argument. If TRUE then percentiles are only used
#' for alpha parameter, all others use medians. If FALSE then percentiles are
#' calculated and returned for all parameters.
#' @param highP,lowP Numerics that represent high and low percentiles that are
#' being estimated.
#' @param stks A vector that is used to identify a subset of stocks for which
#' parameters are returned. Can be either a numeric or character, but should
#' match the format of the first column of the input dataframe.
#' @return Returns a list of length three (\code{pLow, pMed, pHigh}) containing
#' lower percentiles, medians and higher percentiles, respectively. Note that
#' if alphaOnly = TRUE, the same median values will be returned for non-alpha
#' parameters in all list elements.
#'
#' @examples
#' head(rickerParameters)
#' d <- getSRPars(rickerParameters, alphaOnly = FALSE)
#' head(d$pLow) #lower percentile estimates
#' head(d$pMed) #median estimates
#'
#' d <- getSRPars(larkinParameters, alphaOnly = TRUE)
#' head(d$pLow) #lower percentile estimates for alpha only
#' @export

getSRPars <- function(pars, alphaOnly = TRUE, highP = 0.9, lowP = 0.1,
                      stks = NULL) {
  srLow <- NULL
  srMed <- NULL
  srHigh <- NULL
  if (!is.null(stks)) {
    pars <- pars %>%
      filter(stk %in% stks)
  }
  if (!is.null(pars$CU)) {
    pars <- pars %>%
      select(-CU)
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
      select(-c(alpha)) %>%
      group_by(stk) %>%
      summarise_all(median)
    srLow <- merge(perc[, c("stk", "alphaLow")], meds, by = "stk")
    srMed <- merge(perc[, c("stk", "alphaMed")], meds, by = "stk")
    srHigh <- merge(perc[, c("stk", "alphaHigh")], meds, by = "stk")
  }

  if (alphaOnly == FALSE) { #calculate alpha percentiles and sample other pars accordingly
    for (k in seq_along(stkKey)) {
      d <- subset(pars, stk == stkKey[k])
      #ID and pull row for each stock that most closely matches percentile
      dLow <- d[which.min(abs(d$alpha - perc$alphaLow[k])), ]
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

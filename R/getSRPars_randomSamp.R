#' Sample stock-recruit parameters from posterior distribution
#'
#' This function samples relevant stock-recruit parameters from a dataframe
#' containing posterior samples from MCMC sampling. Sampling among MCMC replicates is
#' random with replacement.
#'
#' @importFrom dplyr filter group_by select summarise summarise_all
#'
#' @param pars A dataframe containing stock-recruit parameters. Should be
#' formatted so that stock/population (\code{stk}) is in first column and each
#' subsequent column contains one parameter; each row represents one posterior
#' sample.
#' @param stks A vector that is used to identify a subset of stocks for which
#' parameters are returned. Can be either a numeric or character, but should
#' match the format of the first column of the input dataframe.
#' @return Returns a data frame with SR parameters for each stock that come from a single, randomly
#' sampled MCMC replicate.Each row in the data frame represents a stock.
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

getSRPars_randomSamp <- function(pars,
                                 stks = NULL) {
  srSamp<-NULL
  if (!is.null(stks)) {
    pars <- pars %>%
      filter(stk %in% stks)
  }
  if (!is.null(pars$CU)) {
    pars <- pars %>%
      select(-CU)
  }

  stkKey <- unique(pars$stk)

  # Add column identifying mcmc replicate number
  dum<-pars %>% group_by(stk) %>% mutate(rep=row_number())

  #Select mcmc replicate to sample
  sampID<-round(runif(1,1,max(dum$rep)),0)

  samp<-dum %>% filter(rep==sampID)

  SRout <- data.frame(alpha=samp$alpha, beta=samp$beta, sigma=samp$sigma, gamma=samp$gamma)

  return(SRout)

}

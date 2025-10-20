#' Metrics relating to harvest control rules
#'
#' A list containing example output data from a \code{recoverySimulator}
#' simulation run that is intended to be passed to \code{calcSynchMetrics}.
#'
#' @format a data frame.
#' \describe{
#'   \item{year}{Simulation Year}
#'   \item{iteration}{Simulation iteration}
#'   \item{CU}{Simulated CU number}
#'   \item{aboveLowerBM}{0 and 1, indicating if spawners are above the true lower benchmark }
#'   \item{aboveUpperBM}{0 and 1, indicating if spawners are above the true upper benchmark}
#'   \item{aboveLowerObsBM}{0 and 1, indicating if spawners are above the estimated lower benchmark }
#'   \item{aboveUpperObsBM}{0 and 1, indicating if spawners are above the estimated upper benchmark}
#'   \item{totalCatch}{total catch}
#'   \item{totalCACatch}{Canadian catches}
#'   \item{totalUSCatch}{Us catches}
#'   \item{sMSYEst}{yearly estimates of Smsy}
#'   \item{sGenEst}{yearly estimates of Sgen}
#'   \item{uMSyEst}{yearly estimates of Umsy}
#'   \item{lowerBM}{true lower benchmark}
#'   \item{upperBM}{true upper benchmark}
#'   \item{lowerObsBM}{yearly estimated lower benchmark}
#'   \item{upperObsBM}{yearly estimated upper benchmark}
#'   \item{UmsyBM}{opimum harvest rate, updated with the assessment schedule }
#'   \item{scenario}{string with combination of OM amd MP}
#'   \item{nameMP}{Management procedure string}
#'   \item{aboveLowerassess}{}
#'   \item{aboveUpperassess}{}
#'    \item{HCRtype}{}
#'  \item{wsp.status}{}
#'  \item{status}{}
#'  \item{status_agg}{}
#'  \item{plotOM}{}
#'  \item{plotMP}{}
#'  \item{HCR}{}
#'  \item{rp_type}{}
#'  \item{assess_freq}{}
#'  
#' 



#' Example data for harvest control rule functions
#'
#' A list containing data necessary to run \code{constrain} and \code{calcTAC}
#' functions. Generated with an example run of \code{recoverySim} using Fraser
#' River salmon data.
#'
#' @format A list with nine elements.
#' \describe{
#'   \item{stock}{A character vector of stock names, which are approximately
#'   equivalent to conservation units.}
#'   \item{mu}{A character vector of management unit names representing the MU
#'   that each stock belongs to.}
#'   \item{forecastMU}{A numeric vector representing examples of forecasted
#'   recruit abundance (millions of individuals). Note that
#'   although \code{length(forecastMU)} is equal to the number of stocks, the
#'   number of unique values is equal to the number of MUs because this is the
#'   scale at which test fisheries occur}
#'   \item{lowFRP}{A numeric vector representing MU-specific lower fishery
#'   reference points.}
#'   \item{highFRP}{A numeric vector representing MU-specific upper fishery
#'   reference points.}
#'   \item{minER}{A numeric vector representing MU-specific minimum exploitation
#'   rates. Intended to account for incidental harvest and test fisheries.}
#'   \item{maxER}{A numeric representing the maximum target exploitation rate
#'   allowed under the harvest control rule.}
#'   \item{ppnMixedFisher}{A numeric vector representing the proportion of the
#'   TAC that should be allocated to mixed-stock, as opposed to single-stock,
#'   fisheries.}
#' }
#' @source recoverySim model run with parameters provided by Pacific Salmon
#' Commission and Fraser River Sockeye Spawning Initiative (DFO)
"exampleHCRData"

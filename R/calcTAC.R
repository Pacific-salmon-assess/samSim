#' Calculate total allowable catch
#'
#' This function calculates total allowable catch (TAC) for different MUs using
#' either a) a fixed exploitation rate or b) a simplified version of the total
#' allowable mortality rule currently used to manage the fishery. Total
#' allowable catch is divided between one American and two Canadian fisheries
#' (mixed stock and single stock) based on \code{ppnMix} variable.
#'
#' All values should be passed as single values using apply family or for loops
#' because if statements are common.
#'
#' In the case of the TAM rule TAC is based on  abundance relative to two
#' fishery reference points. This determines whether exploitation is based
#' on a minimum exploitation rate, fixed escapement goal, or maximum
#' exploitation rate. Note that abundance relative to reference points is
#' adjusted downwards to account for anticipated en route mortality
#' (\code{manAdjustment}). TAC will be further reduced based on
#' \code{overlapConstraint} which represents whether other MUs that co-migrate
#' are at sufficiently low abundance to limit a given fishery (see
#' \code{overlapConstraint} for additional details).
#'
#' @param rec A numeric representing MU-specific return abundance.
#' @param canER A numeric representing the target Canadian exploitation rate if
#' \code{harvContRule = fixedER}
#' @param harvContRule A character signifying whether TACs are based on a fixed
#' exploitation rate (\code{fixedER}) or the simplified TAM rule (\code{TAM}).
#' @param amER A numeric representing the target American exploitation rate.
#' @param ppnMix A numeric representing the proportion of the Canadian TAC
#' allocated to mixed stock fisheries.
#' @param species A character that can currently take either \code{chum} or
#' \code{sockeye} values. Determines how American TAC is partitioned relative
#' to Canadian.
#' @param manAdjustment A numeric  representing MU-specific management
#' adjustments. These values are used to adjust forecasted spawner abundance
#' to account for en route mortality (i.e. they increase the target escapement
#' goal).
#' @param lowFRP A numeric representing a MU-specific lower fishery reference
#' point.
#' @param highFRP A numeric vector a MU-specific upper fishery reference point.
#' @param minER A numeric representing minimum exploitation rate (intended to
#' represent mortality due to bycatch or test fisheries).
#' @param maxER A numeric representing maximum exploitation rate that is
#' applied when MU is above its higher FRP.
#' @param overlapConstraint A numeric representing whether a given MU's TAC
#' should be constrained or not.
#' @return Returns a six element list of numeric vectors with length equal to
#' forecast: American TAC, single fishery TAC, mixed fishery TAC, total TAC,
#' unconstrained American TAC, and unconstrained Canadian TAC (latter two for
#' reference purposes only).
#'
#' @examples
#' #Note that the function is intended to receive vectors rather than the DF
#' #used in this example to increase efficiency within the full closed-loop
#' simulation.
#' head(exampleHCRList)
#' names(exampleHCRList)[4] <- "recRYMU"
#'
#' rec <- exampleHCRList$recRYMU
#' lowFRP <- exampleHCRList$lowFRP
#' highFRP <- exampleHCRList$highFRP
#' manAdjustment <- exampleHCRList$adjustment
#' manUnit <- exampleHCRList$mu
#' overlapConstraint <- constrain(rec, highFRP, manAdjustment,
#'                                manUnit)$muConstrained
#'
#' ## Fixed ER version
#' calcTAC(rec, canER = 0.4, harvContRule = "fixedER", amER = 0.1, ppnMix = 1,
#'         species = "sockeye")
#' ## TAM version
#' calcTAC(rec, canER, harvContRule = "TAM", amER = 0.1, ppnMix = 1,
#'         species = "sockeye", manAdjustment = manAdjustment, lowFRP = lowFRP,
#'         highFRP = highFRP,  minER = 0.1, maxER = 0.6,
#'         overlapConstraint = overlapConstraint)
#'
#' @export
calcTAC <- function(rec, canER, harvContRule, amER, ppnMixVec, species = NULL,
                    manAdjustment = NULL, lowFRP = NULL, highFRP = NULL,
                    minER = NULL, maxER = NULL, overlapConstraint = NULL,
                    constrainMix = TRUE) {
  if (species == "sockeye") {
    if (harvContRule == "fixedER") {
      totalTAC <- rec * canER
    }
    if (harvContRule == "TAM") {
      if (length(minER) == 1) { #adjust input data to appropriate vector length
        minER <- rep(minER, length.out = length(rec))
      }
      totalTAC <- rep(NA, length = length(rec))
      for (k in seq_along(rec)) {
        #if recruitment below lower RP, en route mort not accounted for; minEr used
        if (rec[k] < lowFRP[k]) {
          totalTAC[k] <- minER[k] * rec[k]
        }
        #if stock is between ref points, ERs are either scaled to adjusted
        #recruitment or minER
        if ((rec[k] > lowFRP[k]) & (rec[k] < highFRP[k])) {
          escTarget <- lowFRP[k]
          # escapement goal is adjusted up based on pMA and recruitment
          # (equivalent to esc target + MA)
          adjTarget <- escTarget * (1 + manAdjustment[k])
          calcER <- ifelse(rec[k] > adjTarget,
                           (rec[k] - adjTarget) / rec[k],
                           0)
          #if recruitment greater than adjusted target, potential TAC = diff
          #between the two (converted to er for next line)
          tacER <- ifelse(calcER > minER[k],
                          calcER,
                          minER[k])
          totalTAC[k] <- tacER * rec[k]
        }
        #if stock is above upper reference point, ERs set to max ()
        if (rec[k] > highFRP[k]) {
          #escapement target increases w/ abundance (i.e. constant ER)
          escTarget <- ((1 - maxER) * rec[k])
          adjTarget <- escTarget * (1 + manAdjustment[k])
          calcER <- ifelse(rec[k] > adjTarget,
                           (rec[k] - adjTarget) / rec[k],
                           0)
          tacER <- ifelse(calcER > minER[k],
                          calcER,
                          minER[k])
          totalTAC[k] <- tacER * rec[k]
        }
      }
    }
    #adjust total TAC to account for AFE (400k fish) when calculating us TAC
    #divide 400k allocation evenly among all CUs
    afe <- totalTAC / sum(unique(totalTAC)) * 0.4
    amTAC <- amER * (totalTAC - afe)
    canTAC <- totalTAC - amTAC
  } else {
    #sockeye is odd in that Americans get a proportion of total TAC, not total
    #escapement; all other fisheries simplify by calculating simultaneously
    amTAC <- amER * rec
    canTAC <- canER * rec
  }
  mixTAC <- canTAC * ppnMixVec
  unconMixTAC <- mixTAC
  singTAC <- (canTAC * (1 -  ppnMixVec))
  if (harvContRule == "TAM") {
    if (constrainMix == TRUE) {
      for (k in seq_along(rec)) {
        #apply overlap constraints to marine fisheries if above lower FRP
        if (rec[k] > lowFRP[k] & overlapConstraint[k] == 1) {
          mixTAC[k] <- 0.75 * mixTAC[k]
        }
      }
    }
  }
  tacList <- list(amTAC, mixTAC, singTAC, unconMixTAC)
  tacList <- lapply(tacList, function (x){ #replace NAs with 0s
    x[is.na(x)] <- 0
    return(x)
  })
  names(tacList) <- c("amTAC", "mixTAC", "singTAC", "unconMixTAC")
  return(tacList)
}

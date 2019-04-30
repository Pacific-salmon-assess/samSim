#' Calculate total allowable catch
#'
#' This function calculates total allowable catch (TAC) for different MUs using
#' either a) a fixed exploitation rate, b) a simplified version of the total
#' allowable mortality (TAM) rule currently used to manage the Fraser sockeye
#' fishery, or c) a generic HCR that is consistent with the precautionary
#' approach (genPA). TAC is divided between one American and two Canadian
#' fisheries (mixed stock and single stock) based on \code{ppnMix} variable.
#'
#' In the case of the TAM rules TAC is based on abundance relative to two
#' fishery reference points, which determines whether exploitation is based
#' on a minimum exploitation rate, fixed escapement goal, or maximum
#' exploitation rate. Note that abundance relative to reference points is
#' adjusted downwards to account for anticipated en route mortality
#' (\code{manAdjustment}). TAC will be further reduced based on
#' \code{overlapConstraint} which represents whether other MUs that co-migrate
#' are at sufficiently low abundance to limit a given fishery (see
#' \code{overlapConstraint} for additional details).
#'
#' In the case of the genPA rule TAC is also based on abundance relative to two
#' fishery reference points, but there is no management adjustment or overlap
#' constraint. When abundance is below the lower FRP min ER is applied, between
#' the two it increases linearly, and above it is max (typically Umsy).
#'
#' @param rec A numeric representing MU-specific return abundance.
#' @param canER A numeric representing the target Canadian exploitation rate if
#' \code{harvContRule = fixedER}
#' @param harvContRule A character signifying whether TACs are based on a fixed
#' exploitation rate (\code{fixedER}) or the simplified TAM rule (\code{TAM}).
#' @param amER A numeric representing the target American exploitation rate.
#' **Note** that American catch is a function of TAC not escapement for Fraser
#' sockeye and takes into account the Aboriginal Fishery Exclusion (400k fish).
#' Thus their catch/harvest rate is generally below the input parameter.
#' @param ppnMix A numeric representing the proportion of the Canadian TAC
#' allocated to mixed stock fisheries.
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
#' calcTAC(rec, canER = 0.4, harvContRule = "fixedER", amER = 0.1, ppnMix = 1)
#' ## TAM version
#' calcTAC(rec, canER, harvContRule = "TAM", amER = 0.1, ppnMix = 1,
#'         manAdjustment = manAdjustment, lowFRP = lowFRP,
#'         highFRP = highFRP,  minER = 0.1, maxER = 0.6,
#'         overlapConstraint = overlapConstraint)
#'
#' @export
calcTAC <- function(rec, canER, harvContRule, amER, ppnMixVec,
                    manAdjustment = NULL, lowFRP = NULL, highFRP = NULL,
                    minER = NULL, maxER = NULL, overlapConstraint = NULL,
                    constrainMix = TRUE) {
  if (length(minER) == 1) { #adjust input data to appropriate vector length
    minER <- rep(minER, length.out = length(rec))
  }
  if (harvContRule == "TAM") {
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
        escTarget <- ((1 - maxER[k]) * rec[k])
        adjTarget <- escTarget * (1 + manAdjustment[k])
        calcER <- ifelse(rec[k] > adjTarget,
                         (rec[k] - adjTarget) / rec[k],
                         0)
        tacER <- ifelse(calcER > minER[k],
                        calcER,
                        minER[k])
        totalTAC[k] <- tacER * rec[k]
      }
    } #end for k in seq_along
    #adjust total TAC to account for AFE (400k fish) when calculating us TAC
    #divide 400k allocation evenly among all CUs; pmin used to constrain amTAC
    #to positive values (i.e. if afe greater than cu-specific TAC)
    afe <- pmin(totalTAC, totalTAC / sum(unique(totalTAC)) * 0.4)
    amTAC <- amER * (totalTAC - afe)
    canTAC <- totalTAC - amTAC
  } else {
    #Fraser sockeye is odd in that Americans get a proportion of total TAC, not total
    #escapement; all other fisheries simplify by calculating simultaneously
    if (harvContRule == "fixedER") {
      amTAC <- amER * rec
      canTAC <- canER * rec
    }
    if (harvContRule == "genPA") {
      tacER <- rep(NA, length = length(rec))
      if (length(minER) == 1) { #adjust input data to appropriate vector length
        minER <- rep(minER, length.out = length(rec))
      }
      for (k in seq_along(rec)) {
        #if recruitment below lower RP, en route mort not accounted for; minEr used
        if (rec[k] < lowFRP[k]) {
          tacER[k] <- minER[k]
        }
        if ((rec[k] > lowFRP[k]) & (rec[k] < highFRP[k])) {
          calcER <- maxER[k] * ((rec[k] - lowFRP[k]) / (highFRP[k] - lowFRP[k]))
          tacER[k] <- ifelse(calcER > minER[k],
                          calcER,
                          minER[k])
        }
        #if stock is above upper reference point, ERs set to max ()
        if (rec[k] > highFRP[k]) {
          #escapement target increases w/ abundance (i.e. constant ER)
          tacER[k] <- maxER[k]
        }
      } #end for k in seq_along
      amTAC <- amER * rec
      canTAC <- tacER * rec
    } #end harvContRule == "genPA
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
#
#
# rec = recRYManU[y, ]; canER = canER;
# harvContRule = harvContRule; amER = amER;
# ppnMixVec = ppnMixVec;
# species = species; manAdjustment = manAdjustment;
# lowFRP = lowRefPt[y, ]; highFRP = highRefPt[y, ];
# minER = minER; maxER = maxER;
# overlapConstraint = overlapConstraint[y, ];
# constrainMix = constrainMix

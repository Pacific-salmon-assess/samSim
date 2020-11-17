#' Diagnostic plot
#'
#' This function generates diagnostic figures showing CU-specific traits from
#' one Monte Carlo trial within the recoverySimulator function. First figure is
#' a time series of spawner and recruit abundance through the observation and
#' simulation period for one CU (selected at random or with \code{focalCU}
#' argument). Second figure shows the stock-recruit relationship for that same
#' CU separating abundances from process (filled circles) and observation model
#' (closed circles). Remaining figures are time series from up to six CUs (if
#' more are passed, 6 will be selected at random) showing the metrics passed in
#' the input array.
#'
#' @param varArray A numeric array containing CU-specific data.
#' @param varNames A character vector containing the names of each variable
#' passed to the function (i.e. names for the third dimension of
#' \code{varArray}).
#' @param stkName An optional character vector containing stock/CU names.
#' Recommended for plotting.
#' @param model A character vector specifying whether each CU is fit to a
#' Ricker or Larkin model. Length equal to the number of CUs.
#' @param ricB A numeric vector of each CU's Ricker beta parameter.
#' @param larBList A numeric list of four numeric vectors, each representing
#' CU-specific estimates of lagged density dependent effects from Larkin models.
#' @param medAbundance A numeric matrix of each CU's median abundance plus upper
#' and lower quantiles.
#' @param nPrime Numeric representing length of pre-simulation priming period
#' (i.e. observation period) in the Monte Carlo trial.
#' @param extinct A binary vector of extinct (\code{1}) and extant (\code{0})
#' CUs.
#' @param focalCU A numeric identifying the CU that should be selected for
#' presentation. If \code{NULL} (default), a CU is selected at random.
#' @return Returns a multi-page panel, multi-page plot.
#'
#' @export

plotDiagCU <- function(varArray, varNames, stkName = NULL, model, ricB,
                       larBList = NULL, medAbundance = NULL, nPrime, extinct,
                       focalCU = NULL){
  nCU <- dim(varArray)[2]
  if(nCU > 6){ #subset to keep plot readable
    selectCUs <- sample(1:nCU, 6)
    selectCUs <- selectCUs[order(selectCUs)]
  } else{
    selectCUs <- c(1:nCU)
  }
  ## Modifications to allow function to work with simplified PSF model
  if (is.null(stkName)) {
    stkName <- paste("Pop", seq(from = 1, to = nCU, by = 1), sep = " ")
  }
  if (length(model) == 1) {
    model <- rep(model, length.out = nCU)
    ricB <- rep(ricB, length.out = nCU)
  }
  p <- ifelse(is.null(focalCU), sample(selectCUs, 1), focalCU)
  cuName <- stkName[p]

  # Store SR pars for plotting
  #use mean value in case alpha varies over time
  trueA <- mean(varArray[ , p, "Productivity"])
  # note that although this is "true" alpha used to forward simulate, it was
  # drawn from a distribution and may not match historical observations well
  if (model[p] == "ricker" | model[p]=="rickerSurv"){
    trueB <- ricB[p]
  }
  if (model[p] == "larkin"){
    trueB <- c(larBList$lag0[p], larBList$lag1[p], larBList$lag2[p],
               larBList$lag3[p])
  }
  #index value for last year w/ SR estimate
  lastNonNA <- max(which(!is.na(varArray[ , p, "Est Productivity"])))
  modelA <- varArray[lastNonNA, p, "Est Productivity"]
  modelB <- varArray[lastNonNA, p, "Est Beta"]

  par(mfrow = c(3, 2), mar = c(3.5, 4, 0.5, 0.5), oma = c(0, 0, 0, 0),
      cex.lab = 1)
  nYears <- dim(varArray)[1]
  colPal <- data.frame(col = viridis::viridis(length(selectCUs), begin = 0, end = 1),
                       cu = selectCUs)
  colPal$col <- as.character(colPal$col)

  # Plot single CU's abundance relative to median and 25/75 percentiles
  temp1 <- plot(1, type = "n", ylab = "", xlab = "", xlim = c(0.5, nYears),
                ylim = c(0, max(c(varArray[ , p, "Spawners"],
                                  varArray[ , p, "Recruits BY"],
                                  medAbundance[p, 3]), na.rm = TRUE)))
  mtext(side = 2, line = 2.5, "Abundance")
  lines(varArray[ , p, "Spawners"], xaxt = "n", yaxt = "n", type = "l",
        lty = 1, lwd = 1.5, col = colPal[colPal$cu == p, ]$col)
  lines(varArray[ , p, "Recruits BY"], xaxt = "n", yaxt = "n", type = "l",
        lty = 2, lwd = 1.5, col = colPal[colPal$cu == p, ]$col)
  legend("topright", legend = c("S", "Rec"), lty = c(1, 2), cex = 1.1, lwd = 1.25,
         col = colPal[colPal$cu == p, ]$col)
  abline(v = nPrime, lty = 3, lwd = 1.25)
  if (is.null(medAbundance) == FALSE) {
    abline(h = medAbundance[p,1], lwd = 1.25, lty = 1) #median line
    abline(h = medAbundance[p,2], lwd = 1, lty = 3) #25% line
    abline(h = medAbundance[p,3], lwd = 1, lty = 3) #75% line
  }
  text(0.1 * nrow(varArray[ , p, ]),
       0.95 * max(c(varArray[ , p, "Spawners"], varArray[ , p, "Recruits BY"],
                  medAbundance[p, 3]), na.rm = TRUE),
       paste(cuName), cex = 1.4)

  # Plot stock-recruit relationship (observed and predicted) first
  par(mar = c(4, 4, 0.5, 0.5))
  plot(varArray[ , p, "Recruits BY"] ~ varArray[ , p, "Spawners"],
       bg = scales::alpha("black", 0.4), pch = 21,
       ylim = c(0, max(c(varArray[, p, "Recruits BY"],
                         varArray[, p, "Obs Recruits BY"]), na.rm = TRUE)),
       xlim = c(0, max(c(varArray[, p, "Spawners"],
                         varArray[, p, "Obs Spawners"]), na.rm = TRUE)),
       xlab = "", ylab = "")
  mtext(side = 1, line = 2.5, "Spawners")
  mtext(side = 2, line = 2.5, "Recruits")
  #plot Ricker curve if only one beta, plot Larkin if more than one beta
  if (length(trueB) == 1) {
    curve(x * exp(trueA - trueB * x), from = 0,
          to = max(varArray[, p, "Spawners"], na.rm = TRUE),
          add = TRUE, lwd = 2)
  } else { #complicated because requires drawing 4 curves, each based on 4 years
    #of spawner abundance; based on last generation
    yr4 <- sample(10:30, 1) #randomly select one year where extinction unlikely
    for(i in 0:3) {
      dum <- varArray[c((yr4 - 4 - i):(yr4 - i)), p, "Spawners"]
      curve(x * exp(trueA - trueB[1] * x - (trueB[2] * dum[2]) -
                      (trueB[3] * dum[3]) - (trueB[4] * dum[4])),
            from = 0, to = max(varArray[, p, "Spawners"], na.rm = TRUE),
            add = TRUE)
    }
  }
  points(varArray[, p, "Obs Recruits BY"] ~ varArray[, p, "Obs Spawners"],
         bg = scales::alpha("white", 0.4), pch = 21)
  curve(x * exp(modelA - modelB * x), from = 0, #estimated curve
        to = max(varArray[ , p, "Obs Spawners"], na.rm = TRUE),
        add = TRUE, lwd = 2, lty = 2)
  legend("topright", legend = c("True", "Obs"), lty = c(1, 2), pch = c(21, 21),
         bg = "white", pt.bg = c(scales::alpha("black", 0.4), "white"),
         cex = 1.1, lwd = 1.25)
  text(0.1*max(c(varArray[ , p, "Spawners"], varArray[ , p, "Obs Spawners"]),
               na.rm = TRUE),
       0.95*max(c(varArray[ , p, "Recruits BY"],
                  varArray[ , p, "Obs Recruits BY"]), na.rm = TRUE),
       paste(cuName), cex = 1.4)

  # Loop across for remainder of diagnostics
  par(mar = c(3.5, 4, 0.5, 0.5))
  #skip estimated beta
  plotVars <- which(varNames %in% varNames[!varNames %in% c("Est Beta",
                                                            "Obs Spawners",
                                                            "Obs Recruits BY")])
  for (i in plotVars) {
    #do not plot if all NAs in a simulation
    if (all(is.na(varArray[, selectCUs, i]))) {
      temp1 <- plot(1, type = "n", ylab = "", xlab = "", xlim = c(0.5, nYears),
                    ylim = c(0, 1)) #empty plot
      text(nYears / 2, 0.5, "All 0s or NAs")
      mtext(side = 2, line = 2.5, print(varNames[i]))
    } else {
      temp1 <- plot(1, type="n", ylab="", xlab="", xlim = c(0.5, nYears),
                    ylim = c(min(0, min(varArray[ , selectCUs, i],
                                        na.rm = TRUE)),
                             max(varArray[ , selectCUs, i], na.rm = TRUE)))
      mtext(side = 2, line = 2.5, print(varNames[i]))
      abline(v = nPrime, lty = 3) #vert line when management actions kick in
      for(j in selectCUs){ # plot each column
        lines(varArray[, j, i], xaxt = "n", yaxt = "n", type = "l", lty = 1,
              lwd = 1.25, col = colPal[colPal$cu == j, 1])
      }
      if (i > 3) {
        #for spawner abundance plot X when CUs below threshold
        for (j in selectCUs) {
          for (h in 6:length(varArray[, j, i])) {
            if (extinct[h, j] == 1) {
              text(x = h, y = 0, labels = c("X"), cex = 2.5,
                   col = colPal[colPal$cu == j, 1])
              break
            }
          } #end for(h in 1:length(varArray[,j,i]))
        } #end for(j in 1:ncol(varArray[,selectCUs,i]))
      }
      if (i == 4 | i == 8 | i == 14 | i == 20 | i == 26) {
        legend("topleft", legend = c(stkName[selectCUs]), ncol = 3, lty = 1,
               col = colPal$col, cex = 1.1, lwd = 1.25, bg="white")
      }
    }
  } #end time series plots
} #end function


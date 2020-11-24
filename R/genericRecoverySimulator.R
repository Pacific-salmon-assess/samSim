#' Recovery Simulator - Generic Version
#'
#' Closed-loop simulation tool to assess management procedures and inform
#' Pacific salmon rebuilding strategies.
#'
#' Adapted from the original recoverySim() function developed to inform rebuilding
#' strategies for Fraser sockeye salmon and North Coast Chum (see Holt, C.A., Freshwater, C.,
#' Holt, K., and Huang, A.-M. 2020. A quantitative tool for evaluating rebuilding plans for
#' Pacific salmon. Can. Tech. Rep. Fish. Aquat. Sci. 3402: v + 26 p.;
#' https://waves-vagues.dfo-mpo.gc.ca/Library/40889385.pdf)
#'
#' This version has been adapted to easily accomodate alternative life history types.
#' Additionally, the Total Allowable Mortality (TAM) HCR specific to Fraser sockeye fisheries,
#' and assocaited performance measures, has been removed.
#'
#' Simulation runs primed with observed SR data
#' (recDat input). The model includes data generation, variation in age
#' structure, survey design, and variable exploitation rules.
#' @importFrom here here
#' @importFrom dplyr group_by summarise
#' @param y TO BE DEFINED
#' @return TO BE DEFINED
#' @export
#'
#'
genericRecoverySim <- function(simPar, cuPar, catchDat=NULL, srDat=NULL,
                               variableCU=FALSE, makeSubDirs=TRUE, ricPars,
                               larkPars=NULL, cuCustomCorrMat=NULL,
                               erCorrMat=NULL, dirName, nTrials=100, uniqueProd=TRUE,
                               uniqueSurv=FALSE, random=FALSE) {
  # If random = TRUE then each simulation will start at a different point
  # i.e. should ALWAYS be FALSE except for convenience when running independent
  # chains to test convergence
  if (random != TRUE) {
    set.seed(123)
  }

  # Silence warnings present in R 3.5.1
  options(warnPartialMatchArgs = FALSE)
  #_______________________________________________________________________
  ## Set up input arguments

  # Simulation parameters (biological, observation, and management)
  nameOM <- simPar$nameOM
  nameMP <- simPar$nameMP
  plotOrder <- simPar$plotOrder
  prod <- simPar$prodRegime #what is current prodRegime; requires list of posterior estimates for sampling
  harvContRule <- simPar$harvContRule
  bm <- simPar$benchmark
  species <- simPar$species #species (sockeye, chum, pink, coho)
  simYears <- simPar$simYears #total length of simulation period
  nTrials <- nTrials #number of trials to simulate
  canER <- simPar$canER #baseline exploitation rate divided among mixed and single CU fisheries
  ppnMix <- simPar$propMixHigh #ppn of Canadian harvest allocated to mixed stock fisheries when abundance is high (default)
  singleHCR <- ifelse(is.null(simPar$singleHCR), FALSE, simPar$singleHCR) #if TRUE single stock TAC is only harvested when BMs met
  moveTAC <- ifelse(is.null(simPar$moveTAC), FALSE, simPar$moveTAC) #if TRUE single stock TAC for low abundance CUs is re-allocated
  rho <- simPar$rho #autocorrelation coefficient in recruitment residuals
  correlCU <- simPar$correlCU #correlation among CUs in recruitment deviations
  adjSig <- simPar$adjustSig # used to scale CU specific sigma up or down
  tauCatch <- simPar$tauCatch # CU-specific catch assignment error for observation model
  obsSig <- simPar$obsSig #estimated spawner abundance error
  mixOUSig <- simPar$mixOUSig #mixed fisheries outcome uncertainty error (same for canadian and US fisheries)
  singOUSig <- simPar$singOUSig #single fisheries outcome uncertainty error; assumed to be same across CUs but could be varied
  obsMixCatchSig <- simPar$obsMixCatchSig #estimated catch error in mixed-CU fisheries (same for canadian and US fisheries)
  obsSingCatchSig <- simPar$obsSingCatchSig #estimated catch error in single CU fisheries
  ageErr <- simPar$obsAgeErr #error assigning observed recruits to a brood year
  lowCatchThresh <- simPar$lowCatchThresh #area specific-BM representing minimum catch levels
  highCatchThresh <- simPar$highCatchThresh
  agEscTarget <- ifelse(is.null(simPar$agEscTarget), NA, simPar$agEscTarget) #aggregate escapement target
  extinctThresh <- simPar$extinctThresh #minimum number of spawners before set to 0
  preFMigMort <- ifelse(is.null(simPar$preFMigMort), 1,
                        as.numeric(simPar$preFMigMort)) #proportion of en route mortality occurring before single stock fisheries

  # Should BMs be fixed at normative period?; if yes, then BMs aren't updated during sim period
  normPeriod <- ifelse(is.null(simPar$normPeriod), TRUE, simPar$normPeriod)

  ## CU-specific parameters
  cuPar$stkName <- abbreviate(cuPar$stkName, minlength = 4)
  cuPar <- with(cuPar, cuPar[order(as.numeric(stk)), ]) #temporary subset of CUs to examine
  cuNameOM <- cuPar$nameOM
  cuNameMP <- cuPar$nameMP
  nCU <- nrow(cuPar) #number of CUs in analysis
  stkName <- cuPar$stkName
  stkID <- cuPar$stk
  manUnit <- cuPar$manUnit #aggregate of CUs managed coherently
  muName <- unique(manUnit)
  nMU <- length(muName)
  model <- cuPar$model #stock recruit model type (larkin, ricker)
  domCycle <- sapply(seq_along(stkName), function (x) {
    ifelse(is.null(cuPar$domCycle[x]), NA, cuPar$domCycle[x])
  }) #which cycle line is dominant? (for CUs with Larkin model only)

  # Minimum exploitation rate applied
  minER <- cuPar$minER

  # American ER
  if (is.numeric(simPar$usER) == TRUE) {
    amER <- rep(simPar$usER, length.out = nCU)
  } else {
    amER <- cuPar$usER #American exploitation rate shared
  }

  # # En-route mortality
  # enRouteMort <- ifelse(is.null(simPar$enRouteMort), FALSE, simPar$enRouteMort)
  # if (enRouteMort == TRUE) {
  #   # ER mort rate (i.e. between marine and term. fisheries);
  #   # taken from in-river difference between estimates (post-2000)
  #   enRouteMR <- cuPar$meanDBE
  #   enRouteSig <- cuPar$sdDBE
  # } else if (enRouteMort == FALSE) {
  #   enRouteMR <- rep(0, length.out = nrow(cuPar))
  #   enRouteSig <- rep(0, length.out = nrow(cuPar))
  # } #end if is.null(enRouteMort)
  # # adjust en route mortality variation for sensitivity analysis
  # if (!is.null(simPar$adjustEnRouteSig)) {
  #   enRouteSig <- enRouteSig * simPar$adjustEnRouteSig
  # }

  # Forecast parameters
  forecastMean <- cuPar$meanForecast
  forecastSig <- if (is.null(simPar$adjustForecast)) {
    cuPar$sdForecast
  } else {
    cuPar$sdForecast * simPar$adjustForecast
  }

  # Age structure
  ageStruc <- matrix(c(cuPar$meanRec2, cuPar$meanRec3, cuPar$meanRec4, cuPar$meanRec5, cuPar$meanRec6), nrow=nCU, ncol=5) #mean proportion of each age class in returns
  nAges <- ncol(ageStruc) #total number of ages at return in ageStruc matrix (does not mean that modeled populations actually contain 4 ages at maturity)
  tauAge <- cuPar$tauCycAge * simPar$adjustAge #CU-specific variation in age-at-maturity, adjusted by scenario
  medAbundance <- cbind(cuPar$medianRec, cuPar$lowQRec, cuPar$highQRec) #matrix of long term abundances (median, lower and upper quantile)
  gen<-max(cuPar$aveGenTime)  # generation time
  ageFirstRec <- max(cuPar$ageFirstRec) # age at first recruitment
  ageMaxRec<-max(cuPar$ageMaxRec) # maximum age of recruits


  ## Stock-recruitment parameters
  ricA <- cuPar$alpha
  ricB <- cuPar$beta0
  ricSig <- cuPar$sigma
  larA <- cuPar$larkAlpha
  larB <- cuPar$larkBeta0
  larB1 <- cuPar$larkBeta1
  larB2 <- cuPar$larkBeta2
  larB3 <- cuPar$larkBeta3
  larSig <- cuPar$larkSigma

  coef1<-cuPar$coef1
  coVarInit<-cuPar$covarInit
  mu_logCoVar<-cuPar$mu_logCovar
  sig_logCoVar<-cuPar$sig_logCovar
  min_logCoVar<-cuPar$min_logCovar
  max_logCoVar<-cuPar$max_logCovar

  # Coerce all stocks to have the same alpha parameter (regardless of model
  # structure), others will vary
  if (uniqueProd == FALSE) {
    ricA <- rep(mean(cuPar$alpha), length.out = nCU)
    larA <- rep(mean(cuPar$alpha), length.out = nCU)
    ricSig <- rep(mean(cuPar$sigma), length.out = nCU)
    larSig <- rep(mean(cuPar$sigma), length.out = nCU)
  }

  # Coerce all stocks to have the same alpha parameter (regardless of model
  # structure), others will vary
  if (is.na(mu_logCoVar[1]) == FALSE & uniqueSurv == FALSE) {
    mu_logCoVar<-mean(mu_logCoVar)
    sig_logCoVar<-mean(sig_logCoVar)
    min_logCoVar<-mean(min_logCoVar)
    max_logCoVar<-mean(max_logCoVar)
  }


  # If .csv of par dist is not passed and productivity is something other than "med", give error warning
  if (prod != "med" & is.null(ricPars) == TRUE) {
    stop("Full SR parameter dataset necessary to simulate alternative
         productivity scenarios")
  }

  # If .csv of par dist is passed, change from median values to sample from par dist
  if (is.null(ricPars) == FALSE) {
    dum <- getSRPars(pars = ricPars, alphaOnly = FALSE, highP = 0.95,
                     lowP = 0.05, stks = stkID)
    ricA <- dum$pMed[["alpha"]]
    ricB <- dum$pMed[["beta"]]
    ricSig <- dum$pMed[["sigma"]]
    ricGamma <- dum$pMed[["gamma"]]
    if (is.null(larkPars) == FALSE) {
      #getSRPars provides quantiles as well which can be used for high/low prod
      #treatments but this is currently replaced by applying a simple scalar
      #based on mean differences estimated from kalman filter
      dumLark <- getSRPars(pars = larkPars, alphaOnly = TRUE, highP = 0.95,
                           lowP = 0.05, stks = stkID)
      larA <- (dumLark$pMed[["alpha"]])
      larB <- (dumLark$pMed[["beta0"]])
      larB1 <- (dumLark$pMed[["beta1"]])
      larB2 <- (dumLark$pMed[["beta2"]])
      larB3 <- (dumLark$pMed[["beta3"]])
      larSig <- (dumLark$pMed[["sigma"]])
    }
  }

  # Adjust sigma (Ricker only) following Holt and Folkes 2015 if a
  # transformation term is present and TRUE
  if (is.null(simPar$arSigTransform) == FALSE & simPar$arSigTransform == TRUE) {
    ricSig <- ricSig^2 * (1 - rho^2)
  }

  # Save reference alpha values for use when calculating BMs, then adjust alpha
  # based on productivity scenario
  refAlpha <- ifelse(model == "ricker" | model =="rickerSurv", ricA, larA)
  # Then, adjust alpha based on productivity scenario
  if (prod == "low" | prod == "lowStudT") {
    alpha <- 0.65 * refAlpha
  } else if (prod == "high") {
    alpha <- 1.35 * refAlpha
  } else {
    alpha <- refAlpha
  }

  # is the productivity scenario stable
  stable <- ifelse(prod %in% c("decline", "divergent", "divergentSmall",
                               "oneUp", "oneDown"),
                   FALSE,
                   TRUE)

  # For stable trends use as placeholder for subsequent ifelse
  finalAlpha <- alpha
  prodScalars <- rep(1, nCU)
  if (prod == "decline" ) {
    prodScalars <- rep(0.65, nCU)
  } else if (prod == "divergent") {
    prodScalars <- sample(c(0.65, 1, 1.35), nCU, replace = TRUE)
  } else if (prod == "divergentSmall") {
    prodScalars <- ifelse(cuPar$medianRec < median(cuPar$medianRec), 0.65, 1)
  } else if (prod == "oneDown") {
    drawCU <- round(runif(1, min = 0.5, max = nCU))
    prodScalars[drawCU] <- 0.65
  } else if (prod == "oneUp") {
    drawCU <- round(runif(1, min = 0.5, max = nCU))
    prodScalars[drawCU] <- 1.35
  } else if (prod == "scalar") {
    prodScalars <- rep(simPar$prodScalar, nCU)
  }

  finalAlpha <- prodScalars * alpha
  trendLength <- 3 * gen
  trendAlpha <- (finalAlpha - alpha) / trendLength
  cuProdTrends <- dplyr::case_when(
    prodScalars == "0.65" ~ "decline",
    prodScalars == "1" ~ "stable",
    prodScalars == "1.35" ~ "increase"
  )
  beta <- ifelse(model == "ricker" | model == "rickerSurv", ricB, larB)
  if (is.null(simPar$adjustBeta) == FALSE) {
    beta <- beta * simPar$adjustBeta
  }
  # Adjust sigma up or down
  sig <- ifelse(model == "ricker" | model=="rickerSurv", ricSig, larSig) * adjSig

  gamma<-ifelse(model == "rickerSurv", ricGamma, NA)

  #Add correlations in rec deviations
  if (simPar$corrMat == TRUE) { #replace uniform correlation w/ custom matrix
    if (nrow(cuCustomCorrMat) != nCU) {
      stop("Custom correlation matrix does not match number of CUs")
    }
    correlCU <- as.matrix(cuCustomCorrMat)
  }
  #calculate correlations among CUs
  sigMat <- matrix(as.numeric(sig), nrow = 1, ncol = nCU)
  #calculate shared variance and correct based on correlation
  covMat <- (t(sigMat) %*% sigMat) * correlCU
  diag(covMat) <- as.numeric(sig^2) #add variance


  ## Priming period: get stock-recruitment data
  if (!is.null(srDat)) { #transform rec data if available
    if (!is.null(catchDat)) { #add catch data if available
      #catchDat comes first so that R can be pulled regardless of DF width
      recDat <- Reduce(function(x, y) merge(x, y, by = c("stk", "yr")),
                       list(catchDat, srDat))
    } else{
      recDat <- srDat
    }
    # Abbreviate CU names
    if (!is.null(recDat$CU)) {
      recDat$CU <- abbreviate(recDat$CU, minlength = 4)
    }

    # Remove stocks from SR dataset that aren't in CU parameter inputs
    recDat <- recDat %>%
      dplyr::filter(stk %in% cuPar$stk) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(totalRec = sum(rec2, rec3, rec4, rec5, rec6)) %>%
      ungroupRowwiseDF()

    if (length(unique(recDat$stk)) != length(unique(cuPar$stk))) {
      stop("SR input dataset does not match parameter inputs")
    }

    # Trim SR data so that empty rows are excluded; otherwise gaps may be
    # retained when catch data is more up to date
    maxYears <- recDat %>%
      dplyr::group_by(stk) %>%
      dplyr::filter(!is.na(totalRec)) %>%
      dplyr::summarise(maxYr = max(yr))
    recDat <- recDat %>%
      dplyr::filter(!yr > min(maxYears$maxYr))
    recDat <- with(recDat, recDat[order(stk, yr),])
    summRec <- recDat %>%
      dplyr::group_by(stk) %>%
      dplyr::summarise(tsLength = length(ets),
                       maxRec = max(totalRec, na.rm = TRUE))
    # Define nPrime
    nPrime <- max(summRec[, "tsLength"])
    dumFull <- vector("list", nCU)
    #list of stk numbers to pass to following for loop
    stkList <- unique(cuPar$stk)
    #add NAs to front end of shorter TS to ensure all matrices are same length
    for (k in 1:nCU) {
      dum <- recDat %>%
        dplyr::filter(stk == stkList[k])
      if (nrow(dum) < nPrime) {
        empties <- nPrime - nrow(dum)
        emptyMat <- matrix(NA, nrow = empties, ncol = ncol(recDat))
        colnames(emptyMat) <- colnames(dum)
        dum <- rbind(emptyMat, dum)
      }
      if(!"singCatch" %in% colnames(dum) & "frfnCatch" %in% colnames(dum)) {
        dum <- dum %>%
          dplyr::rename(singCatch = frfnCatch)
      }
      if(!"mixCatch" %in% colnames(dum) & "marCatch" %in% colnames(dum)) {
        dum <- dum %>%
          dplyr::rename(mixCatch = marCatch)
      }
      #convert from tbbl to dataframe to silence unknown column warnings
      dumFull[[k]] <- as.data.frame(dum)
    }
    recOut <- dumFull
    #total model run length = TS priming period duration + sim length
    nYears <- nPrime + simYears
    #vector used to orient Larkin BM estimates and TAM rules; note that by
    #default cycle 1 = first year of input data
    cycle <- genCycle(min(recOut[[1]]$yr), nYears)
    residMatrix <- getResiduals(recOut, model) #pull residuals from observed data and save
    #calculate firstYr here because catch and rec data may differ in length
    firstYr <- min(sapply(recOut, function(x) min(x$yr, na.rm = TRUE)))
  } # end of (!is.null(srDat))

  # Extract proportions at age a for each CU (where, a = 2,3,4,5,6)
  ppn2 <- ageStruc[, 1] #proportion at age parameters
  ppn3 <- ageStruc[, 2]
  ppn4 <- ageStruc[, 3]
  ppn5 <- ageStruc[, 4]
  ppn6 <- ageStruc[, 5]


  # DF used to store MU-specific observation errors; can't use matrix because mix of numeric and characters; updated annually
  obsErrDat <- data.frame(mu = manUnit,
                          cu = stkName,
                          spwn = NA, #spawner obs error
                          mixC = NA, #mixed fishery obs catch errror
                          singC = NA, #single CU fishery obs catch error
                          forecast = NA #forecasted recruitment error
  )
  earlyPeriod <- (nPrime + 1):(nPrime + (gen * 3)) #defines years representing first three generations after sim starts; used for some PMs
  endEarly <- max(earlyPeriod)


  #_____________________________________________________________________
  ## Create directories (based on all scenarios in a sim run)
  dir.create(paste(here("outputs/diagnostics"), dirName, sep = "/"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(paste(here("outputs/simData"), dirName, sep = "/"),
             recursive = TRUE, showWarnings = FALSE)

  ## Create subdirectories if multiple OMs and MPs are being run
  if (makeSubDirs == TRUE) {
    subDirName <- simPar$nameOM
    dir.create(paste(here("outputs/diagnostics"), dirName, subDirName,
                     sep = "/"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(paste(here("outputs/simData"), dirName, subDirName,
                     sep = "/"),
               recursive = TRUE, showWarnings = FALSE)
  }
  #use this to generate figs/data in subsequent calls
  dirPath <- ifelse(makeSubDirs == TRUE,
                    paste(dirName, subDirName, sep = "/"),
                    dirName)
  #_____________________________________________________________________


  ## Specify key variable used for outputs
  keyVar <- switch(simPar$keyVar,
                   "prodRegime" = prod,
                   "synch" = correlCU,
                   "expRate" = ifelse(harvContRule == "fixedER", canER,
                                      harvContRule),
                   "ppnMix" = ppnMix,
                   "sigma" = adjSig,
                   "endYear" = endYr,
                   "adjustAge" = simPar$adjustAge,
                   "mixOUSig" = mixOUSig,
                   "adjustForecast" = simPar$adjustForecast,
                   "adjustEnRoute" = simPar$adjustEnRoute,
                   "obsSig" = simPar$obsSig,
                   "obsMixCatch" = simPar$obsMixCatch)
  if (is.null(keyVar)) {
    warning("Key variable misspecified; necessary for plotting")
  }



  #_______________________________________________________________________
  ## Set-up empty vectors and matrices for ALL trials

  #Population dynamics
  sEqVar <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  mSurvAge4<- matrix(NA, nrow = nYears, ncol = nTrials)

  #Harvest and observation
  expFactor <- rep(1, length.out = nTrials)

  #Management and performance
  sAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  obsSAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  recRYAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  obsRecRYAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  recBYAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  obsRecBYAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  catchAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  amCatchAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  mixCatchAg <- matrix(NA, nrow = nYears, ncol = nTrials) #sum of true Canadian mixed-CU catches (excludes CU-specific)
  obsCatchAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  estS25thBY <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  estS50thBY <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  # estS75thBY <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  estS25th <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  estS50th <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  # estS75th <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  s25th <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  s50th <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  #s75th <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  # True benchmark pars
  sMSY <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  sGen <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  # Estimated ricker pars:
  estRicA <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  estRicB <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  estYi <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  # estYi2 <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  estSlope <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  estSMSY <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  estSGen <- array(NA, dim = c(nYears, nCU, nTrials), dimnames = NULL)
  # lowerAgBM <- matrix(0, nrow = nYears, ncol = nTrials)
  # upperAgBM <- matrix(0, nrow = nYears, ncol = nTrials)
  # lowerAgObsBM <- matrix(0, nrow = nYears, ncol = nTrials)
  # upperAgObsBM <- matrix(0, nrow = nYears, ncol = nTrials)
  # catchMetric <- matrix(0, nrow = nYears, ncol = nTrials)
  # lowCatchAgBM <- matrix(0, nrow = nYears, ncol = nTrials)
  # highCatchAgBM <- matrix(0, nrow = nYears, ncol = nTrials)
  # agEscGoal <- matrix(0, nrow = nYears, ncol = nTrials) #aggregate escapement goal (not relevant for FR sockeye)
  # ppnUpperBM <- matrix(0, nrow = nYears, ncol = nTrials)
  # ppnLowerBM <- matrix(0, nrow = nYears, ncol = nTrials)
  ppnCUsUpperBM <- matrix(NA, nrow = nYears, ncol = nTrials) #proportion of CUs within a year above BM
  ppnCUsLowerBM <- matrix(NA, nrow = nYears, ncol = nTrials)
  ppnCUsUpperObsBM <- matrix(NA, nrow = nYears, ncol = nTrials)
  ppnCUsLowerObsBM <- matrix(NA, nrow = nYears, ncol = nTrials)
  ppnCUsExtinct <- matrix(NA, nrow = nYears, ncol = nTrials)
  # ppnConstrained <- matrix(NA, nrow = nYears, ncol = nTrials)
  ppnCUsOpenSingle <- matrix(NA, nrow = nYears, ncol = nTrials)
  meanSingExpRate <- matrix(NA, nrow = nYears, ncol = nTrials)
  # sGeoMean <- matrix(NA, nrow = nYears, ncol = nCU)
  # ppnSLow <- matrix(0, nrow = nTrials, ncol = nCU)
  ppnYrsCOS <- matrix(0, nrow = nTrials, ncol = nCU) #ppn of years where declining ppns are >30%
  ppnYrsWSP <- matrix(0, nrow = nTrials, ncol = nCU) #ppn of years where declining ppns are >25%
  # ppnChangeMat <- array(NA, dim = c(nYears, nCU, nTrials))
  # openFishery <- matrix(0, nrow = nYears, ncol = nMU) #should the fishery open based on species specific ref pts
  # lowRefPtMU <- matrix(NA, nrow = nYears, ncol = nMU) #ref pts that determine whether fishery should be open
  # highRefPtMU <- matrix(NA, nrow = nYears, ncol = nMU)
  # overlapConstraint <- matrix(0, nrow = nYears, ncol = nCU) #vector representing whether MUs were constrained to 75% of TAC in a given year
  targetExpRateAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  expRateAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  obsExpRateAg <- matrix(NA, nrow = nYears, ncol = nTrials)
  spwnrArray <- array(NA, dim = c(nYears, nCU, nTrials))
  recArray <- array(NA, dim = c(nYears, nCU, nTrials))
  alphaArray <- array(NA, dim = c(nYears, nCU, nTrials))
  returnArray <- array(NA, dim = c(nYears, nCU, nTrials))
  logRSArray <- array(NA, dim = c(nYears, nCU, nTrials))
  recDevArray <- array(NA, dim = c(nYears, nCU, nTrials))
  # migMortArray <- array(NA, dim = c(nYears, nCU, nTrials))
  singCatchArray <- array(NA, dim = c(nYears, nCU, nTrials))
  singTACArray <- array(NA, dim = c(nYears, nCU, nTrials))
  totalCatchArray <- array(NA, dim = c(nYears, nCU, nTrials))

  #Plotting matrices and vectors
  # hcr <- matrix(NA, nrow = nTrials, ncol = nCU)
  targetER <- matrix(NA, nrow = nTrials, ncol = nCU)
  medS <- matrix(NA, nrow = nTrials, ncol = nCU)
  varS <- matrix(NA, nrow = nTrials, ncol = nCU)
  medObsS <- matrix(NA, nrow = nTrials, ncol = nCU)
  varObsS <- matrix(NA, nrow = nTrials, ncol = nCU)
  medRecBY <- matrix(NA, nrow = nTrials, ncol = nCU)
  varRecBY <- matrix(NA, nrow = nTrials, ncol = nCU)
  medRecRY <- matrix(NA, nrow = nTrials, ncol = nCU)
  varRecRY <- matrix(NA, nrow = nTrials, ncol = nCU)
  medObsRecRY <- matrix(NA, nrow = nTrials, ncol = nCU)
  varObsRecRY <- matrix(NA, nrow = nTrials, ncol = nCU)
  medObsRecBY <- matrix(NA, nrow = nTrials, ncol = nCU)
  varObsRecBY <- matrix(NA, nrow = nTrials, ncol = nCU)
  medAlpha <- matrix(NA, nrow = nTrials, ncol = nCU)
  varAlpha <- matrix(NA, nrow = nTrials, ncol = nCU)
  medEstAlpha <- matrix(NA, nrow = nTrials, ncol = nCU)
  varEstAlpha <- matrix(NA, nrow = nTrials, ncol = nCU)
  medBeta <- matrix(NA, nrow = nTrials, ncol = nCU)
  medTotalCatchEarly <- matrix(NA, nrow = nTrials, ncol = nCU)
  medTotalCatch <- matrix(NA, nrow = nTrials, ncol = nCU)
  varTotalCatch <- matrix(NA, nrow = nTrials, ncol = nCU)
  stblTotalCatch <- matrix(NA, nrow = nTrials, ncol = nCU)
  stblObsTotalCatch <- matrix(NA, nrow = nTrials, ncol = nCU)
  medObsTotalCatch <- matrix(NA, nrow = nTrials, ncol = nCU)
  varObsTotalCatch <- matrix(NA, nrow = nTrials, ncol = nCU)
  medTotalER <- matrix(NA, nrow = nTrials, ncol = nCU)
  medTotalObsER <- matrix(NA, nrow = nTrials, ncol = nCU)
  # medTAMSingER <- matrix(NA, nrow = nTrials, ncol = nCU)
  # medForgoneCatch <- matrix(NA, nrow = nTrials, ncol = nCU) #diff between constrained and unconstrained TACs
  ppnYrsUpperBM <- matrix(NA, nrow = nTrials, ncol = nCU) #proportion of years where a CU above BM
  ppnYrsLowerBM <- matrix(NA, nrow = nTrials, ncol = nCU)
  ppnYrsUpperObsBM <- matrix(NA, nrow = nTrials, ncol = nCU)
  ppnYrsLowerObsBM <- matrix(NA, nrow = nTrials, ncol = nCU)
  ppnYrsOpenSingle <- matrix(NA, nrow = nTrials, ncol = nCU)
  counterEarlyUpperBM <- matrix(0, nrow = nTrials, ncol = nCU)
  counterEarlyLowerBM <- matrix(0, nrow = nTrials, ncol = nCU)
  counterLateUpperBM <- matrix(0, nrow = nTrials, ncol = nCU)
  counterLateLowerBM <- matrix(0, nrow = nTrials, ncol = nCU)
  counterLateUpperObsBM <- matrix(0, nrow = nTrials, ncol = nCU)
  counterLateLowerObsBM <- matrix(0, nrow = nTrials, ncol = nCU)
  medEarlyS <- matrix(NA, nrow = nTrials, ncol = nCU)
  medEarlyRecRY <- matrix(NA, nrow = nTrials, ncol = nCU)
  medEarlyTotalCatch <- matrix(NA, nrow = nTrials, ncol = nCU)


  # Randomly selects trial to draw and plot
  drawTrial <- round(runif(1, min = 0.5, max = nTrials))


  #_______________________________________________________________________
  ## Simulation model
  for (n in 1:nTrials) {

    #_____________________________________________________________________
    # Set up empty vectors and matrices for each MC trial
    ## Population dynamics
    S <- matrix(NA, nrow = nYears, ncol = nCU)
    alphaMat <- matrix(NA, nrow = nYears, ncol = nCU)
    alphaPrimeMat <- matrix(NA, nrow = nYears, ncol = nCU) # variable used to estimate sEq, sGen and sMsy when spawners generated w/ larkin
    ppnAges <- array(NA, dim=c(nYears, nCU, nAges), dimnames=NULL)
    ppnCatches <- matrix(NA, nrow = nYears, ncol = nCU)
    recBY <- matrix(NA, nrow = nYears, ncol = nCU) #recruits by brood year
    recRY <- matrix(NA, nrow = nYears, ncol = nCU) #recruits by return year
    recRY2 <- matrix(NA, nrow = nYears, ncol = nCU)
    recRY3 <- matrix(NA, nrow = nYears, ncol = nCU)
    recRY4 <- matrix(NA, nrow = nYears, ncol = nCU)
    recRY5 <- matrix(NA, nrow = nYears, ncol = nCU)
    recRY6 <- matrix(NA, nrow = nYears, ncol = nCU)
    logRS <- matrix(NA, nrow = nYears, ncol = nCU) #realized productivity (i.e. recruits/S)
    extinct <- matrix(0, nrow = nYears, ncol = nCU)
    extinctAg <- rep(0, nTrials)
    errorCU <- matrix(NA, nrow = nYears, ncol = nCU)
    migMortRate <- matrix(NA, nrow = nYears, ncol = nCU)
    laggedError <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge2Ret5 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge2Ret4 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge2Ret3 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge2Ret2 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge2Ret1 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge3Ret5 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge3Ret4 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge3Ret3 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge3Ret2 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge3Ret1 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge4Ret5 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge4Ret4 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge4Ret3 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge4Ret2 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge4Ret1 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge5Ret5 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge5Ret4 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge5Ret3 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge5Ret2 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge5Ret1 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge6Ret5 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge6Ret4 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge6Ret3 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge6Ret2 <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnAge6Ret1 <- matrix(NA, nrow = nYears, ncol = nCU)

    ## Observation and harvest
    obsS <- matrix(NA, nrow = nYears, ncol = nCU) #observed spawners
    obsLogRS <- matrix(NA, nrow = nYears, ncol = nCU)
    obsRecRY <- matrix(NA, nrow = nYears, ncol = nCU) #observed recruits by return year
    obsRecBY <- matrix(NA, nrow = nYears, ncol = nCU) #observed recruits by brood year; used as a proxy for forecasts
    obsRecBY_noAgeErr <- matrix(NA, nrow = nYears, ncol = nCU)
    ppnObsSRet5 <- array(NA, dim=c(nYears, nCU, nAges), dimnames=NULL)
    ppnObsSRet4 <- array(NA, dim=c(nYears, nCU, nAges), dimnames=NULL)
    ppnObsSRet3 <- array(NA, dim=c(nYears, nCU, nAges), dimnames=NULL)
    ppnObsSRet2 <- array(NA, dim=c(nYears, nCU, nAges), dimnames=NULL)
    ppnObsSRet1 <- array(NA, dim=c(nYears, nCU, nAges), dimnames=NULL)
    randAges <- matrix(NA, nrow=nYears, ncol=nAges)
    amTAC <- matrix(NA, nrow = nYears, ncol = nCU)
    mixTAC <- matrix(NA, nrow = nYears, ncol = nCU)
    unconMixTAC <- matrix(NA, nrow = nYears, ncol = nCU)
    singTAC <- matrix(NA, nrow = nYears, ncol = nCU)
    canTAC <- matrix(NA, nrow = nYears, ncol = nCU)
    totalTAC <- matrix(NA, nrow = nYears, ncol = nCU)
    amCatch <- matrix(NA, nrow = nYears, ncol = nCU)
    mixCatch <- matrix(NA, nrow = nYears, ncol = nCU)
    migMort <- matrix(NA, nrow = nYears, ncol = nCU)
    singCatch <- matrix(NA, nrow = nYears, ncol = nCU)
    totalCatch <- matrix(NA, nrow = nYears, ncol = nCU)
    forgoneCatch <- matrix(NA, nrow = nYears, ncol = nCU)
    amExpRate <- matrix(NA, nrow = nYears, ncol = nCU)
    mixExpRate <- matrix(NA, nrow = nYears, ncol = nCU)
    singExpRate <- matrix(NA, nrow = nYears, ncol = nCU)
    singCUStatus <- matrix(NA, nrow = nYears, ncol = nCU) #compared to BBs to determine whether singleTAC is taken (unless singleHCR = false)
    obsAmCatch <- matrix(NA, nrow = nYears, ncol = nCU)
    obsMixCatch <- matrix(NA, nrow = nYears, ncol = nCU)
    obsMigMort <- matrix(NA, nrow = nYears, ncol = nCU)
    obsSingCatch <- matrix(NA, nrow = nYears, ncol = nCU)
    obsTotalCatch <- matrix(NA, nrow = nYears, ncol = nCU)
    expRate <- matrix(NA, nrow = nYears, ncol = nCU)
    obsExpRate <- matrix(NA, nrow = nYears, ncol = nCU)
    extYears <- rep(nYears, nCU) #years extinct
    extYearsAg <- nYears
    cycleSGen <- matrix(NA, nrow = nYears, ncol = nCU) #for storing cycle line specific benchmarks in Larkin model
    cycleSMSY <- matrix(NA, nrow = nYears, ncol = nCU)
    upperBM <- matrix(0, nrow = nYears, ncol = nCU)
    lowerBM <- matrix(0, nrow = nYears, ncol = nCU)
    upperObsBM <- matrix(0, nrow = nYears, ncol = nCU)
    lowerObsBM <- matrix(0, nrow = nYears, ncol = nCU)
    counterSingleBMLow <- matrix(0, nrow = nYears, ncol = nCU) #should single stock TAC be fished given secondary HCR
    counterSingleBMHigh <- matrix(0, nrow = nYears, ncol = nCU) #should single stock TAC be re-assigned given secondary HCR
    counterUpperBM <- matrix(0, nrow = nYears, ncol = nCU)
    counterLowerBM <- matrix(0, nrow = nYears, ncol = nCU)
    counterUpperObsBM <- matrix(0, nrow = nYears, ncol = nCU)
    counterLowerObsBM <- matrix(0, nrow = nYears, ncol = nCU)
    # Management
    foreRecRY <- matrix(NA, nrow = nYears, ncol = nCU)
    recRYManU <- matrix(NA, nrow = nYears, ncol = nCU) #recruits by return year summed across MU (ncol = nCU)
    foreRecRYManU <- matrix(NA, nrow = nYears, ncol = nCU)
    # lowRefPt <- matrix(NA, nrow = nYears, ncol = nCU)
    # highRefPt <- matrix(NA, nrow = nYears, ncol = nCU)
    adjForeRec <-  matrix(NA, nrow = nYears, ncol = nCU)
    targetCanER <- matrix(NA, nrow = nYears, ncol = nCU)
    # tamSingER <- matrix(NA, nrow = nYears, ncol = nCU)
    foreRecErr <- matrix(NA, nrow = nYears, ncol = nCU)
    # # Fall-back matrices for diagnostics
    fb1 <- matrix(0, nrow = nYears, ncol = nCU) # fall back matrix for when observed BMs can't be estimated because ricB value > 4*obsS
    fb2 <- matrix(0, nrow = nYears, ncol = nCU) # fall back matrix for when true BMs can't be estimated because alpha prime is negative (only relevant for Larkin CUs)
    fb3 <- matrix(0, nrow = nYears, ncol = nCU) # fall back matrix for when estimated lower BM > higher
    fb4 <- matrix(0, nrow = nYears, ncol = nCU) # flicks on when allocation switches to single CU fishery
    fb5 <- matrix(0, nrow = nYears, ncol = nCU) # fall back matrix for when true lower BM > higher
    fb6 <- matrix(0, nrow = nYears, ncol = nCU) #fall back for when both BMs are NA

    #__________________________________________________________________________________________
    ### LOOP 1: Priming loop
    # - Includes only past data, used to represent both real and observed abundances to "prime" the simulation
    # -	Loops over observed years: 1:nPrime
    # -	Reconstructs return year recruits for historic time series (and catch if applicable)
    # -	Calculate BMs during last 2 generations prior to nPrime (neededy to prime single CU fishery and estimate Larkin BMs which depend on dom cycle line)
    #__________________________________________________________________________________________________

    for (y in 1:nPrime) {

      ## Population model: store SR pars, spawner, and recruit abundances
      alphaMat[y, ] <- refAlpha

      for (k in 1:nCU) {
        S[y, k] <- recOut[[k]]$ets[y]
        #calculate total recruitment as sum of all age classes
        recBY[y, k] <- sum(recOut[[k]][y, c("rec2", "rec3", "rec4", "rec5",
                                            "rec6")])
        #each age is a matrix, columns CUs, rows years
        ppnAges[y, k, ] <- as.matrix(recOut[[k]][y, c("rec2", "rec3", "rec4",
                                                      "rec5", "rec6")] /
                                       recBY[y, k])
        #if ppns can't be estimated due to TS gaps replace with mean values
        for (j in 1:nAges) {
          ppnAges[y, k, j] <- ifelse(is.na(ppnAges[y, k, j]),
                                     ageStruc[k, j],
                                     ppnAges[y, k, j])
        } # end of nAges loop

        if (is.null(catchDat) == FALSE) {
          amCatch[y, k] <- ifelse(is.null(recOut[[k]]$amCatch[y]),
                                  0, #not available for Fraser
                                  recOut[[k]]$amCatch[y])
          mixCatch[y, k] <- ifelse(is.null(recOut[[k]]$mixCatch[y]),
                                   0, #not available for Fraser
                                   recOut[[k]]$mixCatch[y])
          singCatch[y, k] <- ifelse(is.null(recOut[[k]]$singCatch[y]),
                                    0,
                                    recOut[[k]]$singCatch[y])
        } else {
          amCatch[y , k] <- 0
          mixCatch[y, k] <- 0
          singCatch[y, k] <- 0
        }
        expRate[y, k] <- ifelse(is.null(catchDat)==FALSE, recOut[[k]]$totalER[y], 0)
        logRS[y, k] <- log(recBY[y, k] / S[y, k])

      } # end of CU loop


      totalCatch[y, ] <- amCatch[y, ] + mixCatch[y, ] + singCatch[y, ]

      # Aggregated spawners(S), recruitment by BY (recBY), and catches over all CUs for year y of trial n
      sAg[y, n] <- sum(S[y, ])
      recBYAg[y, n] <- sum(recBY[y, ])
      mixCatchAg[y, n] <- sum(mixCatch[y, ])
      amCatchAg[y, n] <- sum(amCatch[y, ])
      catchAg[y, n] <- sum(mixCatch[y, ], singCatch[y, ], amCatch[y, ])

      # Aggregated recruitment by return year (recRY) over all CUs in year y of trial n
      if (y > 6) { # note: 6 is used because max number of age classes is 6
        recRY2[y, ] <- recBY[y - 2, ] * ppnAges[y - 2, , 1]
        recRY3[y, ] <- recBY[y - 3, ] * ppnAges[y - 3, , 2]
        recRY4[y, ] <- recBY[y - 4, ] * ppnAges[y - 4, , 3]
        recRY5[y, ] <- recBY[y - 5, ] * ppnAges[y - 5, , 4]
        recRY6[y, ] <- recBY[y - 6, ] * ppnAges[y - 6, , 5]
        recRY[y, ]<- recRY2[y, ] + recRY3[y, ] + recRY4[y, ] + recRY5[y, ] + recRY6[y, ]
      }
      recRYAg[y, n] <- sum(recRY[y, ], na.rm = TRUE)

      ## Calculate benchmark submodel: calculate BMs during last 2 gen
      # necessary to estimate  to prime single CU fishery and Larkin BMs which
      # depend on dom cycle line
      # (note that DL CUs will still be at 0, realistic for precautionary app)


      ## Management and assessment submodels: calculate BMs during last 2 gen
      # necessary to estimate  to prime single CU fishery and Larkin BMs which
      # depend on dom cycle line
      # (note that DL CUs will still be at 0, realistic for precautionary app)
      if (y > (nPrime - 2 * gen)) {
        for (k in 1:nCU) {
          # calculate percentile BMs
          temp <- S[1:y, k]
          sNoNA <- temp[!is.na(temp)]
          n25th <- round(length(sNoNA) * 0.25, 0)
          n50th <- round(length(sNoNA) * 0.50, 0)
          n75th <- round(length(sNoNA) * 0.75, 0)
          s25th[y, k, n] <- sort(sNoNA)[n25th]
          s50th[y, k, n] <- sort(sNoNA)[n50th]
          #s75th[y, k, n] <- sort(sNoNA)[n75th]
          #calculate SR BMs
          if (model[k] == "ricker") {
            sEqVar[y, k, n] <- refAlpha[k] / beta[k]
            sMSY[y, k, n] <- (1 - gsl::lambert_W0(exp(1 - refAlpha[k]))) /
              beta[k]
            sGen[y, k, n] <- as.numeric(sGenSolver(
              theta = c(refAlpha[k], refAlpha[k] / sEqVar[y, k, n], ricSig[k]),
              sMSY = sMSY[y, k, n]
            ))
          }
          if (model[k] == "rickerSurv") {
            refAlpha_prime<- refAlpha[k] + (gamma[k]*log(coVarInit[k]))
            sEqVar[y, k, n] <- refAlpha_prime / beta[k]
            sMSY[y, k, n] <- (1 - gsl::lambert_W0(exp(1 - refAlpha_prime))) /
              beta[k]
            sGen[y, k, n] <- as.numeric(sGenSolver(
              theta = c(refAlpha_prime, refAlpha_prime / sEqVar[y, k, n], ricSig[k]),
              sMSY = sMSY[y, k, n]
            ))
          }
          if (model[k] == "larkin") {
            #modified alpha used to estimate Larkin BMs
            alphaPrimeMat[y, k] <- refAlpha[k] - (larB1[k] * S[y - 1, k]) -
              (larB2[k] * S[y-2, k]) - (larB3[k] * S[y - 3,  k])
            sEqVar[y, k, n] <- ifelse(alphaPrimeMat[y, k] > 0,
                                      alphaPrimeMat[y, k] / beta[k],
                                      NA)
            cycleSMSY[y, k] <- ifelse(alphaPrimeMat[y, k] > 0,
                                      ((1 - gsl::lambert_W0(exp(
                                        1 - alphaPrimeMat[y, k]))) / beta[k]),
                                      NA)
            # if (is.na(log(cycleSMSY[y, k]))) {
            #   stop("Benchmark calculation is NA")
            # }
            cycleSGen[y, k] <- ifelse(alphaPrimeMat[y, k] > 0,
                                      as.numeric(sGenSolver(
                                        theta = c(alphaPrimeMat[y, k],
                                                  alphaPrimeMat[y, k] /
                                                    sEqVar[y, k, n],
                                                  larSig[k]),
                                        sMSY = cycleSMSY[y, k])),
                                      NA)
            #calculate annual benchmarks as medians within cycle line
            sMSY[y, k, n] <- median(cycleSMSY[seq(cycle[y], y, 4), k],
                                    na.rm = TRUE)
            sGen[y, k, n] <- median(cycleSGen[seq(cycle[y], y, 4), k],
                                    na.rm = TRUE)
          }
          if (is.na(sGen[y, k, n] & sMSY[y, k, n]) == FALSE) {
            if (sGen[y, k, n] > sMSY[y, k, n]) {
              warning("True lower benchmark greater than upper benchmark; set
                      to NA")
              sMSY[y, k, n] <- NA
              sGen[y, k, n] <- NA
              fb5[y, k] <- 1
            }
          } else {
            fb6[y, k] <- 1
          }
        }

        for (k in 1:nCU) {
          if (model[k] == "ricker" | model[k] == "rickerSurv" |
              model[k] == "larkin" & cycle[y] == domCycle[k]) {
            if (bm == "stockRecruit") {
              upperBM[y, k] <- ifelse(!is.na(sMSY[y, k, n]),
                                      0.8 * sMSY[y, k, n],
                                      0)
              lowerBM[y, k] <- ifelse(!is.na(sGen[y, k, n]),
                                      sGen[y, k, n],
                                      0)
            }
            if (bm == "percentile") {
              upperBM[y, k] <- s50th[y, k, n]
              lowerBM[y, k] <- s25th[y, k, n]
            }
          }
          # only save status for Larkin stocks when on dom cycle line otherwise
          # use previous status
          if (model[k] == "larkin" & cycle[y] != domCycle[k]) {
            upperBM[y, k] <- upperBM[y - 1, k]
            lowerBM[y, k] <- lowerBM[y - 1, k]
          } #end if (model[k] == "larkin" & cycle[y] != domCycle[k])
          if (!is.na(S[y, k])) {
            if (!is.na(upperBM[y, k]) & S[y, k] > upperBM[y, k]) {
              counterUpperBM[y, k] <- 1 #is spawner greater than upper BM
            }
            if (!is.na(lowerBM[y, k]) & S[y, k] > lowerBM[y, k]) {
              counterLowerBM[y, k] <- 1 #is spawner greater than lower BM
            }
          }#end if(!is.na(S[y, k]))
        }#end for (k in 1:nCU)
      }#end if (y > (nPrime - 2 * gen))


      ## Observation submodel: to prime simulation assume that observed are
      #equal to true
      obsS[y, ] <- S[y, ]
      obsRecBY[y, ] <- recBY[y, ]
      obsRecBYAg[y, n] <- sum(obsRecBY[y, ])
      obsRecRY[y, ] <- recRY[y, ]
      obsRecRYAg[y, n] <- recRYAg[y, n]
      obsLogRS[y, ] <- log(obsRecBY[y, ] / obsS[y, ])
      obsMixCatch[y, ] <- mixCatch[y, ]
      obsSingCatch[y, ] <- singCatch[y, ]
      obsCatchAg[y, n] <- catchAg[y, n]
      obsExpRate[y, ] <- expRate[y, ]
      obsSAg[y, n] <- sAg[y, n]
      estSMSY[y, , n] <- sMSY[y, , n]
      estSGen[y, , n] <- sGen[y, , n]
      upperObsBM[y, ] <- upperBM[y, ] #obs = true during priming
      lowerObsBM[y, ] <- lowerBM[y, ]

      } # end of year loop 1 (y = 1:nPrime)

    #__________________________________________________________________________________________________
    ### Loop 2: Infilling of last 12 years
    # -	Loops over years: (nPrime - 12):nPrime
    # -	Note: 12 years chosen because it allows default 6-year age structure to occur twice
    # -	Infills missing BY recruitment and spawner abundance for these years
    # -	For last 6 years (one full 6-year age structure): reconstructs RY recruitment-at-age using infilled data series
    #__________________________________________________________________________________________________

    # Pull residuals from observed data to make time series coherent
    errorCU[1:nPrime, ] <- residMatrix
    # Only necessary to infill values in gappy time series
    infillRecBY <- infill(recBY[1:nPrime, ])
    infillS <- infill(S[1:nPrime, ])

    #Default recruitment cap reflecting observed abundance (not quantiles)
    recCap <- 2 * apply(recBY[1:nPrime, ], 2, function(x) max(x, na.rm = TRUE))


    for (y in (nPrime - 12):nPrime) {  # Note that BMs and aggregate PMs are not recalculated after interpolation
      for (k in 1:nCU) {
        if (is.na(recBY[y, k])) {
          recBY[y, k] <- infillRecBY[y, k]
        }
        if (is.na(S[y, k])) {
          S[y, k] <- infillS[y, k]
        }
        logRS[y, k] <- log(recBY[y, k] / S[y, k])
      } #end for k
      if (y > (nPrime - 6)) {
        recRY[y, ] <- recBY[y - 2, ] * ppnAges[y - 2, , 1] +
          recBY[y - 3, ] * ppnAges[y - 3,  , 2] + recBY[y - 4, ] *
          ppnAges[y - 4, , 3] + recBY[y - 5, ] * ppnAges[y - 5, , 4] +
          recBY[y - 6, ] * ppnAges[y - 6, , 5]
        recRY2[y, ] <- recBY[y - 2, ] * ppnAges[y - 2, , 1]
        recRY3[y, ] <- recBY[y - 3, ] * ppnAges[y - 3, , 2]
        recRY4[y, ] <- recBY[y - 4, ] * ppnAges[y - 4, , 3]
        recRY5[y, ] <- recBY[y - 5, ] * ppnAges[y - 5, , 4]
        recRY6[y, ] <- recBY[y - 6, ] * ppnAges[y - 6, , 5]
      } #end if y > nPrime-6
      obsS[y, ] <- S[y, ]
      obsRecBY[y, ] <- recBY[y, ]
      obsRecRY[y, ] <- recRY[y, ]
      obsLogRS[y, ] <- log(obsRecBY[y, ] / obsS[y, ])
    } #end loop 2

    #prime AR error
    laggedError[y, ] <- log(recBY[y, ] / S[y, ]) -
      (alphaMat[y, ] - beta * S[y, ])



    #_____________________________________________________________________
    ### LOOP 3: Forward simulation
    # -	Loops over years: (nPrime + 1):nYears
    #___________________________________________________________________


    for (y in (nPrime + 1):nYears) {
      #________________________________________________________________________
      ### Population dynamics submodel

      # Specify alpha
      #In first year, switch from reference alpha used in priming to testing alpha; add trend for 3 generations by default
      if (y > (nPrime + 1)) {
        if (stable == FALSE & y < (nPrime + trendLength + 1)) {
          alphaMat[y, ] <- alphaMat[y - 1, ] + trendAlpha
        } else {
          alphaMat[y, ] <- finalAlpha
        } #end if stable == FALSE and inside trendPeriod
      } else {
        alphaMat[y, ] <- alphaMat[y - 1, ]
      } #end if y > (nPrime + 1)

      #Estimate BMs if normative period not being used, otherwise assume they are equal to last year of observation
      for (k in 1:nCU) {
        if (model[k] == "ricker") {
          if (normPeriod == TRUE) {
            sMSY[y, k, n] <- sMSY[nPrime, k, n]
            sGen[y, k, n] <- sGen[nPrime, k, n]
          } else if (normPeriod == FALSE) {
            sEqVar[y, k, n] <- refAlpha[k] / beta[k]
            sMSY[y, k, n] <- (1 - gsl::lambert_W0(exp(1 - refAlpha[k]))) / beta[k]
            sGen[y, k, n] <- as.numeric(sGenSolver(
              theta = c(refAlpha[k], refAlpha[k] / sEqVar[y, k, n], ricSig[k]),
              sMSY = sMSY[y, k, n]
            ))
          }
        } #end if model == ricker
        if (model[k] == "rickerSurv") {
          if (normPeriod == TRUE) {
            sMSY[y, k, n] <- sMSY[nPrime, k, n]
            sGen[y, k, n] <- sGen[nPrime, k, n]
          } else if (normPeriod == FALSE) {
            refAlpha_prime<- refAlpha[k] + (gamma[k]*log(coVarInit[k]))
            sEqVar[y, k, n] <- refAlpha_prime / beta[k]
            sMSY[y, k, n] <- (1 - gsl::lambert_W0(exp(1 - refAlpha_prime))) / beta[k]
            sGen[y, k, n] <- as.numeric(sGenSolver(
              theta = c(refAlpha_prime, refAlpha_prime / sEqVar[y, k, n], ricSig[k]),
              sMSY = sMSY[y, k, n]
            ))
          }
        } #end if model == rickerSurv
        if (model[k] == "larkin") {
          if (normPeriod == TRUE) {
            #calculate last observed year on same cycle line as current so that
            #4 different normative BMs are used for Larkin stocks
            larkinBMYear <- max(seq(cycle[y], nPrime, 4))
            sMSY[y, k, n] <- sMSY[larkinBMYear, k, n]
            sGen[y, k, n] <- sGen[larkinBMYear, k, n]
          } else if (normPeriod == FALSE) {
            #modified alpha used to estimate Larkin BMs
            alphaPrimeMat[y, k] <- refAlpha[k] - (larB1[k] * S[y - 1, k]) -
              (larB2[k] * S[y-2, k]) - (larB3[k] * S[y - 3, k])
            sEqVar[y, k, n] <- ifelse(alphaPrimeMat[y, k] > 0,
                                      alphaPrimeMat[y, k] / beta[k],
                                      NA)
            cycleSMSY[y, k] <- ifelse(alphaPrimeMat[y, k] > 0,
                                      (1 - gsl::lambert_W0(exp(
                                        1 - alphaPrimeMat[y, k]))) / beta[k],
                                      NA)
            cycleSGen[y, k] <- ifelse(alphaPrimeMat[y, k] > 0,
                                      as.numeric(sGenSolver(
                                        theta = c(alphaPrimeMat[y, k],
                                                  alphaPrimeMat[y, k] /
                                                    sEqVar[y, k, n],
                                                  larSig[k]),
                                        sMSY = cycleSMSY[y, k])),
                                      NA)
            #calculate annual benchmarks as medians within cycle line
            sMSY[y, k, n] <- median(cycleSMSY[seq(cycle[y], y, 4), k],
                                    na.rm = TRUE)
            sGen[y, k, n] <- median(cycleSGen[seq(cycle[y], y, 4), k],
                                    na.rm = TRUE)
          } #end if normPeriod == FALSE
        }#end if model == Larkin

        if (is.na(sGen[y, k, n] & sMSY[y, k, n]) == FALSE) {
          if (sGen[y, k, n] > sMSY[y, k, n]) {
            warning("True lower benchmark greater than upper benchmark;
                    set to NA")
            sMSY[y, k, n] <- NA
            sGen[y, k, n] <- NA
          }
        } else {
          warning("Neither benchmark could be estimated")
        }

      } #end for k in 1:nCU

      # Calculate recruitment by return year
      recRY2[y, ] <- recBY[y - 2, ] * ppnAges[y - 2, , 1]
      recRY3[y, ] <- recBY[y - 3, ] * ppnAges[y - 3, , 2]
      recRY4[y, ] <- recBY[y - 4, ] * ppnAges[y - 4, , 3]
      recRY5[y, ] <- recBY[y - 5, ] * ppnAges[y - 5, , 4]
      recRY6[y, ] <- recBY[y - 6, ] * ppnAges[y - 6, , 5]
      recRY[y, ] <- recRY2[y, ] + recRY3[y, ] + recRY4[y, ] + recRY5[y, ] + recRY6[y, ]


      recRYAg[y, n] <- sum(recRY[y, ])

      #________________________________________________________________________
      ### Observation submodel 1 (previous years abundance)

      #Calculate observedRecBY from observedRecRY, i.e. make giant brood table...
      recRY[recRY == 0] <- 0.1 * extinctThresh #necessary to avoid NAs from 0/0

      # e.g. 5 years ago, the ppn of age-2 returns in the R.RY
      ppnAge2Ret5[y - 5, ] <- recRY2[y - 5, ] / recRY[y - 5, ]
      ppnAge2Ret4[y - 4, ] <- recRY2[y - 4, ] / recRY[y - 4, ]
      ppnAge2Ret3[y - 3, ] <- recRY2[y - 3, ] / recRY[y - 3, ]
      ppnAge2Ret2[y - 2, ] <- recRY2[y - 2, ] / recRY[y - 2, ]
      ppnAge2Ret1[y - 1, ] <- recRY2[y - 1, ] / recRY[y - 1, ]

      ppnAge3Ret5[y - 5, ] <- recRY3[y - 5, ] / recRY[y - 5, ]
      ppnAge3Ret4[y - 4, ] <- recRY3[y - 4, ] / recRY[y - 4, ]
      ppnAge3Ret3[y - 3, ] <- recRY3[y - 3, ] / recRY[y - 3, ]
      ppnAge3Ret2[y - 2, ] <- recRY3[y - 2, ] / recRY[y - 2, ]
      ppnAge3Ret1[y - 1, ] <- recRY3[y - 1, ] / recRY[y - 1, ]

      ppnAge4Ret5[y - 5, ] <- recRY4[y - 5, ] / recRY[y - 5, ]
      ppnAge4Ret4[y - 4, ] <- recRY4[y - 4, ] / recRY[y - 4, ]
      ppnAge4Ret3[y - 3, ] <- recRY4[y - 3, ] / recRY[y - 3, ]
      ppnAge4Ret2[y - 2, ] <- recRY4[y - 2, ] / recRY[y - 2, ]
      ppnAge4Ret1[y - 1, ] <- recRY4[y - 1, ] / recRY[y - 1, ]

      ppnAge5Ret5[y - 5, ] <- recRY5[y - 5, ] / recRY[y - 5, ]
      ppnAge5Ret4[y - 4, ] <- recRY5[y - 4, ] / recRY[y - 4, ]
      ppnAge5Ret3[y - 3, ] <- recRY5[y - 3, ] / recRY[y - 3, ]
      ppnAge5Ret2[y - 2, ] <- recRY5[y - 2, ] / recRY[y - 2, ]
      ppnAge5Ret1[y - 1, ] <- recRY5[y - 1, ] / recRY[y - 1, ]

      ppnAge6Ret5[y - 5, ] <- recRY6[y - 5, ] / recRY[y - 5, ]
      ppnAge6Ret4[y - 4, ] <- recRY6[y - 4, ] / recRY[y - 4, ]
      ppnAge6Ret3[y - 3, ] <- recRY6[y - 3, ] / recRY[y - 3, ]
      ppnAge6Ret2[y - 2, ] <- recRY6[y - 2, ] / recRY[y - 2, ]
      ppnAge6Ret1[y - 1, ] <- recRY6[y - 1, ] / recRY[y - 1, ]

      recRY[recRY == (0.1 * extinctThresh)] <- 0

      # Get observed recruitment by brood year
      obsBYLag<-ageMaxRec + 1


      # Loop over number of years in which fish that came from spawners in BY will recruit
      obsRecBY_noAgeErr[y-obsBYLag,]<-0
      for (i in 1:(obsBYLag-ageFirstRec)) {
        ppnAgeName<-paste("ppnAge",obsBYLag-i,"Ret",i,sep="")
        obsRecBY_noAgeErr[y-obsBYLag,] <- obsRecBY_noAgeErr[y-obsBYLag,] +
          obsRecRY[y-i,] *eval(as.name( ppnAgeName))[y-i,]
      }

      #Observed age composition of spawners (all ages in the last three years, which are needed
      #for multivariate logistic error in obs ages)
      randAges[y, ] <- runif(nAges, 0.001, 0.999)

      for (k in 1:nCU) {
        # Compile age composition of returns 5 years ago
        ppnSRet5 <- c(ppnAge2Ret5[y - 5, k], ppnAge3Ret5[y - 5, k],
                      ppnAge4Ret5[y - 5, k], ppnAge5Ret5[y - 5, k],
                      ppnAge6Ret5[y - 5, k])
        # Use ifelse statements to test if obs RecYR from 5 years ago are equal
        # to true RecRY values
        #  -- If False, generate observed age composition data with error
        ppnObsSRet5[y - 5, k, ] <- ifelse(rep(obsRecRY[y - 5,  k] ==
                                                recRY[y - 5,  k],
                                              length.out = nAges),
                                          ppnSRet5,
                                          ppnAgeErr(ppnSRet5, ageErr,
                                                    randAges[y - 5, ]))

        # Repeat above for years y - 4, y-3, etc
        ppnSRet4 <- c(ppnAge2Ret4[y - 4, k], ppnAge3Ret4[y - 4, k],
                      ppnAge4Ret4[y - 4, k], ppnAge5Ret4[y - 4, k],
                      ppnAge6Ret4[y - 4, k])
        ppnObsSRet4[y - 4, k, ] <- ifelse(rep(obsRecRY[y - 4, k] ==
                                                recRY[y - 4, k],
                                              length.out = nAges),
                                          ppnSRet4,
                                          ppnAgeErr(ppnSRet4, ageErr,
                                                    randAges[y - 4, ]))
        ppnSRet3 <- c(ppnAge2Ret3[y - 3, k], ppnAge3Ret3[y - 3, k],
                      ppnAge4Ret3[y - 3, k], ppnAge5Ret3[y - 3, k],
                      ppnAge6Ret3[y - 3, k])
        ppnObsSRet3[y - 3, k, ] <- ifelse(rep(obsRecRY[y - 3, k] ==
                                                recRY[y - 3, k],
                                              length.out = nAges),
                                          ppnSRet3,
                                          ppnAgeErr(ppnSRet3, ageErr,
                                                    randAges[y - 3, ]))
        ppnSRet2 <- c(ppnAge2Ret2[y - 2, k], ppnAge3Ret2[y - 2, k],
                      ppnAge4Ret2[y - 2, k], ppnAge5Ret2[y - 2, k],
                      ppnAge6Ret2[y - 2, k])
        ppnObsSRet2[y - 2, k, ] <- ifelse(rep(obsRecRY[y - 2, k] ==
                                                recRY[y - 2, k],
                                              length.out = nAges),
                                          ppnSRet2,
                                          ppnAgeErr(ppnSRet2, ageErr,
                                                    randAges[y - 2, ]))
        ppnSRet1 <- c(ppnAge2Ret1[y - 1, k], ppnAge3Ret1[y - 1, k],
                      ppnAge4Ret1[y - 1, k], ppnAge5Ret1[y - 1, k],
                      ppnAge6Ret1[y - 1, k])
        ppnObsSRet1[y - 1, k, ] <- ifelse(rep(obsRecRY[y - 1,k] ==
                                                recRY[y - 1, k],
                                              length.out = nAges),
                                          ppnSRet1,
                                          ppnAgeErr(ppnSRet1, ageErr,
                                                    randAges[y - 1, ]))
      }

      #Sum returns by age class ppn structured by yr to get recruits by brood yr
      #(i.e. 3 years old 3 yrs prior, 4 yr olds 2 yrs prior)
      #Only calc. obsRecBY after sim has been running for y=obsLag, otherwise
      #overwriting true observations

      if (y > (nPrime + obsBYLag)) {
        #         for (k in 1:nCU) {
        #           if (species == "coho") {
        #             obsRecBY[y - obsBYLag, k] <- obsRecRY[y - 3, k] * ppnObsSRet3[y - 3, k, 1] +
        #             obsRecRY[y - 2, k] * ppnObsSRet2[y - 2, k, 2] +
        #             obsRecRY[y - 1, k] * ppnObsSRet1[y - 2, k, 3]
        #           }
        #
        #           if (species == "sockeye") {
        #             obsRecBY[y - obsBYLag, k] <- obsRecRY[y - 4, k] * ppnObsSRet4[y - 4, k, 1] +
        #                obsRecRY[y - 3, k] * ppnObsSRet3[y - 3, k, 2] +
        #                obsRecRY[y - 2, k] * ppnObsSRet2[y - 2, k, 3] +
        #                obsRecRY[y - 1, k] * ppnObsSRet1[y - 1, k, 4]
        #           }
        #
        #           if (obsS[y - obsBYLag, k] == 0) {
        #             obsRecBY[y - obsBYLag, k] <- 0
        #           }
        #           if (is.na(obsRecBY[y - obsBYLag, k])) {
        #             obsRecBY[y - obsBYLag, k] <- 0
        #           }
        #           obsLogRS[y - obsBYLag, k] <- ifelse(obsRecBY[y - obsBYLag, k] == 0,
        #                                             0,
        #                                             log(obsRecBY[y - obsBYLag, k] /
        #                                                   obsS[y - obsBYLag, k]))
        #         } # end of nCU loop over k's



        ## Non-species-specific alternative to test:
        # # Loop over number of years in which fish that came from spawners in BY will recruit
        obsRecBY[y-obsBYLag,] <- 0
        for (i in 1:(obsBYLag-ageFirstRec)) {
          ppnAgeName<-paste("ppnObsSRet",i,sep="")
          obsRecBY[y-obsBYLag,] <- obsRecBY[y-obsBYLag,] +
            obsRecRY[y-i,] * eval(as.name( ppnAgeName))[y-i,,ageMaxRec-i]
          #note: additional "-1" added to index for colum of ppnAgeName to offset first column starting at age 2
        }

        # Force obsRecBY to zero when S is 0 or NA:
        obsRecBY[y - obsBYLag, which(obsS[y - obsBYLag, ] == 0)] <- 0
        obsRecBY[y - obsBYLag, which(is.na(obsS[y - obsBYLag, ]) == TRUE)] <- 0

        # Calculate obs log R/S:
        obsLogRS[y - obsBYLag, ] <- log(obsRecBY[y - obsBYLag, ] /
                                          obsS[y - obsBYLag, ])
        # -- and, force obs log R/S to 0 when obsS is 0
        obsLogRS[y - obsBYLag, which(obsS[y - obsBYLag, ] == 0)] <- 0

        # calculate aggregate obsRecBY over all CUs
        obsRecBYAg[y - obsBYLag, n] <- sum(obsRecBY[y - obsBYLag, ])


      } # end of if (y > (nPrime + obsBYLag))

      ### Management submodel (i.e. HCRs and catch)

      #Calculate forecasts and benchmarks at CU level
      #If no forecast mean available, assume it's obs perfectly when setting ERs
      # If any CU has a forecast mean of NA, set all obsErr forecast values to 1
      if (sum(is.na(forecastMean)) > 0) {
        obsErrDat[["forecast"]] <- rep(1,nCU)
      } else {
        # parameterize cautiously; forecastSig values should be constrained to
        # smaller values than example input data (e.g. forecastSig = 3.6 can
        # produce absurdly large errors)
        obsErrDat[["forecast"]] <- exp(qnorm(runif(nCU, 0.0001, 0.9999),
                                             log(forecastMean), forecastSig))
      }

      foreRecRY[y, ] <- obsErrDat[["forecast"]] * recRY[y, ]
      foreRecRY[y, ] <- sapply(foreRecRY[y, ],
                               function(x) ifelse(x < extinctThresh,
                                                  extinctThresh,
                                                  x))
      foreRecErr[y, ] <- obsErrDat$forecast

      # Calculate MU-level recruitment and forecasted recruitment associated with each CU
      for (k in 1:nCU) {
        #identify CUs with same MU and sum their forecasts
        MUs <- which(manUnit %in% manUnit[k])
        recRYManU[y, k] <- sum(recRY[y, MUs])
        foreRecRYManU[y, k] <- sum(foreRecRY[y, MUs])
      }

      # Create vector of ppnMix for each CU
      if (ppnMix == "flex") {
        ppnMixVec <- rep(1, length.out = nCU)
      } else {
        ppnMixVec <- rep(as.numeric(ppnMix), length.out = nCU)
      }

      #replace forecasted recruitment with true recruitment (by MU)
      # tacs <- calcTAC(rec = recRYManU[y, ], canER = ER,
      #                 harvContRule = harvContRule, amER = 0,
      #                 ppnMixVec = ppnMixVec,
      #                 manAdjustment = manAdjustment,
      #                 lowFRP = lowRefPt[y, ], highFRP = highRefPt[y, ],
      #                 minER = minER, maxER = maxER,
      #                 overlapConstraint = overlapConstraint[y, ],
      #                 constrainMix = constrainMix)


      if (harvContRule == "fixedER") {
        tacs <- calcTAC_fixedER(rec = recRYManU[y, ], canER=canER, amER = amER, ppnMixVec)

      }


      # Proportion of MU-level RY recruitement that each CU accounts for
      truePpn <- recRY[y, ] / recRYManU[y, ]
      forecastPpn <- foreRecRY[y, ] / foreRecRYManU[y, ]

      #Adjust TAC by proportions so that CU-specific harvest rates can be calc;
      #Adjust mixed fishery TAC by true ppns because these are not influenced by
      #management, single TAC by forecasted ppns
      # Note: with sockeye amTAC often lower than expected (~17%) due to AFE
      amTAC[y, ] <- tacs[['amTAC']] * truePpn
      amTAC[is.na(amTAC)] <- 0
      mixTAC[y, ] <- tacs[['canMixTAC']] * truePpn
      mixTAC[is.na(mixTAC)] <- 0
      singTAC[y, ] <- tacs[['canSingTAC']] * truePpn
      singTAC[is.na(singTAC)] <- 0

      # Apply single CU HCR
      if (singleHCR != FALSE) {

        for (k in 1:nCU) {

          if (singleHCR == "forecast") {
            singCUStatus[y, k] <- max(0,
                                      (foreRecRY[y, k]) -
                                        amTAC[y, k] - mixTAC[y, k])
          }

          if (singleHCR == "retro") {

            if (model[k] == "ricker") {
              #Secondary HCR option 2 uses median obsS abundance over previous gen
              #Needs to be in for loop because calc depends on whether stock is
              #cyclic or not; should eventually be replaced w/ estimated BMs

              singCUStatus[y, k] <- median(obsS[(y - 1):(y - gen), k])


              if (singCUStatus[y, k] >= lowerObsBM[y - 1, k]) {
                counterSingleBMLow[y, k] <- 1
              }
              if (singCUStatus[y, k] >= upperObsBM[y - 1, k]) {
                counterSingleBMHigh[y, k] <- 1
              }
            } # end if (model[k] == "ricker")

            #Larkin HCR only applies to dominant line (no median)
            if (model[k] == "larkin" & cycle[y] == domCycle[k]) {
              singCUStatus[y, k] <- obsS[y - gen, k]

              #NOTE that BMs for Larkin stocks use only dominant cycle line so
              #aligning lowerBM with cycle line is not necessary
              if (singCUStatus[y, k] >= lowerObsBM[y - 1, k]) {
                counterSingleBMLow[y, k] <- 1
              }
              if (singCUStatus[y, k] >= upperObsBM[y - 1, k]) {
                counterSingleBMHigh[y, k] <- 1
              }
            } #end if (model[k] == "larkin")

          } # end of if (singleHCR == "retro")


        } #end for k in 1:nCU
      } #end if singleHCR != FALSE


      # if there is no single stock HCR applied than all CUs assumed to be
      # "above" the lower BM so that TAC is taken
      if (singleHCR == FALSE) {
        counterSingleBMLow[y, ] <- 1
      }

      singTAC[y, ] <- singTAC[y, ] * counterSingleBMLow[y, ]
      singTAC[is.na(singTAC)] <- 0
      canTAC[y, ] <- apply(rbind(mixTAC[y, ], singTAC[y, ]), 2, sum)
      totalTAC[y, ] <- apply(rbind(amTAC[y, ], mixTAC[y, ], singTAC[y, ]), 2,
                             sum)

      # Calculate realized catches (if statement necessary because
      # calcRealCatch has random number generator reset)
      if (random != TRUE) {

        amCatch[y, ] <- calcRealCatch(recRY[y, ], amTAC[y, ], sigma = mixOUSig,
                                      setSeedInput = n * y)
        remRec1 <- pmax(recRY[y, ] - amCatch[y, ], 0)
        mixCatch[y, ] <- calcRealCatch(remRec1, mixTAC[y, ], sigma = mixOUSig,
                                       setSeedInput = n * y)
        remRec2 <- pmax(remRec1 - mixCatch[y, ] - extinctThresh, 0)
        #  migMortRate[y, ] <- enRouteMR * migMortErr
        # migMort1 <- remRec2 * (preFMigMort * migMortRate[y, ])
        #  remRec3 <- pmax(remRec2 - migMort1 - extinctThresh, 0)
        singCatch[y, ] <- calcRealCatch(remRec2, singTAC[y, ],
                                        sigma = singOUSig, setSeedInput = n * y)
      } else {
        amCatch[y, ] <- calcRealCatch(recRY[y, ], amTAC[y, ], sigma = mixOUSig,
                                      random = TRUE)
        remRec1 <- pmax(recRY[y, ] - amCatch[y, ], 0)
        mixCatch[y, ] <- calcRealCatch(remRec1, mixTAC[y, ], sigma = mixOUSig,
                                       random = TRUE)
        remRec2 <- pmax(remRec1 - mixCatch[y, ] - extinctThresh, 0)
        #  migMortRate[y, ] <- enRouteMR * migMortErr
        #  migMort1 <- remRec2 * (preFMigMort * migMortRate[y, ])
        #  remRec3 <- pmax(remRec2 - migMort1 - extinctThresh, 0)
        singCatch[y, ] <- calcRealCatch(remRec2, singTAC[y, ],
                                        sigma = singOUSig, random = TRUE)
      }

      remRec3 <- pmax(remRec2 - singCatch[y, ] - extinctThresh, 0)

      #summary catch calculations
      totalCatch[y, ] <- amCatch[y, ] + mixCatch[y, ] + singCatch[y, ]
      amCatchAg[y, n] <- sum(amCatch[y, ])
      mixCatchAg[y, n] <- sum(mixCatch[y, ])
      catchAg[y, n] <- sum(amCatch[y, ] + mixCatch[y, ] + singCatch[y, ])
      amExpRate[y, ] <- amCatch[y, ] / recRY[y, ]
      mixExpRate[y, ] <- mixCatch[y, ] / recRY[y, ]
      singExpRate[y, ] <- singCatch[y, ] / recRY[y, ]
      singExpRate[y, ][is.na(singExpRate[y, ])] <- 0
      expRate[y, ] <- (amCatch[y, ] + mixCatch[y, ] + singCatch[y, ]) /
        recRY[y, ]

      expRate[recRY == 0] <- 0
      expRateAg[y, n] <- ifelse(recRYAg[y, n] == 0, 0,
                                catchAg[y, n] / recRYAg[y, n])

      # Calculate target harvest rate for Canadian fisheries only
      targetCanER[y, ] <- apply(rbind(mixTAC[y, ], singTAC[y, ]), 2, sum) /
        recRY[y, ]
      targetExpRateAg[y, n] <- ifelse(harvContRule == "fixedER",
                                      canER + amER,
                                      sum(totalTAC[y, ]) / recRYAg[y, n])


      S[y, ] <- recRY[y, ] * (1 - expRate[y, ])

      for (k in 1:nCU) {
        if (S[y, k] < extinctThresh) {
          S[y, k] <- 0
        }
      }
      sAg[y, n] <- sum(S[y, ])


      #___________________________________________________________________
      ### Observation submodel 2 (this years abundance)
      # Generate MU-specific observation error
      for (m in seq_along(muName)) {
        #catch observation error for both Can and US fisheries
        obsErrDat[obsErrDat$mu == muName[m], "mixC"] <- exp(qnorm(
          runif(1, 0.0001, 0.9999), 0, obsMixCatchSig))
      }
      #observed spawner error; also used to generate en route mortality estimate
      obsErrDat[["spwn"]] <- exp(qnorm(runif(nCU, 0.0001, 0.9999), 0, obsSig))
      #CU specific draws from shared distribution
      obsErrDat[["singC"]] <- exp(qnorm(runif(nCU, 0.0001, 0.9999), 0,
                                        obsSingCatchSig))

      obsS[y, ] <- S[y, ] * obsErrDat[["spwn"]]
      obsSAg[y, n] <- sum(obsS[y, ])
      obsAmCatch[y, ] <- calcObsCatch(amCatch[y, ], recRY[y, ], manUnit,
                                      tauCatch, stkID, obsErrDat[["mixC"]],
                                      extinctThresh)
      obsMixCatch[y, ] <- calcObsCatch(mixCatch[y, ], recRY[y, ], manUnit,
                                       tauCatch, stkID, obsErrDat[["mixC"]],
                                       extinctThresh)
      #all fish caught in terminal fishery assumed to originate there
      obsSingCatch[y, ] <- singCatch[y, ] * obsErrDat[["singC"]]
      obsSingCatch[y, ] <- ifelse(obsSingCatch[y, ] < (0.1 * extinctThresh), 0,
                                  obsSingCatch[y, ])
      obsCatchAg[y, n] <- sum(obsAmCatch[y, ] + obsMixCatch[y, ] +
                                obsSingCatch[y, ],
                              na.rm = TRUE)
      obsRecRY[y, ] <- obsS[y, ] +  obsAmCatch[y, ] +
        obsMixCatch[y, ] + obsSingCatch[y, ]
      obsRecRYAg[y, n] <- sum(obsRecRY[y, ], na.rm = TRUE)
      obsExpRate[y, ] <- ifelse((obsAmCatch[y, ] + obsMixCatch[y, ] +
                                   obsSingCatch[y, ]) / obsRecRY[y, ] > 1,
                                1,
                                ifelse(obsRecRY[y, ] == 0,
                                       0,
                                       (obsAmCatch[y, ] + obsMixCatch[y, ] +
                                          obsSingCatch[y, ]) / obsRecRY[y, ]))
      obsExpRate[obsRecRY == 0] <- 0
      obsExpRateAg[y, n] <- mean(obsExpRate[y, ])


      #___________________________________________________________________
      ### Assessment submodel

      #Build SRR (even with normative period useful for diagnostics)
      for (k in 1:nCU) {
        srMod <- quickLm(xVec = obsS[, k], yVec = obsLogRS[, k])
        estYi[y, k, n] <- srMod[[1]]
        estSlope[y, k, n] <- srMod[[2]]
        estRicB[y, k, n] <- ifelse(extinct[y, k] == 1, NA, -estSlope[y, k, n])
        estRicA[y, k, n] <- ifelse(extinct[y, k] == 1, NA, estYi[y, k, n])
      }
      #Benchmark estimates
      # -- If normative period is TRUE than do not estimate BMs (they will diverge
      # from true BMs for reasons other than obs error)
      if (normPeriod == TRUE) {
        s25th[y, , n] <- s25th[nPrime, , n]
        s50th[y, , n] <- s50th[nPrime, , n]
        estS25th[y, , n] <- s25th[nPrime, , n]
        estS50th[y, , n] <- s50th[nPrime, , n]
        estSGen[y, , n] <- sGen[y, , n]
        estSMSY[y, , n] <- sMSY[y, , n]
      } else if (normPeriod == FALSE) {

        if (model[k] == "larkin" | model[k]=="rickerSurv") {
          warning("Normative period is TRUE for Larkin or RickerSurv. Assessment model will not equal true OM for this case. Default assessment model is Ricker")
        }

        for (k in 1:nCU) {
          #Calculate true percentile BMs
          temp <- S[1:y, k]
          sNoNA <- temp[!is.na(temp)]
          n25th <- round(length(sNoNA) * 0.25, 0)
          n50th <- round(length(sNoNA) * 0.50, 0)
          s25th[y, k, n] <- sort(sNoNA)[n25th]
          s50th[y, k, n] <- sort(sNoNA)[n50th]
          #Calculate observed percentile BMs
          temp <- obsS[1:y, k]
          obsSNoNA <- temp[!is.na(temp)]
          obsN25th <- round(length(obsSNoNA) * 0.25, 0)
          obsN50th <- round(length(obsSNoNA) * 0.50, 0)
          estS25th[y, k, n] <- sort(obsSNoNA)[obsN25th]
          estS50th[y, k, n] <- sort(obsSNoNA)[obsN50th]
          #Calculate SR BMs
          estSMSY[y, k, n] <- ifelse(extinct[y, k] == 1, NA,
                                     (1 - gsl::lambert_W0(exp(
                                       1 - estRicA[y, k, n]))) /
                                       estRicB[y, k, n])
          if (is.na(estRicB[y, k, n]) == FALSE) {
            if (estRicB[y, k, n] > 0) {
              if ((1 / estRicB[y, k, n]) <= max(obsS[,k], na.rm = TRUE) * 4) {
                estSGen[y, k, n] <- as.numeric(sGenSolver(
                  theta = c(estRicA[y, k, n], estRicB[y, k, n], ricSig[k]),
                  sMSY = estSMSY[y, k, n]))
              } else {
                #if a BM cannot be estimated set it to the last estimated value
                estSGen[y, k, n] <- estSGen[max(which(!is.na(
                  estSGen[, k, n]))), k, n]
                estSMSY[y, k, n] <- estSMSY[max(which(!is.na(
                  estSMSY[, k, n]))), k, n]
                fb1[y, k] <- 1
              }
            }# End of if(estRicB[y, n]>0)
          } else{
            estSGen[y, k, n] <- estSGen[max(which(!is.na(estSGen[, k, n]))), k, n]
            estSMSY[y, k, n] <- estSMSY[max(which(!is.na(estSMSY[, k, n]))), k, n]
            fb2[y, k] <- 1
          } #end(if(is.na(estRicB)))
        } #end for(k in 1:nCU)
        if (is.na(estSGen[y, k, n]) == FALSE & is.na(estSMSY[y, k, n]) == FALSE) {
          if (estSGen[y, k, n] > estSMSY[y, k, n]) {
            warning("Estimated lower benchmark greater than upper benchmark;
                    set to NA")
            estSGen[y, k, n] <- estSGen[max(which(!is.na(estSGen[, k, n]))), k, n]
            estSMSY[y, k, n] <- estSMSY[max(which(!is.na(estSMSY[, k, n]))), k, n]
            fb3[y, k] <- 1
          }
        }
      } #end if normPeriod = FALSE


      #___________________________________________________________________
      ### Population dynamics submodel


      # Get marine rvival covariate (only used for rickerSurv SR model)
      ## - all CUs have the same marine survival (only option available at present)
      if (y == nPrime+1) mSurvAge4[1:nPrime, n] <- rep(exp(mu_logCoVar),nPrime)
      mSurvAge4[y, n]<-rnorm(1,mu_logCoVar,sig_logCoVar)
      if (mSurvAge4[y, n] > max_logCoVar) mSurvAge4[y, n] <- max_logCoVar
      if (mSurvAge4[y, n] < min_logCoVar) mSurvAge4[y, n] <- min_logCoVar

      for (k in 1:nrow(ageStruc)) {
        ppnAges[y, k, ] <- ppnAgeErr(ageStruc[k, ], tauAge[k],
                                     error = runif(nAges, 0.0001, 0.9999))
      }
      if (prod == "skew") {
        #draw process variance from skewed normal distribution
        errorCU[y, ] <- sn::rmst(n = 1, xi = rep(0, nCU),
                                 alpha = rep(log(0.65), nCU), nu = 1000,
                                 Omega = covMat)
      } else if (prod == "skewT") {
        #draw process variance from skewed T distribution
        errorCU[y, ] <- sn::rmst(n = 1, xi = rep(0, nCU),
                                 alpha = rep(log(0.65), nCU), nu = 2,
                                 Omega = covMat)
      } else if (prod == "studT" | prod == "lowStudT") {
        #draw process variance from student T distribution
        errorCU[y, ] <- sn::rmst(n = 1, xi = rep(0, nCU),
                                 alpha = rep(0, nCU), nu = 2,
                                 Omega = covMat)
      } else { #otherwise draw from normal
        errorCU[y, ] <- sn::rmst(n = 1, xi = rep(0, nCU),
                                 alpha = rep(0, nCU), nu = 10000,
                                 Omega = covMat)
      }


      for (k in 1:nCU) {
        if (S[y, k] > 0) {
          if (model[k] == "ricker") {
            dum <- rickerModel(S[y, k], alphaMat[y, k], beta[k],
                               error = errorCU[y, k], rho = rho,
                               prevErr = laggedError[y - 1, k])
            laggedError[y, k] <- dum[[2]]
            #keep recruitment below CU-specific cap
            recBY[y, k] <- min(dum[[1]], recCap[k])
          }
          if (model[k] == "rickerSurv") {
            # K.Holt - current set-up assumes ageMaxRec == 4 (i.e., IF coho); will need to update for other stocks


            mSurvAtAge<-c(mSurvAge4[y-2,n],mSurvAge4[y-1,n],mSurvAge4[y,n],0,0)
            dum <- rickerSurvModel(S=S[y, k], a=alphaMat[y, k], b=beta[k],
                                   ppnAges = ppnAges[y,k,], gamma=gamma[k],
                                   mSurvAtAge=mSurvAtAge,
                                   error = errorCU[y, k])
            #keep recruitment below CU-specific cap
            recBY[y, k] <- min(dum[[6]], recCap[k])
          }
          if (model[k] == "larkin") {
            dum <- larkinModel(S[y, k], S[y - 1, k], S[y-2, k], S[y - 3, k],
                               alphaMat[y, k], beta[k], larB1[k], larB2[k],
                               larB3[k], error = errorCU[y, k])
            #keep recruitment below CU-specific cap
            recBY[y, k] <- min(dum[[1]], recCap[k])
          }

          logRS[y, k] <- log(recBY[y, k] / S[y, k])
        } #end if(S[y, k]>0)
        if (is.na(laggedError[y, k])) {
          laggedError[y, k] <- 0
        }
        if (S[y, k] == 0) {
          recBY[y, k] <- 0
          logRS[y, k] <- 0
        }
        if (recBY[y, k] <= extinctThresh) {
          recBY[y, k] <- 0
        }
        extinct[y, ] <- extinctionCheck(y = y, gen = gen,
                                        extinctThresh = extinctThresh,
                                        spwnMat = S)
      } #end for(k in 1:nCU)

      recBYAg[y, n] <- sum(recBY[y, ])

      for (k in 1:nCU) {
        if (model[k] == "ricker" | model[k] == "rickerSurv"  |
            model[k] == "larkin" & cycle[y] == domCycle[k]) {
          if (bm == "stockRecruit") {
            #is spawner abundance greater than upper/lower BM
            upperBM[y, k] <- 0.8 * sMSY[y, k, n]
            lowerBM[y, k] <- sGen[y, k, n]
            upperObsBM[y, k] <- 0.8 * estSMSY[y, k, n]
            lowerObsBM[y, k] <- estSGen[y, k, n]
          }
          if (bm == "percentile") {
            upperBM[y, k] <- s50th[y, k, n]
            lowerBM[y, k] <- s25th[y, k, n]
            upperObsBM[y, k] <- estS50th[y, k, n]
            lowerObsBM[y, k] <- estS25th[y, k, n]
          }
        }
        #only save status for Larkin stocks on dom cycle otherwise use prev status
        if (model[k] == "larkin" & cycle[y] != domCycle[k]) {
          upperBM[y, k] <- upperBM[y - 1, k]
          lowerBM[y, k] <- lowerBM[y - 1, k]
          upperObsBM[y, k] <- upperObsBM[y - 1, k]
          lowerObsBM[y, k] <- lowerObsBM[y - 1, k]
        }
        if (!is.na(upperBM[y, k]) & S[y, k] > upperBM[y, k]) {
          counterUpperBM[y, k] <- 1 #is spawner abundance greater than upper BM
        }

        if (!is.na(lowerBM[y, k]) & S[y, k] > lowerBM[y, k]) {
          counterLowerBM[y, k] <- 1 #is spawner abundance greater than lower BM
        }
        if (!is.na(upperObsBM[y, k]) & obsS[y, k] > upperObsBM[y, k]) {
          counterUpperObsBM[y, k] <- 1 #is spawner abundance greater than upper BM
        }
        if (!is.na(lowerObsBM[y, k]) & obsS[y, k] > lowerObsBM[y, k]) {
          counterLowerObsBM[y, k] <- 1 #is spawner abundance greater than lower BM
        }

      } # end of k loop over CUs

      } # end of loop 3

    #__________________________________________________________________________
    ### Draw one trial for plotting
    if (n == drawTrial) {
      # Generate indices of observation error
      obsSErr <- calcErr(obsS, S)
      obsRecBYErr <- calcErr(obsRecBY, recBY)
      obsRecRYErr <- calcErr(obsRecRY, recRY)
      obsMixCatchErr <- calcErr(obsMixCatch, mixCatch)
      obsSingCatchErr <- calcErr(obsSingCatch, singCatch)
      obsAmCatchErr <- calcErr(obsAmCatch, amCatch)
      obsExpErr <- calcErr(obsExpRate, expRate)
      forecastErr <- calcErr(foreRecRY, recRY)

      # Combine relevant data into array passed to plotting function
      varNames <- c("Productivity", "Est Productivity", "Est Beta",
                    "Spawners", "Obs Spawners", "Recruits BY",
                    "Obs Recruits BY", "Recruits RY",
                    "Mix Catch", "Single Catch", "US Catch",
                    "Total Exp Rate", "Mix Exp Rate", "Single Exp Rate",
                    "Smsy", "Sgen", "S 50th Percentile", "S 25th Percentile",
                    "Obs Spawners Err", "Obs RecBY Err",
                    "Obs RecRY Err", "Obs Mix Catch Err",
                    "Obs Single Catch Err",
                    "Obs US Catch Err", "Obs Exp Rate Err", "Forecast Err"
      )

      plotTrialDat <- array(c(alphaMat, estRicA[ , , n], estRicB[ , , n],
                              S, obsS, recBY, obsRecBY, recRY,
                              mixCatch, singCatch, amCatch,
                              expRate, mixExpRate, singExpRate,
                              sMSY[ , , n], sGen[ , , n], s50th[ , , n],
                              s25th[ , , n], obsSErr, obsRecBYErr,
                              obsRecRYErr, obsMixCatchErr, obsSingCatchErr,
                              obsAmCatchErr, obsExpErr, forecastErr),
                            dim = c(nYears, nCU, length(varNames)))
      dimnames(plotTrialDat)[[3]] <- varNames

      # Figure settings
      fileName <- ifelse(variableCU == "TRUE",
                         paste(cuNameOM, cuNameMP, "singleTrialFig.pdf",
                               sep = "_"),
                         paste(nameOM, nameMP, "singleTrialFig.pdf", sep = "_"))
      pdf(file = paste(here("outputs/diagnostics", dirPath, fileName),
                       sep = "/"), height = 6, width = 7)
      if (exists("larB")) { # if larkin terms are present they need to be passed
        larBList <- list(larB, larB1, larB2, larB3)
        names(larBList) <- c("lag0", "lag1", "lag2", "lag3")
        plotDiagCU(plotTrialDat, varNames, stkName, model, ricB,
                           larBList = larBList, medAbundance, nPrime, extinct,
                           focalCU = NULL)
      } else {
        plotDiagCU(plotTrialDat, varNames, stkName, model, ricB,
                           larBList = NULL, medAbundance, nPrime, extinct,
                           focalCU = NULL)
      }
      dev.off()
    }



    #_____________________________________________________________________
    ## Summary calculations with full datasets
    # Aggregate BM status
    for (k in 1:nCU) {
      #did a CU recover early, i.e. above BM year of 4th generation
      if(sum(counterUpperBM[earlyPeriod, k]) == gen) {
        counterEarlyUpperBM[n, k] <- 1
      }
      if(sum(counterLowerBM[earlyPeriod, k]) == gen) {
        counterEarlyLowerBM[n, k] <- 1
      }
      #has a CU "recovered", i.e. above BM every year of last generation
      if(sum(counterUpperBM[(nYears - gen + 1):nYears, k]) == gen) {
        counterLateUpperBM[n, k] <- 1
      }
      if(sum(counterLowerBM[(nYears - gen + 1):nYears, k]) == gen) {
        counterLateLowerBM[n, k] <- 1
      }
      if(sum(counterUpperObsBM[(nYears - gen + 1):nYears, k]) == gen) {
        counterLateUpperObsBM[n, k] <- 1
      }
      if(sum(counterLowerObsBM[(nYears - gen + 1):nYears, k]) == gen) {
        counterLateLowerObsBM[n, k] <- 1
      }
    }


    # Calculate CU-specific trends across geometric means
    ## NOTE: even with quickLm these functions are a big bottleneck
    ## Make an optional output PM (rather than default) based on simPar
    if (!is.null(simPar$statusTrendPM)) {
      if (simPar$statusTrendPM == TRUE) {
        for (k in 1:nCU) {
          S[, k] <- ifelse(S[, k] == 0, 0.00005, S[, k])
          #ask C Holt how this is generated, right now running blind
          nYrsGeoMean <- nYears - nPrime
          #generate geo mean spwnr abundance w.out NAs
          strtGMean <- nPrime - gen
          sGeoMean[strtGMean:nYears, k] <- genMean(S[strtGMean:nYears, k], gen)
          lnSGeoMean <- log(sGeoMean)
          ppnSLow[n, k] <- length(which(sGeoMean[, k] < 0.1)) / nYrsGeoMean
          sl <- NA
          ppnChange <- rep(NA, nYears)
          # calculate slope over 12 year period i.e.,from year j-11 to j
          for (j in (nPrime + (3 * gen - 1)):nYears) {
            if (extinct[j, k] == 0) { #only calculate slopes when non-extinct
              if (is.na(lnSGeoMean[j - (3 * gen - 1), k]) == "FALSE") {
                sl[j] <- quickLm(xVec = c(1:(3 * gen)),
                                 yVec = lnSGeoMean[(j - (3 * gen - 1)):j, k])[[2]]
                ppnChange[j] <- exp(sl[j] * 3 * gen) - 1
              }
            }
            if (extinct[j, k] == 1) {
              ppnChange[j] <- NA
            }
          }

          ppnChangeMat[, k, n] <- ppnChange
          ppnChangeNoNA <- which(is.na(ppnChange)[1:nYears] == "FALSE")
          ppnYrsCOS[n, k] <- ifelse(length(ppnChangeNoNA) == 0,
                                    0,
                                    length(which(ppnChange[1:nYears] > -0.3)) /
                                      length(ppnChangeNoNA))
          ppnYrsWSP[n, k] <-ifelse(length(ppnChangeNoNA) == 0,
                                   0,
                                   length(which(ppnChange[1:nYears] > -0.25)) /
                                     length(ppnChangeNoNA))
          #invert to make positive (i.e. same directionality as above BMs)
          S[, k][which(S[, k] == 0.00005)] <- 0
        }
      }
    }

    #__________________________________________________________________________
    ## Store trial specific outputs
    # Store diagnostic outputs
    spwnrArray[ , , n] <- S # these arrays generated to pass to synch list
    recArray[ , , n] <- recBY
    alphaArray[ , , n] <- alphaMat
    returnArray[ , , n] <- recRY
    logRSArray[ , , n] <- logRS
    obsTotalCatch <- obsAmCatch + obsMixCatch + obsSingCatch
    recDevArray[ , , n] <- errorCU
    singCatchArray[ , , n] <- singCatch
    singTACArray[ , , n] <- singTAC
    totalCatchArray[ , , n] <- totalCatch


    #Store trial and CU specific means, variances, and proportions to add to
    #aggregate data frame
    #NOTE CHANGE IF FIXED ERs VARY AMONG CUs
    targetER[n, ] <- ifelse(harvContRule == "fixedER",
                            rep(canER, nCU),
                            apply(targetCanER[(nPrime+1):nYears, ], 2,
                                  function(x) mean(x, na.rm = TRUE)))
    #data of interest
    yrsSeq <- seq(from = nPrime + 1, to = nYears, by = 1)
    #median and CVs of true or obs PMs through time per trial
    medS[n, ] <- apply(na.omit(S[yrsSeq, ]), 2, median)
    varS[n, ] <- apply(na.omit(S[yrsSeq, ]), 2, cv)
    medObsS[n, ] <- apply(na.omit(obsS[yrsSeq, ]), 2, median)
    varObsS[n, ] <- apply(na.omit(obsS[yrsSeq, ]), 2, cv)
    medRecRY[n, ] <- apply(na.omit(recRY[yrsSeq, ]), 2, median)
    varRecRY[n, ] <- apply(na.omit(recRY[yrsSeq, ]), 2, cv)
    medRecBY[n, ] <- apply(na.omit(recBY[yrsSeq, ]), 2, median)
    varRecBY[n, ] <- apply(na.omit(recBY[yrsSeq, ]), 2, cv)
    medObsRecRY[n, ] <- apply(na.omit(obsRecRY[yrsSeq, ]), 2,
                              median)
    varObsRecRY[n, ] <- apply(na.omit(obsRecRY[yrsSeq, ]), 2, cv)
    medObsRecBY[n, ] <- apply(na.omit(obsRecBY[yrsSeq, ]), 2,
                              median)
    varObsRecBY[n, ] <- apply(na.omit(obsRecBY[yrsSeq, ]), 2, cv)
    medAlpha[n, ] <- apply(na.omit(alphaMat[yrsSeq, ]), 2, median)
    varAlpha[n, ] <- apply(na.omit(alphaMat[yrsSeq, ]), 2, cv)
    medEstAlpha[n, ] <- apply(na.omit(estRicA[yrsSeq, , n]), 2, median)
    varEstAlpha[n, ] <- apply(na.omit(estRicA[yrsSeq, , n]), 2, cv)
    medBeta[n, ] <- if (is.null(nrow(beta))) {
      beta
    } else {
      apply(na.omit(beta), 2, median)
    }#avg capacity through time per trial
    medTotalCatch[n, ] <- apply(totalCatch[yrsSeq, ], 2, median, na.rm = TRUE)
    varTotalCatch[n, ] <- apply(totalCatch[yrsSeq, ], 2, cv)
    #stability in catch
    stblTotalCatch[n, ] <- apply(totalCatch[yrsSeq, ], 2, function(x) 1 / cv(x))
    medObsTotalCatch[n, ] <- apply(obsTotalCatch[yrsSeq, ], 2, median,
                                   na.rm=TRUE)
    varObsTotalCatch[n, ] <- apply(obsTotalCatch[yrsSeq, ], 2, cv)
    stblObsTotalCatch[n, ] <- apply(obsTotalCatch[yrsSeq, ], 2,
                                    function(x) 1 / cv(x))
    #median TAC lost due to overlap constraints
    #medForgoneCatch[n, ] <- apply(forgoneCatch[yrsSeq, ], 2, median)
    medTotalER[n, ] <- apply(expRate[yrsSeq, ], 2, median, na.rm=TRUE)
    medTotalObsER[n, ] <- apply(obsExpRate[yrsSeq, ], 2, median, na.rm=TRUE)
    #median single stock ER relative to escapement across sim period
    #medTAMSingER[n, ] <- apply(tamSingER[yrsSeq, ], 2, median, na.rm=TRUE)
    #ppn of years true abundance above true upper BM per trial
    ppnYrsUpperBM[n, ] <- apply(na.omit(counterUpperBM[yrsSeq, ]), 2, mean)
    ppnYrsLowerBM[n, ] <- apply(na.omit(counterLowerBM[yrsSeq, ]), 2, mean)
    ppnYrsUpperObsBM[n, ] <- apply(na.omit(counterUpperObsBM[yrsSeq, ]), 2,
                                   mean)
    ppnYrsLowerObsBM[n, ] <- apply(na.omit(counterLowerObsBM[yrsSeq, ]), 2,
                                   mean)
    ppnYrsOpenSingle[n, ] <- apply(na.omit(counterSingleBMLow[yrsSeq, ]), 2,
                                   mean)
    #PMs for early period
    medEarlyS[n, ] <- apply(na.omit(S[(nPrime + 1):endEarly, ]), 2, median)
    medEarlyRecRY[n, ] <- apply(na.omit(recRY[(nPrime + 1):endEarly, ]), 2,
                                median)
    medEarlyTotalCatch[n, ] <- apply(na.omit(
      totalCatch[(nPrime + 1):endEarly, ]), 2, median)
    #aggregate data
    meanSingExpRate[yrsSeq, n] <- apply(singExpRate[yrsSeq, ], 1, mean)
    ppnCUsUpperBM[yrsSeq, n] <- apply(counterUpperBM[yrsSeq, ], 1, mean)
    ppnCUsLowerBM[yrsSeq, n] <- apply(counterLowerBM[yrsSeq, ], 1, mean)
    ppnCUsUpperObsBM[yrsSeq, n] <- apply(counterUpperObsBM[yrsSeq, ], 1, mean)
    ppnCUsLowerObsBM[yrsSeq, n] <- apply(counterLowerObsBM[yrsSeq, ], 1, mean)
    ppnCUsExtinct[yrsSeq, n] <- apply(extinct[yrsSeq, ], 1, mean)
    #ppnConstrained[yrsSeq, n] <- apply(overlapConstraint[yrsSeq, ], 1, mean)
    ppnCUsOpenSingle[yrsSeq, n] <- apply(counterSingleBMLow[yrsSeq, ], 1,
                                         mean)

    } # end of nTrials loop

  #____________________________________________________________________________
  ## CU-specific outputs
  # Data
  meanSMSY <- arrayMean(sMSY) # CU's average BMs through time and trials
  meanSGen <- arrayMean(sGen)
  cuList <- list(nameOM, keyVar, plotOrder, nameMP, harvContRule, stkName,
                 stkID, manUnit, targetER, cuProdTrends, meanSMSY, meanSGen,
                 medS, varS,
                 medObsS, varObsS, medRecRY, varRecRY, medRecBY, varRecBY,
                 medObsRecRY, varObsRecRY, medAlpha, varAlpha, medEstAlpha,
                 varEstAlpha, medBeta,
                 medTotalCatch, varTotalCatch, (1 / varTotalCatch),
                 medObsTotalCatch, varObsTotalCatch, (1 / varObsTotalCatch),
                 medTotalER, medTotalObsER,
                 counterEarlyUpperBM, counterEarlyLowerBM, ppnYrsUpperBM,
                 ppnYrsLowerBM, ppnYrsUpperObsBM, ppnYrsLowerObsBM, ppnYrsCOS,
                 ppnYrsWSP, medEarlyS, medEarlyRecRY, medEarlyTotalCatch,
                 ppnYrsOpenSingle)
  names(cuList) <- c("opMod", "keyVar", "plotOrder", "manProc", "hcr",
                     "stkName", "stkNumber", "manUnit", "targetER",
                     "cuProdTrends", "meanSMSY",
                     "meanSGen", "medSpawners", "varSpawners", "medObsSpawners",
                     "varObsSpawners", "medRecRY", "varRecRY", "medRecBY",
                     "varRecBY", "medObsRecRY", "varObsRecRY", "medAlpha",
                     "varAlpha", "medEstAlpha", "varEstAlpha", "medBeta",
                     "medCatch",
                     "varCatch", "stblCatch", "medObsCatch", "varObsCatch",
                     "stblObsCatch", "medTotalER", "medObsTotalER",
                     "counterEarlyUpper",
                     "counterEarlyLower", "ppnYrsUpper", "ppnYrsLower",
                     "ppnYrsEstUpper", "ppnYrsEstLower", "ppnYrsCOS",
                     "ppnYrsWSP", "medEarlyS", "medEarlyRecRY",
                     "medEarlyTotalCatch", "ppnYrsSingleOpen")
  fileName <- ifelse(variableCU == "TRUE",
                     paste(cuNameOM, cuNameMP, "cuDat.RData", sep = "_"),
                     paste(nameOM, nameMP, "cuDat.RData", sep = "_"))

  saveRDS(cuList, file = paste(here("outputs/simData"), dirPath, fileName,
                               sep = "/"), version=3)

  #_____________________________________________________________________
  ## Aggregate outputs
  # Generate array of median, upper and lower quantiles that are passed to
  # plotting function
  agNames <- c("Ag Spawners", "Obs Ag Spawners", "Ag Recruits RY",
               "Obs Ag Recruits RY", "Ag Catch", "Obs Ag Catch", "Exp Rate",
               "Obs Exp Rate", "Change Ag Catch",
               "Prop Above Upper BM", "Prop Above Lower BM",
               "Obs Prop Above Upper BM", "Obs Prop Above Lower BM")
  agDat <- array(c(sAg, obsSAg, recRYAg, obsRecRYAg, catchAg, obsCatchAg,
                   expRateAg, obsExpRateAg, ppnCUsUpperBM,
                   ppnCUsLowerBM, ppnCUsUpperObsBM, ppnCUsLowerObsBM),
                 dim = c(nYears, nTrials, length(agNames)))
  dimnames(agDat)[[3]] <- agNames


  # Save aggregate data as list to create TS plot
  agTSList <- c(list(nameOM, keyVar, plotOrder, nameMP, harvContRule,
                     targetExpRateAg, firstYr, nPrime, nYears),
                plyr::alply(agDat, 3, .dims = TRUE))
  names(agTSList)[1:9] <- c("opMod", "keyVar", "plotOrder", "manProc", "hcr",
                            "targetExpRate", "firstYr", "nPrime", "nYears")
  fileName <- ifelse(variableCU == "TRUE",
                     paste(cuNameOM, cuNameMP, "aggTimeSeries.RData",
                           sep = "_"),
                     paste(nameOM, nameMP, "aggTimeSeries.RData", sep = "_"))
  saveRDS(agTSList, file = paste(here("outputs/simData"), dirPath, fileName,
                                 sep = "/"), version=3)

  # Store aggregate data as data frame; each variable is a vector of single, trial-specific values
  yrsSeq <- (nPrime + 1):nYears
  aggDat <- data.frame(opMod = rep(nameOM, length.out = nTrials),
                       manProc = rep(nameMP, length.out = nTrials),
                       keyVar = rep(keyVar, length.out = nTrials),
                       plotOrder = rep(plotOrder, length.out = nTrials),
                       hcr = rep(harvContRule, length.out = nTrials),
                       trial = seq(from = 1, to = nTrials, length.out = nTrials),
                       targetER = apply(targetExpRateAg[(nPrime + 1):nYears, ], 2, median), #median target ER (stable through time unless TAM rule active)
                       medSpawners = apply(sAg[yrsSeq, ], 2, median), #median aggregate spawner abundance across management period (i.e. loop 3 when status is assessed)
                       varSpawners = apply(sAg[yrsSeq, ], 2, cv), #cv aggregate spawner abundance across management period
                       medObsSpawners = apply(obsSAg[yrsSeq, ], 2, median), #median aggregate estimated spawner abundance across management period (i.e. loop 3 when status is assessed)
                       varObsSpawners = apply(obsSAg[yrsSeq, ], 2, cv), #cv aggregate estimated spawner abundance across management period
                       medRecRY = apply(recRYAg[yrsSeq, ], 2, median), #median aggregate spawner abundance across management period (i.e. loop 3 when status is assessed)
                       varRecRY = apply(recRYAg[yrsSeq, ], 2, cv), #cv aggregate recruit abundance across management period
                       medRecBY = apply(recBYAg[yrsSeq, ], 2, median), #median aggregate recruit abundance across management period (i.e. loop 3 when status is assessed)
                       varRecBY = apply(recBYAg[yrsSeq, ], 2, cv), #cv aggregate recruit abundance across management period
                       medObsRecBY = apply(obsRecBYAg[yrsSeq, ], 2,
                                           function(x) median(x, na.rm = TRUE)), #median aggregate estimated recruit abundance across management period (i.e. loop 3 when status is assessed)
                       varObsRecBY = apply(obsRecBYAg[yrsSeq, ], 2, cv), #cv aggregate estimated recruit abundance across management period
                       medObsRecRY = apply(obsRecRYAg[yrsSeq, ], 2, median), #median aggregate estimated recruit abundance across management period (i.e. loop 3 when status is assessed)
                       varObsRecRY = apply(obsRecRYAg[yrsSeq, ], 2, cv), #cv aggregate estimated recruit abundance across management period
                       medSpawnersLate = apply(sAg[(nYears - (3 * gen)):nYears, ], 2, median), #median aggregate spawner abundance in last 2 generations of management period
                       medCatch = apply(catchAg[yrsSeq, ], 2, median), #median aggregate catch across management period
                       varCatch = apply(catchAg[yrsSeq, ], 2, cv), #cv aggregate catch across management period
                       stabilityCatch = apply(catchAg[yrsSeq, ], 2,
                                              function(x) 1 / cv(x)), #stability of aggregate catch across management period
                       medObsCatch = apply(obsCatchAg[yrsSeq, ], 2, median), #median aggregate catch across management period
                       varObsCatch = apply(obsCatchAg[yrsSeq, ], 2, cv), #cv aggregate catch across management period
                       stabilityObsCatch = apply(obsCatchAg[yrsSeq, ], 2,
                                                 function(x) 1 / cv(x)), #stability of obs aggregate catch across management period
                       medCatchLate = apply(catchAg[(nYears - 3*gen):nYears, ], 2, median), #median aggregate catch in last 2 generations of management period
                       medER = apply(expRateAg[yrsSeq,], 2, median), #median true aggregate ER
                       medObsER = apply(obsExpRateAg[yrsSeq, ], 2, median), #median true aggregate ER
                       #ppnYrsLowCatch = apply(lowCatchAgBM[yrsSeq, ], 2, mean), #proportion of years in management period aggregate catch is above summed catch thresholds
                       #ppnYrsHighCatch = apply(highCatchAgBM[yrsSeq, ], 2, mean), #proportion of years in management period aggregate catch is above summed catch thresholds
                       # ppnYrsEscGoal = apply(agEscGoal[yrsSeq, ], 2, mean), #ppn of years aggregate escapement goal met
                       # ppnYrsCUsLower = apply(ppnLowerBM[yrsSeq, ], 2, mean), #proportion of years at least 50% of CUs are above lower BM
                       #  ppnYrsCUsUpper = apply(ppnUpperBM[yrsSeq, ], 2, mean), #proportion of years at least 50% of CUs are above upper BM
                       #  ppnCUUpper = apply(ppnCUsUpperBM[yrsSeq, ], 2, mean), #mean proportion of CUs above upper benchmark in last generations of management period
                       #  ppnCULower = apply(ppnCUsLowerBM[yrsSeq, ], 2, mean), #mean proportion of CUs above lower benchmark in last generations of management period
                       #   ppnCUEstUpper = apply(na.omit(ppnCUsUpperObsBM), 2, mean), #proportion of CUs estimated above upper benchmark in last 2 generations of management period
                       #    ppnCUEstLower = apply(na.omit(ppnCUsLowerObsBM), 2, mean), #proportion of CUs estimated above lower benchmark in last 2 generations of management period
                       ppnCURecover = apply(na.omit(counterLateUpperBM), 1, mean), #proportion of CUs above upper benchmark in last generations of management period
                       ppnCUStable = apply(na.omit(counterLateLowerBM), 1, mean), #proportion of CUs above lower benchmark in last generations of management period
                       ppnCUExtinct = ppnCUsExtinct[nYears, ], #proportion of CUs extinct at end of simulation period
                       ppnCUExtant = (1 - ppnCUsExtinct[nYears, ]), #proportion of CUs EXTANT at end of simulation period
                       #ppnCUConstrained = apply(na.omit(ppnConstrained), 2,
                       #                         mean),
                       medSpawnersEarly = apply(sAg[(nPrime + 1):endEarly, ], 2, median),
                       medRecRYEarly = apply(recRYAg[(nPrime + 1):endEarly, ], 2, median),
                       medCatchEarly = apply(catchAg[(nPrime + 1):endEarly, ], 2, median) #median aggregate catch in first 2 generations of management period
  )
  fileName <- ifelse(variableCU == "TRUE", paste(cuNameOM, cuNameMP, "aggDat.csv", sep = "_"),
                     paste(nameOM, nameMP, "aggDat.csv", sep = "_"))
  write.csv(aggDat, file = paste(here("outputs/simData"), dirPath, fileName, sep = "/"), row.names = FALSE)


  # Create LRP data for output
  colnames(sAg)<-as.character(1:nTrials)
  sAg.dat<-as_tibble(sAg)
  sAg.dat<-sAg.dat %>% add_column(year=1:nYears)
  sAg.dat<-sAg.dat %>% pivot_longer(as.character(1:nTrials), names_to="iteration", values_to="sAg")

  colnames(ppnCUsLowerBM)<-as.character(1:nTrials)
  ppnCUs.dat<-as_tibble(ppnCUsLowerBM)
  ppnCUs.dat<-ppnCUs.dat %>% add_column(year=1:nYears)
  ppnCUs.dat<-ppnCUs.dat %>% pivot_longer(as.character(1:nTrials), names_to="iteration", values_to="ppnCUsLowerBM")

  LRP.dat <- sAg.dat %>% left_join(ppnCUs.dat)

  fileName <- ifelse(variableCU == "TRUE", paste(cuNameOM, cuNameMP, "lrpDat.csv", sep = "_"),
                     paste(nameOM, nameMP, "lrpDat.csv", sep = "_"))

  write.csv(LRP.dat, file = paste(here("outputs/simData"), dirPath, fileName, sep = "/"),
            row.names = FALSE)



  # Create CU spawner abundance data for output

  for(i in 1:nTrials) {

    spnDat.i<-as_tibble(spwnrArray[,,i])

    spnDat.i<-spnDat.i %>% add_column(year=1:nrow(spwnrArray)) %>% add_column(iteration=rep(i,nrow(spwnrArray)))


    spnDat_long.i <- spnDat.i %>% select(starts_with("V"),iteration, year) %>% pivot_longer(starts_with("V"),names_to="CU", values_to="spawners")
    spnDat_long.i$CU<-rep(1:nCU,length=nrow(spnDat_long.i))

    if (i == 1) spnDat<-spnDat_long.i
    if (i > 1) {
      spnDat<-bind_rows(spnDat,spnDat_long.i)
    }

  }

  fileName <- ifelse(variableCU == "TRUE", paste(cuNameOM, cuNameMP, "CUspwnDat.csv", sep = "_"),
                     paste(nameOM, nameMP, "CUspwnDat.csv", sep = "_"))

  write.csv(spnDat, file = paste(here("outputs/simData"), dirPath, fileName, sep = "/"),
            row.names = FALSE)


  } # end of genericRecoverySim()

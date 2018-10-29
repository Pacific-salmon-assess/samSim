#*************************************************************************************
# runModelRobustOM_synchSensitivity.R
# Date revised: ONGOING
# Inputs: recoverySim.R
# Outputs: pdf plots
# Explainer: Equivalent to runModelRobustOM_synch.R, but focused on sensitivity analyses
# so smaller range of OMs examined (i.e. only median levels sigma and correlations); 
# plots are stripped down as a result
#*************************************************************************************


# Check if required packages are installed and run
listOfPackages <- c("plyr", "here", "parallel", "doParallel", "foreach", 
                    "reshape2", "tidyr", "gsl", "tictoc", "stringr", "dplyr",
                    "synchrony", "zoo", "Rcpp", "RcppArmadillo", "sn",
                    "sensitivity", "mvtnorm", "forcats", "ggpubr")

newPackages <- listOfPackages[!(listOfPackages %in% 
                                  installed.packages()[ , "Package"])]
if(length(newPackages)) install.packages(newPackages)
lapply(listOfPackages, require, character.only = TRUE)

source(here("scripts/func/postProcessing.R"))
source(here("scripts/recoverySim.R"))
source(here("scripts/func/simUtilityFunc.R"))

simPar <- read.csv(here("data/opModelScenarios/fraserOMInputs_varyCorrSensitivity.csv"), 
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/fraserCUpars.csv"), stringsAsFactors = F)
srDat <- read.csv(here("data/fraserRecDatTrim.csv"), stringsAsFactors = F)
catchDat <- read.csv(here("data/fraserCatchDatTrim.csv"), stringsAsFactors = F)
ricPars <- read.csv(here("data/fraserDat/rickerMCMCPars.csv"), 
                    stringsAsFactors = F)
larkPars <- read.csv(here("data/fraserDat/larkinMCMCPars.csv"), 
                     stringsAsFactors = F)
tamFRP <- read.csv(here("data/fraserDat/tamRefPts.csv"), 
                   stringsAsFactors = F)
cuCustomCorrMat <- read.csv(here("data/fraserDat/prodCorrMatrix.csv"), 
                            stringsAsFactors=F)


### SET UP MODEL RUN -----------------------------------------------------

## Define simulations to be run
nTrials <- 200

## General robustness runs
simParTrim <- subset(simPar,
                     scenario == "refSensitivity" | scenario == "ageTau" |
                       scenario == "ouSig" |  scenario == "forecastSig" |
                       scenario == "enRouteSig"
)
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), sep = "_"))
# 
# recoverySim(simParTrim[1, ], cuPar, catchDat = catchDat, srDat = srDat, variableCU = FALSE,
#                   ricPars, larkPars = larkPars, tamFRP = tamFRP, cuCustomCorrMat = cuCustomCorrMat,
#                   dirName = "test", nTrials = 5, multipleMPs = FALSE)

for (i in seq_along(dirNames)) {
  dirName <- dirNames[i]
  d <- subset(simParTrim, scenario == scenNames[i])
  simsToRun <- split(d, seq(nrow(d)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 3) #save two cores
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(MASS),
                     library(here),
                     library(sensitivity),
                     library(mvtnorm),
                     library(scales), #shaded colors for figs
                     library(here), #use this package so wd are common across different computers
                     library(synchrony),
                     library(zoo), #synch and zoo used to calculate rolling estimates of synchrony
                     library(viridis), #color blind gradient palette
                     library(ggplot2),
                     library(gsl), #to calculate exact estimate of MSY following Scheuerell 2016 PeerJ
                     library(dplyr),
                     library(Rcpp),
                     library(RcppArmadillo),
                     library(sn)))
  if (simsToRun[[1]]$species == "sockeye") {
    clusterExport(cl, c("simsToRun", "recoverySim", "cuPar", "dirName", "nTrials", "catchDat", "srDat",
                        "ricPars", "dirName", "larkPars", "tamFRP", "cuCustomCorrMat"), envir = environment()) #export custom function and objects
    tic("run in parallel")
    parLapply(cl, simsToRun, function(x) {
      recoverySim(x, cuPar, catchDat = catchDat, srDat = srDat, variableCU = FALSE, 
                  ricPars, larkPars = larkPars, tamFRP = tamFRP, cuCustomCorrMat = cuCustomCorrMat, 
                  dirName = dirName, nTrials = nTrials, multipleMPs = FALSE)
    })
    stopCluster(cl) #end cluster
    toc()
  }
  if (simsToRun[[1]]$species == "chum") {
    clusterExport(cl, c("simsToRun","recoverySim","cuPar","dirName","nTrials","catchDat","srDat",
                        "ricPars","dirName"), envir = environment()) #export custom function and objects
    tic("run in parallel")
    parLapply(cl, simsToRun, function(x) {
      recoverySim(x, cuPar, catchDat = catchDat, srDat = srDat, variableCU = FALSE, ricPars, 
                  larkPars = NULL, tamFRP = NULL, dirName = dirName, nTrials = nTrials, 
                  multipleMPs = FALSE)
    })
    stopCluster(cl) #end cluster
    toc()
  }
}


#_________________________________________________________________________
## Modified version of multi-OM grouped box plots
axisSize = 15; dotSize = 3.5; lineSize = 0.8; legendSize = 14

vars <- c("medRecRY", "ppnCULower", "ppnCUUpper", "ppnCUExtant",
          "medCatch", "ppnFisheriesOpen", "ppnYrsHighCatch", "stabilityCatch")
plotDat = NULL
for (h in seq_along(dirNames)) {
  agList <- genOutputList(dirNames[h], agg = TRUE)
  singleScen = NULL
  for (i in seq_along(vars)) {
    dum <- data.frame(scen = as.factor(rep(scenNames[h], length.out =  length(agList))),
                      var = rep(vars[i], length.out = length(agList)),
                      om = as.factor(sapply(agList, function(x) unique(x$opMod))),
                      avg = sapply(agList, function(x) median(x[,vars[i]])),
                      lowQ = sapply(agList, function(x) qLow(x[,vars[i]])),
                      highQ = sapply(agList, function(x) qHigh(x[,vars[i]])),
                      row.names = NULL
    )
    singleScen <- rbind(singleScen, dum)
  }
  plotDat <- rbind(plotDat, singleScen) #merge multiple scenarios into one dataframe
}
colPal <- c("black", "orange", "blue")
names(colPal) <- levels(plotDat$om)

# Plot
consVars <- c("medRecRY", "ppnCULower", "ppnCUUpper", "ppnCUExtant") 
consYLabs <- c("Recruit\nAbundance", "Prop. CUs\nLower", "Prop. CUs\nUpper", "Prop. CUs\nExtant")
consPlots <- lapply(seq_along(consVars), function(i) {
  temp <- plotDat %>% 
    filter(var == consVars[i])
  q <- ggplot(temp, aes(x = scen, y = avg, ymin = lowQ, ymax = highQ, 
                        color = om)) +
    labs(x = "Scenario", y = consYLabs[i], color = "Operating\nModel") +
    geom_pointrange(fatten = dotSize, size = lineSize, position = 
                      position_dodge(width = 0.5)) +
    scale_x_discrete(labels = c("refSensitivity" = "Reference",
                                "ageTau" = "Maturation\nAge", "ouSig" = 
                                  "Outcome\nUncertainty", "forecastSig" = 
                                  "Forecast\nError", "enRouteSig" = "En Route\nMortality")) +
    scale_colour_manual(values = colPal) 
  if (i == 1) {
    q <- q + theme_sleekX(position = "top") 
  } 
  if (i == 2 | i == 3) {
    q <- q + theme_sleekX(position = "mid")
  }
  if (i == 4) {
    q <- q + theme_sleekX(position = "bottom")
  }
  return(q)
})

catchVars <- c("medCatch", "stabilityCatch", "ppnFisheriesOpen", "ppnYrsHighCatch")
catchYLabs <- c("Catch\nAbundance", "Catch Stability", "Prop.\nFisheries Open", "Prop. Years\nHigher Catch")
catchPlots <- lapply(seq_along(catchVars), function(i) {
  temp <- plotDat %>% 
    filter(var == catchVars[i])
  q <- ggplot(temp, aes(x = scen, y = avg, ymin = lowQ, ymax = highQ, 
                        color = om)) +
    labs(x = "Scenario", y = catchYLabs[i], color = "Operating\nModel") +
    geom_pointrange(fatten = dotSize, size = lineSize, position = 
                      position_dodge(width = 0.5)) +
    scale_x_discrete(labels = c("refSensitivity" = "Reference",
                                "ageTau" = "Maturation\nAge", "ouSig" = 
                                  "Outcome\nUncertainty", "forecastSig" = 
                                  "Forecast\nError", "enRouteSig" = "En Route\nMortality")) +
    scale_colour_manual(values = colPal) 
  if (i == 1) {
    q <- q + theme_sleekX(position = "top")
  }
  if (i == 2 | i == 3) {
    q <- q + theme_sleekX(position = "mid")
  }
  if (i == 4) {
    q <- q + theme_sleekX(position = "bottom")
  }
  return(q)
})

png(file = paste(here("/outputs/summaryFigs/synchTrials/sens/consGroupedPlots.png"),
                 sep = ""), 
    height = 6.5, width = 7, units = "in", res = 150)
ggarrange(consPlots[[1]], consPlots[[2]], consPlots[[3]], consPlots[[4]], 
          ncol = 1, nrow = 4, common.legend = TRUE, legend = "right", 
          align = "v", heights = c(1,1,1,1.3))
dev.off()
png(file = paste(here("/outputs/summaryFigs/synchTrials/sens/catchGroupedPlots.png"),
                 sep = ""), 
    height = 6.5, width = 7, units = "in", res = 150)
ggarrange(catchPlots[[1]], catchPlots[[2]], catchPlots[[3]], catchPlots[[4]],
          ncol = 1, nrow = 4, common.legend = TRUE, legend = "right", 
          align = "v", heights = c(1,1,1,1.3))
dev.off()


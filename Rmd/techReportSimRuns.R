#******************************************************************************
# techReportSimRuns.R
# Date revised: June 19, 2019
# Explainer: Contains simulation runs to generate figures for SPERA technical
# report
#******************************************************************************


listOfPackages <- c("here", "parallel", "doParallel", "foreach", "viridis",
                    "tidyverse", "tictoc", "samSim")
lapply(listOfPackages, library, character.only = TRUE)

simPar <- read.csv(here("data", "manProcScenarios",
                        "fraserMPInputs_varyAllocation.csv"),
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/fraserDat/fraserCUpars.csv"), stringsAsFactors = F)
srDat <- read.csv(here("data/fraserDat/fraserRecDatTrim.csv"), stringsAsFactors = F)
catchDat <- read.csv(here("data/fraserDat/fraserCatchDatTrim.csv"), stringsAsFactors = F)
ricPars <- read.csv(here("data/trimRecursiveRickerMCMCPars.csv"), stringsAsFactors = F)
larkPars <- read.csv(here("data/trimRecursiveLarkinMCMCPars.csv"),
                     stringsAsFactors = F)
tamFRP <- read.csv(here("data/tamRefPts.csv"), stringsAsFactors = F)

## Define simulations to be run
nTrials <- 150

for(k in 1:nrow(simPar)) {
  if(is.na(simPar$nameMP[k])) {
    simPar$nameMP[k] <- paste(simPar$propMixHigh[k], simPar$singleHCR[k], "_",
                              simPar$canER[k], simPar$harvContRule[k], sep = "")
  }
}

simParTrim <- simPar
# %>%
#     filter(keyVar == "ppnMix",
#            scenario == "fixed0.7Retro_enRoute")
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species),
                                                sep = "_"))
agVarsToPlot <-  c("medRecRY", "medSpawners", "ppnCUUpper", "ppnCUExtinct",
                   "medCatch")

#------------------------------------------------------------------------------

## Run model
for (i in seq_along(dirNames)) {
  dirName <- dirNames[i]
  d <- subset(simParTrim, scenario == scenNames[i])
  simsToRun <- split(d, seq(nrow(d)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 1) #save two cores
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(MASS),
                     library(here),
                     library(sensitivity),
                     library(mvtnorm),
                     library(scales), #shaded colors for figs
                     library(here),
                     library(synchrony),
                     library(zoo),
                     library(viridis), #color blind gradient palette
                     library(ggplot2),
                     library(gsl), #to calculate exact estimate of MS
                     library(dplyr),
                     library(Rcpp),
                     library(RcppArmadillo),
                     library(sn),
                     library(samSim)))
  clusterExport(cl, c("simsToRun","cuPar","dirName","nTrials",
                      "catchDat","srDat","ricPars","larkPars",
                      "tamFRP"), envir=environment())
  tic("run in parallel")
  parLapply(cl, simsToRun, function(x) {
    recoverySim(x, cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE,
                ricPars, larkPars=larkPars, tamFRP=tamFRP, dirName=dirName,
                nTrials=nTrials, makeSubDirs=TRUE, random = FALSE)
  })
  stopCluster(cl) #end cluster
  toc()
}

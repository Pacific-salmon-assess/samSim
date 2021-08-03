## Check if necessary packages are available and install if necessary
listOfPackages <- c("here", "parallel", "doParallel", "foreach",
                    "tidyverse", "dplyr", "tictoc", "samSim", "tibble")
newPackages <- listOfPackages[!(listOfPackages %in%
                                  installed.packages()[ , "Package"])]
if (length(newPackages)) {
  install.packages(newPackages)
}
lapply(listOfPackages, require, character.only = TRUE)


# When de-bugging, uncomment out file that you want to get into
#source("R/genericRecoverySimulator.r")
#source("R/stockRecruitModels.r")
#source("R/getSRPars_randomSamp.r")
#source("R/calcTAC_fixedER.r")


# Coho ====================================================================

## (1) Load relevant input data

  # Simulation run parameters describing different scenarios
  simPar <- read.csv(here("data", "IFCohoPars",
                        "cohoSimPars.csv"), stringsAsFactors = F)

  # CU-specific parameters
  cuPars <- read.csv(here("data", "IFCohoPars", "cohoCUPars.csv"),
                  stringsAsFactors=F)

  # Stock-recruit and catch data that are used to populate the simulation priming period
  srDat <- read.csv(here("data", "IFCohoPars", "cohoRecDatTrim.csv"),
                  stringsAsFactors=F)

  # Posterior MCMC samples for SR pars
  ricPars<-read.csv(here("data", "IFCohoPars", "cohoRickerSurv_mcmc.csv"),
                  stringsAsFactors=F)

  # Correlation matrix for between-CU correlation in recruitment residuals
  corMatrix <- read.csv(here("data", "IFCohoPars", "cohoCorrMat.csv"),
                      stringsAsFactors=F, header=F)

  ## Store relevant object names to help run simulation
  scenNames <- unique(simPar$scenario)
  dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),
                                                sep = "_"))

## (2) Check to ensure that package works with only a small number of scenarios ======
genericRecoverySim(simPar[1, ], cuPar=cuPars, srDat=srDat,
                   variableCU=FALSE, ricPars=ricPars,
                   cuCustomCorrMat = corMatrix,
                   nTrials=5, makeSubDirs=FALSE,
                   random=TRUE, outDir="outDir")


## (3) Run with more scenarios
genericRecoverySim(simPar[1, ], cuPar=cuPars, srDat=srDat,
                   variableCU=FALSE, ricPars=ricPars,
                   cuCustomCorrMat = corMatrix,
                   nTrials=5000, makeSubDirs=FALSE,
                   random=TRUE, outDir="outDir")






## Old code to run in parallel =============================================
# simsToRun <- split(simPar, seq(nrow(simPar)))
# Ncores <- detectCores()
# cl <- makeCluster(Ncores - 1) #save one core
# registerDoParallel(cl)
# clusterEvalQ(cl, c(library(samSim)))
#
# clusterExport(cl, c("simsToRun", "cuPars",
#                     "srDat", "corMatrix", "ricPars"), envir=environment())
#
#
# tic("run in parallel")
# parLapply(cl, simsToRun, function(x) {
#
#   genericRecoverySim(x, cuPar=cuPars, srDat=srDat,
#                      variableCU=FALSE, ricPars=ricPars,
#                      cuCustomCorrMat = corMatrix,
#                      nTrials=100, makeSubDirs=FALSE,
#                      random=TRUE, outDir="outDir")
# })
# stopCluster(cl) #end cluster
# toc()

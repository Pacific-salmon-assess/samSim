## Check if necessary packages are available and install if necessary
listOfPackages <- c("here", "parallel", "doParallel", "foreach",
                    "tidyverse", "dplyr", "tictoc", "samSim")
newPackages <- listOfPackages[!(listOfPackages %in%
                                  installed.packages()[ , "Package"])]
if (length(newPackages)) {
  install.packages(newPackages)
}
lapply(listOfPackages, require, character.only = TRUE)


source("R/genericRecoverySimulator.r")

# Coho ====================================================================

## Load relevant input data

# Simulation run parameters describing different scenarios
simPar <- read.csv(here("data", "IFCohoPars",
                        "cohoSimPar.csv"),
                   stringsAsFactors = F)
# CU-specific parameters
cuPar <- read.csv(here("data", "IFCohoPars", "cohoCUPars_rickerSurv.csv"),
                  stringsAsFactors=F)

# Stock-recruit and catch data that are used to populate the simulation priming
# period
srDat <- read.csv(here("data", "IFCohoPars", "cohoRecDatTrim.csv"),
                  stringsAsFactors=F)


ricPars<-read.csv(here("data", "IFCohoPars", "cohoRickerSurv_mcmc.csv"),
                  stringsAsFactors=F)


corMatrix <- read.csv(here("data", "IFCohoPars", "cohoCorrMat.csv"),
                      stringsAsFactors=F, header=F)

## Store relevant object names to help run simulation
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),
                                                sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)

genericRecoverySim(simPar[1, ], cuPar=cuPar, srDat=srDat,
                 variableCU=FALSE, ricPars=ricPars, cuCustomCorrMat = corMatrix,
                 dirName="test.co", nTrials=100, makeSubDirs=FALSE, random=FALSE)



genericRecoverySim(simPar[2, ], cuPar=cuPar, srDat=srDat,
                   variableCU=FALSE, ricPars=ricPars, cuCustomCorrMat = corMatrix,
                   dirName="test.co", nTrials=100, makeSubDirs=FALSE, random=FALSE)

genericRecoverySim(simPar[3, ], cuPar=cuPar, srDat=srDat,
                   variableCU=FALSE, ricPars=ricPars, cuCustomCorrMat = corMatrix,
                   dirName="test.co", nTrials=100, makeSubDirs=FALSE, random=FALSE)




## Define a larger number of simulations to be run (note still well below
## suggested number for stability)
nTrials <- 100

## Divide each scenario into a list element to pass to parLapply()
simsToRun <- split(simPar, seq(nrow(simPar)))
dirName <- "test.co.cluster"
Ncores <- detectCores()
cl <- makeCluster(Ncores - 1) #save one core
registerDoParallel(cl)
clusterEvalQ(cl, c(library(samSim)))

clusterExport(cl, c("simsToRun", "cuPar", "nTrials", "dirName",
                     "srDat", "ricPars", "corMatrix"), envir=environment())

tic("run in parallel")

parLapply(cl, simsToRun, function(x) {
  genericRecoverySim(x, cuPar=cuPar, srDat=srDat, ricPars=ricPars, variableCU=FALSE,
                     cuCustomCorrMat = corMatrix, dirName=dirName, nTrials=nTrials, makeSubDirs=FALSE,
                      random=FALSE)
})
stopCluster(cl) #end cluster
toc()



# Sockeye ===================


  # Read-in sockeye pars

  # Simulation run parameters describing different scenarios
  simPar.sk <- read.csv(here("data", "sockeyeTestPars",
                             "fraserMPInputs_exampleSimPar.csv"),
                        stringsAsFactors = F)
  # CU-specific parameters
  cuPar.sk <- read.csv(here("data", "sockeyeTestPars", "fraserCUpars.csv"),
                       stringsAsFactors=F)
  # Stock-recruit and catch data that are used to populate the simulation priming
  # period
  srDat.sk <- read.csv(here("data", "sockeyeTestPars", "fraserRecDatTrim.csv"),
                       stringsAsFactors=F)
  catchDat.sk <- read.csv(here("data", "sockeyeTestPars", "fraserCatchDatTrim.csv"),
                          stringsAsFactors=F)
  # Posterior values of  CU-specific stock recruitment parameters for Ricker and
  # Larkin models; when available, passed and used to calculate alpha, beta and
  # sigma parameters rather than passing point values
  ricPars.sk <- read.csv(here("data", "sockeyeTestPars", "pooledRickerMCMCPars.csv"),
                         stringsAsFactors=F)
  larkPars.sk <- read.csv(here("data", "sockeyeTestPars", "pooledLarkinMCMCPars.csv"),
                          stringsAsFactors=F)



## Store relevant object names to help run simulation
scenNames.sk <- unique(simPar.sk$scenario)
dirNames.sk <- sapply(scenNames.sk, function(x) paste(x, unique(simPar.sk$species),
                                                      sep = "_"))

simpleRecoverySim(simPar.sk[1, ], cuPar=cuPar.sk, catchDat=catchDat.sk, srDat=srDat.sk,
                  variableCU=FALSE, ricPars=ricPars.sk,larkPars=larkPars.sk,
                  dirName="test.sk", nTrials=2, makeSubDirs=FALSE, random=FALSE)







#### Rest of example ============================

## Define a larger number of simulations to be run (note still well below
## suggested number for stability)
nTrials <- 50

## Divide each scenario into a list element to pass to parLapply()
simsToRun <- split(simPar, seq(nrow(simPar)))
dirName <- "example"
Ncores <- detectCores()
cl <- makeCluster(Ncores - 1) #save one core
registerDoParallel(cl)
clusterEvalQ(cl, c(library(samSim)))
# clusterExport(cl, c("simsToRun", "cuPar", "nTrials", "dirName",
#                     "catchDat", "srDat", "ricPars", "larkPars",
#                     "tamFRP"), envir=environment())
clusterExport(cl, c("simsToRun", "cuPar", "nTrials", "dirName",
                    "catchDat", "srDat", "ricPars",
                    "tamFRP"), envir=environment())

tic("run in parallel")
# parLapply(cl, simsToRun, function(x) {
#   recoverySim(x, cuPar=cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE,
#               ricPars=ricPars, larkPars=larkPars, tamFRP=tamFRP,
#               dirName=dirName, nTrials=nTrials, makeSubDirs=FALSE,
#               random=FALSE)
parLapply(cl, simsToRun, function(x) {
  recoverySim(x, cuPar=cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE,
              ricPars=ricPars, tamFRP=tamFRP,
              dirName=dirName, nTrials=nTrials, makeSubDirs=FALSE,
              random=FALSE)
})
stopCluster(cl) #end cluster
toc()

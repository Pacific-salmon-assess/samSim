#******************************************************************************
# runModelVaryMP_ppnMix.R
# Date revised: ONGOING
# Inputs: recoverySim.R
# Outputs: pdf plots, matrices for plotting
# Explainer: Runs closed loop simulation model with different OMs and MPs based
# on input csv files; structured so that MPs vary within a scenario, not OMs
#******************************************************************************


# Check if required packages are installed and run
listOfPackages <- c("plyr", "here", "parallel", "doParallel", "foreach", 
                    "reshape2", "tidyr", "gsl", "tictoc", "stringr", "dplyr",
                    "synchrony", "zoo", "Rcpp", "RcppArmadillo", "sn", 
                    "sensitivity", "mvtnorm", "devtools")
newPackages <- listOfPackages[!(listOfPackages %in% 
                                  installed.packages()[ , "Package"])]
if (length(newPackages)) {
  install.packages(newPackages)
}
lapply(listOfPackages, require, character.only = TRUE)

here <- here::here()

source(here("scripts/func/postProcessing.R"))
source(here("scripts/func/simUtilityFunc.R"))
source(here("scripts/recoverySim.R"))


simPar <- read.csv(here("data/manProcScenarios/fraserMPInputs_varyMixPpn.csv"), 
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/fraserCUpars.csv"), stringsAsFactors=F)
srDat <- read.csv(here("data/fraserRecDatTrim.csv"), stringsAsFactors=F)
catchDat <- read.csv(here("data/fraserCatchDatTrim.csv"), stringsAsFactors=F)
ricPars <- read.csv(here("data/fraserDat/rickerMCMCPars.csv"), stringsAsFactors=F)
larkPars <- read.csv(here("data/fraserDat/larkinMCMCPars.csv"), stringsAsFactors=F)
tamFRP <- read.csv(here("data/fraserDat/tamRefPts.csv"), stringsAsFactors=F)


### SET UP MODEL RUN ----------------------------------------------------------

## Define simulations to be run
nTrials <- 75

simParTrim <- simPar
# simParTrim <- subset(simPar, scenario == "ref" | scenario == "constrain" |
#                        scenario == "singleHCR" | scenario == "singleHCR_constrain" |
#                        scenario == "moveTAC")
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), 
                                                sep = "_"))
# recoverySim(simParTrim[7,], cuPar, catchDat = catchDat, srDat = srDat,
#             variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP,
#             cuCustomCorrMat = cuCustomCorrMat, dirName = "Test", nTrials = 2,
#             multipleMPs = FALSE)
# d <- subset(simParTrim, scenario == scenNames[1])
# simsToRun <- split(d, seq(nrow(d)))
# lapply(simsToRun, function(x) {
#   recoverySim(x, cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE,
#               ricPars, larkPars=larkPars, tamFRP=tamFRP, dirName=dirNames[1],
#               nTrials=2, multipleMPs=TRUE)
# })

for (i in seq_along(dirNames)) {
  dirName <- dirNames[i]
  d <- subset(simParTrim, scenario == scenNames[i])
  simsToRun <- split(d, seq(nrow(d)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 2) #save two cores
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
                  library(sn)))
  clusterExport(cl, c("simsToRun","recoverySim","cuPar","dirName","nTrials",
                      "catchDat","srDat","ricPars","dirName","larkPars",
                      "tamFRP"), envir=environment()) 
  tic("run in parallel")
  parLapply(cl, simsToRun, function(x) {
    recoverySim(x, cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE, 
                ricPars, larkPars=larkPars, tamFRP=tamFRP, dirName=dirName, 
                nTrials=nTrials, multipleMPs=TRUE)
    })
  stopCluster(cl) #end cluster
  toc()
}


### GENERATE OUTPUT FIGS ------------------------------------------------------
source(here("scripts/func/postProcessing.R"))

## Read in data from directories and convert CU-lists into dataframes then plot
cuDatF <- buildDataCU(dirNames = dirNames, 
                      cuVars = c("medSpawners", "medCatch", "ppnYrsUpper", 
                                 "ppnYrsLower", "medTotalER", 
                                 "medEarlyS", "medEarlyTotalCatch"),
                      keyVarName = "ppnMixed", 
                      selectedCUs = c("E.St", "Bwrn", "Chlk", "L.Sh", "Clts", 
                                      "Hrrs"))
cuDatF$muName <- factor(cuDatF$muName, levels(cuDatF$muName)[c(1, 2, 4, 3)])
# cuDatF$om <- factor(cuDatF$om, levels(cuDatF$om)[c(2, 1)])
cuDat <- cuDatF 
# %>% 
#   mutate(byShape = case_when( #add argument for multiple allocation strategies 
#     ppnMixed == "flex" ~ "flex",
#     ppnMixed != "flex" ~ "static"
#   ))

# trade off plots
consVec <- c("medSpawners", "medEarlyS")
yLabVec <- c("Median Spawner Abundance", "Median Early Spawner Abundance")
catchVec <- c("medCatch", "medEarlyTotalCatch")
xLabVec <- c("Median Catch", "Median Early Catch")

toPlotList <- lapply(seq_along(consVec), function (h) {
  p <- plotCUTradeoff(cuDat, consVar = consVec[h], catchVar = catchVec[h], 
                      facet = "mu", panel = "om", showUncertainty = TRUE,
                      xLab = xLabVec[h], yLab = yLabVec[h], 
                      legendLab = "Proportion Mixed\nStock TAC",
                      lineSize = 0.3)
  return(p)
})

pdf(here("outputs/varyAllocation/allSingleScenarios/Fig1_CUTradeoffs.pdf"), 
    height = 6.5, width = 8.5)
sapply(toPlotList, function (h) 
  sapply(h, function(i) print(i))
)
dev.off()

# dot plots
varVec <- c("medSpawners", "medCatch", "ppnYrsLower")
yLabVec <- c("Median Spawner Abundance", "Median Catch Abundance", 
             "Proportion Above Lower BM")
cuDatTrim <- cuDat %>%
  filter(om %in% c("ref", "singleHCR_constrain", "moveTAC", "singConsERBefore")) %>% 
  mutate(om = factor(om))
dotPlotList <- lapply(seq_along(varVec), function(x) {
  p <- plotCUDot(cuDatTrim, plotVar = varVec[x], group = "om", 
                 legendLab = "Operating\nModel", 
                 xLab = "Propotion of TAC in Mixed Fishery", 
                 yLab = yLabVec[x], axisSize = 14, dotSize = 2.5, 
                 lineSize = 0.9, legendSize = 14)
  return(p)
})

pdf(here("outputs/varyAllocation/allSingleScenarios/Fig2_CUDotplots.pdf"), 
    height = 6.5,
    width = 9.5)
sapply(dotPlotList, function(x) print(x))
dev.off()


## Read in data from directories and convert aggregate lists into df then plot
agDatF <- buildDataAgg(dirNames, agVars =  c("medSpawners", "medCatch", 
                                             "ppnCULower", "ppnCUExtinct", 
                                             "ppnCUStable", "ppnYrsHighCatch", 
                                             "ppnFisheriesOpen", "medCatchEarly",
                                             "medSpawnersEarly", "medRecRYEarly"), 
                       keyVarName = "ppnMixed")
# agDatF$om <- factor(agDatF$om, levels(agDatF$om)[c(2, 1, 3)])
agDat <- agDatF #%>% 
  # filter(mp == "TAM", !plotOrder == "6")

# trade off plots
consVec <- c("medSpawners", "medSpawnersEarly")
yLabVec <- c("Median Spawner Abundance", "Median Early Spawner Abundance")
catchVec <- c("medCatch", "medCatchEarly")
xLabVec <- c("Median Catch", "Median Early Catch")
toAgPlotList <- lapply(seq_along(consVec), function(x) {
  p <- plotAgTradeoff(agDat, consVar = consVec[x], catchVar = catchVec[x],
                      facet = "om", showUncertainty = TRUE,
                      xLab = xLabVec[x], yLab = yLabVec[x], 
                      legendLab = "Propotion\nTAC Mix\nFishery",
                      axisSize = 14, dotSize = 4, legendSize = 14,
                      lineSize = 0.5, freeY = FALSE)
  return(p)
})


pdf(here("outputs/varyAllocation/allSingleScenarios/Fig3_AgTradeoffs.pdf"), 
    height = 4.5, width = 6.5)
sapply(toAgPlotList, function(x) print(x))
dev.off()

agDatTrim <- agDat %>%
  filter(om %in% c("ref", "singleHCR_constrain", "moveTAC", "singConsERBefore")) %>% 
  mutate(om = factor(om))

pdf(here("outputs/varyAllocation/allSingleScenarios/Fig4_AgDotPlot.pdf"), 
    height = 4.5, width = 8.5)
plotAgDot(agDatTrim, group = "om", legendLab = "Operating\nModel",
          xLab = "Proportion of TAC in Mixed Fishery", yLab = "", axisSize = 14,
          dotSize = 4, lineSize = 1, legendSize = 14)
dev.off()


### Radar plots
dat <- buildDataAgg(dirNames, agVars =  c("ppnCUExtant", "ppnYrsHighCatch", 
                                          "ppnCULower", "ppnFisheriesOpen",
                                             "ppnCUStable"), 
                       keyVarName = "ppnMixed") %>% 
  filter(om == "staticTAC", mp == "TAM")

pdf(here("outputs/varyAllocation/shiftTAC/Fig5b_AgRadarPlot_staticTAC.pdf"), 
    height = 6, width = 7.5)
plotRadar(dat, xLab = c("Extant", "Years\nHigh Catch", "CUs\nLower BM", 
                        "Fisheries\n Open", "CUs\nStable"),
          plotVars = NULL, groupingVar = NULL, cu = FALSE, 
          legendLab = "Proportion\nTAC in\nMixed Catch", axisSize = 13)
dev.off()

#______________________________________________________________________________
### Pull trial data and look at relationship between median mix catch and 
## years constrained
subDirs <- list.dirs(path = paste(here("outputs/simData"), dirNames[i], sep = "/"), full.names = FALSE,
                     recursive = FALSE) #alternatively ID OMs based on prespecified directory

for (j in seq_along(subDirs)) {
  cuList <- genOutputList(dirNames[i], subDirs[j], agg = FALSE)
}
dum <- cuList[[3]]
plot(as.vector(dum$medForegoneCatch) ~ as.vector(dum$medCatch))


### Pull trial data and look at catch trends through time
agList <- genOutputList(dirNames[i], subDirs[1], agg = TRUE, aggTS = TRUE)

singDat2 <- agList[[2]][["Ag Catch"]]
mixDat2 <- agList[[6]][["Ag Catch"]]
par(mfrow = c(2, 4)) 
for (i in c(3, 27, 91, 18)) {
  plot(singDat[ , i], type = "l", ylab = "S catch")
  plot(mixDat[ , i], type = "l", ylab = "M catch")
}

s <- apply(singDat, 2, median)
m <- apply(mixDat, 2, median)
plot(s ~ m)

singDat <- agList[[2]][["Ag Spawners"]]
mixDat <- agList[[6]][["Ag Spawners"]]
par(mfrow = c(2, 4)) 
for (i in c(3, 27, 91, 18)) {
  plot(singDat[60:100, i], type = "l", ylab = "S catch")
  lines(singDat2[60:100, i], col = "red")
  plot(mixDat[60:100, i], type = "l", ylab = "M catch")
  lines(singDat2[60:100, i], col = "red")
}
par(mfrow = c(1,1))
s <- apply(singDat[60:100,], 2, sd)
m <- apply(mixDat[60:100,], 2, sd)
plot(s ~ m)
abline(1,1)

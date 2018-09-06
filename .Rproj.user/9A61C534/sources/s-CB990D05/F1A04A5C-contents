#' Make a list of output data
#'
#' This function is used to generate lists of summary data based on the outputs of recoverySim.
#'
#' This is a function specific to the samSim package. Generally a directory represents either a
#' single operating model/management procedure with subdirectories representing representing OM/MPs
#' nested within creating multiple unique scenarios.
#'
#' @importFrom here here
#'
#' @param dirName A character vector describing directory where output data saved.
#' @param subDirName A character vector describing subdirectory where output data saved.
#'  Defaults to NULL.
#' @param selectedCUs A character vector of CUs that will be subsetted from aggregate. Defaults
#'  to NULL.
#' @param agg Is the generated list summarizing aggregate or CU-specific data? Defaults to TRUE.
#' @param aggTS If agg == TRUE, is the generated list summarizing aggregate time series data?
#'  Defaults to FALSE.
#'
#' @return Returns a list of subdirectory outputs within the directory.
#'
#' @examples
#'
#' @export

genOutputList <- function(dirName, subDirName = NULL, selectedCUs = NULL, agg = TRUE,
                          aggTS = FALSE) {
  dirPath <- ifelse(is.null(subDirName), dirName, paste(dirName, subDirName, sep="/"))

  if (agg == TRUE) {
    if (aggTS == TRUE) {
      arrayNames <- list.files(paste(here("outputs/simData"), dirPath, sep="/"), pattern="\\Series.RData$")
      aggList <- list()
      for(i in 1:length(arrayNames)){ #make list of lists!
        aggList[[i]] <- readRDS(paste(here("outputs/simData"), dirPath, arrayNames[i], sep="/"))
      }
      names(aggList) <- arrayNames
      return(aggList)
    } else {
      listNames <- list.files(paste(here("outputs/simData"), dirPath, sep="/"), pattern="*aggDat.csv")
      aggList <- list()
      for(i in 1:length(listNames)){
        aggList[[i]] <- read.csv(paste(here("outputs/simData"), dirPath, listNames[i], sep="/"))
      }
      names(aggList) <- listNames
      return(aggList)
    }
  }

  if (agg == FALSE) {
    dfNames <- list.files(paste(here("outputs/simData"), dirPath, sep="/"), pattern="\\cuDat.RData$")
    cuList <- list()
    newNames <- NULL
    for(i in 1:length(dfNames)){ #make list of lists!
      cuList[[i]] <- readRDS(paste(here("outputs/simData"), dirPath, dfNames[i], sep="/"))
      newName <- unlist(strsplit(dfNames[i], "_cuDat"))[1] #identify and assign truncated name to CU-specific list
      # newNames <- c(newNames, paste("fixed", newName, sep = ""))
      newNames <- c(newNames, newName)
    }
    names(cuList) <- newNames

    if (is.null(selectedCUs) == FALSE) { #subset CU list based on input vector
      cuNumbers <- which(cuList[[1]][["stkName"]] %in% selectedCUs)
      cuList <- lapply(cuList, function(lst) {
        tempList <- vector("list", length = length(lst))
        tempList[1:4] <- lst[c("opMod", "keyVar", "plotOrder", "hcr")]
        for (i in 5:length(lst)) {
          if (is.matrix(lst[[i]]) == TRUE){
            temp <- lst[[i]][, cuNumbers]
          }
          if (is.vector(lst[[i]]) == TRUE){
            temp <- lst[[i]][cuNumbers]
          }
          tempList[[i]] <- temp
        }
        names(tempList) <- names(lst)
        return(tempList)
      })
    }
    return(cuList)
  }
}

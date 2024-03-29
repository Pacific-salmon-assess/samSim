% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDblAxis.R
\name{plotContTradeOffs}
\alias{plotContTradeOffs}
\title{Double y-axis Walters plot}
\usage{
plotContTradeOffs(agDat, keyVar = "expRate", double = TRUE)
}
\arguments{
\item{agDat}{Dataframe generated by \code{buildDataAgg}.}

\item{keyVar}{A character value corresponding to the harvest control
rule-esqure variable that is on the x-axis (defaults to "expRate") that is
used to label x-axis and subset data input data Should correspond to the
\code{keyVarName} in \code{buildDataAgg}.}

\item{double}{A logical that determines whether a double or single y-axis
plot is generated.}
}
\value{
Returns a base plot object.
}
\description{
This function generates a double y-axis plot of the variety popularized by C.
Walters to visualize tradeoffs between multiple performance metrics
simultaneously. It may be untoward to try place so many different variables
with divergents units in the same figure, but hey, people seem to want to see
them...
}

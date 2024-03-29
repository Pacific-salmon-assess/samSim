% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcBenchmarks.R
\name{sGenSolver}
\alias{sGenSolver}
\title{Solve for sGen}
\usage{
sGenSolver(theta, sMSY)
}
\arguments{
\item{theta}{A numeric vector of Ricker stock recruit parameters: alpha,
beta, and sigma.}

\item{sMSY}{A numeric of spawner abundance estimated to result in maximum
sustainable yield.}
}
\value{
Returns a numeric that is the spawner abundance that minimizes the
the log likelihood.
}
\description{
This function solves for sGen based on sMSY and the log-likelihood estimated
in sGenOptimum.
}
\examples{
Stock-recruit parameters approximate those of Fraser River sockeye salmon
Chilko CU.
alpha = 1.2
beta = 1.5
sigma = 0.8
theta = c(alpha, beta, sigma)
sMSY = 0.3
sGenSolver(theta, sMSY)
}

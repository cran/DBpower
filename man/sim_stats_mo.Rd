% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_power_indiv.R
\name{sim_stats_mo}
\alias{sim_stats_mo}
\title{Simulate power starting from individual-level data for multiple
outcomes setting.}
\usage{
sim_stats_mo(B, covY, xMat, gVec, alphaMat, gammaVec, checkpoint = FALSE)
}
\arguments{
\item{B}{Number of simulations.}

\item{covY}{Covariance matrix of outcomes.}

\item{xMat}{Design matrix of non-genetic covariates, n*p.}

\item{gVec}{n*1 vector of genotypes.}

\item{alphaMat}{p*K vector of regression coefficients for xMat.}

\item{gammaVec}{K*1 vector of regression coefficients for each outcome.}

\item{checkpoint}{Boolean, if true then print message every 50 simulations.}
}
\value{
A list with the elements:
\item{zMat}{Matrix of test statistics Z.}
\item{zVecGBJ}{Check on Z statistics, vector should match first row of zMat.}
\item{iMat}{Innovated statistics using correlation matrix under the null.}
}
\description{
Simulate power starting from individual-level data for multiple
outcomes setting.
}
\examples{
\dontrun{
covY <- matrix(data=0.3, nrow=10, ncol=10); diag(covY) <- 1
xMat <- cbind(1, rnorm(n = 1000), rbinom(n = 1000, size=1, prob=0.5))
gVec <- rbinom(n= 1000, size = 2, prob=0.3)
alphaMat <-matrix(data = 1, nrow=3, ncol=10)
gammaVec <- rep(0, 10)
sim_stats_mo(B=10000, covY = covY, xMat = xMat, gVec = gVec,
alphaMat = alphaMat, gammaVec = gammaVec)
}
}

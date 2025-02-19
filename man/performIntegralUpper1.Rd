% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_b1.R
\name{performIntegralUpper1}
\alias{performIntegralUpper1}
\title{Apply this function over 1:J to calculate each portion of the integral
we need for the upper bound.}
\usage{
performIntegralUpper1(j, muVec, sigMat, lBounds1, uBounds1, lBounds2, uBounds2)
}
\arguments{
\item{j}{Apply over this integer, the element that will be the largest in magnitude.}

\item{muVec}{Mean vector of test statistics under the alternative (assuming it's MVN).}

\item{sigMat}{Covariance matrix of test statistics under the alternative (assuming it's MVN).}

\item{lBounds1}{A 3J-2 vector of lower bounds for the first integral (see paper), bounds will be longer than for performIntegralLower1.}

\item{uBounds1}{A 3J-2 vector of upper bounds for the second integral (see paper).}

\item{lBounds2}{A 3J-2 vector of lower bounds for the first integral (see paper).}

\item{uBounds2}{A 3J-2 vector of upper bounds for the second integral (see paper).}
}
\value{
The value of the integration.
}
\description{
Apply this function over 1:J to calculate each portion of the integral
we need for the upper bound.
}

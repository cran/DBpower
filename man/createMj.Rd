% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_b1.R
\name{createMj}
\alias{createMj}
\title{Create the matrix that linearly transforms the vector of test statistics
into a quantity amenable for pmvnorm.}
\usage{
createMj(j, size)
}
\arguments{
\item{j}{The element of the vector that is the largest.}

\item{size}{The length of the set.}
}
\value{
The transformation matrix of dimension (2J-1)*(2J-1)
}
\description{
Create the matrix that linearly transforms the vector of test statistics
into a quantity amenable for pmvnorm.
}
\examples{
createMj(j=3, size=5)
}

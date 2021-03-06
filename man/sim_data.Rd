% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_data.R
\name{sim_data}
\alias{sim_data}
\title{Simulating data}
\usage{
sim_data(ntot, missing = FALSE, incorrect = FALSE, b_G = 0.6,
  b_X = c(0.02, 0.3, -0.12), b_X_Gsize = c(0.04, -2, -3),
  b_X_Gprob = c(-0.01, 0.5, 1), Sigma = diag(4), cond_cov = FALSE)
}
\arguments{
\item{ntot}{the number of observation to simulate}

\item{missing}{logical flag whether missing indicator should 
be used. Default is \code{FALSE}.}

\item{b_G}{strength parameter for linear association between the covariate "G"
and the outcome of interest "Y".}
}
\value{
simulated data matrix of size \code{ntot x 4}
}
\description{
Simulating data
}
\keyword{internal}

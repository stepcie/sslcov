% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smooth_sslCPP.R
\name{smooth_sslCPP}
\alias{smooth_sslCPP}
\title{Computes a smoothed projection using C++}
\usage{
smooth_sslCPP(ri, fi, fnew, wgt = NULL, bw = NULL, cdf_trans = TRUE, rsup)
}
\arguments{
\item{ri}{label data correlation}

\item{fi}{the labeled data}

\item{fnew}{the new data to be predicted}

\item{wgt}{optionnal weights (used for perturbations). Default is \code{NULL} in which case
no weighting is performed.}

\item{bw}{kernel bandwith. Default is \code{NULL} in which case the sd of the new data divided
by the total number of observation to the power 0.3 is used.}

\item{cdf_trans}{a logical flag indicating wether the smoothing should be
performed on the data transformed with their cdf. Default is \code{TRUE}.}

\item{rsup}{the supervised estimate of rhat}
}
\value{
a vector of length 3 containing:\itemize{
 \item rhat.ssl the semi-supervised estimation of rhat
 \item rhat.ssl.bc the semi-supervised estimation of rhat accounting for smoothing bias
 \item bw the value of the bandwith actually used
}
}
\description{
Smoothing predictor using Gaussian kernel using a C++ implementation
}
\details{
Smoothing over the CDF transformed data preven_learns some tail estimation issues
when the new data are subsequen_learnly large.
}
\keyword{internal}

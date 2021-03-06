% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extremeSampling.R
\name{extremeSampling}
\alias{extremeSampling}
\title{Extreme sampling function}
\usage{
extremeSampling(data, nn, surrogate_name = NULL)
}
\arguments{
\item{data}{the data}

\item{nn}{the number of observations to sample}

\item{surrogate_name}{a character string vector containing the name of the
columns from data containing the surrogate variables (at least 2)}
}
\value{
a list :
\code{extreme_index}
\code{weights}: the sampling weights associated with each sampled observation (inverse of their 
respective sampling probabilities). Can take only two values: 1 for extremes or wi0 for random
observations (wi0 > 1).
}
\description{
Extreme sampling function
}

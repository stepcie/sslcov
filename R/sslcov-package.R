#' sslcov
#'
#' Semi-Supervised Estimation Of Covariance
#'
#' This package implements semi-supervised estimation and testing of covariance 
#' and conditional covariance. The main function of the package is 
#' \code{\link{sslcov_test}}.
#'
#' \tabular{ll}{
#' Package: \tab sslcov\cr
#' Type: \tab Package\cr
#' Version: \tab sslcov.0.1.5\cr
#' Date: \tab 2017-08-03\cr
#' License: \tab \href{http://www.gnu.org/licenses/lgpl.txt}{MIT}\cr
#' }
#'
#' @author Boris P. Hejblum, Stephanie Chan, Tianxi Cai
#' --- Maintainer: Stephanie Chan
#'
#' @docType package
#' @name sslcov-package
#' @aliases sslcov
#' 
#' @references S Chan, BP Hejblum, A Chakrabortty, T Cai, Semi-Supervised 
#' Estimation of Covariance with Application to Phenome-wide Association Studies 
#' with Electronic Medical Records Data, 2017, \emph{submitted}.
#'
#' @useDynLib sslcov, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
NULL

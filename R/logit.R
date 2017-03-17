#'expit
#'
#'expit function
#'
#'@param x a real number
#'
#'@return a real number on the probabilty scale (between \code{0} and \code{1})
#'
#'@export
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}


#'logit
#'
#'logit function
#'
#'@param p a number between \code{0} and \code{1} (e.g. a probability)
#'
#'@return a real number between \code{-Inf} and \code{Inf}
#'
#'@export
logit <- function(p){
  return(log(p/(1-p)))
}
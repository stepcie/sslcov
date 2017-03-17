#'expit
#'
#'@export
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}


#'logit
#'
#'@export
logit <- function(p){
  return(log(p/(1-p)))
}
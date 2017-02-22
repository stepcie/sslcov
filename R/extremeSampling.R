#' Extreme sampling function
#'
#' @param data the data
#'
#' @param nn the number of observations to sample
#'
#' @param surrogate_name a character string vector containing the name of the
#'columns from data containing the surrogate variables (at least 2)
#'
#'@export
extremeSampling <- function(data, nn, surrogate_name=NULL){

  ntot <- nrow(data) #TOCHECK
  nsamp <- nn

  stopifnot(length(surrogate_name)>1)
  stopifnot(floor(nsamp/8)==nsamp/8)

  Si <- data[, surrogate_name]

  rank1 <- rank(Si[, 1],na.last=NA)
  rank2 <- rank(rowSums(Si[, -1]),na.last=NA)

  ind_lowlow <-  order(rank1 + rank2, decreasing = FALSE,na.last=NA)[1:(nsamp/4)] #1/4
  ind_highhigh <- order(rank1 + rank2, decreasing = TRUE,na.last=NA)[1:(nsamp/4)] #1/4
  ind_lowhigh <- order(rank1 - rank2, decreasing = FALSE,na.last=NA)[1:(nsamp/8)] #1/8
  ind_highlow <- order(rank1 - rank2, decreasing = TRUE,na.last=NA)[1:(nsamp/8)] #1/8

  ind_extreme <- c(ind_lowlow, ind_highhigh, ind_lowhigh, ind_highlow)
  ind_notextreme <- setdiff(1:ntot, ind_extreme)
  ind_random <- sample(ind_notextreme, nsamp/4, replace=FALSE)
  ind_samp <-  c(ind_extreme, ind_random) # ask Stephanie about the union

  delta <-  1*(ind_samp %in% ind_extreme)
  vri <-  1*(ind_samp %in% ind_random)
  pi <-  (nsamp/4)/ntot
  wi  <-  delta + (1-delta)*vri/pi
  wgt <- c(wi, rep(1/pi, ntot-nsamp))

  wi <- wi/max(wi)
  wgt <- wgt/max(wgt)

  return(list("extreme_index"=ind_samp, "weights" = wi))
}

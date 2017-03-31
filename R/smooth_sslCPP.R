#' Computes a smoothed projection using C++
#'
#' Smoothing predictor using Gaussian kernel using a C++ implementation
#'
#' Smoothing over the CDF transformed data preven_learns some tail estimation issues
#' when the new data are subsequen_learnly large.
#'
#' @param ri label data correlation
#'
#' @param fi the labeled data
#'
#' @param fnew the new data to be predicted
#'
#' @param wgt optionnal weights (used for perturbations). Default is \code{NULL} in which case
#' no weighting is performed.
#'
#' @param bw kernel bandwith. Default is \code{NULL} in which case the sd of the new data divided
#' by the total number of observation to the power 0.3 is used.
#'
#' @param cdf_trans a logical flag indicating wether the smoothing should be
#' performed on the data transformed with their cdf. Default is \code{TRUE}.
#'
#' @param rsup the supervised estimate of rhat
#'
#' @importFrom stats sd
#'
#' @return a vector of length 3 containing:\itemize{
#'  \item rhat.ssl the semi-supervised estimation of rhat
#'  \item rhat.ssl.bc the semi-supervised estimation of rhat accounting for smoothing bias
#'  \item bw the value of the bandwith actually used
#' }
#'
#' @keywords internal
#'
#' @export
#'
smooth_sslCPP <- function(ri, fi, fnew, wgt = NULL, bw = NULL, cdf_trans=TRUE, rsup){

  n_learn <- length(ri)
  n_new <- length(fnew)
  nn <- length(fi)

  if(is.null(wgt)){
    wgti <- rep(1,n_learn)
    #wgt_j <- rep(1,n_new)
  }else{
    wgti <- wgt[(1:n_learn)]
    #wgt_j <- wgt[-(1:n_learn)]
  }

  if(cdf_trans){
    pall <- ecdf_cpp(c(fi, fnew), c(fi, fnew))
    #pi <- sum.I(c(fi, fnew), "<=", c(fi, fnew), Vi=wgt)
    pnew <-  pall[-(1:n_learn)]
    pi <-  pall[1:n_learn]
  }else{
    pi <-  fi
    pall <- c(fi,fnew)
  }

  if(is.null(bw)){
    bw <-  sd(pall)/nn^0.3
  }


  kij <- dnormC_multi(pi, pall, bw, Log = FALSE)
  rhat.num <- c(t(wgti*ri)%*%kij)
  rhat.den <- c(t(wgti)%*%kij)

  rhat.ssl <- mean(rhat.num/rhat.den)

  rhat.ssl.bc <- rhat.ssl - (sum(wgti*rhat.num[1:n_learn]/rhat.den[1:n_learn])/sum(wgti) - rsup)

  return(c(rhat.ssl, rhat.ssl.bc, bw, rhat.num/rhat.den))
}

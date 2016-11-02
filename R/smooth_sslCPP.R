#' Computes a smoothed projection
#'
#' Smoothing predictor using Gaussian kernel
#'
#' Smoothing over the CDF transformed data preven_learns some tail estimation issues
#' when the new data are subsequen_learnly large.
#'
#' @param ri label data correlation
#' @param fi the odl data
#' @param fnew the new data to be
#' @param wgt facultative weights (used for perturbations)
#' @param bw bandwith
#' @param cdf_trans a logical flag indicating wether the smoothing should be
#' performed on the data transformed with their cdf. Default is \code{TRUE}.
#' See Details.
#'
#'@return a
#'
#'@importFrom landpred VTM
#'
#'@export
smooth_sslCPP <- function(ri, fi, fnew, rsup, wgt=NULL, bw=NULL, cdf_trans=TRUE, condi=FALSE){

  n_learn <- length(ri)
  n_new <- length(fnew)
  nn <- length(fi)

  if(is.null(wgt)){
    wgt_i <- rep(1,n_learn)
    #wgt_j <- rep(1,n_new)
  }else{
    wgt_i <- wgt[(1:n_learn)]
    #wgt_j <- wgt[-(1:n_learn)]
  }

  if(cdf_trans){
    pi <- ecdf_cpp(c(fi, fnew), c(fi, fnew))
    #pi <- sum.I(c(fi, fnew),"<=" ,c(fi, fnew),Vi=wgt)
    pnew <-  pi[-(1:n_learn)]
    pi <-  pi[1:n_learn]
  }else{
    pi <-  fi
    pnew <- fnew
  }

  if(is.null(bw)){
    bw <-  sd(c(pi,pnew))/nn^0.3
  }

  kij <- dnormC_multi(pi, c(pi,pnew), bw, Log = FALSE)
  #ki <- dnormC_multi(pi, pnew, bw, Log = FALSE)
  #smthed <- c(t(wgt_i*ri)%*%ki)/c(t(wgt_i)%*%ki)
  rhat.num = c(t(wgt_i*ri)%*%kij); rhat.den = c(t(wgt_i)%*%kij)
  rhat.ssl = mean(rhat.num/rhat.den)
  #rhat.ssl.bc = rhat.ssl - (mean(rhat.num[1:n_learn]/rhat.den[1:n_learn])-rsup)
  rhat.ssl.bc = rhat.ssl - (mean(rhat.num[1:n_learn]/rhat.den[1:n_learn]*wgt_i)/mean(wgt_i)-rsup)

  #return(c("rhat"=mean(smthed*wgt_j,na.rm=T)/mean(wgt_j), bw))
  #return(c("rhat"=mean(smthed, na.rm=T), bw))
  return(c(rhat.ssl,rhat.ssl.bc,bw))
}

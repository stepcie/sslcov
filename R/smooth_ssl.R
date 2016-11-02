#' Computes a smoothed projection
#'
#' Smoothing predictor using Gaussian kernel
#'
#' Smoothing over the CDF transformed data preven_learns some tail estimation issues
#' when the new data are subsequen_learnly large.
#'
#' @param ri label data correlation
#' @param fi the labeled data
#' @param fnew the new data to be predicted
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
# smooth_ssl <- function(ri, fi, fnew, rsup, wgt=NULL, bw=NULL, cdf_trans=TRUE){
#
#   n_learn <- length(ri)
#   n_new <- length(fnew)
#   nn <- length(fi)
#
#   if(is.null(wgt)){
#     wgt_i <- rep(1,n_learn)
#   }else{
#     wgt_i <- wgt[1:n_learn]
#   }
#
#   if(cdf_trans){
#     pall <- ecdf(c(fi, fnew))(c(fi, fnew))
#     #pi <- sum.I(c(fi, fnew), FUN=">=", c(fi, fnew), wgt)/sum(wgt) #compute the quantiles => cdf
#     pnew <-  pall[-(1:n_learn)]
#     pi <-  pall[1:n_learn]
#   }else{
#     pi <-  fi
#     pall <- c(fi,fnew)
#   }
#
#   if(is.null(bw)){
#     bw <-  sd(pall)/nn^0.3
#   }
#
#   kij <- dnorm(pi-landpred::VTM(pall, n_learn), sd=bw)
#   rhat.num = c(t(wgt_i*ri)%*%kij); rhat.den = c(t(wgt_i)%*%kij)
#   rhat.ssl = mean(rhat.num/rhat.den)
#   rhat.ssl.bc = rhat.ssl - (mean(rhat.num[1:n_learn]/rhat.den[1:n_learn])-rsup)
#   return(c(rhat.ssl,rhat.ssl.bc,bw))
#
#   #return(c("rhat" = mean(smthed*wgt_j, na.rm=TRUE)/mean(wgt_j), bw))
#   #return(c("rhat"=mean(smthed, na.rm=T), bw))
# }
smooth_ssl <- function(ri, fi, fnew, wgti=NULL, bw=NULL, cdf_trans=TRUE, rsup){
  n_learn <- length(ri); n_new <- length(fnew); nn <- length(fi)
  if(is.null(wgti)){wgti = rep(1,nn)}
  if(cdf_trans){
    pall <- ecdf(c(fi, fnew))(c(fi, fnew)); pi <-  pall[1:n_learn]
  }else{
    pi <-  fi; pall <- c(fi,fnew)
  }
  if(is.null(bw)){bw <-  sd(pall)/nn^0.3}
  kij <- dnorm(pi-VTM(pall, n_learn), sd=bw); ## n_learn x n_all
  rhat.num = c(t(wgti*ri)%*%kij); rhat.den = c(t(wgti)%*%kij)
  rhat.ssl = mean(rhat.num/rhat.den)
  rhat.ssl.bc = rhat.ssl - (mean(rhat.num[1:n_learn]/rhat.den[1:n_learn]*wgti)/mean(wgti)-rsup)
  c(rhat.ssl,rhat.ssl.bc,bw)
}

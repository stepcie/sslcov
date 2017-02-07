#'Semi-supervised correlation estimation
#'
#'This function estimates the correlation between an outcome available only
#'for a small subset of the data and a covariate. The outcome is imputed to all
#'the data using a smoothed predictor learned thanks to a set of surrogate variables,
#'available for all the data.
#'
#'Smoothing over the CDF transformed data prevents some tail estimation issues when the new data are subsequently large.
#'
#'@param data the data. The first \code{nn} rows should be the labeled data, the
#'remaining rows should be the unlabeled data.
#'
#'@param nn the number of labeled data
#'
#'@param outcome_name a character string containing the name of the
#'column from data containing the partly missing outcome of interest
#'
#'@param covariate_name a character string containing the name of the
#'column from data containing the covariate to be related to the outcome of
#'interest
#'
#'@param surrogate_name a character string vector containing the name of the
#'column(s) from data containing the surrogate variable(s)
#'
#'@param bw the bandwidth to use
#'
#'@param cdf_trans a logical flag indicating wether the smoothing should be
#'performed on the data transformed with their cdf. Default is TRUE. See Details.
#'
#'@param ptb_nolabel a logical flag indicating whether accounting for the variation
#'due to unlabeled data. Default is \code{FALSE}
#'
#'@param weights a vector of weights in case a weighted version of the
#'correlation has to be computed. Default is \code{NULL}, in which case, no
#'additional weighting is done and regular perturbation is performed.
#'
#'@param ptb_beta logical flag indicating whether beta coefficient should be perturbed.
#'Dafault is \code{TRUE}.
#'
#'@param adjust_covariates_name optional vector of names of the covariates to adjust on during imputation and smoothing.
#'Default is \code{NULL}.
#'
#'@param do_interact logical flag indicating whether interactins between \code{x} and
#'covariates should be taken into account when imputing \code{y}. Default is \code{TRUE}.
#'
#'@return a list with the following elements:\itemize{
#'    \item rhat
#'    \item bw the bandwith used
#'    \item data_sup
#'    \item W_unlabel
#'    }
#'
#'
#'@seealso \code{\link{smooth_ssl}} \code{\link{r_hat}}
#'
#'@export
rhat_ptb_cond <- function(data, nn, outcome_name=NULL, covariate_name=NULL,
                     surrogate_name=NULL, bw, cdf_trans = TRUE, ptb_nolabel = FALSE,
                     weights=NULL, ptb_beta=TRUE, adjust_covariates_name=NULL,
                     do_interact=TRUE){

  NN <- nrow(data)
  stopifnot(length(weights)==nn)

  outcome_colnum <- which(colnames(data)==outcome_name)

  #make sure the data are ordered correctly:

  # compute the perturbation weights:
  if(ptb_nolabel){
    Vij <-  rbeta(NN, 0.5, 1.5)*4
  }else{
    Vij <- c(rbeta(nn, 0.5, 1.5)*4, rep(1,NN-nn))
  }

  if(is.null(weights)){
    weights <- rep(1, NN)
  }else{
    weights <- c(weights, rep(1, NN-nn))
  }
  Vij_w <- Vij*weights
  Vi <- Vij_w[1:nn]

  data_centered_ptb <- data[, covariate_name, drop=FALSE]# - mean(data[, covariate_name]*Vij)/mean(Vij)
  data_all_ptb <- cbind(data[, outcome_colnum, drop=FALSE], data_centered_ptb, data[, surrogate_name, drop=FALSE])
  if(!is.null(adjust_covariates_name)){
    data_all_ptb <- cbind(data_all_ptb, data[, adjust_covariates_name, drop=FALSE])
  }
  if(do_interact){
    data_interact_ptb <- data_centered_ptb[, covariate_name]*data[, surrogate_name, drop=FALSE]
    data_all_ptb <- cbind(data_all_ptb, data_interact_ptb)
  }
  outcome_colnum <- 1
  covariate_colnum <- 2

  ncoef <- ncol(data_all_ptb)
  W_unlabel_ptb <- data_all_ptb[-(1:nn), -outcome_colnum]

  data_sup_ptb <- data_all_ptb[1:nn, ]

  covariate_colnum <- which(colnames(data_all_ptb)==covariate_name)
  adjust_covariates_colnums <- which(colnames(data_all_ptb) %in% adjust_covariates_name)
  surrogate_colnums <- which(colnames(data_all_ptb) %in% surrogate_name)

  if(all.equal(floor(data_all_ptb[, covariate_name]), data_all_ptb[, covariate_name])!=TRUE){
    #cat("Covariates are not counts. They should be log(counts + 1) then.\n")
    covariate_counts <- exp(data_all_ptb[, covariate_name]) - 1
  }else{
    covariate_counts <- data_all_ptb[, covariate_name]
  }
  alpha_0 <- MASS::glm.nb(covariate_counts~data_all_ptb[, -c(outcome_colnum, covariate_colnum, surrogate_colnums)], weights=Vij)$coef
  pred_G <- exp(cbind(1, data_all_ptb[, -c(outcome_colnum, covariate_colnum, surrogate_colnums)])%*%matrix(alpha_0, ncol=1))
  cond_G_res <- (data_all_ptb[, covariate_name] - pred_G)

  #yi_cen_ptb <- data_sup_ptb[,1] - mean(data_sup_ptb[,1])
  ri_ptb <- lm(data_sup_ptb[, outcome_colnum]~data_sup_ptb[, adjust_covariates_colnums], weights=Vi)$residuals*cond_G_res[1:nn]
  rhat_ptb_sup <- mean(ri_ptb*Vi)/mean(Vi)

  yi_cen_ptb <- lm(data_sup_ptb[, outcome_colnum]~data_sup_ptb[, adjust_covariates_colnums], weights=Vi)$residuals
  if(ptb_beta){
    betptb <- lm(yi_cen_ptb ~ data_sup_ptb[, -outcome_colnum], weights=Vi)$coef[1:ncoef]
  }
  else{
    betptb <- lm(yi_cen_ptb ~ data_sup_ptb[,-outcome_colnum])$coef[1:ncoef]
  }
  if(length(which(is.na(betptb)))>0){
    betptb[which(is.na(betptb))] <- 0 #TODO
  }
  fi_ptb <- c(cbind(1, data_sup_ptb[, -outcome_colnum])[, 1:ncoef]%*%betptb)*cond_G_res[1:nn]
  fj_ptb <- c(cbind(1, W_unlabel_ptb)%*%betptb)*cond_G_res[-c(1:nn)]
  ptb_ssl <- smooth_sslCPP(ri=ri_ptb, fi=fi_ptb, fnew=fj_ptb, rsup=rhat_ptb_sup, wgt=Vij_w, bw=bw, cdf_trans=cdf_trans)
  #rptb.ssl = mean(Vj*npreg(bws=bw,txdat=fi_ptb,tydat=ri_ptb*Vi,exdat=fj_ptb)$mean/
  #                  npreg(bws=bw,txdat=fi_ptb,tydat=Vi,exdat=fj_ptb)$mean,na.rm=T)/mean(Vj)
  c("rhat_sup"=rhat_ptb_sup,"rhat_ns"=mean(c(fi_ptb,fj_ptb)), "rhat_ssl"=ptb_ssl[1], "rhat_ssl_bc"=ptb_ssl[2], "beta_lm_ptb"=betptb)
}

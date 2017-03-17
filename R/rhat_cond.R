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
#'column from data containing the covariate of interest to be related to the outcome of
#'interest
#'
#'@param surrogate_name a character string vector containing the name of the
#'column(s) from data containing the surrogate variable(s)
#'
#'@param bw the bandwidth to use
#'
#'@param cdf_trans a logical flag indicating wether the smoothing should be performed on the data
#'transformed with their cdf. Default is TRUE. See Details.
#'
#'@param weights a vector of weights in case a weighted version of the
#'correlation has to be computed. Default is \code{NULL}, in which case, no
#'weighting is done.
#'
#'@param adjust_covariates_name optional vector of names of the covariates to adjust on during imputation and smoothing.
#'Default is \code{NULL}.
#'
#'@param do_interact logical flag indicating whether interactins between \code{x} and
#'covariates should be taken into account when imputing \code{y}. Default is \code{TRUE}.
#'
#'@importFrom stats lm
#'
#'
#'@return a list with the following elements:\itemize{
#'    \item rhat
#'    \item bw the bandwith used
#'    \item data_sup
#'    \item W_unlabel
#'    }
#'
#'@seealso \code{\link{smooth_ssl}}
#'
#'@importFrom MASS glm.nb
#'@export
rhat_cond <- function(data, nn, outcome_name=NULL, covariate_name=NULL,
                 surrogate_name=NULL, bw=NULL, cdf_trans=TRUE,
                 weights=NULL, adjust_covariates_name=NULL, do_interact=TRUE){


  outcome_colnum <- which(colnames(data)==outcome_name)
  NN <- nrow(data)
  stopifnot(length(weights)==nn)

  if(is.null(weights)){
    Vij <- rep(1,NN)
  }else{
    Vij <-  c(weights, rep(1,NN-nn))
  }
  Vi <- Vij[1:nn]

  data_centered <- data[, covariate_name, drop=FALSE]# - mean(data[, covariate_name]*Vij, na.rm = TRUE)/mean(Vij) # center G with mean from the entire dataset
  data_all <- cbind(data[, outcome_colnum, drop=FALSE], data_centered, data[, surrogate_name, drop=FALSE])
  if(!is.null(adjust_covariates_name)){
    data_all <- cbind(data_all, data[, adjust_covariates_name, drop=FALSE])
  }
  if(do_interact){
    data_interact <- data_centered[, covariate_name]*data[, surrogate_name, drop=FALSE]
    data_all <- cbind(data_all, data_interact)
  }
  outcome_colnum <- 1

  ncoef <- ncol(data_all)

  W_label <- data_all[(1:nn), -outcome_colnum]
  W_unlabel <- data_all[-(1:nn), -outcome_colnum]
  data_sup <- data_all[1:nn,]


  #new for cond
  covariate_colnum <- which(colnames(data_all)==covariate_name)
  adjust_covariates_colnums <- which(colnames(data_all) %in% adjust_covariates_name)
  surrogate_colnums <- which(colnames(data_all) %in% surrogate_name)

  if(all.equal(floor(data_all[, covariate_name]), data_all[, covariate_name])!=TRUE){
    #cat("Covariates are not counts. They should be log(counts + 1) then.\n")
    covariate_counts <- exp(data_all[, covariate_name]) - 1
  }else{
    covariate_counts <- data_all[, covariate_name]
  }

  #SUP
  gamma_tilde <- MASS::glm.nb(covariate_counts[1:nn]~data_all[1:nn, -c(outcome_colnum, covariate_colnum, surrogate_colnums), drop=FALSE])$coef
  #TODO weights=Vi ?
  linear_predictor <- cbind(1, data_all[1:nn, -c(outcome_colnum, covariate_colnum, surrogate_colnums)])%*%matrix(gamma_tilde, ncol=1)
  # residuals are computed as exp(log(Y)-(linear_predictor +1))...
  pred_G_sup <- exp(linear_predictor)
  cond_G_res_sup <- (covariate_counts[1:nn] - pred_G_sup)

  #yi_cen <- data_sup[, 1] - mean(data_sup[, 1]*Vi)/mean(Vi)
  linearmodel_y_sup <- lm(data_sup[, outcome_colnum]~data_sup[, adjust_covariates_colnums], weights=Vi)
  yi_cen <- linearmodel_y_sup$residuals

  mu_y_tilde_i <- cbind(1, data_all[, adjust_covariates_colnums]) %*% linearmodel_y_sup$coef


  ri_hat <- yi_cen*cond_G_res_sup[1:nn]
  rhat_sup <- mean(ri_hat*Vi)/mean(Vi)


  #SSL
  beta_hat <- lm(data_sup[, outcome_colnum]~data_sup[, -outcome_colnum], weights = Vi)$coef
  if(length(which(is.na(beta_hat)))>0){
    warning("some betahat are NA...")
    beta_hat[which(is.na(beta_hat))] <- 0 #TODO
  }
  gamma_hat <- MASS::glm.nb(covariate_counts~data_all[, -c(outcome_colnum, covariate_colnum, surrogate_colnums), drop=FALSE])$coef
  pred_G <- exp(cbind(1, data_all[, -c(outcome_colnum,covariate_colnum, surrogate_colnums)])%*%matrix(gamma_hat, ncol=1))
  cond_G_res <- (covariate_counts - pred_G)
  #print(beta_hat)
  #betax <- cbind(1, data_sup[, -outcome_colnum])[, 1:ncoef]%*%beta_hat
  #plot(density(betax))
  fi_hat <- (c(cbind(1, W_label)%*%beta_hat) - mu_y_tilde_i[1:nn])*cond_G_res[1:nn]
  fj_hat <- (c(cbind(1,W_unlabel)%*%beta_hat) - mu_y_tilde_i[-c(1:nn)])*cond_G_res[-c(1:nn)]
  ##rhat_ssl = mean(npreg(bws=bw,txdat=fi_hat,tydat=ri_hat, exdat = fj_hat)$mean,na.rm=T)
  rhat_ssl_smres <- smooth_sslCPP(ri = ri_hat, fi = fi_hat, fnew = fj_hat, rsup = rhat_sup,
                            wgt = weights, bw = bw, cdf_trans = cdf_trans)

  mij_hat <- rhat_ssl_smres[4:length(rhat_ssl_smres)]
  bw <- rhat_ssl_smres[3]
  rhat_ssl_bc <- rhat_ssl_smres[2]
  rhat_ssl <- rhat_ssl_smres[1]
  return(list("rhat" = c("Supervised"=rhat_sup,"NoSmooth"=mean(c(fi_hat,fj_hat)), "SemiSupervised"=rhat_ssl,
                         "SemiSupervisedBC"=rhat_ssl_bc),
              "var" = c("Supervised"=sum(Vi*(ri_hat-rhat_sup)^2)/(nn*sum(Vi)),
                        "NoSmooth"=sum(Vi*(ri_hat-fi_hat)^2)/(nn*sum(Vi)),
                        "SemiSupervised"=sum(Vi*(ri_hat-mij_hat[1:nn])^2)/(nn*sum(Vi)),
                        "SemiSupervisedBC"=sum(Vi*(ri_hat-mij_hat[1:nn])^2)/(nn*sum(Vi))),
              "bw" = bw,
              "data_sup" = data_sup,
              "W_unlabel" = W_unlabel,
              "beta_lm" = beta_hat)
  )
}



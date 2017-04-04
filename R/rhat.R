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
#'@param cdf_trans a logical flag indicating wether the smoothing should be performed on the data
#'transformed with their cdf. Default is TRUE. See Details.
#'
#'@param weights a weighting vector of length \code{nn} in case a weighted version of the
#'correlation has to be computed. Default is \code{NULL}, in which case, no
#'weighting is done.
#'
#'@param adjust_covariates_name optional vector of names of the covariates to adjust on during imputation and smoothing.
#'Default is \code{NULL}
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
#'@export
rhat <- function(data, nn, outcome_name=NULL, covariate_name=NULL,
                 surrogate_name=NULL, bw=NULL, cdf_trans=TRUE,
                 weights=NULL, adjust_covariates_name=NULL, do_interact=TRUE){
  
  
  outcome_colnum <- which(colnames(data)==outcome_name)
  NN <- nrow(data)
  stopifnot(length(weights)==nn)
  
  if(is.null(weights)){
    wi0 <- NN/nn # sampling weight of a random obs 
    weights <- rep(wi0, nn)
  }else{
    wi0 <- max(weights) # sampling weight of a random obs (extremes sampling weights are 1)
  }
  Vij <- c(weights, rep(wi0, NN-nn))
  Vi <- Vij[1:nn]
  n0 <- sum(Vi!=1)
  
  # variance weights from sampling probabilities :
  pi <- rep(nn/n0, nn)
  pi[Vi==1] <- nn/NN
  
  # data processing
  data_centered <- data[, covariate_name, drop=FALSE] - mean(data[, covariate_name], na.rm = TRUE)#mean(data[, covariate_name]*Vij, na.rm = TRUE)/mean(Vij) # center G with mean from the entire dataset
  data_all <- cbind(data[, outcome_colnum], data_centered, data[, surrogate_name])
  if(!is.null(adjust_covariates_name)){
    data_all <- cbind(data_all, data[, adjust_covariates_name])
  }
  if(do_interact){
    data_interact <- data_centered[, covariate_name]*data[, surrogate_name, drop=FALSE]
    data_all <- cbind(data_all, data_interact)
  }
  outcome_colnum <- 1
  
  ncoef <- ncol(data_all)
  W_unlabel <- data_all[-(1:nn), -outcome_colnum]
  
  data_sup <- data_all[1:nn,]
  
  yi_cen <- data_sup[, 1] - mean(data_sup[, 1]*Vi)/mean(Vi)
  
  ri_hat <- yi_cen*data_sup[, covariate_name]
  rhat_sup <- mean(Vi*ri_hat)/mean(Vi)
  
  bethat <- lm(yi_cen~data_sup[, -outcome_colnum], weights = Vi)$coef[1:ncoef]
  if(length(which(is.na(bethat)))>0){
    bethat[which(is.na(bethat))] <- 0 #TODO
  }
  #print(bethat)
  #betax <- cbind(1, data_sup[, -outcome_colnum])[, 1:ncoef]%*%bethat
  #plot(density(betax))
  fi_hat <- c(cbind(1, data_sup[, -outcome_colnum])[, 1:ncoef]%*%bethat)*data_sup[, covariate_name]
  fj_hat <- (cbind(1, W_unlabel)%*%bethat)*W_unlabel[, covariate_name]
  ##rhat_ssl = mean(npreg(bws=bw,txdat=fi_hat,tydat=ri_hat, exdat = fj_hat)$mean,na.rm=T)
  rhat_ssl <- smooth_sslCPP(ri = ri_hat, fi = fi_hat, fnew = fj_hat, rsup = rhat_sup,
                            wgt = weights, bw = bw, cdf_trans = cdf_trans)
  
  mij_hat <- rhat_ssl[4:length(rhat_ssl)]
  bw <- rhat_ssl[3]
  rhat_ssl_bc <- rhat_ssl[2]
  rhat_ssl <- rhat_ssl[1]
  
  return(list("rhat" = c("Supervised" = rhat_sup,"NoSmooth" = mean(c(fi_hat,fj_hat)), "SemiSupervised" = rhat_ssl,
                         "SemiSupervisedBC" = rhat_ssl_bc),
              "var" = c("Supervised" = mean(pi^2*(ri_hat-rhat_sup)^2/nn), #sum(Vi^2*(ri_hat-rhat_sup)^2/sum(Vi)^2)
                        "NoSmooth" = mean(pi^2*(ri_hat-fi_hat)^2/nn), 
                        "SemiSupervised" = mean(pi^2*(ri_hat-mij_hat[1:nn])^2/nn),
                        "SemiSupervisedBC" = mean(pi^2*(ri_hat-mij_hat[1:nn])^2/nn)),
              "bw" = bw,
              "data_sup" = data_sup,
              "W_unlabel" = W_unlabel)
  )
}



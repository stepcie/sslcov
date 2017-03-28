#'Wrapper function for testing TODO
#'
#'
#'@param y vector of length \code{N} of the outcome
#'@param x vector of length \code{N} containing the covariate of interest
#'@param index_sup vector of length \code{n} indicating which observations whave supervised info
#'@param adjust_covariates optional matrix with \code{N} rows containing covariates to adjust on in
#'the semi-supervised estimation. Default is \code{NULL}
#'@param surrogate matrix with \code{N} rows containing the surrogate information
#'@param sampling_weights optional vector of length \code{n} containing the weights with which the
#'\code{n} supervised observation were sampled from the population. Default is random sampling (all the weigths equal to 1).
#'@param nperturb The number of perturbations to be run. Default is \code{500}.
#'
#'@param do_interact logical flag indicating whether interactins between \code{x} and
#'covariates should be taken into account when imputing \code{y}. Default is \code{TRUE}.
#'
#'@param condi logical flag indicating whether the covariance estimated should be condition on the
#'covariates indicated in . Default is \code{TRUE}.
#'@param do_ptb logical flag indicating whether to use perturbation to calculate the variance instead of using the asymptotic variance
#'
#'@importFrom stats na.omit sd quantile qnorm
#'
#'@examples
#'\dontrun{
#'#rm(list=ls())
#'
#'#Simulate data
#'nn_divide <- 10
#'NN <- 2000
#'nn <- NN/nn_divide
#'mySigma <- matrix(rep(0.6,16), 4, 4) + 0.4*diag(4)
#'mySigma[3,4] <- mySigma[3,4] + 0.2
#'mySigma[4,3] <- mySigma[4,3] + 0.2
#'beta <- 0 #0.6 #
#'beta_X <- c(0, 0, 0) # c(0.02, 0.3, -0.12) #
#'#set.seed(1234)
#'data_sim <- sim_data(ntot = NN, Sigma = 3*mySigma, b_G = beta, b_X = beta_X, cond_cov = TRUE)
#'cov_sim <- data_sim$cov_cond
#'data_sim <- data_sim$data
#'cov_sim
#'cov(data_sim[,"Y"], data_sim[,"G"])
#'es <- extremeSampling(data_sim, nn=nn, surrogate_name=c("S1", "S2", "S3"))
#'data_sampled <- rbind(data_sim[es$extreme_index,], data_sim[-es$extreme_index,])
#'
#'#True Covariance:
#'cov(data_sim[,"Y"], log(1+data_sim[,"G"]))
#'cov(data_sampled[1:nn,"Y"], log(1+data_sampled[1:nn,"G"]))
#'cov(data_sim[1:nn,"Y"], log(1+data_sim[1:nn,"G"]))
#'
#'res_ssl_randomsampling <- sslcov_test(y = data_sim[,"Y"], x = log(1 + data_sim[,"G"]),
#'                                      index_sup = 1:nn,
#'                                      surrogate = data_sim[,c("S1", "S2", "S3")],
#'                                      do_interact=FALSE, condi = FALSE, do_ptb=FALSE)
#'res_ssl <- sslcov_test(y = data_sampled[,"Y"], x = log(1 + data_sampled[,"G"]),
#'                       index_sup = 1:nn,
#'                       surrogate = data_sampled[,c("S1", "S2", "S3")],
#'                       sampling_weights = es$weights,
#'                       do_interact=FALSE, condi = FALSE, do_ptb=FALSE)
#'res_ssl_noWeights <- sslcov_test(y = data_sampled[,"Y"], x = log(1 + data_sampled[,"G"]),
#'                                 index_sup = 1:nn,
#'                                 surrogate = data_sampled[,c("S1", "S2", "S3")],
#'                                 do_interact=FALSE, condi = FALSE, do_ptb=FALSE)
#'
#'
#'# Conditional:
#'cov_sim
#'cov(data_sim[,"Y"], data_sim[,"G"])
#'cov(data_sim[,"Y"], log(1+data_sim[,"G"]))
#'cov(lm(data_sim[,"Y"]~data_sim[, c("Age", "Race", "Gender")])$residuals, data_sim[,"G"] -
#'    exp(MASS::glm.nb(data_sim[,"G"]~data_sim[, c("Age", "Race", "Gender")])$linear.predictors))
#'#library(profvis)
#'#profvis(
#'res_ssl_random_condi <- sslcov_test(y = data_sim[,"Y"], x = data_sim[,"G"], index_sup = 1:nn,
#'                          surrogate = data_sim[,c("S1", "S2", "S3")],
#'                          adjust_covariates = data_sim[, c("Age", "Race", "Gender"), 
#'                                                        drop=FALSE],
#'                          do_interact=FALSE, condi = TRUE, do_ptb=FALSE)
#'#)
#'res_ssl_condi <- sslcov_test(y = data_sampled[,"Y"], x = data_sampled[,"G"], index_sup = 1:nn,
#'                          surrogate = data_sampled[,c("S1", "S2", "S3")],
#'                          adjust_covariates = data_sampled[, c("Age", "Race", "Gender"),
#'                                                            drop=FALSE],
#'                          sampling_weights = es$weights,
#'                          do_interact=FALSE, condi = TRUE, do_ptb=FALSE)
#'#
#'res_ssl_noWeights_condi <- sslcov_test(y = data_sampled[,"Y"], x = data_sampled[,"G"],
#'                            index_sup = 1:nn, surrogate = data_sampled[,c("S1", "S2", "S3")],
#'                            adjust_covariates = data_sampled[, c("Age", "Race", "Gender"), 
#'                                                              drop=FALSE],
#'                            do_interact=FALSE, condi = TRUE, do_ptb=FALSE)
#'}
#'
#'@export
sslcov_test <- function(y, x, index_sup, surrogate, adjust_covariates=NULL,
                        sampling_weights = rep(1/length(y), length(index_sup)),
                        nperturb = 500, do_interact=TRUE, condi=FALSE, do_ptb=TRUE){
  
  if(condi & is.null(adjust_covariates)){
    stop("no covariates to condition on")
  }
  
  n <- length(index_sup)
  N <- length(y)
  
  stopifnot(nrow(x)==N)
  stopifnot(nrow(surrogate)==N)
  stopifnot(nrow(adjust_covariates)==N)
  stopifnot(length(sampling_weights)==n)
  
  ri0 <-(y-mean(y, na.rm = TRUE))*(x-mean(x))
  ri_nowgt <- (y[index_sup]-mean(y[index_sup], na.rm = TRUE))*(x[index_sup]-mean(x[index_sup]))
  
  
  data <- cbind(y, x, surrogate)
  surrogate_names <- paste0("S", 1:ncol(surrogate))
  if(!is.null(adjust_covariates)){
    data <- cbind(data, adjust_covariates)
    adjust_covariates_name <- paste0("AC", 1:ncol(adjust_covariates))
    colnames(data) <- c("y", "x", surrogate_names, adjust_covariates_name)
  }else{
    adjust_covariates_name <- NULL
    colnames(data) <- c("y", "x", surrogate_names)
  }
  data <- as.matrix(data)
  data_ord <- rbind(data[index_sup,], data[-index_sup,])
  weights_ord <- sampling_weights
  
  if(condi){
    rhat_out <- rhat_cond(data=data_ord, nn = n, outcome_name = "y", covariate_name = "x",
                          surrogate_name = surrogate_names, weights = weights_ord, cdf_trans = TRUE,
                          adjust_covariates_name = adjust_covariates_name, do_interact = do_interact)
  }else{
    rhat_out <- rhat(data=data_ord, nn = n, outcome_name = "y", covariate_name = "x",
                     surrogate_name = surrogate_names, weights = weights_ord,cdf_trans = TRUE,
                     adjust_covariates_name = adjust_covariates_name, do_interact = do_interact)
  }
  rhat_sup_w <-  rhat_out$rhat["Supervised"]
  rhat_ns_w <-  rhat_out$rhat["NoSmooth"]
  rhat_ssl_w <-  rhat_out$rhat["SemiSupervised"]
  rhat_ssl_bc_w <-  rhat_out$rhat["SemiSupervisedBC"]
  
  bw <- rhat_out$bw
  
  if(do_ptb){
    if(condi){
      res_ptb <-  na.omit(t(sapply(X = 1:nperturb, FUN = rhat_ptb_cond, data = data_ord, nn = n,
                                   outcome_name = "y", covariate_name = "x",
                                   surrogate_name = surrogate_names,
                                   bw = bw, cdf_trans = TRUE,
                                   weights = weights_ord, ptb_beta = TRUE,
                                   adjust_covariates_name = adjust_covariates_name,
                                   do_interact = do_interact)))
    }else{
      res_ptb <-  na.omit(t(sapply(X = 1:nperturb, FUN = rhat_ptb, data = data_ord, nn = n,
                                   outcome_name = "y", covariate_name = "x",
                                   surrogate_name = surrogate_names,
                                   bw = bw, cdf_trans = TRUE,
                                   weights = weights_ord, ptb_beta = TRUE,
                                   adjust_covariates_name = adjust_covariates_name,
                                   do_interact = do_interact)))
    }
    res_ptb_finite <- res_ptb[apply(abs(res_ptb), 1, max) < Inf,]
    
    rhat_sup_w_ptb <- res_ptb_finite[,"Supervised"]
    rhat_ns_w_ptb <- res_ptb_finite[,"NoSmooth"]
    rhat_ssl_w_ptb <- res_ptb_finite[,"SemiSupervised"]
    rhat_ssl_bc_w_ptb <- res_ptb_finite[,"SemiSupervisedBC"]
    
    
    #beta_bias <- rhat_out$beta_lm - colMeans(res_ptb_finite[, -c(1:3)])
    
    #sup
    sigma_sup_w <- sd(rhat_sup_w_ptb)
    # pval_sup_w <- 2*(1-pnorm(abs(rhat_sup_w/sigma_sup_w)))
    # pval_sup_w_dist = 2*(1-pnorm(abs(rhat_sup_w_ptb-rhat_sup_w)/sigma_sup_w))
    #sup bias correction
    se_sup_w <- sqrt(mean((rhat_sup_w - rhat_sup_w_ptb)^2))
    # bias_sup <- mean(rhat_sup_w - rhat_sup_w_ptb)
    # rhat_sup_w_unbias <- rhat_sup_w - bias_sup
    # pval_sup_w_unbias <- 2*(1-pnorm(abs(rhat_sup_w-bias_sup)/se_sup_w))
    pval_sup_w <- 2*(1-pnorm(abs(rhat_sup_w)/se_sup_w))
    
    #ssl
    sigma_ssl_w <-  sd(rhat_ssl_w_ptb)
    # pval_ssl_w <-  2*(1-pnorm(abs(rhat_ssl_w/sigma_ssl_w)))
    # pval_ssl_w_dist = 2*(1-pnorm(abs(rhat_ssl_w_ptb-rhat_ssl_w)/sigma_ssl_w))
    #ssl bias correction
    se_ssl_w <- sqrt(mean((rhat_ssl_w - rhat_ssl_w_ptb)^2))
    # bias_ssl <- mean(rhat_ssl_w - rhat_ssl_w_ptb)
    # rhat_ssl_w_unbias <- rhat_ssl_w - bias_ssl
    # pval_ssl_w_unbias <- 2*(1-pnorm(abs(rhat_ssl_w-bias_ssl)/se_ssl_w))
    pval_ssl_w <- 2*(1-pnorm(abs(rhat_ssl_w)/se_ssl_w))
    
    #ssl bc
    sigma_ssl_bc_w <-  sd(rhat_ssl_bc_w_ptb)
    se_ssl_bc_w <- sqrt(mean((rhat_ssl_bc_w - rhat_ssl_bc_w_ptb)^2))
    # bias_ssl_bc <- mean(rhat_ssl_bc_w - rhat_ssl_bc_w_ptb)
    # rhat_ssl_bc_w_unbias <- rhat_ssl_bc_w - bias_ssl
    # pval_ssl_bc_w_unbias <- 2*(1-pnorm(abs(rhat_ssl_bc_w-bias_ssl_bc)/se_ssl_bc_w))
    pval_ssl_bc_w <- 2*(1-pnorm(abs(rhat_ssl_bc_w)/se_ssl_w))
    pval_ssl_bc_w_bc <- 2*(1-pnorm(abs(rhat_ssl_bc_w)/se_ssl_bc_w))
    
    #no smoothing
    sigma_ns_w <- sd(rhat_ns_w_ptb)
    # pval_ns_w <- 2*(1-pnorm(abs(rhat_ns_w/sigma_ns_w)))
    # pval_ns_w_dist = 2*(1-pnorm(abs(rhat_ns_w_ptb-rhat_ns_w)/sigma_ns_w))
    #sup bias correction
    se_ns_w <- sqrt(mean((rhat_ns_w - rhat_ns_w_ptb)^2))
    # bias_ns <- mean(rhat_ns_w - rhat_ns_w_ptb)
    # rhat_ns_w_unbias <- rhat_ns_w - bias_ns
    # pval_ns_w_unbias <- 2*(1-pnorm(abs(rhat_ns_w-bias_ns)/se_ns_w))
    pval_ns_w <- 2*(1-pnorm(abs(rhat_ns_w)/se_ns_w))
    #
    # #combined
    # t_combined_w = -2*(log(pval_sup_w)+log(pval_ssl_w))
    # t.combined_w_dist = -2*(log(pval_sup_w_dist)+log(pval_ssl_w_dist))
    # pval_combined_w = mean(t_combined_w < t.combined_w_dist)
    
    
    ans <- list("r0"=mean(ri0, na.rm = TRUE),
                #"rhat_sup_unweighted"=mean(ri_nowgt, na.rm = TRUE),
                "rhat_sup_w" = rhat_sup_w, "sigma_sup_w" = sigma_sup_w, "pval_sup_w" = pval_sup_w,
                #"rhat_sup_w_unbias" = rhat_sup_w_unbias,
                "se_sup_w" = se_sup_w, #"pval_sup_w_unbias" = pval_sup_w_unbias,
                "rhat_ssl_w" = rhat_ssl_w, "sigma_ssl_w" = sigma_ssl_w, "pval_ssl_w" = pval_ssl_w,
                #"rhat_ssl_w_unbias" = rhat_ssl_w_unbias,
                "se_ssl_w"=se_ssl_w, #"pval_ssl_w_unbias" = pval_ssl_w_unbias,
                "rhat_ssl_w_avgPtb"=mean(rhat_ssl_w_ptb),
                "rhat_ssl_bc_w" = rhat_ssl_bc_w,"sigma_ssl_bc_w" = sigma_ssl_bc_w,"se_ssl_bc_w"=se_ssl_bc_w,
                "pval_ssl_bc_w" = pval_ssl_bc_w,"pval_ssl_bc_w_bc" = pval_ssl_bc_w_bc,
                "rhat_ns_w" = rhat_ns_w, "sigma_ns_w" = sigma_ns_w, "pval_ns_w" = pval_ns_w,
                #"rhat_ns_w_unbias" = rhat_ns_w_unbias,
                "se_ns_w"=se_ns_w, #"pval_ns_w_unbias" = pval_ns_w_unbias,
                #"rhat_ns_w_avgPtb"=mean(rhat_ns_w_ptb),
                #"pval_combined_w" = pval_combined_w,
                "q025_sup" = quantile(rhat_sup_w_ptb, 0.025), "q975_sup" = quantile(rhat_sup_w_ptb, 0.975),
                "q025_ssl" = quantile(rhat_ssl_w_ptb, 0.025), "q975_ssl" = quantile(rhat_ssl_w_ptb, 0.975),
                "q025_ssl_bc" = quantile(rhat_ssl_bc_w_ptb, 0.025), "q975_ssl_bc" = quantile(rhat_ssl_bc_w_ptb, 0.975),
                "q025_ns" = quantile(rhat_ns_w_ptb, 0.025), "q975_ns" = quantile(rhat_ns_w_ptb, 0.975)
                #"beta_bias"=beta_bias,"beta"=rhat_out$beta_lm,"beta_se"=apply(beta_ptb,2,sd),
                #"q025_beta"=apply(beta_ptb,2,function(x) quantile(x,0.025)),"q975_beta"=apply(beta_ptb,2,function(x) quantile(x,0.975))
    )
  } else {
    #sup
    sigma_sup_w <- sqrt(rhat_out$var["Supervised"])
    se_sup_w <- sqrt(rhat_out$var["Supervised"])
    pval_sup_w <- 2*(1-pnorm(abs(rhat_sup_w)/se_sup_w))
    
    #ssl
    sigma_ssl_w <-  sqrt(rhat_out$var["SemiSupervised"])
    se_ssl_w <- sqrt(rhat_out$var["SemiSupervised"])
    pval_ssl_w <- 2*(1-pnorm(abs(rhat_ssl_w)/se_ssl_w))
    
    #ssl bc
    sigma_ssl_bc_w <-  sqrt(rhat_out$var["SemiSupervisedBC"])
    se_ssl_bc_w <- sqrt(rhat_out$var["SemiSupervisedBC"])
    pval_ssl_bc_w <- 2*(1-pnorm(abs(rhat_ssl_bc_w)/se_ssl_w))
    pval_ssl_bc_w_bc <- 2*(1-pnorm(abs(rhat_ssl_bc_w)/se_ssl_bc_w))
    
    #no smoothing
    sigma_ns_w <- sqrt(rhat_out$var["NoSmooth"])
    se_ns_w <- sqrt(rhat_out$var["NoSmooth"])
    pval_ns_w <- 2*(1-pnorm(abs(rhat_ns_w)/se_ns_w))
    
    ans <- list("r0"=mean(ri0, na.rm = TRUE),
                #"rhat_sup_unweighted"=mean(ri_nowgt, na.rm = TRUE),
                "rhat_sup_w" = rhat_sup_w, "sigma_sup_w" = sigma_sup_w, "pval_sup_w" = pval_sup_w,
                #"rhat_sup_w_unbias" = rhat_sup_w_unbias,
                "se_sup_w" = se_sup_w, #"pval_sup_w_unbias" = pval_sup_w_unbias,
                "rhat_ssl_w" = rhat_ssl_w, "sigma_ssl_w" = sigma_ssl_w, "pval_ssl_w" = pval_ssl_w,
                #"rhat_ssl_w_unbias" = rhat_ssl_w_unbias,
                "se_ssl_w"=se_ssl_w, #"pval_ssl_w_unbias" = pval_ssl_w_unbias,
                "rhat_ssl_w_avgPtb"=NULL,
                "rhat_ssl_bc_w" = rhat_ssl_bc_w,"sigma_ssl_bc_w" = sigma_ssl_bc_w,"se_ssl_bc_w"=se_ssl_bc_w,
                "pval_ssl_bc_w" = pval_ssl_bc_w,"pval_ssl_bc_w_bc" = pval_ssl_bc_w_bc,
                "rhat_ns_w" = rhat_ns_w, "sigma_ns_w" = sigma_ns_w, "pval_ns_w" = pval_ns_w,
                #"rhat_ns_w_unbias" = rhat_ns_w_unbias,
                "se_ns_w"=se_ns_w, #"pval_ns_w_unbias" = pval_ns_w_unbias,
                #"rhat_ns_w_avgPtb"=mean(rhat_ns_w_ptb),
                #"pval_combined_w" = pval_combined_w,
                "q025_sup" = rhat_sup_w+qnorm(0.025)*se_sup_w, "q975_sup" = rhat_sup_w+qnorm(0.975)*se_sup_w,
                "q025_ssl" = rhat_ssl_w+qnorm(0.025)*se_ssl_w, "q975_ssl" = rhat_ssl_w+qnorm(0.975)*se_ssl_w,
                "q025_ssl_bc" = rhat_ssl_bc_w+qnorm(0.025)*se_ssl_bc_w, "q975_ssl_bc" = rhat_ssl_bc_w+qnorm(0.975)*se_ssl_bc_w,
                "q025_ns" = rhat_ns_w+qnorm(0.025)*se_ns_w, "q975_ns" = rhat_ns_w+qnorm(0.975)*se_ns_w
                #"beta_bias"=beta_bias,"beta"=rhat_out$beta_lm,"beta_se"=apply(beta_ptb,2,sd),
                #"q025_beta"=apply(beta_ptb,2,function(x) quantile(x,0.025)),"q975_beta"=apply(beta_ptb,2,function(x) quantile(x,0.975))
    )
    
  }
  class(ans) <- "ssl_test"
  print(ans)
  return(ans)
}



#'Print methods for ssl_test results
#'@method print ssl_test
#'
#'@param x a \code{ssl_test} object
#'@param ... further arguments passed to or from other methods
#'
#'@export
print.ssl_test <- function(x, ...){
  sup <- x[grep("sup", names(x))]
  ssl <- x[grep("ssl", names(x))]
  ssl_bc <- x[grep("ssl_bc", names(x))]
  ns <- x[grep("ns", names(x))]
  tab <- rbind(sup[c(grep("rhat_sup_w", names(sup)),grep("se_sup_w", names(sup)), grep("pval_sup_w", names(sup)))],
               ssl[c(which(names(ssl)=="rhat_ssl_w"),grep("se_ssl_w", names(ssl)), grep("pval_ssl_w", names(ssl)))],
               ssl_bc[c(grep("rhat_ssl_bc_w", names(ssl_bc)),grep("se_ssl_bc_w", names(ssl_bc)), grep("pval_ssl_bc_w_bc", names(ssl_bc)))],
               ns[c(grep("rhat_ns_w", names(ns)),grep("se_ns_w", names(ns)), grep("pval_ns_w", names(ns)))])
  colnames(tab) <- gsub("_sup", "", names(x)[grep("sup", names(x))])[c(1:3)]
  rownames(tab) <- c("Supervised", "Semi-supervised" , "Semi-supervised unbiased","Parametric (no smoothing)")
  print(tab)
}


#'Print method for ssl_test results
#'@method summary ssl_test
#'
#'@param object a \code{ssl_test} object
#'@param ... further arguments passed to or from other methods
#'
#'@export
summary.ssl_test <- function(object, ...){
  sup <- object[grep("sup", names(object))]
  ssl <- object[grep("ssl", names(object))]
  ssl_bc <- object[grep("ssl_bc", names(object))]
  ns <- object[grep("ns", names(object))]
  tab <- rbind(unlist(sup[c(grep("rhat_sup_w", names(sup)),grep("se_sup_w", names(sup)), grep("pval_sup_w", names(sup)))]),
               unlist(ssl[c(which(names(ssl)=="rhat_ssl_w"),grep("se_ssl_w", names(ssl)), grep("pval_ssl_w", names(ssl)))]),
               unlist(ssl_bc[c(grep("rhat_ssl_bc_w", names(ssl_bc)),grep("se_ssl_bc_w", names(ssl_bc)), grep("pval_ssl_bc_w_bc", names(ssl_bc)))]),
               unlist(ns[c(grep("rhat_ns_w", names(ns)),grep("se_ns_w", names(ns)), grep("pval_ns_w", names(ns)))]))
  colnames(tab) <- gsub("_sup", "", names(object)[grep("sup", names(object))])[c(1:3)]
  rownames(tab) <- c("Supervised", "Semi-supervised" , "Semi-supervised unbiased","Parametric (no smoothing)")
  tab.df <- as.data.frame(tab)
  return(tab.df)
}
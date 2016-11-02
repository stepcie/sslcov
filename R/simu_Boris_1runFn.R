#' Simulation function computing 1 run implemented by Boris
#'
#'@examples
#'
#'sims1 <- simu_Boris_1runFn(nperturb=2)
#'
#' @export
simu_Boris_1runFn <- function(NN = 200000, nn_divide = 1000, nperturb = 500,
                              mySigma = matrix(rep(0.3,16),4,4) + 0.7*diag(4),
                              underTheNull = TRUE,
                              verbose = FALSE, i = NULL, extreme = TRUE,
                              beta_2bePtb=TRUE, cond=FALSE){

  if(verbose){cat(i, "/100...", sep="")}

  if(underTheNull){
    beta <- 0
  }else{
    beta <- 0.6
  }

  if(!is.null(i)){
    set.seed(i*147852)
  }


  nn <- NN/nn_divide

  data_sim <- sim_data(ntot = NN, Sigma = mySigma, b_G = beta)

  ri0 <- (data_sim[,"Y"]-mean(data_sim[,"Y"]))*(data_sim[,"G"]-mean(data_sim[,"G"]))

  if(extreme){
    es <- extremeSampling(data_sim, nn, surrogate_name=c("S1", "S2", "S3"))
    data_sampled <- rbind(data_sim[es$extreme_index,], data_sim[-es$extreme_index,])
    rhat_out <- rhat(data=data_sampled, nn=nn, outcome_name="Y", covariate_name="G",
                     surrogate_name=c("S1", "S2", "S3"), weights = es$weights)
  }else{
    ranbdom_index <- sample(1:NN, size = nn, replace = FALSE)
    data_sampled <- rbind(data_sim[ranbdom_index, ], data_sim[-ranbdom_index, ])
    rhat_out <- rhat(data=data_sampled, nn=nn, outcome_name="Y", covariate_name="G",
                     surrogate_name=c("S1", "S2", "S3"))
  }

  rhat_sup_w <-  rhat_out$rhat["Supervised"]
  rhat_ssl_w <-  rhat_out$rhat["SemiSupervised.rhat"]
  bw <- rhat_out$bw

  if(extreme){
    res_ptb <-  na.omit(t(sapply(1:nperturb, rhat_ptb, data=data_sampled, nn=nn, bw=bw,
                                 ptb_nolab=FALSE, outcome_name="Y",
                                 covariate_name="G",
                                 surrogate_name=c("S1", "S2", "S3"),
                                 weights = es$weights,
                                 beta_2bePtb=beta_2bePtb)))
  }else{
    res_ptb <-  na.omit(t(sapply(1:nperturb, rhat_ptb, data=data_sampled, nn=nn, bw=bw,
                                 ptb_nolab=FALSE, outcome_name="Y",
                                 covariate_name="G",
                                 surrogate_name=c("S1", "S2", "S3"),
                                 beta_2bePtb=beta_2bePtb)))
  }
  res_ptb_finite <- res_ptb[apply(abs(res_ptb),1,max) < Inf,]


  rhat_sup_w_ptb <- res_ptb_finite[,"rhat_sup"]
  rhat_ssl_w_ptb <- res_ptb_finite[,"rhat_ssl.rhat"]

  sigma_sup_w <- sd(rhat_sup_w_ptb)
  pval_sup_w <- 2*(1-pnorm(abs(rhat_sup_w/sigma_sup_w)))
  pval_sup_w_dist = 2*(1-pnorm(abs((rhat_sup_w_ptb-rhat_sup_w)/sigma_sup_w)))

  se_sup_w <- sqrt(mean((rhat_sup_w - rhat_sup_w_ptb)^2))
  bias_sup <- mean(rhat_sup_w - rhat_sup_w_ptb)
  pval_sup_w_bias <- 2*(1-pnorm(abs(rhat_sup_w-bias_sup)/se_sup_w))

  sigma_ssl_w <-  sd(rhat_ssl_w_ptb)
  pval_ssl_w <-  2*(1-pnorm(abs(rhat_ssl_w/sigma_ssl_w)))
  pval_ssl_w_dist = 2*(1-pnorm(abs((rhat_ssl_w_ptb-rhat_ssl_w)/sigma_ssl_w)))

  se_ssl_w <- sqrt(mean((rhat_ssl_w - rhat_ssl_w_ptb)^2))
  bias_ssl <- mean(rhat_ssl_w - rhat_ssl_w_ptb)
  pval_ssl_w_bias <- 2*(1-pnorm(abs(rhat_ssl_w-bias_ssl)/se_ssl_w))

  t_combined_w = -2*(log(pval_sup_w)+log(pval_ssl_w))
  t.combined_w_dist = -2*(log(pval_sup_w_dist)+log(pval_ssl_w_dist))
  pval_combined_w = mean(t_combined_w < t.combined_w_dist)

  res <- data.frame("r0"=mean(ri0),
                    "rhat_sup_w" = rhat_sup_w, "bias_sup" = bias_sup,
                    "sigma_sup_w" = sigma_sup_w, "se_sup_w" = se_sup_w,
                    "pval_sup_w" = pval_sup_w, "pval_sup_w_bias" = pval_sup_w_bias,
                    "rhat_ssl_w" = rhat_ssl_w, "rhat_ssl_w_avgPtb"=mean(rhat_ssl_w_ptb),
                    "bias_ssl" = bias_ssl, "sigma_ssl_w" = sigma_ssl_w, "se_ssl_w" = se_ssl_w,
                    "pval_ssl_w" = pval_ssl_w, "pval_ssl_w_bias" = pval_ssl_w_bias,
                    "pval_combined_w" = pval_combined_w,
                    "q025_sup" = quantile(rhat_sup_w_ptb, 0.025), "q975_sup" = quantile(rhat_sup_w_ptb, 0.975),
                    "q025_ssl" = quantile(rhat_ssl_w_ptb, 0.025), "q975_ssl" = quantile(rhat_sup_w_ptb, 0.975)
  )

  if(verbose){cat(" DONE !\n")}

  return(res)

}

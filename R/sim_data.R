#' Simulating data
#'
#' @param ntot the number of observation to simulate
#'
#' @param missing logical flag whether missing indicator should 
#' be used. Default is \code{FALSE}.
#'
#' @param b_G strength parameter for linear association between the covariate "G"
#' and the outcome of interest "Y".
#'
#' @return simulated data matrix of size \code{ntot x 4}
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rnbinom rbinom pnorm rnorm cov
#'
#' @keywords internal
#' @export


sim_data <-  function(ntot, missing=FALSE,incorrect=FALSE,
                      b_G=0.6, 
                      b_X=c(0.02, 0.3, -0.12), 
                      b_X_Gsize = c(0.04, -2, -3), b_X_Gprob=c(-0.01, 0.5, 1),
                      Sigma=diag(4), cond_cov=FALSE){
  
  Xi <- cbind("Age"=stats::rnorm(ntot, m=50, sd=7),
              "Race"=stats::rbinom(ntot, size=1, prob=0.7),
              "Gender"= stats::rbinom(ntot, size=1, prob=0.5))
  
  # #Pi <- rbinom(ntot, size=1, prob=0.25)
  # #Gi <- rnbinom(ntot, size=1, mu = 0.3 + 15*Pi + 0*5*Xi*(b_G==0))
  # Gi <- rnbinom(ntot, size=0.1, mu = exp(1.4 + 0.3*Xi*(b_G!=0)))
  # crp <- MASS::mvrnorm(ntot, c(0,0,0,0), Sigma) + b_X*b_G*Xi + b_G*Gi + rnorm(ntot,0,0.3)
  # #mtemp <- MASS::glm.nb(Gi~Xi)
  # #ist(Gi, n=100)
  
  
  # Gi <- rnbinom(ntot, size=1, 0.3)
  # Gi <- log(Gi+1)
  # #Gi <- rbinom(ntot, size=3, 0.1)
  # #crp <- 0.5*MASS::mvrnorm(ntot, c(-1,-1,-1,-1), Sigma) +0.5*MASS::mvrnorm(ntot, c(1,1,1,1), Sigma)  + b_G*Gi + rnorm(ntot,0,0.3)
  # crp <- MASS::mvrnorm(ntot, c(0,0,0,0), Sigma) + b_G*Gi + rnorm(ntot,0,0.3)
  #
  # Yi <- crp[, 1]
  # Si <- crp[, 2:4]
  # if(missing){
  #   pmiss<-max(0,0.5-0.01*(Gi+1))
  #   Imiss <-  rbinom(ntot, 1, pmiss)
  #   Si[Imiss==1, ] <-  NA
  # }
  # if(incorrect){
  #   #Si <- exp(Si + Gi*rnorm(ntot,0,0.03))
  #   Si <- exp(Si)
  #
  #   # delta <- rbinom(ntot,1,0.5)
  #   # crp <- delta*MASS::mvrnorm(ntot, c(-1,-1,-1,-1), Sigma) +(1-delta)*MASS::mvrnorm(ntot, c(1,1,1,1), Sigma)  + b_G*Gi + rnorm(ntot,0,0.3)
  #   # Yi <- crp[, 1]
  #   # Si <- crp[, 2:4]
  #
  # }
  # data <- cbind(Yi, Gi, Si)
  # colnames(data) <- c("Y", "G", paste("S",1:ncol(Si), sep=""))
  # return(data)
  
  size_G=1
  prob_G=0.2
  mG=size_G*(1-prob_G)/prob_G
  
  #E(G|X)=E(G)
  cov_cond_est <- NULL
  if(cond_cov){
    Gi_raw <- rnbinom(ntot, size=pmax(floor(size_G + Xi%*%b_X_Gsize), 1), 
                      prob=expit(logit(prob_G) + Xi%*%b_X_Gprob))
    
    nMC <- 1000
    cov_X <- numeric(nMC)
    covlog_X <- numeric(nMC)
    #stats::cov(Yi_MC, Gi_raw)
    
    for (j in 1:nMC){
      Gi_raw_MC <- floor(stats::rnbinom(ntot, size=pmax(floor(size_G + Xi%*%b_X_Gsize), 1),
                                 prob=expit(logit(prob_G) + Xi%*%b_X_Gprob)))
      Gi_MC <- log(1 + Gi_raw_MC)
      p.S_MC = nrow(Sigma)-1
      Yi_MC <- b_G*Gi_MC + as.vector(Xi%*%b_X) + MASS::mvrnorm(ntot, mu=rep(0, p.S_MC+1), Sigma)
      Yi_MC = Yi_MC[,1]
      cov_X[j] <- stats::cov(Yi_MC-mean(Yi_MC), Gi_raw_MC-mean(Gi_raw_MC))
      covlog_X[j] <- stats::cov(Yi_MC-mean(Yi_MC), Gi_MC-mean(Gi_MC))
    }
    cov_cond_est <- mean(cov_X)
    cov_cond_est_log <- mean(covlog_X)
  }else{
    Gi_raw <- stats::rnbinom(ntot, size=size_G, prob=prob_G)
  }
  
  Gi <- log(1 + Gi_raw)
  p.S = nrow(Sigma)-1
  Yi <- b_G*Gi + as.vector(Xi%*%b_X) + MASS::mvrnorm(ntot, mu=rep(0, p.S+1), Sigma)
  Si = Yi[,-1]
  Yi = Yi[,1]
  
  
  
  colnames(Si) = paste("S",1:ncol(Si),sep="")
  if(incorrect){
    Yi =  Yi + b_G*(Gi^2 - Gi)
    Si = Si - b_G*Gi^2
  }
  if(missing){
    Imiss <-  as.matrix(stats::rbinom(ntot, 1, stats::pnorm(1-Si[,1]-Si[,2])))
    colnames(Imiss) = paste("I_miss",1:ncol(Imiss),sep="")
    Si[Imiss==1, ] <-  landpred::VTM(apply(Si[Imiss==0,],2,mean),sum(Imiss))
  }else{
    Imiss = NULL
  }
  colnames(Si) = paste("S",1:ncol(Si),sep="")
  
  data <- cbind(Yi, Gi_raw, Xi, Si, Imiss)
  colnames(data)[1:2] <- c("Y", "G")
  
  return(list("data"=data, "cov_cond_est"=c(cov_cond_est, cov_cond_est_log)))
  
}

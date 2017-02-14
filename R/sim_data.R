#' Simulating data
#'
#' @param ntot the number of observation to simulate
#'
#' @param missing logical flag wether
#' . Default is FALSE.
#'
#' @param b_G strength parameter for linear association between the covariate "G"
#' and the outcome of interest "Y".
#'
#' @return simulated data matrix of size \code{ntot x 4}
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rnbinom rbinom pnorm
#'
#' @keywords internal
#' @export


sim_data <-  function(ntot, missing=FALSE,incorrect=FALSE,
                      b_G=0.6,
                      Sigma=diag(4)){
  # Xi <- rnorm(ntot)
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
  Gi <- log(1+rnbinom(ntot, size=2, 0.1))
  p.S = nrow(Sigma)-1
  Yi <- b_G*Gi + mvrnorm(ntot,mu=rep(0,p.S+1),Sigma)
  Si = Yi[,-1]
  Yi = Yi[,1]
  colnames(Si) = paste("S",1:ncol(Si),sep="")
  if(incorrect){
    Yi = Gi^2*b_G + Yi - b_G*Gi
    Si= Si-b_G*Gi^2
  }
  if(missing){
    Imiss <-  as.matrix(rbinom(ntot, 1, pnorm(1-Si[,1]-Si[,2])))
    colnames(Imiss) = paste("I_miss",1:ncol(Imiss),sep="")
    Si[Imiss==1, ] <-  landpred::VTM(apply(Si[Imiss==0,],2,mean),sum(Imiss))
  }else{
    Imiss = NULL
  }
  colnames(Si) = paste("S",1:ncol(Si),sep="")

  data <- cbind(Yi, Gi, Si, Imiss)
  colnames(data)[1:2] <- c("Y", "G")
  return(data)

}

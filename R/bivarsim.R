#' @name bivar.sim
#' @title Simulating an Interacted Vector Autoregression
#' @description This function is used to produce simulated realizations which follow a Vector Autorgression (GVAR). It will also automatically simulate coefficients. All parameters can also be set by the user.
#' @usage bvar.sim(len, M, plag=1, cons=FALSE, trend=FALSE, SV=FALSE)
#' @details For testing purposes, this function enables to simulate time series processes which can be described by a Global Vector Autoregression. Since stability conditions are not checked, it is only implemented for \code{M=3}.
#' @param len length of the simulated time series.
#' @param M number of endogenous variables.
#' @param plag number of lags.
#' @param cons logical indicating whether to include an intercept. Default set to \code{FALSE}.
#' @param trend logical indicating whether to include an intercept. Default set to \code{FALSE}.
#' @param SV logical indicating whether the process should be simulated with or without stochastic volatility. Default set to \code{FALSE}.
#' @return Returns a list with the following elements
#' @author Maximilian Boeck
#' @examples
#' library(BTSM)
#' sim  <- bivar.sim(len=200, M=3, plag=1, cons=TRUE, trend=FALSE, SV=FALSE)
#' Data <- sim$obs$Yraw
#' @importFrom stats rnorm
#' @importFrom stochvol svsim
#' @export
bivar.sim <- function(len, M=3, plag=1, cons=TRUE, SV=FALSE, type="dummy", share=0.1, sd.1=0.1){
  idi.para <- t(matrix(rep(c(-5,0.9,0.01),M),3,M))
  # -----------------------------------------------------------------------------------------
  burn <- 50
  # simluate shocks
  vol.true <- matrix(NA, len+plag+burn, M)
  if(SV){
    for(mm in 1:M) {
      temp <- svsim(len=len+plag+burn, mu=idi.para[mm,1], phi=idi.para[mm,2], sigma=idi.para[mm,3])
      vol.true[,mm] <- temp$vol
    }
  } else {
    for(mm in 1:M) {
      vol.true[,mm] <- exp(-5)
    }
  }
  Yraw <- matrix(0,len+plag+burn,M)
  colnames(Yraw) <- paste0("var.",seq(1,M))
  # interaction
  if(type=="dummy"){
    Draw <- matrix(0,len+plag+burn,1)
    Draw[sample(1:len, size=floor(len*share), replace=FALSE),] <- 1
  }else if(type=="continuous"){
    Draw <- matrix(rnorm(len+plag+burn,0,1),len+plag+burn,1)
  }
  #Draw[1:floor(len/3),] <- 1

  # state 0/1
  crit_eig0 <- crit_eig1 <- 1.1; icounter<-0
  while(crit_eig0>0.99 & crit_eig1>0.99 & icounter<1e7){
    icounter<-icounter+1
    # state 0
    A0.true   <- diag(M)*runif(1,0.7,0.9)
    A0.true[upper.tri(A0.true)] <- runif((M-1)*M/2,-0.2,0.2)
    A0.true[lower.tri(A0.true)] <- runif((M-1)*M/2,-0.2,0.2)
    # state 1
    A1.true   <- matrix(rnorm(M^2,0,sd.1),M,M)
    crit_eig0 <- max(abs(Re(eigen(A0.true)$values)))
    crit_eig1 <- max(abs(Re(eigen(A0.true+A1.true)$values)))
    #cat(paste0("\n Round: ",icounter," - State 0: ", crit_eig0," - State 1: ",crit_eig1))
  }
  if(crit_eig0>0.99) stop("Couldn't find suitable coefficient matrix, maybe try a smaller-dimensional system..")
  if(crit_eig1>0.99) stop("Couldn't find suitable coefficient matrix, maybe try a smaller-dimensional system..")

  # J
  Jstar0.true <- diag(M); Jstar0.true[lower.tri(Jstar0.true)] <- runif((M-1)*M/2,-1,1)# rnorm((M-1)*M/2,0,sd.1)
  Jstar1.true <- diag(M); #Jstar1.true[lower.tri(Jstar1.true)] <- rnorm((M-1)*M/2,0,sd.1)

  a0.true <- rep(0,M)
  a1.true <- rep(0,M)
  if(cons) a0.true <- runif(M,-0.2,0.5)
  if(cons) a1.true <- runif(M,-0.5,0.2)
  # -----------------------------------------------------------------------------------------
  for(pp in 1:plag) Yraw[pp,] <- a0.true + a1.true*Draw[pp,] + rnorm(M,0,vol.true[pp,])
  for(tt in (plag+1):(len+plag+burn)){
    xlag <- NULL; Abig <- NULL
    for(pp in 1:plag) {
      xlag <- cbind(xlag,Yraw[tt-pp,,drop=FALSE],Draw[tt,,drop=FALSE]%*%Yraw[tt-pp,,drop=FALSE])
      Abig <- rbind(Abig,A0.true,A1.true)
      rownames(Abig) <- colnames(xlag) <- paste0(rep(colnames(Yraw),2),rep(c("",".interaction"),each=M))
      colnames(Abig) <- colnames(Yraw)
    }
    if(cons){
      xlag <- cbind(xlag,1,1)
      Abig <- rbind(Abig,a0.true,a1.true)
      colnames(xlag)[(M*plag*2+1):ncol(xlag)] <- rownames(Abig)[(M*plag*2+1):nrow(Abig)] <- c("cons","cons.interaction")
    }
    for(mm in 1:M){
      Ause <- Abig[,mm,drop=FALSE]
      xuse <- xlag
      if(mm>1){
        xuse <- cbind(Yraw[tt,1:(mm-1),drop=FALSE],Draw[tt,]%*%Yraw[tt,1:(mm-1),drop=FALSE],xuse)
        Ause <- rbind(t(Jstar0.true[mm,1:(mm-1),drop=FALSE]),t(Jstar1.true[mm,1:(mm-1),drop=FALSE]),Ause)
      }
      Yraw[tt,mm] <- xuse%*%Ause + rnorm(1,0,vol.true[tt,mm])
    }
  }
  Yraw <- Yraw[(plag+1+burn):nrow(Yraw),,drop=FALSE]
  Draw <- Draw[(plag+1+burn):nrow(Draw),,drop=FALSE]
  vol.true <- vol.true[(plag+1+burn):nrow(vol.true),,drop=FALSE]
  # change sign of Jstar
  Jstar0.true <- -Jstar0.true; diag(Jstar0.true) <- 1
  Jstar1.true <- -Jstar1.true; diag(Jstar1.true) <- 1
  #-----------------------------------------------------------------------------------------
  nhor    <- 60
  impresp <- array(0, c(M, M, nhor, 2))
  # state 0
  S0.true <- solve(Jstar0.true)%*%diag(apply(vol.true,2,median))%*%t(solve(Jstar0.true))
  S1.true <- solve(Jstar0.true+Jstar1.true)%*%diag(apply(vol.true,2,median))%*%t(solve(Jstar0.true+Jstar1.true))
  temp0   <- .gen_compMat(A0.true, M, plag)
  temp1   <- .gen_compMat(A0.true+A1.true, M, plag)
  Jm0     <- temp0$Jm
  Cm0     <- temp0$Cm
  Jm1     <- temp1$Jm
  Cm1     <- temp1$Cm

  # identification
  shock0   <- t(chol(apply(S0.true,c(1,2),median)))
  diagonal <- diag(diag(shock0))
  shock0   <- solve(diagonal)%*%shock0 # unit initial shock

  shock1   <- t(chol(apply(S1.true,c(1,2),median)))
  diagonal <- diag(diag(shock1))
  shock1   <- solve(diagonal)%*%shock1 # unit initial shock

  impresp[,,1,1] <- shock0
  impresp[,,1,2] <- shock1
  compMati0 <- Cm0
  compMati1 <- Cm1
  for(j in 2:nhor) {
    temp0 <- t(Jm0) %*% compMati0 %*% Jm0 %*% shock0
    temp1 <- t(Jm1) %*% compMati1 %*% Jm1 %*% shock1
    compMati0 <- compMati0 %*% Cm0
    compMati1 <- compMati1 %*% Cm1
    impresp[,,j,1] <- temp0
    impresp[,,j,2] <- temp1
  }

  true.list <- list(A0.true=A0.true, A1.true=A1.true, a0.true=a0.true, a1.true=a1.true, Jstar0.true=Jstar0.true, Jstar1.true=Jstar1.true, vol.true=vol.true)
  return(list(Yraw=Yraw,Draw=Draw,true.list=true.list,impresp.true=impresp))
}

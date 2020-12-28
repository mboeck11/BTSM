#' @name tvpbvar.sim
#' @title Simulating a Time-varying Parameter Bayesian Vector Autoregression
#' @description This function is used to produce simulated realizations which follow a Vector Autorgression (GVAR). It will also automatically simulate coefficients. All parameters can also be set by the user.
#' @usage tvpbvar.sim(len, M, plag=1, cons=FALSE, trend=FALSE, SV=FALSE)
#' @details For testing purposes, this function enables to simulate time series processes which can be described by a Global Vector Autoregression. Since stability conditions are not checked, it is only implemented for \code{M=3}.
#' @param len length of the simulated time series.
#' @param M number of endogenous variables.
#' @param plag number of lags.
#' @param cons logical indicating whether to include an intercept. Default set to \code{FALSE}.
#' @param trend logical indicating whether to include an intercept. Default set to \code{FALSE}.
#' @param SV logical indicating whether the process should be simulated with or without stochastic volatility. Default set to \code{FALSE}.
#' @param sparse.coef sparsification of coefficients, has to be provided as percentage between zero and one.
#' @param sparse.tvp sparsification of time-variation, has to be provided as percentage between zero and one.
#' @param tvp.var variance of state process governing the time-variation in the coefficients
#' @return Returns a list with the following elements
#' @author Maximilian Boeck
#' @examples
#' library(BTSM)
#' sim <- tvpbvar.sim(len=200, M=3, plag=1, cons=TRUE, trend=FALSE, SV=FALSE)
#' Data = sim$obs$xglobal
#' W    = sim$obs$W
#' @importFrom stats rnorm
#' @importFrom stochvol svsim
#' @export
tvpbvar.sim <- function(len, M=3, plag=1, cons=FALSE, trend=FALSE, SV=FALSE,
                        sparse.coef=0, sparse.tvp=0, tvp.var=0.002, snr=0.2){
  if(sparse.coef>1 || sparse.coef<0){
    stop("Please provide argument 'sparse.coef' as percentage between zero and one.")
  }
  if(sparse.tvp>1 || sparse.tvp<0){
    stop("Please provide argument 'sparse.tvp' as percentage between zero and one.")
  }
  vol.mean <- log(tvp.var/snr)
  idi.para <- t(matrix(rep(c(vol.mean,0.9,0.01),M),3,M))
  # -----------------------------------------------------------------------------------------
  # simluate shocks
  vol.true <- matrix(NA, len, M)
  if(SV){
    for(mm in 1:M) {
      temp <- svsim(len=len, mu=idi.para[mm,1], phi=idi.para[mm,2], sigma=idi.para[mm,3])
      vol.true[,mm] <- temp$vol
    }
  } else {
    for(mm in 1:M) {
      vol.true[,mm] <- exp(vol.mean)
    }
  }
  Yraw <- matrix(0,len+plag,M)

  # state 0
  a0.true  <- rep(0,M)
  a1.true  <- rep(0,M)
  indic.coef <- sample(1:(M^2*plag), size=floor(M^2*plag*sparse.coef))
  crit_eig <- 1.05; icounter<-0
  while(crit_eig>1.0 & icounter<1e4){
    icounter<-icounter+1
    A.true   <- diag(M)*runif(1,0.7,0.9)
    A.true[upper.tri(A.true)] <- runif((M-1)*M/2,-0.2,0.2)
    A.true[lower.tri(A.true)] <- runif((M-1)*M/2,-0.2,0.2)
    avec <- as.vector(A.true)
    if(length(indic.coef)>0){
      for(kk in 1:length(indic.coef)){
        avec[indic.coef[kk]] <- 0
      }
    }
    A.true <- matrix(avec,M*plag,M)
    crit_eig <- max(abs(Re(eigen(A.true)$values)))
  }
  if(crit_eig>1.00) stop("Couldn't find suitable coefficient matrix, maybe try a smaller-dimensional system..")
  At.true <- array(NA, c(len,M,M))
  Qt.true <- diag(M^2*plag)*tvp.var
  indic.tvp<-sample(1:(M^2*plag), size=floor(M^2*plag*sparse.tvp))
  if(length(indic.tvp)>0){
    for(kk in 1:length(indic.tvp)){
      Qt.true[indic.tvp[kk],indic.tvp[kk]] <- 0
    }
  }
  At.true[1,,] <- A.true
  tt<-2
  while(tt <= len){
    temp <- At.true[tt-1,,] + matrix(rnorm(M^2*plag,rep(0,M^2*plag),diag(Qt.true)),M*plag,M)
    if(max(abs(Re(eigen(.gen_compMat(temp,M,plag)$Cm)$values)))<crit_eig){
      At.true[tt,,] <- temp
      tt<-tt+1
    }
  }
  # par(mfrow=c(M,M),mar=c(2,2,1,1))
  # for(mm in 1:M) for(mmm in 1:M) plot.ts(At.true[,mm,mmm])

  L.true   <- diag(M)
  L.true[lower.tri(L.true)] <- runif((M-1)*M/2,-0.2,0.2)
  S.true   <- array(NA,c(M,M,len))

  if(cons) a0.true <- runif(M,-1,1)
  if(trend) a1.true <- runif(M,-1,1)
  # -----------------------------------------------------------------------------------------
  S.true[,,1] <- L.true%*%diag(vol.true[1,])%*%t(L.true)
  for(pp in 1:plag) Yraw[pp,] <- a0.true
  for (tt in 1:len){
    S.true[,,tt] <- L.true%*%diag(vol.true[tt,])%*%t(L.true)
    xlag <- NULL; Abig <- At.true[tt,,]
    for(pp in 1:plag) xlag <- cbind(xlag,Yraw[tt+plag-1,,drop=FALSE])
    if(cons){
      xlag <- cbind(xlag,1)
      Abig <- rbind(Abig,a0.true)
    }
    if(trend){
      xlag <- cbind(xlag,tt)
      Abig <- rbind(Abig,a1.true)
    }
    Yraw[tt+plag,] <- xlag%*%Abig + t(t(chol(S.true[,,tt]))%*%rnorm(M))
  }
  Yraw <- Yraw[(plag+1):(len+plag),,drop=FALSE]
  #-----------------------------------------------------------------------------------------
  temp <- .gen_compMat(A.true, M, plag)
  Jm   <- temp$Jm
  Cm   <- temp$Cm

  # identification
  shock <- t(chol(apply(S.true,c(1,2),median)))
  diagonal <- diag(diag(shock))
  shock <- solve(diagonal)%*%shock # unit initial shock

  nhor    <- 60
  impresp <- array(0, c(M, M, nhor))
  impresp[,,1]  <- t(shock)
  compMati <- Cm
  for(j in 2:nhor) {
    temp <- t(Jm) %*% compMati %*% Jm %*% shock
    compMati <- compMati %*% Cm
    impresp[,,j] <- temp
  }

  true.list <- list(A.true=A.true, At.true=At.true, a0.true=a0.true, a1.true=a1.true, L.true=L.true, vol.true=vol.true, S.true=S.true)

  return(list(Yraw=Yraw,true.list=true.list,impresp.true=impresp))
}

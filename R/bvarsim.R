#' @name bvar.sim
#' @title Simulating a Vector Autoregression
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
#' sim <- bvar.sim(len=200, M=3, plag=1, cons=TRUE, trend=FALSE, SV=FALSE)
#' Data = sim$obs$xglobal
#' W    = sim$obs$W
#' @importFrom stats rnorm
#' @importFrom stochvol svsim
#' @export
bvar.sim <- function(len, M=3, plag=1, cons=FALSE, trend=FALSE, SV=FALSE){
  idi.para <- t(matrix(rep(c(-5,0.9,0.01),M),3,M))
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
      vol.true[,mm] <- exp(-5)
    }
  }
  Yraw <- matrix(0,len,M)

  # state 0
  a0.true  <- rep(0,M)
  a1.true  <- rep(0,M)
  crit_eig <- 1.1; icounter<-0
  while(crit_eig>1.0 & icounter<1e4){
    icounter<-icounter+1
    A.true   <- diag(M)*runif(1,0.7,0.9)
    A.true[upper.tri(A.true)] <- runif((M-1)*M/2,-0.2,0.2)
    A.true[lower.tri(A.true)] <- runif((M-1)*M/2,-0.2,0.2)
    crit_eig <- max(abs(Re(eigen(A.true)$values)))
  }
  if(crit_eig>1.00) stop("Couldn't find suitable coefficient matrix, maybe try a smaller-dimensional system..")
  #A.true   <- matrix(c(0.7,-0.03,0.08,-0.11,0.4,0,0.1,0.03,0.7),M,M)
  L.true   <- diag(M)
  L.true[lower.tri(L.true)] <- runif((M-1)*M/2,-1,1)
  S.true   <- array(NA,c(M,M,len))

  if(cons) a0.true <- runif(M,-1,1)
  if(trend) a1.true <- runif(M,-1,1)
  # -----------------------------------------------------------------------------------------
  S.true[,,1] <- L.true%*%diag(vol.true[1,])%*%t(L.true)
  for(pp in 1:plag) Yraw[pp,] <- a0.true + pp*a1.true + t(chol(S.true[,,1]))%*%rnorm(M)
  for (tt in (plag+1):len){
    S.true[,,tt] <- L.true%*%diag(vol.true[tt,])%*%t(L.true)
    xlag <- NULL; Abig <- A.true
    for(pp in 1:plag) xlag <- cbind(xlag,Yraw[tt-1,,drop=FALSE])
    if(cons){
      xlag <- cbind(xlag,1)
      Abig <- rbind(Abig,a0.true)
    }
    if(trend){
      xlag <- cbind(xlag,tt)
      Abig <- rbind(Abig,a1.true)
    }
    Yraw[tt,] <- xlag%*%Abig + t(t(chol(S.true[,,tt]))%*%rnorm(M))
  }
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

  true.list <- list(A.true=A.true, a0.true=a0.true, a1.true=a1.true, L.true=L.true, vol.true=vol.true, S.true=S.true)

  return(list(Yraw=Yraw,true.list=true.list,impresp.true=impresp))
}

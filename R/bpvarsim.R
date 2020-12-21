#' @name bpvar.sim
#' @title Simulating a Panel Vector Autoregression
#' @description This function is used to produce simulated realizations which follow a Vector Autorgression (GVAR). It will also automatically simulate coefficients. All parameters can also be set by the user.
#' @usage bpvar.sim(len, M, N, plag=1, cons=FALSE, trend=FALSE, SV=FALSE)
#' @details For testing purposes, this function enables to simulate time series processes which can be described by a Global Vector Autoregression. Since stability conditions are not checked, it is only implemented for \code{M=3}.
#' @param len length of the simulated time series.
#' @param M number of endogenous variables.
#' @param N number of countries.
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
bpvar.sim <- function(len, M=3, N=5, plag=1, cons=FALSE, trend=FALSE, SV=FALSE, ex=FALSE){
  idi.para <- t(matrix(c(-10,0.95,0.01),M))
  # -----------------------------------------------------------------------------------------
  # simluate shocks
  shock <- matrix(NA, len+plag, M*N)
  vol_true <- matrix(NA, len+plag, M*N)
  if(SV) {
    for(mm in 1:(M*N)) {
      temp <- svsim(len=len+plag, mu  = idi.para[1,1],
                    phi = idi.para[1,2],
                    sigma = idi.para[1,3])
      shock[,mm] <- temp$y
      vol_true[,mm] <- temp$vol
    }
  } else {
    for(mm in 1:(M*N)) {
      vol_true[,mm] <- rexp(1, rate = 2)
      shock[,mm]    <- rnorm(len+plag, 0, vol_true[,mm])
    }
  }

  par(mfrow=c(1,1))
  ts.plot(vol_true, col = palette()[1:(M*N)])
  ts.plot(shock, col = palette()[1:(M*N)])

  Yraw <- matrix(0,len+plag,M*N)
  if(ex) Wraw <- matrix(rnorm(2*(len+plag)),len+plag,2) else Wraw <- NULL

  ident.country <- t(matrix(seq(1,M*N),M,N))
  cN <- paste0(letters[1:N],letters[1:N],sep="")

  A1.mean <- matrix(c(0.5,-0.2,0.3,
                      0.5,0.8,-0.3,
                      -0.2,0.1,0.9),M,M)
  if(ex) B1.mean <- matrix(c(0.05,0.03,-0.1,0.2,-0.02,0.2),2,M) else B1.mean <- NULL
  max(abs(Re(eigen(A1.mean)$values)))

  V.A <- diag(M) * 0.001
  # -----------------------------------------------------------------------------------------
  # simulate VAR coefficients
  A <- array(NA, dim=c(M, M, N))
  if(ex) B <- array(NA, dim=c(2, M, N)) else B <- NULL
  SIGMA <- array(NA, dim=c(M, M, N))
  for(cc in 1:N) {
    sl.country <- ident.country[cc,]

    for(mm in 1:M) {
      A[,mm,cc] <- rnorm(M*plag, A1.mean[,mm], sd = sqrt(diag(V.A)))
      if(ex) B[,mm,cc] <- rnorm(2,   B1.mean[,mm], sd = sqrt(diag(V.A)))
    }

    SIGMA[,,cc] <- var(shock[,sl.country])
  }

  # -----------------------------------------------------------------------------------------
  # simulate data
  for(cc in 1:N){
    sl.country <- ident.country[cc,]
    if(ex) Yraw[1:plag,] <- Wraw[1:plag,]%*%B[,,cc] + shock[1:plag,sl.country] else Yraw[1:plag,] <- shock[1:plag,sl.country]
    for (tt in (plag+1):(len+plag)){
      if(ex) Yraw[tt,sl.country] <- Yraw[tt-1,sl.country] %*% A[,,cc] + Wraw[tt,] %*% B[,,cc] + shock[tt,sl.country] else Yraw[tt,sl.country] <- Yraw[tt-1,sl.country] %*% A[,,cc] + shock[tt,sl.country]
    }
  }

  Yraw <- Yraw[(plag+1):nrow(Yraw),,drop=FALSE]
  if(ex) Wraw <- Wraw[(plag+1):nrow(Wraw),,drop=FALSE]

  # par(mfrow=c(1,1))
  # ts.plot(Yraw, col = palette()[1:M])
  #
  # colnames(Yraw) <- paste0(rep(cN,each=M),".",c("y","inf","stir"))
  # if(ex) colnames(Wraw) <- c("oilprice","comprice")

  # MEAN BETA
  A_draw_mean <- apply(A, c(1,2), mean)
  # MEAN SIGMA
  SIGMA_mean <- apply(SIGMA, c(1,2), mean)

  # IMPULSE RESPONSES - KOOP/KOROBILIS a la Primiceri (2005)
  J <- matrix(0, M*plag, M)
  J[1:M,1:M] <- diag(M)

  compMat <- matrix(0, M*plag, M*plag)
  if(plag==1) compMat <- t(A_draw_mean) else {
    for(j in 1:(p-1)) {
      compMat[(j*M+1):(M*(j+1)),(M*(j-1)+1):(j*M)] <- diag(M)
    }
  }
  bbtemp <- A_draw_mean[1:(M*plag),] # kick intercept
  splace <- 0
  for(ii in 1:plag) {
    for(iii in 1:M) {
      compMat[iii,((ii-1)*M+1):(ii*M)] <- t(bbtemp[(splace+1):(splace+M),iii])
    }
    splace <- splace+M
  }

  # identification
  shock <- t(chol(SIGMA_mean))
  diagonal <- diag(diag(shock))
  shock <- solve(diagonal)%*%shock # unit initial shock

  nhor <- 60
  impresp <- array(0, c(M, M, nhor))
  impresp[1:M,1:M, 1] <- shock
  compMati <- compMat
  for(j in 2:nhor) {
    impresp[,,j] <- t(J) %*% compMati %*% J %*% shock
    compMati <- compMati %*% compMat
  }

  true.list <- list(A1.true=A1.mean, B1.true=B1.mean, A.true=A, B.true=B, vol.true=vol_true, SIGMA1.true=SIGMA_mean, SIGMA.true=SIGMA)

  return(list(Yraw=Yraw,Wraw=Wraw,true.list=true.list,impresp.true=impresp))
}

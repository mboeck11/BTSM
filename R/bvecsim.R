#' @name bvec.sim
#' @title Simulating a Vector Error Correction Model
#' @description This function is used to produce simulated realizations which follow a Vector error correction model. It will also automatically simulate coefficients. All parameters can also be set by the user.
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
#' len <- 1000
#' M <- 3
#' plag <- 1
#' include <- "none"
#' B <- rbind(c(0.7, 0, 0, 0.4), c(0.6, -0.5, 0.2, 0), c(0.2, 0, 0.4, 0.9))
#' beta <- c(1,-1,-1)
#' yb <- bvec.sim(len=len, M=M, B=B, beta=beta, plag=plag, include=include)
#' @importFrom stats rnorm
#' @importFrom stochvol svsim
#' @importFrom mnormt rmnorm
#' @export
bvec.sim <- function(len, M, B, beta, plag=1, include="none"){
  if(!is.matrix(B)){
    stop("Please provide with argument 'B' a matrix.")
  }
  ninc<- switch(include, "none"=0, "const"=1, "trend"=1, "both"=2)
  incVal<- switch(include, "none"=NULL, "const"="const", "trend"="trend", "both"=c("const","trend"))
  innov <- rbind(matrix(0, plag+1, M), rmnorm(len, varcov=diag(M)))
  rownames(B) <- paste("Equ x", 1:M, ":",sep="")

  smallT <- len-plag-1
  if(is.vector(beta)){
    if(length(beta)==M-1) beta <- c(1, -beta)
    tBETA<-matrix(beta, nrow=1)
    r <- 1
  } else {
    if(nrow(beta)!=M) stop("beta should have k rows and r cols")
    r <- ncol(beta)
    tBETA <- t(beta)
  }
  if(dim(B)[1]!=M || dim(B)[2]!=r+ninc+M*plag){
    stop("Please provide matrix B with dimensions M times (r+cons+trend+M*p).")
  }
  delta <- 1
  tBETAo <- tBETA
  if(include=="const") tBETAo <- cbind(tBETA,delta)

  colnames(B) <- c(paste0("cointegration vector ",seq(1,r)),incVal,paste0("coef ",seq(1,M)))
  Bmat <- matrix(0,nrow=M,ncol=r+2+plag*M)
  cols <- c(paste0("cointegration vector ",seq(1,r)),"const","trend",paste0("coef ",seq(1,M)))
  colnames(Bmat) <- cols
  for(c in cols) if(c%in%colnames(B)) Bmat[,c] <- B[,c]

  Yb<-matrix(0, len, M)
  trend<-c(rep(NA, len-smallT),1:smallT)
  Eb<-matrix(NA, len, r)

  for(tt in (plag+2):len){
    ECT<-Bmat[,1:r,drop=FALSE]%*%tBETA%*%matrix(Yb[tt-1,], ncol=1)
    Yb[tt,]<-rowSums(cbind(t(Yb[tt-1,,drop=FALSE]),
                           Bmat[,r+1,drop=FALSE],
                           Bmat[,r+2,drop=FALSE]*trend[tt],
                           ECT,
                           Bmat[,-c(1:(r+2))]%*%matrix(t(Yb[tt-c(1:plag),,drop=FALSE]-Yb[tt-c(2:(plag+1)),,drop=FALSE]), ncol=1),
                           t(innov[tt,,drop=FALSE])))
    # Ybo <- matrix(Yb[tt-1,],ncol=1)
    # if(include=="const") Ybo <- rbind(Ybo,1)
    # Eb[i,]<-tBETAo%*%Ybo
  }

  return(Yb)
}

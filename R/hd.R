#' @export
"hd" <- function(x, R=NULL, verbose=TRUE){
  UseMethod("hd", x)
}

#' @name hd
#' @title Historical Decomposition
#' @description A function that calculates historical decomposition (HD) of the time series and the structural error.
#' @method hd bvar.irf
#' @usage hd(x, R=NULL, verbose=TRUE)
#' @param x an item fitted by \code{IRF}.
#' @param R If \code{NULL} and the \code{irf.bgvar} object has been fitted via sign restrictions, the rotation matrix is used that minimizes the distance to the median impulse responses at the posterior median.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @details To save computational time as well as due to storage limits, both functions are based on the posterior median (as opposed to calculating HDs and the structural error for each draw of the MCMC chain). In case the shock has been identified via sign restrictions, a rotation matrix has to be selected to calculate both statistics. If not specified otherwise (via \code{R}), the algorithm searches for 50 rotation matrices that fulfill the sign restrictions at the \emph{posterior median} of the coefficients and then singles out the rotation matrix that minimizes the distance to the median of the impulse responses as suggested in Fry and Pagan (2011).
#' @export
hd.bvar.irf<-function(x, R=NULL, verbose=TRUE){
  start.hd <- Sys.time()
  if(verbose) cat("\nStart computing historical decomposition of Bayesian Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  xglobal <- x$model.obj$xglobal
  plag    <- x$model.obj$plag
  ident   <- x$ident
  Traw    <- nrow(xglobal)
  bigK    <- ncol(xglobal)
  xdat    <- xglobal[(plag+1):Traw,,drop=FALSE]
  bigT    <- nrow(xdat)
  Amat    <- x$struc.obj$Amat
  Smat    <- x$struc.obj$Smat
  varNames<- colnames(xglobal)
  cons    <- x$model.obj$cons
  trend   <- x$model.obj$trend
  if(!is.null(R)){
    R<-x$struc.obj$Rmed
  }else{
    R<-diag(bigK)
  }
  rownames(R) <- colnames(R) <- varNames
  #------------------------checks-------------------------------------------------------------------#
  if(ident=="girf"){
    message("Historical decomposition of the time series not implemented for GIRFs since cross-correlation is unequal to zero (and hence decompositions do not sum up to original time series).")
    return(list(hd_array=NA,struc.shock=vv,xglobal=xglobal) )
  }
  #------ initialize objects -----------------------------------------------------------------------#
  struc_post <- array(NA,dim=c(bigT,bigK))
  hd_array   <- array(0,c(bigK,bigT,(bigK+2+cons+trend)))
  if(cons & !trend){
    dimnames(hd_array)<-list(rownames(x), NULL, c(paste("contribution of shock to", c(varNames)),"constant","initial cond.","residual"))
  }else if(!cons & trend){
    dimnames(hd_array)<-list(rownames(x), NULL, c(paste("contribution of shock to", c(varNames)),"trend","initial cond.","residual"))
  }else if(cons & trend){
    dimnames(hd_array)<-list(rownames(x), NULL, c(paste("contribution of shock to", c(varNames)),"constant","trend","initial cond.","residual"))
  }else{
    dimnames(hd_array)<-list(rownames(x), NULL, c(paste("contribution of shock to", c(varNames)),"initial cond.","residual"))
  }
  #------------------------------------------------------------------------------------------------#
  Rinv       <- solve(R)
  Sigchol_u  <- t(chol(Smat))
  Sigcholinv <- solve(Sigchol_u)

  vv <- matrix(0,bigT,bigK,dimnames=list(NULL,varNames))
  YY <- xglobal[(plag+1):Traw,]
  XX <- .mlag(xglobal,plag)
  XX <- XX[(plag+1):nrow(XX),]
  if(cons)  XX <- cbind(XX,1)
  if(trend) XX <- cbind(XX,seq(1,bigT))

  strMat <- Rinv%*%Sigcholinv

  for (tt in 1:bigT){
    Yhat <- strMat%*%YY[tt,]
    Xhat <- strMat%*%t(XX[tt,]%*%Amat)

    #PHI <- Rinv%*%Sigcholinv%*%t(Amat)
    vv[tt,] <- Yhat-Xhat
  }
  #Start historical decompositions -------------------------------------------------------------------------------#
  if(verbose) cat("Start computing HDs...\n")
  HDshock_big <- array(0,c(plag*bigK,bigT,bigK))
  HDinit_big  <- matrix(0,plag*bigK,bigT)
  HDshock     <- array(0,c(bigK,bigT,bigK))
  HDinit      <- matrix(0,bigK,bigT)
  if(cons){
    HDconst_big <- matrix(0,plag*bigK,bigT)
    HDconst     <- matrix(0,bigK,bigT)
  }
  if(trend){
    HDtrend_big <- matrix(0,plag*bigK,bigT)
    HDtrend     <- matrix(0,bigK,bigT)
  }

  #NOTE FOR MARTIN: IF SIGN RESTRICTIONS: R = ROTATION ELSE R = diag(M)
  solveA <- (Sigchol_u%*%R) #Depends on identification, if Cholesky then solveA = t(chol(SIGMA)), where SIGMA is the VC of the global model
  eps <- (YY-XX%*%Amat)%*%t(solve(solveA)) #Atilda is the matrix of autoregressive coefficients of the global model
  Fcomp <- .gen_compMat(Amat, bigK, plag)$Cm

  invA_big <- matrix(0,bigK*plag,bigK)  #M is the number of endogenous variables ; p is the number of lags
  invA_big[1:bigK,] <- solveA
  Icomp <- cbind(diag(bigK),matrix(0,bigK,(plag-1)*bigK))
  for (nn in 2:bigT){
    for (jj in 1:bigK){
      eps_big <- matrix(0,bigK,1)
      eps_big[jj,] <- eps[nn,jj]
      HDshock_big[,nn,jj] <- (invA_big)%*%eps_big+Fcomp%*%HDshock_big[,nn-1,jj]
      HDshock[,nn,jj] <- Icomp%*%HDshock_big[,nn,jj]
    }
    #Initial value
    HDinit_big[,1] <- XX[1,1:(plag*bigK)]
    HDinit[,1] <- Icomp%*%HDinit_big[,1]
    HDinit_big[,nn] <- Fcomp%*%HDinit_big[,nn-1]
    HDinit[,nn] <- Icomp%*%HDinit_big[,nn]

    #Constant
    if(cons){
      CC <- matrix(0,bigK*plag,1)
      CC[1:bigK] <- Amat[(bigK*plag)+1,]
      HDconst_big[,nn] <- CC+Fcomp%*%HDconst_big[,nn-1]
      HDconst[,nn] <- Icomp%*%HDconst_big[,nn]
    }

    # Trend
    if(trend){
      TT <- matrix(0,bigK*plag,1)
      TT[1:bigK] <- t(Amat)[(bigK*plag)+2,]
      HDtrend_big[,nn] <- TT+Fcomp%*%HDtrend_big[,nn-1]
      HDtrend[,nn] <- Icomp%*%HDtrend_big[,nn]
    }
  }
  hd_array[,,1:bigK]   <- HDshock_big
  if(cons)  hd_array[,,(bigK+1)] <- HDconst_big
  if(trend) hd_array[,,(bigK+1+trend)] <- HDtrend_big
  hd_array[,,(bigK+2+trend)] <- HDinit_big
  hd_array[,,(bigK+3+trend)] <- (t(xdat)-apply(hd_array,c(1,2),sum)) # residual part
  #----------------------------------------------------------------------------------#
  hd_array <- aperm(hd_array,c(2,1,3))
  out      <- structure(list(hd_array=hd_array,struc_shock=vv,xglobal=xdat, R=NULL), class="bvar.hd")
  if(verbose) cat(paste("Size of object:", format(object.size(out),unit="MB")))
  end.hd <- Sys.time()
  diff.hd <- difftime(end.hd,start.hd,units="mins")
  mins.hd <- round(diff.hd,0); secs.hd <- round((diff.hd-floor(diff.hd))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.hd," ",ifelse(mins.hd==1,"min","mins")," ",secs.hd, " ",ifelse(secs.hd==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @method print bvar.hd
#' @export
print.bvar.hd <- function(x, ...){
  cat("---------------------------------------------------------------------------------------")
  cat("\n")
  cat("Object contains historical decomposition of object estimated with 'bgvar':")
  cat("\n")
  cat(paste0("Size of hd_array containing historical decompositions: ",dim(x$hd_array)[[1]]," x ",dim(x$hd_array)[[2]]," x ",dim(x$hd_array)[[3]],"."))
  cat("\n")
  cat(paste0("Size of struc_shock containing structural errors: ",dim(x$struc_shock)[[1]]," x ",dim(x$struc_shock)[[2]],"."))
  cat("\n")
  cat("Identification scheme: ")
  if(is.null(x$R)){
    cat("Short-run restrictions via Cholesky decomposition.")
  }else{
    cat("Sign-restrictions.")
  }
  cat("\n")
  cat(paste0("Size ob object: ",format(object.size(x),unit="MB")))
  cat("\n")
  cat("---------------------------------------------------------------------------------------")

  return(invisible(x))
}

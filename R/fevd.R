#' @export
"fevd" <- function(x, R=NULL, var.slct=NULL, verbose=TRUE){
  UseMethod("fevd", x)
}

#' @name fevd
#' @title Forecast Error Variance Decomposition
#' @description This function calculates the forecast error variance decomposition (FEVDs) for Cholesky and sign-identified shocks.
#' @usage fevd(x, R=NULL, var.slct=NULL, verbose=TRUE)
#' @details Since the calculations are very time consuming, the FEVDs are based on the posterior median only (as opposed to calculating FEVDs for each MCMC sweep). In case the underlying shock has been identified via sign restrictions, the rotation matrix corresponds to the one that fulfills the sign restrictions at the posterior median of the estimated coefficients. More precisely, the algorithm searches for 50 rotation matrices that fulfill the sign restrictions at the \emph{posterior median} of the coefficients and then singles out the rotation matrix that minimizes the distance to the median of the impulse responses as suggested in Fry and Pagan (2011).
#' @param x an object of class \code{bvar.irf}.
#' @param R If \code{NULL} and the \code{x} has been fitted via sign restrictions, the rotation matrix is used that minimizes the distance to the median impulse responses at the posterior median.
#' @param var.slct character vector that contains the variables for which forecast error variance decomposition should be performed. If \code{NULL} the FEVD is computed for the whole system, which is very time consuming.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @return Returns a list with two elements \itemize{
#' \item{\code{FEVD}}{  an array of size (K times horizon times N), where K are all variables in the system, horizon is the specified impulse response horizon and N is the size of the decomposed structural variables (if \code{var.slct=NULL} then K=N).}
#' \item{\code{xglobal}}{ used data of the model.}
#' }
#' @export
fevd.bvar.irf <- function(x, R=NULL, var.slct=NULL, verbose=TRUE){
  start.fevd <- Sys.time()
  if(verbose) cat("\nStart computing forecast error variance decomposition of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  xglobal <- x$model.obj$xglobal
  plag    <- x$model.obj$plag
  ident   <- x$ident
  Traw    <- nrow(xglobal)
  bigK    <- ncol(xglobal)
  xdat    <- xglobal[(plag+1):Traw,,drop=FALSE]
  bigT    <- nrow(x)
  Amat    <- x$struc.obj$Amat
  Smat    <- x$struc.obj$Smat
  horizon <- dim(x$posterior)[1]
  varNames<- colnames(xglobal)
  shock   <- x$shock
  sign.constr <- x$sign.constr
  if(ident=="sign"){
    if(verbose) cat("Identification scheme: Sign-restrictions provided.\n")
  }else if(ident=="chol"){
    if(verbose) cat("Identification scheme: Short-run restrictions via Cholesky decomposition.\n")
  }else if(ident=="proxy"){
    if(verbose)
      cat("Identification schem: Identification via proxy variable.\n")
  }
  #-------------------- some checks ------------------------------------------------------------------#
  if(!ident%in%c("sign","chol","proxy")){
    stop("FEVD implemented for shocks identified via cholesky ordering, sign restrictions and external instruments only.")
  }
  if(!is.null(var.slct)){
    if(!all(var.slct%in%varNames)){
      stop("One of the variables you want to decompose is not contained in the system. Please re-specify!")
    }
  }
  if(is.null(var.slct)){
    if(verbose) cat("FEVD computed for all variables.\n\n")
    var.slct<-varNames
  }else{
    var.print <- var.slct[1]
    if(length(var.slct)>1) for(kk in 2:length(var.slct)) var.print <- paste(var.print,", ",var.slct[kk],sep="")
    if(verbose) cat(paste("FEVD computed for the following variables: ",var.print,".\n",sep=""))
  }
  if(ident=="sign" && is.null(R)){
    R <- x$struc.obj$Rmed
  }else{
    R<-diag(bigK)
  }
  rownames(R) <- colnames(R) <- varNames
  #----------------------------------------------------------------------------------------------------#
  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+horizon+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  temp <- .gen_compMat(Amat, bigK, plag)
  Cm   <- temp$Cm
  Jm   <- temp$Jm
  Cmat <- diag(bigK*plag)
  for (ihor in (plag+2):(plag+horizon+1)){
    Cmat  <- Cmat%*%Cm
    PHIx[,,ihor]  <- t(Jm)%*%Cmat%*%Jm
  }
  PHI  <-  PHIx[,,(plag+1):(plag+horizon+1)]
  #----------------------------------------------------------------------------------------------------#
  if(verbose) cat("Start computing FEVDs...\n")
  vslct <- diag(bigK)
  P0G   <- t(chol(Smat))
  if(ident == "proxy"){
    P0G <- x$struc.obj$Rmed
  }

  FEVDres  <-  array(0,dim=c(bigK,length(var.slct),horizon+1))
  dimnames(FEVDres) <- list(varNames,paste("Decomp. of",var.slct),0:horizon)

  for(zz in 1:length(var.slct)){
    eslct <-matrix(0,bigK,1);rownames(eslct) <- varNames
    eslct[var.slct[zz],1] <- 1

    num  <-  matrix(0,bigK,horizon+1)
    den  <-  matrix(0,bigK,horizon+1)

    N <- 1
    while (N<=horizon+1){
      for (l in 1:N){
        acc1  <-  t((t(eslct)%*%R%*%PHI[,,l]%*%P0G%*%vslct)^2)
        num[,N]  <-  num[,N] + acc1
        acc2  <-  (t(eslct)%*%R%*%PHI[,,l]%*%Smat%*%t(R%*%PHI[,,l])%*%eslct)
        den[,N]  <-  den[,N] + matrix(1,bigK,1)*as.numeric(acc2)
      }
      FEVDres[,paste("Decomp. of",var.slct[zz]),N]  <-  (num[,N])/den[,N]
      N <- N+1
    }
  }
  #------------------------------------------------------------------------------------------------------
  out <- structure(list(FEVD=FEVDres,
                        xglobal=xglobal,
                        R=R),
                   class="bvar.fevd", type="fevd")
  if(verbose) cat(paste("\nSize of FEVD object: ", format(object.size(FEVDres),unit="MB")))
  end.fevd <- Sys.time()
  diff.fevd <- difftime(end.fevd,start.fevd,units="mins")
  mins.fevd <- round(diff.fevd,0); secs.fevd <- round((diff.fevd-floor(diff.fevd))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.fevd," ",ifelse(mins.fevd==1,"min","mins")," ",secs.fevd, " ",ifelse(secs.fevd==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @method print bvar.fevd
#' @export
print.bvar.fevd <- function(x, ...){
  cat("---------------------------------------------------------------------------------------")
  cat("\n")
  cat("Object contains forecast error variance decomposition of object estimated with 'bvar':")
  cat("\n")
  cat(paste0("Size of FEVD containing forecast error variance decompositions: ",dim(x$FEVD)[[1]]," x ",dim(x$FEVD)[[2]]," x ",dim(x$FEVD)[[3]],"."))
  cat("\n")
  if(attributes(x)$type=="fevd"){
    cat("Identification scheme: ")
    if(is.null(x$R)){
      cat("Short-run restrictions via Cholesky decomposition.")
    }else{
      cat("Sign-restrictions.")
    }
  }else if(attributes(x)$type=="gfevd"){
    cat("Identification scheme: Generalized - no identification scheme employed.")
  }
  cat("\n")
  cat(paste0("Size ob object: ",format(object.size(x),unit="MB")))
  cat("\n")
  cat("---------------------------------------------------------------------------------------")

  return(invisible(x))
}

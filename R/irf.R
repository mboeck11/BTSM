#' @name irf
#' @title Impulse Response Function
#' @usage irf(x, n.ahead=24, ident=NULL, shockinfo=NULL, save.store=FALSE,
#'    applyfun=NULL, cores=NULL, verbose=TRUE,
#'    quantiles=c(.05,.10,.16,.50,.84,.90,.95), ...)
#' @param x object of class \code{bvar}.
#' @param n.ahead Forecasting horizon.
#' @param ident Preferred identification scheme.
#' @param shockinfo Dataframe with specified details on the shock.
#' @param proxy Matrix with proxy variables.
#' @param save.store If set to \code{TRUE} the full posterior is returned. Default is set to \code{FALSE} in order to save storage.
#' @param applyfun Allows for user-specific apply function, which has to have the same interface than \code{lapply}. If \code{cores=NULL} then \code{lapply} is used, if set to a numeric either \code{parallel::parLapply()} is used on Windows platforms and \code{parallel::mclapply()} on non-Windows platforms.
#' @param cores Specifies the number of cores which should be used. Default is set to \code{NULL} and \code{applyfun} is used.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @export
"irf" <- function(x, n.ahead=24, ident=NULL, shockinfo=NULL, proxy=NULL, save.store=FALSE, applyfun=NULL, cores=NULL, verbose=TRUE,
                  quantiles=c(.05,.10,.16,.50,.84,.90,.95), ...){
  #------------------------------ do checks ---------------------------------------------------------------#
  .irf.checks(x=x,n.ahead=n.ahead,ident=ident,shockinfo=shockinfo,proxy=proxy,save.store=save.store,quantiles=quantiles)
  #------------------------------ use irf method ----------------------------------------------------------#
  UseMethod("irf", x)
}

#' @export
irf.bvar <- function(x, n.ahead=24, ident=NULL, shockinfo=NULL, proxy=NULL, save.store=FALSE, applyfun=NULL, cores=NULL, verbose=TRUE,
                     quantiles=c(.05,.10,.16,.50,.84,.90,.95), ...){
  start.irf <- Sys.time()
  cat("\nStart computing impulse response functions of Bayesian Vector Autoregression.\n\n")
  if(ident=="chol-shortrun"){
    if(verbose)
      cat("Identification scheme: Short-run identification via Cholesky decomposition.\n")
  }else if(ident=="chol-longrun"){
    if(verbose)
      cat("Identification schem: Long-run identification via Cholesky decomposition.\n")
  }else if(ident=="girf"){
    if(verbose)
      cat("Identification scheme: Generalized impulse responses.\n")
  }else if(ident=="sign"){
    if(verbose)
      cat("Identification scheme: identification via sign-restrictions.\n")
  }else if(ident=="proxy"){
    if(verbose)
      cat("Identification schem: Identification via proxy variable.\n")
  }
  if(verbose) cat(paste("Start impulse response analysis on ", ifelse(is.null(cores),1,cores), " cores", " (",x$args$thindraws," stable draws in total).",sep=""),"\n")
  out <- .irf.generator(x,n.ahead=n.ahead,ident=ident,shockinfo=shockinfo,proxy=proxy,save.store=save.store,applyfun=applyfun,quantiles=quantiles,verbose=verbose)
  cat(paste("\nSize of irf object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @export
irf.tvpbvar <- function(x, n.ahead=24, ident=NULL, shockinfo=NULL, proxy=NULL, save.store=FALSE, applyfun=NULL, cores=NULL, verbose=TRUE,
                        quantiles=c(.05,.10,.16,.50,.84,.90,.95), period="med", ...){
  start.irf <- Sys.time()
  cat("\nStart computing impulse response functions of Time-varying Parameter Bayesian Vector Autoregression.\n\n")
  if(ident=="chol-shortrun"){
    if(verbose)
      cat("Identification scheme: Short-run identification via Cholesky decomposition.\n")
  }else if(ident=="chol-longrun"){
    if(verbose)
      cat("Identification schem: Long-run identification via Cholesky decomposition.\n")
  }else if(ident=="girf"){
    if(verbose)
      cat("Identification scheme: Generalized impulse responses.\n")
  }else if(ident=="sign"){
    if(verbose)
      cat("Identification scheme: identification via sign-restrictions.\n")
  }else if(ident=="proxy"){
    if(verbose)
      cat("Identification schem: Identification via proxy variable.\n")
  }
  #------------ check ----------------------------#
  if(length(period)!=1) stop("Please provide argument 'period' with length 1.")
  #------------ get data -------------------------#
  Y    <- x$args$Y
  bigT <- nrow(Y)
  M    <- ncol(Y)
  thindraws <- x$args$thindraws
  #-------------------------------------------------------------------------------------------------------#
  if(verbose) cat(paste("Start impulse response analysis on ", cores, " cores", " (",thindraws," stable draws in total).",sep=""),"\n")
  # median response
  if(period=="med"){
    x.med <- x
    x.med$store$A_store <- apply(x.med$store$A_store,c(1,3,4),median)
    out <- .irf.generator(x.med,n.ahead=n.ahead,ident=ident,shockinfo=shockinfo,proxy=proxy,save.store=save.store,applyfun=applyfun,cores=cores,verbose=verbose)
  }else if(period%in%seq(bigT)){
    tt <- as.numeric(period)
    if(verbose) cat(paste0("Time point: ", tt, " of ", bigT,".\n"))
    x.t <- x
    x.t$store$A_store <- x.t$store$A_store[,tt,,]
    x.t$store$Smed_store <- x.t$store$S_store[,tt,,]
    out <- .irf.generator(x.t,n.ahead=n.ahead,ident=ident,shockinfo=shockinfo,proxy=proxy,save.store=save.store,applyfun=applyfun,cores=cores,quantiles=quantiles,verbose=verbose)$posterior
  }else if(period=="full"){
    x.med <- x
    x.med$store$A_store <- apply(x.med$store$A_store,c(1,3,4),median)
    out <- .irf.generator(x.med,n.ahead=n.ahead,ident=ident,shockinfo=shockinfo,proxy=proxy,save.store=save.store,applyfun=applyfun,cores=cores,quantiles=quantiles,verbose=verbose)
    out$posterior.full <- array(NA,c(bigT,n.ahead,M,M,length(quantiles)),dimnames=c(list(NULL),dimnames(out$posterior)))
    for(tt in 1:bigT){
      x.t <- x
      x.t$store$A_store <- x.t$store$A_store[,tt,,]
      x.t$store$Smed_store <- x.t$store$S_store[,tt,,]
      out$posterior.full[tt,,,,] <- .irf.generator(x.t,n.ahead=n.ahead,ident=ident,shockinfo=shockinfo,proxy=proxy,save.store=save.store,applyfun=applyfun,quantiles=quantiles,cores=cores,verbose=FALSE)$posterior
    }
  }
  ## bind together somehow
  cat(paste("\nSize of irf object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @export
irf.bvec <- function(x, n.ahead=24, ident=NULL, shockinfo=NULL, proxy=NULL, save.store=FALSE, applyfun=NULL, cores=NULL, verbose=TRUE,
                     quantiles=c(.05,.10,.16,.50,.84,.90,.95), ...){
  start.irf <- Sys.time()
  cat("\nStart computing impulse response functions of Bayesian Vector Error Correction Model.\n\n")
  if(ident=="chol-shortrun"){
    if(verbose)
      cat("Identification scheme: Short-run identification via Cholesky decomposition.\n")
  }else if(ident=="chol-longrun"){
    if(verbose)
      cat("Identification schem: Long-run identification via Cholesky decomposition.\n")
  }else if(ident=="girf"){
    if(verbose)
      cat("Identification scheme: Generalized impulse responses.\n")
  }else if(ident=="sign"){
    if(verbose)
      cat("Identification scheme: identification via sign-restrictions.\n")
  }else if(ident=="proxy"){
    if(verbose)
      cat("Identification schem: Identification via proxy variable.\n")
  }
  if(verbose) cat(paste("Start impulse response analysis on ", cores, " cores", " (",x$args$thindraws," stable draws in total).",sep=""),"\n")
  out <- .irf.generator(x,n.ahead=n.ahead,ident=ident,shockinfo=shockinfo,proxy=proxy,save.store=save.store,applyfun=applyfun,quantiles=quantiles,verbose=verbose)
  cat(paste("\nSize of irf object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  if(verbose) cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @export
irf.bivar <- function(x, n.ahead=24, ident=NULL, shockinfo=NULL, proxy=NULL, save.store=FALSE, applyfun=NULL, cores=NULL, verbose=TRUE,
                      quantiles=c(.05,.10,.16,.50,.84,.90,.95), eval.q=NULL, ...){
  start.irf <- Sys.time()
  cat("\nStart computing impulse response functions of Bayesian Interacted Vector Autoregression.\n\n")
  if(ident=="chol-shortrun"){
    if(verbose){
      cat("Identification scheme: Short-run identification via Cholesky decomposition.\n")
    }
  }else if(ident=="chol-longrun"){
    if(verbose){
      cat("Identification schem: Long-run identification via Cholesky decomposition.\n")
    }
  }else if(ident=="girf"){
    if(verbose){
      cat("Identification scheme: Generalized impulse responses.\n")
    }
  }else if(ident=="sign"){
    if(verbose){
      cat("Identification scheme: identification via sign-restrictions.\n")
    }
  }else if(ident=="proxy"){
    if(verbose){
      cat("Identification schem: Identification via proxy variable.\n")
    }
  }
  #------------ get data -------------------------#
  Y         <- x$args$Y
  D         <- coredata(x$xint)
  Jt        <- x$store$J_store
  M         <- ncol(x$args$Data)
  Ki        <- ncol(D)
  Ki1       <- Ki+1
  thindraws <- x$args$thindraws
  k         <- dim(x$store$Atilde_store)[[2]]
  plag      <- x$args$plag
  bigT      <- nrow(x$xglobal)-plag
  cons      <- x$args$cons
  trend     <- x$args$trend
  X         <- .mlag(x$args$Data,plag)
  X         <- X[(plag+1):nrow(X),,drop=FALSE]
  if(cons) X <- cbind(X,1)
  if(trend) X <- cbind(X,seq(1,bigT))
  #---------- extra checks------------------------#
  if(length(eval.q)!=Ki)
    stop("Please provide 'eval.q' with same length as number of columns in interaction matrix D.")
  #------ transform to reduced-form model --------#
  Dval <- c(1,apply(D,2,quantile,eval.q))

  # old stuff
  a0tilde_store   <- x$store$a0tilde_store
  a1tilde_store   <- x$store$a1tilde_store
  Phitilde_store  <- x$store$Phitilde_store
  volatilde_store <- x$store$volatilde_store

  # new stuff
  A_store    <- array(NA, c(thindraws,M*plag+cons+trend,M))
  S_store    <- array(NA, c(thindraws,bigT,M,M))
  res_store  <- array(NA, c(thindraws,bigT,M))

  # loop over draws
  for(irep in 1:thindraws){
    J <- diag(M)
    for(mm in 2:M){
      for(jj in 1:(mm-1)){
        for(kk in 1:Ki1) J[mm,jj] <- J[mm,jj]+Jt[irep,mm,jj,kk]*Dval[kk]
      }
    }
    Jinv <- solve(J)
    a0 <- a1 <- matrix(0,1,M)
    Phi <- matrix(0,M*plag,M)
    S <- array(NA,c(bigT,M,M))
    if(cons)
      for(kk in 1:Ki1) a0 <- a0+a0tilde_store[irep,kk,]*Dval[kk] else a0 <- NULL
    if(trend)
      for(kk in 1:Ki1) a1 <- a1+a1tilde_store[irep,kk,]*Dval[kk] else a1 <- NULL
    for(pp in 1:plag){
      for(kk in 1:Ki1) Phi[((pp-1)*M+1):(pp*M),] <- Phi[((pp-1)*M+1):(pp*M),] + Phitilde_store[[pp]][irep,,,kk]*Dval[kk]
    }
    # multiply with Jinv
    if(cons) a0 <- t(Jinv%*%t(a0))
    if(trend) a1 <- t(Jinv%*%t(a1))
    Phi <- t(Jinv%*%t(Phi))
    for(tt in 1:bigT){
      S[tt,,] <- Jinv%*%diag(volatilde_store[irep,tt,])%*%t(Jinv)
    }
    A_store[irep,,]   <- rbind(Phi,a0,a1)
    S_store[irep,,,]  <- S
    res_store[irep,,] <- Y - X%*%A_store[irep,,]
  }

  x$store$A_store    <- A_store
  x$store$S_store    <- S_store
  x$store$res_store  <- res_store
  x$store$Smed_store <- apply(S_store,c(1,3,4),median)
  x$store$res_store  <- res_store
  if(verbose) cat(paste("Start impulse response analysis on ", cores, " cores", " (",thindraws," stable draws in total).",sep=""),"\n")
  out <- .irf.generator(x=x,n.ahead=n.ahead,ident=ident,shockinfo=shockinfo,proxy=proxy,save.store=save.store,applyfun=applyfun,verbose=verbose)
  if(verbose) cat(paste("\nSize of irf object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  if(verbose) cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @name .irf.generator
#' @noRd
#' @importFrom stats median
#' @importFrom stringr str_pad
#' @importFrom utils object.size
.irf.generator <- function(x,n.ahead=24,ident=NULL,shockinfo=NULL,proxy=NULL,save.store=FALSE,applyfun=NULL,cores=NULL,quantiles=c(.05,.10,.16,.50,.84,.90,.95),verbose=TRUE){
  #------------------------------ get stuff -------------------------------------------------------#
  plag        <- x$args$plag
  xglobal     <- x$args$Data
  Traw        <- nrow(xglobal)
  bigK        <- ncol(xglobal)
  bigT        <- Traw-plag
  bigQ        <- length(quantiles)
  A_large     <- x$store$A_store
  S_large     <- x$store$Smed_store
  E_large     <- x$store$res_store
  xdat        <- xglobal[(plag+1):Traw,,drop=FALSE]
  thindraws   <- x$args$thindraws
  varNames    <- colnames(xglobal)
  MaxTries    <- 7500
  Rmed        <- NULL
  epsmed      <- NULL
  type        <- NULL
  rot.nr      <- NULL
  #------------------------------ prepare applyfun --------------------------------------------------------#
  if(is.null(applyfun)) {
    applyfun <- if(is.null(cores)) {
      lapply
    } else {
      if(.Platform$OS.type == "windows") {
        cl_cores <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl_cores))
        function(X, FUN, ...) parallel::parLapply(cl = cl_cores, X, FUN, ...)
      } else {
        function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores =
                                                   cores)
      }
    }
  }
  if(is.null(cores)) {cores <- 1}
  #------------------------------ assign irf function  ----------------------------------------------------#
  if(ident=="sign"){
    irf.fun<-.irf.sign.zero
    select_shocks <- which(varNames%in%shockinfo$shock)
    scale <- shockinfo$scale
    shock.nr <- length(select_shocks)
    if(shock.nr == 0)
      stop("Please provide shock in the dataset. Respecify 'shockinfo' argument.")
    # adjust for rationality conditions
    if(any(shockinfo$sign=="ratio.H")){
      idx <- which(shockinfo$sign=="ratio.H")
      for(ii in idx){
        Kshock <- nrow(shockinfo)
        Mshock <- as.numeric(shockinfo$horizon[ii])
        shockinfo[(Kshock+1):(Kshock+2),] <- NA
        shockinfo$shock[(Kshock+1):nrow(shockinfo)] <- rep(shockinfo$shock[ii],2)
        shockinfo$restrictions[(Kshock+1):nrow(shockinfo)] <- c(shockinfo$restrictions[ii], strsplit(shockinfo$restrictions[ii],"_")[[1]][1])
        shockinfo$sign[(Kshock+1):nrow(shockinfo)] <- c("0","-1")
        shockinfo$horizon[(Kshock+1):nrow(shockinfo)] <- c(1,Mshock)
        shockinfo$scale[(Kshock+1):nrow(shockinfo)] <- rep(shockinfo$scale[ii],2)
      }
      shockinfo <- shockinfo[-idx,]
      rownames(shockinfo)<-seq(1,nrow(shockinfo))
    }
    if(any(shockinfo$sign=="ratio.avg")){
      idx <- which(shockinfo$sign=="ratio.avg")
      for(ii in idx){
        Kshock <- nrow(shockinfo)
        Mshock <- as.numeric(shockinfo$horizon[ii])
        shockinfo[(Kshock+1):(Kshock+Mshock),] <- NA
        shockinfo$shock[(Kshock+1):nrow(shockinfo)] <- rep(shockinfo$shock[ii],Mshock)
        shockinfo$restrictions[(Kshock+1):nrow(shockinfo)] <- c(shockinfo$restrictions[ii],rep(strsplit(shockinfo$restrictions[ii],"_")[[1]][1],Mshock-1))
        shockinfo$sign[(Kshock+1):nrow(shockinfo)] <- c("0",rep(-1/(Mshock-1),Mshock-1))
        shockinfo$horizon[(Kshock+1):nrow(shockinfo)] <- seq(1,Mshock)
        shockinfo$scale[(Kshock+1):nrow(shockinfo)] <- rep(shockinfo$scale[ii],Mshock)
      }
      shockinfo <- shockinfo[-idx,]
      rownames(shockinfo)<-seq(1,nrow(shockinfo))
    }
  }else if(ident%in%c("chol","chol-shortrun")){
    select_shocks <- which(varNames%in%shockinfo$shock)
    shock.nr <- length(select_shocks)
    if(shock.nr == 0)
      stop("Please provide shock in the dataset. Respecify 'shockinfo' argument.")
    scale <- shockinfo$scale
    irf.fun <- .irf.chol
    type="short-run"
  }else if(ident=="girf"){
    irf.fun <- .irf.girf
  }else if(ident=="chol-longrun"){
    select_shocks <- which(varNames%in%shockinfo$shock)
    shock.nr <- length(select_shocks)
    if(shock.nr == 0)
      stop("Please provide shock in the dataset. Respecify 'shockinfo' argument.")
    scale <- shockinfo$scale
    irf.fun <- .irf.chol
    type="short-run"
  }else if(ident=="proxy"){
    select_shocks <- which(varNames%in%shockinfo$shock)
    shock.nr <- length(select_shocks)
    if(shock.nr == 0)
      stop("Please provide shock in the dataset. Respecify 'shockinfo' argument.")
    scale <- shockinfo$scale
    irf.fun <- .irf.proxy
    if(nrow(proxy)==Traw){
      proxy <- proxy[(plag+1):Traw,,drop=FALSE]
    }else if(nrow(proxy)==bigT){
      proxy <- proxy
    }else{
      stop("Please provide proxy of appropriate length!")
    }
  }
  #--------------------------------------------------------------------------------------------------------#
  # initialize objects to save IRFs, HDs, etc.
  R_store       <- array(NA, dim=c(thindraws,bigK,bigK))
  IRF_store     <- array(NA, dim=c(thindraws,n.ahead,bigK,bigK));dimnames(IRF_store)[[3]] <- varNames
  eps_store     <- array(NA, dim=c(thindraws,bigT,bigK));dimnames(eps_store)[[3]]<-varNames
  #------------------------------ start computing irfs  ---------------------------------------------------#
  start.comp <- Sys.time()
  imp.obj <- applyfun(1:thindraws,function(irep){
    Amat <- A_large[irep,,]
    Smat <- S_large[irep,,]
    Emat <- E_large[irep,,]
    imp.obj    <- irf.fun(xdat=xdat,plag=plag,n.ahead=n.ahead,Amat=Amat,Smat=Smat,shockinfo=shockinfo,MaxTries=MaxTries,type=type,Emat=Emat,proxy=proxy)
    if(verbose && ident=="sign"){
      if(!any(is.null(imp.obj$rot))){
        cat("\n",as.character(Sys.time()), "MCMC draw", irep, ": rotation found after ",imp.obj$icounter," tries", "\n")
      }else{
        cat("\n",as.character(Sys.time()), "MCMC draw", irep, ": no rotation found", "\n")
      }
    }
    return(list(impl=imp.obj$impl,rot=imp.obj$rot,eps=imp.obj$eps))
  })
  for(irep in 1:thindraws){
    if(all(is.na(imp.obj[[irep]]$impl))) next
    IRF_store[irep,,,] <- imp.obj[[irep]]$impl
    R_store[irep,,]    <- imp.obj[[irep]]$rot
    eps_store[irep,,]  <- imp.obj[[irep]]$eps
  }
  end.comp <- Sys.time()
  diff.comp <- difftime(end.comp,start.comp,units="mins")
  mins <- round(diff.comp,0); secs <- round((diff.comp-floor(diff.comp))*60,0)
  if(verbose) cat(paste("\nImpulse response analysis took ",mins," ",ifelse(mins==1,"min","mins")," ",secs, " ",ifelse(secs==1,"second.\n","seconds.\n"),sep=""))
  #------------------------------ post processing  ---------------------------------------------------#
  # re-set IRF object in case we have found only a few rotation matrices
  if(ident=="sign"){
    idx<-which(!is.na(apply(IRF_store,1,sum)))
    rot.nr<-paste("For ", length(idx), " draws out of ", thindraws, " draws, a rotation matrix has been found.")
    if(length(idx)==0){
      stop("No rotation matrix found with imposed sign restrictions. Please respecify.")
    }
    if(verbose) cat(rot.nr)
    # subset posterior draws
    IRF_store <- IRF_store[idx,,,,drop=FALSE]
    A_large   <- A_large[idx,,,drop=FALSE]
    S_large   <- S_large[idx,,,drop=FALSE]
    R_store   <- R_store[idx,,,drop=FALSE]
    eps_store <- eps_store[idx,,,drop=FALSE]
    thindraws <- length(idx)
  }
  # re-set IRF object in case we have found only a few rotation matrices
  if(ident=="proxy"){
    idx<-which(!is.na(apply(IRF_store,1,sum)))
    rot.nr<-paste("For ", length(idx), " draws out of ", thindraws, " draws, a rotation matrix has been found.")
    if(length(idx)==0){
      stop("No rotation matrix found with imposed sign restrictions. Please respecify.")
    }
    if(verbose) cat(rot.nr)
    # subset posterior draws
    IRF_store <- IRF_store[idx,,,,drop=FALSE]
    A_large   <- A_large[idx,,,drop=FALSE]
    S_large   <- S_large[idx,,,drop=FALSE]
    R_store   <- R_store[idx,,,drop=FALSE]
    eps_store <- eps_store[idx,,,drop=FALSE]
    thindraws <- length(idx)
  }
  # Subset to shocks under consideration
  IRF_store <- IRF_store[,,,select_shocks,drop=FALSE]
  imp_posterior <- array(NA, dim=c(n.ahead,bigK,shock.nr,bigQ),
                         dimnames=list(1:n.ahead, colnames(xglobal),paste("shock",colnames(xglobal)[select_shocks],sep="_"),
                                       paste("Q",str_pad(gsub("0\\.","",quantiles),width=2,side="right",pad="0"),sep=".")))
  # Normalization
  for(z in 1:shock.nr)
  {
    Mean<-IRF_store[,1,select_shocks[z],z]
    for(irep in 1:thindraws){
      IRF_store[irep,,,z]<-(IRF_store[irep,,,z]/Mean[irep])*scale[z]
    }
    for(qq in 1:bigQ){
      imp_posterior[,,z,qq] <- apply(IRF_store[,,,z],c(2,3),quantile,quantiles[qq],na.rm=TRUE)
    }
  }
  # calculate objects needed for HD and struc shock functions later---------------------------------------------
  # median quantitities
  Amat <- apply(A_large,c(2,3),median)
  Smat <- apply(S_large,c(2,3),median)
  Emat <- apply(E_large,c(2,3),median)
  if(ident=="sign"){
    imp.obj    <- try(irf(xdat=xdat,plag=plag,n.ahead=n.ahead,Amat=Amat,Smat=Smat,Emat=Emat,shockinfo=shockinfo,MaxTries=MaxTries),silent=TRUE)
    if(!is(imp.obj,"try-error")){
      Rmed<-imp.obj$rot
      epsmed<-imp.obj$eps
    }else{
      Rmed<-epsmed<-NULL
    }
  }
  struc.obj <- list(Amat=Amat,Smat=Smat,Emat=Emat,Rmed=Rmed,epsmed=epsmed)
  model.obj <- list(xglobal=xglobal,plag=plag,cons=x$args$cons,trend=x$args$trend)
  #--------------------------------- prepare output----------------------------------------------------------------------#
  out <- structure(list("posterior"   = imp_posterior,
                        "ident"       = ident,
                        "rot.nr"      = rot.nr,
                        "shockinfo"   = shockinfo,
                        "struc.obj"   = struc.obj,
                        "model.obj"   = model.obj),
                   class="bvar.irf")
  if(save.store){
    out$IRF_store = IRF_store
    out$eps_store = eps_store
    out$R_store   = R_store
  }
  return(out)
}

#' @name .irf.checks
#' @noRd
.irf.checks <- function(x=x, n.ahead=n.ahead, ident=ident, shockinfo=shockinfo, proxy=proxy, save.store=save.store, quantiles=quantiles){
  #-----------------------------------------------------------------------------------------------------#
  # check arguments
  if(!is.numeric(n.ahead) || n.ahead<1){
    stop("Please provide 'n.ahead' as numeric bigger than zero.")
  }
  if(is.null(ident)){
    stop("Please provide preferred identification scheme.")
  }
  if(!(ident%in%c("chol","chol-shortrun","chol-longrun","sign","proxy","girf"))){
    stop("Chosen identification scheme not available. Please respecify.")
  }
  if(!is.logical(save.store)){
    stop("Please provide argument 'save.store' as logical. Respecify.")
  }
  if(any(quantiles<0||quantiles>1)){
    stop("Please specify quantiles within 0 and 1.")
  }
  #-----------------------------------------------------------------------------------------------------#
  # check shockinfo argument for different identification schemes
  if(ident=="chol-shortrun"){
    # no checks
  }else if(ident=="chol-longrun"){
    # no checks
  }else if(ident=="girf"){
    # no checks
  }else if(ident=="sign"){
    if(!all(c("shock","restrictions","sign","horizon","scale")%in%colnames(shockinfo))){
      stop("Please provide columns 'shock', 'restrictions', 'sign', 'horizon' and 'scal' in dataframe 'shockinfo'.")
    }
    vars<-colnames(x$args$Data)
    if(!(all(shockinfo$shock%in%vars) && all(shockinfo$restrictions%in%vars))){
      stop("Please provide in columns 'shock' and 'restrictions' of 'shockinfo' only variable names available in dataset used for estimation.")
    }
    if(!any(shockinfo$sign%in%c(">","<","0","ratio.H","ratio.avg"))){
      stop("Misspecification in 'sign'. Only the following is allowed: <, >, 0, ratio.H, ratio.avg")
    }
  }else if(ident=="proxy"){
    if(!all(c("shock","instr","scale")%in%colnames(shockinfo)))
    if(!is.matrix(proxy)){
      stop("Please provide with argument 'shockinfo' matrix with instruments.")
    }
    if(nrow(proxy)!=nrow(x$xglobal) || nrow(proxy)!=nrow(x$xglobal-x$args$plag)){
      stop("Provide argument 'shockinfo' containing a matrix with same length as dataset used for estimation.")
    }
  }
}


#' @name get_shockinfo
#' @title Create \code{shockinfo} argument
#' @description Creates dummy \code{shockinfo} argument for appropriate use in  \code{irf} function.
#' @param ident Definition of identification scheme, either \code{chol} or \code{sign}.
#' @details Depending on the identification scheme a different \code{shockinfo} argument in the \code{irf} function is needed. To handle this convenient, an appropriate data.frame with is created with this function.
#' @usage get_shockinfo(ident="chol", nr_rows=1)
#' @seealso \code{\link{irf}}
#' @export
get_shockinfo <- function(ident="chol", nr_rows=1){
  if(ident=="chol")
    return(data.frame(shock=rep(NA,nr_rows),scale=rep(NA,nr_rows)))
  if(ident=="sign")
    return(data.frame(shock=rep(NA,nr_rows),restrictions=rep(NA,nr_rows),sign=rep(NA,nr_rows),
                      horizon=rep(NA,nr_rows),scale=rep(NA,nr_rows),prob=rep(NA,nr_rows),info=rep(NA,nr_rows)))
  if(ident=="proxy")
    return(data.frame(shock=rep(NA,nr_rows),instr=rep(NA,nr_rows),scale=rep(NA,nr_rows)))
}

#' @name add_shockinfo
#' @title Adding shocks to 'shockinfo' argument
#' @description Adds automatically rows to 'shockinfo' data.frame for appropriate use in \code{irf}.
#' @usage add_shockinfo(shockinfo=NULL, shock=NULL, restrictions=NULL, sign=NULL, horizon=NULL,
#' prob=NULL, scale=NULL, horizon.fillup=TRUE)
#' @param shockinfo Dataframe to append shocks. If \code{shockinfo=NULL} appropriate dataframe for sign-restrictions will be created.
#' @param shock String element. Variable of interest for structural shock. Only possible to add restrictions to one structural shock at a time.
#' @param restrictions Character vector with variables that are supposed to be sign restricted.
#' @param sign Character vector with signs.
#' @param horizon Numeric vector with horizons to which restriction should hold. Set \code{horizon.fillup} to \code{FALSE} to just restrict one specific horizon.
#' @param prob Number between zero and one determining the probability with which restriction is supposed to hold.
#' @param scale Scaling parameter.
#' @param horizon.fillup Default set to \code{TRUE}, horizon specified up to given horizon. Otherwise just one specific horizon is restricted.
#' @details This is only possible for sign restriction, hence if \code{ident="sign"} in \code{get_shockinfo()}.
#' @seealso \code{\link{irf}}
#' @export
add_shockinfo <- function(shockinfo=NULL, shock=NULL, restrictions=NULL, sign=NULL, horizon=NULL, prob=NULL, scale=NULL, horizon.fillup=TRUE){
  if(is.null(shockinfo)){
    shockinfo <- get_shockinfo(ident="sign")
  }
  if(is.null(shock)){
    stop("Please specify structural shock. This corresponds to the variable the shock is originating from.")
  }
  if(length(shock)>1){
    stop("Please only specify one structural shock at once.")
  }
  if(is.null(restrictions) || is.null(sign)){
    stop("Please specify 'restrictions' together with 'sign'.")
  }
  if(length(restrictions)!=length(sign) || length(restrictions)!=length(horizon) || length(sign)!=length(horizon)){
    if(length(horizon)!=1) stop("Please provide the arguments 'restrictions' and 'sign' with equal length. Please respecify.")
  }
  nr <- length(sign)
  if(!(is.null(restrictions) && is.null(sign)) && is.null(horizon)){
    warning("No horizon specified, is set to one, i.e., a shock restriction on impact.")
    horizon <- rep(1,nr)
  }
  if(!any(sign%in%c(">","<","0","ratio.H","ratio.avg"))){
    stop("Misspecification in 'sign'. Only the following is allowed: <, >, 0, ratio.H, ratio.avg")
  }
  if(is.null(scale)){
    warning("Scaling is not specified, set positive.")
    scale <- rep(1,nr)
  }
  if(length(scale)==1) scale <- rep(scale,nr)
  scale <- sign(scale)
  if(length(unique(scale))>1){
    warning("Different scaling supplied. Set to default value: positive.")
    scale <- rep(1,nr)
  }
  if(is.null(prob)){
    warning("Restriction proabilities not specified, set to one.")
    prob <- rep(1,nr)
  }
  if(length(prob)==1) prob <- rep(prob,nr)
  if(length(prob)!=nr || length(scale)!=nr){
    stop("Please specify 'prob' or 'scale' with unit length for all restrictions or equal length than restriction.")
  }
  if(length(horizon)==1 && length(horizon)<nr){
    warning("Only one horizon specified, is used for all horizons.")
    horizon <- rep(horizon,nr)
  }
  # if horizon is bigger than one
  idx_nr <- which(!(sign%in%c("ratio.H","ratio.avg")))
  idx_r  <- which(sign%in%c("ratio.H","ratio.avg"))
  if(any(horizon[idx_nr]>1) && horizon.fillup){
    repetition <- c(horizon[idx_nr],rep(1,length(idx_r)))
    restrictions <- rep(restrictions, repetition)
    sign <- rep(sign, repetition)
    prob <- rep(prob, repetition)
    scale <- rep(scale, repetition)
    horizon <- c(unlist(sapply(horizon[idx_nr],seq)),horizon[idx_r])
  }
  nr <- length(sign)
  # add to shockinfo
  nt<-ifelse(all(is.na(shockinfo)),0,nrow(shockinfo))
  for(nn in 1:nr){
    shockinfo[nt+nn,] <- NA
    shockinfo$shock[nt+nn] <- shock
    shockinfo$restrictions[nt+nn] <- restrictions[nn]
    shockinfo$sign[nt+nn] <- sign[nn]
    shockinfo$horizon[nt+nn] <- horizon[nn]
    shockinfo$prob[nt+nn] <- prob[nn]
    shockinfo$scale[nt+nn] <- scale[nn]
  }
  # delete duplicate lines
  shockinfo<-shockinfo[!duplicated(shockinfo),]
  return(shockinfo)
}


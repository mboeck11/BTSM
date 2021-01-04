#' @name irf
#' @title Impulse Response Function
#' @usage irf(x, n.ahead=24, ident=NULL, scal=1, sign.constr=NULL, proxy=NULL, save.store=FALSE,
#'            applyfun=NULL, cores=NULL, verbose=TRUE)
#' @param x object of class \code{bvar}.
#' @param n.ahead forecasting horizon.
#' @param ident preferred identification scheme.
#' @param scal scaling factor.
#' @param sign.constr the user should submit a list containing the following entries \itemize{
#' \item{\code{shock1}}{ is a list object that defines sign restrictions for a particular shock.}
#' \itemize{
#' \item{\code{shockvar}}{ is a character vector containing the variable to shock.}
#' \item{\code{restrictions}}{ is a list containing the variables to restrict.}
#' \item{\code{sign}}{ is a character vector of length of set of restrictions + 1, specifying the signs to impose. Use either \code{>}, \code{<} or \code{0}. The latter implements zero restrictions according to Arias et al. (2019). First entry is for the shock, say \code{AT.ltir} should go up, the following entries refer to the restrictions. \code{sign.constr$shock1$sign=c(">", "<", "<")} would impose \code{AT.ltir} to increase, and variables specified in \code{sign.constr$shock1$restricionts$rest1} and \code{sign.constr$shock1$restricionts$rest2} to decrease.}
#' \item{\code{rest.horz}}{ is a vector with same length as slot \code{sign} above and specifies the length of periods the restrictions are imposed. If \code{rest.horz} is 1, only impact restrictions are considered.}
#' \item{\code{constr}}{ is a vector with same length as slot \code{sign} above with elements lying between \code{0} and \code{1}. It specifies the percentage of countries for which cross-country restrictions have to hold. If no cross-country restrictions are supplied, set all elements of \code{constr} to 1.}
#' \item{\code{scal}}{ optional numeric in case impact normalization is desired.}
#' }
#' \item{\code{MaxTries}}{ Optional numeric corresponding to the maximum tries to search for a rotation matrix that fulfills the user-specified restrictions. Default is set to 7500. After \code{MaxTries} unsuccessful tries the algorithm sets the impulse response for that specific posterior draw to \code{NA}.}
#' \item{\code{shock2}}{ define a second list with the same arguments as \code{shock1} to identify a second shock. Can be used iteratively to identify multiple shocks.}
#' }
#' @param proxy in case of identification via proxy.
#' @param save.store If set to \code{TRUE} the full posterior is returned. Default is set to \code{FALSE} in order to save storage.
#' @param applyfun Allows for user-specific apply function, which has to have the same interface than \code{lapply}. If \code{cores=NULL} then \code{lapply} is used, if set to a numeric either \code{parallel::parLapply()} is used on Windows platforms and \code{parallel::mclapply()} on non-Windows platforms.
#' @param cores Specifies the number of cores which should be used. Default is set to \code{NULL} and \code{applyfun} is used.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @export
"irf" <- function(x, n.ahead=24, ident=NULL, scal=1, sign.constr=NULL, proxy=NULL, save.store=FALSE, applyfun=NULL, cores=NULL, verbose=TRUE,
                  quantiles=c(.05,.10,.16,.50,.84,.90,.95), ...){
  #------------------------------ message to console -------------------------------------------------------#
  if(verbose){
    if(class(x)=="bvar")
      cat("\nStart computing impulse response functions of Bayesian Vector Autoregression.\n\n")
    if(class(x)=="bvec")
      cat("\nStart computing impulse response functions of Bayesian Vector Error Correction Model.\n\n")
    if(class(x)=="bivar")
      cat("\nStart computing impulse response functions of Bayesian Interacted Vector Autoregression.\n\n")
    if(class(x)=="tvpbvar")
      cat("\nStart computing impulse response functions of Time-varying Parameter Bayesian Vector Autoregression.\n\n")
  }
  #------------------------------ do checks ---------------------------------------------------------------#
  .irf.checks(x=x, n.ahead=n.ahead, ident=ident, scal=scal, sign.constr=sign.constr, proxy=proxy)
  UseMethod("irf", x)
}

#' @export
irf.bvar <- function(x,n.ahead=24,ident=NULL,scal=NULL,sign.constr=NULL,proxy=NULL,save.store=FALSE,applyfun=NULL,cores=NULL,verbose=TRUE,
                     quantiles=c(.05,.10,.16,.50,.84,.90,.95), ...){
  start.irf <- Sys.time()
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
  out <- .irf.generator(x,n.ahead=n.ahead,ident=ident,scal=1,sign.constr=sign.constr,proxy=proxy,save.store=save.store,applyfun=applyfun,verbose=verbose)
  cat(paste("\nSize of irf object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @export
irf.tvpbvar <- function(x,n.ahead=24,ident=NULL,scal=NULL,sign.constr=NULL,proxy=NULL,save.store=FALSE,applyfun=NULL,cores=NULL,verbose=TRUE,
                        quantiles=c(.05,.10,.16,.50,.84,.90,.95), ...){
  start.irf <- Sys.time()
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
  #------------ get data -------------------------#
  Y    <- x$args$Y
  bigT <- nrow(Y)
  M    <- ncol(Y)
  thindraws <- x$args$thindraws
  #-------------------------------------------------------------------------------------------------------#
  if(verbose) cat(paste("Start impulse response analysis on ", cores, " cores", " (",thindraws," stable draws in total).",sep=""),"\n")
  # median response
  x.med <- x
  x.med$store$A_store <- apply(x.med$store$A_store,c(1,3,4),median)
  out <- .irf.generator(x.med,n.ahead=n.ahead,ident=ident,scal=1,sign.constr=sign.constr,proxy=proxy,save.store=save.store,applyfun=applyfun,verbose=FALSE)
  out$posterior.full <- array(NA,c(bigT,n.ahead,M,M,length(quantiles)),dimnames=c(list(NULL),dimnames(out$posterior)))
  for(tt in 1:bigT){
    x.t <- x
    x.t$store$A_store <- x.t$store$A_store[,tt,,]
    x.t$store$Smed_store <- x.t$store$S_store[,tt,,]
    out$posterior.full[tt,,,,] <- .irf.generator(x.t,n.ahead=n.ahead,ident=ident,scal=1,sign.constr=sign.constr,proxy=proxy,save.store=save.store,applyfun=applyfun,verbose=FALSE)$posterior
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
irf.bvec <- function(x,n.ahead=24,ident=NULL,scal=1,sign.constr=NULL,proxy=NULL,save.store=FALSE,applyfun=NULL,cores=NULL,verbose=TRUE,
                     quantiles=c(.05,.10,.16,.50,.84,.90,.95), ...){
  start.irf <- Sys.time()
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
  out <- .irf.generator(x,n.ahead=n.ahead,ident=ident,scal=1,sign.constr=sign.constr,proxy=proxy,save.store=save.store,applyfun=applyfun,verbose=verbose)
  cat(paste("\nSize of irf object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  if(verbose) cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @export
irf.bivar <- function(x,n.ahead=24,ident=NULL,scal=1,sign.constr=NULL,proxy=NULL,save.store=FALSE,applyfun=NULL,cores=NULL,verbose=TRUE,
                      eval.q=NULL){
  start.irf <- Sys.time()
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
  x$store$Smed_store <- apply(S_store,c(1,3,4),median)
  x$store$res_store  <- res_store

  out <- .irf.generator(x=x,n.ahead=n.ahead,ident=ident,scal=1,sign.constr=sign.constr,proxy=proxy,save.store=save.store,applyfun=applyfun,verbose=verbose)
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
.irf.generator <- function(x,n.ahead=24,ident=NULL,scal=NULL,sign.constr=NULL,proxy=NULL,save.store=FALSE,applyfun=NULL,
                     quantiles=c(.05,.10,.16,.50,.84,.90,.95), verbose=TRUE){
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
  #------------------------------ assign irf function  ---------------------------------------------------#
  if(ident=="sign"){
    # check MaxTries, if no MaxTries, set it to 7500
    MaxTries <- 7500
    if(!is.null(sign.constr$MaxTries)){
      MaxTries <- sign.constr$MaxTries
      sign.constr$MaxTries <- NULL
    }
    shock.nr<-length(sign.constr)
    irf<-.irf.sign.zero
    # first check whether all elements of sign restriction list have been specified
    res_len<-unlist(lapply(sign.constr, function(l){
      if(is.list(l$restrictions)){return(length(l$restrictions))}else{return(0)}
    }))
    varNames <- vars <- colnames(xglobal)
    shockvar <- unlist(lapply(sign.constr,function(l)l$shock))
    for(m in 1:shock.nr){
      aux<-sign.constr[[m]]
      if(!all(c("shock","restrictions","sign","rest.horz")%in%names(aux))){
        stop("Please specify for the object 'sign.contr' slots labeled 'shock', 'restrictions', 'sign', and 'rest.horz'. A slot 'scal' is optional. See the help files for more information.")
      }
      mm<-1+length(aux$restrictions) # one for the shock, rest for the restrictions
      if(length(aux$sign)!=mm || length(aux$rest.horz)!=mm){
        stop("Please specify M +1 signs, rest.horz and constr with M denoting the nr. of restrictions.")
      }
      if(!any(unlist(aux$restrictions)%in%varNames) && !is.null(aux$restrictions)){
        stop("Please restrict variables available in the dataset. Respecify.")
      }
      #-------------------------------------------------------
      sign.constr[[m]]$restrictions <- as.character(unlist(sign.constr[[m]]$restrictions))
    }
    names(sign.constr) <- paste("shock",seq(1,length(sign.constr)),sep="")
    for(m in 1:length(sign.constr)){
      if(any(sign.constr[[m]]$sign=="ratio.H"|sign.constr[[m]]$sign=="ratio.avg")){
        ratios <- which(sign.constr[[m]]$sign=="ratio.H"|sign.constr[[m]]$sign=="ratio.avg")
        temp.restr <- sign.constr[[m]]$restrictions[ratios-1]
        temp.signs <- sign.constr[[m]]$sign[ratios]
        temp.horiz <- sign.constr[[m]]$rest.horz[ratios]
        temp.list  <- list(restrictions=c(),
                           sign=c(),
                           rest.horz=c())
        for(ii in 1:length(ratios)){
          var_exp  <- temp.restr[ii]
          var_real <- strsplit(temp.restr[ii],"_")[[1]][1]
          if(temp.signs[ii]=="ratio.avg"){
            horizons <- seq(2,temp.horiz[ii]) # not on impact
          }else if(temp.signs[ii]=="ratio.H"){
            horizons <- temp.horiz[ii]
          }
          #### 0 restriction on expectations
          temp.list$restrictions <- c(temp.list$restrictions,var_exp)
          temp.list$sign         <- c(temp.list$sign,"0")
          temp.list$rest.horz    <- c(temp.list$rest.horz,1)
          #### -1 restriction indicates sum
          if(length(horizons)==1){
            temp.list$restrictions <- c(temp.list$restrictions,var_real)
            temp.list$sign         <- c(temp.list$sign,"-1")
            temp.list$rest.horz    <- c(temp.list$rest.horz,horizons)
          } else if(length(horizons)>1){ ### -1/H restriction indicates average
            temp.list$restrictions <- c(temp.list$restr,rep(var_real,length(horizons)))
            share <- as.character(-1/length(horizons))
            temp.list$sign         <- c(temp.list$sign,rep(share,length(horizons)))
            temp.list$rest.horz    <- c(temp.list$rest.horz,horizons)
          }
        }
        # kick old stuff
        sign.constr[[m]]$restrictions <- sign.constr[[m]]$restrictions[-c(ratios-1)]
        sign.constr[[m]]$sign <- sign.constr[[m]]$sign[-ratios]
        sign.constr[[m]]$rest.horz <- sign.constr[[m]]$rest.horz[-ratios]
        # get new stuff in
        sign.constr[[m]]$restrictions <- c(sign.constr[[m]]$restrictions,temp.list$restrictions)
        sign.constr[[m]]$sign <- c(sign.constr[[m]]$sign,temp.list$sign)
        sign.constr[[m]]$rest.horz <- c(sign.constr[[m]]$rest.horz,temp.list$rest.horz)
      }else if(any(sign.constr[[m]]$sign=="0")){
        zeros <- which(sign.constr[[m]]$sign=="0")
        temp.restr  <- sign.constr[[m]]$restrictions[zeros-1]
        temp.signs  <- sign.constr[[m]]$sign[zeros]
        temp.horiz  <- sign.constr[[m]]$rest.horz[zeros]
        temp.constr <- sign.constr[[m]]$constr[zeros]
        temp.list   <- list(restrictions=c(),
                            sign=c(),
                            rest.horz=c())
        for(ii in 1:length(zeros)){
          var_zero <- temp.restr[ii]
          horizons <- seq(1,temp.horiz[ii])
          #### 0 restriction on expectations
          temp.list$restrictions <- c(temp.list$restrictions,rep(var_zero,length(horizons)))
          temp.list$sign         <- c(temp.list$sign,rep("0",length(horizons)))
          temp.list$rest.horz    <- c(temp.list$rest.horz,horizons)
        }
        # kick old stuff
        sign.constr[[m]]$restrictions <- sign.constr[[m]]$restrictions[-c(zeros-1)]
        sign.constr[[m]]$sign <- sign.constr[[m]]$sign[-zeros]
        sign.constr[[m]]$rest.horz <- sign.constr[[m]]$rest.horz[-zeros]
        # get new stuff in
        sign.constr[[m]]$restrictions <- c(sign.constr[[m]]$restrictions,temp.list$restrictions)
        sign.constr[[m]]$sign <- c(sign.constr[[m]]$sign,temp.list$sign)
        sign.constr[[m]]$rest.horz <- c(sign.constr[[m]]$rest.horz,temp.list$rest.horz)
      }else if(any(sign.constr[[m]]$sign==">"|sign.constr[[m]]$sign=="<")){
        signs <- which(sign.constr[[m]]$sign==">"|sign.constr[[m]]$sign=="<")
        if(any(signs==1)){
          temp.restr <- c(sign.constr[[m]]$shock,sign.constr[[m]]$restrictions[signs-1])
          first <- TRUE
        }else{
          temp.restr <- sign.constr[[m]]$restrictions[signs-1]
        }
        temp.signs <- sign.constr[[m]]$sign[signs]
        temp.horiz <- sign.constr[[m]]$rest.horz[signs]
        temp.list  <- list(restrictions=c(),
                           sign=c(),
                           rest.horz=c())
        for(ii in 1:length(signs)){
          var_horz <- seq(1,temp.horiz[ii])
          temp.list$restrictions <- c(temp.list$restrictions,
                                      rep(temp.restr[ii],length(var_horz)))
          temp.list$sign         <- c(temp.list$sign,
                                      rep(temp.signs[ii],length(var_horz)))
          temp.list$rest.horz    <- c(temp.list$rest.horz,var_horz)
        }
        # kick old stuff
        sign.constr[[m]]$restrictions <- sign.constr[[m]]$restrictions[-c(signs-1)]
        sign.constr[[m]]$sign         <- sign.constr[[m]]$sign[-signs]
        sign.constr[[m]]$rest.horz    <- sign.constr[[m]]$rest.horz[-signs]
        # get new stuff in
        if(first){
          sign.constr[[m]]$restrictions <- c(sign.constr[[m]]$restrictions,
                                             temp.list$restrictions[-1])
        }else{
          sign.constr[[m]]$restrictions <- c(sign.constr[[m]]$restrictions,
                                             temp.list$restrictions)
        }
        sign.constr[[m]]$sign         <- c(sign.constr[[m]]$sign,
                                           temp.list$sign)
        sign.constr[[m]]$rest.horz    <- c(sign.constr[[m]]$rest.horz,
                                           temp.list$rest.horz)
      }
    }
    strg.list <- NULL
  }else if(ident=="chol-shortrun"){
    irf <- .irf.chol
    type <- "short-run"
    if(is.null(scal)) scal <- 1
    if(length(scal)==1) scal <- rep(scal,bigK)
    MaxTries<-str<-sign.constr<-rot.nr<-Rmed<-NULL
  }else if(ident=="girf"){
    irf <- .irf.girf
    if(is.null(scal)) scal <- 1
    if(length(scal)==1) scal <- rep(scal,bigK)
    MaxTries<-str<-sign.constr<-rot.nr<-Rmed<-NULL
  }else if(ident=="chol-longrun"){
    irf <- .irf.chol
    type <- "long-run"
    if(is.null(scal)) scal <- 1
    if(length(scal)==1) scal <- rep(scal,bigK)
    MaxTries<-str<-sign.constr<-rot.nr<-Rmed<-NULL
  }else if(ident=="proxy"){
    irf <- .irf.proxy
    if(is.null(scal)) scal <- 1
    if(length(scal)==1) scal <- rep(scal,bigK)
    MaxTries<-str<-sign.constr<-rot.nr<-Rmed<-type<-NULL
    if(nrow(proxy)==Traw) proxy <- proxy[(plag+1):Traw,,drop=FALSE]
  }
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
  #--------------------------------------------------------------------------------------------------------#
  # initialize objects to save IRFs, HDs, etc.
  R_store       <- array(NA, dim=c(thindraws,bigK,bigK))
  IRF_store     <- array(NA, dim=c(thindraws,n.ahead,bigK,bigK));dimnames(IRF_store)[[3]] <- varNames
  imp_posterior <- array(NA, dim=c(n.ahead,bigK,bigK,bigQ))
  dimnames(imp_posterior)[[1]] <- 1:n.ahead
  dimnames(imp_posterior)[[2]] <- colnames(xglobal)
  dimnames(imp_posterior)[[3]] <- paste("shock",colnames(xglobal),sep="_")
  dimnames(imp_posterior)[[4]] <- paste("Q",str_pad(gsub("0\\.","",quantiles),width=2,side="right",pad="0"),sep=".")
  #------------------------------ start computing irfs  ---------------------------------------------------#
  start.comp <- Sys.time()
  if(verbose) cat(paste("Start impulse response analysis on ", cores, " cores", " (",thindraws," stable draws in total).",sep=""),"\n")
  imp.obj <- applyfun(1:thindraws,function(irep){
    Amat <- A_large[irep,,]
    Smat <- S_large[irep,,]
    Emat <- E_large[irep,,]
    imp.obj    <- irf(xdat=xdat,plag=plag,n.ahead=n.ahead,Amat=Amat,Smat=Smat,sign.constr=sign.constr,
                      MaxTries=MaxTries,type=type,Emat=Emat,proxy=proxy)
    if(verbose){
      if(!is.null(sign.constr)){
        if(!any(is.null(imp.obj$rot))){
          cat("\n",as.character(Sys.time()), "MCMC draw", irep, ": rotation found after ",imp.obj$icounter," tries", "\n")
        }else{
          cat("\n",as.character(Sys.time()), "MCMC draw", irep, ": no rotation found", "\n")
        }
      }
    }
    return(list(impl=imp.obj$impl,rot=imp.obj$rot))
  })
  for(irep in 1:thindraws){
    if(all(is.na(imp.obj[[irep]]$impl))) next
    IRF_store[irep,,,] <- imp.obj[[irep]]$impl
    if(ident=="sign"){
      R_store[irep,,] <- imp.obj[[irep]]$rot
    }
  }
  end.comp <- Sys.time()
  diff.comp <- difftime(end.comp,start.comp,units="mins")
  mins <- round(diff.comp,0); secs <- round((diff.comp-floor(diff.comp))*60,0)
  if(verbose) cat(paste("\nImpulse response analysis took ",mins," ",ifelse(mins==1,"min","mins")," ",secs, " ",ifelse(secs==1,"second.","seconds.\n"),sep=""))
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
    thindraws <- length(idx)
  }
  # Normalization
  if(thindraws>0){
    for(z in 1:bigK){
      Mean<-IRF_store[,1,z,z]
      for(irep in 1:thindraws){
        IRF_store[irep,,,z]<-(IRF_store[irep,,,z]/Mean[irep])*scal[z]
      }
    }
    for(i in 1:bigK){
      for(q in 1:bigQ){
        imp_posterior[,,i,q]  <- apply(IRF_store[,,,i],c(2,3),quantile,quantiles[q],na.rm=TRUE)
      }
    }
  }
  # calculate objects needed for HD and struc shock functions later---------------------------------------------
  # median quantitities
  Amat    <- apply(A_large,c(2,3),median)
  Smat    <- apply(S_large,c(2,3),median)
  if(ident=="sign"){
    imp.obj    <- try(irf(xdat=xdat,plag=plag,n.ahead=n.ahead,Amat=Amat,Smat=Smat,
                          sign.constr=sign.constr,MaxTries=MaxTries,shock.nr=shock.nr),silent=TRUE)
    if(!is(imp.obj,"try-error")){
      Rmed<-imp.obj$rot
    }else{
      Rmed<-NULL
    }
  }
  struc.obj <- list(Amat=Amat,Smat=Smat,Rmed=Rmed)
  model.obj <- list(xglobal=xglobal,plag=plag,cons=x$args$cons,trend=x$args$trend)
  #--------------------------------- prepare output----------------------------------------------------------------------#
  out <- structure(list("posterior"   = imp_posterior,
                        "ident"       = ident,
                        "rot.nr"      = rot.nr,
                        "sign.constr" = sign.constr,
                        "struc.obj"   = struc.obj,
                        "model.obj"   = model.obj),
                   class="bvar.irf")
  if(save.store){
    out$IRF_store = IRF_store
  }
  return(out)
}

#' @name .irf.checks
#' @noRd
.irf.checks <- function(x, n.ahead, ident, scal, sign.constr, proxy){
  # checks general
  if(is.null(ident)){
    stop("Please provide preferred identification scheme.")
  }
  # checks identification via sign restrictions
  if(!is.null(sign.constr)){
    if(ident!="sign"){
      stop("Please select 'sign' as identification scheme when providing a list of sign restrictions.")
    }
    # check whether sign restriction list is correctly specified, for each shock have to specify
    tt<-all(sapply(lapply(sign.constr,names),function(x) all(c("shock","sign","restrictions")%in%x)))
    if(!tt){
      stop("For each shock (i.e., first layer in the list), please specify lists named shock, sign and restrictions. See the details and examples in the manual.")
    }
    # check scaling, if no scaling, set it to 1
    if(is.null(scal)) scal <- 1
    if(length(scal)!=bigK) scal <- rep(scal,bigK)
    type <- sign.constr$type
    if(is.null(type)) type <- "short-run"
    # check signs horizons
    hh<-unlist(lapply(sign.constr,function(l)l$rest.horz))
    if(any(hh==0)){
      stop("Please provide contemporanous signs as horizon==1.")
    }
    if(any(hh>n.ahead)){
      stop("Do not restrict anything after the impulse response horizon given by 'n.ahead'.")
    }
  }
  # checks for proxy
  if(ident=="proxy"){
    if(nrow(proxy)!=nrow(x$xglobal)){
      stop("Provide 'proxy' with same length as dataset.")
    }
  }
}


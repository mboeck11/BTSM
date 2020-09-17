#' @name irf
#' @title Impulse Response Function
#' @usage irf(x, n.ahead=24, shock=NULL, sign.constr=NULL, save.store=FALSE,
#'            applyfun=NULL, cores=NULL, verbose=TRUE)
#' @param x object of class \code{bvar}.
#' @param n.ahead forecasting horizon.
#' @param shock This is a list object. It should contain an entry labeled \code{var} that contains the name of the variable to be shocked. Finally it has to contain an entry labeled \code{ident} that is either \code{chol} if the shock is based on a short-run identification scheme done with the Cholesky decomposition or \code{girf} if generalized impulse responses should be calculated. In case impulses should be normalized (e.g., a +100bp increase in the interest rate), add another entry \code{scal} that contains a numeric value of the desired impact normalization.
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
#' @param save.store If set to \code{TRUE} the full posterior is returned. Default is set to \code{FALSE} in order to save storage.
#' @param applyfun Allows for user-specific apply function, which has to have the same interface than \code{lapply}. If \code{cores=NULL} then \code{lapply} is used, if set to a numeric either \code{parallel::parLapply()} is used on Windows platforms and \code{parallel::mclapply()} on non-Windows platforms.
#' @param cores Specifies the number of cores which should be used. Default is set to \code{NULL} and \code{applyfun} is used.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @export
"irf" <- function(x, n.ahead=24, shock=NULL, sign.constr=NULL, save.store=FALSE, applyfun=NULL, cores=NULL, verbose=TRUE){
  UseMethod("irf", x)
}

#' @export
#' @importFrom parallel parLapply mclapply
#' @importFrom stats median
#' @importFrom utils object.size
irf.bvar <- function(x,n.ahead=24,shock=NULL,sign.constr=NULL,save.store=FALSE,applyfun=NULL,cores=NULL,verbose=TRUE){
  start.irf <- Sys.time()
  if(verbose) cat("\nStart computing impulse response functions of Bayesian Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  plag        <- x$args$plag
  xglobal     <- x$args$Yraw
  Traw        <- nrow(xglobal)
  bigK        <- ncol(xglobal)
  bigT        <- Traw-plag
  A_large     <- x$store$A
  S_large     <- x$store$SIGMA
  xdat        <- xglobal[(plag+1):Traw,,drop=FALSE]
  thindraws   <- x$args$draws
  varNames    <- colnames(xglobal)
  #------------------------------ user checks  ---------------------------------------------------#
  # checks general
  if(is.null(shock)&&is.null(sign.constr)){
    stop("Please provide either the list objects shock OR sign.constr to run GIRF / Cholesky or sign restriction based impulse response analysis.")
  }
  if(!is.null(shock)&&!is.null(sign.constr)){
    stop("Please provide either the list objects shock OR sign.constr to run GIRF / Cholesky or sign restriction based impulse response analysis.")
  }
  # checks cholesky/girf
  if(!is.null(shock)){
    if(!all(c("var","ident")%in%names(shock))){
      stop("Please provide list entries shock, ccode and ident for the object shock. See the help files for more details.")
    }
    if(!shock$ident%in%c("girf","chol")){
      stop("Please provide identification scheme.")
    }
    shockvar <- shock$var
    ident    <- shock$ident
    scal     <- shock$scal
    type     <- shock$type
    if(is.null(type)) type <- "short-run"
    if(is.null(scal)) scal <- rep(1,bigK)
    if(length(scal)!=bigK) scal <- rep(scal,bigK)
    Rmed <- NULL
    shock.nr <- length(shock$var)
  }
  # checks identification via sign restrictions
  if(!is.null(sign.constr)){
    ident<-"sign"
    # check MaxTries, if no MaxTries, set it to 7500
    MaxTries <- 7500
    shock.nr <- length(sign.constr)
    if(!is.null(sign.constr$MaxTries)){
      MaxTries <- sign.constr$MaxTries
      sign.constr$MaxTries <- NULL
    }
    shock.nr<-length(sign.constr)
    # check whether sign restriction list is correctly specified, for each shock have to specify
    tt<-all(sapply(lapply(sign.constr,names),function(x) all(c("shock","sign","restrictions")%in%x)))
    if(!tt){
      stop("For each shock (i.e., first layer in the list), please specify lists named shock, sign and restrictions. See the details and examples in the manual.")
    }
    # check scaling, if no scaling, set it to 1
    for(kk in 1:shock.nr){
      sign.constr[[kk]]$scal<-ifelse(is.null(sign.constr[[kk]]$scal),1,sign.constr[[kk]]$scal)
    }
    scal <- unlist(sapply(sign.constr,function(x) x$scal))
    type <- sign.constr$type
    if(is.null(type)) type <- "short-run"
  }
  #------------------------------ assign irf function  ---------------------------------------------------#
  if(ident=="sign"){
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
      if(!any(unlist(aux$restrictions)%in%varNames)){
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
  }else if(ident=="chol"){
    if(verbose) {
      cat("Identification scheme: Short-run identification via Cholesky decomposition.\n")
      cat(paste("Structural shock of interest: ", paste(shock$var,collapse=", "),".\n",sep=""))
    }
    irf <- .irf.chol
    MaxTries<-str<-sign.constr<-Srots<-rot.nr<-NULL
  }else{
    if(verbose){
      cat("Identification scheme: Generalized impulse responses.\n")
      cat(paste("Structural shock of interest: ", paste(shock$var,collapse=", "),".\n",sep=""))
    }
    irf <- .irf.girf
    MaxTries<-str<-sign.constr<-rot.nr<-NULL
  }

  # initialize objects to save IRFs, HDs, etc.
  R_store       <- array(NA, dim=c(thindraws,bigK,bigK))
  IRF_store     <- array(NA, dim=c(thindraws,n.ahead,bigK,shock.nr));dimnames(IRF_store)[[3]] <- varNames
  imp_posterior <- array(NA, dim=c(n.ahead,bigK,shock.nr,7))
  dimnames(imp_posterior)[[1]] <- 1:n.ahead
  dimnames(imp_posterior)[[2]] <- colnames(xglobal)
  dimnames(imp_posterior)[[3]] <- paste("shock",shockvar,sep="_")
  dimnames(imp_posterior)[[4]] <- c("low05","low10","low16","median","high84","high90","high95")
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
  #------------------------------ start computing irfs  ---------------------------------------------------#
  start.comp <- Sys.time()
  if(verbose) cat(paste("Start impulse response analysis on ", cores, " cores", " (",thindraws," stable draws in total).",sep=""),"\n")
  imp.obj <- applyfun(1:thindraws,function(irep){
    Amat <- A_large[irep,,]
    Smat <- S_large[irep,,]
    imp.obj    <- irf(xdat=xdat,plag=plag,n.ahead=n.ahead,Amat=Amat,Smat=Smat,sign.constr=sign.constr,
                      MaxTries=MaxTries,shock.nr=shock.nr,type=type)
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
    IRF_store[irep,,,] <- imp.obj[[irep]]$impl[,,shockvar,drop=FALSE]
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
    if(!is.null(scal)){
      for(z in 1:shock.nr){
        Mean<-IRF_store[,1,z,z]
        for(irep in 1:thindraws){
          IRF_store[irep,,,z]<-(IRF_store[irep,,,z]/Mean[irep])*scal[z]
        }
      }
    }

    for(i in 1:shock.nr){
      imp_posterior[,,i,"low16"]  <- apply(IRF_store[,,,i],c(2,3),quantile,0.16,na.rm=TRUE)
      imp_posterior[,,i,"low10"]  <- apply(IRF_store[,,,i],c(2,3),quantile,0.10,na.rm=TRUE)
      imp_posterior[,,i,"low05"]  <- apply(IRF_store[,,,i],c(2,3),quantile,0.05,na.rm=TRUE)
      imp_posterior[,,i,"median"] <- apply(IRF_store[,,,i],c(2,3),median,na.rm=TRUE)
      imp_posterior[,,i,"high84"] <- apply(IRF_store[,,,i],c(2,3),quantile,0.84,na.rm=TRUE)
      imp_posterior[,,i,"high90"] <- apply(IRF_store[,,,i],c(2,3),quantile,0.90,na.rm=TRUE)
      imp_posterior[,,i,"high95"] <- apply(IRF_store[,,,i],c(2,3),quantile,0.95,na.rm=TRUE)
    }
  }
  # calculate objects needed for HD and struc shock functions later---------------------------------------------
  # median quantitities
  Amat    <- apply(A_large,c(2,3),median)
  Smat    <- apply(S_large,c(2,3),median)
  if(ident=="sign"){
    imp.obj    <- try(irf(xdat=xdat,plag=plag,n.ahead=n.ahead,Amat=Amat,Smat=Smat,shock=shock,
                          sign.constr=sign.constr,MaxTries=MaxTries,shock.nr=shock.nr),silent=TRUE)
    if(!is(imp.obj,"try-error")){
      Rmed<-imp.obj$rot
    }else{
      Rmed<-NULL
    }
  }
  struc.obj <- list(Amat=Amat,Smat=Smat,Rmed=Rmed)
  model.obj <- list(xglobal=xglobal,plag=plag)
  #--------------------------------- prepare output----------------------------------------------------------------------#
  out <- structure(list("posterior"   = imp_posterior,
                        "rot.nr"      = rot.nr,
                        "shock"       = shock,
                        "sign.constr" = sign.constr,
                        "ident"       = ident,
                        "struc.obj"   = struc.obj,
                        "model.obj"   = model.obj),
                   class="bvar.irf")
  if(save.store){
    out$IRF_store = IRF_store
  }
  if(verbose) cat(paste("\nSize of irf object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  if(verbose) cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

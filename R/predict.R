#' @name predict
#' @title Predictions
#' @param object an object of class \code{bvar}.
#' @param ... additional arguments.
#' @param n.ahead the forecast horizon.
#' @param quantiles posterior quantiles to be computed.
#' @param applyfun parallelization
#' @param cores number of cores
#' @param save.store If set to \code{TRUE} the full distribution is returned. Default is set to \code{FALSE} in order to save storage.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @importFrom stats rnorm tsp sd
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
"predict" <- function(object, ..., n.ahead=4, quantiles=c(.05,.10,.16,.50,.84,.90,.95), applyfun=NULL, cores=NULL, save.store=FALSE, verbose=TRUE){
  #------------------------------ message to console -------------------------------------------------------#
  if(verbose){
    if(class(object)=="bvar")
      cat("\nStart doing predictions of Bayesian Vector Autoregression.\n\n")
    if(class(object)=="bvec")
      cat("\nStart doing predictions of Bayesian Vector Error Correction Model.\n\n")
    if(class(object)=="bivar")
      cat("\nStart doing predictions of Bayesian Interacted Vector Autoregression.\n\n")
    if(class(object)=="tvpbvar")
      cat("\nStart doing predictions of Time-varying Parameter Bayesian Vector Autoregression.\n\n")
  }
  UseMethod("predict", object)
}

#' @export
predict.bvar <- function(object, ..., n.ahead=4, quantiles=c(.05,.10,.16,.50,.84,.90,.95), applyfun=NULL, cores=NULL, save.store=FALSE, verbose=TRUE){
  start.pred <- Sys.time()
  if(verbose) cat("\nStart computing predictions of Bayesian Vector Autoregression.\n\n")
  if(verbose) cat("Start computing...\n")
  out <- predict.generator(object=object, n.ahead=n.ahead, quantiles=quantiles, applyfun=applyfun, cores=cores, save.store=save.store, TVP=FALSE)
  if(verbose) cat(paste("\n\nSize of object:", format(object.size(out),unit="MB")))
  end.pred <- Sys.time()
  diff.pred <- difftime(end.pred,start.pred,units="mins")
  mins.pred <- round(diff.pred,0); secs.pred <- round((diff.pred-floor(diff.pred))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.pred," ",ifelse(mins.pred==1,"min","mins")," ",secs.pred, " ",ifelse(secs.pred==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @export
predict.tvpbvar <- function(object, ..., n.ahead=4, quantiles=c(.05,.10,.16,.50,.84,.90,.95), applyfun=NULL, cores=NULL, save.store=FALSE, verbose=TRUE){
  start.pred <- Sys.time()
  if(verbose) cat("\nStart computing predictions of Time-varying parameter Bayesian Vector Autoregression.\n\n")
  if(verbose) cat("Start computing...\n")
  out <- predict.generator(object=object, n.ahead=n.ahead, quantiles=quantiles, applyfun=applyfun, cores=cores, save.store=save.store, TVP=TRUE)
  if(verbose) cat(paste("\n\nSize of object:", format(object.size(out),unit="MB")))
  end.pred <- Sys.time()
  diff.pred <- difftime(end.pred,start.pred,units="mins")
  mins.pred <- round(diff.pred,0); secs.pred <- round((diff.pred-floor(diff.pred))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.pred," ",ifelse(mins.pred==1,"min","mins")," ",secs.pred, " ",ifelse(secs.pred==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @name .predict.generator
#' @importFrom stringr str_pad
#' @noRd
predict.generator <- function(object, n.ahead, quantiles=c(.05,.10,.16,.50,.84,.90,.95), applyfun=NULL, cores=NULL, save.store=FALSE, TVP=FALSE){
  thindraws  <- object$args$thindraws
  plag       <- object$args$plag
  prior      <- object$args$prior
  xglobal    <- object$xglobal
  x          <- xglobal[(plag+1):nrow(xglobal),]
  S_store    <- object$store$S_store
  A_store    <- object$store$A_store
  pars_store <- object$store$pars_store
  vola_store <- object$store$vola_store
  L_store    <- object$store$L_store
  thetasqrt_store<-object$store$thetasqrt_store
  Lthetasqrt_store<-object$store$Lthetasqrt_store
  Omega_store<-object$store$Omega_store
  LDOmega_store<-object$store$LOmega_store
  varNames   <- colnames(xglobal)
  Traw       <- nrow(xglobal)
  bigT       <- nrow(x)
  M          <- ncol(xglobal)
  h          <- object$args$h
  cons       <- ifelse(object$args$cons,1,0)
  trend      <- ifelse(object$args$trend,1,0)
  K          <- M*plag+cons+trend
  SV         <- object$args$SV
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
  #------------------------------ looping --------------------------------------------------------#
  # start loop here
  pred_store <- array(NA, dim=c(thindraws,M,n.ahead), dimnames=list(NULL,varNames,1:n.ahead))
  pred.obj <- applyfun(1:thindraws,function(irep){
    #Step I: get coefficients and varianes
    if(TVP){
      A_t <- A_store[irep,bigT,,]
      if(prior=="TVP" || prior == "TVP-NG") Q_t <- as.vector(thetasqrt_store[irep,,]^2)
      if(prior=="TTVP") Q_t <- as.vector(Omega_store[irep,bigT,,])
      L_t <- L_store[irep,bigT,,]
      if(prior=="TVP" || prior == "TVP-NG") LQ_t <- lapply(Lthetasqrt_store,function(l)l[irep,]^2)
      if(prior=="TTVP") LQ_t <- lapply(LDOmega_store,function(l)l[irep,bigT,])
    } else {
      A_t <- A_store[irep,,]
      L_t <- L_store[irep,,]
    }
    Htt <- vola_store[irep,bigT,]
    Stt <- S_store[irep,bigT,,]
    if(SV){
      pars_var <- pars_store[irep,,]
    }
    # Step II: get last data point
    Xpred <- .mlag_pred(xglobal,plag)
    Mean00  <- t(Xpred[bigT,,drop=FALSE])
    Sigma00 <- matrix(0,M*plag,M*plag)
    y2      <- NULL
    #gets companion form
    aux   <- .gen_compMat(A_t,M,plag)
    Mm    <- aux$Cm
    Jm    <- aux$Jm
    if(!is.null(Htt)) Sig_t <- L_t %*% diag(exp(Htt)) %*% t(L_t) else Sig_t <- Stt
    Jsigt <- Jm%*%Sig_t%*%t(Jm)
    if(cons==1) consf <- rbind(t(A_t["cons",,drop=FALSE]),matrix(0,(plag-1)*M,1)) else consf <- matrix(0,M*plag,1)
    if(trend==1) trendf <- rbind(t(A_t["trend",,drop=FALSE]),matrix(0,(plag-1)*M,1)) else trendf <- matrix(0,M*plag,1)
    # this is the forecast loop
    for (ih in 1:n.ahead){
      Mean00  <- Mm%*%Mean00 + consf + trendf*(bigT+ih)
      Sigma00 <- Mm%*%Sigma00%*%t(Mm) + Jsigt
      chol_varyt <- try(t(chol(Sigma00[1:M,1:M])),silent=TRUE)
      if(is(chol_varyt,"try-error")){
        yf <- mvrnorm(1,mu=Mean00[1:M],Sigma00[1:M,1:M])
      }else{
        yf <- Mean00[1:M]+chol_varyt%*%rnorm(M,0,1)
      }
      # if(TVP){
      #   A_t <- matrix(rnorm(K*M,as.vector(A_t),Q_t),K,M, dimnames=dimnames(A_t))
      #   Q_t <- Q_t + Q_t
      #   for(mm in 2:M){
      #     L_t[mm,1:(mm-1)] <- rnorm(mm-1,L_t[mm,1:(mm-1)],LQ_t[[mm-1]])
      #     LQ_t[[mm-1]] <- LQ_t[[mm-1]] + LQ_t[[mm-1]]
      #   }
      # }
      if(SV){
        Htt   <- pars_var["mu",]+pars_var["phi",]*(Htt-pars_var["mu",])+rnorm(M,0,sqrt(pars_var["sigma",]))
        Jsigt <- L_t %*% diag(exp(Htt))%*%t(L_t)
      }
      y2 <- cbind(y2,yf)
    }
    return(y2)
  })
  for(irep in 1:thindraws){
    pred_store[irep,,] <- pred.obj[[irep]]
  }
  pred_post <- array(NA, dim=c(length(quantiles),M,n.ahead),
                     dimnames=list(paste("Q",str_pad(gsub("0\\.","",quantiles),width=2,side="right",pad="0"),sep="."),varNames,1:n.ahead))
  for(qq in 1:length(quantiles)){
    pred_post[qq,,] <- apply(pred_store,c(2,3),quantile,quantiles[qq])
  }

  if(h>n.ahead) {
    h <- n.ahead
    warning("Argument 'h' bigger than 'n.ahead'. Evaluation of forecasts only 'n.ahead' periods, hence 'h' is set to 'n.ahead'.")
  }
  yfull <- object$args$yfull
  if(h>0){
    lps.stats                <- array(0,dim=c(M,2,h))
    dimnames(lps.stats)[[1]] <- colnames(xglobal)
    dimnames(lps.stats)[[2]] <- c("mean","sd")
    dimnames(lps.stats)[[3]] <- 1:h
    lps.stats[,"mean",]      <- apply(pred_store[,,1:h,drop=FALSE],c(2,3),mean)
    lps.stats[,"sd",]        <- apply(pred_store[,,1:h,drop=FALSE],c(2,3),sd)
    hold.out<-yfull[(nrow(yfull)+1-h):nrow(yfull),,drop=FALSE]
  }else{
    lps.stats<-NULL
    hold.out<-NULL
  }
  out <- structure(list(fcast=pred_post,
                        xglobal=xglobal,
                        n.ahead=n.ahead,
                        lps.stats=lps.stats,
                        hold.out=hold.out),
                   class="bvar.pred")
  if(save.store){
    out$pred_store = pred_store
  }
  return(out)
}

#' @name cond.predict
#' @title Conditional Forecasts
#' @usage cond.predict(constr, bvar.obj, pred.obj, constr_sd=NULL, verbose=TRUE)
#' @details Conditional forecasts need a fully identified system. Therefore this function utilizes short-run restrictions via the Cholesky decomposition on the global solution of the variance-covariance matrix of the Bayesian VAR.
#' @param constr a matrix containing the conditional forecasts of size horizon times K, where horizon corresponds to the forecast horizon specified in \code{pred.obj}, while K is the number of variables in the system. The ordering of the variables have to correspond the ordering of the variables in the system. Rest is just set to NA.
#' @param bvar.obj an item fitted by \code{bvar}.
#' @param pred.obj an item fitted by \code{predict}. Note that \code{save.store=TRUE} is required as argument!
#' @param constr_sd a matrix containing the standard deviations around the conditional forecasts. Must have the same size as \code{constr}.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @importFrom abind adrop
#' @importFrom stats rnorm
#' @export
cond.predict <- function(constr, bvar.obj, pred.obj, constr_sd=NULL, verbose=TRUE){
  start.cond <- Sys.time()
  if(verbose) cat("\nStart conditional forecasts of Bayesian Global Vector Autoregression.\n\n")
  #----------------get stuff-------------------------------------------------------#
  plag        <- bvar.obj$args$plag
  xglobal     <- pred.obj$xglobal
  Traw        <- nrow(xglobal)
  bigK        <- ncol(xglobal)
  bigT        <- Traw-plag
  A_large     <- bvar.obj$stacked.results$A_large
  F_large     <- bvar.obj$stacked.results$F_large
  S_large     <- bvar.obj$stacked.results$S_large
  Ginv_large  <- bvar.obj$stacked.results$Ginv_large
  F.eigen     <- bvar.obj$stacked.results$F.eigen
  thindraws   <- length(F.eigen)
  x           <- xglobal[(plag+1):Traw,,drop=FALSE]
  horizon     <- pred.obj$n.ahead
  varNames    <- colnames(xglobal)
  cN          <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[1]))
  var         <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[2]))
  #---------------------checks------------------------------------------------------#
  if(is.null(pred.obj$pred_store)){
    stop("Please set 'save.store=TRUE' when computing predictions.")
  }
  if(!all(dim(constr)==c(horizon,bigK))){
    stop("Please respecify dimensions of 'constr'.")
  }
  if(!is.null(constr_sd)){
    if(!all(dim(constr_sd)==c(horizon,bigK))){
      stop("Please respecify dimensions of 'constr_sd'.")
    }
    constr_sd[is.na(constr_sd)] <- 0
  }else{
    constr_sd <- matrix(0,horizon,bigK)
  }
  pred_array <- pred.obj$pred_store
  #---------------container---------------------------------------------------------#
  cond_pred <- array(NA, c(thindraws, bigK, horizon))
  dimnames(cond_pred)[[2]] <- varNames
  #----------do conditional forecasting -------------------------------------------#
  if(verbose) cat("Start computing...\n")
  if(verbose) pb <- txtProgressBar(min = 0, max = thindraws, style = 3)
  for(irep in 1:thindraws){
    pred    <- pred_array[irep,,]
    Sigma_u <- Ginv_large[irep,,]%*%S_large[irep,,]%*%t(Ginv_large[irep,,])
    irf     <- .impulsdtrf(B=adrop(F_large[irep,,,,drop=FALSE],drop=1),
                           smat=t(chol(Sigma_u)),nstep=horizon)

    temp <- as.vector(constr) + rnorm(bigK*horizon,0,as.vector(constr_sd))
    constr_use <- matrix(temp,horizon,bigK)

    v <- sum(!is.na(constr))
    s <- bigK * horizon
    r <- c(rep(0, v))
    R <- matrix(0, v, s)
    pos <- 1
    for(i in 1:horizon) {
      for(j in 1:bigK) {
        if(is.na(constr_use[i, j])) {next}
        r[pos] <- constr_use[i, j] - pred[j, i]
        for(k in 1:i) {
          R[pos, ((k - 1) * bigK + 1):(k * bigK)] <- irf[j,,(i - k + 1)]
        }
        pos <- pos + 1
      }
    }

    R_svd <- svd(R, nu=nrow(R), nv=ncol(R))
    U     <- R_svd[["u"]]
    P_inv <- diag(1/R_svd[["d"]])
    V1    <- R_svd[["v"]][,1:v]
    V2    <- R_svd[["v"]][,(v+1):s]
    eta   <- V1 %*% P_inv %*% t(U) %*% r + V2 %*% rnorm(s-v)
    eta   <- matrix(eta, horizon, bigK, byrow=TRUE)

    for(h in 1:horizon) {
      temp <- matrix(0, bigK, 1)
      for(k in 1:h) {
        temp <- temp + irf[, , (h - k + 1)] %*% t(eta[k , , drop=FALSE])
      }
      cond_pred[irep,,h] <- pred[,h,drop=FALSE] + temp
    }
    if(verbose) setTxtProgressBar(pb, irep)
  }
  #------------compute posteriors----------------------------------------------#
  imp_posterior<-array(NA,dim=c(bigK,horizon,5))
  dimnames(imp_posterior)[[1]] <- varNames
  dimnames(imp_posterior)[[2]] <- 1:horizon
  dimnames(imp_posterior)[[3]] <- c("low25","low16","median","high75","high84")

  imp_posterior[,,"low25"]  <- apply(cond_pred,c(2,3),quantile,0.25,na.rm=TRUE)
  imp_posterior[,,"low16"]  <- apply(cond_pred,c(2,3),quantile,0.16,na.rm=TRUE)
  imp_posterior[,,"median"] <- apply(cond_pred,c(2,3),quantile,0.50,na.rm=TRUE)
  imp_posterior[,,"high75"] <- apply(cond_pred,c(2,3),quantile,0.75,na.rm=TRUE)
  imp_posterior[,,"high84"] <- apply(cond_pred,c(2,3),quantile,0.84,na.rm=TRUE)

  #----------------------------------------------------------------------------------#
  out <- structure(list(fcast=imp_posterior,
                        xglobal=xglobal,
                        n.ahead=horizon),
                   class="bgvar.pred")
  if(verbose) cat(paste("\n\nSize of object:", format(object.size(out),unit="MB")))
  end.cond <- Sys.time()
  diff.cond <- difftime(end.cond,start.cond,units="mins")
  mins.cond <- round(diff.cond,0); secs.cond <- round((diff.cond-floor(diff.cond))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.cond," ",ifelse(mins.cond==1,"min","mins")," ",secs.cond, " ",ifelse(secs.cond==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @export
"lps" <- function(object){
  UseMethod("lps", object)
}

#' @name lps
#' @title Compute Log-predictive Scores
#' @method lps bvar.pred
#' @description  Computes and prints log-predictive score of an object of class \code{bvar.predict}.
#' @param object an object of class \code{bvar.predict}.
#' @param ... additional arguments.
#' @return Returns an object of class \code{bvar.lps}, which is a matrix of dimension h times K, whereas h is the forecasting horizon and K is the number of variables in the system.
#' @author Maximilian Boeck
#' @importFrom stats dnorm
#' @export
lps.bvar.pred <- function(object, ...){
  hold.out <- object$hold.out
  h        <- nrow(hold.out)
  K        <- ncol(hold.out)
  if(is.null(hold.out)){
    stop("Please submit a forecast object that includes a hold out sample for evaluation (set h>0 when estimating the model with bgvar)!")
  }
  lps.stats  <- object$lps.stats
  lps.scores <- matrix(NA,h,K)
  for(i in 1:K){
    lps.scores[,i]<-dnorm(hold.out[,i],mean=lps.stats[i,"mean",],sd=lps.stats[i,"sd",],log=TRUE)
  }
  colnames(lps.scores)<-dimnames(lps.stats)[[1]]
  out <- structure(lps.scores, class="bgvar.lps")
  return(out)
}

#' @export
"rmse" <- function(object){
  UseMethod("rmse", object)
}

#' @name rmse
#' @title Compute Root Mean Squared Errors
#' @method rmse bvar.pred
#' @description  Computes and prints root mean squared errors (RMSEs) of an object of class \code{bvar.predict}.
#' @param object an object of class \code{bvar.predict}.
#' @param ... additional arguments.
#' @return Returns an object of class \code{bvar.rmse}, which is a matrix of dimension h times K, whereas h is the forecasting horizon and K is the number of variables in the system.
#' @author Maximilian Boeck
#' @importFrom knitr kable
#' @importFrom stats dnorm
#' @export
rmse.bvar.pred <- function(object, ...){
  hold.out <- object$hold.out
  h        <- nrow(hold.out)
  K        <- ncol(hold.out)
  if(is.null(hold.out)){
    stop("Please submit a forecast object that includes a hold out sample for evaluation (set h>0 in fcast)!")
  }
  lps.stats   <- object$lps.stats
  rmse.scores <- matrix(NA,h,K)
  for(i in 1:K){
    rmse.scores[,i]<-sqrt((hold.out[,i]-lps.stats[i,"mean",])^2)
  }
  colnames(rmse.scores)<-dimnames(lps.stats)[[1]]
  out <- structure(rmse.scores, class="bgvar.rmse")
  return(out)
}

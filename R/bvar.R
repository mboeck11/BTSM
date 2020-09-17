#' @name bvar
#' @title Bayesian Vector Autoregression
#' @usage bvar(Yraw, plag, Wraw=NULL, draws=5000, burnin=5000, cons=FALSE, h=0,
#'             sv=TRUE, shrink.type="lagwise", sample_theta=FALSE,
#'             cfit=FALSE, crit_eig=1.00, prmean=1)
#' @param Data Data in matrix form
#' @param plag number of lags
#' @param draws number of saved draws.
#' @param burnin number of burn-ins.
#' @param prior which prior
#' @param SV If set to \code{TRUE} stochastic volatility is enabled.
#' @param h holdout-sample.
#' @param thin thinning factor
#' @param hyperpara hyperparameter set
#' @param eigen should eigenvalues be computed?
#' @param Ex exogenous variables to add to the model
#' @param cons If set to \code{TRUE} a constant is included.
#' @param trend If set to \code{TRUE} a trend is included.
#' @param applyfun parallelization
#' @param cores number of cores
#' @param verbose verbosity option
#' @export
#' @importFrom MASS ginv
#' @importFrom methods is
#' @importFrom stats is.ts median time ts
#' @importFrom xts is.xts
#' @importFrom zoo coredata
bvar<-function(Data,W,plag=1,draws=5000,burnin=5000,prior="NG",SV=TRUE,h=0,thin=1,hyperpara=NULL,eigen=FALSE,Ex=NULL,cons=FALSE,trend=FALSE,applyfun=NULL,cores=NULL,verbose=TRUE){
  start.bvar <- Sys.time()
  #--------------------------------- checks  ------------------------------------------------------#
  if(!is.matrix(Data)){
    stop("Please provide the argument 'Data' either as 'matrix' object.")
  }
  if(!is.null(Ex)){
    if(!is.list(Ex) & !is.matrix(Ex)){
      stop("Please provide the argument 'Ex' either as 'list' or as 'matrix' object.")
    }
  }
  if(!is.numeric(plag)){
    stop("Please specify number of lags as numeric.")
  }
  if(any(is.na(plag))){
    stop("Please specify number of lags.")
  }
  if(length(plag)>1 || plag<1){
    stop("Please specify number of lags accordingly. One lag length parameter for the whole model.")
  }
  if(!is.numeric(draws) | !is.numeric(burnin)){
    stop("Please specify number of draws and burnin as numeric.")
  }
  if(length(draws)>1 || draws<0 || length(burnin)>1 || burnin<0){
    stop("Please specify number of draws and burnin accordingly. One draws and burnin parameter for the whole model.")
  }
  if(prior%in%c("NC","MN","SSVS","NG")){
    stop("Please choose an available prior specification.")
  }
  #-------------------------- construct arglist ----------------------------------------------------#
  args <- .construct.arglist(bgvar)
  if(verbose){
    cat("\nStart estimation of Bayesian Vector Autoregression.\n\n")
    cat(paste("Prior: ",ifelse(prior=="MN","Minnesota prior",ifelse(prior=="SSVS","Stochastic Search Variable Selection prior","Normal-Gamma prior")),".\n",sep=""))
    cat(paste("Lag order: ",plag,"\n",sep=""))
    cat(paste("Stochastic volatility: ", ifelse(SV,"enabled","disabled"),".\n",sep=""))
  }
  #------------------------------ user checks  ---------------------------------------------------#
  # check Data
  if(is.matrix(Data)){
    if(any(is.na(Data))){
      stop("The data you have submitted contains NAs. Please check the data.")
    }
    isTS  <- is.ts(Data)
    isXTS <- is.xts(Data)
    Traw  <- nrow(Data)
    if(isTS || isXTS){
      temp       <- as.character(time(Data))
      years      <- unique(regmatches(temp,regexpr("^[0-9]{4}",temp)))
      months     <- temp
      for(kk in 1:length(years)) months <- gsub(paste(years[kk],"(\\.)?",sep=""),"",months)
      freq       <- length(unique(months))
      months     <- strtrim(months,3)
      startmonth <- ifelse(months[1]=="","01",ifelse(months[1]=="083","02",ifelse(months[1]=="166","03",ifelse(months[1]=="25","04",
                                                                                                               ifelse(months[1]=="333","05",ifelse(months[1]=="416","06",ifelse(months[1]=="5","07",ifelse(months[1]=="583","08",
                                                                                                                                                                                                           ifelse(months[1]=="666","09",ifelse(months[1]=="75","10",ifelse(months[1]=="833","11","12")))))))))))
      timeindex  <- seq.Date(from=as.Date(paste(years[1],"-",startmonth,"-01",sep=""), format="%Y-%m-%d"),
                             by=ifelse(freq==12,"months","quarter"), length.out = Traw)
      Data       <- ts(coredata(Data), start=c(as.numeric(years[1]),as.numeric(startmonth)),frequency=freq)
    }else{
      timeindex  <- seq.Date(from=as.Date("1830-08-01", format="%Y-%m-%d"), by="month", length.out = Traw)
      temp       <- coredata(Data)
      Data       <- ts(temp, start=c(1830,8), frequency=12)
    }
    args$time <- timeindex
    args$Traw <- length(timeindex)
  }
  args$Data <- Data
  # check truly exogenous variables
  if(!is.null(Ex)){
    if(is.matrix(Ex)){
      if(any(is.na(Ex))){
        stop("The data for exogenous variables you have submitted contains NAs. Please check the data.")
      }
      if(nrow(Ex)!=args$Traw){
        stop("Provided data and truly exogenous data not equally long. Please check.")
      }
    }
  }
  args$Ex <- Ex
  # check thinning factor
  if(thin<1){
    thin_mess <- paste("Thinning factor of ",thin," not possible. Adjusted to ",round(1/thin,2),".\n",sep="")
    thin <- round(1/thin,2)
  }
  if(draws%%thin!=0){
    thin_mess <- paste("Thinning factor of ",thin," no divisor of ",draws," (number of draws to save for posterior analysis).\n",sep="")
    div <- .divisors(draws,thin)
    thin <- min(div[which(abs(div-thin)==min(abs(div-thin)))])
    thin_mess <- paste(thin_mess,"New thinning factor: ", thin,". This means every", ifelse(thin==1,"",ifelse(thin==2,paste(" ",thin,"nd ",sep=""), ifelse(thin==3,paste(" ",thin,"rd ",sep=""),paste(" ",thin,"th ",sep="")))), "draw is saved.\n",sep="")
  }else{
    thin_mess <- paste("Thinning factor: ", thin,". This means every ",ifelse(thin==1,"",ifelse(thin==2,paste(thin,"nd ",sep=""),ifelse(thin==3,paste(thin,"rd ",sep=""),paste(thin,"th ",sep="")))),"draw is saved.\n",sep="")
  }
  if(verbose) cat(thin_mess)
  args$thindraws <- draws/thin
  # set default
  if(verbose) cat("Hyperparameter setup: \n")
  default_hyperpara <- list(c=10, # hyperparameter setup for natural conjugate case
                            a_1=0.01,b_1=0.01, prmean=0,# Gamma hyperparameter SIGMA (homoskedastic case) and mean
                            Bsigma=1, a0=25, b0=1.5, bmu=0, Bmu=100^2, # SV hyper parameter
                            shrink1=0.1,shrink2=0.2,shrink3=10^2,shrink4=0.1, # MN
                            tau0=.1,tau1=3,kappa0=0.1,kappa1=7,p_i=0.5,q_ij=0.5,   # SSVS
                            e_lambda=0.01,d_lambda=0.01,a_start=0.7,sample_A=FALSE,a_log=TRUE) # NG
  paras     <- c("c","a_1","b_1","prmean","Bsigma_sv","a0","b0","bmu","Bmu","shrink1","shrink2","shrink3",
                 "shrink4","tau0","tau1","kappa0","kappa1","p_i","q_ij","e_lambda","d_lambda","a_start","sample_A")
  if(is.null(hyperpara)){
    if(verbose) cat("\t No hyperparameters are chosen, default setting applied.\n")
  }
  if(!is.null(hyperpara)){
    for(para in names(hyperpara)){
      if(!para%in%paras){
        warning(paste0(para," no valid hyperparameter. Please check.\n"))
        next
      }
      default_hyperpara[para] <- hyperpara[para]
      if(para=="a_start") a_log <- FALSE
    }
    if(verbose) cat("Default values for chosen hyperparamters overwritten.\n")
  }
  #----------------------------------transform to matirx--------------------------------------------------------#
  Yraw <- as.matrix(Data)
  #---------------------------------hold out sample------------------------------------------------------------#
  args$yfull <- Yraw
  xglobal    <- Yraw[1:(nrow(Yraw)-h),,drop=FALSE]
  args$time  <- args$time[1:(length(args$time)-h)]
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
  #------------------------------ estimate BVAR ---------------------------------------------------------------#
  if(verbose) cat("\nEstimation of model starts... ")
  start.estim <- Sys.time()
  if(prior=="NC"){
    globalpost <- .BVAR_natural_conjugate()
  }else{
    globalpost <- .BVAR_linear_wrapper(Yraw=Yraw,prior=prior,plag=plag,draws=draws,burnin=burnin,cons=cons,trend=trend,SV=SV,thin=thin,default_hyperpara=default_hyperpara,Ex=Ex)
  }
  names(globalpost) <- cN
  end.estim <- Sys.time()
  diff.estim <- difftime(end.estim,start.estim,units="mins")
  mins <- round(diff.estim,0); secs <- round((diff.estim-floor(diff.estim))*60,0)
  if(verbose) cat(paste(" took ",mins," ",ifelse(mins==1,"min","mins")," ",secs, " ",ifelse(secs==1,"second.","seconds.\n"),sep=""))
  #--------------------------- stacking part for global model -----------------------------------------------------#
  if(is.logical(eigen)){
    if(eigen){trim<-1.05}else{trim<-NULL}
  }else{
    trim<-eigen;eigen<-TRUE
  }
  if(verbose) cat("Start stacking: \n")
  # insert stacking function here
  stacked.results <- .gvar.stacking.wrapper(xglobal=xglobal,plag=plag,globalpost=globalpost,draws=draws,thin=thin,trend=trend,eigen=eigen,trim=trim,verbose=verbose)
  if(!is.null(trim)) {args$thindraws <- length(stacked.results$F.eigen)}
  if(verbose) cat("\nStacking finished.\n")
  #---------------------- return output ---------------------------------------------------------------------------#
  out  <- structure(list("args"=args,
                         "xglobal"=xglobal,
                         "gW"=gW,
                         "stacked.results"=stacked.results,
                         "cc.results"=cc.results), class = "bvar")
  end.bvar <- Sys.time()
  diff.bvar <- difftime(end.bvar,start.bvar,units="mins")
  mins.bvar <- round(diff.bvar,0); secs.bvar <- round((diff.bvar-floor(diff.bvar))*60,0)
  if(verbose) cat(paste("\n Needed time for estimation of bgvar: ",mins.bvar," ",ifelse(mins.bvar==1,"min","mins")," ",secs.bvar, " ",ifelse(secs.bvar==1,"second.","seconds.\n"),sep=""))
  return(out)
}


#' @name bvar_old
#' @noRd
#' @export
bvar_old <- function(Yraw, plag, Wraw=NULL, draws=5000, burnin=5000, cons=FALSE, h=0, sv=TRUE, shrink.type = "lagwise",
                 sample_theta = FALSE, cfit=FALSE, crit_eig=1.00, prmean=1) {
  args <- c(as.list(environment()))
  varnames<-colnames(Yraw)
  Multiplier <- 5
  # if(!is.null(fixvar)) fixvar_pointer <- which(varnames == fixvar)
  if(!is.null(Wraw)) varnamesw <- colnames(Wraw) else varnamesw <- NULL
  #------------------------------Setup----------------------------------------------#
  ntot <- burnin+draws
  M <- ncol(Yraw)
  Xraw <- cbind(.mlag(Yraw,plag))
  Y <- Yraw[(plag+1):(nrow(Yraw)-h),]
  X <- Xraw[(plag+1):(nrow(Xraw)-h),]
  if(!is.null(Wraw)){
    Wex <- Wraw[(plag+1):(nrow(Wraw)),,drop=FALSE]
    X <- cbind(X,Wex)
    w<-ncol(Wex)
  }else{w<-0}
  if(cons) {X <- cbind(X,1);c<-1}else{c<-0}
  bigT <- nrow(X)
  K <- M*plag
  k <- ncol(X)
  v <- (M*(M-1))/2
  #--------------------------Initialize Gibbs sampler--------------------------------#
  A_OLS  <- try(solve(crossprod(X)) %*% crossprod(X,Y),silent=TRUE)
  if(is(A_OLS,"try-error")) A_OLS <- ginv(crossprod(X)) %*% crossprod(X,Y)
  Em    <- Em.str    <- Y - X %*% A_OLS
  SIGMA <- SIGMA_OLS <- crossprod(Y - X %*% A_OLS)/(bigT-k)
  #----------------------------PRIORS-----------------------------------------------#
  A_draw <- A_OLS
  # prior mean for autoregressive parameters
  A_prior <- matrix(0,k,M)
  A_prior[1:M,1:M] <- diag(M)*prmean

  # Normal-Gamma prior stuff
  tau2_draw <- matrix(1000, k, M)
  #if(cons) tau2_draw[k,] <- 1000
  lambda2_draw <- matrix(10^2,1,1)
  theta_draw <- 0.3
  cl0 <- 0.01
  dl0 <- 0.01
  tau2_scal_draw <- matrix(.1, M, M)
  lambda2_scal_draw <- matrix(10^2,1,1)
  theta_scal_draw <- 0.3
  # MH zeugs
  scale_coef <- 0.43
  accept_coef <- 0
  scale_scal <- 0.43
  accept_scal <- 0

  # variances
  # lower triangular cholesky matrix
  H_draw <- matrix(0, c(M, M))
  H_draw <- lower.tri(SIGMA_OLS)*1
  diag(H_draw) <- 1

  #SV
  svdraw <- list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,bigT))
  svl <- list()
  for (jj in 1:M) svl[[jj]] <- svdraw

  hv <- svdraw$latent
  para <- list(mu=-3,phi=.9,sigma=.2)
  Sv_draw <- matrix(-3,bigT,M) # log volatilties
  eta_list <- list()

  # no sv
  a_i=0.01
  b_i=0.01

  #---------------------------------Create storage matrices--------------------------#
  A_store                <- array(NA, c(draws, k, M))
  dimnames(A_store)[[3]] <- varnames
  Em_store               <- array(NA, c(draws, bigT, M))
  SIGMA_store            <- array(NA, c(draws, M, M))
  tau2_store             <- array(NA, c(draws, k, M))
  tau2_scal_store        <- array(NA, c(draws, M, M))
  lambda2_store          <- array(NA, c(draws, 1))
  lambda2_scal_store     <- array(NA, c(draws, 1))
  if(sample_theta) {
    theta_coef_store     <- array(NA, c(draws, 1))
    theta_scal_store     <- array(NA, c(draws, 1))
  }
  H_store                <- array(NA, c(draws, M, M))
  if(sv) pars_store      <- array(NA, c(draws, 3, M))
  Sv_store               <- array(NA, c(draws, bigT, M))
  if(cfit) {
    fit_store            <- array(NA, c(draws, bigT, M))
    Lik_store            <- array(NA, c(draws, 1))
  }
  #---------------------------Gibbs loop---------------------------------------#
  counter<-0
  irep<-1
  while(irep < (ntot+1)) {
    #----------------------------------------------------------------------------------------
    # Step I: Sample autoregressive parameters
    coefs <- .drawVARcoef(Y=Y, X=X, Sv=Sv_draw, aprior=A_prior,
                          Vprior=tau2_draw, Hprior=tau2_scal_draw)

    A_draw                     <- coefs$A
    H_draw                     <- coefs$H
    eta_list                   <- coefs$eta
    Em                         <- coefs$Em
    Em.str                     <- coefs$Em.str
    #----------------------------------------------------------------------------------------
    # Step II: Normal-Gamma prior on T_draw
    if(irep==1&shrink.type=="lagwise") {
      lambda2_draw  <- matrix(0.01,plag,1)
      lambda2_store <- matrix(0,draws,plag)
      theta_draw    <- matrix(theta_draw,plag,1)
      if(sample_theta) {
        theta_store   <- matrix(0,draws,plag)
        scale_coef    <- rep(scale_coef,plag)
        accept_coef   <- rep(accept_coef,plag)
      }
    }
    ng <- .drawNG(coef=A_draw[1:K,], tau2=tau2_draw[1:K,], lambda2=lambda2_draw,
                  theta=theta_draw, cl0=cl0, dl0=dl0, shrink.type=shrink.type,
                  sample_theta = sample_theta, scale = scale_coef, accept=accept_coef,
                  irep = irep, burnin = burnin, c=0, w=w)
    tau2_draw[1:K,]   <- ng$tau2
    lambda2_draw      <- ng$lambda2
    if(sample_theta) {
      theta_draw      <- ng$theta
      scale_coef      <- ng$scale
      accept_coef     <- ng$accept
    }
    #----------------------------------------------------------------------------------------
    # Step III: Sample covariance parameters
    if(sv) {
      pars_var <- matrix(0,3,M)
      for (jj in 1:M){
        svdraw <- svsample2(Em.str[,jj],
                            startpara=para(svl[[jj]]),
                            startlatent=Sv_draw[,jj],
                            priorphi = c(25, 1.5),
                            priorsigma=1)
        svl[[jj]] <- svdraw
        S_ <- exp(svdraw$latent)
        pars_var[,jj] <- para(svl[[jj]])
        Sv_draw[,jj] <- log(S_)
      }
    }else{
      for (mm in 1:M){
        S_1 <- a_i/2+bigT/2
        S_2 <- b_i/2+crossprod(Em.str[,mm])/2

        sig_eta <- 1/rgamma(1,S_1,S_2)
        Sv_draw[,mm] <- log(sig_eta)
      }
    }
    # compute SIGMA
    SIGMA <- H_draw %*% diag(exp(Sv_draw[bigT,])) %*% t(H_draw)
    # Step IV: Draw prior scaling factors from GIG for covariances
    ng <- .drawNG(coef=eta_list, tau2=tau2_scal_draw,
                  lambda2=lambda2_scal_draw,
                  theta=theta_scal_draw,
                  cl0=cl0, dl0=dl0, shrink.type = "covariance",
                  sample_theta=sample_theta, scale=scale_scal,
                  accept=accept_scal, irep=irep, burnin=burnin)
    tau2_scal_draw    <- ng$tau2
    lambda2_scal_draw <- ng$lambda2
    if(sample_theta) {
      theta_scal_draw <- ng$theta
      accept_scal     <- ng$accept
      scale_scal      <- ng$scale
    }
    # Step V: Check Stationarity
    Cm <- .gen_compMat(A_draw[1:K,],M,plag)$Cm
    if(max(abs(Re(eigen(Cm)$values)))>crit_eig && irep > burnin && counter < (burnin+Multiplier*draws)){
      irep    <- irep
      counter <- counter+1
      next
    }
    # Step VI: Store draws after burn-in/ Compute forecasts/ Impulse responses etc.
    if(irep > burnin) {
      if(cfit) {
        fit<-NULL;Lik<-0
        for(t in 1:bigT){
          fit <- rbind(fit,X[t,]%*%A_draw)
          # use faster version of multivariate normal from mvnfast package
          #Lik <- Lik + dmvnorm(Y[t,sl.country],fitt[t,],SIGMA[,,cc],log=TRUE)
          Lik <- Lik + dmvn(Y[t,],fit[t,],SIGMA,log=TRUE)
        }
      }

      A_store[irep-burnin,,]           <- A_draw
      SIGMA_store[irep-burnin,,]       <- SIGMA
      Em_store[irep-burnin,,]          <- Em
      tau2_store[irep-burnin,,]        <- tau2_draw
      tau2_scal_store[irep-burnin,,]   <- tau2_scal_draw
      lambda2_store[irep-burnin,]      <- lambda2_draw
      lambda2_scal_store[irep-burnin,] <- lambda2_scal_draw
      if(sample_theta) {
        theta_store[irep-burnin,]      <- theta_draw
        theta_scal_store[irep-burnin,] <- theta_scal_draw
      }
      H_store[irep-burnin,,]           <- solve(H_draw)
      if(sv) pars_store[irep-burnin,,] <- pars_var
      Sv_store[irep-burnin,,]          <- Sv_draw
      if(cfit) {
        fit_store[irep-burnin,,]       <- fit
        Lik_store[irep-burnin,]        <- Lik
      }
    } # END OF STEP VI
    irep    <- irep+1
    counter <- counter+1
    if(irep%%50==0) print(paste0("Round: ",irep))
  }
  print(paste0("Needed rounds: ", counter, " for ", irep, " total rounds."))
  #------------------------EX POST STUFF-------------------------------------#
  store <- post <-  NULL
  namespace    <- ls()
  #---------------------store draws-------------------------------------------#
  namestore    <- namespace[grepl("_store",namespace)]
  store        <- lapply(namestore, get, envir=sys.frame(sys.parent(0)))
  names(store) <- gsub("_store","",namestore)
  #---------------------compute posteriors-------------------------------------#
  quantile_set <- c(0.05,0.10,0.16,0.5,0.84,0.9,0.95)
  quantile_nam <- c("q05","q10","q16","q50","q84","q90","q95")

  for(nn in 1:length(namestore)){
    temp        <- get(namestore[nn], envir=sys.frame(sys.parent(0)))
    dims        <- 2:length(dim(temp))
    post        <- lapply(quantile_set, function(q) apply(temp, dims, quantile, q,na.rm=TRUE))
    names(post) <- quantile_nam
    assign(gsub("store","post",namestore[nn]),post)
  }

  namepost     <- ls()[grepl("_post",ls())]
  post         <- lapply(namepost, get, envir=sys.frame(sys.parent(0)))
  names(post)  <- gsub("_post","",namepost)

  # define output
  out <- structure(list(post=post,
                        store=store,
                        args=args), class="bvar")
  return(out)
}

#' @name bvar_natconj
#' @title Bayesian Vector Autoregression with Natural Conjugate Prior setup
#' @usage bvar(Yraw, plag, Wraw=NULL, draws=5000, burnin=5000, cons=FALSE, h=0,
#'             sv=TRUE, shrink.type="lagwise", sample_theta=FALSE,
#'             cfit=FALSE, crit_eig=1.00, prmean=1)
#' @param Yraw Data in matrix form
#' @param plag number of lags
#' @param Wraw exogenous variables.
#' @param draws number of saved draws.
#' @param burnin number of burn-ins.
#' @param cons If set to \code{TRUE} a constant is included.
#' @param h holdout-sample.
#' @param sv If set to \code{TRUE} stochastic volatility is enabled.
#' @param shrink.type which kind of shrinkage should be applied. Either \code{global} or \code{lagwise}.
#' @param sample_theta whether \code{theta} should be sampled
#' @param cfit whether log-likelihood should be computed
#' @param crit_eig critical eigenvalue
#' @param prmean prior mean.
#' @export
#' @importFrom stats rgamma rnorm quantile
#' @importFrom mvnfast dmvn
#' @importFrom mvtnorm rmvnorm dmvnorm rmvt
#' @importFrom stochvol svsample2
#' @importFrom MASS ginv
#' @importFrom methods is
bvar_natconj <- function(Yraw, plag, Wraw=NULL, draws=5000, cons=FALSE, h=0, cfit=FALSE, crit_eig=1.00, prmean=1) {
  args <- c(as.list(environment()))
  varnames<-colnames(Yraw)
  Multiplier <- 5
  c <- 10 # bigger => less information
  # if(!is.null(fixvar)) fixvar_pointer <- which(varnames == fixvar)
  if(!is.null(Wraw)) varnamesw <- colnames(Wraw) else varnamesw <- NULL
  #------------------------------Setup----------------------------------------------#
  ntot <- burnin+draws
  M <- ncol(Yraw)
  Xraw <- cbind(.mlag(Yraw,plag))
  Y <- Yraw[(plag+1):(nrow(Yraw)-h),]
  X <- Xraw[(plag+1):(nrow(Xraw)-h),]
  if(!is.null(Wraw)){
    Wex <- Wraw[(plag+1):(nrow(Wraw)),,drop=FALSE]
    X <- cbind(X,Wex)
    w<-ncol(Wex)
  }else{w<-0}
  if(cons) {X <- cbind(X,1);c<-1}else{c<-0}
  bigT <- nrow(X)
  K <- M*plag
  k <- ncol(X)
  v <- (M*(M-1))/2
  #--------------------------Initialize Gibbs sampler--------------------------------#
  A_OLS  <- try(solve(crossprod(X)) %*% crossprod(X,Y),silent=TRUE)
  if(is(A_OLS,"try-error")) A_OLS <- ginv(crossprod(X)) %*% crossprod(X,Y)
  S_OLS <- crossprod(Y - X %*% A_OLS)
  #----------------------------PRIORS------------------------------------------------#
  # prior mean for autoregressive parameters
  A_prior <- matrix(0,k,M)
  A_prior[1:M,1:M] <- diag(M)*prmean
  # prior variance for autoregressive coefficients
  V_prior    <- diag(k)
  V_priorinv <- diag(1/diag(V_prior))
  # prior degrees of scaling
  v_prior <- 1
  # prior scaling matrix
  S_prior <- (1/c)*diag(M)
  #---------------------POSTERIOR MOMENTS--------------------------------------------#
  # posterior of coefficients
  V_post <- solve(crossprod(X) + V_priorinv)
  # A_post <- V_post %*% (crossprod(X)%*%A_OLS + V_priorinv%*%A_prior)
  A_post <- V_post %*% (crossprod(X,Y) + V_priorinv%*%A_prior)
  # posterior of variance
  S_post <- S_OLS + S_prior + t(X%*%A_OLS)%*%X%*%A_OLS + t(A_prior)%*%V_priorinv%*%A_prior - t(A_post)%*%(V_priorinv + crossprod(X))%*%A_post
  v_post <- v_prior + bigT
  # posterior of coefficient variance for t-distribution
  bigVpost  <- kronecker(S_post, V_post)/(v_post-M-1)
  #---------------------------STORAGE-----------------------------------------------#
  A_store                <- array(NA, c(draws, k, M))
  dimnames(A_store)[[3]] <- varnames
  Em_store               <- array(NA, c(draws, bigT, M))
  S_store                <- array(NA, c(draws, M, M))
  if(cfit) {
    fit_store            <- array(NA, c(draws, bigT, M))
    Lik_store            <- array(NA, c(draws, 1))
  }
  #-------------------MONTE CARLO SIMULATION----------------------------------------#
  counter<-0
  irep<-1
  while(irep < (draws+1)){
    # draw coefficients
    Sinv_draw <- matrix(rWishart(1,v_post,solve(S_post)),M,M)
    S_draw    <- solve(Sinv_draw)
    A_draw    <- matrix(mvtnorm::rmvt(1, sigma=bigVpost, df=v_post, delta=as.vector(A_post)),k,M)

    # check stationarity
    Cm <- .gen_compMat(A_draw[1:K,],M,plag)$Cm
    if(max(abs(Re(eigen(Cm)$values)))>crit_eig && irep > burnin && counter < (burnin+Multiplier*draws)){
      irep    <- irep
      counter <- counter+1
      next
    }
    # save everything
    A_store[irep,,]  <- A_draw
    S_store[irep,,]  <- S_draw
    Em_store[irep,,] <- Y - X%*%A_draw

    if(cfit) {
      fit<-NULL;Lik<-0
      for(tt in 1:bigT){
        fit <- rbind(fit,X[tt,]%*%A_draw)
        # use faster version of multivariate normal from mvnfast package
        #Lik <- Lik + dmvnorm(Y[t,sl.country],fitt[t,],SIGMA[,,cc],log=TRUE)
        Lik <- Lik + mvnfast::dmvn(Y[tt,],fit[tt,],S_draw,log=TRUE)
      }
    }
    irep    <- irep+1
    counter <- counter+1
    if(irep%%50==0) print(paste0("Round: ",irep))
  }
  print(paste0("Needed rounds: ", counter, " for ", irep, " total rounds."))
  #------------------------EX POST STUFF-------------------------------------#
  store <- post <-  NULL
  namespace    <- ls()
  #---------------------store draws-------------------------------------------#
  namestore    <- namespace[grepl("_store",namespace)]
  store        <- lapply(namestore, get, envir=sys.frame(sys.parent(0)))
  names(store) <- gsub("_store","",namestore)
  #---------------------compute posteriors-------------------------------------#
  quantile_set <- c(0.05,0.10,0.16,0.5,0.84,0.9,0.95)
  quantile_nam <- c("q05","q10","q16","q50","q84","q90","q95")

  for(nn in 1:length(namestore)){
    temp        <- get(namestore[nn], envir=sys.frame(sys.parent(0)))
    dims        <- 2:length(dim(temp))
    post        <- lapply(quantile_set, function(q) apply(temp, dims, quantile, q,na.rm=TRUE))
    names(post) <- quantile_nam
    assign(gsub("store","post",namestore[nn]),post)
  }

  namepost     <- ls()[grepl("_post",ls())]
  post         <- lapply(namepost, get, envir=sys.frame(sys.parent(0)))
  names(post)  <- gsub("_post","",namepost)

  # define output
  out <- structure(list(post=post,
                        store=store,
                        args=args), class="bvar")
  return(out)
}

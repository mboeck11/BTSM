#' @name bvec
#' @title Bayesian Vector Error Correction Model
#' @usage bvec(Data, plag=1, r=1, beta=NULL, draws=5000,burnin=5000,prior="NG",SV=TRUE,h=0,thin=1,
#'             hyperpara=NULL,eigen=FALSE,
#'             Ex=NULL,cons=FALSE,trend=FALSE,applyfun=NULL,cores=NULL,verbose=TRUE)
#' @param Data Data in matrix form
#' @param plag number of lags
#' @param r number of cointegration relationships.
#' @param beta cointegration vector can be supplied, default is set to \code{NULL}.
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
#' @importFrom GIGrvg rgig
#' @importFrom parallel parLapply mclapply
#' @importFrom Rcpp evalCpp
#' @importFrom stats is.ts median time ts
#' @importFrom xts is.xts
#' @importFrom zoo coredata
bvec<-function(Data,plag=1,r=1,beta=NULL,draws=5000,burnin=5000,prior="NG",SV=TRUE,h=0,thin=1,hyperpara=NULL,eigen=FALSE,Ex=NULL,cons=FALSE,trend=FALSE,applyfun=NULL,cores=NULL,verbose=TRUE){
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
  if(!prior%in%c("NG")){
    stop("Please choose an available prior specification.")
  }
  #-------------------------- construct arglist ----------------------------------------------------#
  args <- .construct.arglist(bvar)
  if(verbose){
    cat("\nStart estimation of Bayesian Vector Error Correction Model.\n\n")
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
  default_hyperpara <- list(a_1=0.01,b_1=0.01, prmean=1,# Gamma hyperparameter SIGMA (homoskedastic case) and mean
                            Bsigma=1, a0=25, b0=1.5, bmu=0, Bmu=100^2, # SV hyper parameter
                            shrink1=0.1,shrink2=0.2,shrink3=10^2,shrink4=0.1, # MN
                            tau0=.1,tau1=3,kappa0=0.1,kappa1=7,p_i=0.5,q_ij=0.5,   # SSVS
                            e_lambda=0.01,d_lambda=0.01,a_start=0.7,sample_A=FALSE,a_log=TRUE) # NG
  paras     <- c("a_1","b_1","prmean","Bsigma_sv","a0","b0","bmu","Bmu","shrink1","shrink2","shrink3",
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
      if(para=="a_start") default_hyperpara["a_log"] <- FALSE
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
  if(verbose) cat("\nEstimation of model starts...\n")
  globalpost <- .BVEC_linear_wrapper(Yraw=Yraw,r=r,beta=beta,prior=prior,plag=plag,draws=draws,burnin=burnin,cons=cons,trend=trend,SV=SV,thin=thin,default_hyperpara=default_hyperpara,Ex=Ex)
  #--------------------------- checking eigenvalues ----------------------------------------------------------#
  if(is.logical(eigen)){
    if(eigen){trim<-1.00}else{trim<-NULL}
  }else{
    trim<-eigen;eigen<-TRUE
  }
  if(eigen){
    A.eigen <- applyfun(1:args$thindraws,function(irep){
      Cm <- .gen_compMat(globalpost$store$A_store[irep,,],ncol(Yraw),plag+1)$Cm
      return(max(abs(Re(eigen(Cm)$values))))
    })
    globalpost$post$A.eigen <- unlist(A.eigen)
  }
  #---------------------- return output ---------------------------------------------------------------------------#
  out  <- structure(list("args"=args,
                         "xglobal"=xglobal,
                         "post"=globalpost$post,
                         "store"=globalpost$store), class = "bvec")
  end.bvar <- Sys.time()
  diff.bvar <- difftime(end.bvar,start.bvar,units="mins")
  mins.bvar <- round(diff.bvar,0); secs.bvar <- round((diff.bvar-floor(diff.bvar))*60,0)
  if(verbose) cat(paste("\n Needed time for estimation of bvar: ",mins.bvar," ",ifelse(mins.bvar==1,"min","mins")," ",secs.bvar, " ",ifelse(secs.bvar==1,"second.","seconds.\n"),sep=""))
  return(out)
}

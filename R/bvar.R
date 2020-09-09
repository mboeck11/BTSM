#' @name bvar
#' @title Bayesian Vector Autoregression
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom GIGrvg rgig

bvar <- function(Yraw, Wraw=NULL, p, nsave=5000, nburn=5000, cons=FALSE, fcst = FALSE,
                 fhorz=0, h=0, sv=TRUE, shrink.type = "lagwise",
                 sample_theta = FALSE, cfit=FALSE, store.post=TRUE, crit_eig=1.00,
                 Multiplier=5,prmean=1) {
  args <- c(as.list(environment()))
  #---------------------------load Packages------------------------------------------#
  # require(mvtnorm, quietly=TRUE)
  # # require(Rcpp, quietly=TRUE)
  # # require(RcppArmadillo, quietly=TRUE)
  # require(MASS, quietly=TRUE)
  # require(MCMCpack, quietly=TRUE)
  # require(magic, quietly=TRUE)
  # require(stochvol,quietly=TRUE)
  # require(mvnfast, quietly=TRUE)
  # require(stringr, quietly=TRUE)
  # # sourceCpp("./functions/comp_irf2.cpp")
  # # sourceCpp("./functions/do_rgig.cpp")
  # require(GIGrvg, quietly=TRUE)
  #------------------------------checks----------------------------------------------#
  #if(is.null(cN)) cN <- unique(str_extract(colnames(Yraw), "^[a-z,A-Z]+"))
  #varnames <- unique(str_extract(colnames(Yraw),"[a-z,A-Z]+$")) # this does not work always
  # variables are separated by ".", first comes the country name, then the variable name, e.g., AT.eq
  # important, don't use "." as a separator in the variable names, only to separate the country name from the
  # variable name (i.e., not AT.eq.2 but AT.eq2)
  #max.char<-max(nchar(colnames(Yraw)))
  varnames<-colnames(Yraw)
  # if(!is.null(fixvar)) fixvar_pointer <- which(varnames == fixvar)
  if(!is.null(Wraw)) varnamesw <- colnames(Wraw) else varnamesw <- NULL
  #------------------------------Setup----------------------------------------------#
  ntot <- nburn+nsave
  M <- ncol(Yraw)
  Xraw <- cbind(mlag(Yraw,p))
  Y <- Yraw[(p+1):(nrow(Yraw)-h),]
  X <- Xraw[(p+1):(nrow(Xraw)-h),]
  if(!is.null(Wraw)){
    Wex <- Wraw[(p+1):(nrow(Wraw)),,drop=FALSE]
    X <- cbind(X,Wex)
    w<-ncol(Wex)
  }else{w<-0}
  if(cons) {X <- cbind(X,1);c<-1}else{c<-0}
  bigT <- nrow(X)
  K <- M*p
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
  A_store                <- array(NA, c(nsave, k, M))
  dimnames(A_store)[[3]] <- varnames
  Em_store               <- array(NA, c(nsave, bigT, M))
  SIGMA_store            <- array(NA, c(nsave, M, M))
  tau2_store             <- array(NA, c(nsave, k, M))
  tau2_scal_store        <- array(NA, c(nsave, M, M))
  lambda2_store          <- array(NA, c(nsave, 1))
  lambda2_scal_store     <- array(NA, c(nsave, 1))
  if(sample_theta) {
    theta_coef_store     <- array(NA, c(nsave, 1))
    theta_scal_store     <- array(NA, c(nsave, 1))
  }
  H_store                <- array(NA, c(nsave, M, M))
  if(sv) pars_store      <- array(NA, c(nsave, 3, M))
  Sv_store               <- array(NA, c(nsave, bigT, M))
  if(cfit) {
    fit_store            <- array(NA, c(nsave, bigT, M))
    Lik_store            <- array(NA, c(nsave, 1))
  }
  if(fcst) {
    pred_fcst            <- array(NA, c(nsave, fhorz, M))
    dimnames(pred_fcst)[[3]] <- varnames
    rmse_fcst            <- array(NA, c(nsave, h))
    lps_fcst             <- array(NA, c(nsave, h))
  }
  #---------------------------Gibbs loop---------------------------------------#
  counter<-0
  irep<-1
  while(irep < (ntot+1)) {
    #----------------------------------------------------------------------------------------
    # Step I: Sample autoregressive parameters
    coefs <- drawVARcoef(Y=Y, X=X, Sv=Sv_draw, aprior=A_prior,
                         Vprior=tau2_draw, Hprior=tau2_scal_draw)

    A_draw                     <- coefs$A
    H_draw                     <- coefs$H
    eta_list                   <- coefs$eta
    Em                         <- coefs$Em
    Em.str                     <- coefs$Em.str
    #----------------------------------------------------------------------------------------
    # Step II: Normal-Gamma prior on T_draw
    if(irep==1&shrink.type=="lagwise") {
      lambda2_draw  <- matrix(0.01,p,1)
      lambda2_store <- matrix(0,nsave,p)
      theta_draw    <- matrix(theta_draw,p,1)
      if(sample_theta) {
        theta_store   <- matrix(0,nsave,p)
        scale_coef    <- rep(scale_coef,p)
        accept_coef   <- rep(accept_coef,p)
      }
    }
    ng <- NormalGammaPrior(coef=A_draw[1:K,], tau2=tau2_draw[1:K,], lambda2=lambda2_draw,
                           theta=theta_draw, cl0=cl0, dl0=dl0, shrink.type=shrink.type,
                           sample_theta = sample_theta, scale = scale_coef, accept=accept_coef,
                           irep = irep, nburn = nburn, c=0, w=w)
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
    ng <- NormalGammaPrior(coef=eta_list, tau2=tau2_scal_draw,
                           lambda2=lambda2_scal_draw,
                           theta=theta_scal_draw,
                           cl0=cl0, dl0=dl0, shrink.type = "covariance",
                           sample_theta=sample_theta, scale=scale_scal,
                           accept=accept_scal, irep=irep, nburn=nburn)
    tau2_scal_draw    <- ng$tau2
    lambda2_scal_draw <- ng$lambda2
    if(sample_theta) {
      theta_scal_draw <- ng$theta
      accept_scal     <- ng$accept
      scale_scal      <- ng$scale
    }
    # Step V: Check Stationarity
    Cm <- gen_compMat(A_draw[1:K,],M,p)$Cm
    if(max(abs(Re(eigen(Cm)$values)))>crit_eig && irep > nburn && counter < (nburn+Multiplier*nsave)){
      irep    <- irep
      counter <- counter+1
      next
    }
    # Step VI: Store draws after burn-in/ Compute forecasts/ Impulse responses etc.
    if(irep > nburn) {
      if(cfit) {
        fit<-NULL;Lik<-0
        for(t in 1:bigT){
          fit <- rbind(fit,X[t,]%*%A_draw)
          # use faster version of multivariate normal from mvnfast package
          #Lik <- Lik + dmvnorm(Y[t,sl.country],fitt[t,],SIGMA[,,cc],log=TRUE)
          Lik <- Lik + dmvn(Y[t,],fit[t,],SIGMA,log=TRUE)
        }
      }

      A_store[irep-nburn,,]           <- A_draw
      SIGMA_store[irep-nburn,,]       <- SIGMA
      Em_store[irep-nburn,,]          <- Em
      tau2_store[irep-nburn,,]        <- tau2_draw
      tau2_scal_store[irep-nburn,,]   <- tau2_scal_draw
      lambda2_store[irep-nburn,]      <- lambda2_draw
      lambda2_scal_store[irep-nburn,] <- lambda2_scal_draw
      if(sample_theta) {
        theta_store[irep-nburn,]      <- theta_draw
        theta_scal_store[irep-nburn,] <- theta_scal_draw
      }
      H_store[irep-nburn,,]           <- solve(H_draw)
      if(sv) pars_store[irep-nburn,,] <- pars_var
      Sv_store[irep-nburn,,]          <- Sv_draw
      if(cfit) {
        fit_store[irep-nburn,,]       <- fit
        Lik_store[irep-nburn,]        <- Lik
      }

      if(fcst){
        #Compute draw from the predictive density
        xt <- Xraw[(nrow(Xraw)-h):nrow(Xraw),,drop=FALSE]
        yt <- Yraw[(nrow(Yraw)-h):nrow(Yraw),,drop=FALSE]

        if(h==0) {
          if(p==1) Mean00 <- t(yt) else Mean00 <- t(cbind(yt,xt[,1:(M*(p-1)),drop=FALSE]))
        }else {
          Mean00 <- t(xt)
        }
        Sigma00 <- matrix(0,K,K)
        Sigma00[1:M,1:M] <- SIGMA
        Htt <- Sv_draw[bigT,]
        temp <- gen_compMat(A_draw,M,p)
        Jm <- temp$Jm
        Cm <- temp$Cm
        for(ih in 1:fhorz) {
          # first and second moments
          Mean00 <- Cm%*%Mean00
          Sigma00 <- Cm%*%Sigma00%*%t(Cm)+Jm%*%SIGMA%*%t(Jm)
          # add stuff
          meanvec <- Mean00[1:M,,drop=FALSE]
          if(cons) meanvec <- meanvec + A_draw[K+1,]
          if(!is.null(Wraw)) meanvec <- meanvec + Wraw%*%A_draw[(K+c):(K+c+w),]
          cholSig <- try(t(chol(Sigma00[1:M,1:M])), silent=TRUE)
          if(is(cholSig,"try-error")){
            yf <- rmvnorm(1,meanvec,Sigma00[1:M,1:M])
          } else {
            yf <- meanvec + cholSig%*%rnorm(M,0,1)
          }
          pred_fcst[irep-nburn,ih,] <- yf
          if(sv) {
            Htt <- pars_var[1,]+pars_var[2,]*(Htt-pars_var[1,])+rnorm(M,0,sqrt(pars_var[3,]))
            Sigmat <- H_draw %*% diag(exp(Htt))%*%t(H_draw)
          }
          # only if we use a holdout sample
          if((h-ih+1)>0){
            # endpoint forecast evaluation
            rmse_fcst[irep-nburn,ih] <- sum(sqrt((yf-yt[ih,])^2))
            lps_fcst[irep-nburn,ih]  <- dmvnorm(c(yf), mean = yt[ih,],
                                                sigma=Sigma00[1:M,1:M], log=TRUE)
          }
        }
      }
    } # END OF STEP VI
    irep    <- irep+1
    counter <- counter+1
    if(irep%%50==0) print(paste0("Round: ",irep))
  }
  print(paste0("Needed rounds: ", counter, " for ", irep, " total rounds."))
  #------------------------EX POST STUFF-------------------------------------#
  store <- fcststore <- post <-  NULL
  namespace    <- ls()
  #---------------------store draws-------------------------------------------#
  namestore    <- namespace[grepl("_store",namespace)]
  store        <- lapply(namestore, get, envir=sys.frame(sys.parent(0)))
  names(store) <- gsub("_store","",namestore)
  #-----------------------forecasts--------------------------------------------#
  if(fcst){
    namefcst     <- namespace[grepl("_fcst",namespace)]
    fcststore    <- lapply(namefcst, get, envir=sys.frame(sys.parent(0)))
    names(fcststore) <- gsub("_fcst","",namefcst)
  }
  #---------------------compute posteriors-------------------------------------#
  if(store.post) {
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
  }

  out <- structure(list(post=post,
                        store=store,
                        fcststore=fcststore,
                        args=args), class="bvar")

  return(out)
}

#' @name bvar
#' @title Bayesian Vector Autoregression
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
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom stochvol svsample2
#' @importFrom MASS ginv
#' @importFrom methods is
bvar <- function(Yraw, plag, Wraw=NULL, draws=5000, burnin=5000, cons=FALSE, h=0, sv=TRUE, shrink.type = "lagwise",
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

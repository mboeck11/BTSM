#' @name .TVPBVAR_linear_wrapper
#' @noRd
#' @importFrom abind adrop
#' @importFrom utils capture.output
.TVPBVAR_linear_wrapper <- function(Yraw, prior, plag, draws, burnin, cons, trend, SV, thin, default_hyperpara, Ex, applyfun, cores, eigen, trim){
  class(Yraw) <- "numeric"
  prior_in <- prior
  if(default_hyperpara[["a_log"]]) default_hyperpara["a_start"] <- 1/log(ncol(Yraw))
  if(prior=="TVP" || prior=="TVP-NG"){
    prior_in <- ifelse(prior=="TVP",1,2)
    post_draws <- applyfun(1:ncol(Yraw), function(nr){
      .TVPBVAR_noncentered_R(nr=nr,Y_in=Yraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)
    })
    tvpbvar <- .var_posterior(post_draws, prior, draws/thin, applyfun, cores)
  }else if(prior=="TTVP"){
    prior_in <- 3
    post_draws <- applyfun(1:ncol(Yraw), function(nr){
      .TVPBVAR_centered_R(nr=nr,Y_in=Yraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)
    })
    tvpbvar <- .var_posterior(post_draws, prior, draws/thin, applyfun, cores)
  }
  #------------------------------------------------ get data ----------------------------------------#
  Y <- tvpbvar$Y; colnames(Y) <- colnames(Yraw); X <- tvpbvar$X
  M <- ncol(Y); bigT <- nrow(Y); K <- ncol(X)
  if(!is.null(Ex)) Mex <- ncol(Ex)
  names <- colnames(Yraw)
  if(is.null(names)) names <- rep("Y",M)
  xnames <-NULL
  for(ii in 1:plag) xnames <- c(xnames,paste0(names,".lag",ii))
  if(!is.null(Ex)) enames <- c(enames,paste(rep("Tex",Mex))) else enames <- NULL
  if(cons)  cnames <- "cons" else cnames <- NULL
  if(trend) tnames <- "trend" else tnames <- NULL
  colnames(X) <- c(xnames,enames,cnames,tnames)
  #-----------------------------------------get containers ------------------------------------------#
  A_store <- tvpbvar$A_store; dimnames(A_store)[[2]] <- paste("t",seq(1,bigT),sep="."); dimnames(A_store)[[3]] <- colnames(X); dimnames(A_store)[[4]] <- colnames(Y)
  # splitting up stores
  dims  <- dimnames(A_store)[[3]]
  a0_store <- a1_store <- Ex_store <- Phi_store <- NULL
  if(cons) a0_store <- A_store[,,which(dims%in%cnames),]
  if(trend) a1_store <- A_store[,,which(dims%in%tnames),]
  if(!is.null(Ex)) Ex_store <- A_store[,which(dims%in%enames),,drop=FALSE]
  for(jj in 1:plag) {
    xnames.jj <- xnames[grepl(paste0("lag",jj),xnames)]
    Phi_store[[jj]] <- A_store[,,which(dims%in%xnames.jj),]
    dimnames(Phi_store[[jj]]) <- list(NULL,paste("t",seq(1,bigT),sep="."),xnames.jj,names)
  }
  L_store <- tvpbvar$L_store
  S_store <- tvpbvar$S_store
  Smed_store <- tvpbvar$Smed_store
  vola_store <- tvpbvar$Sv_store; dimnames(vola_store) <- list(NULL,NULL,colnames(Y))
  if(SV){
    pars_store <- tvpbvar$pars_store; dimnames(pars_store) <- list(NULL,c("mu","phi","sigma","latent0"),colnames(Y))
  }else pars_store <- NULL
  res_store <- tvpbvar$res_store; dimnames(res_store) <- list(NULL,NULL,colnames(Y))
  # NG
  if(prior=="TVP"){
    thetasqrt_store<- tvpbvar$thetasqrt_store
    Lthetasqrt_store<-tvpbvar$Lthetasqrt_store
    tau2_store<-xi2_store<-lambda2_store<-kappa2_store<-a_tau_store<-a_xi_store<-Ltau2_store<-Lxi2_store <- NULL
    D_store<-Omega_store<-thrsh_store<-kappa_store<-V0_store<-LD_store<-LOmega_store<-Lthrsh_store<-LV0_store <- NULL
  }else if(prior=="TVP-NG"){
    D_store<-Omega_store<-thrsh_store<-kappa_store<-V0_store<-LD_store<-LOmega_store<-Lthrsh_store<-LV0_store<-NULL
    thetasqrt_store <- tvpbvar$thetasqrt_store
    tau2_store      <- tvpbvar$tau2_store
    xi2_store       <- tvpbvar$xi2_store
    lambda2_store   <- tvpbvar$lambda2_store
    kappa2_store    <- tvpbvar$kappa2_store
    a_tau_store     <- tvpbvar$a_tau_store
    a_xi_store      <- tvpbvar$a_xi_store
    Lthetasqrt_store<-tvpbvar$Lthetasqrt_store
    Ltau2_store     <- tvpbvar$Ltau2_store
    Lxi2_store      <- tvpbvar$Lxi2_store
  }else if(prior=="TTVP"){
    thetasqrt_store<-Lthetasqrt_store<-tau2_store<-xi2_store<-lambda2_store<-kappa2_store<-a_tau_store<-a_xi_store<-Ltau2_store<-Lxi2_store<-NULL
    D_store       <- tvpbvar$D_store
    Omega_store   <- tvpbvar$Omega_store
    thrsh_store   <- tvpbvar$thrsh_store
    kappa_store   <- tvpbvar$kappa_store
    V0_store      <- tvpbvar$V0_store
    LD_store      <- tvpbvar$LD_store
    LOmega_store  <- tvpbvar$LOmega_store
    Lthrsh_store  <- tvpbvar$Lthrsh_store
    LV0_store     <- tvpbvar$LV0_store
  }
  if(eigen){
    # check medians: could be done more carefully
    A.eigen <- unlist(applyfun(1:(draws/thin),function(irep){
      Cm <- .gen_compMat(apply(A_store[irep,,,],c(2,3),median),ncol(Yraw),plag)$Cm
      return(max(abs(Re(eigen(Cm)$values))))
    }))
    trim_eigen <- which(A.eigen<trim)
    if(length(trim_eigen)==0) stop("No stable draws found. Either increase number of draws or trimming factor.")
    A_store<-A_store[trim_eigen,,,,drop=FALSE]
    if(cons) a0_store <- a0_store[trim_eigen,,,drop=FALSE]
    if(trend) a1_store <- a1_store[trim_eigen,,,drop=FALSE]
    if(!is.null(Ex)) Ex_store <- Ex_store[trim_eigen,,,,drop=FALSE]
    Phi_store<-lapply(Phi_store,function(l)l[trim_eigen,,,,drop=FALSE])
    L_store<-L_store[trim_eigen,,,,drop=FALSE]
    S_store<-S_store[trim_eigen,,,,drop=FALSE]
    Smed_store<-Smed_store[trim_eigen,,,drop=FALSE]
    vola_store<-vola_store[trim_eigen,,,drop=FALSE]
    if(SV) pars_store<-pars_store[trim_eigen,,,drop=FALSE]
    res_store<-res_store[trim_eigen,,,drop=FALSE]
    if(prior=="TVP"){
      thetasqrt_store<-thetasqrt_store[trim_eigen,,,drop=FALSE]
      Lthetasqrt_store<-lapply(Lthetasqrt_store,function(l)l[trim_eigen,,drop=FALSE])
    }
    if(prior=="TVP-NG"){
      thetasqrt_store<-thetasqrt_store[trim_eigen,,,drop=FALSE]
      tau2_store<-tau2_store[trim_eigen,,,drop=FALSE]
      xi2_store<-xi2_store[trim_eigen,,,drop=FALSE]
      lambda2_store<-lambda2_store[trim_eigen,,,drop=FALSE]
      kappa2_store<-kappa2_store[trim_eigen,,,drop=FALSE]
      a_tau_store<-a_tau_store[trim_eigen,,,drop=FALSE]
      a_xi_store<-a_xi_store[trim_eigen,,,drop=FALSE]
      Lthetasqrt_store<-lapply(Lthetasqrt_store,function(l)l[trim_eigen,,drop=FALSE])
      Ltau2_store<-lapply(Ltau2_store,function(l)l[trim_eigen,,drop=FALSE])
      Lxi2_store<-lapply(Lxi2_store,function(l)l[trim_eigen,,drop=FALSE])
    }else if(prior=="TTVP"){
      D_store<-D_store[trim_eigen,,,,drop=FALSE]
      Omega_store<-Omega_store[trim_eigen,,,,drop=FALSE]
      thrsh_store<-thrsh_store[trim_eigen,,,drop=FALSE]
      kappa_store<-kappa_store[trim_eigen,,,drop=FALSE]
      V0_store<-V0_store[trim_eigen,,,drop=FALSE]
      LD_store<-lapply(LD_store,function(l)l[trim_eigen,,,drop=FALSE])
      LOmega_store<-lapply(LOmega_store,function(l)l[trim_eigen,,,drop=FALSE])
      Lthrsh_store<-lapply(Lthrsh_store,function(l)l[trim_eigen,,drop=FALSE])
      LV0_store<-lapply(LV0_store,function(l)l[trim_eigen,,drop=FALSE])
    }
  }else{A.eigen<-NULL}
  store <- list(A_store=A_store,a0_store=a0_store,a1_store=a1_store,Phi_store=Phi_store,Ex_store=Ex_store,S_store=S_store,Smed_store=Smed_store,L_store=L_store,Lthetasqrt_store=Lthetasqrt_store,
                vola_store=vola_store,pars_store=pars_store,res_store=res_store,thetasqrt_store=thetasqrt_store,tau2_store=tau2_store,xi2_store=xi2_store,lambda2_store=lambda2_store,
                kappa2_store=kappa2_store,a_tau_store=a_tau_store,a_xi_store=a_xi_store,Ltau2_store=Ltau2_store,Lxi2_store=Lxi2_store,D_store=D_store,
                Omega_store=Omega_store,thrsh_store=thrsh_store,kappa_store=kappa_store,V0_store=V0_store,LD_store=LD_store,LOmega_store=LOmega_store,
                Lthrsh_store=Lthrsh_store,LV0_store=LV0_store,A.eigen=A.eigen)
  #------------------------------------ compute posteriors -------------------------------------------#
  A_post      <- apply(A_store,c(2,3,4),median)
  L_post      <- apply(L_store,c(2,3,4),median)
  S_post      <- apply(S_store,c(2,3,4),median)
  Smed_post   <- apply(Smed_store,c(2,3),median)
  Sig         <- apply(S_post,c(2,3),mean)/(bigT-K)
  res_post    <- apply(res_store,c(2,3),median)
  # splitting up posteriors
  a0_post <- a1_post <- Ex_post <- NULL
  if(cons)  a0_post <- A_post[,which(dims=="cons"),,drop=FALSE]
  if(trend) a1_post <- A_post[,which(dims=="trend"),,drop=FALSE]
  if(!is.null(Ex)) Ex_post <- A_post[,which(dims=="Tex"),,drop=FALSE]
  Phi_post<- NULL
  for(jj in 1:plag){
    Phi_post[[jj]]    <- A_post[,which(dims==paste("Ylag",jj,sep="")),,drop=FALSE]
  }
  vola_post <- apply(vola_store,c(2,3),median); dimnames(vola_post) <- list(NULL,colnames(Y))
  if(SV){
    pars_post <- apply(pars_store,c(2,3),median); dimnames(pars_post) <- list(c("mu","phi","sigma","latent0"),colnames(Y))
  }else pars_post <- NULL
  if(prior=="TVP"){
    thetasqrt_post<-apply(thetasqrt_store,c(2,3),median)
    Lthetasqrt_post<-lapply(Lthetasqrt_store,function(l)apply(l,2,median))
    tau2_post<-xi2_post<-lambda2_post<-kappa2_post<-a_tau_post<-a_xi_post<-Ltau2_post<-Lxi2_post<-NULL
    D_post<-Omega_post<-thrsh_post<-kappa_post<-V0_post<-LD_post<-LOmega_post<-Lthrsh_post<-LV0_post<-NULL
  }else if(prior=="TVP-NG"){
    D_post<-Omega_post<-thrsh_post<-kappa_post<-V0_post<-LD_post<-LOmega_post<-Lthrsh_post<-LV0_post<-NULL
    thetasqrt_post<-apply(thetasqrt_store,c(2,3),median)
    tau2_post <- apply(tau2_store,c(2,3),median)
    xi2_post  <- apply(xi2_store,c(2,3),median)
    lambda2_post <- apply(lambda2_store,c(2,3),median)
    kappa2_post <- apply(kappa2_store,c(2,3),median)
    a_tau_post <- apply(a_tau_store,c(2,3),median)
    a_xi_post <- apply(a_xi_store,c(2,3),median)
    Lthetasqrt_post<-lapply(Lthetasqrt_store,function(l)apply(l,2,median))
    Ltau2_post <- lapply(Ltau2_store,function(l)apply(l,c(2),median))
    Lxi2_post <- lapply(Lxi2_store,function(l)apply(l,c(2),median))
  }else if(prior=="TTVP"){
    thetasqrt_post<-Lthetasqrt_post<-tau2_post<-xi2_post<-lambda2_post<-kappa2_post<-a_tau_post<-a_xi_post<-Ltau2_post<-Lxi2_post<-NULL
    D_post <- apply(D_store,c(2,3.4),median)
    Omega_post <- apply(Omega_store,c(2,3,4),median)
    thrsh_post <- apply(thrsh_store,c(2,3),median)
    kappa_post <- apply(kappa_store,c(2,3),median)
    V0_post <- apply(V0_store,c(2,3),median)
    LD_post <- lapply(LD_store,function(l)apply(l,c(2,3),median))
    LOmega_post <- lapply(LOmega_store,function(l)apply(l,c(2,3),median))
    Lthrsh_post <- lapply(Lthrsh_store,function(l)apply(l,2,median))
    LV0_post <- lapply(LV0_store,function(l)apply(l,2,median))
  }
  post <- list(A_post=A_post,a0_post=a0_post,a1_post=a1_post,Phi_post=Phi_post,Ex_post=Ex_post,S_post=S_post,Smed_post=Smed_post,L_post=L_post,Lthetasqrt_post=Lthetasqrt_post,
                vola_post=vola_post,pars_post=pars_post,res_post=res_post,tau2_post=tau2_post,thetasqrt_post=thetasqrt_post,xi2_post=xi2_post,lambda2_post=lambda2_post,
                kappa2_post=kappa2_post,a_tau_post=a_tau_post,a_xi_post=a_xi_post,Ltau2_post=Ltau2_post,Lxi2_post=Lxi2_post,D_post=D_post,
                Omega_post=Omega_post,thrsh_post=thrsh_post,kappa_post=kappa_post,V0_post=V0_post,LD_post=LD_post,LOmega_post=LOmega_post,
                Lthrsh_post=Lthrsh_post,LV0_post=LV0_post)
  return(list(Y=Y,X=X,store=store,post=post))
}

#' @name .TVPBVAR_noncentered_R.m
#' @importFrom stochvol svsample_fast_cpp specify_priors default_fast_sv sv_normal sv_beta sv_gamma
#' @importFrom MASS ginv mvrnorm
#' @importFrom matrixcalc hadamard.prod
#' @importFrom methods is
#' @importFrom stats rnorm rgamma runif dnorm
#' @noRd
.TVPBVAR_noncentered_R <- function(nr,Y_in,p_in,draws_in,burnin_in,cons_in,trend_in,sv_in,thin_in,quiet_in,prior_in,hyperparam_in,Ex_in){
  #----------------------------------------INPUTS----------------------------------------------------#
  Yraw  <- Y_in
  p     <- p_in
  Traw  <- nrow(Yraw)
  M     <- ncol(Yraw)
  K     <- M*p
  Ylag  <- .mlag(Yraw,p)
  names <- colnames(Yraw)
  if(is.null(names)) names <- rep("Y",M)
  colnames(Yraw) <- names
  nameslags <- NULL
  for(ii in 1:p) nameslags <- c(nameslags,paste0(names,".lag",ii))
  colnames(Ylag) <- nameslags

  texo <- FALSE; Mex <- 0; Exraw <- NULL; enames <- NULL
  if(!is.null(Ex_in)){
    Exraw <- Ex_in; Mex <- ncol(Exraw); texo <- TRUE
    enames <- colnames(Exraw)
    if(is.null(enames)) enames <- rep("Tex",Mex)
    colnames(Exraw) <- enames
  }

  if(nr==1) slct <- NULL else slct <- 1:(nr-1)

  Xraw  <- cbind(Yraw[,slct],Ylag,Exraw)
  colnames(Xraw) <- c(colnames(Yraw)[slct],nameslags,enames)
  X     <- Xraw[(p+1):nrow(Xraw),,drop=FALSE]
  y     <- Yraw[(p+1):Traw,nr,drop=FALSE]
  bigT  <- nrow(X)
  M_    <- M-length(slct)

  cons  <- cons_in
  if(cons){
    X <- cbind(X,1)
    colnames(X)[ncol(X)] <- "cons"
  }
  trend <- trend_in
  if(trend){
    X <- cbind(X,seq(1,bigT))
    colnames(X)[ncol(X)] <- "trend"
  }

  d <- ncol(X)
  n <- d*M
  v <- (M*(M-1))/2
  #---------------------------------------------------------------------------------------------------------
  # HYPERPARAMETERS
  #---------------------------------------------------------------------------------------------------------
  hyperpara <- hyperparam_in
  prior     <- prior_in
  sv        <- sv_in
  prmean    <- hyperpara$prmean
  # non-SV
  c0        <- hyperpara$c0
  g0        <- hyperpara$g0
  # SV
  bmu       <- hyperpara$bmu
  Bmu       <- hyperpara$Bmu
  a0        <- hyperpara$a0
  b0        <- hyperpara$b0
  Bsigma    <- hyperpara$Bsigma
  # TVP-NG
  d1        <- hyperpara$d1
  d2        <- hyperpara$d2
  e1        <- hyperpara$e1
  e2        <- hyperpara$e2
  b_xi      <- hyperpara$b_xi
  b_tau     <- hyperpara$b_tau
  nu_xi     <- hyperpara$nu_xi
  nu_tau    <- hyperpara$nu_tau
  a_start   <- hyperpara$a_start
  sample_A  <- hyperpara$sample_A
  #---------------------------------------------------------------------------------------------------------
  # OLS Quantitites
  #---------------------------------------------------------------------------------------------------------
  XtXinv <- try(solve(crossprod(X)),silent=TRUE)
  if(is(XtXinv,"try-error")) XtXinv <- ginv(crossprod(X))
  A_OLS  <- XtXinv%*%(t(X)%*%y)
  E_OLS  <- y - X%*%A_OLS
  S_OLS  <- crossprod(E_OLS)/(bigT-d)
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  A_draw  <- matrix(A_OLS, bigT+1, d, byrow=TRUE, dimnames=list(NULL,colnames(X)))
  S_draw  <- matrix(S_OLS, bigT, 1)

  # time-varying stuff
  Am_draw    <- A_OLS
  At_draw    <- matrix(0, bigT+1, d)
  theta_draw <- rep(1,d)
  theta_sqrt <- sqrt(theta_draw)
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  # Priors on VAR coefs
  #-----------------------------
  # prior mean
  A_prior <- matrix(0,2*d, 1)
  A_prior[2*nr-1,1] <- prmean
  # prior variance
  tau2_draw <- rep(10,d)
  xi2_draw  <- rep(10,d)

  # NG stuff
  lambda2      <- 10
  a_tau        <- a_start
  scale_tau    <- .43
  acc_tau      <- 0

  kappa2       <- 10
  a_xi         <- a_start
  scale_xi     <- .43
  acc_xi       <- 0
  #------------------------------------
  # SV quantities
  #------------------------------------
  svdraw  <- list(para=c(mu=-10,phi=.9,sigma=.2,latent0=-3),latent=rep(-3,bigT))
  Sv_draw <- svdraw$latent
  pars_var <- matrix(c(-3,.9,.2,-3),4,1,dimnames=list(c("mu","phi","sigma","latent0"),NULL))
  Sv_priors <- specify_priors(mu=sv_normal(mean=bmu, sd=Bmu), phi=sv_beta(a0,b0), sigma2=sv_gamma(shape=0.5,rate=1/(2*Bsigma)))
  #-----------------------------------
  # non-SV quantities
  #-----------------------------------
  sig_eta <- exp(-3)
  G0 <- g0/S_OLS*(c0-1)
  #---------------------------------------------------------------------------------------------------------
  # SAMPLER MISCELLANEOUS
  #---------------------------------------------------------------------------------------------------------
  nsave <- draws_in
  nburn <- burnin_in
  ntot  <- nsave+nburn

  # thinning
  thin         <- thin_in
  count <- 0
  thindraws    <- nsave/thin
  thin.draws   <- seq(nburn+1,ntot,by=thin)
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store      <- array(NA,c(thindraws,bigT+1,d))
  Am_store     <- array(NA,c(thindraws,d,1))
  At_store     <- array(NA,c(thindraws,bigT+1,d))
  res_store    <- array(NA,c(thindraws,bigT,1))
  Sv_store     <- array(NA,c(thindraws,bigT,1))
  pars_store   <- array(NA,c(thindraws,4,1))
  # state variances
  thetasqrt_store <- array(NA,c(thindraws,d,1))
  # TVP-NG
  tau2_store   <- array(NA,c(thindraws,d,1))
  xi2_store    <- array(NA,c(thindraws,d,1))
  lambda2_store<- array(NA,c(thindraws,p,1))
  kappa2_store <- array(NA,c(thindraws,1,1))
  a_xi_store   <- array(NA,c(thindraws,1,1))
  a_tau_store  <- array(NA,c(thindraws,p,1))
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 0: Normalize data
    Xt <- apply(X,2,function(x)x*exp(-0.5*Sv_draw))
    yt <- y*exp(-0.5*Sv_draw)
    #----------------------------------------------------------------------------
    # Step 1: Sample coefficients
    Zt <- cbind(Xt,hadamard.prod(Xt,At_draw[2:(bigT+1),]))

    Vpriorinv <- diag(1/c(tau2_draw,xi2_draw))
    # V_post <- try(chol2inv(chol(crossprod(Zt)+Vpriorinv)),silent=TRUE)
    # if (is(V_post,"try-error")) V_post <- ginv(crossprod(Zt)+Vpriorinv)
    # alternative a la supplementary bitto/sfs s.3
    Vpriorsqrt <- diag(c(sqrt(tau2_draw),sqrt(xi2_draw)))
    V_poststar <- solve(Vpriorsqrt%*%crossprod(Zt)%*%Vpriorsqrt + diag(2*d))
    V_post <- Vpriorsqrt%*%V_poststar%*%Vpriorsqrt

    A_post <- V_post%*%(crossprod(Zt,yt)+Vpriorinv%*%A_prior)
    alph_draw <- try(A_post+t(chol(V_post))%*%rnorm(ncol(Zt)),silent=TRUE)
    if (is(alph_draw,"try-error")) alph_draw <- matrix(mvrnorm(1,A_post,V_post),ncol(Zt),1)

    Am_draw    <- alph_draw[1:d,,drop=FALSE]
    theta_sqrt <- alph_draw[(d+1):(2*d),,drop=TRUE]
    theta_draw <- theta_sqrt^2
    #----------------------------------------------------------------------------
    # Step 2: Sample TVP-coef
    ystar <- yt - Xt%*%Am_draw
    Fstar <- Xt%*%diag(theta_sqrt)

    At_draw <- sample_McCausland(ystar, Fstar)
    #----------------------------------------------------------------------------
    # Step 3: Interweaving
    theta_sign <- sign(theta_sqrt)
    A_draw     <- matrix(Am_draw,bigT+1,d,byrow=TRUE) + At_draw%*%diag(theta_sqrt)
    A_diff     <- diff(At_draw%*%diag(theta_sqrt))
    #A_diff     <- diff(A_draw) # same as line above
    for(dd in 1:d){
      # theta.new
      res <- do_rgig1(lambda=-bigT/2,
                      chi=sum(A_diff[,dd]^2)+(A_draw[1,dd]-Am_draw[dd,1])^2,
                      psi=1/xi2_draw[dd])
      theta_draw[dd] <- res
      theta_sqrt[dd] <- sqrt(res)*theta_sign[dd]
      # betam.new
      sigma2_A_mean <- 1/((1/tau2_draw[dd]) + (1/theta_draw[dd]))
      mu_A_mean     <- A_draw[1,dd]*tau2_draw[dd]/(tau2_draw[dd] + theta_draw[dd])
      Am_draw[dd,1] <- rnorm(1, mu_A_mean, sqrt(sigma2_A_mean))
    }
    At_draw <- sapply(1:d,function(dd)A_draw[,dd]-Am_draw[dd,1])%*%diag(1/theta_sqrt)
    #----------------------------------------------------------------------------
    # Step 4: Prior choice
    if(prior==1){ # TVP
      ### no hierarchical priors
    }else if(prior==2){ # TVP-NG
      kappa2  <- rgamma(1, d1+a_xi*d,  d2+0.5*a_xi*mean(xi2_draw)*d)
      lambda2 <- rgamma(1, e1+a_tau*d, e2+0.5*a_tau*mean(tau2_draw)*d)
      for(dd in 1:d){
        xi2_draw[dd]  <- do_rgig1(lambda=a_xi-0.5,  chi=theta_draw[dd], psi=a_xi*kappa2)
        tau2_draw[dd] <- do_rgig1(lambda=a_tau-0.5, chi=(Am_draw[dd,1]-A_prior[dd,1])^2, psi=a_tau*lambda2)
      }
      xi2_draw[xi2_draw<1e-7] <- 1e-7
      tau2_draw[tau2_draw<1e-7] <- 1e-7
      if(sample_A){
        before <- a_xi
        a_xi   <- MH_step(a_xi, scale_xi, d, kappa2, theta_sqrt, b_xi, nu_xi, d1, d2)
        if(before!=a_xi){
          acc_xi <- acc_xi + 1
        }
        before <- a_tau
        a_tau  <- MH_step(a_tau, scale_xi, d, lambda2, Am_draw, b_tau, nu_tau, e1, e2)
        if(before!=a_tau){
          acc_tau <- acc_tau + 1
        }
        # scale MH proposal during the first 50% of the burn-in stage
        if(irep<(0.5*burnin)){
          if((acc_xi/irep)>0.30){scale_xi <- 1.01*scale_xi}
          if((acc_xi/irep)<0.15){scale_xi <- 0.99*scale_xi}
          if((acc_tau/irep)>0.30){scale_xi <- 1.01*scale_xi}
          if((acc_tau/irep)<0.15){scale_xi <- 0.99*scale_xi}
        }
      }
    } # END PRIOR QUERY
    #----------------------------------------------------------------------------
    # Step 5: Sample variances
    eps <- y - cbind(Xt,hadamard.prod(Xt,At_draw[2:(bigT+1),]))%*%alph_draw
    if(sv){
      para <- as.list(pars_var); names(para) <- c("mu","phi","sigma","latent0")
      para$nu = Inf; para$rho=0; para$beta<-0
      svdraw <- svsample_fast_cpp(y=eps, draws=1, burnin=0, designmatrix=matrix(NA_real_),
                                  priorspec=Sv_priors, thinpara=1, thinlatent=1, keeptime="all",
                                  startpara=para, startlatent=Sv_draw,
                                  keeptau=FALSE, print_settings=list(quiet=TRUE, n_chains=1, chain=1),
                                  correct_model_misspecification=FALSE, interweave=TRUE, myoffset=0,
                                  fast_sv=default_fast_sv)
      h_           <- exp(svdraw$latent[1,])
      para$mu      <- svdraw$para[1,"mu"]
      para$phi     <- svdraw$para[1,"phi"]
      para$sigma   <- svdraw$para[1,"sigma"]
      para$latent0 <- svdraw$latent0[1,"h_0"]
      pars_var     <- unlist(para[c("mu","phi","sigma","latent0")])
      Sv_draw       <- log(h_)
    }else{
      C0  <- rgamma(1, g0+c0, G0+sig_eta)
      S_1 <- c0+bigT/2
      S_2 <- C0+crossprod(eps)/2

      sig_eta <- 1/rgamma(1,S_1,S_2)
      Sv_draw <- matrix(log(sig_eta),bigT,1)
    }
    #-------------------------------------------------------------------------#
    # STEP 6: RANDOM SIGN SWITCH
    for(dd in 1:d){
      if(runif(1,0,1)>0.5){
        theta_sqrt[dd] <- -theta_sqrt[dd]
      }
    }
    #----------------------------------------------------------------------------
    # Step 7: store draws
    if(irep %in% thin.draws){
      count <- count+1
      A_store[count,,]<- A_draw
      res_store[count,,]<- eps
      # SV
      Sv_store[count,,] <- Sv_draw
      pars_store[count,,] <- pars_var
      # NG
      thetasqrt_store[count,,] <- theta_sqrt
      tau2_store[count,,]<- tau2_draw
      xi2_store[count,,] <- xi2_draw
      lambda2_store[count,,] <- lambda2
      kappa2_store[count,,] <- kappa2
      a_xi_store[count,,] <- a_xi
      a_tau_store[count,,]<- a_tau
    }
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,paste("t",seq(0,bigT),sep="."),colnames(X))
  ret <- list(Y=y,X=X,A_store=A_store,Sv_store=Sv_store,pars_store=pars_store,res_store=res_store,
              thetasqrt_store=thetasqrt_store,tau2_store=tau2_store,xi2_store=xi2_store,lambda2_store=lambda2_store,kappa2_store=kappa2_store,a_xi_store=a_xi_store,a_tau_store=a_tau_store)
  return(ret)
}

#' @name .TVPBVAR_centered_R.m
#' @importFrom stochvol svsample_fast_cpp specify_priors default_fast_sv sv_normal sv_beta sv_gamma
#' @importFrom dlm dlmModReg dlmMLE dlmSmooth
#' @importFrom MASS ginv mvrnorm
#' @importFrom matrixcalc hadamard.prod
#' @importFrom methods is
#' @importFrom stats rnorm rgamma runif dnorm
#' @noRd
.TVPBVAR_centered_R <- function(nr,Y_in,p_in,draws_in,burnin_in,cons_in,trend_in,sv_in,thin_in,quiet_in,prior_in,hyperparam_in,Ex_in){
  #----------------------------------------INPUTS----------------------------------------------------#
  Yraw  <- Y_in
  p     <- p_in
  Traw  <- nrow(Yraw)
  M     <- ncol(Yraw)
  K     <- M*p
  Ylag  <- .mlag(Yraw,p)
  names <- colnames(Yraw)
  if(is.null(names)) names <- rep("Y",M)
  colnames(Yraw) <- names
  nameslags <- NULL
  for(ii in 1:p) nameslags <- c(nameslags,paste0(names,".lag",ii))
  colnames(Ylag) <- nameslags

  texo <- FALSE; Mex <- 0; Exraw <- NULL; enames <- NULL
  if(!is.null(Ex_in)){
    Exraw <- Ex_in; Mex <- ncol(Exraw); texo <- TRUE
    enames <- colnames(Exraw)
    if(is.null(enames)) enames <- rep("Tex",Mex)
    colnames(Exraw) <- enames
  }

  if(nr==1) slct <- NULL else slct <- 1:(nr-1)

  Xraw  <- cbind(Yraw[,slct],Ylag,Exraw)
  colnames(Xraw) <- c(colnames(Yraw)[slct],nameslags,enames)
  X     <- Xraw[(p+1):nrow(Xraw),,drop=FALSE]
  y     <- Yraw[(p+1):Traw,nr,drop=FALSE]
  bigT  <- nrow(X)
  M_    <- M-length(slct)

  cons  <- cons_in
  if(cons){
    X <- cbind(X,1)
    colnames(X)[ncol(X)] <- "cons"
  }
  trend <- trend_in
  if(trend){
    X <- cbind(X,seq(1,bigT))
    colnames(X)[ncol(X)] <- "trend"
  }

  d <- ncol(X)
  n <- d*M
  v <- (M*(M-1))/2
  #---------------------------------------------------------------------------------------------------------
  # HYPERPARAMETERS
  #---------------------------------------------------------------------------------------------------------
  hyperpara      <- hyperparam_in
  prior          <- prior_in
  sv             <- sv_in
  prmean         <- hyperpara$prmean
  # non-SV
  c0             <- hyperpara$c0
  g0             <- hyperpara$g0
  # SV
  bmu            <- hyperpara$bmu
  Bmu            <- hyperpara$Bmu
  a0             <- hyperpara$a0
  b0             <- hyperpara$b0
  Bsigma         <- hyperpara$Bsigma
  # TTVP
  B_1            <- hyperpara$B_1
  B_2            <- hyperpara$B_2
  kappa0         <- hyperpara$kappa0
  a_tau          <- hyperpara$a_tau
  c_tau          <- hyperpara$c_tau
  d_tau          <- hyperpara$d_tau
  h0prior        <- hyperpara$h0prior
  grid.length    <- hyperpara$grid.length
  thrsh.pct      <- hyperpara$thrsh.pct
  thrsh.pct.high <- hyperpara$thres.pct.high
  TVS            <- hyperpara$TVS
  a.approx       <- hyperpara$a.approx
  sim.kappa      <- hyperpara$sim.kappa
  kappa.grid     <- hyperpara$kappa.grid
  MaxTrys        <- hyperpara$MaxTrys
  #---------------------------------------------------------------------------------------------------------
  # OLS Quantitites
  #---------------------------------------------------------------------------------------------------------
  XtXinv <- try(solve(crossprod(X)),silent=TRUE)
  if(is(XtXinv,"try-error")) XtXinv <- ginv(crossprod(X))
  A_OLS  <- XtXinv%*%(t(X)%*%y)
  E_OLS  <- y - X%*%A_OLS
  S_OLS  <- crossprod(E_OLS)/(bigT-d)
  V_OLS  <- as.numeric(S_OLS)*XtXinv
  sd_OLS <- sqrt(diag(V_OLS))
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  A_draw  <- matrix(A_OLS, bigT+1, d, byrow=TRUE, dimnames=list(NULL,colnames(X)))
  S_draw  <- matrix(S_OLS, bigT,1)

  # state variances
  Omega_t <- matrix(1,bigT,d)
  # state indicator
  D_t <- matrix(1,bigT,d)
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  # Priors on VAR coefs
  #-----------------------------
  # prior mean
  A_prior <- matrix(0,2*d, 1)
  A_prior[2*nr-1,1] <- prmean
  # prior variance
  sqrttheta1 <- diag(d)*0.1
  sqrttheta2 <-diag(d)*0.01
  Omega_t <- D_t%*%sqrttheta1+(1-D_t)%*%sqrttheta2
  kappa00 <- kappa0
  if(kappa0<0) kappa00 <- -kappa0 * sd.OLS else kappa00 <- matrix(kappa0,d,1)

  if (a.approx){
    buildCapm <- function(u){
      dlm::dlmModReg(X, dV = exp(u[1]), dW = exp(u[2:(d+1)]),addInt = FALSE)
    }
    outMLE <- dlm::dlmMLE(y, parm = rep(0,d+1), buildCapm)
    mod <- buildCapm(outMLE$par)
    outS <- dlm::dlmSmooth(y, mod)
    states.OLS <- t(matrix(outS$s,bigT+1,d))
    Achg.OLS <- t(diff(t(states.OLS)))#t(as.numeric(ALPHA0)+ALPHA2)
  }
  # threshold
  thrsh <- matrix(0,d,1)

  # priors on initial state
  B0prior <- matrix(0,d,1)
  V0prior <- rep(4,d)
  #------------------------------------
  # SV quantities
  #------------------------------------
  svdraw  <- list(para=c(mu=-10,phi=.9,sigma=.2,latent0=-3),latent=rep(-3,bigT))
  Sv_draw <- svdraw$latent
  pars_var <- matrix(c(-3,.9,.2,-3),4,1,dimnames=list(c("mu","phi","sigma","latent0"),NULL))
  Sv_priors <- specify_priors(mu=sv_normal(mean=bmu, sd=Bmu), phi=sv_beta(a0,b0), sigma2=sv_gamma(shape=0.5,rate=1/(2*Bsigma)))
  #-----------------------------------
  # non-SV quantities
  #-----------------------------------
  sig_eta <- exp(-3)
  G0 <- g0/S_OLS*(c0-1)
  #---------------------------------------------------------------------------------------------------------
  # SAMPLER MISCELLANEOUS
  #---------------------------------------------------------------------------------------------------------
  nsave <- draws_in
  nburn <- burnin_in
  ntot  <- nsave+nburn

  # thinning
  thin         <- thin_in
  count <- 0
  thindraws    <- nsave/thin
  thin.draws   <- seq(nburn+1,ntot,by=thin)
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store      <- array(NA,c(thindraws,bigT,d))
  res_store    <- array(NA,c(thindraws,bigT,1))
  Sv_store     <- array(NA,c(thindraws,bigT,1))
  pars_store   <- array(NA,c(thindraws,4,1))
  # TTVP
  D_store      <- array(NA,c(thindraws,bigT,d,1))
  Omega_store  <- array(NA,c(thindraws,bigT,d,1))
  thrsh_store  <- array(NA,c(thindraws,d,1))
  kappa_store  <- array(NA,c(thindraws,1))
  V0_store     <- array(NA,c(thindraws,d,1))
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Draw A
    invisible(capture.output(
        A_draw1 <- try(KF_fast(t(as.matrix(y)), X,as.matrix(exp(Sv_draw)),Omega_t,d, 1, bigT, B0prior, diag(V0prior)),silent=TRUE),
    type="message"))
    if (is(A_draw1,"try-error")){
      invisible(capture.output(
        A_draw1 <- KF(t(as.matrix(y)), X,as.matrix(exp(Sv_draw)),Omega_t,d, 1, bigT, B0prior, diag(V0prior)),
      type="message"))
      try0 <- 0
      while (any(abs(A_draw1$bdraw)>1e+10) && try0<MaxTrys){ #This block resamples if the draw from the state vector is not well behaved
        invisible(capture.output(
          A_draw1 <- try(KF(t(as.matrix(y)), X,as.matrix(exp(Sv_draw)),Omega_t,d, 1, bigT, B0prior, diag(V0prior)),silent=TRUE),
        type="message"))
        try0 <- try0+1
      }
    }
    A_draw <- t(A_draw1$bdraw)
    VCOV   <- A_draw1$Vcov
    #----------------------------------------------------------------------------
    # Step 2: Prior choice
    if(prior==3){
      #------------------------------------------
      # Step 2a: Sample variances
      A_diff <- diff(A_draw)
      for(dd in 1:d){
        sig_q <- sqrttheta1[dd,dd]

        if (!a.approx){
          si <- (abs(A_diff[,dd])>thrsh[dd,1])*1
        }else{
          si <- (abs(Achg.OLS[,dd])>thrsh[dd,1])*1
        }

        si <- D_t[2:bigT,dd]
        s_1 <- B_1 + sum(si)/2 + 0.5
        s_2 <- B_2 + 0.5*crossprod(A_diff[si==1,dd,drop=FALSE])
        sig_q <- 1/rgamma(1,s_1,s_2)

        sqrttheta1[dd,dd] <- sig_q
        sqrttheta2[dd,dd] <- kappa00[dd,1]^2
      }
      #------------------------------------------
      # sample indicator
      if(TVS){
        #Check whether coefficient is time-varying or constant at each point in time
        Achg <- t(A_diff)
        Achg <- cbind(matrix(0,d,1),Achg) #we simply assume that the parameters stayed constant between t=0 and t=1
        if(a.approx) Achg.approx <- Achg.OLS else Achg.approx <- Achg
        grid.mat <- matrix(unlist(lapply(1:d,function(x) .get_grid(Achg[x,],sqrt(sqrttheta1[x,x]),grid.length=grid.length,thrsh.pct=thrsh.pct,thrsh.pct.high=thrsh.pct.high))),ncol = d)
        probs    <- get_threshold(Achg, sqrttheta1, sqrttheta2, grid.mat, Achg.approx)
        for(dd in 1:d){
          post1 <- probs[,dd]
          probs1 <- exp(post1-max(post1))/sum(exp(post1-max(post1)))
          thrsh[dd,] <- sample(grid.mat[,dd],1,prob=probs1)
          if (!a.approx){
            D_t[,dd] <- (abs(Achg[dd,])>thrsh[dd,])*1 #change 2:T usw. here
          }else{
            D_t[,dd] <- (abs(Achg.OLS[dd,])>thrsh[dd,])*1 #change 2:T usw. here
          }
        }
        if (sim.kappa){
          grid.kappa <- kappa.grid
          Lik.kappa <- matrix(0,length(grid.kappa),1)
          count <- 0
          for (grid.i in grid.kappa){
            count <- count+1
            sqrttheta.prop <- (grid.i*sd_OLS)^2
            cov.prop <- sqrt(D_t*diag(sqrttheta1)+(1-D_t)*sqrttheta.prop)

            Lik.kappa[count,1] <-sum(dnorm(t(Achg),matrix(0,bigT,2),cov.prop,log=TRUE))
          }
          Lik.kappa.norm <- exp(Lik.kappa-max(Lik.kappa))
          probs.kappa <- Lik.kappa.norm/sum(Lik.kappa.norm)
          kappa0 <- sample(grid.kappa,size=1, prob=probs.kappa)
          kappa00 <- kappa0*sd_OLS
        }
      }
      Omega_t <- D_t%*%sqrttheta1+(1-D_t)%*%sqrttheta2
      #------------------------------------------
      # Step 2b: Draw variance of initial state
      lambda2_tau <- rgamma(1,c_tau+a_tau*d,d_tau+a_tau/2*sum(V0prior)) # global component
      # local component
      for(dd in 1:d){
        res <- try(do_rgig1(lambda=a_tau-0.5,
                            chi=A_draw[1,dd]^2,
                            psi=a_tau*lambda2_tau), silent=TRUE)
        V0prior[dd] <- ifelse(is(res,"try-error"),next,res)
      }
    } # END PRIOR QUERY
    #----------------------------------------------------------------------------
    # Step 3: Sample variances
    eps <- y - rowSums(hadamard.prod(X,A_draw))
    if(sv){
      para <- as.list(pars_var); names(para) <- c("mu","phi","sigma","latent0")
      para$nu = Inf; para$rho=0; para$beta<-0
      svdraw <- svsample_fast_cpp(y=eps, draws=1, burnin=0, designmatrix=matrix(NA_real_),
                                  priorspec=Sv_priors, thinpara=1, thinlatent=1, keeptime="all",
                                  startpara=para, startlatent=Sv_draw,
                                  keeptau=FALSE, print_settings=list(quiet=TRUE, n_chains=1, chain=1),
                                  correct_model_misspecification=FALSE, interweave=TRUE, myoffset=0,
                                  fast_sv=default_fast_sv)
      h_           <- exp(svdraw$latent[1,])
      para$mu      <- svdraw$para[1,"mu"]
      para$phi     <- svdraw$para[1,"phi"]
      para$sigma   <- svdraw$para[1,"sigma"]
      para$latent0 <- svdraw$latent0[1,"h_0"]
      pars_var     <- unlist(para[c("mu","phi","sigma","latent0")])
      Sv_draw       <- log(h_)
    }else{
      C0  <- rgamma(1, g0+c0, G0+sig_eta)
      S_1 <- c0+bigT/2
      S_2 <- C0+crossprod(eps)/2

      sig_eta <- 1/rgamma(1,S_1,S_2)
      Sv_draw <- matrix(log(sig_eta),bigT,1)
    }
    #----------------------------------------------------------------------------
    # Step 4: store draws
    if(irep %in% thin.draws){
      count <- count+1
      A_store[count,,]<- A_draw
      res_store[count,,]<- eps
      # SV
      Sv_store[count,,] <- Sv_draw
      pars_store[count,,] <- pars_var
      # TTVP
      D_store[count,,,] <- D_t
      Omega_store[count,,,] <- Omega_t
      thrsh_store[count,,] <- thrsh
      kappa_store[count,] <- kappa0
      V0_store[count,,] <- V0prior
    }
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,paste("t",seq(1,bigT),sep="."),colnames(X))
  ret <- list(Y=y,X=X,A_store=A_store,Sv_store=Sv_store,pars_store=pars_store,res_store=res_store,
              D_store=D_store,Omega_store=Omega_store,thrsh_store=thrsh_store,kappa_store=kappa_store,V0_store=V0_store)
  return(ret)
}

#' @name .get_grid
#' @noRd
.get_grid <- function(Achg,sd.state,grid.length=150,thrsh.pct=0.1,thrsh.pct.high=0.9){
  d_prop <- seq(thrsh.pct*sd.state,thrsh.pct.high*sd.state,length.out=grid.length)
  return(d_prop)
}

#' @name .gck
#' @noRd
.gck <- function(yg,gg,hh,capg,f,capf,sigv,kold,t,ex0,vx0,nvalk,kprior,kvals,p,kstate){
  # GCK's Step 1 on page 821
  lpy2n=0;
  mu = matrix(0,t*kstate,1);
  omega = matrix(0,t*kstate,kstate);
  for (i in seq(t-1,1,by=-1)){
    gatplus1 = sigv%*%kold[i+1]
    ftplus1 = capf[(kstate*i+1):(kstate*(i+1)),]
    cgtplus1 = capg[(i*p+1):((i+1)*p),]
    htplus1 = t(hh[(i*p+1):((i+1)*p),])

    htt1 <- crossprod(htplus1,gatplus1)
    rtplus1 = tcrossprod(htt1)+tcrossprod(cgtplus1,cgtplus1)
    rtinv = solve(rtplus1)
    btplus1 = tcrossprod(gatplus1)%*%htplus1%*%rtinv
    atplus1 = (diag(kstate)-tcrossprod(btplus1,htplus1))%*%ftplus1

    if (kold[i+1] == 0){
      ctplus1 = matrix(0,kstate,kstate)
    }else{
      cct = gatplus1%*%(diag(kstate)-crossprod(gatplus1,htplus1)%*%tcrossprod(rtinv,htplus1)%*%gatplus1)%*%t(gatplus1)
      ctplus1 = t(chol(cct))
    }
    otplus1 = omega[(kstate*i+1):(kstate*(i+1)),]

    dtplus1 = crossprod(ctplus1,otplus1)%*%ctplus1+diag(kstate)
    omega[(kstate*(i-1)+1):(kstate*i),] = crossprod(atplus1,(otplus1 - otplus1%*%ctplus1%*%solve(dtplus1)%*%t(ctplus1)%*%otplus1))%*%atplus1+t(ftplus1)%*%htplus1%*%rtinv%*%t(htplus1)%*%ftplus1
    satplus1 = (diag(kstate)-tcrossprod(btplus1,htplus1))%*%(f[,i+1]-btplus1%*%gg[,i+1]) #CHCKCHCKCHKC
    mutplus1 = mu[(kstate*i+1):(kstate*(i+1)),]
    mu[(kstate*(i-1)+1):(kstate*i),] = crossprod(atplus1,(diag(kstate)-otplus1%*%ctplus1%*%solve(dtplus1)%*%t(ctplus1)))%*%(mutplus1-otplus1%*%(satplus1+btplus1%*%yg[i+1]))+t(ftplus1)%*%htplus1%*%rtinv%*%(yg[i+1]-gg[,i+1]-t(htplus1)%*%f[,i+1])
  }

  # GCKs Step 2 on pages 821-822
  kdraw = kold;
  ht = t(hh[1:p,])
  ft = capf[1:kstate,]
  gat = matrix(0,kstate,kstate)
  # Note: this specification implies no shift in first period -- sensible
  rt = t(ht)%*%ft%*%vx0%*%t(ft)%*%ht + crossprod(ht,gat)%*%crossprod(gat,ht)+ tcrossprod(capg[1:p,])
  rtinv = solve(rt)
  jt = (ft%*%vx0%*%t(ft)%*%ht + tcrossprod(gat)%*%ht)%*%rtinv
  mtm1 = (diag(kstate) - tcrossprod(jt,ht))%*%(f[,1] + ft%*%ex0) + jt%*%(yg[1] - gg[,1])
  vtm1 <- ft%*%tcrossprod(vx0,ft)+tcrossprod(gat)-jt%*%tcrossprod(rt,jt)
  lprob <- matrix(0,nvalk,1)

  for (i in 2:t){
    ht <-  t(hh[((i-1)*p+1):(i*p),])
    ft <- capf[(kstate*(i-1)+1):(kstate*i),]
    for (j in 1:nvalk){
      gat <- kvals[j,1]%*%sigv
      rt <- crossprod(ht,ft)%*%tcrossprod(vtm1,ft)%*%ht+crossprod(ht,gat)%*%crossprod(gat,ht)+tcrossprod(capg[((i-1)*p+1):(i*p),])
      rtinv <- solve(rt)
      jt <- (ft%*%tcrossprod(vtm1,ft)%*%ht+tcrossprod(gat)%*%ht)%*%rtinv
      mt <- (diag(kstate)-tcrossprod(jt,ht))%*%(f[,i]+ft%*%mtm1)+jt%*%(yg[i]-gg[,i])
      vt <- ft%*%tcrossprod(vtm1,ft)+tcrossprod(gat)-jt%*%tcrossprod(rt,jt)

      lpyt = -.5*log(det(rt)) - .5*t(yg[i] - gg[,i] - t(ht)%*%t(f[,i] + ft%*%mtm1))%*%rtinv%*%(yg[i] - gg[,i] - t(ht)%*%(f[,i] + ft%*%mtm1))

      if (det(vt)<=0){
        tt <- matrix(0,kstate,kstate)
      }else{
        tt <- t(chol(vt))
      }
      ot = omega[(kstate*(i-1)+1):(kstate*i),]
      mut = mu[(kstate*(i-1)+1):(kstate*i),]
      tempv = diag(kstate) + crossprod(tt,ot)%*%tt
      lpyt1n = -.5*log(det(tempv)) -.5*(crossprod(mt,ot)%*%mt-2*crossprod(mut,mt)-t(mut-ot%*%mt)%*%tt%*%solve(tempv)%*%t(tt)%*%(mut-ot%*%mt))
      lprob[j,1] <- log(kprior[j,1])+lpyt1n+lpyt
      if (i==2){
        lpy2n <- lpyt1n+lpyt
      }
    }
    pprob = exp(lprob-max(lprob))/sum(exp(lprob-max(lprob)))
    tempv = runif(1)
    tempu = 0
    for (j in 1:nvalk){
      tempu <- tempu+pprob[j,1]
      if (tempu> tempv){
        kdraw[i] <- kvals[j,1]
        break
      }
    }
    gat = kdraw[i]%*%sigv
    rt = crossprod(ht,ft)%*%tcrossprod(vtm1,ft)%*%ht+t(ht)%*%tcrossprod(gat)%*%ht+tcrossprod(capg[((i-1)*p+1):(i*p)])
    rtinv = solve(rt)
    jt = (ft%*%tcrossprod(vtm1,ft)%*%ht+tcrossprod(gat)%*%ht)%*%rtinv
    mtm1 <- (diag(kstate)-tcrossprod(jt,ht))%*%(f[,i]+ft%*%mtm1)+jt%*%(yg[i]-gg[,i])
    vtm1 = ft%*%tcrossprod(vtm1,ft)+tcrossprod(gat)-jt%*%tcrossprod(rt,jt)
  }
  return(kdraw)
}

#' @name .var_posterior
#' @importFrom MASS ginv
#' @importFrom abind adrop abind
#' @noRd
.var_posterior <- function(post_draws, prior, draws, applyfun, cores){
  M <- length(post_draws)
  bigT <- nrow(post_draws[[1]]$Y)
  bigK <- ncol(post_draws[[1]]$X)
  K    <- unlist(lapply(post_draws,function(l)ncol(l$X)))
  # bind data
  Y <- do.call("cbind",lapply(1:M,function(mm)post_draws[[mm]]$Y))
  X <- post_draws[[1]]$X
  # general stuff
  res_store  <- abind(lapply(1:M,function(mm)post_draws[[mm]]$res_store),along=3)
  Sv_store   <- abind(lapply(1:M,function(mm)post_draws[[mm]]$Sv_store),along=3)
  pars_store <- abind(lapply(1:M,function(mm)post_draws[[mm]]$pars_store),along=3)
  # container
  A_store        <- array(NA,c(draws,bigT,bigK,M))
  L_store        <- array(NA,c(draws,bigT,M,M))
  S_store        <- array(NA,c(draws,bigT,M,M))
  timepoints     <- paste("t",seq(1,bigT),sep=".")
  store.obj <- applyfun(1:draws,function(irep){
    At_store <- array(NA,c(bigT,bigK,M))
    Lt_store <- array(NA,c(bigT,M,M))
    St_store <- array(NA,c(bigT,M,M))
    for(tt in 1:bigT){
      A0 <- diag(M)
      for(mm in 2:M){
        A0[mm,1:(mm-1)] <- -post_draws[[mm]]$A_store[irep,timepoints[tt],1:(mm-1)]
      }
      A0inv <- try(solve(A0),silent=TRUE)
      if(is(A0inv,"try-error")) A0inv <- ginv(A0)
      Lt_store[tt,,] <- A0inv
      St_store[tt,,] <- A0inv%*%diag(exp(Sv_store[irep,tt,]))%*%t(A0inv)
      Atilde <- NULL
      for(mm in 1:M) Atilde <- cbind(Atilde,post_draws[[mm]]$A_store[irep,timepoints[tt],mm:K[mm]])
      At_store[tt,,] <- t(A0inv%*%t(Atilde))
    }
    return(list(At_store=At_store,Lt_store=Lt_store,St_store=St_store))
  })
  for(irep in 1:draws){
    A_store[irep,,,] <- store.obj[[irep]]$At_store
    L_store[irep,,,] <- store.obj[[irep]]$Lt_store
    S_store[irep,,,] <- store.obj[[irep]]$St_store
  }
  dimnames(A_store) <- list(NULL,timepoints,colnames(X),colnames(Y))
  Smed_store <- apply(S_store,c(1,3,4),median)
  if(prior=="TVP"){
    thetasqrt_store <- abind(lapply(1:M,function(mm)post_draws[[mm]]$thetasqrt_store[,mm:K[mm],]),along=3)
    Lthetasqrt_store <- lapply(2:M,function(mm)adrop(post_draws[[mm]]$thetasqrt_store[,1:(mm-1),,drop=FALSE],drop=3))
    tau2_store<-xi2_store<-Ltau2_store<-Lxi2_store<-lambda2_store<-kappa2_store<-a_xi_store<-a_tau_store<-D_store<-Omega_store<-thrsh_store<-kappa_store<-V0_store<-LD_store<-LOmega_store<-Lthrsh_store<-LV0_store<-NULL
  }else if(prior=="TVP-NG"){
    D_store<-Omega_store<-thrsh_store<-kappa_store<-V0_store<-LD_store<-LOmega_store<-Lthrsh_store<-LV0_store<-NULL
    # general stuff
    thetasqrt_store <- abind(lapply(1:M,function(mm)post_draws[[mm]]$thetasqrt_store[,mm:K[mm],]),along=3)
    lambda2_store <- abind(lapply(1:M,function(mm)post_draws[[mm]]$lambda2_store),along=3)
    kappa2_store  <- abind(lapply(1:M,function(mm)post_draws[[mm]]$kappa2_store),along=3)
    a_xi_store    <- abind(lapply(1:M,function(mm)post_draws[[mm]]$a_xi_store),along=3)
    a_tau_store   <- abind(lapply(1:M,function(mm)post_draws[[mm]]$a_tau_store),along=3)
    tau2_store    <- abind(lapply(1:M,function(mm)post_draws[[mm]]$tau2_store[,mm:K[mm],]),along=3)
    xi2_store     <- abind(lapply(1:M,function(mm)post_draws[[mm]]$xi2_store[,mm:K[mm],]),along=3)
    ## ATTENTION: variances of L just as list !!
    Lthetasqrt_store <- lapply(2:M,function(mm)adrop(post_draws[[mm]]$thetasqrt_store[,1:(mm-1),,drop=FALSE],drop=3))
    Ltau2_store   <- lapply(2:M,function(mm)adrop(post_draws[[mm]]$tau2_store[,1:(mm-1),,drop=FALSE],drop=3))
    Lxi2_store    <- lapply(2:M,function(mm)adrop(post_draws[[mm]]$xi2_store[,1:(mm-1),,drop=FALSE],drop=3))
  }else if(prior=="TTVP"){
    thetasqrt_store<-Lthetasqrt_store<-tau2_store<-xi2_store<-Ltau2_store<-Lxi2_store<-lambda2_store<-kappa2_store<-a_xi_store<-a_tau_store<-NULL
    # general stuff
    kappa_store <- abind(lapply(1:M,function(mm)post_draws[[mm]]$kappa_store),along=3)
    D_store     <- abind(lapply(1:M,function(mm)post_draws[[mm]]$D_store[,,mm:K[mm],]),along=4)
    Omega_store <- abind(lapply(1:M,function(mm)post_draws[[mm]]$Omega_store[,,mm:K[mm],]),along=4)
    thrsh_store <- abind(lapply(1:M,function(mm)post_draws[[mm]]$thrsh_store[,mm:K[mm],]),along=3)
    V0_store    <- abind(lapply(1:M,function(mm)post_draws[[mm]]$V0_store[,mm:K[mm],]),along=3)
    ## ATTENTION: variances of L just as list !!
    LD_store    <- lapply(2:M,function(mm)adrop(post_draws[[mm]]$D_store[,,1:(mm-1),,drop=FALSE],drop=4))
    LOmega_store<- lapply(2:M,function(mm)adrop(post_draws[[mm]]$Omega_store[,,1:(mm-1),,drop=FALSE],drop=4))
    Lthrsh_store<- lapply(2:M,function(mm)adrop(post_draws[[mm]]$thrsh_store[,1:(mm-1),,drop=FALSE],drop=3))
    LV0_store   <- lapply(2:M,function(mm)adrop(post_draws[[mm]]$V0_store[,1:(mm-1),,drop=FALSE],drop=3))
  }
  ret <- list(Y=Y,X=X,A_store=A_store,L_store=L_store,Sv_store=Sv_store,S_store=S_store,Smed_store=Smed_store,pars_store=pars_store,res_store=res_store,thetasqrt_store=thetasqrt_store,Lthetasqrt_store=Lthetasqrt_store,
              tau2_store=tau2_store,xi2_store=xi2_store,Ltau2_store=Ltau2_store,Lxi2_store=Lxi2_store,lambda2_store=lambda2_store,kappa2_store=kappa2_store,a_xi_store=a_xi_store,a_tau_store=a_tau_store,
              D_store=D_store,Omega_store=Omega_store,thrsh_store=thrsh_store,kappa_store=kappa_store,V0_store=V0_store,LD_store=LD_store,LOmega_store=LOmega_store,Lthrsh_store=Lthrsh_store,LV0_store=LV0_store)
  return(ret)
}

#' @name .TVPBVAR_linear_R
#' @importFrom stochvol svsample_fast_cpp specify_priors default_fast_sv sv_normal sv_beta sv_gamma
#' @importFrom MASS ginv mvrnorm
#' @importFrom matrixcalc hadamard.prod
#' @importFrom methods is
#' @importFrom stats rnorm rgamma runif dnorm
#' @noRd
.TVPBVAR_linear_R <- function(Y_in,p_in,draws_in,burnin_in,cons_in,trend_in,sv_in,thin_in,quiet_in,prior_in,hyperparam_in,Ex_in){
  #----------------------------------------INPUTS----------------------------------------------------#
  Yraw  <- Y_in
  p     <- p_in
  Traw  <- nrow(Yraw)
  M     <- ncol(Yraw)
  K     <- M*p
  Ylag  <- .mlag(Yraw,p)
  nameslags <- NULL
  for (ii in 1:p) nameslags <- c(nameslags,rep(paste("Ylag",ii,sep=""),M))
  colnames(Ylag) <- nameslags

  texo <- FALSE; Mex <- 0; Exraw <- NULL
  if(!is.null(Ex_in)){
    Exraw <- Ex_in; Mex <- ncol(Exraw)
    texo <- TRUE
    colnames(Exraw) <- rep("Tex",Mex)
  }

  X <- cbind(Ylag,Exraw)
  X <- X[(p+1):nrow(X),,drop=FALSE]
  Y <- Yraw[(p+1):Traw,,drop=FALSE]
  bigT  <- nrow(X)

  cons  <- cons_in
  if(cons){
    X <- cbind(X,1)
    colnames(X)[ncol(X)] <- "cons"
  }
  trend <- trend_in
  if(trend){
    X <- cbind(X,seq(1,bigT))
    colnames(X)[ncol(X)] <- "trend"
  }

  k     <- ncol(X)
  n <- k*M
  v <- (M*(M-1))/2
  #---------------------------------------------------------------------------------------------------------
  # HYPERPARAMETERS
  #---------------------------------------------------------------------------------------------------------
  hyperpara <- hyperparam_in
  prior     <- prior_in
  sv        <- sv_in
  prmean    <- hyperpara$prmean
  a_1       <- hyperpara$a_1
  b_1       <- hyperpara$b_1
  # SV
  Bsigma    <- hyperpara$Bsigma
  a0        <- hyperpara$a0
  b0        <- hyperpara$b0
  bmu       <- hyperpara$bmu
  Bmu       <- hyperpara$Bmu
  # other stuff
  d1        <- hyperpara$d1
  d2        <- hyperpara$d2
  e1        <- hyperpara$e1
  e2        <- hyperpara$e2
  b_xi      <- hyperpara$b_xi
  b_tau     <- hyperpara$b_tau
  nu_xi     <- hyperpara$nu_xi
  nu_tau    <- hyperpara$nu_tau
  a_start   <- hyperpara$a_start
  sample_A  <- hyperpara$sample_A
  #---------------------------------------------------------------------------------------------------------
  # OLS Quantitites
  #---------------------------------------------------------------------------------------------------------
  XtXinv <- try(solve(crossprod(X)),silent=TRUE)
  if(is(XtXinv,"try-error")) XtXinv <- ginv(crossprod(X))
  A_OLS  <- XtXinv%*%(t(X)%*%Y)
  E_OLS  <- Y - X%*%A_OLS
  S_OLS  <- crossprod(E_OLS)/(bigT-k)
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  A_draw  <- array(A_OLS, c(bigT+1,k,M))
  S_draw  <- array(S_OLS, c(M,M,bigT))
  Em_draw <- Em_str <- E_OLS
  L_draw  <- diag(M)

  # time-varying stuff
  Am_draw    <- A_OLS
  At_draw    <- array(0, c(bigT+1, k, M))
  theta_draw <- matrix(1, k, M)
  theta_sqrt <- sqrt(theta_draw)
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  # Priors on VAR coefs
  #-----------------------------
  # prior mean
  A_prior <- matrix(0,2*k,M)
  diag(A_prior) <- prmean
  # prior variance
  tau2.draw <- matrix(10,k,M)
  xi2.draw  <- matrix(10,k,M)

  # NG stuff
  lambda2      <- matrix(10,p,1)
  a_tau        <- matrix(a_start,p,1)
  scale_tau    <- rep(.43,p)
  acc_tau      <- rep(0,p)

  kappa2       <- 10
  a_xi         <- a_start
  scale_xi     <- .43
  acc_xi       <- 0
  #------------------------------------
  # Priors on coefs in H matrix of VCV
  #------------------------------------
  # prior mean
  l_prior <- matrix(0,M,M)
  # prior variance
  L_prior <- matrix(10,M,M)
  L_prior[upper.tri(L_prior)] <- 0; diag(L_prior) <- 0

  # NG
  lambda2_L   <- 10
  a_L_tau     <- a_start
  scale_L_tau <- .43
  acc_L_tau   <- 0
  #------------------------------------
  # SV quantities
  #------------------------------------
  Sv_draw <- matrix(-3,bigT,M)
  svdraw <- list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,bigT))
  svl <- list()
  for (jj in 1:M) svl[[jj]] <- svdraw
  pars_var <- matrix(c(-3,.9,.2,-3),4,M,dimnames=list(c("mu","phi","sigma","latent0"),NULL))

  hv <- svdraw$latent
  para <- list(mu=-3,phi=.9,sigma=.2)
  Sv_priors <- specify_priors(mu=sv_normal(mean=bmu, sd=Bmu), phi=sv_beta(a0,b0), sigma2=sv_gamma(shape=0.5,rate=1/(2*Bsigma)))
  eta <- list()
  #---------------------------------------------------------------------------------------------------------
  # SAMPLER MISCELLANEOUS
  #---------------------------------------------------------------------------------------------------------
  nsave <- draws_in
  nburn <- burnin_in
  ntot  <- nsave+nburn

  # thinning
  thin         <- thin_in
  count <- 0
  thindraws    <- nsave/thin
  thin.draws   <- seq(nburn+1,ntot,by=thin)
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store      <- array(NA,c(thindraws,bigT+1,k,M))
  Am_store     <- array(NA,c(thindraws,k,M))
  At_store     <- array(NA,c(thindraws,bigT+1,k,M))
  L_store      <- array(NA,c(thindraws,M,M))
  res_store    <- array(NA,c(thindraws,bigT,M))
  Sv_store     <- array(NA,c(thindraws,bigT,M))
  pars_store   <- array(NA,c(thindraws,4,M))
  # # NG
  tau2_store   <- array(NA,c(thindraws,k,M))
  xi2_store    <- array(NA,c(thindraws,k,M))
  lambda2_store<- array(NA,c(thindraws,p,1))
  kappa2_store <- array(NA,c(thindraws,1,1))
  a_xi_store   <- array(NA,c(thindraws,1,1))
  a_tau_store  <- array(NA,c(thindraws,p,1))
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Sample coefficients
    for (mm in 1:M){
      if(mm==1){
        Ystar <- (Y[,mm]-X%*%Am_draw)*exp(-0.5*Sv_draw[,mm])
        Fstar <- (X%*%diag(theta_sqrt[,mm]))*exp(-0.5*Sv_draw[,mm])

        At_draw[,,mm] <- sample_McCausland(Ystar, Fstar)

        Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
        Z.i <- cbind(X,hadamard.prod(X,At_draw[2:(bigT+1),,mm]))*exp(-0.5*Sv_draw[,mm])

        Vpriorinv <- diag(1/c(tau2.draw[,mm],xi2.draw[,mm]))

        V_post <- try(chol2inv(chol(crossprod(Z.i)+Vpriorinv)),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv(crossprod(Z.i)+Vpriorinv)
        A_post <- V_post%*%(crossprod(Z.i,Y.i)+Vpriorinv%*%A_prior[,mm])

        A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(Z.i)),silent=TRUE)
        if (is(A.draw.i,"try-error")) A.draw.i <- matrix(mvrnorm(1,A_post,V_post),ncol(Z.i),1)

        Am_draw[,mm]    <- A.draw.i[1:k,,drop=FALSE]
        theta_sqrt[,mm] <- A.draw.i[(k+1):(2*k),,drop=FALSE]
        # compute errors
        Em_draw[,mm] <- Em_str[,mm] <- Y[,mm] - X%*%Am_draw[,mm] - apply(hadamard.prod(X,At_draw[2:(bigT+1),,mm]%*%diag(theta_sqrt[,mm])),1,sum)
      }else{
        Ystar <- (Y[,mm]-X%*%Am_draw)*exp(-0.5*Sv_draw[,mm])
        Fstar <- (X%*%diag(theta_sqrt[,mm]))*exp(-0.5*Sv_draw[,mm])

        At_draw[,,mm] <- sample_McCausland(Ystar, Fstar)

        Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
        Z.i <- cbind(X,hadamard.prod(X,At_draw[2:(bigT+1),,mm]),Em_draw[,1:(mm-1)])*exp(-0.5*Sv_draw[,mm])

        Vpriorinv <- diag(1/c(tau2.draw[,mm],xi2.draw[,mm],L_prior[mm,1:(mm-1)]))

        V_post <- try(chol2inv(chol((crossprod(Z.i)+Vpriorinv))),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv((crossprod(Z.i)+Vpriorinv))
        A_post <- V_post%*%(crossprod(Z.i,Y.i)+Vpriorinv%*%c(A_prior[,mm],l_prior[mm,1:(mm-1)]))

        A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(Z.i)),silent=TRUE)
        if (is(A.draw.i,"try-error")) A.draw.i <- matrix(mvrnorm(1,A_post,V_post),ncol(Z.i),1)

        Am_draw[,mm]        <- A.draw.i[1:k,,drop=FALSE]
        theta_sqrt[,mm]     <- A.draw.i[(k+1):(2*k),,drop=FALSE]
        L_draw[mm,1:(mm-1)] <- A.draw.i[(2*k+1):ncol(Z.i),,drop=FALSE]
        # compute errors
        Em_draw[,mm] <- Y[,mm]-X%*%Am_draw[,mm]-apply(hadamard.prod(X,At_draw[2:(bigT+1),,mm]%*%diag(theta_sqrt[,mm])),1,sum)
        Em_str[,mm]  <- Y[,mm]-X%*%Am_draw[,mm]-apply(hadamard.prod(X,At_draw[2:(bigT+1),,mm]%*%diag(theta_sqrt[,mm])),1,sum)-Em_draw[,1:(mm-1),drop=FALSE]%*%t(L_draw[mm,1:(mm-1),drop=FALSE])
      }
    }
    rownames(Am_draw)    <- colnames(X)
    theta_draw <- theta_sqrt^2
    #----------------------------------------------------------------------------
    # Step 3: Interweaving
    theta_sign <- sign(theta_sqrt)
    for(mm in 1:M){
      A_draw[,,mm] <- matrix(Am_draw[,mm],bigT+1,k,byrow=TRUE) + At_draw[,,mm]%*%diag(theta_sqrt[,mm])
      A_diff <- diff(At_draw[,,mm]%*%diag(theta_sqrt[,mm]))
      for(kk in 1:k){
        #theta.new
        res <- do_rgig1(lambda=-bigT/2,
                        chi=sum(A_diff[,kk]^2)+(A_draw[1,kk,mm]-Am_draw[kk,mm])^2,
                        psi=1/xi2.draw[kk,mm])
        theta_draw[kk,mm] <- res
        theta_sqrt[kk,mm] <- sqrt(res)*theta_sign[kk,mm]
        # Am_new
        sigma2_A_mean  <- 1/((1/tau2.draw[kk,mm]) + (1/theta_draw[kk,mm]))
        mu_A_mean      <- A_draw[1,kk,mm]*tau2.draw[kk,mm]/(tau2.draw[kk,mm] + theta_draw[kk,mm])
        Am_draw[kk,mm] <- rnorm(1, mu_A_mean, sqrt(sigma2_A_mean))
      }
      At_draw[,,mm] <- sapply(1:k,function(kk)A_draw[,kk,mm]-Am_draw[kk,mm])%*%diag(1/theta_sqrt[,mm])
    }
    #----------------------------------------------------------------------------
    # Step 4a: Shrinkage priors on state variances
    kappa2       <- rgamma(1, d1+a_xi*k,  d2+0.5*k*a_xi*mean(xi2.draw))
    for(ii in 1:k){
      for(jj in 1:M){
        xi2.draw[ii,jj]  <- do_rgig1(lambda=a_xi-0.5, chi=theta_draw[ii,jj], psi=a_xi*kappa2)
      }
    }
    xi2.draw[xi2.draw<1e-7] <- 1e-7
    if(sample_A){
      before <- a_xi
      a_xi   <- MH_step(a_xi, scale_xi, k, kappa2, as.vector(theta_sqrt), b_xi, nu_xi, d1, d2)
      if(before!=a_xi){
        acc_xi <- acc_xi + 1
      }
      # scale MH proposal during the first 50% of the burn-in stage
      if(irep<(0.5*burnin)){
        if((acc_xi/irep)>0.30){scale_xi <- 1.01*scale_xi}
        if((acc_xi/irep)<0.15){scale_xi <- 0.99*scale_xi}
      }
    }
    # Step 4b: Shrinkage prior on mean (multiplicative Gamma prior)
    for(pp in 1:p){
      slct.i       <- which(rownames(Am_draw)==paste("Ylag",pp,sep=""))
      if(pp==1 & cons) slct.i <- c(slct.i,which(rownames(Am_draw)=="cons"))
      if(pp==1 & trend) slct.i <- c(slct.i,which(rownames(Am_draw)=="trend"))
      Am_lag.i     <- Am_draw[slct.i,,drop=FALSE]
      A_prior.i    <- A_prior[slct.i,,drop=FALSE]
      tau2.i       <- tau2.draw[slct.i,,drop=FALSE]

      if(pp==1){
        lambda2[pp,1] <- rgamma(1, e1+a_tau[pp,1]*M^2, e2+0.5*a_tau[pp,1]*mean(tau2.i))
      }else{
        lambda2[pp,1] <- rgamma(1, e1+a_tau[pp,1]*M^2, e2+0.5*a_tau[pp,1]*prod(lambda2[1:(pp-1)])*mean(tau2.i))
      }
      Mend <- M + ifelse(pp==1&cons,1,0) + ifelse(pp==1&trend,1,0)
      for(ii in 1:Mend){
        for(jj in 1:M){
          tau2.i[ii,jj] <- do_rgig1(lambda=a_tau[pp,1]-0.5, chi=(Am_lag.i[ii,jj]-A_prior.i[ii,jj])^2, psi=a_tau[pp,1]*prod(lambda2[1:pp,1]))
        }
      }
      tau2.i[tau2.i<1e-7] <- 1e-7
      if(sample_A){
        before <- a_tau[pp,1]
        a_tau[pp,1]   <- MH_step(a_tau[pp,1], scale_xi[pp], M^2, lambda2[pp,1], as.vector(Am_lag.i), b_tau, nu_tau, e1, e2)
        if(before!=a_tau[pp,1]){
          acc_tau[pp] <- acc_tau[pp] + 1
        }
        # scale MH proposal during the first 50% of the burn-in stage
        if(irep<(0.5*burnin)){
          if((acc_tau[pp]/irep)>0.30){scale_xi[pp] <- 1.01*scale_xi[pp]}
          if((acc_tau[pp]/irep)<0.15){scale_xi[pp] <- 0.99*scale_xi[pp]}
        }
      }
      tau2.draw[slct.i,] <- tau2.i
    }
    # Step 4c: Shrinkage prior on covariances
    lambda2_L <- rgamma(1, e1+a_L_tau*v,  e2+0.5*v*a_L_tau*mean(L_prior[lower.tri(L_prior)]))
    for(ii in 2:M){
      for(jj in 1:(ii-1)){
        res  <- do_rgig1(lambda=a_L_tau-0.5, chi=(L_draw[mm,ii]-l_prior[mm,ii])^2, psi=a_L_tau*lambda2_L)
        L_prior[ii,jj] <- ifelse(res<1e-7,1e-7,res)
      }
    }
    if(sample_A){
      before  <- a_L_tau
      a_L_tau <- MH_step(a_L_tau, scale_L_tau, v, lambda2_L, L_draw[lower.tri(L_draw)], b_tau, nu_tau, e1, e2)
      if(before!=a_L_tau){
        acc_L_tau <- acc_L_tau + 1
      }
      # scale MH proposal during the first 50% of the burn-in stage
      if(irep<(0.5*burnin)){
        if((acc_L_tau/irep)>0.30){scale_L_tau <- 1.01*scale_L_tau}
        if((acc_L_tau/irep)<0.15){scale_L_tau <- 0.99*scale_L_tau}
      }
    }
    #----------------------------------------------------------------------------
    # Step 5: Sample variances
    if (sv){
      for (jj in 1:M){
        para   <- as.list(pars_var[,jj])
        para$nu = Inf; para$rho=0; para$beta<-0
        svdraw <- svsample_fast_cpp(y=Em_str[,jj], draws=1, burnin=0, designmatrix=matrix(NA_real_),
                                    priorspec=Sv_priors, thinpara=1, thinlatent=1, keeptime="all",
                                    startpara=para, startlatent=Sv_draw[,jj],
                                    keeptau=FALSE, print_settings=list(quiet=TRUE, n_chains=1, chain=1),
                                    correct_model_misspecification=FALSE, interweave=TRUE, myoffset=0,
                                    fast_sv=default_fast_sv)
        svl[[jj]]     <- svdraw
        h_            <- exp(svdraw$latent[1,])
        para$mu       <- svdraw$para[1,"mu"]
        para$phi      <- svdraw$para[1,"phi"]
        para$sigma    <- svdraw$para[1,"sigma"]
        para$latent0  <- svdraw$latent0[1,"h_0"]
        pars_var[,jj] <- unlist(para[c("mu","phi","sigma","latent0")])
        Sv_draw[,jj]  <- log(h_)
      }
    }else{
      for (jj in 1:M){
        S_1 <- a_1+bigT/2
        S_2 <- b_1+crossprod(Em_str[,jj])/2

        sig_eta <- 1/rgamma(1,S_1,S_2)
        Sv_draw[,jj] <- log(sig_eta)
      }
    }
    #-------------------------------------------------------------------------#
    # STEP 6: RANDOM SIGN SWITCH
    for(mm in 1:M){
      for(kk in 1:k){
        if(runif(1,0,1)>0.5){
          theta_sqrt[kk,mm] <- -theta_sqrt[kk,mm]
        }
      }
    }
    #----------------------------------------------------------------------------
    # Step 7: store draws
    if(irep %in% thin.draws){
      count <- count+1
      A_store[count,,,]      <- A_draw
      L_store[count,,]       <- L_draw
      res_store[count,,]     <- Em_draw
      # SV
      Sv_store[count,,]      <- Sv_draw
      pars_store[count,,]    <- pars_var
      # NG
      tau2_store[count,,]    <- tau2.draw
      xi2_store[count,,]     <- xi2.draw
      lambda2_store[count,,] <- lambda2
      kappa2_store[count,,]  <- kappa2
      a_xi_store[count,,]    <- a_xi
      a_tau_store[count,,]   <- a_tau
    }
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,paste("t",seq(0,bigT),sep="."),colnames(X),colnames(A_OLS))
  ret <- list(Y=Y,X=X,A_store=A_store,L_store=L_store,Sv_store=Sv_store,pars_store=pars_store,res_store=res_store,
              tau2_store=tau2_store,xi2_store=xi2_store,lambda2_store=lambda2_store,kappa2_store=kappa2_store,a_xi_store=a_xi_store,a_tau_store=a_tau_store)
  return(ret)
}


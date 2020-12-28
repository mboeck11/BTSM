#' @name .TVPBVAR_linear_wrapper
#' @noRd
#' @importFrom abind adrop
#' @importFrom utils capture.output
.TVPBVAR_linear_wrapper <- function(Yraw, prior, plag, draws, burnin, cons, trend, SV, thin, default_hyperpara, Ex, applyfun, cores){
  class(Yraw) <- "numeric"
  prior_in <- prior
  if(!is.null(default_hyperpara[["a_log"]])){
    default_hyperpara["a_start"] <- 1/log(ncol(Yraw))
  }
  nr <- 1
  Y_in=Yraw
  p_in=plag
  draws_in=draws
  burnin_in=burnin
  cons_in=cons
  trend_in=trend
  sv_in=SV
  thin_in=thin
  prior_in=prior_in
  hyperparam_in=default_hyperpara
  Ex_in=Ex
  if(prior=="TVP" || prior=="TVP-NG"){
    prior_in <- ifelse(prior=="TVP",1,2)
    post_draws <- applyfun(1:ncol(Yraw), function(nr){
      .TVPBVAR_noncentered_R(nr=nr,Y_in=Yraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)
    })
  }else if(prior=="TTVP"){
    prior_in <- 3
    post_draws <- applyfun(1:ncol(Yraw), function(nr){
      .TVPBVAR_centered_R(nr=nr,Y_in=Yraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)
    })
  }

  #tvpbvar<-.TVPBVAR_linear_R(Y_in=Yraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)
  #------------------------------------------------ get data ----------------------------------------#
  Y <- tvpbvar$Y; colnames(Y) <- colnames(Yraw); X <- tvpbvar$X
  M <- ncol(Y); bigT <- nrow(Y); K <- ncol(X)
  if(!is.null(Ex)) Mex <- ncol(Ex)
  xnames <- paste(rep("Ylag",M),rep(seq(1,plag),each=M),sep="")
  if(!is.null(Ex)) xnames <- c(xnames,paste(rep("Tex",Mex)))
  if(cons)  xnames <- c(xnames,"cons")
  if(trend) xnames <- c(xnames,"trend")
  colnames(X) <- xnames
  #-----------------------------------------get containers ------------------------------------------#
  A_store <- tvpbvar$A_store; dimnames(A_store)[[3]] <- colnames(X); dimnames(A_store)[[4]] <- colnames(Y)
  # splitting up stores
  dims          <- dimnames(A_store)[[3]]
  a0store <- a1store <- Exstore <- NULL
  if(cons) {
    a0store       <- adrop(A_store[,,which(dims=="cons"),,drop=FALSE],drop=3)
  }
  if(trend){
    a1store     <- adrop(A_store[,,which(dims=="trend"),,drop=FALSE],drop=3)
  }
  if(!is.null(Ex)){
    Exstore     <- A_store[,,which(dims=="Tex"),,drop=FALSE]
  }
  Phistore    <- NULL
  for(jj in 1:plag){
    Phistore[[jj]]  <- A_store[,,which(dims==paste("Ylag",jj,sep="")),,drop=FALSE]
  }
  S_store <- array(NA, c(draws/thin,bigT,M,M)); dimnames(S_store) <- list(NULL,NULL,colnames(Y),colnames(Y))
  if(prior%in%c("TVP-NG")){
    L_store <- tvpbvar$L_store
    for(irep in 1:(draws/thin)){
      for(tt in 1:bigT){
        if(M>1){
          S_store[irep,tt,,] <- L_store[irep,,]%*%diag(exp(tvpbvar$Sv_store[irep,tt,]))%*%t(L_store[irep,,])
        }else{
          S_store[irep,tt,,] <- L_store[irep,,]%*%exp(tvpbvar$Sv_store[irep,tt,])%*%t(L_store[irep,,])
        }
      }
    }
    Smed_store <- apply(S_store,c(1,3,4),median)
    if(SV){
      vola_store  <- tvpbvar$Sv_store; dimnames(vola_store) <- list(NULL,NULL,colnames(Y))
      pars_store  <- tvpbvar$pars_store
      vola_post   <- apply(vola_store,c(2,3),median)
      pars_post   <- apply(pars_store,c(2,3),median)
    }else{
      vola_store  <- tvpbvar$Sv_store; pars_store <- NULL;
      vola_post   <- apply(vola_store,c(2,3),median); pars_post <- NULL
    }
  }
  res_store       <- tvpbvar$res_store; dimnames(res_store) <- list(NULL,NULL,colnames(Y))
  # NG
  if(prior=="TVP-NG"){
    tau2_store    <- tvpbvar$tau2_store
    xi2_store     <- tvpbvar$xi2_store
    lambda2_store <- tvpbvar$lambda2_store
    kappa2_store  <- tvpbvar$kappa2_store
    a_tau_store   <- tvpbvar$a_tau_store
    a_xi_store    <- tvpbvar$a_xi_store
    tau2_post     <- apply(tau2_store,c(2,3),median)
    xi2_post      <- apply(xi2_store,c(2,3),median)
    lambda2_post  <- apply(lambda2_store,c(2,3),median)
    kappa2_post   <- apply(kappa2_store,c(2,3),median)
    a_tau_post    <- apply(a_tau_store,c(2,3),median)
    a_xi_post     <- apply(a_xi_store,c(2,3),median)
  }
  store <- list(A_store=A_store,a0store=a0store,a1store=a1store,Phistore=Phistore,Exstore=Exstore,S_store=S_store,Smed_store=Smed_store,
                L_store=L_store,vola_store=vola_store,pars_store=pars_store,res_store=res_store,
                tau2_store=tau2_store,xi2_store=xi2_store,lambda2_store=lambda2_store,kappa2_store=kappa2_store,a_tau_store=a_tau_store,a_xi_store=a_xi_store)
  #------------------------------------ compute posteriors -------------------------------------------#
  A_post      <- apply(A_store,c(2,3,4),median)
  S_post      <- apply(S_store,c(2,3,4),median)
  Sig         <- apply(S_post,c(2,3),mean)/(bigT-K)
  res_post    <- apply(res_store,c(2,3),median)
  # splitting up posteriors
  a0post <- a1post <- Expost <- NULL
  if(cons)  a0post <- A_post[,which(dims=="cons"),,drop=FALSE]
  if(trend) a1post <- A_post[,which(dims=="trend"),,drop=FALSE]
  if(!is.null(Ex)) Expost <- A_post[,which(dims=="Tex"),,drop=FALSE]
  Phipost     <- NULL
  for(jj in 1:plag){
    Phipost[[jj]]    <- A_post[,which(dims==paste("Ylag",jj,sep="")),,drop=FALSE]
  }
  post <- list(A_post=A_post,a0post=a0post,a1post=a1post,Phipost=Phipost,Expost=Expost,S_post=S_post,Sig=Sig,vola_post=vola_post,pars_post=pars_post,res_post=res_post,
               tau2_post=tau2_post,xi2_post=xi2_post,lambda2_post=lambda2_post,kappa2_post=kappa2_post,a_tau_post=a_tau_post,a_xi_post=a_xi_post)
  return(list(Y=Y,X=X,store=store,post=post))
}

#' @name .TVPBVAR_noncentered_R.m
#' @importFrom stochvol svsample_fast_cpp specify_priors default_fast_sv
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

  texo <- FALSE; Mex <- 0; Exraw <- NULL
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
  Sv_store     <- array(NA,c(thindraws,bigT,1))
  pars_store   <- array(NA,c(thindraws,4,1))
  tau2_store   <- array(NA,c(thindraws,d))
  xi2_store    <- array(NA,c(thindraws,d))
  # TVP-NG
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
      res <- GIGrvg::rgig(1,
                          lambda=-bigT/2,
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
      kappa2  <- rgamma(1, d1+a_xi*M^2,  d2+0.5*a_xi*mean(xi2_draw))
      lambda2 <- rgamma(1, e1+a_tau*M^2, e2+0.5*a_tau*mean(tau2_draw))
      for(dd in 1:d){
        xi2_draw[dd]  <- do_rgig1(lambda=a_xi-0.5,  chi=theta_draw[dd], psi=a_xi*kappa2)
        tau2_draw[dd] <- do_rgig1(lambda=a_tau-0.5, chi=(Am_draw[dd,1]-A_prior[dd,1])^2, psi=a_tau*lambda2)
      }
      xi2_draw[xi2_draw<1e-7] <- 1e-7
      tau2_draw[tau2_draw<1e-7] <- 1e-7
      if(sample_A){
        before <- a_xi
        a_xi   <- MH_step(a_xi, scale_xi, d*M, kappa2, theta_sqrt, b_xi, nu_xi, d1, d2)
        if(before!=a_xi){
          acc_xi <- acc_xi + 1
        }
        before <- a_tau
        a_tau  <- MH_step(a_tau, scale_xi, M^2, lambda2, Am_draw, b_tau, nu_tau, e1, e2)
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
      tau2_store[count,]<- tau2_draw
      xi2_store[count,] <- xi2_draw
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
  ret <- list(Y=Y,X=X,A_store=A_store,Sv_store=Sv_store,pars_store=pars_store,res_store=res_store,
              tau2_store=tau2_store,xi2_store=xi2_store,lambda2_store=lambda2_store,kappa2_store=kappa2_store,a_xi_store=a_xi_store,a_tau_store=a_tau_store)
  return(ret)
}

#' @name .TVPBVAR_centered_R.m
#' @importFrom stochvol svsample_fast_cpp specify_priors default_fast_sv
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

  texo <- FALSE; Mex <- 0; Exraw <- NULL
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
  cons.mod       <- hyperpara$cons.mod
  robust         <- hyperpara$robust
  a.approx       <- hyperpara$a.approx
  sim.kappa      <- hyperpara$sim.kappa
  kappa.grid     <- hyperpara$kappa.grid
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
  Sv_store     <- array(NA,c(thindraws,bigT,1))
  pars_store   <- array(NA,c(thindraws,4,1))
  tau2_store   <- array(NA,c(thindraws,d))
  xi2_store    <- array(NA,c(thindraws,d))
  # TVP-NG
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
      res <- GIGrvg::rgig(1,
                          lambda=-bigT/2,
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
      kappa2  <- rgamma(1, d1+a_xi*M^2,  d2+0.5*a_xi*mean(xi2_draw))
      lambda2 <- rgamma(1, e1+a_tau*M^2, e2+0.5*a_tau*mean(tau2_draw))
      for(dd in 1:d){
        xi2_draw[dd]  <- do_rgig1(lambda=a_xi-0.5,  chi=theta_draw[dd], psi=a_xi*kappa2)
        tau2_draw[dd] <- do_rgig1(lambda=a_tau-0.5, chi=(Am_draw[dd,1]-A_prior[dd,1])^2, psi=a_tau*lambda2)
      }
      xi2_draw[xi2_draw<1e-7] <- 1e-7
      tau2_draw[tau2_draw<1e-7] <- 1e-7
      if(sample_A){
        before <- a_xi
        a_xi   <- MH_step(a_xi, scale_xi, d*M, kappa2, theta_sqrt, b_xi, nu_xi, d1, d2)
        if(before!=a_xi){
          acc_xi <- acc_xi + 1
        }
        before <- a_tau
        a_tau  <- MH_step(a_tau, scale_xi, M^2, lambda2, Am_draw, b_tau, nu_tau, e1, e2)
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
      tau2_store[count,]<- tau2_draw
      xi2_store[count,] <- xi2_draw
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
  ret <- list(Y=Y,X=X,A_store=A_store,Sv_store=Sv_store,pars_store=pars_store,res_store=res_store,
              tau2_store=tau2_store,xi2_store=xi2_store,lambda2_store=lambda2_store,kappa2_store=kappa2_store,a_xi_store=a_xi_store,a_tau_store=a_tau_store)
  return(ret)
}

#' @name .TVPBVAR_linear_R
#' @importFrom stochvol svsample_fast_cpp specify_priors default_fast_sv
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
    kappa2       <- rgamma(1, d1+a_xi*M^2,  d2+0.5*a_xi*mean(xi2.draw))
    for(ii in 1:k){
      for(jj in 1:M){
        xi2.draw[ii,jj]  <- do_rgig1(lambda=a_xi-0.5, chi=theta_draw[ii,jj], psi=a_xi*kappa2)
      }
    }
    xi2.draw[xi2.draw<1e-7] <- 1e-7
    if(sample_A){
      before <- a_xi
      a_xi   <- MH_step(a_xi, scale_xi, k*M, kappa2, as.vector(theta_sqrt), b_xi, nu_xi, d1, d2)
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
    lambda2_L <- rgamma(1, e1+a_L_tau*v,  e2+0.5*a_L_tau*mean(L_prior[lower.tri(L_prior)]))
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

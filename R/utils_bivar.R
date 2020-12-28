#' @name .BVAR_linear_wrapper
#' @noRd
#' @importFrom abind adrop
#' @importFrom utils capture.output
.BIVAR_linear_wrapper <- function(Yraw, Draw, prior, plag, draws, burnin, cons, trend, SV, thin, default_hyperpara, Ex, applyfun, cores, eigen, trim){
  class(Yraw) <- class(Draw) <- "numeric"
  prior_in <- ifelse(prior=="MN",1,ifelse(prior=="SSVS",2,3))
  if(default_hyperpara[["a_log"]]){
    default_hyperpara["a_start"] <- 1/log(ncol(Yraw))
  }
  bivar<-.BIVAR_linear_R(Y_in=Yraw,D_in=Draw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)
  #------------------------------------------------ get data ----------------------------------------#
  Y <- bivar$Y; colnames(Y) <- colnames(Yraw); X <- bivar$X; D <- bivar$D; colnames(D) <- colnames(Draw)
  M <- ncol(Y); bigT <- nrow(Y); K <- ncol(X); Ki <- ncol(D); Ki1 <- Ki+1
  if(!is.null(Ex)) Mex <- ncol(Ex)
  names <- colnames(Yraw)
  if(is.null(names)) names <- rep("Y",M)
  xnames <- dnames <- NULL
  for(ii in 1:plag) xnames <- c(xnames,paste0(names,".lag",ii))
  for(ii in 1:plag) dnames <- c(dnames,rep(paste0("D",names,".lag",ii),Ki))
  if(!is.null(Ex)) enames <- paste(rep("Tex",Mex)) else enames <- NULL
  if(cons)  cnames <- c("cons",paste0("Dcons",seq(1,Ki))) else cnames <- NULL
  if(trend) tnames <- c("trend",paste0("Dtrend",seq(1,Ki))) else tnames <- NULL
  colnames(X) <- c(xnames,dnames,enames,cnames,tnames)
  #-----------------------------------------get containers ------------------------------------------#
  Atilde_store <- bivar$A_store; dimnames(Atilde_store)[[2]] <- colnames(X); dimnames(Atilde_store)[[3]] <- colnames(Y)
  # splitting up stores
  dims  <- dimnames(Atilde_store)[[2]]
  a0tilde_store <- a1tilde_store <- Extilde_store <- Phitilde_store <- NULL
  if(cons) a0tilde_store <- Atilde_store[,which(dims%in%cnames),]
  if(trend) a1tilde_store <- Atilde_store[,which(dims%in%tnames),]
  if(!is.null(Ex)) Extilde_store <- Atilde_store[,which(dims%in%enames),,drop=FALSE]
  for(jj in 1:plag) {
    xnames.jj <- xnames[grepl(paste0("lag",jj),xnames)]
    dnames.jj <- dnames[grepl(paste0("lag",jj),dnames)]
    Phitilde_store[[jj]] <- abind(Atilde_store[,which(dims%in%xnames.jj),],Atilde_store[,which(dims%in%dnames.jj),],along=4)
    dimnames(Phitilde_store[[jj]]) <- list(NULL,xnames.jj,names,c("",rep("D",Ki)))
  }
  J_store <- bivar$J_store
  if(SV){
    volatilde_store <- bivar$Sv_store; dimnames(vola_store) <- list(NULL,NULL,colnames(Y))
    parstilde_store <- bivar$pars_store
  }else{
    volatilde_store <- bivar$Sv_store; parstilde_store <- NULL;
  }
  theta_store   <- bivar$theta_store; dimnames(theta_store)[[2]] <- colnames(X); dimnames(theta_store)[[3]] <- colnames(Y)
  zeta_store    <- bivar$zeta_store
  # MN
  if(prior=="MN"){
    shrink_store  <- bivar$shrink_store; dimnames(shrink_store) <- list(NULL,c("shrink1","shrink2","shrink4"))
    shrink_post   <- apply(shrink_store,2,median)
  }else{
    shrink_store  <- NULL
  }
  # SSVS
  if(prior=="SSVS"){
    gamma_store <- bivar$gamma_store; dimnames(gamma_store) <- list(NULL,colnames(X),colnames(Y))
    omega_store <- bivar$omega_store; dimnames(omega_store) <- list(NULL,colnames(Y),colnames(Y))
    PIP         <- apply(gamma_store,c(2,3),mean)
    PIP_omega   <- apply(omega_store,c(2,3,4),mean)
  }else{
    gamma_store <- omega_store <- NULL
  }
  # NG
  if(prior=="NG"){
    lambda2_store <- bivar$lambda2_store
    tau_store     <- bivar$tau_store
    dimnames(lambda2_store) <- list(NULL,paste("lag",1:plag,sep="_"),c("endogenous","covariance"))
    dimnames(lambda2_store) <- list(NULL,paste("lag",1:plag,sep="_"),c("endogenous","covariance"))
    lambda2_post  <- apply(lambda2_store,c(2,3),median)
    tau_post      <- apply(tau_store,c(2,3),median)
  }else{
    lambda2_store <- tau_store <- NULL
  }
  #-----------------------------------------check eigenvalues ---------------------------------------#
  if(eigen){
    D.quantile <- cbind(1,matrix(apply(D,2,quantile,seq(.1,.9,by=.1)),9,Ki))
    A.eigen <- NULL
    for(dd in 1:9){
      A.eigen.list <- applyfun(1:draws,function(irep){
        J <- diag(M)
        for(mm in 2:M){
          for(jj in 1:(mm-1)){
            for(kk in 1:Ki1) J[mm,jj] <- J[mm,jj]+J_store[irep,mm,jj,kk]*D.quantile[dd,kk]
          }
        }
        Jinv <- solve(J)
        Phi <- matrix(0,M*plag,M)
        for(pp in 1:plag){
          for(kk in 1:Ki1) Phi[((pp-1)*M+1):(pp*M),] <- Phi[((pp-1)*M+1):(pp*M),] + Phitilde_store[[pp]][irep,,,kk]*D.quantile[dd,kk]
        }
        Cm <- .gen_compMat(Phi,M,plag)$Cm
        return(max(abs(Re(eigen(Cm)$values))))
      })
      A.eigen <- cbind(A.eigen,unlist(A.eigen.list))
    }
    A.eigen <- apply(A.eigen,1,max)
    trim_eigen <- which(A.eigen<trim)
    Atilde_store<-Atilde_store[trim_eigen,,]
    if(cons) a0tilde_store<-a0tilde_store[trim_eigen,,]
    if(trend) a1tilde_store<-a1tilde_store[trim_eigen,,]
    Phitilde_store<-lapply(Phitilde_store,function(l)l[trim_eigen,,,])
    if(!is.null(Ex)) Extilde_store<-Extilde_store[trim_eigen,]
    J_store<-J_store[trim_eigen,,,]
    theta_store<-theta_store[trim_eigen,,]
    zeta_store<-zeta_store[trim_eigen,,,]
    volatilde_store<-volatilde_store[trim_eigen,,]
    parstilde_store<-parstilde_store[trim_eigen,,]
    if(prior=="MN"){
      shrink_store<-shrink_store[trim_eigen,,drop=FALSE]
    }
    if(prior=="SSVS"){
      gamma_store<-gamma_store[trim_eigen,,,drop=FALSE]
      omega_store<-omega_store[trim_eigen,,,,drop=FALSE]
    }
    if(prior=="NG"){
      lambda2_store<-lambda2_store[trim_eigen,,,drop=FALSE]
      tau_store<-tau_store[trim_eigen,,,drop=FALSE]
    }
  }else{A.eigen<-NULL}
  #-----------------------------------------get containers ------------------------------------------#
  store <- list(Atilde_store=Atilde_store,a0tilde_store=a0tilde_store,a1tilde_store=a1tilde_store,Phitilde_store=Phitilde_store,Extilde_store=Extilde_store,
                J_store=J_store,theta_store=theta_store,zeta_store=zeta_store,volatilde_store=volatilde_store,parstilde_store=parstilde_store,shrink_store=shrink_store,
                gamma_store=gamma_store,omega_store=omega_store,lambda2_store=lambda2_store,tau_store=tau_store,A.eigen=A.eigen)
  #------------------------------------ compute posteriors -------------------------------------------#
  Atilde_post <- apply(Atilde_store,c(2,3),median)
  J_post <- apply(J_store,c(2,3,4),median)
  # splitting up posteriors
  a0tilde_post <- a1tilde_post <- Extilde_post <- Phitilde_post <- NULL
  if(cons)  a0tilde_post <- Atilde_post[which(dims%in%cnames),,drop=FALSE]
  if(trend) a1tilde_post <- Atilde_post[which(dims%in%tnames),,drop=FALSE]
  if(!is.null(Ex)) Extilde_post <- Atilde_post[which(dims%in%enames),,drop=FALSE]
  for(jj in 1:plag){
    xnames.jj <- xnames[grepl(paste0("lag",jj),xnames)]
    dnames.jj <- dnames[grepl(paste0("lag",jj),dnames)]
    Phitilde_post[[jj]] <- abind(Atilde_post[which(dims%in%xnames.jj),],Atilde_post[which(dims%in%dnames.jj),],along=3)
    dimnames(Phitilde_post[[jj]]) <- list(xnames.jj,names,c("",rep("D",Ki)))
  }
  theta_post  <- apply(theta_store,c(2,3),median)
  zeta_post   <- apply(zeta_store,c(2,3,4),median)
  if(SV){
    volatilde_post  <- apply(volatilde_store,c(2,3),median)
    parstilde_post  <- apply(parstilde_store,c(2,3),median)
  }else{
    volatilde_post  <- apply(volatilde_store,c(2,3),median); parstilde_post <- NULL
  }
  # MN
  if(prior=="MN"){
    shrink_post   <- apply(shrink_store,2,median)
  }else{
    shrink_post <- NULL
  }
  # SSVS
  if(prior=="SSVS"){
    PIP         <- apply(gamma_store,c(2,3),mean)
    PIP_omega   <- apply(omega_store,c(2,3,4),mean)
  }else{
    PIP <- PIP_omega <- NULL
  }
  # NG
  if(prior=="NG"){
    lambda2_post  <- apply(lambda2_store,c(2,3),median)
    tau_post      <- apply(tau_store,c(2,3),median)
  }else{
    lambda2_post <- tau_post <- NULL
  }
  post <- list(Atilde_post=Atilde_post,a0tilde_post=a0tilde_post,a1tilde_post=a1tilde_post,Phitilde_post=Phitilde_post,J_post=J_post,Extilde_post=Extilde_post,theta_post=theta_post,
               volatilde_post=volatilde_post,parstilde_post=parstilde_post,shrink_post=shrink_post,PIP=PIP,PIP_omega=PIP_omega,zeta_post=zeta_post,
               lambda2_post=lambda2_post,tau_post=tau_post)
  return(list(Y=Y,X=X,D=D,store=store,post=post))
}

#' @name .BIVAR_linear_R
#' @importFrom stochvol svsample_fast_cpp specify_priors default_fast_sv
#' @importFrom MASS ginv mvrnorm
#' @importFrom methods is
#' @importFrom matrixcalc hadamard.prod
#' @importFrom stats rnorm rgamma runif dnorm
#' @noRd
.BIVAR_linear_R <- function(Y_in,D_in,p_in,draws_in,burnin_in,cons_in,trend_in,sv_in,thin_in,quiet_in,prior_in,hyperparam_in,Ex_in){
  #----------------------------------------INPUTS----------------------------------------------------#
  Yraw  <- Y_in
  Draw  <- D_in
  p     <- p_in
  Traw  <- nrow(Yraw)
  M     <- ncol(Yraw)
  K     <- M*p
  Ki    <- ncol(Draw)
  Ki1   <- Ki+1
  Ylag  <- .mlag(Yraw,p)
  names <- colnames(Yraw)
  if(is.null(names)) names <- rep("Y",M)
  nameslags <- NULL
  for(ii in 1:p) nameslags <- c(nameslags,paste0(names,".lag",ii))
  colnames(Ylag) <- nameslags

  DYlag <- hadamard.prod(matrix(Ylag,Traw,Ki*M*plag),matrix(Draw,Traw,Ki*M*plag))
  colnames(DYlag) <- paste0("D",nameslags)

  DYraw <- matrix(NA,Traw,Ki1*M)
  for(mm in 1:M){
    idx <- seq((mm-1)*Ki1+1,mm*Ki1)
    DYraw[,idx] <- cbind(Yraw[,mm],hadamard.prod(matrix(Yraw[,mm],Traw,Ki),matrix(Draw,Traw,1)))
  }
  colnames(DYraw) <- paste0(rep(c("","D"),M),rep(names,each=Ki1))
  DYmat <- matrix(seq(1,M*Ki1),Ki1,M)

  texo <- FALSE; Mex <- 0; Exraw <- NULL
  if(!is.null(Ex_in)){
    Exraw <- Ex_in; Mex <- ncol(Exraw)
    texo <- TRUE
    colnames(Exraw) <- rep("Tex",Mex)
  }

  X  <- cbind(Ylag,DYlag,Exraw)
  X  <- X[(p+1):nrow(X),,drop=FALSE]
  Y  <- Yraw[(p+1):Traw,,drop=FALSE]
  D  <- Draw[(p+1):Traw,,drop=FALSE]
  DY <- DYraw[(p+1):Traw,,drop=FALSE]
  bigT  <- nrow(X)

  cons  <- cons_in
  if(cons){
    X <- cbind(X,1,D)
    colnames(X)[(ncol(X)-Ki):ncol(X)] <- c("cons",paste0("Dcons",seq(1,Ki)))
  }
  trend <- trend_in
  if(trend){
    X <- cbind(X,seq(1,bigT),seq(1,bigT)*D)
    colnames(X)[(ncol(X)-Ki):ncol(X)] <- c("trend",paste0("Dtrend",seq(1,Ki)))
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
  # prior == 1: MN
  shrink1   <- hyperpara$shrink1
  shrink2   <- hyperpara$shrink2
  shrink3   <- hyperpara$shrink3
  shrink4   <- hyperpara$shrink4
  # prior == 2: SSVS
  tau00     <- hyperpara$tau0
  tau11     <- hyperpara$tau1
  p_i       <- hyperpara$p_i
  kappa00   <- hyperpara$kappa0
  kappa11   <- hyperpara$kappa1
  q_ij      <- hyperpara$q_ij
  # prior == 3: NG
  d_lambda  <- hyperpara$d_lambda
  e_lambda  <- hyperpara$e_lambda
  a_start   <- hyperpara$a_start
  sample_A  <- hyperpara$sample_A
  #---------------------------------------------------------------------------------------------------------
  # OLS Quantitites
  #---------------------------------------------------------------------------------------------------------
  XtXinv <- try(solve(crossprod(X)),silent=TRUE)
  if(is(XtXinv,"try-error")) XtXinv <- ginv(crossprod(X))
  A_OLS <- matrix(0, k, M, dimnames=list(colnames(X),names))
  J_OLS <- array(0, c(M,M,Ki1)); dimnames(J_OLS)[[3]] <- c("NoI",colnames(D)); for(kk in 1:Ki1) J_OLS[,,kk] <- diag(M)
  S_OLS <- matrix(0,M,M)
  E_OLS <- matrix(NA,bigT,M)
  V_OLS <- matrix(0, k, M)
  Z_OLS <- array(0, c(M,M,Ki1))
  X1 <- X
  for(mm in 1:M){
    if(mm>1) X1 <- cbind(DY[,c(DYmat[,1:(mm-1)])],X) else X1 <- X
    XtXinv1 <- solve(crossprod(X1))
    temp    <- XtXinv1%*%t(X1)%*%Y[,mm]
    A_OLS[,mm]  <- temp[((mm-1)*Ki1+1):nrow(temp),]
    if(mm>1) for(kk in 1:Ki1) J_OLS[mm,1:(mm-1),kk] <- temp[(mm-2)*Ki1+kk,]
    E_OLS[,mm]   <- Y[,mm] - X1%*%temp
    S_OLS[mm,mm] <- crossprod(E_OLS[,mm])/(bigT-ncol(X1))
    temp         <- diag(XtXinv1*S_OLS[mm,mm])
    V_OLS[,mm] <- temp[((mm-1)*Ki1+1):length(temp)]
    if(mm>1) for(kk in 1:Ki1) Z_OLS[mm,1:(mm-1),kk] <- temp[(mm-2)*Ki1+kk]
  }
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  A_draw <- A_OLS
  J_draw <- J_OLS
  SIGMA  <- array(S_OLS, c(M,M,bigT))
  Em_str <- E_OLS
  Em     <- E_OLS
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  # Priors on VAR coefs
  #-----------------------------
  # prior mean A_draw
  A_prior <- matrix(0,k,M)
  diag(A_prior) <- prmean
  a_prior  <-  as.vector(A_prior)
  # prior variance A_draw
  theta <- V_OLS

  # prior mean J_draw
  J_prior <- array(0,c(M,M,Ki1))
  # prior variance
  zeta <- Z_OLS

  # MN A_draw
  accept1 <- accept2 <- accept4 <- 0
  scale1  <- scale2  <- scale4 <- .43
  sigma_sq  <- matrix(0,M,1) #vector which stores the residual variance
  sigma_co <- matrix(1,M,1)
  for (i in 1:M){
    Ylag_i         <- .mlag(Yraw[,i],p)
    Ylag_i         <- Ylag_i[(p+1):nrow(Ylag_i),,drop=FALSE]
    Y_i            <- Yraw[(p+1):nrow(Yraw),i,drop=FALSE]
    Ylag_i         <- cbind(Ylag_i,seq(1,nrow(Y_i)))
    alpha_i        <- solve(crossprod(Ylag_i))%*%crossprod(Ylag_i,Y_i)
    sigma_sq[i,1]  <- (1/(nrow(Y_i)-p-1))*t(Y_i-Ylag_i%*%alpha_i)%*%(Y_i-Ylag_i%*%alpha_i)
  }
  for(mm in 2:M){
    Ylag_i         <- cbind(Yraw[,1:(mm-1)],Draw)
    Y_i            <- Yraw[,mm,drop=FALSE]
    alpha_i        <- solve(crossprod(Ylag_i))%*%crossprod(Ylag_i,Y_i)
    sigma_co[mm,1] <- (1/(nrow(Y_i)-(mm-1)-Ki))*crossprod(Y_i-Ylag_i%*%alpha_i)
  }
  if(prior==1){
    temp <- .get_V2(M=M,p=p,Ki1=Ki1,a_bar_1=shrink1,a_bar_2=shrink2,a_bar_3=shrink3,a_bar_4=shrink4,
                     sigma_sq=sigma_sq,sigma_co=sigma_co,cons=cons,trend=trend)
    theta <- temp$theta
    zeta <- temp$zeta
  }

  # SSVS stuff
  gamma  <-  matrix(1,k,M)
  tau0 <- matrix(NA, k, M); tau1 <- matrix(NA, k, M)
  for(mm in 1:M){
    for(kk in 1:k){
      tau0[kk,mm] <- tau00*sqrt(V_OLS[kk,mm])
      tau1[kk,mm] <- tau11*sqrt(V_OLS[kk,mm])
    }
  }
  omega <- array(1,c(M,M,Ki1))
  kappa0 <- array(NA, c(M,M,Ki1)); kappa1 <- array(NA, c(M, M,Ki1))
  for(kk in 1:Ki1){
    for(mm in 2:M){
      for(mmm in 1:(mm-1)){
        kappa0[mm,mmm,kk] <- kappa00*sqrt(Z_OLS[mm,mmm,kk])
        kappa1[mm,mmm,kk] <- kappa11*sqrt(Z_OLS[mm,mmm,kk])
      }
    }
  }

  # NG A_draw
  lambda2_A    <- matrix(0.01,p,1)
  A_tau        <- matrix(a_start,p,1)
  colnames(A_tau) <- colnames(lambda2_A) <- "endo"
  rownames(A_tau) <- rownames(lambda2_A) <- paste("lag.",seq(1,p),sep="")
  A_tuning     <- matrix(.43,p,1)
  A_accept     <- matrix(0,p,1)
  # NG J_draw
  lambda2_J <- 0.01
  J_tau     <- a_start
  J_accept  <- 0
  J_tuning  <- .43
  #------------------------------------
  # SV quantities
  #------------------------------------
  Sv_draw <- matrix(-3,bigT,M)
  svdraw <- list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,bigT))
  svl <- list()
  for (jj in 1:M) svl[[jj]] <- svdraw
  pars_var <- matrix(c(-3,.9,.2,-3),4,M,dimnames=list(c("mu","phi","sigma","latent0"),NULL))
  Sv_priors <- specify_priors(mu=sv_normal(mean=bmu, sd=Bmu), phi=sv_beta(a0,b0), sigma2=sv_gamma(shape=0.5,rate=1/(2*Bsigma)))
  hv <- svdraw$latent
  para <- list(mu=-3,phi=.9,sigma=.2)
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
  A_store      <- array(NA,c(thindraws,k,M))
  J_store      <- array(NA,c(thindraws,M,M,Ki1))
  # SV
  Sv_store     <- array(NA,c(thindraws,bigT,M))
  pars_store   <- array(NA,c(thindraws,4,M))
  # MN
  shrink_store <- array(NA,c(thindraws,3))
  # SSVS
  gamma_store  <- array(NA,c(thindraws,k,M))
  omega_store  <- array(NA,c(thindraws,M,M,Ki1))
  # NG
  theta_store  <- array(NA,c(thindraws,k,M))
  lambda2_store<- array(NA,c(thindraws,p,2))
  tau_store    <- array(NA,c(thindraws,p,2))
  zeta_store   <- array(NA,c(thindraws,M,M,Ki1))
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Sample coefficients
    for(mm in 1:M){
      Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
      if(mm>1) X1 <- cbind(DY[,c(DYmat[,1:(mm-1)])],X) else X1 <- X
      X.i <- X1*exp(-0.5*Sv_draw[,mm])

      Aprior    <- A_prior[,mm,drop=FALSE]
      Vprior    <- theta[,mm,drop=FALSE]
      if(mm>1) {
        for(kk in 1:Ki1) {
          Aprior <- rbind(matrix(J_prior[mm,1:(mm-1),kk],mm-1,1),Aprior)
          Vprior <- rbind(matrix(zeta[mm,1:(mm-1),kk],mm-1,1),Vprior)
        }
      }
      Vpriorinv <- diag(1/c(Vprior))

      V_post <- try(chol2inv(chol(crossprod(X.i)+Vpriorinv)),silent=TRUE)
      if (is(V_post,"try-error")) V_post <- ginv(crossprod(X.i)+Vpriorinv)
      A_post <- V_post%*%(crossprod(X.i,Y.i)+Vpriorinv%*%Aprior)

      A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
      if (is(A.draw.i,"try-error")) A.draw.i <- mvrnorm(1,A_post,V_post)
      A_draw[,mm] <- A.draw.i[((mm-1)*(1+Ki)+1):nrow(A.draw.i),]
      if(mm>1) for(kk in 1:Ki1) J_draw[mm,1:(mm-1),kk] <- -A.draw.i[(mm-2)*Ki1+kk,]
      Em_str[,mm] <- Y[,mm]-X1%*%A.draw.i
    }
    #----------------------------------------------------------------------------
    # Step 2: different shrinkage prior setups
    # MN
    if(prior==1){
      #Step for the first shrinkage parameter (own lags)
      shrink1.prop <- exp(rnorm(1,0,scale1))*shrink1
      if(shrink1.prop<1e-17) shrink1.prop <- 1e-17
      if(shrink1.prop>1e+17) shrink1.prop <- 1e+17
      theta1.prop   <- .get_V2(M=M,p=p,Ki1=Ki1,a_bar_1=shrink1,a_bar_2=shrink2,a_bar_3=shrink3,a_bar_4=shrink4,
                               sigma_sq=sigma_sq,sigma_co=sigma_co,cons=cons,trend=trend)$theta
      post1.prop<-sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta1.prop)),log=TRUE))+dgamma(shrink1.prop,0.01,0.01,log=TRUE)
      post1.prop<-post1.prop+log(shrink1.prop) # correction term
      post1 <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta)),log=TRUE))+dgamma(shrink1,0.01,0.01,log=TRUE)
      post1 <- post1+log(shrink1) # correction term
      if ((post1.prop-post1)>log(runif(1,0,1))){
        shrink1 <- shrink1.prop
        theta   <- theta1.prop
        accept1 <- accept1+1
      }

      #Step for the second shrinkage parameter (cross equation)
      shrink2.prop <- exp(rnorm(1,0,scale2))*shrink2
      if(shrink2.prop<1e-17) shrink2.prop <- 1e-17
      if(shrink2.prop>1e+17) shrink2.prop <- 1e+17
      theta2.prop   <- .get_V2(M=M,p=p,Ki1=Ki1,a_bar_1=shrink1,a_bar_2=shrink2,a_bar_3=shrink3,a_bar_4=shrink4,
                               sigma_sq=sigma_sq,sigma_co=sigma_co,cons=cons,trend=trend)$theta
      post2.prop <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta2.prop)),log=TRUE))+dgamma(shrink2.prop,0.01,0.01,log=TRUE)
      post2.prop <- post2.prop + log(shrink2.prop) # correction term
      post2 <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta)),log=TRUE))+dgamma(shrink2,0.01,0.01,log=TRUE)
      post2 <- post2 + log(shrink2) # correction term
      if ((post2.prop-post2)>log(runif(1,0,1))){
        shrink2 <- shrink2.prop
        theta   <- theta2.prop
        accept2 <- accept2+1
      }

      # #Step for the fourth shrinkage parameter (covariances)
      # shrink4.prop <- exp(rnorm(1,0,scale4))*shrink4
      # if(shrink4.prop<1e-17) shrink4.prop <- 1e-17
      # if(shrink4.prop>1e+17) shrink4.prop <- 1e+17
      # zeta4.prop <- .get_V2(M=M,p=p,Ki1=Ki1,a_bar_1=shrink1,a_bar_2=shrink2,a_bar_3=shrink3,a_bar_4=shrink4,
      #                       sigma_sq=sigma_sq,sigma_co=sigma_co,cons=cons,trend=trend)$zeta
      # post4.prop <- post4 <- 0
      # for(kk in 1:Ki1){
      #   J_draw.kk  <- J_draw[,,kk]
      #   J_prior.kk <- J_prior[,,kk]
      #   zeta.kk    <- zeta[,,kk]
      #   zeta4.kk   <- zeta4.prop[,,kk]
      #   post4.prop <- post4.prop + sum(dnorm(as.vector(J_draw.kk[lower.tri(J_draw.kk)]),as.vector(J_prior.kk[lower.tri(J_prior.kk)]),sqrt(as.vector(zeta4.kk[lower.tri(zeta4.kk)])),log=TRUE)) +
      #     dgamma(shrink4.prop,0.01,0.01,log=TRUE)
      #   post4      <- post4 + sum(dnorm(as.vector(J_draw.kk[lower.tri(J_draw.kk)]),as.vector(J_prior.kk[lower.tri(J_prior.kk)]),sqrt(as.vector(zeta.kk[lower.tri(zeta.kk)])),log=TRUE)) +
      #     dgamma(shrink4.prop,0.01,0.01,log=TRUE)
      # }
      # # correction term
      # post4.prop <- post4.prop + log(shrink4.prop)
      # post4 <- post4 + log(shrink4)
      # if ((post4.prop-post4)>log(runif(1,0,1))){
      #   shrink4 <- shrink4.prop
      #   zeta    <- zeta4.prop
      #   accept4 <- accept4+1
      # }

      if (irep<(0.5*nburn)){
        if ((accept1/irep)<0.15) scale1 <- 0.99*scale1
        if ((accept1/irep)>0.3)  scale1 <- 1.01*scale1
        if ((accept2/irep)<0.15) scale2 <- 0.99*scale2
        if ((accept2/irep)>0.3)  scale2 <- 1.01*scale2
        if ((accept4/irep)<0.15) scale4 <- 0.99*scale4
        if ((accept4/irep)>0.3)  scale4 <- 1.01*scale4
      }
    }
    # SSVS
    if(prior==2){
      for(mm in 1:M){
        for(kk in 1:k){
          u_i1  <-  dnorm(A_draw[kk,mm],A_prior[kk,mm],tau0[kk,mm]) * p_i
          u_i2  <-  dnorm(A_draw[kk,mm],A_prior[kk,mm],tau1[kk,mm]) * (1-p_i)
          gst  <-  u_i1/(u_i1 + u_i2)
          if(gst=="NaN") gst <- 0
          gamma[kk,mm]  <-  .bernoulli(gst)
          gamma[is.na(gamma)] <- 1
          if (gamma[kk,mm] == 0){
            theta[kk,mm]  <-  tau0[kk,mm]^2
          }else if (gamma[kk,mm] == 1){
            theta[kk,mm]  <-  tau1[kk,mm]^2
          }
        }
      }
      for(kk in 1:Ki1){
        for(mm in 2:M){
          for(ii in 1:(mm-1)){
            u_ij1  <-  dnorm(J_draw[mm,ii,kk],J_prior[mm,ii,kk],kappa0[mm,ii,kk]) * q_ij
            u_ij2  <-  dnorm(J_draw[mm,ii,kk],J_prior[mm,ii,kk],kappa1[mm,ii,kk]) * (1-q_ij)
            ost  <-  u_ij1/(u_ij1 + u_ij2)
            if(is.na(ost)) ost <- 1
            omega[mm,ii,kk] <-  .bernoulli(ost)
            if (is.na(omega[mm,ii,kk])) omega[mm,ii,kk] <- 1
            if(omega[mm,ii,kk]==1){
              zeta[mm,ii,kk] <- kappa1[mm,ii,kk]^2
            }else{
              zeta[mm,ii,kk] <- kappa0[mm,ii,kk]^2
            }
          }
        }
      }
    }
    # NG
    if(prior==3){
      # Normal-Gamma for Covariances
      summand <- 0
      for(kk in 1:Ki1) for(mm in 2:M) for(ii in 1:(mm-1)) summand<-summand+zeta[mm,ii,kk]
      lambda2_J    <- rgamma(1,d_lambda+J_tau*v,e_lambda+J_tau/2*summand)
      #Step VI: Sample the prior scaling factors for covariances from GIG
      zetas <- NULL
      for(kk in 1:Ki1){
        for(mm in 2:M){
          for(ii in 1:(mm-1)){
            temp           <- do_rgig1(lambda=J_tau-0.5, chi=(J_draw[mm,ii,kk]-J_prior[mm,ii,kk])^2, psi=J_tau*lambda2_J)
            zetas          <- c(zetas,ifelse(temp<1e-8,1e-8,temp))
            zeta[mm,ii,kk] <- temp
          }
        }
      }
      if(sample_A){
        #Sample L_tau through a simple RWMH step
        J_tau_prop       <- exp(rnorm(1,0,J_tuning))*J_tau
        post_J_tau_prop  <- .atau_post(atau=J_tau_prop, thetas=zetas, k=Ki1*v, lambda2=lambda2_J)
        post_J_tau_old   <- .atau_post(atau=J_tau,      thetas=zetas, k=Ki1*v, lambda2=lambda2_J)
        post.diff    <- post_J_tau_prop-post_J_tau_old
        post.diff    <- ifelse(is.nan(post.diff),-Inf,post.diff)
        if (post.diff > log(runif(1,0,1))){
          J_tau      <- J_tau_prop
          J_accept   <- J_accept+1
        }
        if (irep<(0.5*nburn)){
          if ((J_accept/irep)>0.3)  J_tuning <- 1.01*J_tuning
          if ((J_accept/irep)<0.15) J_tuning <- 0.99*J_tuning
        }
      }
      # Normal-Gamma for endogenous variables
      for (ss in 1:p){
        slct.i    <- grepl(paste0(".lag",ss),rownames(A_draw))
        A.lag     <- A_draw[slct.i,,drop=FALSE]
        A.prior   <- A_prior[slct.i,,drop=FALSE]
        theta.lag <- theta[slct.i,,drop=FALSE]

        M.end <- nrow(A.lag)
        if (ss==1){
          lambda2_A[ss,1] <- rgamma(1,d_lambda+A_tau[ss,1]*M.end^2,e_lambda+A_tau[ss,1]/2*sum(theta.lag))
        }else{
          lambda2_A[ss,1] <- rgamma(1,d_lambda+A_tau[ss,1]*M.end^2,e_lambda+A_tau[ss,1]/2*prod(lambda2_A[1:(ss-1),1])*sum(theta.lag))
        }
        for(jj in 1:M.end){
          for (ii in 1:M){
            theta.lag[jj,ii] <- do_rgig1(lambda=A_tau[ss,1]-0.5,
                                         chi=(A.lag[jj,ii]-A.prior[jj,ii])^2,
                                         psi=A_tau[ss,1]*prod(lambda2_A[1:ss,1]))
          }
        }
        theta.lag[theta.lag<1e-8] <- 1e-8
        theta[slct.i,] <- theta.lag
        #TO BE MODIFIED
        if (sample_A){
          #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
          A_tau_prop <- exp(rnorm(1,0,A_tuning[ss,1]))*A_tau[ss,1]
          post_A_tau_prop <- .atau_post(atau=A_tau_prop,  thetas=as.vector(theta.lag), lambda2=prod(lambda2_A[1:ss,1]), k=length(theta.lag))
          post_A_tau_old  <- .atau_post(atau=A_tau[ss,1], thetas=as.vector(theta.lag), lambda2=prod(lambda2_A[1:ss,1]), k=length(theta.lag))
          post.diff <- post_A_tau_prop-post_A_tau_old
          post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)

          if (post.diff > log(runif(1,0,1))){
            A_tau[ss,1]    <- A_tau_prop
            A_accept[ss,1] <- A_accept[ss,1]+1
          }
          if (irep<(0.5*nburn)){
            if ((A_accept[ss,1]/irep)>0.3)  A_tuning[ss,1] <- 1.01*A_tuning[ss,1]
            if ((A_accept[ss,1]/irep)<0.15) A_tuning[ss,1] <- 0.99*A_tuning[ss,1]
          }
        }
      }
    }
    #----------------------------------------------------------------------------
    # Step 3: Sample variances
    if (sv){
      for (jj in 1:M){
        para   <- as.list(pars_var[,jj])
        para$nu = Inf; para$rho=0; para$beta<-0
        svdraw <- svsample_fast_cpp(y=Em[,jj], draws=1, burnin=0, designmatrix=matrix(NA_real_),
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
    #----------------------------------------------------------------------------
    # Step 4: store draws
    if(irep %in% thin.draws){
      count <- count+1
      A_store[count,,]      <- A_draw
      J_store[count,,,]     <- J_draw
      # SV
      Sv_store[count,,]     <- Sv_draw
      pars_store[count,,]   <- pars_var
      # MN
      shrink_store[count,]  <- c(shrink1,shrink2,shrink4)
      # SSVS
      gamma_store[count,,]  <- gamma
      omega_store[count,,,] <- omega
      # NG
      theta_store[count,,]     <- theta
      zeta_store[count,,,]     <- zeta
      lambda2_store[count,1,2] <- lambda2_J
      lambda2_store[count,,1]  <- lambda2_A
      tau_store[count,1,2]     <- J_tau
      tau_store[count,,1]      <- A_tau
    }
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,colnames(X),colnames(A_OLS))
  ret <- list(Y=Y,X=X,D=D,A_store=A_store,J_store=J_store,Sv_store=Sv_store,shrink_store=shrink_store,gamma_store=gamma_store,omega_store=omega_store,theta_store=theta_store,zeta_store=zeta_store,lambda2_store=lambda2_store,tau_store=tau_store,pars_store=pars_store)
  return(ret)
}

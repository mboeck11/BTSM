#' @name .BVAR_linear_wrapper
#' @noRd
#' @importFrom abind adrop
#' @importFrom utils capture.output
.BVAR_linear_wrapper <- function(Yraw, prior, plag, draws, burnin, cons, trend, SV, thin, default_hyperpara, Ex, applyfun, cores){
  class(Yraw) <- "numeric"
  prior_in <- ifelse(prior=="MN",1,ifelse(prior=="SSVS",2,ifelse(prior=="NG",3,NA)))
  if(default_hyperpara[["a_log"]]){
    default_hyperpara["a_start"] <- 1/log(ncol(Yraw))
  }
  if(!is.na(prior_in)){
    invisible(capture.output(
      bvar<-BVAR_linear(Y_in=Yraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)
      ,type="message"))
    if(is(bvar,"try-error")){
      bvar<-.BVAR_linear_R(Y_in=Yraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)
    }
  }else if(prior=="NC"){
    bvar <- bvar_natural_conjugate(Y_in=Yraw,p_in=plag,draws_in=draws,cons_in=cons,trend_in=trend,thin_in=thin,hyperparam_in=default_hyperpara,Ex_in=Ex,applyfun=applyfun,cores=cores)
  }
  #------------------------------------------------ get data ----------------------------------------#
  Y <- bvar$Y; colnames(Y) <- colnames(Yraw); X <- bvar$X
  M <- ncol(Y); bigT <- nrow(Y); K <- ncol(X)
  if(!is.null(Ex)) Mex <- ncol(Ex)
  xnames <- paste(rep("Ylag",M),rep(seq(1,plag),each=M),sep="")
  if(!is.null(Ex)) xnames <- c(xnames,paste(rep("Tex",Mex)))
  if(cons)  xnames <- c(xnames,"cons")
  if(trend) xnames <- c(xnames,"trend")
  colnames(X) <- xnames
  #-----------------------------------------get containers ------------------------------------------#
  A_store <- bvar$A_store; dimnames(A_store)[[2]] <- colnames(X); dimnames(A_store)[[3]] <- colnames(Y)
  # splitting up stores
  dims          <- dimnames(A_store)[[2]]
  a0store <- a1store <- Exstore <- NULL
  if(cons) {
    a0store       <- adrop(A_store[,which(dims=="cons"),,drop=FALSE],drop=2)
  }
  if(trend){
    a1store     <- adrop(A_store[,which(dims=="trend"),,drop=FALSE],drop=2)
  }
  if(!is.null(Ex)){
    Exstore     <- A_store[,which(dims=="Tex"),,drop=FALSE]
  }
  Phistore    <- NULL
  for(jj in 1:plag){
    Phistore[[jj]]  <- A_store[,which(dims==paste("Ylag",jj,sep="")),,drop=FALSE]
  }
  S_store <- array(NA, c(draws/thin,bigT,M,M)); dimnames(S_store) <- list(NULL,NULL,colnames(Y),colnames(Y))
  if(prior%in%c("MN","SSVS","NG")){
    L_store <- bvar$L_store
    for(irep in 1:(draws/thin)){
      for(tt in 1:bigT){
        if(M>1){
          S_store[irep,tt,,] <- L_store[irep,,]%*%diag(exp(bvar$Sv_store[irep,tt,]))%*%t(L_store[irep,,])
        }else{
          S_store[irep,tt,,] <- L_store[irep,,]%*%exp(bvar$Sv_store[irep,tt,])%*%t(L_store[irep,,])
        }
      }
    }
    Smed_store <- apply(S_store,c(1,3,4),median)
    if(SV){
      vola_store  <- bvar$Sv_store; dimnames(vola_store) <- list(NULL,NULL,colnames(Y))
      pars_store  <- bvar$pars_store
      vola_post   <- apply(vola_store,c(2,3),median)
      pars_post   <- apply(pars_store,c(2,3),median)
    }else{
      vola_store  <- bvar$Sv_store; pars_store <- NULL;
      vola_post   <- apply(vola_store,c(2,3),median); pars_post <- NULL
    }
  }else if(prior=="NC"){
    Smed_store  <- bvar$S_store
    for(irep in 1:(draws/thin)){
      for(tt in 1:bigT){
        S_store[irep,tt,,] <- bvar$S_store[irep,,]
      }
    }
    L_store     <- NULL
    theta_store <- NULL
    vola_store  <- NULL
    pars_store  <- NULL
    vola_post   <- NULL
    pars_post   <- NULL
  }
  theta_store   <- bvar$theta_store; dimnames(theta_store)[[2]] <- colnames(X); dimnames(theta_store)[[3]] <- colnames(Y)
  res_store     <- bvar$res_store; dimnames(res_store) <- list(NULL,NULL,colnames(Y))
  # MN
  if(prior=="MN"){
    shrink_store  <- bvar$shrink_store; dimnames(shrink_store) <- list(NULL,c("shrink1","shrink2","shrink4"))
    shrink_post   <- apply(shrink_store,2,median)
  }else{
    shrink_store  <- shrink_post <- NULL
  }
  # SSVS
  if(prior=="SSVS"){
    gamma_store <- bvar$gamma_store; dimnames(gamma_store) <- list(NULL,colnames(X),colnames(Y))
    omega_store <- bvar$omega_store; dimnames(omega_store) <- list(NULL,colnames(Y),colnames(Y))
    PIP         <- apply(gamma_store,c(2,3),mean)
    PIP_omega   <- apply(omega_store,c(2,3),mean)
  }else{
    gamma_store <- omega_store <- PIP <- PIP_omega <- NULL
  }
  # NG
  if(prior=="NG"){
    lambda2_store <- bvar$lambda2_store
    tau_store     <- bvar$tau_store
    dimnames(lambda2_store) <- list(NULL,paste("lag",1:plag,sep="_"),c("endogenous","covariance"))
    dimnames(lambda2_store) <- list(NULL,paste("lag",1:plag,sep="_"),c("endogenous","covariance"))
    lambda2_post  <- apply(lambda2_store,c(2,3),median)
    tau_post      <- apply(tau_store,c(2,3),median)
  }else{
    lambda2_store <- tau_store <- lambda2_post <- tau_post <- NULL
  }
  store <- list(A_store=A_store,a0store=a0store,a1store=a1store,Phistore=Phistore,Exstore=Exstore,S_store=S_store,Smed_store=Smed_store,
                L_store=L_store,theta_store=theta_store,vola_store=vola_store,pars_store=pars_store,res_store=res_store,
                shrink_store=shrink_store,gamma_store=gamma_store,omega_store=omega_store,lambda2_store=lambda2_store,tau_store=tau_store)
  #------------------------------------ compute posteriors -------------------------------------------#
  A_post      <- apply(A_store,c(2,3),median)
  S_post      <- apply(S_store,c(2,3,4),median)
  Sig         <- apply(S_post,c(2,3),mean)/(bigT-K)
  theta_post  <- apply(theta_store,c(2,3),median)
  res_post    <- apply(res_store,c(2,3),median)
  # splitting up posteriors
  a0post <- a1post <- Expost <- NULL
  if(cons)  a0post <- A_post[which(dims=="cons"),,drop=FALSE]
  if(trend) a1post <- A_post[which(dims=="trend"),,drop=FALSE]
  if(!is.null(Ex)) Expost <- A_post[which(dims=="Tex"),,drop=FALSE]
  Phipost     <- NULL
  for(jj in 1:plag){
    Phipost    <- rbind(Phipost,A_post[which(dims==paste("Ylag",jj,sep="")),,drop=FALSE])
  }
  post <- list(A_post=A_post,a0post=a0post,a1post=a1post,Phipost=Phipost,Expost=Expost,S_post=S_post,Sig=Sig,theta_post=theta_post,
               vola_post=vola_post,pars_post=pars_post,res_post=res_post,shrink_post=shrink_post,PIP=PIP,PIP_omega=PIP_omega,
               lambda2_post=lambda2_post,tau_post=tau_post)
  return(list(Y=Y,X=X,store=store,post=post))
}

#' @name .BVAR_linear_R
#' @importFrom stochvol svsample_fast_cpp specify_priors default_fast_sv
#' @importFrom MASS ginv mvrnorm
#' @importFrom methods is
#' @importFrom stats rnorm rgamma runif dnorm
#' @noRd
.BVAR_linear_R <- function(Y_in,p_in,draws_in,burnin_in,cons_in,trend_in,sv_in,thin_in,quiet_in,prior_in,hyperparam_in,Ex_in){
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
  crit_eig  <- hyperpara$crit_eig
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
  kappa0    <- hyperpara$kappa0
  kappa1    <- hyperpara$kappa1
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
  A_OLS  <- XtXinv%*%(t(X)%*%Y)
  E_OLS  <- Y - X%*%A_OLS
  #a_OLS <- as.vector(A_OLS)
  #SSE  <-  t((Y - X%*%A_OLS))%*%(Y - X%*%A_OLS)
  SIGMA_OLS  <- crossprod(E_OLS)/(bigT-k)
  #IXY  <-   kronecker(diag(M),(t(X)%*%Y))
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  A_draw <- A_OLS
  SIGMA  <- array(SIGMA_OLS, c(M,M,bigT))
  Em     <- Em_str <- E_OLS
  L_draw <- diag(M)
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  # Priors on VAR coefs
  #-----------------------------
  # prior mean
  A_prior <- matrix(0,k,M)
  diag(A_prior) <- prmean
  a_prior  <-  as.vector(A_prior)
  # prior variance
  theta <- matrix(10,k,M)

  # MN stuff
  accept1 <- 0
  accept2 <- 0
  accept4 <- 0
  scale1  <- .43
  scale2  <- .43
  scale4  <- .43
  sigma_sq  <- matrix(0,M,1) #vector which stores the residual variance
  for (i in 1:M){
    Ylag_i         <- .mlag(Yraw[,i],p)
    Ylag_i         <- Ylag_i[(p+1):nrow(Ylag_i),,drop=FALSE]
    Y_i            <- Yraw[(p+1):nrow(Yraw),i,drop=FALSE]
    Ylag_i         <- cbind(Ylag_i,seq(1,nrow(Y_i)))
    alpha_i        <- solve(crossprod(Ylag_i))%*%crossprod(Ylag_i,Y_i)
    sigma_sq[i,1]  <- (1/(nrow(Y_i)-p-1))*t(Y_i-Ylag_i%*%alpha_i)%*%(Y_i-Ylag_i%*%alpha_i)
  }
  if(prior==1){
    theta <- .get_V(k=k,M=M,p=p,a_bar_1=shrink1,a_bar_2=shrink2,a_bar_3=shrink3,
                    a_bar_4=shrink4,sigma_sq=sigma_sq,trend=trend)
  }

  # SSVS stuff
  gamma  <-  matrix(1,k,M)
  sigma_alpha  <-  sqrt(diag(kronecker(SIGMA_OLS,XtXinv)))
  tau0 <- matrix(NA, k, M); tau1 <- matrix(NA, k, M)
  ii <- 1
  for(mm in 1:M){
    for(kk in 1:k){
      tau0[kk,mm] <- tau00*sigma_alpha[ii]
      tau1[kk,mm] <- tau11*sigma_alpha[ii]
      ii <- ii+1
    }
  }
  # NG stuff
  lambda2_A    <- matrix(0.01,p,1)
  A_tau        <- matrix(a_start,p,1)
  colnames(A_tau) <- colnames(lambda2_A) <- "endo"
  rownames(A_tau) <- rownames(lambda2_A) <- paste("lag.",seq(1,p),sep="")
  A_tuning     <- matrix(.43,p,1)
  A_accept     <- matrix(0,p,1)
  #------------------------------------
  # Priors on coefs in H matrix of VCV
  #------------------------------------
  # prior mean
  l_prior <- matrix(0,M,M)

  # prior variance
  L_prior <- matrix(10,M,M)
  L_prior[upper.tri(L_prior)] <- 0; diag(L_prior) <- 0

  # SSVS
  omega <- matrix(1,M,M)
  omega[upper.tri(omega)] <- 0; diag(omega) <- 0

  # NG
  lambda2_L <- 0.01
  L_tau     <- a_start
  L_accept <- 0
  L_tuning  <- .43
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
  A_store      <- array(NA,c(thindraws,k,M))
  L_store      <- array(NA,c(thindraws,M,M))
  res_store    <- array(NA,c(thindraws,bigT,M))
  # SV
  Sv_store     <- array(NA,c(thindraws,bigT,M))
  pars_store   <- array(NA,c(thindraws,4,M))
  # MN
  shrink_store <- array(NA,c(thindraws,3))
  # SSVS
  gamma_store  <- array(NA,c(thindraws,k,M))
  omega_store  <- array(NA,c(thindraws,M,M))
  # NG
  theta_store  <- array(NA,c(thindraws,k,M))
  lambda2_store<- array(NA,c(thindraws,p,2))
  tau_store    <- array(NA,c(thindraws,p,2))
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Sample coefficients
    for (mm in 1:M){
      if (mm==1){
        Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
        X.i <- X*exp(-0.5*Sv_draw[,mm])

        V_post <- try(chol2inv(chol(crossprod(X.i)+diag(1/theta[,mm]))),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv(crossprod(X.i)+diag(1/theta[,mm]))
        A_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/theta[,mm])%*%A_prior[,mm])

        A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
        if (is(A.draw.i,"try-error")) A.draw.i <- mvrnorm(1,A_post,V_post)
        A_draw[,mm] <- A.draw.i
        Em[,mm] <-  Em_str[,mm] <- Y[,mm]-X%*%A.draw.i
      }else{
        Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
        X.i <- cbind(X,Em[,1:(mm-1)])*exp(-0.5*Sv_draw[,mm])

        V_post <- try(chol2inv(chol((crossprod(X.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))))),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv((crossprod(X.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))))
        A_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))%*%c(A_prior[,mm],l_prior[mm,1:(mm-1)]))

        A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
        if (is(A.draw.i,"try-error")) A.draw.i <- mvrnorm(1,A_post,V_post)

        A_draw[,mm] <- A.draw.i[1:ncol(X)]
        Em[,mm] <- Y[,mm]-X%*%A.draw.i[1:ncol(X)]
        Em_str[,mm] <- Y[,mm]-X%*%A.draw.i[1:ncol(X)]-Em[,1:(mm-1),drop=FALSE]%*%A.draw.i[(ncol(X)+1):ncol(X.i),drop=FALSE]
        L_draw[mm,1:(mm-1)] <- A.draw.i[(ncol(X)+1):ncol(X.i)]
      }
    }
    rownames(A_draw) <- colnames(X)
    #----------------------------------------------------------------------------
    # Step 2: different shrinkage prior setups
    # MN
    if(prior==1){
      #Step for the first shrinkage parameter (own lags)
      shrink1.prop <- exp(rnorm(1,0,scale1))*shrink1
      if(shrink1.prop<1e-17) shrink1.prop <- 1e-17
      if(shrink1.prop>1e+17) shrink1.prop <- 1e+17
      theta1.prop   <- .get_V(k=k,M=M,p=p,a_bar_1=shrink1.prop,a_bar_2=shrink2,a_bar_3=shrink3,a_bar_4=shrink4,sigma_sq=sigma_sq)
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
      theta2.prop   <- .get_V(k=k,M=M,p=p,a_bar_1=shrink1,a_bar_2=shrink2.prop,a_bar_3=shrink3,a_bar_4=shrink4,sigma_sq=sigma_sq)
      post2.prop <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta2.prop)),log=TRUE))+dgamma(shrink2.prop,0.01,0.01,log=TRUE)
      post2.prop <- post2.prop + log(shrink2.prop) # correction term
      post2 <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta)),log=TRUE))+dgamma(shrink2,0.01,0.01,log=TRUE)
      post2 <- post2 + log(shrink2) # correction term
      if ((post2.prop-post2)>log(runif(1,0,1))){
        shrink2 <- shrink2.prop
        theta   <- theta2.prop
        accept2 <- accept2+1
      }

      #Step for the final shrinkage parameter (weakly exogenous)
      shrink4.prop <- exp(rnorm(1,0,scale4))*shrink4
      if(shrink4.prop<1e-17) shrink4.prop <- 1e-17
      if(shrink4.prop>1e+17) shrink4.prop <- 1e+17
      theta4.prop   <- .get_V(k=k,M=M,p=p,a_bar_1=shrink1,a_bar_2=shrink2,a_bar_3=shrink3,a_bar_4=shrink4.prop,sigma_sq=sigma_sq)
      post4.prop <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta4.prop)),log=TRUE))+dgamma(shrink4.prop,0.01,0.01,log=TRUE)
      post4.prop <- post4.prop + log(shrink4.prop)
      post4 <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta)),log=TRUE))+dgamma(shrink4,0.01,0.01,log=TRUE)
      post4 <- post4 + log(shrink4)
      if ((post4.prop-post4)>log(runif(1,0,1))){
        shrink4  <- shrink4.prop
        theta    <- theta4.prop
        accept4  <- accept4+1
      }

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
      for(mm in 2:M){
        for(ii in 1:(mm-1)){
          u_ij1  <-  dnorm(L_draw[mm,ii],l_prior[mm,ii],kappa0) * q_ij
          u_ij2  <-  dnorm(L_draw[mm,ii],l_prior[mm,ii],kappa1) * (1-q_ij)
          ost  <-  u_ij1/(u_ij1 + u_ij2)
          if(is.na(ost)) ost <- 1
          omega[mm,ii] <-  .bernoulli(ost)
          if (is.na(omega[mm,ii])) omega[mm,ii] <- 1
          if(omega[mm,ii]==1){
            L_prior[mm,ii] <- kappa1^2
          }else{
            L_prior[mm,ii] <- kappa0^2
          }
        }
      }
    }
    # NG
    if(prior==3){
      # Normal-Gamma for Covariances
      lambda2_L    <- rgamma(1,d_lambda+L_tau*v,e_lambda+L_tau/2*sum(L_prior[lower.tri(L_prior)]))
      #Step VI: Sample the prior scaling factors for covariances from GIG
      for(mm in 2:M){
        for(ii in 1:(mm-1)){
          L_prior[mm,ii] <- do_rgig1(lambda=L_tau-0.5, chi=(L_draw[mm,ii]-l_prior[mm,ii])^2, psi=L_tau*lambda2_L)
        }
      }
      if(sample_A){
        #Sample L_tau through a simple RWMH step
        L_tau_prop       <- exp(rnorm(1,0,L_tuning))*L_tau
        post_L_tau_prop  <- .atau_post(atau=L_tau_prop, thetas=L_prior[lower.tri(L_prior)], k=v, lambda2=lambda2_L)
        post_L_tau_old   <- .atau_post(atau=L_tau,      thetas=L_prior[lower.tri(L_prior)], k=v, lambda2=lambda2_L)
        post.diff    <- post_L_tau_prop-post_L_tau_old
        post.diff    <- ifelse(is.nan(post.diff),-Inf,post.diff)
        if (post.diff > log(runif(1,0,1))){
          L_tau      <- L_tau_prop
          L_accept   <- L_accept+1
        }
        if (irep<(0.5*nburn)){
          if ((L_accept/irep)>0.3)  L_tuning <- 1.01*L_tuning
          if ((L_accept/irep)<0.15) L_tuning <- 0.99*L_tuning
        }
      }
      # Normal-Gamma for endogenous variables
      for (ss in 1:p){
        slct.i    <- which(rownames(A_draw)==paste("Ylag",ss,sep=""))
        A.lag     <- A_draw[slct.i,,drop=FALSE]
        A.prior   <- A_prior[slct.i,,drop=FALSE]
        theta.lag <- theta[slct.i,,drop=FALSE]

        M.end <- nrow(A.lag)
        if (ss==1){
          lambda2_A[ss,1] <- rgamma(1,d_lambda+A_tau[ss,1]*M^2,e_lambda+A_tau[ss,1]/2*sum(theta.lag))
        }else{
          lambda2_A[ss,1] <- rgamma(1,d_lambda+A_tau[ss,1]*M^2,e_lambda+A_tau[ss,1]/2*prod(lambda2_A[1:(ss-1),1])*sum(theta.lag))
        }
        for (jj in 1:M){
          for (ii in 1:M){
            theta.lag[jj,ii] <- do_rgig1(lambda=A_tau[ss,1]-0.5,
                                         chi=(A.lag[jj,ii]-A.prior[jj,ii])^2,
                                         psi=A_tau[ss,1]*prod(lambda2_A[1:ss,1]))
          }
        }
        theta[slct.i,] <- theta.lag
        theta[theta<1e-8] <- 1e-8
        #TO BE MODIFIED
        if (sample_A){
          #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
          A_tau_prop <- exp(rnorm(1,0,A_tuning[ss,1]))*A_tau[ss,1]
          post_A_tau_prop <- .atau_post(atau=A_tau_prop,  thetas=as.vector(theta.lag), lambda2=prod(lambda2_A[1:ss,1]), k=length(theta.lag))
          post_A_tau_old  <- .atau_post(atau=A_tau[ss,1], thetas=as.vector(theta.lag), lambda2=prod(lambda2_A[1:ss,1]), k=length(theta.lag))
          post.diff <- post_A_tau_prop-post_A_tau_old
          post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)

          if (post.diff > log(runif(1,0,1))){
            A_tau[ss,1] <- A_tau_prop
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
        svdraw <- svsample_fast_cpp(y=Em_str[,jj], draws=1, burnin=0, designmatrix=matrix(NA_real_),
                                    priorspec=specify_priors(), thinpara=1, thinlatent=1, keeptime="all",
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
      A_store[count,,] <- A_draw
      L_store[count,,] <- L_draw
      res_store[count,,] <- Y-X%*%A_draw
      # SV
      Sv_store[count,,] <- Sv_draw
      pars_store[count,,] <- pars_var
      # MN
      shrink_store[count,] <- c(shrink1,shrink2,shrink4)
      # SSVS
      gamma_store[count,,] <- gamma
      omega_store[count,,] <- omega
      # NG
      theta_store[count,,]     <- theta
      lambda2_store[count,1,2] <- lambda2_L
      lambda2_store[count,,1]  <- lambda2_A
      tau_store[count,1,2]     <- L_tau
      tau_store[count,,1]      <- A_tau
    }
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,colnames(X),colnames(A_OLS))
  ret <- list(Y=Y,X=X,A_store=A_store,L_store=L_store,Sv_store=Sv_store,shrink_store=shrink_store,gamma_store=gamma_store,omega_store=omega_store,theta_store=theta_store,lambda2_store=lambda2_store,tau_store=tau_store,pars_store=pars_store,res_store=res_store)
  return(ret)
}

#' @name .BVAR_natural_conjugate
#' @importFrom stats rgamma rnorm quantile rWishart
#' @importFrom mvnfast dmvn
#' @importFrom mvtnorm rmvnorm dmvnorm rmvt
#' @importFrom MASS ginv
#' @importFrom methods is
#' @noRd
bvar_natural_conjugate <- function(Y_in,p_in,draws_in,cons_in,trend_in,thin_in,quiet_in,hyperparam_in,Ex_in,applyfun,cores) {
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
  c         <- hyperpara$c
  prmean    <- hyperpara$prmean
  Multiplier<- hyperpara$Multiplier
  #--------------------------Initialize Gibbs sampler--------------------------------#
  A_OLS  <- try(solve(crossprod(X)) %*% crossprod(X,Y),silent=TRUE)
  if(is(A_OLS,"try-error")) A_OLS <- ginv(crossprod(X)) %*% crossprod(X,Y)
  S_OLS <- crossprod(Y - X %*% A_OLS)
  #----------------------------PRIORS------------------------------------------------#
  # prior mean for autoregressive parameters
  A_prior <- matrix(0,k,M)
  A_prior[1:M,1:M] <- diag(M)*prmean
  # prior variance for autoregressive coefficients
  theta_prior <- diag(k)*c
  theta_priorinv <- diag(1/diag(theta_prior))
  # prior degrees of scaling
  v_prior <- 1
  # prior scaling matrix
  S_prior <- (1/c)*diag(M)
  #---------------------POSTERIOR MOMENTS--------------------------------------------#
  # posterior of coefficients
  V_post <- solve(crossprod(X) + theta_priorinv)
  # A_post <- V_post %*% (crossprod(X)%*%A_OLS + V_priorinv%*%A_prior)
  A_post <- V_post %*% (crossprod(X,Y) + theta_priorinv%*%A_prior)
  # posterior of variance
  S_post <- S_OLS + S_prior + t(X%*%A_OLS)%*%X%*%A_OLS + t(A_prior)%*%theta_priorinv%*%A_prior - t(A_post)%*%(theta_priorinv + crossprod(X))%*%A_post
  v_post <- v_prior + bigT
  # posterior of coefficient variance for t-distribution
  bigVpost  <- kronecker(S_post, V_post)/(v_post-M-1)
  #---------------------------------------------------------------------------------------------------------
  # SAMPLER MISCELLANEOUS
  #---------------------------------------------------------------------------------------------------------
  draws        <- draws_in
  # thinning
  thin         <- thin_in
  count        <- 0
  thindraws    <- draws/thin
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store      <- array(NA,c(thindraws,k,M))
  theta_store  <- array(NA,c(thindraws,k,M))
  res_store    <- array(NA,c(thindraws,bigT,M))
  S_store      <- array(NA,c(thindraws,M,M))
  #-------------------MONTE CARLO SIMULATION----------------------------------------#
  storage <- applyfun(1:thindraws,function(irep){
    # draw coefficients
    Sinv_draw <- matrix(rWishart(1,v_post,solve(S_post)),M,M)
    S_draw    <- solve(Sinv_draw)
    A_draw    <- matrix(mvtnorm::rmvt(1, sigma=bigVpost, df=v_post, delta=as.vector(A_post)),k,M)
    res_draw  <- Y - X%*%A_draw
    return(list(A_draw=A_draw,S_draw=S_draw,res_draw=res_draw))
  })
  # save everything
  for(irep in 1:thindraws){
    A_store[irep,,]     <- storage[[irep]]$A_draw
    S_store[irep,,]     <- storage[[irep]]$S_draw
    res_store[irep,,]   <- storage[[irep]]$res
    theta_store[irep,,] <- matrix(rep(diag(theta_prior),M),k,M)
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,colnames(X),colnames(A_OLS))
  ret <- list(Y=Y,X=X,A_store=A_store,S_store=S_store,theta_store=theta_store,res_store=res_store)
  return(ret)
}

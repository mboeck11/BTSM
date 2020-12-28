#' @name .BVAR_linear_wrapper
#' @noRd
#' @importFrom abind adrop
#' @importFrom utils capture.output
.BVEC_linear_wrapper <- function(Yraw, r, beta, prior, plag, draws, burnin, cons, trend, SV, thin, default_hyperpara, Ex){
  class(Yraw) <- "numeric"
  prior_in <- ifelse(prior=="MN",1,ifelse(prior=="SSVS",2,ifelse(prior=="NG",3,NA)))
  if(default_hyperpara[["a_log"]]){
    default_hyperpara["a_start"] <- 1/log(ncol(Yraw))
  }
  if(!is.na(prior_in)){
    bvec<-.BVEC_linear_R(Y_in=Yraw,r_in=r,beta_in=beta,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)
  }
  #------------------------------------------------ get data ----------------------------------------#
  Y <- bvec$Y; X <- bvec$X; Z <- bvec$Z
  M <- ncol(Y); bigT <- nrow(Y); K <- ncol(X)
  if(!is.null(Ex)) Mex <- ncol(Ex)
  xnames <- paste(rep("Ylag",M),rep(seq(1,plag),each=M),sep="")
  if(!is.null(Ex)) xnames <- c(xnames,paste(rep("Tex",Mex)))
  if(cons)  xnames <- c(xnames,"cons")
  if(trend) xnames <- c(xnames,"trend")
  colnames(X) <- xnames
  znames <- paste("ECM",seq(1,r),sep="")
  #-----------------------------------------get containers ------------------------------------------#
  PHI_store <- bvec$PHI_store; dimnames(PHI_store)[[2]] <- c(znames,xnames); dimnames(PHI_store)[[3]] <- colnames(Y)
  beta_store <- bvec$beta_store; dimnames(beta_store)[[2]] <- colnames(Y); dimnames(beta_store)[[3]] <- paste0("ECM",seq(1,r))
  # splitting up stores
  dims          <- dimnames(PHI_store)[[2]]
  alpha_store   <- PHI_store[,which(dims==paste0("ECM",seq(1,r))),,drop=FALSE]
  a0store <- a1store <- Exstore <- NULL
  if(cons) {
    a0store       <- adrop(PHI_store[,which(dims=="cons"),,drop=FALSE],drop=2)
  }
  if(trend){
    a1store     <- adrop(PHI_store[,which(dims=="trend"),,drop=FALSE],drop=2)
  }
  if(!is.null(Ex)){
    Exstore     <- PHI_store[,which(dims=="Tex"),,drop=FALSE]
  }
  Phistore    <- NULL
  for(jj in 1:plag){
    Phistore[[jj]]  <- PHI_store[,which(dims==paste("Ylag",jj,sep="")),,drop=FALSE]
  }
  S_store <- array(NA, c(draws/thin,bigT,M,M)); dimnames(S_store) <- list(NULL,NULL,colnames(Y),colnames(Y))
  if(prior%in%c("MN","SSVS","NG")){
    L_store <- bvec$L_store
    for(irep in 1:(draws/thin)){
      for(tt in 1:bigT){
        if(M>1){
          S_store[irep,tt,,] <- L_store[irep,,]%*%diag(exp(bvec$Sv_store[irep,tt,]))%*%t(L_store[irep,,])
        }else{
          S_store[irep,tt,,] <- L_store[irep,,]%*%exp(bvec$Sv_store[irep,tt,])%*%t(L_store[irep,,])
        }
      }
    }
    Smed_store <- apply(S_store,c(1,3,4),median)
    if(SV){
      vola_store  <- bvec$Sv_store; dimnames(vola_store) <- list(NULL,NULL,colnames(Y))
      pars_store  <- bvec$pars_store
      vola_post   <- apply(vola_store,c(2,3),median)
      pars_post   <- apply(pars_store,c(2,3),median)
    }else{
      vola_store  <- bvec$Sv_store; pars_store <- NULL;
      vola_post   <- apply(vola_store,c(2,3),median); pars_post <- NULL
    }
  }
  theta_store   <- bvec$theta_store; dimnames(theta_store)[[2]] <- c(znames,xnames); dimnames(theta_store)[[3]] <- colnames(Y)
  res_store     <- bvec$res_store; dimnames(res_store) <- list(NULL,NULL,colnames(Y))
  # MN
  if(prior=="MN"){
    shrink_store  <- bvec$shrink_store; dimnames(shrink_store) <- list(NULL,c("shrink1","shrink2","shrink4"))
    shrink_post   <- apply(shrink_store,2,median)
  }else{
    shrink_store  <- shrink_post <- NULL
  }
  # SSVS
  if(prior=="SSVS"){
    gamma_store <- bvec$gamma_store; dimnames(gamma_store) <- list(NULL,c(znames,xnames),colnames(Y))
    omega_store <- bvec$omega_store; dimnames(omega_store) <- list(NULL,colnames(Y),colnames(Y))
    PIP         <- apply(gamma_store,c(2,3),mean)
    PIP_omega   <- apply(omega_store,c(2,3),mean)
  }else{
    gamma_store <- omega_store <- PIP <- PIP_omega <- NULL
  }
  # NG
  if(prior=="NG"){
    lambda2_store <- bvec$lambda2_store
    tau_store     <- bvec$tau_store
    dimnames(lambda2_store) <- list(NULL,paste("lag",1:plag,sep="_"),c("endogenous","covariance","beta"))
    dimnames(lambda2_store) <- list(NULL,paste("lag",1:plag,sep="_"),c("endogenous","covariance","beta"))
    lambda2_post  <- apply(lambda2_store,c(2,3),median)
    tau_post      <- apply(tau_store,c(2,3),median)
  }else{
    lambda2_store <- tau_store <- lambda2_post <- tau_post <- NULL
  }
  Astore <- list()
  for(jj in (plag+1):1){
    Astore[[jj]] <- array(NA,c(draws,ncol(Y),ncol(Y)))
    for(irep in 1:draws){
      if(jj == (plag+1)){
        Astore[[jj]][irep,,] <- -t(Phistore[[jj-1]][irep,,])
      }else if(jj <= plag && jj > 1){
        Astore[[jj]][irep,,] <- t(Phistore[[jj]][irep,,]) - t(Phistore[[jj-1]][irep,,])
      }else if(jj == 1){
        PItemp <- adrop(beta_store[irep,,,drop=FALSE],drop=1)%*%adrop(alpha_store[irep,,,drop=FALSE],drop=1)
        Astore[[jj]][irep,,] <- diag(M) + PItemp - t(Phistore[[jj]][irep,,])
      }
    }
  }
  A_store <- array(NA,c(draws,ncol(Y)*(plag+1),ncol(Y)))
  for(irep in 1:draws){
    temp <- NULL
    for(jj in 1:(plag+1)){
      temp <- rbind(temp,t(Astore[[jj]][irep,,]))
    }
    A_store[irep,,] <- temp
  }
  store <- list(PHI_store=PHI_store,A_store=A_store,alpha_store=alpha_store,beta_store=beta_store,a0store=a0store,a1store=a1store,Phistore=Phistore,Astore=Astore,Exstore=Exstore,S_store=S_store,Smed_store=Smed_store,
                L_store=L_store,theta_store=theta_store,vola_store=vola_store,pars_store=pars_store,res_store=res_store,
                shrink_store=shrink_store,gamma_store=gamma_store,omega_store=omega_store,lambda2_store=lambda2_store,tau_store=tau_store)
  #------------------------------------ compute posteriors -------------------------------------------#
  PHI_post    <- apply(PHI_store,c(2,3),median)
  alpha_post  <- apply(alpha_store,c(2,3),median)
  beta_post   <- apply(beta_store,c(2,3),median)
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
  Gammapost     <- NULL
  for(jj in 1:plag){
    Gammapost    <- rbind(Gammapost,A_post[which(dims==paste("Ylag",jj,sep="")),,drop=FALSE])
  }
  post <- list(PHI_post=PHI_post,A_post=A_post,alpha_post=alpha_post,beta_post=beta_post,a0post=a0post,a1post=a1post,Gammapost=Gammapost,Expost=Expost,S_post=S_post,Sig=Sig,theta_post=theta_post,
               vola_post=vola_post,pars_post=pars_post,res_post=res_post,shrink_post=shrink_post,PIP=PIP,PIP_omega=PIP_omega,
               lambda2_post=lambda2_post,tau_post=tau_post)
  return(list(Y=Y,X=X,store=store,post=post))
}

#' @name .BVEC_linear_R
#' @importFrom MASS ginv mvrnorm
#' @importFrom Matrix bdiag
#' @importFrom methods is
#' @importFrom stats rnorm rgamma runif dnorm
#' @noRd
.BVEC_linear_R <- function(Y_in,r_in,beta_in,p_in,draws_in,burnin_in,cons_in,trend_in,sv_in,thin_in,quiet_in,prior_in,hyperparam_in,Ex_in){
  #----------------------------------------INPUTS----------------------------------------------------#
  Yraw  <- Y_in
  p     <- p_in
  r     <- r_in
  Traw  <- nrow(Yraw)
  M     <- ncol(Yraw)
  K     <- M*p
  Zraw  <- .mlag(Yraw,1)
  Ydiff <- diff(Yraw)
  Ylag  <- .mlag(Ydiff,p)
  nameslags <- NULL
  for (ii in 1:p) nameslags <- c(nameslags,rep(paste("Ylag",ii,sep=""),M))
  colnames(Ylag) <- nameslags
  colnames(Zraw) <- rep("Yt-1",M)

  texo <- FALSE; Mex <- 0; Exraw <- NULL
  if(!is.null(Ex_in)){
    Exraw <- Ex_in; Mex <- ncol(Exraw)
    texo <- TRUE
    colnames(Exraw) <- rep("Tex",Mex)
    Exraw <- Exraw[-1]
  }

  X <- cbind(Ylag,Exraw)
  Z <- Zraw[(p+2):Traw,,drop=FALSE]
  X <- X[(p+1):(Traw-1),,drop=FALSE]
  Y <- Ydiff[(p+1):(Traw-1),,drop=FALSE]
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
  kr    <- k+r
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
  if(!is.null(beta_in)){
    beta <- matrix(beta_in,M,r)
    draw_beta <- FALSE
  }else{
    beta <- matrix(0,M,r)
    diag(beta) <- 1
    draw_beta <- TRUE
  }
  D      <- Z%*%beta
  colnames(D) <- paste0("ECM",seq(1,r))
  Xtilde <- cbind(D,X)

  XtXinv <- try(solve(crossprod(Xtilde)),silent=TRUE)
  if(is(XtXinv,"try-error")) XtXinv <- ginv(crossprod(Xtilde))
  PHI_OLS  <- XtXinv%*%crossprod(Xtilde,Y)
  E_OLS    <- Y - Xtilde%*%PHI_OLS
  S_OLS    <- crossprod(E_OLS)/(bigT-k)
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  beta_draw <- beta
  PHI_draw  <- PHI_OLS
  S_draw    <- array(S_OLS, c(M,M,bigT))
  Smed_draw <- S_OLS
  Smedinv   <- solve(S_OLS)
  Em        <- Em_str <- E_OLS
  L_draw    <- diag(M)
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  # Priors on VAR coefs
  #-----------------------------
  # prior mean
  PHI_prior <- matrix(0,kr,M)
  PHI_prior[(r+1):(r+M),] <- diag(M)*prmean
  phi_prior  <-  as.vector(PHI_prior)
  # prior variance
  theta <- matrix(10,kr,M)

  # MN stuff
  accept1 <- 0
  accept2 <- 0
  accept4 <- 0
  scale1  <- .43
  scale2  <- .43
  scale4  <- .43
  sigma_sq  <- matrix(0,M,1) #vector which stores the residual variance
  for (i in 1:M){
    Ylag_i         <- .mlag(Ydiff[,i],p)
    Ylag_i         <- Ylag_i[(p+1):nrow(Ylag_i),,drop=FALSE]
    Y_i            <- Ydiff[(p+1):nrow(Ydiff),i,drop=FALSE]
    Ylag_i         <- cbind(Ylag_i,seq(1,nrow(Y_i)))
    alpha_i        <- solve(crossprod(Ylag_i))%*%crossprod(Ylag_i,Y_i)
    sigma_sq[i,1]  <- (1/(nrow(Y_i)-p-1))*t(Y_i-Ylag_i%*%alpha_i)%*%(Y_i-Ylag_i%*%alpha_i)
  }
  if(prior==1){
    theta <- .get_V(k=k,M=M,p=p,a_bar_1=shrink1,a_bar_2=shrink2,a_bar_3=shrink3,
                    a_bar_4=shrink4,sigma_sq=sigma_sq,trend=trend)
  }

  # SSVS stuff
  gamma  <-  matrix(1,kr,M)
  sigma_alpha  <-  sqrt(diag(kronecker(S_OLS,XtXinv)))
  tau0 <- matrix(NA, kr, M); tau1 <- matrix(NA, kr, M)
  ii <- 1
  for(mm in 1:M){
    for(kk in 1:kr){
      tau0[kk,mm] <- tau00*sigma_alpha[ii]
      tau1[kk,mm] <- tau11*sigma_alpha[ii]
      ii <- ii+1
    }
  }
  # NG stuff
  lambda2_PHI    <- matrix(0.01,p,1)
  PHI_tau        <- matrix(a_start,p,1)
  colnames(PHI_tau) <- colnames(lambda2_PHI) <- "endo"
  rownames(PHI_tau) <- rownames(lambda2_PHI) <- paste("lag.",seq(1,p),sep="")
  PHI_tuning     <- matrix(.43,p,1)
  PHI_accept     <- matrix(0,p,1)
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
  Sv_priors <- specify_priors(mu=sv_normal(mean=bmu, sd=Bmu), phi=sv_beta(a0,b0), sigma2=sv_gamma(shape=0.5,rate=1/(2*Bsigma)))
  eta <- list()
  #------------------------------------
  # prior for beta
  #------------------------------------
  beta_prior <- matrix(0,M,r)
  bS  <- as.vector(diag(M)[,1:r])
  bH  <- list()
  H.i <- t(cbind(matrix(0,M-r,r),diag(M-r)))
  for(jj in 1:r) bH[[jj]] <- H.i
  H   <- as.matrix(bdiag(bH))
  P_mean <- beta_draw*1
  diag(P_mean) <- 1
  P_mat <- matrix(1,M,r)
  P_mat[upper.tri(P_mat)] <- 0
  P_mat[1:r,][lower.tri(P_mat[1:r,])] <- 0

  # NG
  lambda2_B <- 10^2
  B_tau     <- a_start
  B_accept  <- 0
  B_tuning   <- .43
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
  PHI_store    <- array(NA,c(thindraws,kr,M))
  L_store      <- array(NA,c(thindraws,M,M))
  res_store    <- array(NA,c(thindraws,bigT,M))
  beta_store   <- array(NA,c(thindraws,M,r))
  P_store      <- array(NA,c(thindraws,M,r))
  # SV
  Sv_store     <- array(NA,c(thindraws,bigT,M))
  pars_store   <- array(NA,c(thindraws,4,M))
  # MN
  shrink_store <- array(NA,c(thindraws,3))
  # SSVS
  gamma_store  <- array(NA,c(thindraws,kr,M))
  omega_store  <- array(NA,c(thindraws,M,M))
  # NG
  theta_store  <- array(NA,c(thindraws,kr,M))
  lambda2_store<- array(NA,c(thindraws,p,3))
  tau_store    <- array(NA,c(thindraws,p,3))
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Sample coefficients
    D      <- Z%*%beta_draw
    colnames(D) <- paste0("ECM",seq(1,r))
    Xtilde <- cbind(D,X)

    for (mm in 1:M){
      if (mm==1){
        Y.i <- Y[,mm,drop=FALSE]*exp(-0.5*Sv_draw[,mm])
        X.i <- Xtilde*exp(-0.5*Sv_draw[,mm])

        V_post <- try(chol2inv(chol(crossprod(X.i)+diag(1/theta[,mm]))),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv(crossprod(X.i)+diag(1/theta[,mm]))
        PHI_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/theta[,mm])%*%PHI_prior[,mm])

        PHI.draw.i <- try(PHI_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
        if (is(PHI.draw.i,"try-error")) PHI.draw.i <- mvrnorm(1,PHI_post,V_post)
        PHI_draw[,mm] <- PHI.draw.i
        Em[,mm] <-  Em_str[,mm] <- Y[,mm]-Xtilde%*%PHI.draw.i
      }else{
        Y.i <- Y[,mm,drop=FALSE]*exp(-0.5*Sv_draw[,mm])
        X.i <- cbind(Xtilde,Em[,1:(mm-1)])*exp(-0.5*Sv_draw[,mm])

        V_post <- try(chol2inv(chol((crossprod(X.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))))),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv((crossprod(X.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))))
        PHI_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))%*%c(PHI_prior[,mm],l_prior[mm,1:(mm-1)]))

        PHI.draw.i <- try(PHI_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
        if (is(PHI.draw.i,"try-error")) PHI.draw.i <- mvrnorm(1,PHI_post,V_post)

        PHI_draw[,mm] <- PHI.draw.i[1:ncol(Xtilde)]
        Em[,mm] <- Y[,mm]-X%*%PHI.draw.i[1:ncol(X)]
        Em_str[,mm] <- Y[,mm]-Xtilde%*%PHI.draw.i[1:ncol(Xtilde)]-Em[,1:(mm-1),drop=FALSE]%*%PHI.draw.i[(ncol(Xtilde)+1):ncol(X.i),drop=FALSE]
        L_draw[mm,1:(mm-1)] <- PHI.draw.i[(ncol(Xtilde)+1):ncol(X.i)]
      }
    }
    rownames(PHI_draw) <- colnames(Xtilde)
    #----------------------------------------------------------------------------
    # Step 2: Sample beta
    if(draw_beta){
      alph  <- PHI_draw[1:r,,drop=FALSE]
      AA    <- PHI_draw[(r+1):nrow(PHI_draw),,drop=FALSE]
      Yhat  <- as.vector(Y-X%*%AA)-kronecker(t(alph),Z)%*%bS
      Bmean <- crossprod(H,kronecker(alph%*%Smedinv,t(Z))%*%as.vector(Yhat))
      Bvar  <- t(H)%*%kronecker(alph%*%Smedinv%*%t(alph),crossprod(Z))%*%H

      Pmat  <- diag(as.vector(P_mat[(r+1):M,]))
      Pinv  <- diag(1/diag(Pmat))
      Bpost <- solve(Bvar + Pinv)
      Bmean <- Bpost %*% (Bmean + Pinv%*%as.vector(P_mean[(r+1):M,]))
      Bdraw <- Bmean + t(chol(Bpost))%*%rnorm(M-r)

      beta_draw <- matrix(bS+H%*%Bdraw,M,r)
    }
    #----------------------------------------------------------------------------
    # Step 3: different shrinkage prior setups
    # MN
    if(prior==1){
      # NOT YET IMPLEMENTED
      stop("This prior is not implemented.")
    }
    # SSVS
    if(prior==2){
      for(mm in 1:M){
        for(kk in 1:kr){
          u_i1  <-  dnorm(PHI_draw[kk,mm],PHI_prior[kk,mm],tau0[kk,mm]) * p_i
          u_i2  <-  dnorm(PHI_draw[kk,mm],PHI_prior[kk,mm],tau1[kk,mm]) * (1-p_i)
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
        slct.i    <- which(rownames(PHI_draw)==paste("Ylag",ss,sep=""))
        if(ss==1){
          slct.i <- c(which(rownames(PHI_draw)==paste("ECM",seq(1,r),sep="")),slct.i)
          if(cons) slct.i <- c(slct.i,which(rownames(PHI_draw)=="cons"))
          if(trend) slct.i <- c(slct.i,which(rownames(PHI_draw)=="trend"))
        }
        PHI.lag     <- PHI_draw[slct.i,,drop=FALSE]
        PHI.prior   <- PHI_prior[slct.i,,drop=FALSE]
        theta.lag   <- theta[slct.i,,drop=FALSE]

        M.end <- nrow(PHI.lag)
        if (ss==1){
          lambda2_PHI[ss,1] <- rgamma(1,d_lambda+PHI_tau[ss,1]*M.end^2,e_lambda+PHI_tau[ss,1]/2*sum(theta.lag))
        }else{
          lambda2_PHI[ss,1] <- rgamma(1,d_lambda+PHI_tau[ss,1]*M.end^2,e_lambda+PHI_tau[ss,1]/2*prod(lambda2_PHI[1:(ss-1),1])*sum(theta.lag))
        }
        for (jj in 1:M){
          for (ii in 1:M){
            theta.lag[jj,ii] <- do_rgig1(lambda=PHI_tau[ss,1]-0.5,
                                         chi=(PHI.lag[jj,ii]-PHI.prior[jj,ii])^2,
                                         psi=PHI_tau[ss,1]*prod(lambda2_PHI[1:ss,1]))
          }
        }
        theta[slct.i,] <- theta.lag
        theta[theta<1e-8] <- 1e-8
        #TO BE MODIFIED
        if (sample_A){
          #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
          PHI_tau_prop <- exp(rnorm(1,0,PHI_tuning[ss,1]))*PHI_tau[ss,1]
          post_PHI_tau_prop <- .atau_post(atau=PHI_tau_prop,  thetas=as.vector(theta.lag), lambda2=prod(lambda2_PHI[1:ss,1]), k=length(theta.lag))
          post_PHI_tau_old  <- .atau_post(atau=PHI_tau[ss,1], thetas=as.vector(theta.lag), lambda2=prod(lambda2_PHI[1:ss,1]), k=length(theta.lag))
          post.diff <- post_PHI_tau_prop-post_PHI_tau_old
          post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)

          if (post.diff > log(runif(1,0,1))){
            PHI_tau[ss,1]    <- PHI_tau_prop
            PHI_accept[ss,1] <- PHI_accept[ss,1]+1
          }
          if (irep<(0.5*nburn)){
            if ((PHI_accept[ss,1]/irep)>0.3)  PHI_tuning[ss,1] <- 1.01*PHI_tuning[ss,1]
            if ((PHI_accept[ss,1]/irep)<0.15) PHI_tuning[ss,1] <- 0.99*PHI_tuning[ss,1]
          }
        }
      }
      # Normal-Gamma for beta
      if(draw_beta){
        Pmat.i <- as.vector(P_mat[(r+1):M,,drop=FALSE])
        beta.i <- as.vector(beta_draw[(r+1):M,,drop=FALSE])
        beta.prior.i <- as.vector(beta_prior[(r+1):M,,drop=FALSE])

        M.end <- length(beta.i)
        lambda2_B <- rgamma(1,d_lambda+B_tau*M.end,e_lambda+B_tau/2*sum(Pmat.i))
        for(mm in 1:M.end){
          Pmat.i[mm] <- rgig(n=1,lambda=B_tau-0.5,(beta.i[mm]-beta.prior.i[mm])^2,B_tau*lambda2_B)
        }
        Pmat.i[Pmat.i<1e-7] <- 1e-7
        if(sample_A) {
          # Sample theta through a simple RWMH step
          # (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
          B_tau_prop      <- exp(rnorm(1,0,B_tuning))*B_tau
          post_B_tau_prop <- .atau_post(atau=B_tau_prop, thetas=Pmat.i, lambda2=lambda2_B, k=length(Pmat.i))
          post_B_tau_old  <- .atau_post(atau=B_tau,      thetas=Pmat.i, lambda2=lambda2_B, k=length(Pmat.i))
          post.diff       <- post_B_tau_prop-post_B_tau_old
          post.diff       <- ifelse(is.nan(post.diff),-Inf,post.diff)
          if (post.diff > log(runif(1,0,1))){
            B_tau    <- B_tau_prop
            B_accept <- B_accept+1
          }
          # Scale MH proposal during the first 50% of the burn-in stage
          if (irep<(0.5*nburn)){
            if ((B_accept/irep)>0.3)  scale <- 1.01*scale
            if ((B_accept/irep)<0.15) scale <- 0.99*scale
          }
        }
        P_mat[(r+1):M,] <- Pmat.i
      }
    }
    #----------------------------------------------------------------------------
    # Step 4: Sample variances
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
        svl[[jj]] <- svdraw
        h_ <- exp(svdraw$latent[1,])
        para$mu      <- svdraw$para[1,"mu"]
        para$phi     <- svdraw$para[1,"phi"]
        para$sigma   <- svdraw$para[1,"sigma"]
        para$latent0 <- svdraw$latent0
        Sv_draw[,jj] <- log(h_)
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
      PHI_store[count,,] <- PHI_draw
      beta_store[count,,] <- beta_draw
      L_store[count,,] <- L_draw
      res_store[count,,] <- Y-Xtilde%*%PHI_draw
      P_store[count,,] <- P_mat
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
      lambda2_store[count,,1]  <- lambda2_PHI
      lambda2_store[count,1,3] <- lambda2_B
      tau_store[count,1,2]     <- L_tau
      tau_store[count,,1]      <- PHI_tau
      tau_store[count,1,3]     <- B_tau
    }
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(PHI_store)=list(NULL,colnames(Xtilde),colnames(PHI_OLS))
  ret <- list(Y=Y,Z=Z,X=X,PHI_store=PHI_store,L_store=L_store,beta_store=beta_store,P_store=P_store,Sv_store=Sv_store,shrink_store=shrink_store,gamma_store=gamma_store,omega_store=omega_store,theta_store=theta_store,lambda2_store=lambda2_store,tau_store=tau_store,pars_store=pars_store,res_store=res_store)
  return(ret)
}

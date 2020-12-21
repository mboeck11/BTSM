#' @name .BVAR_linear_wrapper
#' @noRd
#' @importFrom abind adrop
#' @importFrom utils capture.output
.BPVAR_linear_wrapper <- function(Yraw, prior, plag, draws, burnin, cons, trend, SV, thin, default_hyperpara, Ex, applyfun, cores){
  class(Yraw) <- "numeric"
  prior_in <- ifelse(prior=="NG",3,NA)
  if(default_hyperpara[["a_log"]]){
    default_hyperpara["a_start"] <- 1/log(ncol(Yraw))
  }
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
  bpvar<-.BPVAR_linear_R(Y_in=Yraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)
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
    shrink_store  <- bvar$shrink_store; dimnames(shrink_store) <- list(NULL,c("shrink1","shrink2"))
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

#' @name .BPVAR_linear_R
#' @importFrom stochvol svsample_fast_cpp specify_priors default_fast_sv
.BPVAR_linear_R <- function(Yraw, Wraw=NULL, p, q, nsave=5000, nburn=5000, cons=FALSE, fcst = FALSE,
                            fhorz=2, sv=TRUE, ng_mean = FALSE, shrink.type = "none", store.draws = FALSE,
                            sample_theta = FALSE, cfit=FALSE, store.post=TRUE, crit_eig=1.00,
                            Multiplier=5,thin=5) {
  #------------------------------checks----------------------------------------------#
  cN        <- names(Yraw)
  N         <- length(cN)
  varnames  <- colnames(Yraw[[1]])
  Traw      <- unlist(lapply(Yraw,nrow))
  timeraw   <- lapply(Yraw, rownames)
  if(!is.null(Wraw)){
    Wraw.l <- list()
    for(cc in 1:N){
      idx <- which(rownames(Wraw)%in%timeraw[[cc]])
      Wraw.l[[cN[cc]]] <- Wraw[idx,,drop=FALSE]
      if(q > 0) {
        temp <- mlag(Wraw.l[[cN[cc]]],q)
        Wraw.l[[cN[cc]]] <- cbind(Wraw.l[[cN[cc]]], temp)
      }
    }
    varnamesw <- colnames(Wraw.l[[cN[cc]]]);w<-ncol(Wraw.l[[cN[cc]]])
  }else{varnamesw<-NULL;w<-0}
  #------------------------------Setup----------------------------------------------#
  ntot <- nburn+nsave
  Xraw <- list()
  for(cc in 1:N){
    temp <- mlag(Yraw[[cc]],p)
    if(!is.null(Wraw)){temp<-cbind(temp,Wraw.l[[cc]])}
    if(cons){temp <- cbind(temp,1);c<-1;colnames(temp)[ncol(temp)]<-"cons"}else{c<-0}
    Xraw[[cN[cc]]] <- temp
  }
  Y <- lapply(Yraw,function(l)l[(p+1):nrow(l),,drop=FALSE])
  X <- lapply(Xraw,function(l)l[(p+1):nrow(l),,drop=FALSE])
  bigT <- unlist(lapply(Y,nrow))
  M    <- length(varnames)
  K    <- M*p
  k    <- M*p+c+w
  q    <- K*M
  h    <- q*N
  MN   <- M*N
  v    <- (M*(M-1))/2
  #--------------------------OLS estimates------------------------------------------#
  A_OLS <- array(NA, c(k, M, N))
  E_OLS <- list()
  S_OLS <- array(NA, c(M, M, N))
  S_inv <- array(NA, c(M, M, N))
  for(cc in 1:N){
    Y.c <- as.matrix(Y[[cc]])
    X.c <- as.matrix(X[[cc]])
    temp        <- try(solve(crossprod(X.c))%*%t(X.c)%*%Y.c, silent=TRUE)
    if(is(temp,"try-error")) temp <- ginv(crossprod(X.c))%*%t(X.c)%*%Y.c
    A_OLS[,,cc] <- temp
    E_OLS[[cc]] <- Y.c - X.c%*%A_OLS[,,cc]
    S_OLS[,,cc] <- crossprod(E_OLS[[cc]])/(bigT[cc]-k)
    temp <- try(solve(S_OLS[,,cc]), silent=TRUE)
    if(is(temp,"try-error")) temp <- ginv(S_OLS[,,cc])
    S_inv[,,cc] <- temp
  }
  #--------------------------Initialize Gibbs sampler--------------------------------#
  Y.big <- do.call("rbind",Y)
  X.big <- do.call("rbind",X)
  arvar <- arloop(Y.big,X.big,p)
  A_draw <- A_OLS;dimnames(A_draw)[[1]]<-colnames(X.big);dimnames(A_draw)[[2]]<-colnames(Y.big);dimnames(A_draw)[[3]]<-cN
  S_draw <- S_OLS
  Em <- Em.str <- lapply(bigT,function(l)matrix(NA,l,M))
  #----------------------------PRIORS-----------------------------------------------#
  lambda1_draw <- .1
  lambda2 <- .5
  lambda3 <- 1
  lambda4 <- 10^2
  Omegab_prior <- array(0, c(k,k,M))
  for(mm in 1:M) {
    diags <- seq(1,K+w)
    ondiag <- seq(mm,M*p,by=M)
    for(pp in 1:p){
      Omegab_prior[ondiag[pp],ondiag[pp],mm] <- (1/(pp^lambda3))^2
      for(mmm in 1:M) {
        z <- (pp-1)*M+mmm
        if(z != ondiag[pp]) Omegab_prior[diags[z],diags[z],mm] <-
            (arvar[mm,1]/arvar[mmm,1])*(lambda2/(pp^lambda3))^2
      }
    }
    if(!is.null(Wraw)){
      for(ww in 1:w) {
        z <- diags[K+ww]
        Omegab_prior[z,z,mm] <- arvar[mm,1]*(lambda4)^2
      }
    }
    if(cons) Omegab_prior[K+w+c,K+w+c,mm] <- arvar[mm,1]*(lambda4)^2
  }
  # make list
  Omegab_prior <- lapply(seq(dim(Omegab_prior)[3]), function(x) Omegab_prior[ , , x])
  Omegab_prior <- Reduce(adiag,Omegab_prior)
  Omegab_priorinv <- diag(1/diag(Omegab_prior))
  V_prior <- lambda1_draw * Omegab_prior
  Vinv_prior <- diag(1/diag(V_prior))

  # prior mean for autoregressive parameters
  alpha_draw <- apply(A_OLS, c(1,2), mean)

  # prior for lambda1_draw, Gamma here
  s0 <- 0.01 # 0.5
  v0 <- 0.01 # 0.5

  # Normal-Gamma prior stuff
  tau2_coef_draw    <- matrix(.1, k, M)
  lambda2_coef_draw <- 10^2
  theta_coef_draw   <- 0.3
  cl0 <- 0.01
  dl0 <- 0.01
  tau2_scal_draw <- array(0.1, c(M, M, N))
  lambda2_scal_draw <- rep(10^2,N)
  theta_scal_draw <- rep(.3,N)
  # MH zeugs
  scale_coef <- 0.43
  accept_coef <- 0
  scale_scal <- rep(0.43,N)
  accept_scal <- rep(0,N)

  # variances

  sigma2.scale <- array(0,c(M,1,N)) # individual per country or arvar from above?
  for(cc in 1:N) {
    for (mm in 1:M){
      temp0 <- lm(Y[[cc]][-1,mm]~Y[[cc]][-nrow(Y[[cc]]),mm])
      sigma2.scale[mm,,cc] <- summary(temp0)$sigma
    }
  }

  m0 <- 2.5 + (M - 1) / 2
  n0 <- 0.5 + (M - 1) / 2
  Q0 <- 100 * n0 / m0 * diag(as.numeric(arvar))
  # Q0 <- array(0, c(M, M, N))
  # for(cc in 1:N) {
  #   Q0[,,cc] <- 100 * n0 / m0 * diag(as.numeric(sigma2.scale[,,cc]))
  # }

  #---------------------------------Create storage matrices--------------------------#
  A_store              <- array(NA, c(nsave, k, M, N))
  alpha_store          <- array(NA, c(nsave, k, M))
  S_store              <- array(NA, c(nsave, M, M, N))
  C0_store             <- array(NA, c(nsave, M, M))
  tau2_coef_store      <- array(NA, c(nsave, k, M))
  tau2_scal_store      <- array(NA, c(nsave, M, M, N))
  lambda1_store        <- array(NA, c(nsave, 1))
  lambda2_coef_store   <- array(NA, c(nsave, 1))
  lambda2_scal_store   <- array(NA, c(nsave, N))
  if(sample_theta){
    theta_coef_store   <- array(NA, c(nsave, 1))
    theta_scal_store   <- array(NA, c(nsave, N))
  }
  if(sv) pars_store    <- array(NA, c(nsave, 3, MN))
  if(cfit) {
    # fit_store          <- lapply(bigT,function(l)array(NA, c(nsave, l, M)))
    Lik_store          <- array(NA, c(nsave, 1))
    DIC_store          <- array(NA, c(1,1))
  }
  if(fcst){
    pred_store         <- array(NA, c(nsave, M, fhorz, N))
    rmse_store         <- array(NA, c(nsave, fhorz, N))
  }
  dimnames(A_store)[[2]]<-dimnames(alpha_store)[[2]]<-colnames(X.big)
  dimnames(A_store)[[3]]<-dimnames(alpha_store)[[3]]<-varnames
  dimnames(A_store)[[4]]<-cN
  #---------------------------Gibbs loop---------------------------------------#
  ntot <- nburn + nsave * thin
  counter<-0
  irep<-1
  n_thin <- 1
  while(irep < (ntot+1)) {
    #----------------------------------------------------------------------------------------
    # Step I: Sample autoregressive parameters per country
    for(cc in 1:N) {
      Y.c       <- Y[[cc]]
      X.c       <- X[[cc]]
      S_inv.c <- S_inv[,,cc]

      coefs <- drawVARcoef(Y=Y.c, X=X.c, aprior=alpha_draw,
                           Vprior=V_prior, Sigma_inv = S_inv.c)
      A_draw[,,cc]    <- coefs$A
      Em[[cc]]        <- coefs$Em
    }
    #----------------------------------------------------------------------------------------
    #Step II: Pooling prior
    alpha_draw  <- drawPOOLcoef(mean=A_draw, var=V_prior, priorvar=tau2_coef_draw)
    # Step IIb: update heterogeneity coefficient with Gamma prior -> GIG posterior
    dev <- apply(A_draw, 3, function(a) t(c(a)-c(alpha_draw))%*%Omegab_priorinv%*%(c(a)-c(alpha_draw)))
    lambda1_draw <- rgig(1, -h/2 + v0, sum(dev), 2 * s0)
    V_prior <- lambda1_draw * Omegab_prior
    Vinv_prior <- diag(1/diag(V_prior))

    #----------------------------------------------------------------------------------------
    # Step III: Normal-Gamma prior on alpha_draw (optional)
    if(shrink.type!="none"){
      if(irep==1&shrink.type=="lagwise") {
        lambda2_coef_draw  <- matrix(0.01,p,1)
        lambda2_coef_store <- matrix(0,nsave,p)
        theta_coef_draw    <- matrix(theta_coef_draw,p,1)
        if(sample_theta) {
          theta_coef_store   <- matrix(0,nsave,p)
          scale_coef         <- rep(scale_coef,p)
          accept_coef        <- rep(accept_coef,p)
        }
      }
      ng <- NormalGammaPrior(coef=alpha_draw, tau2=tau2_coef_draw, lambda2=lambda2_coef_draw,
                             theta=theta_coef_draw, cl0=cl0, dl0=dl0, shrink.type=shrink.type,
                             sample_theta = sample_theta, scale = scale_coef, accept=accept_coef,
                             irep = irep, nburn = nburn, c=c, w=w)
      tau2_coef_draw    <- ng$tau2
      lambda2_coef_draw <- ng$lambda2
      if(sample_theta) {
        theta_coef_draw <- ng$theta
        scale_coef      <- ng$scale
        accept_coef     <- ng$accept
      }
    }

    #----------------------------------------------------------------------------------------
    # Step IV: Sample Sigma from hierarchical Wishart setup
    C0_j <- bayesm::rwishart(N * (n0 + m0 * N), # N *
                             1/N * chol2inv(chol(Q0 + apply(S_inv, c(1,2), sum))))$W # Flo W, 1/N *
    for(cc in 1:N) {
      # following code in MS_VAR
      scale0 <- crossprod(Em[[cc]])/2 + C0_j
      v_post <- bigT[[cc]] / 2 + m0
      S_draw[,,cc] <- bayesm::rwishart(N * v_post, 1/N * chol2inv(chol(scale0)))$IW # Flo IW, N *, 1/N *

      S_inv[,,cc] <- solve(S_draw[,,cc])
    }

    # Step VI: Check Stationarity
    Cm <- gen_compMat(alpha_draw,M,p)$Cm
    if(max(abs(Re(eigen(Cm)$values)))>crit_eig && irep > nburn && counter < (nburn+Multiplier*nsave)){
      irep    <- irep
      counter <- counter+1
      next
    }
    #----------------------------------------------------------------------------------------
    # Step VII: Store draws after burn-in/ Compute forecasts/ Impulse responses etc.
    if(irep > nburn & (irep %% thin == 0)) {
      # Step VIIa: compute fit
      if(cfit){
        fit<-list()
        logLik<-0
        for(cc in 1:N){
          Y.c <- Y[[cc]]; X.c <- X[[cc]]
          A.c <- A_draw[,,cc]; S.c <- S_draw[,,cc]
          # fit[[cN[cc]]] <- X.c%*%A.c
          for(tt in 1:bigT[cc]) logLik <- logLik + dmvnorm(Y.c[tt,], X.c[tt,]%*%A.c, S.c, log=TRUE)
        }
      }
      # Step VIIb: save in containers
      A_store[n_thin,,,]                          <- A_draw
      alpha_store[n_thin,,]                       <- alpha_draw
      lambda1_store[n_thin,]                      <- lambda1_draw
      S_store[n_thin,,,]                          <- S_draw
      C0_store[n_thin,,]                          <- C0_j
      tau2_coef_store[n_thin,,]                   <- tau2_coef_draw
      tau2_scal_store[n_thin,,,]                  <- tau2_scal_draw
      lambda2_coef_store[n_thin,]                 <- lambda2_coef_draw
      lambda2_scal_store[n_thin,]                 <- lambda2_scal_draw
      if(sample_theta) theta_coef_store[n_thin,]  <- theta_coef_draw
      if(sample_theta) theta_scal_store[n_thin,]  <- theta_scal_draw

      if(cfit) {
        # fit_store[[cc]][n_thin,,] <- fit[[cc]]
        Lik_store[n_thin,]    <- logLik
        if(n_thin == nsave) {
          avgLik <- 0
          A.mean <- apply(A_store, c(2,3,4), mean)
          S.mean <- apply(S_store, c(2,3,4), mean)
          for(cc in 1:N){
            Y.c <- Y[[cc]]; X.c <- X[[cc]]
            A.c <- A.mean[,,cc]; S.c <- S.mean[,,cc]
            # fit[[cN[cc]]] <- X.c%*%A.c
            for(tt in 1:bigT[cc]) avgLik <- avgLik + dmvnorm(Y.c[tt,], X.c[tt,]%*%A.c, S.c, log=TRUE)
          }
          pD <- 2 * (avgLik - mean(Lik_store))
          Dbar <- -2 * mean(Lik_store)
          DIC_store[1,1] <- pD + Dbar
        }
      }

      n_thin <- n_thin + 1
    } # END OF STEP V
    irep    <- irep+1
    counter <- counter+1
    if(irep%%50==0) print(paste0("Round: ",irep))
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,colnames(X),colnames(A_OLS))
  ret <- list(Y=Y,X=X,A_store=A_store,L_store=L_store,Sv_store=Sv_store,shrink_store=shrink_store,gamma_store=gamma_store,omega_store=omega_store,theta_store=theta_store,lambda2_store=lambda2_store,tau_store=tau_store,pars_store=pars_store,res_store=res_store)
  return(ret)
}

#----------------------------------------------------------------------------------------------------
drawVARcoef <- function(Y,X,aprior,Vprior,Sigma_inv) {
  M <- ncol(Y)
  T <- nrow(Y)
  K <- round(ncol(X))
  if(all(dim(Vprior)==c(K,M))){
    Vinvprior <- matrix(0, K*M, K*M)
    for(mm in 1:M) {
      Vinvprior[((mm-1)*K+1):(mm*K),((mm-1)*K+1):(mm*K)] <- diag(1/Vprior[,mm])
    }
  }else{
    # Vinvprior <- solve(Vprior)
    Vinvprior <- diag(1/diag(Vprior))
  }

  # container
  A   <- matrix(NA, K, M)
  Em  <- matrix(NA, T, M)
  rownames(Em) <- rownames(Y)
  colnames(Em) <- colnames(Y)

  psi_xx  <-  kronecker(Sigma_inv,crossprod(X))

  V_post  <-  try(solve(psi_xx + Vinvprior),silent=TRUE)
  if (is(V_post,"try-error")) V_post <- MASS::ginv(psi_xx + Vinvprior)

  IXY  <-   kronecker(diag(M),crossprod(X,Y))
  visig <- as.vector(Sigma_inv)
  a_post  <-  V_post%*%(IXY%*%visig + Vinvprior%*%as.vector(aprior))


  alpha  <-  try(a_post + t(chol(V_post))%*%rnorm(M*ncol(X),0,1),silent=TRUE) # Draw alpha
  if (is(alpha,"try-error")) alpha <- t(mvtnorm::rmvnorm(1, a_post, V_post))

  A <- matrix(alpha, K, M)

  Em <- Y - X %*% A

  return(list(A=A,
              Em=Em))
}

#-------------------------------------------------------------------------------------------------------------
drawPOOLcoef <- function(mean, var, priormean=NULL, priorvar=NULL, ng = FALSE){
  K <- dim(mean)[[1]]
  M <- dim(mean)[[2]]
  N <- dim(mean)[[3]]
  if(!is.null(priorvar)){
    if(all(dim(priorvar)==c(K,M))){
      priorvarinv <- matrix(0, K*M, K*M)
      for(mm in 1:M) {
        priorvarinv[((mm-1)*K+1):(mm*K),((mm-1)*K+1):(mm*K)] <- diag(1/priorvar[,mm])
      }
    }else{
      priorvarinv <- diag(1/diag(priorvar))
    }
  }else{
    priorvarinv <- diag(K*M)/10
  }
  if(is.null(priormean)) priormean <- matrix(0, K, M)


  if(ng){
    coef <- matrix(NA, K, M)
    for(mm in 1:M){
      mean.i        <- mean[,mm,]
      varinv.i      <- diag(1/diag(var[((mm-1)*K+1):(mm*K),((mm-1)*K+1):(mm*K)]))
      priormean.i   <- priormean[,mm]
      priorvarinv.i <- priorvarinv[((mm-1)*K+1):(mm*K),((mm-1)*K+1):(mm*K)]

      # posterior para
      S_post <- try(chol2inv(chol(N * varinv.i + priorvarinv.i)), silent=TRUE)
      if(is(S_post,"try-error")) S_post <- solve(N * varinv.i + priorvarinv.i)
      mu_post <- S_post %*% (varinv.i %*% apply(mean.i, 1, sum) + priorvarinv.i %*% priormean.i)

      # posterior draw
      temp <- try(mu_post + t(chol(S_post))%*%rnorm(K), silent=TRUE)
      if(is(temp,"try_error")) temp <- rmvnorm(1, mu_post, S_post)
      coef[,mm] <- temp
    }
  } else {
    mean.all <- apply(mean, c(1, 2), mean)
    mean.all <- as.vector(mean.all)
    mean.var <- var/N
    coef <- matrix(rmvnorm(1, mean.all, mean.var), K, M)
  }

  return(coef)
}

#------------------------------------------------------------------------------------------------------------

theta_post <- function(theta=theta,lambda2=lambda2,tau2=tau2,k=length(tau2),rat=1){
  logpost <- sum(dgamma(tau2,theta,(theta*lambda2/2),log=TRUE))+dexp(theta,rate=rat,log=TRUE)
  return(logpost)
}

NormalGammaPrior <- function(coef, tau2, lambda2, theta, cl0, dl0, shrink.type = "global",
                             prior=NULL, sample_theta, scale, accept,
                             irep, nburn, c=0, w=0) {
  require(GIGrvg)
  if(is.matrix(coef)) {
    K <- nrow(coef)-w-c
    M <- ncol(coef)
    p <- K/M
    if(is.null(prior)) prior <- matrix(0,K+w+c,M)
  } else {
    M <- length(coef)
    v <- (M*(M-1))/2
    if(is.null(prior)) prior <- matrix(0,M,M)
  }

  if(shrink.type == "lagwise") {
    if(nrow(lambda2)!=p) {
      print("lambda2 has wrong format, check!")
      return(NULL)
    }
  }
  if(is.null(scale)) scale <- 0.5
  if(is.list(coef)) shrink.type = "covariance"

  if(shrink.type == "global") {
    # draw lambda
    cl <- cl0 + (K+w+c)*theta
    dl <- dl0 + .5*theta*sum(tau2)
    lambda2 <- rgamma(1,cl,dl)
    # draw tau
    for(kk in 1:(K+w+c)) {
      for(mm in 1:M) {
        tau2[kk,mm] <- GIGrvg::rgig(1,
                                    lambda = theta-0.5,
                                    chi    = (coef[kk,mm]-prior[kk,mm])^2,
                                    psi    = lambda2 * theta)
      }
    }
    tau2[tau2<1e-7] <- 1e-7
    if(sample_theta) {
      # Sample theta through a simple RWMH step
      # (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
      theta_prop      <- exp(rnorm(1,0,scale))*theta
      post_theta_prop <- theta_post(theta = theta_prop,
                                    tau2  = as.vector(tau2),
                                    lambda2 = lambda2)
      post_theta_old  <- theta_post(theta = theta,
                                    tau2  = as.vector(tau2),
                                    lambda2 = lambda2)
      post.diff       <- post_theta_prop-post_theta_old
      post.diff       <- ifelse(is.nan(post.diff),-Inf,post.diff)
      if (post.diff > log(runif(1,0,1))){
        theta <- theta_prop
        accept <- accept+1
      }
      # Scale MH proposal during the first 50% of the burn-in stage
      if (irep<(0.5*nburn)){
        if ((accept/irep)>0.3)  scale <- 1.01*scale
        if ((accept/irep)<0.15) scale <- 0.99*scale
      }
    }
  }
  if(shrink.type == "lagwise") {
    for(pp in 1:p) {
      if(pp==1&c>0&w>0){
        coef.lag  <- coef[c(((pp-1)*M+1):(pp*M),(K+1):(K+w+c)),]
        prior.lag <- prior[c(((pp-1)*M+1):(pp*M),(K+1):(K+w+c)),]
        tau2.lag  <- tau2[c(((pp-1)*M+1):(pp*M),(K+1):(K+w+c)),]
      }else if(pp==1&c>0&w==0){
        coef.lag  <- coef[c(((pp-1)*M+1):(pp*M),(K+1)),]
        prior.lag <- prior[c(((pp-1)*M+1):(pp*M),(K+1)),]
        tau2.lag  <- tau2[c(((pp-1)*M+1):(pp*M),(K+1)),]
      }else if(pp==1&c==0&w>0){
        coef.lag  <- coef[c(((pp-1)*M+1):(pp*M),(K+1):(K+w)),]
        prior.lag <- prior[c(((pp-1)*M+1):(pp*M),(K+1):(K+w)),]
        tau2.lag  <- tau2[c(((pp-1)*M+1):(pp*M),(K+1):(K+w)),]
      }else{
        coef.lag  <- coef[((pp-1)*M+1):(pp*M),]
        prior.lag <- prior[((pp-1)*M+1):(pp*M),]
        tau2.lag  <- tau2[((pp-1)*M+1):(pp*M),]
      }
      # draw lambda
      if(pp==1) {
        cl <- cl0 + theta*M^2
        dl <- dl0 + .5*theta*sum(tau2.lag)
        lambda2[pp,1] <- rgamma(1,cl,dl)
      } else {
        cl <- cl0 + theta*M^2
        dl <- dl0 + .5*theta*prod(lambda2[1:(pp-1),1])*sum(tau2.lag)
        lambda2[pp,1] <- rgamma(1,cl,dl)
      }
      # draw tau
      Mt <- nrow(coef.lag)
      for(kk in 1:Mt) {
        for(mm in 1:M) {
          xi.scale <- ifelse(prod(lambda2[1:pp,1])==0,1e-8,prod(lambda2[1:pp,1]))
          tau2.lag[kk,mm] <- GIGrvg::rgig(1,
                                          lambda = theta[pp]-0.5,
                                          chi    = (coef.lag[kk,mm]-prior.lag[kk,mm])^2,
                                          psi    = xi.scale*theta[pp])
        }
      }
      tau2.lag[tau2.lag<1e-7]    <- 1e-7
      if(pp==1&c>0&w>0){
        tau2[c(((pp-1)*M+1):(pp*M),(K+1):(K+w+c)),] <- tau2.lag
      }else if(pp==1&c>0&w==0){
        tau2[c(((pp-1)*M+1):(pp*M),(K+1)),]         <- tau2.lag
      }else if(pp==1&c==0&w>0){
        tau2[c(((pp-1)*M+1):(pp*M),(K+1):(K+w)),]   <- tau2.lag
      }else{
        tau2[((pp-1)*M+1):(pp*M),]                  <- tau2.lag
      }
      if(sample_theta){
        #Sample theta through a simple RWMH step
        # (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
        theta_prop <- exp(rnorm(1,0,scale[pp]))*theta[pp]
        post_theta_prop <- theta_post(theta = theta_prop,
                                      tau2  = as.vector(tau2.lag),
                                      lambda2 = prod(lambda2[1:(pp),1]))
        post_theta_old  <- theta_post(theta = theta[pp],
                                      tau2  = as.vector(tau2.lag),
                                      lambda2 = prod(lambda2[1:(pp),1]))
        post.diff <- post_theta_prop-post_theta_old
        post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)

        if (post.diff > log(runif(1,0,1))){
          theta[pp] <- theta_prop
          accept[pp] <- accept[pp]+1
        }
        # Scale MH proposal during the first 50% of the burn-in stage
        if (irep<(0.5*nburn)){
          if ((accept[pp]/irep)>0.3)  scale[pp] <- 1.01*scale[pp]
          if ((accept[pp]/irep)<0.15) scale[pp] <- 0.99*scale[pp]
        }
      }
    }
  }
  if(shrink.type == "covariance") {
    # draw lambda
    cl <- cl0 + theta*v
    cl <- cl0 + theta*v
    dl <- dl0 + .5*theta*sum(tau2[lower.tri(tau2)])

    lambda2 <- rgamma(1,cl,dl)
    # Step Vb: sample the prior scaling factors for covariances from GIG
    for(mm in 2:M) {
      et_0 <- coef[[mm]]
      for(jj in 1:length(et_0)) {
        tau2[mm,jj] <- rgig(1,
                            lambda = theta-0.5,
                            chi = et_0[jj]^2,
                            psi = theta*lambda2)
      }
    }
    tau2[tau2<1e-7 && tau2!=0] <- 1e-7
    if(sample_theta) {
      # Sample theta through a simple RWMH step
      # (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
      theta_prop      <- exp(rnorm(1,0,scale))*theta
      post_theta_prop <- theta_post(theta = theta_prop,
                                    tau2  = as.vector(tau2),
                                    lambda2 = lambda2)
      post_theta_old  <- theta_post(theta = theta,
                                    tau2  = as.vector(tau2),
                                    lambda2 = lambda2)
      post.diff       <- post_theta_prop-post_theta_old
      post.diff       <- ifelse(is.nan(post.diff),-Inf,post.diff)
      if (post.diff > log(runif(1,0,1))){
        theta <- theta_prop
        accept <- accept+1
      }
      # Scale MH proposal during the first 50% of the burn-in stage
      if (irep<(0.5*nburn)){
        if ((accept/irep)>0.3)  scale <- 1.01*scale
        if ((accept/irep)<0.15) scale <- 0.99*scale
      }
    }
  }

  if(sample_theta) {
    return(list(tau2    = tau2,
                lambda2 = lambda2,
                theta   = theta,
                scale  = scale,
                accept  = accept))
  } else {
    return(list(tau2    = tau2,
                lambda2 = lambda2,
                theta   = theta))
  }
}

#---------------------------------------------------------------------------------------------------------
panelarloop <- function(M,N,p,T,Y,X,ident.country,ident.country.X,W=NULL) {
  arvar <- matrix(0, M, 1)

  if(!is.null(W)) w <- ncol(W)

  for(mm in 1:M) {
    Ystack <- matrix(0, T*N, 1)
    Xstack <- matrix(0, T*N, M*p)
    if(!is.null(W)) Wstack <- matrix(0, T*N, w)

    for(cc in 1:N) {
      sl.country <- ident.country[cc,mm]
      Ystack[((cc-1)*T+1):(T*cc),] <- Y[,sl.country]
      for(k in 1:M) {
        sl.country.X <- ident.country.X[cc,seq(k,M*p,by=M)]
        Xstack[((cc-1)*T+1):(T*cc),seq(k,M*p,by=M)] <- X[,sl.country.X]
      }
      if(!is.null(W)) Wstack[((cc-1)*T+1):(T*cc),] <- W
    }
    if(!is.null(W)) Xstack <- cbind(Xstack,Wstack)
    B <- try(solve(crossprod(Xstack))%*%t(Xstack)%*%Ystack,silent=TRUE)
    if(is(B,"try-error")) B <- ginv(crossprod(Xstack))%*%t(Xstack)%*%Ystack
    eps <- Ystack - Xstack%*%B
    K <- ncol(Xstack)

    arvar[mm,] <- crossprod(eps)*(1/(N*T-(K+1)))
  }
  return(arvar)
}

arloop <- function(Y,X,p,W=NULL) {
  M     <- ncol(Y)
  bigT  <- nrow(Y)
  arvar <- matrix(0, M, 1)
  k     <- ncol(X)

  for(mm in 1:M) {
    Y.i <- Y[,mm]
    X.i <- X
    if(!is.null(W)) X.i <- cbind(X.i,W)

    B <- try(solve(crossprod(X.i))%*%t(X.i)%*%Y.i,silent=TRUE)
    if(is(B,"try-error")) B <- ginv(crossprod(X.i))%*%t(X.i)%*%Y.i
    eps <- Y.i - X.i%*%B

    arvar[mm,] <- crossprod(eps)/(bigT-(k))
  }
  return(arvar)
}

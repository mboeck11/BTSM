#' @name .bpvar_linear_wrapper
#' @noRd
#' @importFrom abind adrop
#' @importFrom utils capture.output
.BPVAR_linear_wrapper <- function(Yraw, prior, plag, draws, burnin, cons, trend, SV, thin, default_hyperpara, Ex, applyfun, cores){
  prior_in <- ifelse(prior=="NG",3,NA)
  if(default_hyperpara[["a_log"]]) default_hyperpara["a_start"] <- 1/log(ncol(Yraw))
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
  Y <- bpvar$Y; colnames(Y) <- colnames(Yraw); X <- bpvar$X
  Y.big
  M <- ncol(Y); bigT <- nrow(Y); K <- ncol(X)
  if(!is.null(Ex)) Mex <- ncol(Ex)
  xnames <- paste(rep("Ylag",M),rep(seq(1,plag),each=M),sep="")
  if(!is.null(Ex)) xnames <- c(xnames,paste(rep("Tex",Mex)))
  if(cons)  xnames <- c(xnames,"cons")
  if(trend) xnames <- c(xnames,"trend")
  colnames(X) <- xnames
  #-----------------------------------------get containers ------------------------------------------#
  A_store <- bpvar$A_store; dimnames(A_store)[[2]] <- colnames(X); dimnames(A_store)[[3]] <- colnames(Y)
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
    L_store <- bpvar$L_store
    for(irep in 1:(draws/thin)){
      for(tt in 1:bigT){
        if(M>1){
          S_store[irep,tt,,] <- L_store[irep,,]%*%diag(exp(bpvar$Sv_store[irep,tt,]))%*%t(L_store[irep,,])
        }else{
          S_store[irep,tt,,] <- L_store[irep,,]%*%exp(bpvar$Sv_store[irep,tt,])%*%t(L_store[irep,,])
        }
      }
    }
    Smed_store <- apply(S_store,c(1,3,4),median)
    if(SV){
      vola_store  <- bpvar$Sv_store; dimnames(vola_store) <- list(NULL,NULL,colnames(Y))
      pars_store  <- bpvar$pars_store
      vola_post   <- apply(vola_store,c(2,3),median)
      pars_post   <- apply(pars_store,c(2,3),median)
    }else{
      vola_store  <- bpvar$Sv_store; pars_store <- NULL;
      vola_post   <- apply(vola_store,c(2,3),median); pars_post <- NULL
    }
  }else if(prior=="NC"){
    Smed_store  <- bpvar$S_store
    for(irep in 1:(draws/thin)){
      for(tt in 1:bigT){
        S_store[irep,tt,,] <- bpvar$S_store[irep,,]
      }
    }
    L_store     <- NULL
    theta_store <- NULL
    vola_store  <- NULL
    pars_store  <- NULL
    vola_post   <- NULL
    pars_post   <- NULL
  }
  theta_store   <- bpvar$theta_store; dimnames(theta_store)[[2]] <- colnames(X); dimnames(theta_store)[[3]] <- colnames(Y)
  res_store     <- bpvar$res_store; dimnames(res_store) <- list(NULL,NULL,colnames(Y))
  # MN
  if(prior=="MN"){
    shrink_store  <- bpvar$shrink_store; dimnames(shrink_store) <- list(NULL,c("shrink1","shrink2"))
    shrink_post   <- apply(shrink_store,2,median)
  }else{
    shrink_store  <- shrink_post <- NULL
  }
  # SSVS
  if(prior=="SSVS"){
    gamma_store <- bpvar$gamma_store; dimnames(gamma_store) <- list(NULL,colnames(X),colnames(Y))
    omega_store <- bpvar$omega_store; dimnames(omega_store) <- list(NULL,colnames(Y),colnames(Y))
    PIP         <- apply(gamma_store,c(2,3),mean)
    PIP_omega   <- apply(omega_store,c(2,3),mean)
  }else{
    gamma_store <- omega_store <- PIP <- PIP_omega <- NULL
  }
  # NG
  if(prior=="NG"){
    lambda2_store <- bpvar$lambda2_store
    tau_store     <- bpvar$tau_store
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
#' @noRd
#' @importFrom magic adiag
#' @importFrom mvtnorm rmvnorm
#' @importFrom stochvol svsample_fast_cpp specify_priors default_fast_sv
#' @importFrom bayesm rwishart
.BPVAR_linear_R <- function(Y_in,p_in,draws_in,burnin_in,cons_in,trend_in,sv_in,thin_in,prior_in,hyperparam_in,Ex_in){
  #------------------------------checks----------------------------------------------#
  Yraw    <- Y_in
  p       <- p_in
  cN      <- names(Yraw)
  N       <- length(cN)
  Traw    <- unlist(lapply(Yraw,nrow))
  timeraw <- lapply(Yraw, rownames)
  names <- colnames(Yraw[[1]])
  if(is.null(names)) names <- rep("Y",M)
  M <- length(names)
  K <- M*p
  nameslags <- NULL
  for(ii in 1:p) nameslags <- c(nameslags,paste0(names,".lag",ii))
  Ylag <- lapply(Yraw,function(y){
    y<-.mlag(y,p)
    colnames(y)<-nameslags
    return(y)
  })

  texo <- FALSE; Mex <- 0; Exraw <- NULL; enames <- NULL
  if(!is.null(Ex_in)){
    Exraw <- Ex_in
    enames <- paste0("Tex.",colnames(Exraw[[1]]))
    if(is.null(enames)) enames <- rep("Tex",Mex)
    Exraw <- lapply(Exraw,function(e){
      colnames(e)<-enames
      return(e)
    })
    Mex <- length(enames)
    texo <- TRUE
  }else{
    Exraw <- vector(mode="list", length=N)
  }

  X  <- lapply(1:N,function(cc) cbind(Ylag[[cc]],Exraw[[cc]]))
  names(X) <- cN
  X  <- lapply(1:N,function(cc) X[[cc]][(p+1):nrow(X[[cc]]),,drop=FALSE])
  Y  <- lapply(1:N,function(cc) Yraw[[cc]][(p+1):Traw[cc],,drop=FALSE])
  bigT <- sapply(Y,nrow)

  cons  <- cons_in
  if(cons){
    X <- lapply(X,function(x){
      x<-cbind(x,1)
      colnames(x)[ncol(x)] <- "cons"
      return(x)
    })
  }
  trend <- trend_in
  if(trend){
    X <- lapply(X,function(x){
      x<-cbind(x,seq(1,nrow(x)))
      colnames(x)[ncol(x)] <- "trend"
      return(x)
    })
  }
  k  <- K+Mex+cons+trend
  q  <- k*M
  h  <- q*N
  MN <- M*N
  v  <- (M*(M-1))/2
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
  s0        <- hyperpara$s0
  v0        <- hyperpara$v0
  # prior == 3: NG
  d_lambda  <- hyperpara$d_lambda
  e_lambda  <- hyperpara$e_lambda
  a_start   <- hyperpara$a_start
  sample_A  <- hyperpara$sample_A
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
  sigma_sq  <- matrix(0,M,1) #vector which stores the residual variance
  for(mm in 1:M){
    Y_m            <- Y.big[,mm,drop=FALSE]
    X_m            <- X.big # check!!!!
    alpha_i        <- solve(crossprod(X_m))%*%crossprod(X_m,Y_m)
    sigma_sq[mm,1]  <- (1/(nrow(Y_m)-ncol(X_m)))*crossprod(Y_m-X_m%*%alpha_i)
  }
  A_draw <- A_OLS;dimnames(A_draw)[[1]]<-colnames(X.big);dimnames(A_draw)[[2]]<-colnames(Y.big);dimnames(A_draw)[[3]]<-cN
  alpha_draw <- apply(A_OLS, c(1,2), mean);dimnames(alpha_draw) <- list(colnames(X.big),colnames(Y.big))
  S_draw <- S_OLS
  Em <- Em.str <- lapply(bigT,function(l)matrix(NA,l,M))
  #----------------------------PRIORS-----------------------------------------------#
  Omegab_prior <- array(0, c(k,k,M))
  for(mm in 1:M) {
    diags <- seq(1,K+Mex)
    ondiag <- seq(mm,M*p,by=M)
    for(pp in 1:p){
      Omegab_prior[ondiag[pp],ondiag[pp],mm] <- (1/(pp^shrink3))^2
      for(mmm in 1:M) {
        z <- (pp-1)*M+mmm
        if(z != ondiag[pp]) Omegab_prior[diags[z],diags[z],mm] <-
            (sigma_sq[mm,1]/sigma_sq[mmm,1])*(shrink2/(pp^shrink3))^2
      }
    }
    if(texo){
      for(ww in 1:Mex) {
        z <- diags[K+ww]
        Omegab_prior[z,z,mm] <- sigma_sq[mm,1]*(shrink4)^2
      }
    }
    if(cons) Omegab_prior[K+Mex+1,K+Mex+1,mm] <- sigma_sq[mm,1]*shrink4^2
    if(trend) Omegab_prior[K+Mex+2,K+Mex+2,mm] <- sigma_sq[mm,1]*shrink4^2
  }
  # make list
  Omegab_prior <- lapply(seq(dim(Omegab_prior)[3]), function(x) Omegab_prior[,,x])
  Omegab_prior <- Reduce(adiag,Omegab_prior)
  Omegab_priorinv <- diag(1/diag(Omegab_prior))
  V_prior <- shrink1 * Omegab_prior
  Vinv_prior <- diag(1/diag(V_prior))

  # prior mean for common mean
  alpha_prior <- matrix(0,k,M)
  # prior variance for common mean
  theta <- matrix(.1,k,M)

  # NG A_draw
  lambda2_A    <- matrix(0.01,p,1)
  A_tau        <- matrix(a_start,p,1)
  colnames(A_tau) <- colnames(lambda2_A) <- "endo"
  rownames(A_tau) <- rownames(lambda2_A) <- paste("lag.",seq(1,p),sep="")
  A_tuning     <- matrix(.43,p,1)
  A_accept     <- matrix(0,p,1)

  #------------------------------------
  # non-SV quantities
  #------------------------------------
  # # variances
  # sigma2.scale <- array(0,c(M,1,N)) # individual per country or arvar from above?
  # for(cc in 1:N) {
  #   for (mm in 1:M){
  #     temp0 <- lm(Y[[cc]][-1,mm]~Y[[cc]][-nrow(Y[[cc]]),mm])
  #     sigma2.scale[mm,,cc] <- summary(temp0)$sigma
  #   }
  # }
  m0 <- 2.5 + (M - 1) / 2
  n0 <- 0.5 + (M - 1) / 2
  Q0 <- 100 * n0 / m0 * diag(as.numeric(sigma_sq))
  # Q0 <- array(0, c(M, M, N))
  # for(cc in 1:N) {
  #   Q0[,,cc] <- 100 * n0 / m0 * diag(as.numeric(sigma2.scale[,,cc]))
  # }

  #------------------------------------
  # SV quantities
  #------------------------------------
  Sv_draw <- lapply(bigT, function(tt) matrix(-3,tt,M)); names(Sv_draw) <- cN
  svdraw  <- lapply(bigT, function(tt) list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,tt))); names(svdraw) <- cN
  pars_var <- array(c(-3,.9,.2,-3), c(4,M,N),
                     dimnames=list(c("mu","phi","sigma","latent0"),names,cN))
  hv <- svdraw[[1]]$latent
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
  A_store       <- array(NA, c(nsave, k, M, N))
  alpha_store   <- array(NA, c(nsave, k, M))
  S_store       <- array(NA, c(nsave, M, M, N))
  C0_store      <- array(NA, c(nsave, M, M))
  theta_store   <- array(NA, c(nsave, k, M))
  shrink_store  <- array(NA, c(nsave, 1))
  lambda2_store <- array(NA, c(nsave, p))
  tau_store     <- array(NA, c(nsave, p))
  pars_store    <- array(NA, c(nsave, 3, M, N))
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------------------
    # Step I: Sample autoregressive parameters per country
    for(cc in 1:N) {
      Y.c     <- Y[[cc]]
      X.c     <- X[[cc]]
      S_inv.c <- S_inv[,,cc]
      Vinvprior <- diag(1/diag(V_prior))

      psi_xx  <-  kronecker(S_inv.c,crossprod(X.c))
      V_post  <-  try(solve(psi_xx + Vinvprior),silent=TRUE)
      if(is(V_post,"try-error")) V_post <- MASS::ginv(psi_xx + Vinvprior)

      IXY  <-   kronecker(diag(M),crossprod(X.c,Y.c))
      visig <- as.vector(S_inv.c)
      a_post  <-  V_post%*%(IXY%*%visig + Vinvprior%*%as.vector(alpha_draw))
      a_draw  <-  try(a_post + t(chol(V_post))%*%rnorm(M*k,0,1),silent=TRUE) # Draw alpha
      if (is(a_draw,"try-error")) a_draw <- t(rmvnorm(1, a_post, V_post))

      A_draw[,,cc]    <- matrix(a_draw, k, M)
      Em[[cc]]        <- Y.c - X.c %*% matrix(a_draw, k, M)
    }
    #----------------------------------------------------------------------------------------
    #Step II: Pooling prior
    ng<-TRUE
    if(ng){
      for(mm in 1:M){
        mean.i        <- A_draw[,mm,]
        varinv.i      <- diag(1/diag(V_prior[((mm-1)*k+1):(mm*k),((mm-1)*k+1):(mm*k)]))
        priormean.i   <- alpha_prior[,mm]
        priorvarinv.i <- diag(1/c(theta[,mm]))

        # posterior para
        S_post <- try(chol2inv(chol(N * varinv.i + priorvarinv.i)), silent=TRUE)
        if(is(S_post,"try-error")) S_post <- solve(N * varinv.i + priorvarinv.i)
        mu_post <- S_post %*% (varinv.i %*% apply(mean.i, 1, sum) + priorvarinv.i %*% priormean.i)

        # posterior draw
        temp <- try(mu_post + t(chol(S_post))%*%rnorm(k), silent=TRUE)
        if(is(temp,"try_error")) temp <- rmvnorm(1, mu_post, S_post)
        alpha_draw[,mm] <- temp
      }
    } else {
      mean.all <- apply(A_draw, c(1, 2), mean)
      mean.all <- as.vector(mean.all)
      mean.var <- as.vector(tau2_coef_draw)/N
      alpha_draw <- matrix(rmvnorm(1, mean.all, diag(mean.var)), k, M)
    }
    #----------------------------------------------------------------------------------------
    # Step III: update heterogeneity coefficient with Gamma prior -> GIG posterior
    dev <- apply(A_draw, 3, function(a) t(c(a)-c(alpha_draw))%*%Omegab_priorinv%*%(c(a)-c(alpha_draw)))
    shrink1 <- rgig(1, -h/2 + v0, sum(dev), 2 * s0)
    V_prior <- shrink1 * Omegab_prior
    Vinv_prior <- diag(1/diag(V_prior))
    #----------------------------------------------------------------------------------------
    # Step III: Normal-Gamma prior on alpha_draw
    for (ss in 1:p){
      slct.i <- which(grepl(paste0(".lag",ss),rownames(alpha_draw)))
      if(ss==1) slct.i <- c(slct.i,which(grepl("Tex",rownames(alpha_draw))))
      alpha.lag   <- alpha_draw[slct.i,,drop=FALSE]
      alpha.prior <- alpha_prior[slct.i,,drop=FALSE]
      theta.lag   <- theta[slct.i,,drop=FALSE]

      M.end <- nrow(alpha.lag)
      if (ss==1){
        lambda2_A[ss,1] <- rgamma(1,d_lambda+A_tau[ss,1]*M.end^2,e_lambda+A_tau[ss,1]/2*sum(theta.lag))
      }else{
        lambda2_A[ss,1] <- rgamma(1,d_lambda+A_tau[ss,1]*M.end^2,e_lambda+A_tau[ss,1]/2*prod(lambda2_A[1:(ss-1),1])*sum(theta.lag))
      }
      for(jj in 1:M.end){
        for (ii in 1:M){
          theta.lag[jj,ii] <- do_rgig1(lambda=A_tau[ss,1]-0.5,
                                       chi=(alpha.lag[jj,ii]-alpha.prior[jj,ii])^2,
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
    #----------------------------------------------------------------------------------------
    # Step IV: Sample Sigma from hierarchical Wishart setup
    if(sv){
      stop("SV currently not implemented.")
    }else{
      C0_j <- rwishart(N * (n0 + m0 * N), # N *
                       1/N * chol2inv(chol(Q0 + apply(S_inv, c(1,2), sum))))$W # Flo W, 1/N *
      for(cc in 1:N) {
        # following code in MS_VAR
        scale0 <- crossprod(Em[[cc]])/2 + C0_j
        v_post <- bigT[[cc]] / 2 + m0
        S_draw[,,cc] <- rwishart(N * v_post, 1/N * chol2inv(chol(scale0)))$IW # Flo IW, N *, 1/N *

        S_inv[,,cc] <- solve(S_draw[,,cc])
      }
    }
    #----------------------------------------------------------------------------------------
    # Step VII: Store draws after burn-in/ Compute forecasts/ Impulse responses etc.
    if(irep %in% thin.draws){
      count <- count+1
      A_store[count,,,]  <- A_draw
      alpha_store[count,,] <- alpha_draw
      S_store[count,,,] <- S_draw
      C0_store[count,,] <- C0_j
      shrink_store[count,] <- shrink1
      theta_store[count,,] <- theta
      lambda2_store[count,]<- lambda2_A
      tau_store[count,]  <- A_tau
    } # END OF STEP V
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,colnames(X.big),colnames(Y.big),cN)
  dimnames(alpha_store)=list(NULL,colnames(X.big),colnames(Y.big))
  ret <- list(Y=Y,X=X,Y.big=Y.big,X.big=X.big,A_store=A_store,alpha_store=alpha_store,S_store=S_store,C0_store=C0_store,shrink_store=shrink_store,theta_store=theta_store,lambda2_store=lambda2_store,tau_store=tau_store,pars_store=pars_store)
  return(ret)
}


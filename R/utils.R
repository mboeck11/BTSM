#' @name .mlag
#' @noRd
.mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)] <- X[(p+1-ii):(Traw-ii),(1:N)]
  }
  colnames(Xlag) <- paste0(colnames(X),".lag",rep(seq(p),each=N))
  return(Xlag)
}

#' @name .gen_compMat
#' @export
#' @noRd
.gen_compMat <- function(A, M, p){
  Jm          <- matrix(0, M*p, M)
  Jm[1:M,1:M] <- diag(M)

  A   <- A[1:(M*p),,drop=FALSE]
  Cm  <- matrix(0, M*p, M*p)
  if(p==1) Cm <- t(A) else {
    for(j in 1:(p-1)){
      Cm[(j*M+1):(M*(j+1)),(M*(j-1)+1):(j*M)] <- diag(M)
    }
  }
  bbtemp <- A[1:(M*p),]
  splace <- 0
  for(ii in 1:p){
    for(iii in 1:M) {
      Cm[iii,((ii-1)*M+1):(ii*M)] <- t(bbtemp[(splace+1):(splace+M),iii])
    }
    splace <- splace+M
  }
  return(list(Cm=Cm,
              Jm=Jm))
}

#' @name .divisors
#' @noRd
.divisors <- function (n,div) {
  div <- round(div)
  for(dd in div:1){
    if(n%%div==0) break else div<-div-1
  }
  return(div)
}

#' @name .construct.arglist
#' @noRd
.construct.arglist = function (funobj, envir = NULL){
  namedlist = formals(funobj)
  argnames = names(namedlist)
  if (!is.environment(envir))
    envir = sys.frame(-1)
  for (argn in 1:length(namedlist)) {
    testval = as.logical(try(exists(argnames[argn], envir = envir),
                             silent = TRUE))
    if (is.na(testval))
      testval = FALSE
    if (testval) {
      testout = try(get(argnames[argn], envir = envir),silent = TRUE)
      if (is.null(testout)) {
        namedlist[[argn]] = "list(NULL)blabla"
      } else {
        namedlist[[argn]] = testout
      }
    }
  }
  namedlist = lapply(namedlist,function(x) if (any(x=="list(NULL)blabla")) NULL else x)
  lapply(namedlist, function(l) if(any(l=="list(NULL)blabla")){NULL}else{l})
  return(namedlist)
}

#' @name .theta_post
#' @noRd
#' @importFrom stats dgamma dexp
.theta_post <- function(theta=theta,lambda2=lambda2,tau2=tau2,k=length(tau2),rat=1){
  logpost <- sum(dgamma(tau2,theta,(theta*lambda2/2),log=TRUE))+dexp(theta,rate=rat,log=TRUE)
  return(logpost)
}

#' @name .atau_post
#' @importFrom stats dgamma dexp
#' @noRd
.atau_post <- function(atau,lambda2,thetas,k,rat=1){
  logpost <- sum(dgamma(thetas,atau,(atau*lambda2/2),log=TRUE))+dexp(atau,rate=rat,log=TRUE)
  return(logpost)
}

#' @name .bernoulli
#' @importFrom stats runif
#' @noRd
.bernoulli <- function(p){
  u <- runif(1)
  if (u<p){
    x=0
  }else{
    x=1
  }
  return(x)
}

#' @name .get_V
#' @noRd
.get_V <- function(k=k,M=M,p=p,a_bar_1,a_bar_2,a_bar_3,a_bar_4,sigma_sq,cons=FALSE,trend=FALSE){
  V_i <- matrix(0,k,M)
  # endogenous part
  for(i in 1:M){
    for(pp in 1:p){
      for(j in 1:M){
        if(i==j){
          #V_i[j+M*(pp-1),i] <- a_bar_1/(pp^2) ######
          V_i[j+M*(pp-1),i] <- (a_bar_1/pp)^2
        }else{
          #V_i[j+M*(pp-1),i] <- (a_bar_2 * sigma_sq[i])/(pp^2*sigma_sq[j]) #####
          V_i[j+M*(pp-1),i] <- (a_bar_2/pp)^2 * (sigma_sq[i]/sigma_sq[j])
        }
      }
    }
  }
  # deterministics
  if(cons || trend){
    for(i in 1:M){
      V_i[(k-cons-trend+1),i] <- a_bar_3 * sigma_sq[i]
    }
  }
  return(V_i)
}

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
    invisible(capture.output(bvar<-try(BVAR_linear(Y_in=Yraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=cons,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Ex)),type="message"))
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
#' @importFrom stochvol svsample2
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
  L_prior <- matrix(kappa1,M,M)
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
  pars_var <- matrix(0,3,M)

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
  pars_store   <- array(NA,c(thindraws,3,M))
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
      pars_var <- matrix(0,3,M)
      for (jj in 1:M){
        svdraw <- svsample2(Em_str[,jj],startpara=svl[[jj]]$para,startlatent=Sv_draw[,jj],
                            priormu=c(bmu,sqrt(Bmu)),priorphi=c(a0, b0), priorsigma=Bsigma)
        svl[[jj]] <- svdraw
        h_ <- exp(svdraw$latent)
        pars_var[,jj] <- svl[[jj]]$para
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
#' @importFrom stochvol svsample2
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
#' @importFrom stochvol svsample2
#' @importFrom MASS ginv mvrnorm
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
  pars_var <- matrix(0,3,M)

  hv <- svdraw$latent
  para <- list(mu=-3,phi=.9,sigma=.2)

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
  pars_store   <- array(NA,c(thindraws,3,M))
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
            if ((accept/irep)>0.3)  scale <- 1.01*scale
            if ((accept/irep)<0.15) scale <- 0.99*scale
          }
        }
        P_mat[(r+1):M,] <- Pmat.i
      }
    }
    #----------------------------------------------------------------------------
    # Step 4: Sample variances
    if (sv){
      pars_var <- matrix(0,3,M)
      for (jj in 1:M){
        svdraw <- svsample2(Em_str[,jj],startpara=svl[[jj]]$para,startlatent=Sv_draw[,jj],
                            priormu=c(bmu,sqrt(Bmu)),priorphi=c(a0, b0), priorsigma=Bsigma)
        svl[[jj]] <- svdraw
        h_ <- exp(svdraw$latent)
        pars_var[,jj] <- svl[[jj]]$para
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

#' @name .irf.sign.zero
#' @noRd
#' @importFrom MASS Null
.irf.sign.zero <- function(xdat,plag,n.ahead,Amat,Smat,shock,sign.constr,MaxTries,...){
  bigT     <- nrow(xdat)
  bigK     <- ncol(xdat)
  varNames <- colnames(xdat)

  P0G <- try(t(chol(Smat[,,drop=FALSE])),silent=TRUE)
  if(is(P0G,"try-error")) suppressWarnings(P0G <- t(chol(Smat,pivot=TRUE)))
  colnames(P0G) <- rownames(P0G) <- varNames

  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+n.ahead+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  temp <- .gen_compMat(Amat, bigK, plag)
  Cm   <- temp$Cm
  Jm   <- temp$Jm
  Cmat <- diag(bigK*plag)
  for (ihor in (plag+2):(plag+n.ahead+1)){
    Cmat  <- Cmat%*%Cm
    PHIx[,,ihor]  <- t(Jm)%*%Cmat%*%Jm
  }
  PHI  <-  PHIx[,,(plag+1):(plag+n.ahead+1)]
  #-----------------------------------------------------------------------------
  sign.horizon   <- unique(unlist(lapply(sign.constr, function(l) l$rest.horz)))-1 # zero impact is coded as 1
  sign.horizon   <- sort(sign.horizon, decreasing=FALSE)
  sign.shockvars <- unlist(lapply(sign.constr, function(l) l$shock))
  H.restr        <- length(sign.horizon)
  N.restr        <- bigK*H.restr
  S.cube         <- array(NA, c(N.restr, N.restr, bigK)) # sign restrictions
  Z.cube         <- array(NA, c(N.restr, N.restr, bigK)) # zero restrictions
  dimnames(S.cube)[[1]] <- dimnames(Z.cube)[[1]] <-
    dimnames(S.cube)[[2]] <- dimnames(Z.cube)[[2]] <- paste(rep(varNames,H.restr),".",
                                                            rep(sign.horizon,each=bigK),sep="")
  dimnames(S.cube)[[3]] <- dimnames(Z.cube)[[3]] <- varNames
  for(vv in 1:length(varNames)){
    S.temp <- matrix(0, N.restr, N.restr)
    Z.temp <- matrix(0, N.restr, N.restr)
    if(varNames[vv]%in%sign.shockvars){
      sign.shock <- sign.constr[[which(sign.shockvars == varNames[vv])]]$shock
      sign.restr <- sign.constr[[which(sign.shockvars == varNames[vv])]]$restrictions
      sign.signs <- sign.constr[[which(sign.shockvars == varNames[vv])]]$sign
      sign.horiz <- sign.constr[[which(sign.shockvars == varNames[vv])]]$rest.horz-1

      sign.restr <- c(sign.shock,sign.restr) ## append positive shock on shock variable

      s.point <- which(sign.signs=="<"|sign.signs==">")
      z.point <- seq(1,length(sign.signs))[-s.point]

      if(length(s.point)>0){
        for(ss in 1:length(s.point)){
          grp <- which(sign.horiz[s.point[ss]] == sign.horizon)
          col <- seq(which(sign.restr[s.point[ss]]==varNames),bigK*grp,by=bigK)
          for(ii in 1:length(col)){
            S.temp[col[ii],col[ii]] <- ifelse(sign.signs[s.point[ss]]=="<",-1,1)
          }
        }
      }
      if(length(z.point)>0){
        for(zz in 1:length(z.point)){
          if(sign.signs[z.point[zz]]=="0"){
            grp <- which(sign.horiz[z.point[zz]] == sign.horizon)
            row <- (grp-1)*bigK+which(sign.restr[z.point[zz]]==varNames)
            Z.temp[row,row] <- 1
          }else{ # take row from above
            grp <- which(sign.horiz[z.point[zz]] == sign.horizon)
            col <- (grp-1)*bigK+which(sign.restr[z.point[zz]]==varNames)
            Z.temp[row,col] <- as.numeric(sign.signs[z.point[zz]])
          }
        }
      }
    }
    S.cube[,,vv] <- S.temp
    Z.cube[,,vv] <- Z.temp
  }

  no.zero.restr <- ifelse(base::sum(abs(Z.cube))>0,FALSE,TRUE)
  shock.order   <- rep(NA, bigK)
  search.Znum   <- apply(Z.cube, 3, function(x) base::sum(abs(x)))
  search.Snum   <- apply(S.cube, 3, function(x) base::sum(abs(x)))
  for(mm in 1:bigK){
    if(!no.zero.restr){
      max.Z <- which(search.Znum==max(search.Znum))
      if(length(max.Z)==1){
        shock.order[mm] <- max.Z
      }else{
        shock.order[mm] <- sample(max.Z,1)
      }
    } else {
      shock.order[mm] <- mm
    }
    search.Znum[shock.order[mm]] <- -1
    search.Snum[shock.order[mm]] <- -1
  }
  shock.order <- varNames[shock.order]

  irf.restr         <- matrix(NA, N.restr, bigK)
  for(hh in 1:H.restr){
    # ARRW: Definition 1
    if(sign.horizon[hh]!=Inf) irf.hh<-PHI[,,sign.horizon[hh]+1]%*%P0G
    # ARRW: Definition 2
    #if(sign.horizon[hh]==Inf) irf.hh <- solve(A0-A0%*%Cm[1:M,]%*%do.call("rbind",rep(list(diag(M)),p)))
    irf.restr[((hh-1)*bigK+1):(bigK*hh),1:bigK] <- irf.hh
  }
  colnames(irf.restr) <- varNames
  rownames(irf.restr) <- paste(rep(varNames,H.restr),".",
                               rep(sign.horizon,each=bigK),sep="")

  Z.cube <- Z.cube[,,shock.order]

  # draw rotation matrix here
  icounter <- 0
  condall <- 0
  max.counter <- MaxTries
  impresp<-Q_bar<-NA
  while(condall == 0 && icounter < max.counter){
    signCheck <- matrix(NA, bigK, 1)
    randMat <- matrix(rnorm(bigK*bigK,0,1),bigK,bigK)
    Q <- matrix(0, bigK, bigK)
    if(no.zero.restr){
      Q <- qr(randMat)
      Q <- qr.Q(Q)
    }else{
      for(mm in 1:bigK){
        Z.temp <- Z.cube[,,mm]
        Z.temp <- Z.temp[rowSums(abs(Z.temp))!=0,,drop=F]
        if(nrow(Z.temp)==0){
          Z.temp <- matrix(0, 1, N.restr)
        }
        if(all(Z.temp==0) && mm>1){
          R <- c()
        }else{
          R <- Z.temp%*%irf.restr
        }
        if(mm > 1){R <- rbind(R, t(Q[,(1:(mm-1)), drop=FALSE]))}

        NU  <- Null(t(R))
        x_j <- randMat[,mm,drop =FALSE]

        q_j <- NU%*%(t(NU)%*%x_j/sqrt(as.numeric(crossprod(t(NU)%*%x_j))))
        Q[,mm] <- q_j
      }
    }
    colnames(Q) <- shock.order; rownames(Q) <- varNames
    Q <- Q[,varNames]
    Q_bar <- Q%*%diag(((diag(Q)>0)-(diag(Q)<0)))

    irf.check <- irf.restr%*%Q_bar
    colnames(irf.check) <- varNames
    rownames(irf.check) <- paste(rep(varNames,H.restr),".",rep(sign.horizon,each=bigK),sep="")

    for(ss in 1:bigK){
      STemp <- S.cube[,,ss]
      IrfCheckTemp <- irf.check[,ss,drop = FALSE]
      signCheckVec <- matrix(NA, N.restr, 1)
      rownames(signCheckVec) <- paste(rep(varNames,H.restr),".",rep(sign.horizon,each=bigK),sep="")
      for(kk in 1:N.restr){
        STempRow <- STemp[kk,]
        emptyCheck <- sum(STempRow)
        if(emptyCheck == 0){
          signCheckVec[kk,1] <- 1;
        }else{
          signCheckVec[kk,1] <- as.numeric(STempRow%*%IrfCheckTemp)
        }
      }
      signCheck[ss,] <- prod((signCheckVec > 0)*(signCheckVec > 0))
    }
    condall <- prod(signCheck)
    icounter <- icounter + 1
  }

  shock <- P0G%*%Q_bar

  irfa  <- array(0,c(n.ahead,bigK,bigK)); dimnames(irfa)[[2]] <- dimnames(irfa)[[3]] <- varNames
  for (ihor in 1:n.ahead){
    irfa[ihor,,] <- PHI[,,ihor]%*%shock
  }

  if(icounter==MaxTries){
    irfa <- NA
  }
  # end rotation matrix loop ----------------------------------------------------------------------------
  return(list(impl=irfa,rot=Q_bar,icounter=icounter))
}

#' @name .irf.chol
#' @noRd
.irf.chol <- function(xdat,plag,n.ahead,Amat,Smat,type,...){
  bigT      <- nrow(xdat)
  bigK      <- ncol(xdat)
  varNames  <- colnames(xdat)
  cons      <- FALSE
  if(nrow(Amat)%%plag!=0) cons <- TRUE

  if(type=="long-run"){
    Asum <- matrix(0,bigK,bigK)
    for(pp in 1:plag){
      Asum <- Asum + Amat[((pp-1)*bigK+1):(pp*bigK),]
    }
    A1inv <- solve(diag(bigK)-t(Asum))
    C1 <- try(t(chol(A1inv%*%Smat%*%t(A1inv))),silent=TRUE)
    if(is(C1,"try-error")) suppressWarnings(C1 <- t(chol(A1inv%*%Smat%*%t(A1inv),pivot=TRUE)))
    shock <- solve(A1inv)%*%C1
  }else if(type=="short-run"){
    shock <- try(t(chol(Smat)),silent=TRUE)
    if(is(shock,"try-error")) suppressWarnings(shock <- t(chol(Smat,pivot=TRUE)))
  }
  colnames(shock) <- rownames(shock) <- varNames

  diagonal <- diag(diag(shock))
  shock    <- solve(diagonal)%*%shock

  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+n.ahead+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  temp <- .gen_compMat(Amat, bigK, plag)
  Cm   <- temp$Cm
  Jm   <- temp$Jm
  Cmat <- diag(bigK*plag)
  for (ihor in (plag+2):(plag+n.ahead+1)){
    Cmat  <- Cmat%*%Cm
    PHIx[,,ihor]  <- t(Jm)%*%Cmat%*%Jm
  }
  PHI  <-  PHIx[,,(plag+1):(plag+n.ahead+1)]

  # computing impulse response function
  irfa  <- array(0,c(n.ahead,bigK,bigK)); dimnames(irfa)[[2]] <- dimnames(irfa)[[3]] <- varNames
  for (ihor in 1:n.ahead){
    irfa[ihor,,] <- PHI[,,ihor]%*%shock
  }

  out <- list(impl=irfa,rot=NULL)
  return(out)
}

#' @name .irf.proxy
#' @noRd
.irf.proxy <- function(xdat,plag,n.ahead,Amat,Smat,Emat,proxy,...){
  bigT      <- nrow(xdat)
  bigK      <- ncol(xdat)
  varNames  <- colnames(xdat)

  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+n.ahead+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  temp <- .gen_compMat(Amat, bigK, plag)
  Cm   <- temp$Cm
  Jm   <- temp$Jm
  Cmat <- diag(bigK*plag)
  for (ihor in (plag+2):(plag+n.ahead+1)){
    Cmat  <- Cmat%*%Cm
    PHIx[,,ihor]  <- t(Jm)%*%Cmat%*%Jm
  }
  PHI  <-  PHIx[,,(plag+1):(plag+n.ahead+1)]

  # identification step
  fitted.err <- lm(Emat[,1] ~ proxy)$fitted
  b21ib11    <- t(lm(Emat[,-1] ~ fitted.err-1)$coef)
  Sig11      <- matrix(Smat[1,1], 1, 1)
  Sig21      <- matrix(Smat[2:bigK,1],bigK-1,1)
  Sig12      <- matrix(Smat[1,2:bigK],1,bigK-1)
  Sig22      <- matrix(Smat[2:bigK,2:bigK],bigK-1,bigK-1)
  ZZp        <- b21ib11%*%Sig11%*%t(b21ib11) - Sig21%*%t(b21ib11)+b21ib11%*%t(Sig21)+Sig22
  b12b12p    <- t(Sig21-b21ib11%*%Sig11)%*%solve(ZZp)%*%(Sig21-b21ib11%*%Sig11)
  b11b11p    <- Sig11 - b12b12p
  b11        <- sqrt(b11b11p)
  if(is.nan(b11)) b11 <- 1e-10
  shock      <- c(b11, b21ib11*c(b11))
  shock      <- shock/shock[1] # normalize to unit shock

  randMat <- matrix(rnorm(bigK*bigK,0,1),bigK,bigK)
  Q <- qr(randMat)
  Q <- qr.Q(Q)
  Q[,1] <- shock

  # computing impulse response function
  irfa  <- array(0,c(n.ahead,bigK,bigK)); dimnames(irfa)[[2]] <- dimnames(irfa)[[3]] <- varNames
  for (ihor in 1:n.ahead){
    irfa[ihor,,] <- PHI[,,ihor]%*%Q
  }

  return(list(impl=irfa,rot=NULL))
}

#' @name .irf.girf
#' @noRd
.irf.girf <- function(xdat,plag,n.ahead,Amat,Smat, ...){
  bigT      <- nrow(xdat)
  bigK      <- ncol(xdat)
  varNames  <- colnames(xdat)

  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+n.ahead+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  temp <- .gen_compMat(Amat, bigK, plag)
  Cm   <- temp$Cm
  Jm   <- temp$Jm
  Cmat <- diag(bigK*plag)
  for (ihor in (plag+2):(plag+n.ahead+1)){
    Cmat  <- Cmat%*%Cm
    PHIx[,,ihor]  <- t(Jm)%*%Cmat%*%Jm
  }
  PHI  <-  PHIx[,,(plag+1):(plag+n.ahead+1)]

  # computing impulse response function
  irfa  <- array(0,c(n.ahead,bigK,bigK)); dimnames(irfa)[[2]] <- dimnames(irfa)[[3]] <- varNames
  for (ihor in 1:n.ahead){
    irfa[ihor,,] <- PHI[,,ihor]%*%Smat
  }

  return(list(impl=irfa,rot=NULL))
}

#' @name .impulsdtrf
#' @noRd
.impulsdtrf <- function(B,smat,nstep)
  ### By:             As emerges from rfvar, neqn x nvar x lags array of rf VAR coefficients.
  ### smat:           nshock x nvar matrix of initial shock vectors.  To produce "orthogonalized
  ###                 impulse responses" it should have the property that crossprod(t(smat))=sigma,
  ###                 where sigma is the Var(u(t)) matrix and u(t) is the rf residual vector.  One
  ###                 way to get such a smat is to set smat=t(chol(sigma)).  To get the smat
  ###                 corresponding to a different ordering, use
  ###                 smat = t(chol(P %*% Sigma %*% t(P)) %*% P), where P is a permutation matrix.
  ###                 To get impulse responses for a structural VAR in the form A(L)y=eps, with
  ###                 Var(eps)=I, use B(L)=-A_0^(-1)A_+(L) (where A_+ is the coefficients on strictly
  ###                 positive powers of L in A), smat=A_0^(-1).
  ###                 In general, though, it is not required that smat be invertible.
### response:       nvar x nshocks x nstep array of impulse responses.
###
### Code written by Christopher Sims,mat based on 6/03 matlab code.  This version 3/27/04.
### Added dimension labeling, 8/02/04.
{

  neq <- dim(B)[1]
  nvar <- dim(B)[2]
  lags <- dim(B)[3]
  dimnB <- dimnames(B)
  if(dim(smat)[2] != dim(B)[2]) stop("B and smat conflict on # of variables")
  response <- array(0,dim=c(neq,nvar,nstep+lags-1));
  response[ , , lags] <- smat
  response <- aperm(response, c(1,3,2))
  irhs <- 1:(lags*nvar)
  ilhs <- lags * nvar + (1:nvar)
  response <- matrix(response, ncol=neq)
  B <- B[, , seq(from=lags, to=1, by=-1)]  #reverse time index to allow matrix mult instead of loop
  B <- matrix(B,nrow=nvar)
  for (it in 1:(nstep-1)) {
    response[ilhs, ] <- B %*% response[irhs, ]
    irhs <- irhs + nvar
    ilhs <- ilhs + nvar
  }
  dim(response) <- c(nvar, nstep + lags - 1, nvar)
  #drop the zero initial conditions; array in usual format
  if(lags>1){
    response<-response[,-(1:(lags-1)),]
  }
  response <- aperm(response, c(1, 3, 2))
  dimnames(response) <- list(dimnB[[1]], dimnames(smat)[[2]], NULL)
  ## dimnames(response)[2] <- dimnames(smat)[1]
  ## dimnames(response)[1] <- dimnames(B)[2]
  return(response)
}

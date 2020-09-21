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
    for(irep in 1:(draws/thin)){
      for(tt in 1:bigT){
        S_store[irep,tt,,] <- bvar$S_store[irep,,]
      }
    }
    L_store     <- NULL
    Smed_store  <- NULL
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
          post_A_tau_prop <- .atau_post(atau=A_tau_prop, thetas=as.vector(theta.lag),lambda2=prod(lambda2_A[1:ss,1]))
          post_A_tau_old  <- .atau_post(atau=A_tau[ss,1], thetas=as.vector(theta.lag),lambda2=prod(lambda2_A[1:ss,1]))
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

#' @name .drawNG
#' @noRd
#' @importFrom GIGrvg rgig
#' @importFrom stats rnorm runif
.drawNG <- function(coef, tau2, lambda2, theta, cl0, dl0, shrink.type = "global",
                    prior=NULL, sample_theta, scale, accept,
                    irep, burnin, c=0, w=0) {
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
    for(kk in 1:K) {
      for(mm in 1:M) {
        tau2[kk,mm] <- GIGrvg::rgig(1,
                                    lambda = theta-0.5,
                                    chi    = (coef[kk,mm]-prior[kk,mm])^2,
                                    psi    = lambda2 * theta)
        # tau2[kk,mm] <- do_rgig1(lambda = theta-0.5,
        #                         chi    = (coef[kk,mm]-prior[kk,mm])^2,
        #                         psi    = lambda2 * theta)
      }
    }
    tau2[tau2<1e-7] <- 1e-7
    if(sample_theta) {
      # Sample theta through a simple RWMH step
      # (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
      theta_prop      <- exp(rnorm(1,0,scale))*theta
      post_theta_prop <- .theta_post(theta = theta_prop,
                                     tau2  = as.vector(tau2),
                                     lambda2 = lambda2)
      post_theta_old  <- .theta_post(theta = theta,
                                     tau2  = as.vector(tau2),
                                     lambda2 = lambda2)
      post.diff       <- post_theta_prop-post_theta_old
      post.diff       <- ifelse(is.nan(post.diff),-Inf,post.diff)
      if (post.diff > log(runif(1,0,1))){
        theta <- theta_prop
        accept <- accept+1
      }
      # Scale MH proposal during the first 50% of the burn-in stage
      if (irep<(0.5*burnin)){
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
          # tau2.lag[kk,mm] <- do_rgig1(lambda = theta[pp]-0.5,
          #                             chi    = (coef.lag[kk,mm]-prior.lag[kk,mm])^2,
          #                             psi    = xi.scale*theta[pp])
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
        post_theta_prop <- .theta_post(theta = theta_prop,
                                       tau2  = as.vector(tau2.lag),
                                       lambda2 = prod(lambda2[1:(pp),1]))
        post_theta_old  <- .theta_post(theta = theta[pp],
                                       tau2  = as.vector(tau2.lag),
                                       lambda2 = prod(lambda2[1:(pp),1]))
        post.diff <- post_theta_prop-post_theta_old
        post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)

        if (post.diff > log(runif(1,0,1))){
          theta[pp] <- theta_prop
          accept[pp] <- accept[pp]+1
        }
        # Scale MH proposal during the first 50% of the burn-in stage
        if (irep<(0.5*burnin)){
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
        # tau2[mm,jj] <- do_rgig1(lambda = theta-0.5,
        #                         chi    = et_0[jj]^2,
        #                         psi    = theta*lambda2)
      }
    }
    tau2[tau2<1e-7 && tau2!=0] <- 1e-7
    if(sample_theta) {
      # Sample theta through a simple RWMH step
      # (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
      theta_prop      <- exp(rnorm(1,0,scale))*theta
      post_theta_prop <- .theta_post(theta = theta_prop,
                                     tau2  = as.vector(tau2),
                                     lambda2 = lambda2)
      post_theta_old  <- .theta_post(theta = theta,
                                     tau2  = as.vector(tau2),
                                     lambda2 = lambda2)
      post.diff       <- post_theta_prop-post_theta_old
      post.diff       <- ifelse(is.nan(post.diff),-Inf,post.diff)
      if (post.diff > log(runif(1,0,1))){
        theta <- theta_prop
        accept <- accept+1
      }
      # Scale MH proposal during the first 50% of the burn-in stage
      if (irep<(0.5*burnin)){
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

#' @name .drawVARcoef
#' @noRd
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
.drawVARcoef <- function(Y,X,Sv,aprior,Vprior,Hprior) {
  M <- ncol(Y)
  bigT <- nrow(Y)
  K <- round(ncol(X))
  if(all(dim(Vprior)==c(K,M))){
    Vinvprior <- matrix(0, K*M, K*M)
    for(mm in 1:M) {
      Vinvprior[((mm-1)*K+1):(mm*K),((mm-1)*K+1):(mm*K)] <- diag(1/Vprior[,mm]) # diag(1/Vprior[,mm])
    }
  }else{
    Vinvprior <- diag(1/diag(Vprior))
  }

  # container
  A   <- matrix(NA, K, M)
  H   <- diag(M)
  eta <- list()
  Em  <- Em.str <- matrix(NA, bigT, M)

  for(mm in 1:M) {
    if(mm==1) {
      Y.i      <- Y[,mm]*exp(-0.5*Sv[,mm])
      X.i      <- X*exp(-0.5*Sv[,mm])
      Vinv.i   <- Vinvprior[1:K,1:K]

      V_post   <- try(chol2inv(chol(crossprod(X.i) + Vinv.i)),silent=TRUE)
      if (is(V_post,"try-error")) V_post <- ginv(crossprod(X.i) + Vinv.i)
      A_post   <- V_post %*% (crossprod(X.i,Y.i) + Vinv.i %*% aprior[,mm])

      A_draw.i <- try(A_post + t(chol(V_post))%*%rnorm(K),silent=T)
      if(is(A_draw.i,"try-error")) A_draw.i <- rmvnorm(1, A_post, V_post)

      A[,mm]   <- A_draw.i
      Em[,mm]  <- Em.str[,mm] <- Y[,mm,drop=F] - X %*% A_draw.i
    }else{
      Y.i      <- Y[,mm]*exp(-0.5*Sv[,mm])
      X.i      <- cbind(X, Em[,1:(mm-1)])*exp(-0.5*Sv[,mm])
      Vinv.i   <- matrix(0, K+(mm-1), K+(mm-1))
      Vinv.i[1:K,1:K] <- Vinvprior[((mm-1)*K+1):(mm*K),((mm-1)*K+1):(mm*K)]
      for(mmm in 1:(mm-1)) Vinv.i[K+mmm,K+mmm] <- 1/Hprior[mm,mmm]

      V_post   <- try(chol2inv(chol(crossprod(X.i) + Vinv.i)), silent=TRUE)
      if(is(V_post,"try-error")) V_post <- ginv(crossprod(X.i) + Vinv.i)
      A_post   <- V_post %*% (crossprod(X.i,Y.i) + Vinv.i %*% c(aprior[,mm], rep(0,mm-1)))

      A_draw.i <- try(A_post + t(chol(V_post)) %*% rnorm(K+(mm-1)), silent=TRUE)
      if(is(A_draw.i,"try-error")) A_draw.i <- rmvnorm(1, A_post, V_post)

      A[,mm]   <- A_draw.i[1:K]
      Em[,mm]  <- Y[,mm,drop=FALSE] - X %*% A_draw.i[1:K]
      Em.str[,mm] <- Y[,mm,drop=FALSE] - X %*% A_draw.i[1:K] -
        Em[,1:(mm-1),drop=FALSE] %*% A_draw.i[(K+1):ncol(X.i),drop=FALSE]
      H[mm,1:(mm-1)] <- eta[[mm]] <- A_draw.i[(K+1):ncol(X.i)]
    }
  }

  return(list(A=A,
              H=H,
              eta=eta,
              Em=Em,
              Em.str=Em.str))
}

#' @name .irf.sign.zero
#' @noRd
#' @importFrom MASS Null
.irf.sign.zero <- function(xdat,plag,n.ahead,Amat,Smat,shock,sign.constr,MaxTries,shock.nr,...){
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
      if(sum(abs(STemp))>0){
        shock.nr <- which(unlist(lapply(sign.constr,function(l)l$shock==dimnames(S.cube)[3][[1]][ss])))
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

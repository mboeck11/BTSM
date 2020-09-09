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

#' @name .theta_post
#' @noRd
#' @importFrom stats dgamma dexp
.theta_post <- function(theta=theta,lambda2=lambda2,tau2=tau2,k=length(tau2),rat=1){
  logpost <- sum(dgamma(tau2,theta,(theta*lambda2/2),log=TRUE))+dexp(theta,rate=rat,log=TRUE)
  return(logpost)
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
.irf.sign.zero <- function(xdat,plag,nhor,Amat,Smat,shock,sign.constr,MaxTries,shock.nr,...){
  bigT     <- nrow(xdat)
  bigK     <- ncol(xdat)
  varNames <- colnames(xdat)

  P0G <- try(t(chol(Smat[,,drop=FALSE])),silent=TRUE)
  if(is(P0G,"try-error")) suppressWarnings(P0G <- t(chol(Smat,pivot=TRUE)))
  colnames(P0G) <- rownames(P0G) <- varNames

  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+nhor+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  temp <- .gen_compMat(Amat, bigK, plag)
  Cm   <- temp$Cm
  Jm   <- temp$Jm
  Cmat <- diag(bigK*plag)
  for (ihor in (plag+2):(plag+nhor+1)){
    Cmat  <- Cmat%*%Cm
    PHIx[,,ihor]  <- t(Jm)%*%Cmat%*%Jm
  }
  PHI  <-  PHIx[,,(plag+1):(plag+nhor+1)]
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

  irfa  <- array(0,c(nhor,bigK,bigK)); dimnames(irfa)[[2]] <- dimnames(irfa)[[3]] <- varNames
  for (ihor in 1:nhor){
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
.irf.chol <- function(xdat,plag,nhor,Amat,Smat,...){
  bigT      <- nrow(xdat)
  bigK      <- ncol(xdat)
  varNames  <- colnames(xdat)

  shock <- try(t(chol(Smat)),silent=TRUE)
  if(is(shock,"try-error")) suppressWarnings(shock <- t(chol(Smat,pivot=TRUE)))
  colnames(shock) <- rownames(shock) <- varNames

  diagonal <- diag(diag(shock))
  shock    <- solve(diagonal)%*%shock

  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+nhor+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  temp <- .gen_compMat(Amat, bigK, plag)
  Cm   <- temp$Cm
  Jm   <- temp$Jm
  Cmat <- diag(bigK*plag)
  for (ihor in (plag+2):(plag+nhor+1)){
    Cmat  <- Cmat%*%Cm
    PHIx[,,ihor]  <- t(Jm)%*%Cmat%*%Jm
  }
  PHI  <-  PHIx[,,(plag+1):(plag+nhor+1)]

  # computing impulse response function
  irfa  <- array(0,c(nhor,bigK,bigK)); dimnames(irfa)[[2]] <- dimnames(irfa)[[3]] <- varNames
  for (ihor in 1:nhor){
    irfa[ihor,,] <- PHI[,,ihor]%*%shock
  }

  out <- list(impl=irfa,rot=NULL)
  return(out)
}

#' @name .irf.girf
#' @noRd
.irf.girf <- function(xdat,plag,nhor,Amat,Smat, ...){
  bigT      <- nrow(xdat)
  bigK      <- ncol(xdat)
  varNames  <- colnames(xdat)

  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+nhor+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  temp <- .gen_compMat(Amat, bigK, plag)
  Cm   <- temp$Cm
  Jm   <- temp$Jm
  Cmat <- diag(bigK*plag)
  for (ihor in (plag+2):(plag+nhor+1)){
    Cmat  <- Cmat%*%Cm
    PHIx[,,ihor]  <- t(Jm)%*%Cmat%*%Jm
  }
  PHI  <-  PHIx[,,(plag+1):(plag+nhor+1)]

  # computing impulse response function
  irfa  <- array(0,c(nhor,bigK,bigK)); dimnames(irfa)[[2]] <- dimnames(irfa)[[3]] <- varNames
  for (ihor in 1:nhor){
    irfa[ihor,,] <- PHI[,,ihor]%*%Smat
  }

  return(list(impl=irfa,rot=NULL))
}

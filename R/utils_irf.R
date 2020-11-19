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
#' @importFrom stats lm
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

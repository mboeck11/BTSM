#' @name .mlag
#' @noRd
.mlag <- function(X,plag){
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,plag*N)
  for (ii in 1:plag){
    Xlag[(plag+1):Traw,(N*(ii-1)+1):(N*ii)] <- X[(plag+1-ii):(Traw-ii),(1:N)]
  }
  colnames(Xlag) <- paste0(colnames(X),".lag",rep(seq(plag),each=N))
  return(Xlag)
}

.mlag_pred <- function(X,plag){
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  if(plag>1){
    Xlag <- matrix(0,Traw,(plag-1)*N)
    for(ii in 1:(plag-1)) Xlag[(plag+1):Traw,(N*(ii-1)+1):(N*ii)] <- X[(plag+1-ii):(Traw-ii),(1:N)]
  }else Xlag <- NULL
  Xpred <- cbind(X,Xlag)
  if(plag>1) colnames(Xpred) <- c(colnames(X),paste0(colnames(X),".lag",rep(seq(plag-1),each=N)))
  return(Xpred)
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
  for(i in 1:M){ # for each equation i
    for(pp in 1:p){ # for each lag pp
      for(j in 1:M){ # for each variable j
        if(i==j){
          V_i[j+M*(pp-1),i] <- (a_bar_1/pp^2) ### variance on own lags
        }else{
          V_i[j+M*(pp-1),i] <- (a_bar_2/pp)^2 * (sigma_sq[i]/sigma_sq[j]) # variance on other lags
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

#' @name .get_V
#' @noRd
.get_V2 <- function(M=M,p=p,Ki1=Ki1,a_bar_1,a_bar_2,a_bar_3,a_bar_4,sigma_sq,sigma_co,cons=FALSE,trend=FALSE){
  V_it <- array(0, c(M*p,M,Ki1))
  Z_i  <- array(0, c(M,M,Ki1))
  if(cons)  C_it <- array(0, c(1,M,Ki1))
  if(trend) T_it <- array(0, c(1,M,Ki1))
  for(kk in 1:Ki1){
    # endogenous part
    for(i in 1:M){ # for each equation i
      for(pp in 1:p){ # for each lag pp
        for(j in 1:M){ # for each variable j
          if(i==j){
            V_it[j+M*(pp-1),i,kk] <- (a_bar_1/pp)^2 ### variance on own lags
          }else{
            V_it[j+M*(pp-1),i,kk] <- (a_bar_2/pp)^2 * (sigma_sq[i]/sigma_sq[j]) # variance on other lags
          }
        }
      }
    }
    # constant
    if(cons){
      for(i in 1:M){
        C_it[1,i,kk] <- a_bar_3 * sigma_sq[i]
      }
    }
    # trend
    if(trend){
      for(i in 1:M){
        T_it[1,i,kk] <- a_bar_3 * sigma_sq[i]
      }
    }
    # zeta
    for(mm in 2:M){
      for(mmm in 1:(mm-1)){
        Z_i[mm,mmm,kk] <- a_bar_4^2 * sigma_sq[i] #(sigma_sq[i]/sigma_co[mm])
      }
    }
  }
  V_i <- lapply(seq(1,Ki1), function(x) V_it[,,x])
  V_i <- Reduce(rbind,V_i)
  if(cons) {
    C_i <- lapply(seq(1,Ki1), function(x) C_it[,,x])
    C_i <- Reduce(rbind,C_i)
  } else C_i <- NULL
  if(trend) {
    T_i <- lapply(seq(1,Ki1), function(x) T_it[,,x])
    T_i <- Reduce(rbind,T_i)
  } else T_i <- NULL

  V_i <- rbind(V_i,C_i,T_i)
  rownames(V_i) <- colnames(V_i) <- NULL

  return(list(theta=V_i,zeta=Z_i))
}



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace RcppArmadillo;

// non-normalized log-density for N(mu, sigma^2)
inline double logdnorm(double x, double mu = 0, double sigma = 1) {
  return -log(sigma)-((x-mu)*(x-mu)/(2*sigma*sigma));
}

//this function implements the griddy gibbs sampler for the thresholds
// [[Rcpp::export]]
NumericVector get_threshold(NumericMatrix Achg1, NumericMatrix SIGMAS1, NumericMatrix SIGMAS2,NumericMatrix threshgrid, NumericMatrix Aapprox){
  int K=Achg1.nrow();
  int grid_length=threshgrid.nrow();
  int Time=Achg1.ncol();
  IntegerVector frame = seq_len(grid_length);
  IntegerVector draw_ind(2);
  arma::vec condpost(Time);
  arma::vec threshold(K);
  arma::mat probs_full(grid_length,K);

  for (int jj=0;jj<K;++jj){
    NumericVector gridj=threshgrid(_,jj);
    arma::vec post1(grid_length,1);
    arma::vec Achg_j=abs(Aapprox(jj,_)); //this is based on the approximation
    arma::vec Achg_actual = abs(Achg1(jj,_));
    for (int kk=0;kk<grid_length;++kk){
      arma::vec condpost(Time);
      double dd=threshgrid(kk,jj); //kk,jj
      for (int tt=0;tt<Time;++tt){ //chg tt=1
        if (Achg_j(tt)>dd){ //chg tt-1
          condpost.row(tt)=logdnorm(Achg_actual(tt),0,sqrt(SIGMAS1(jj,jj)));
        }else{
          condpost.row(tt)=logdnorm(Achg_actual(tt),0,sqrt(SIGMAS2(jj,jj)));
        }
      }
      post1(kk)=sum(condpost);
    }
    probs_full.col(jj)=post1;
    //     NumericVector prob2=Rcpp::as<Rcpp::NumericVector>(wrap(probs1));
    // //
    //      draw_ind=RcppArmadillo::sample(frame,2, FALSE, prob2);//  arma::vec ind = arma::find( den <= rnd);//find(den <= rnd);
    //      int indicator=draw_ind(1);
    //      threshold(jj)=threshgrid(indicator,jj);
  }
  return wrap(probs_full);
}

// [[Rcpp::export]]
double dinvgamma(const double x, const double a, const double b){
  return a * log( b ) - lgamma( a ) - ( a + 1 ) * log( x ) - b / x ;
}
// [[Rcpp::export]]
List KF(NumericMatrix y, NumericMatrix Z,NumericMatrix Ht, NumericMatrix Qtt,int m, int p, int t, NumericVector B0, NumericMatrix V0) {
  static Rcpp::Function asVector("as.vector");

  //Define everything calculated down there
  arma::vec  bp = Rcpp::as<arma::vec>(B0);
  arma::mat Vp = Rcpp::as<arma::mat>(V0);
  arma::mat Qt =Rcpp::as<arma::mat>(Qtt);
  arma::mat bt(t,m);
  arma::mat Vt(pow(m,2),t);
  arma::mat R(p,m);
  arma::mat H(t*m,p);
  arma::mat cfe(p,1);
  arma::mat yt(p,1);
  arma::mat f(p,p);
  arma::mat inv_f(p,p);
  arma::mat btt(m,1);
  arma::mat Vtt(m,m);

  arma::cube test(m,m,t);
  arma::cube Vcov(m,m,t);

  //double log_lik = 0;
  for (int i=1;i<(t+1);++i){
    arma::mat Qt =diag(Qtt(i-1,_));
    //int i=t;
    R= Ht(Range((i-1)*p,(i*p)-1),_);
    H = Z(Range((i-1)*p,(i*p)-1),_);
    yt=y(_,i-1);
    cfe= yt-H*bp;
    f = H*Vp*trans(H)+R;
    inv_f =  solve(f.t(),H).t();//inv(f);
    //btt = bp+Vp*trans(H)*inv_f*cfe; // bp=btt-1
    btt = bp+Vp*inv_f*cfe;
    //Vtt = Vp-Vp*trans(H)*inv_f*H*Vp;
    Vtt = Vp-Vp*inv_f*H*Vp;
    //log_lik=log_lik+
    if ((i-1)<t){
      bp=btt;
      Vp=Vtt+Qt;
    }
    bt.row(i-1) = trans(btt);
    test.slice(i-1)=Vtt;
  }
  //draw S(T|T) ~ N(S(T|T),P(T|T))
  arma::mat bdraw(t,m);

  arma::mat Y = arma::randn(m, 1);
  arma::mat bmean(m,1);
  arma::mat bvar(m,m);

  arma::vec eigval;
  arma::mat eigvec;

  //eig_sym(eigval,eigvec,Vtt,"std");
  arma::eig_sym(eigval,eigvec,Vtt);

  //eigvec*diagmat(sqrt(pmax(eigvec,0)))
  NumericVector epsilon(m);
  arma::vec arma_epsilon(epsilon.begin(),epsilon.length(),false);
  epsilon = rnorm(m);

  bdraw.row(t-1)=(btt+eigvec*diagmat(arma::sqrt(arma::max(eigval,arma::zeros<arma::vec>(m))))*arma_epsilon).t();
  //backward recurssions
  arma::mat bf(1,m);
  for (int i=1;i<t;++i){
    arma::mat Qt =diag(Qtt(t-i,_));
    bf=trans(bdraw.row(t-i));
    btt=trans(bt.row(t-i-1));
    Vtt=test.slice(t-i-1);
    f=Vtt+Qt;
    inv_f=arma::solve(f.t(),Vtt.t()).t();

    cfe=bf-btt;
    //bmean=btt+Vtt*inv_f*cfe;
    bmean=btt+inv_f*cfe;
    //bvar=Vtt-Vtt*inv_f*Vtt;
    bvar=Vtt-inv_f*Vtt;
    Vcov.slice(t-i-1)=bvar;

    //eig_sym(eigval,eigvec,bvar,"std");
    arma::eig_sym(eigval,eigvec,bvar);

    //eigvec*diagmat(sqrt(pmax(eigvec,0)))
    epsilon = rnorm(m);

    bdraw.row(t-i-1)= (bmean+eigvec*diagmat(arma::sqrt(arma::max(eigval,arma::zeros<arma::vec>(m))))*arma_epsilon).t();
  }
  bdraw=trans(bdraw);

  return List::create(
    _["bdraw"]  = wrap(bdraw),
    _["Vcov"] = wrap(Vcov)
  ) ;
}


// [[Rcpp::export]]
double get_lik(NumericMatrix y, NumericMatrix Z,NumericMatrix Ht, NumericMatrix Qtt,int m, int p, int t, NumericVector B0, NumericMatrix V0) {
  static Rcpp::Function asVector("as.vector");

  //Define everything calculated down there
  arma::vec  bp = Rcpp::as<arma::vec>(B0);
  arma::mat Vp = Rcpp::as<arma::mat>(V0);
  arma::mat Qt =Rcpp::as<arma::mat>(Qtt);
  arma::mat bt(t,m);
  arma::mat Vt(pow(m,2),t);
  arma::mat R(p,m);
  arma::mat H(t*m,p);
  arma::mat cfe(p,1);
  arma::mat yt(p,1);
  arma::mat f(p,p);
  arma::mat inv_f(p,p);
  arma::mat btt(m,1);
  arma::mat Vtt(m,m);
  arma::vec Lik(t);


  arma::cube test(m,m,t);

  //double log_lik = 0;
  for (int i=1;i<(t+1);++i){
    arma::mat Qt =diag(Qtt(i-1,_));
    //int i=t;
    R= Ht(Range((i-1)*p,(i*p)-1),_);
    H = Z(Range((i-1)*p,(i*p)-1),_);
    yt=y(_,i-1);
    cfe= yt-H*bp;
    f = H*Vp*trans(H)+R;
    inv_f =inv(f);
    btt = bp+Vp*trans(H)*inv_f*cfe; // bp=btt-1
    Vtt = Vp-Vp*trans(H)*inv_f*H*Vp;
    Lik.row(i-1)=log(arma::det(f))+trans(cfe)*inv_f*cfe;
    if ((i-1)<t){
      bp=btt;
      Vp=Vtt+Qt;
    }
    bt.row(i-1) = trans(bp);
    test.slice(i-1)=Vp;
  }

  return sum(Lik);
}


// [[Rcpp::export]]
List KF_fast(NumericMatrix y, NumericMatrix Z,NumericMatrix Ht, NumericMatrix Qtt,int m, int p, int t, NumericVector B0, NumericMatrix V0) {
  static Rcpp::Function asVector("as.vector");

  //Define everything calculated down there
  arma::vec  bp = Rcpp::as<arma::vec>(B0);
  arma::mat Vp = Rcpp::as<arma::mat>(V0);
  arma::mat Qt =Rcpp::as<arma::mat>(Qtt);
  arma::mat bt(t,m);
  arma::mat Vt(pow(m,2),t);
  arma::mat R(p,m);
  arma::mat H(t*m,p);
  arma::mat cfe(p,1);
  arma::mat yt(p,1);
  arma::mat f(p,p);
  arma::mat inv_f(p,p);
  arma::mat btt(m,1);
  arma::mat Vtt(m,m);

  arma::cube test(m,m,t);
  arma::cube Vcov(m,m,t);

  //double log_lik = 0;
  for (int i=1;i<(t+1);++i){
    arma::mat Qt =diag(Qtt(i-1,_));
    //int i=t;
    R= Ht(Range((i-1)*p,(i*p)-1),_);
    H = Z(Range((i-1)*p,(i*p)-1),_);
    yt=y(_,i-1);
    cfe= yt-H*bp;
    f = H*Vp*trans(H)+R;
    inv_f =  inv(f);
    btt = bp+Vp*trans(H)*inv_f*cfe; // bp=btt-1
    //btt = bp+Vp*inv_f*cfe;
    Vtt = Vp-Vp*trans(H)*inv_f*H*Vp;
    //Vtt = Vp-Vp*trans(H)*inv_f*H*Vp;
    //log_lik=log_lik+
    if ((i-1)<t){
      bp=btt;
      Vp=Vtt+Qt;
    }
    bt.row(i-1) = trans(btt);
    test.slice(i-1)=Vtt;
  }
  //draw S(T|T) ~ N(S(T|T),P(T|T))
  arma::mat bdraw(t,m);

  arma::mat Y = arma::randn(m, 1);
  arma::mat bmean(m,1);
  arma::mat bvar(m,m);

  bdraw.row(t-1)=trans(btt+ trans(arma::chol(Vtt))*Y);
  //backward recurssions
  arma::mat bf(1,m);
  for (int i=1;i<t;++i){
    arma::mat Qt =diag(Qtt(t-i,_));
    bf=trans(bdraw.row(t-i));
    btt=trans(bt.row(t-i-1));
    Vtt=test.slice(t-i-1);
    f=Vtt+Qt;
    inv_f=inv(f);

    cfe=bf-btt;
    bmean=btt+Vtt*inv_f*cfe;
    //bmean=btt+inv_f*cfe;
    bvar=Vtt-Vtt*inv_f*Vtt;
    //bvar=Vtt-inv_f*Vtt;
    Vcov.slice(t-i-1)=bvar;
    arma::mat Y = arma::randn(m, 1);
    bdraw.row(t-i-1)=trans(bmean+trans(arma::chol(bvar))*Y);
  }
  bdraw=trans(bdraw);

  return List::create(
    _["bdraw"]  = wrap(bdraw),
    _["Vcov"] = wrap(Vcov)
  ) ;
}


// [[Rcpp::export]]
List KF_MH(NumericMatrix y, NumericMatrix Z,NumericMatrix Ht, NumericMatrix Qtt,int m, int p, int t, NumericVector B0, NumericMatrix V0) {
  static Rcpp::Function asVector("as.vector");

  //Define everything calculated down there
  arma::vec  bp = Rcpp::as<arma::vec>(B0);
  arma::mat Vp = Rcpp::as<arma::mat>(V0);
  arma::mat Qt =Rcpp::as<arma::mat>(Qtt);
  arma::mat bt(t,m);
  arma::mat Vt(pow(m,2),t);
  arma::mat R(p,m);
  arma::mat H(t*m,p);
  arma::mat cfe(p,1);
  arma::mat yt(p,1);
  arma::mat f(p,p);
  arma::mat inv_f(p,p);
  arma::mat btt(m,1);
  arma::mat Vtt(m,m);

  arma::cube test(m,m,t);
  arma::cube Vcov(m,m,t);

  //double log_lik = 0;
  for (int i=1;i<(t+1);++i){
    arma::mat Qt =diag(Qtt(i-1,_));
    //int i=t;
    R= Ht(Range((i-1)*p,(i*p)-1),_);
    H = Z(Range((i-1)*p,(i*p)-1),_);
    yt=y(_,i-1);
    cfe= yt-H*bp;
    f = H*Vp*trans(H)+R;
    inv_f =  inv(f);
    btt = bp+Vp*trans(H)*inv_f*cfe; // bp=btt-1
    //btt = bp+Vp*inv_f*cfe;
    Vtt = Vp-Vp*trans(H)*inv_f*H*Vp;
    //Vtt = Vp-Vp*trans(H)*inv_f*H*Vp;
    //log_lik=log_lik+
    if ((i-1)<t){
      bp=btt;
      Vp=Vtt+Qt;
    }
    bt.row(i-1) = trans(btt);
    test.slice(i-1)=Vtt;
  }

  return List::create(
    _["bdraw"]  = wrap(bt),
    _["Vcov"] = wrap(test)
  ) ;
}

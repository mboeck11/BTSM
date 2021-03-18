// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @name sample_McCausland
//' @noRd
// [[Rcpp::export]]
arma::mat sample_McCausland(arma::vec ystar, arma::mat Fstar){
  // Import Rs chol function
  Environment base = Environment("package:base");
  Function Rchol = base["chol"];

  // parameters
  int d = Fstar.n_cols;
  int bigT = Fstar.n_rows;

  arma::mat betatilde(bigT+1, d);

  arma::mat I_d = arma::eye(d, d);
  arma::mat Omega_offdiag = -1 * I_d;
  arma::cube Omega_ondiag(d, d, bigT+1);

  Omega_ondiag.slice(0) = 2*I_d;
  for(int tt=1; tt < bigT; tt++){
    Omega_ondiag.slice(tt) = Fstar.row(tt-1).t() * Fstar.row(tt-1) + 2*I_d;
  }
  Omega_ondiag.slice(bigT) = Fstar.row(bigT-1).t() * Fstar.row(bigT-1) + I_d;

  arma::cube cvett(d, 1, bigT+1);
  cvett.slice(0) = arma::vec(d, arma::fill::zeros);
  for(int tt=1; tt < bigT+1; tt++){
    cvett.slice(tt) = Fstar.row(tt-1).t() * ystar(tt-1);
  }

  arma::cube L_upper_arr(d, d, bigT+1);
  arma::cube L_lower_arr(d, d, bigT+1);
  arma::cube steptwo_arr(d, d, bigT+1);
  arma::cube m_list(d, 1, bigT+1);

  arma::mat Sigma_inv = Omega_ondiag.slice(0);
  arma::mat L_upper;
  bool chol_success = chol(L_upper, Sigma_inv);
  // Fall back on Rs chol if armadillo fails (it suppports pivoting)
  if (chol_success == false){
    Rcpp::NumericMatrix tmp = Rchol(Sigma_inv, true, false, -1);
    arma::uvec piv = arma::sort_index(as<arma::vec>(tmp.attr("pivot")));
    arma::mat L_upper_tmp = arma::mat(tmp.begin(), d, d, false);
    L_upper = L_upper_tmp.cols(piv);
  }
  arma::mat L_lower = L_upper.t();
  arma::mat steptwo = solve(arma::trimatl(L_lower), Omega_offdiag);

  L_upper_arr.slice(0) = L_upper;
  L_lower_arr.slice(0) = L_lower;
  steptwo_arr.slice(0) = steptwo;

  arma::mat stepthree = steptwo.t() * steptwo;
  arma::mat a0 = solve(arma::trimatl(L_lower), cvett.slice(0));
  arma::mat m = solve(arma::trimatu(L_upper), a0);
  m_list.slice(0) = m;

  for(int tt = 1; tt < bigT+1; tt++){
    Sigma_inv = Omega_ondiag.slice(tt) - stepthree;
    chol_success = chol(L_upper, Sigma_inv);
    // Fall back on Rs chol if armadillo fails (it suppports pivoting)
    if (chol_success == false){
      Rcpp::NumericMatrix tmp = Rchol(Sigma_inv, true, false, -1);
      arma::uvec piv = arma::sort_index(as<arma::vec>(tmp.attr("pivot")));
      arma::mat L_upper_tmp = arma::mat(tmp.begin(), d, d, false);
      L_upper = L_upper_tmp.cols(piv);
    }
    L_lower = L_upper.t();
    steptwo= solve(arma::trimatl(L_lower), Omega_offdiag);

    L_upper_arr.slice(tt) = L_upper;
    L_lower_arr.slice(tt) = L_lower;
    steptwo_arr.slice(tt) = steptwo;

    stepthree = steptwo.t() * steptwo;
    arma::mat l = cvett.slice(tt) - Omega_offdiag.t() * m_list.slice(tt-1);
    arma::mat at = arma::solve(arma::trimatl(L_lower), l);
    m = arma::solve(arma::trimatu(L_upper), at);
    m_list.slice(tt) = m;
  }

  arma::vec eps = Rcpp::rnorm(d, 0, 1);
  arma::mat l = arma::solve(arma::trimatu(L_upper_arr.slice(bigT)), eps);
  betatilde.row(bigT) = (l + m_list.slice(bigT)).t();

  for (int tt = bigT-1; tt >= 0; tt--){
    eps = Rcpp::rnorm(d, 0, 1);
    arma::vec q = eps - steptwo_arr.slice(tt) * betatilde.row(tt+1).t();
    l = arma::solve(arma::trimatu(L_upper_arr.slice(tt)), q);
    betatilde.row(tt) = (l + m_list.slice(tt)).t();
  }

  return(betatilde);
}

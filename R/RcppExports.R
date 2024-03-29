# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @name BVAR_linear
#' @noRd
BVAR_linear <- function(Y_in, p_in, draws_in, burnin_in, cons_in, trend_in, sv_in, thin_in, prior_in, hyperparam_in, Ex_in) {
    .Call(`_BTSM_BVAR_linear`, Y_in, p_in, draws_in, burnin_in, cons_in, trend_in, sv_in, thin_in, prior_in, hyperparam_in, Ex_in)
}

do_rgig1 <- function(lambda, chi, psi) {
    .Call(`_BTSM_do_rgig1`, lambda, chi, psi)
}

#' @name dmvnrm_arma_fast
#' @noRd
dmvnrm_arma_fast <- function(x, mean, sigma, logd = FALSE) {
    .Call(`_BTSM_dmvnrm_arma_fast`, x, mean, sigma, logd)
}

#' @name loglik_C
#' @noRd
loglik_C <- function(Y_in, X_in, A_in, S_in, thindraws_in) {
    .Call(`_BTSM_loglik_C`, Y_in, X_in, A_in, S_in, thindraws_in)
}

#' @name sample_McCausland
#' @noRd
sample_McCausland <- function(ystar, Fstar) {
    .Call(`_BTSM_sample_McCausland`, ystar, Fstar)
}

#' @name MH_step
#' @noRd
MH_step <- function(current_val, c_tuning_par, d, scale_par, param_vec, b, nu, hyp1, hyp2) {
    .Call(`_BTSM_MH_step`, current_val, c_tuning_par, d, scale_par, param_vec, b, nu, hyp1, hyp2)
}

get_threshold <- function(Achg1, SIGMAS1, SIGMAS2, threshgrid, Aapprox) {
    .Call(`_BTSM_get_threshold`, Achg1, SIGMAS1, SIGMAS2, threshgrid, Aapprox)
}

dinvgamma <- function(x, a, b) {
    .Call(`_BTSM_dinvgamma`, x, a, b)
}

KF <- function(y, Z, Ht, Qtt, m, p, t, B0, V0) {
    .Call(`_BTSM_KF`, y, Z, Ht, Qtt, m, p, t, B0, V0)
}

get_lik <- function(y, Z, Ht, Qtt, m, p, t, B0, V0) {
    .Call(`_BTSM_get_lik`, y, Z, Ht, Qtt, m, p, t, B0, V0)
}

KF_fast <- function(y, Z, Ht, Qtt, m, p, t, B0, V0) {
    .Call(`_BTSM_KF_fast`, y, Z, Ht, Qtt, m, p, t, B0, V0)
}

KF_MH <- function(y, Z, Ht, Qtt, m, p, t, B0, V0) {
    .Call(`_BTSM_KF_MH`, y, Z, Ht, Qtt, m, p, t, B0, V0)
}





#' The mean resultant length of a von Mises distribution.
#'
#' @param kp The concentration parameter.
#'
#' @return The mean resultant length.
#' @export
#'
computeMeanResultantLengthVM <- function(kp) {
  besselI(kp, 1) / besselI(kp, 0)
}

#' The mean resultant length of a Batschelet-type distribution.
#'
#' @param kp The concentration parameter.
#' @param lam The peakedness parameter.
#' @param dbat_fun The pdf function of the required Batschelet distribution.
#'
#' @return The mean resultant length.
#' @export
computeMeanResultantLengthBat <- function(kp, lam, dbat_fun = dinvbat) {

  # R = E[cos(theta)], so we need a function f(theta) = cos(theta)
  # p(theta, kp, lam) to integrate over.
  cos_fun <- function(theta) cos(theta) * dbat_fun(theta, 0, kp, lam, log = FALSE)

  integrate(cos_fun, -pi, pi)$value
}

computeCircSD <- function(R_bar) sqrt(-2*log(R_bar))
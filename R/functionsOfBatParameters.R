


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
computeMeanResultantLengthBat <- function(kp, lam, bat_type = "inverse") {

  # Choose the appropriate functions for the bat_type.
  if (bat_type == "inverse") {
    dbatkernfun <- dinvbatkern
    nc_batfun <- K_kplam
  } else if (bat_type == "power") {
    dbatkernfun <- dpowbatkern
    nc_batfun <- powbat_nc
  } else {
    stop("Unknown Batschelet-type.")
  }

  # Use the von Mises resultant length when possible.
  if (lam == 0) computeMeanResultantLengthVM(kp)

  # R = E[cos(theta)], so we need a function f(theta) = cos(theta) p(theta, kp,
  # lam) to integrate over. Note that the normalizing constant is removed for
  # computational efficiency.
  cos_fun <- function(theta) {
    cos(theta) * dbatkernfun(theta, 0, kp, lam, log = FALSE)
  }

  stats::integrate(cos_fun, -pi, pi)$value / nc_batfun(kp, lam)
}



computeCircSD <- function(R_bar) sqrt(-2*log(R_bar))








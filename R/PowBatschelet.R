
# The original transformation by Batschelet 1981
# Note that this version does not include peaked distributions.
tau_lam <- function(x, lam) x + lam * sin(x)

# New power transformations
tpow_lam     <- function(x, lam) sign(x) * pi * (abs(x) / pi)^exp(lam)
tpow_lam_inv <- function(x, lam) sign(x) * pi * (abs(x) / pi)^exp(-lam)


#' Kernel of the von-Mises based symmetric power Batschelet distribution
#'
#' @return The unnormalized density value of the power Batschelet distribution.
#'
dpowbatkern <- function(x, mu = 0, kp = 1, lam = 0, log = FALSE) {
  dvm(tpow_lam(x - mu, lam), mu = 0, kp = kp, log = log)
}


powbat_nc <- function(kp, lam) {
  integrate(function(x) dpowbatkern(x, mu = 0, kp, lam, log = FALSE), -pi, pi)$value
}



#' The von-Mises based symmetric power Batschelet distribution
#'
#' The power Batschelet distribution, with mean direction \code{mu}, concentration parameter \code{kp}, and
#' shape (peakedness) parameter \code{lam}. This is the von Mises based version, without a skewness
#' parameter.
#'
#' @param x An angle in radians.
#' @param mu A mean direction, in radians.
#' @param kp Numeric, \eqn{> 0,}the concentration parameter.
#' @param lam The shape parameter (peakedness), -Inf < \code{lam} < Inf
#' @param log Logical; whether to return the log of the probability or not.
#'
#' @return
#' @export
#'
#' @examples
#' dinvbat(3)
#'
#' # Peaked distribution
#' curve(dinvbat(x, lam = .8), -pi, pi)
#'
#' # Flat-topped distribution
#' curve(dinvbat(x, lam = -.8), -pi, pi)
#'
dpowbat <- function(x, mu = 0, kp = 1, lam = 0, log = FALSE) {
  if (kp < 0) return(NA)

  if (log) {
    dpowbatkern(x, mu = mu, kp = kp, lam = lam, log = TRUE) /
  } else {
    dpowbatkern(x, mu = mu, kp = kp, lam = lam, log = FALSE)
  }
}



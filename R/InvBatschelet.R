#' Compute a partial normalizing constant of a vM-based inverse Batschelet distribution
#'
#' The value givin is only the additional part of the inverse Batschelet. The von Mises normalizing
#' constant must still be added.
#'
#' @inheritParams dinvbat
#'
#' @return Numeric.
#'
K_kplam <- function(kp, lam) {
    (1 + lam) / (1 - lam) - (2 * lam / (1 - lam)) *
      integrate(function(x) dvm(x - 0.5 * (1 - lam) * sin(x), 0, kp), -pi, pi)$value
}

#' Set an angle to have its bounds in (-pi, pi)
#'
#' This function is used to make sure that an angle is represented as being between -pi and pi.
#'
#' @param x Numeric, representing an angle.
#'
#' @return  Numeric, representing an angle, between -pi and pi.
#'
force_neg_pi_pi <- function(x) {
  ((x + pi) %% (2*pi)) - pi
}

#' Inverse Batschelet subfunction
#'
#' @inheritParams t_lam
#'
#' @return Numeric.
#'
s_lam <- function(x, lam) {
  x - 0.5 * (1 + lam) * sin(x)
}

s_lam_inv <- function(x, lam) {
  # Compute the root to obtain the inverse.
  uniroot(function(y) s_lam(y, lam) - x, lower = -pi, upper = pi)$root
}

#' Inverse Batschelet function
#'
#' Function used to generate different shapes for the inverse Batschelet distribution.
#'
#' @param x An angle in radians.
#' @param lam The shape parameter.
#'
#' @return The transformed angle x.
#' @export
#'
#' @examples
#'
#' # Reproduce Figure 5 in Jones & Pewsey
#' curve(t_lam(x, 0), -pi, pi, asp = 1)
#' for (i in seq(-1, 1, by = 1/3)) {
#'   curve(t_lam(x, i), -pi, pi, add = TRUE, lty = "dotted")
#' }
t_lam <- Vectorize(function(x, lam) {

  # Make sure the the input x is in the correct bounds.
  if (x > pi || x < -pi) x <- force_neg_pi_pi(x)

  if (lam != -1) {
    x * (1 - lam) / (1 + lam) + s_lam_inv(x, lam) * 2 * lam / (1 + lam)
  } else {
    x - sin(x)
  }
})


#' Kernel of the von-Mises based symmetric inverse Batschelet distribution
#'
#' @inheritParams dinvbat
#'
#' @return The unnormalized density value of the inverse Batschelet distribution.
#'
dinvbatkern <- function(x, mu = 0, kp = 1, lam = 0, log = FALSE) {
  dvm(t_lam(x - mu, lam), mu = 0, kp = kp, log = log)
}


#' The von-Mises based symmetric inverse Batschelet distribution
#'
#' The inverse Batschelet distribution, with mean direction \code{mu}, concentration parameter \code{kp}, and
#' shape (peakedness) parameter \code{lam}. This is the von Mises based version, without a skewness
#' parameter.
#'
#' @param x An angle in radians.
#' @param mu A mean direction, in radians.
#' @param kp Numeric, \eqn{> 0,}the concentration parameter.
#' @param lam The shape parameter (peakedness), -1 < \code{lam} < 1.
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
dinvbat <- function(x, mu = 0, kp = 1, lam = 0, log = FALSE) {
  dinvbatkern(x, mu = mu, kp = kp, lam = lam, log = log) / K_kplam(kp = kp, lam = lam)
}




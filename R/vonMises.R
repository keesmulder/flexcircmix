#' Log of the Bessel function
#'
#' This can be more numerically stable than taking the log directly.
#'
#' @inheritParams base::besselI
#'
#' @return The logarithm of the bessel funciton
#'
#' @examples
logBesselI <- function(x, nu) {
  x + log(besselI(x, nu, expon.scaled = TRUE))
}


#' @describeIn dvm The kernel of the von Mises distribution.
dvmkern <- function(x, mu = 0, kp = 1, log = FALSE) {
  logp <- kp * cos(x - mu)

  if (log) {
    return(logp)
  } else {
    return(exp(logp))
  }
}

#' Compute the density of the von Mises distribution
#'
#' Compute the probability density function of the von Mises distribution.
#'
#' @param x An angle in radians for which to compute the pdf.
#' @param mu The mean direction.
#' @param kp The concentration parameter, \eqn{kappa}.
#' @param log Logical; If TRUE, return the log of the probability of x.
#'
#' @return The probability density of the von Mises distribution at x.
#' @export
#'
#' @examples
#' dvm(3)
#' dvm(3, -2.2, 10)
dvm <- function(x, mu = 0, kp = 1, log = FALSE) {
  logp <- kp * cos(x - mu) - log(2 * pi) - logBesselI(kp, 0)

  if (log) {
    return(logp)
  } else {
    return(exp(logp))
  }
}







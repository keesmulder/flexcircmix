#' Log of the Bessel function
#'
#' This can be more numerically stable than taking the log directly.
#'
#' @inheritParams base::besselI
#'
#' @return The logarithm of the bessel function
#'
logBesselI <- function(x, nu) {
  x + log(besselI(x, nu, expon.scaled = TRUE))
}


#' Kernel of the von Mises distribution.
#'
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


#' Bessel exponential distribution
#'
#' Functions for the Bessel exponential distribution. This is the conditional
#' posterior distribution of the concentration parameter \code{kappa} of a von
#' Mises distribution with conjugate prior. The random generation algorithm is
#' due to Forbes and Mardia (2015).
#'
#' @param kp Numeric; value of kappa to evaluate.
#' @param n Integer; number of generated samples.
#' @param log Logical; Whether to return the log of the result.
#' @param eta Integer; This is the posterior sample size, which is n + c where c
#'   is the number of observations contained in the conjugate prior. For
#'   uninformative, \code{c = 0} and \code{eta = n}.
#' @param g Numeric; Should be \code{-R*cos(mu-theta_bar)/eta}, where \code{R}
#'   is the posterior mean resultant length, and \code{theta_bar} is the
#'   posterior mean, while \code{mu} is the current value of the mean. Note that
#'   this parameter is called \eqn{\beta_0} in Forbes and Mardia (2015).
#'
#' @name besselexp
#'
#' @return For \code{dbesselexp} and \code{dbesselexpkern}, a scalar. For
#'   \code{rbesselexp}, a vector of random variates from the distribution.
#'
#' @examples
#' dbesselexp(2, 20, -.5)
#' plot(density(rbesselexp(100, 20, -.5)), xlim = c(0, 5))
#'
#' # Plot probability density function
#' dbesexpfun <- Vectorize(function(x) dbesselexp(x, 20, -.5))
#' curve(dbesexpfun, 0, 5)
#'
#' # Compare with density of random draws
#' plot(density(rbesselexp(1000, 20, -.5)), xlim = c(0, 5))
#'
NULL

#' @describeIn besselexp Probability density function.
#' @export
dbesselexp <- function(kp, eta = 1, g = -.5, log = FALSE) {
  nc <- stats::integrate(function(x) dbesselexpkern(x, eta, g, log = log),
                  0, Inf)$value
  dbesselexpkern(kp, eta, g) / nc
}


#' @describeIn besselexp Kernel (unnormalized version) of the pdf.
#' @export
dbesselexpkern <- Vectorize(function(kp, eta, g, log = FALSE) {
  logprob <- -eta * g * kp - eta * logBesselI(kp, 0)
  ifelse(log, logprob, exp(logprob))
}, "kp")


#' @describeIn besselexp Random generation.
#' @export
rbesselexp <- function(n, eta, g) {
  replicate(n, circglmbayes::sampleKappa(eta * g, eta)[1])
}


















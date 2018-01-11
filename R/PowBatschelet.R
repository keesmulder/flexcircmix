
# The original transformation by Batschelet 1981
# Note that this version does not include peaked distributions.
tau_lam <- function(x, lam) x + lam * sin(x)

# Logistic functions
expit <- function(x) exp(x) / (1 + exp(x)) # R -> [0, 1]
logit <- function(p) log(p) - log(1 - p)   # [0, 1] -> R


# New power transformations
tpow_lam     <- function(x, lam) {

  x <- force_neg_pi_pi(x)

  # Reparametrize from (-1, 1) to a power for the power batschelet function.
  pwr <- (1 - 0.4052284*lam) / (1 + 0.4052284*lam)

  sign(x) * pi * (abs(x) / pi)^pwr
}

tpow_lam_inv <- function(x, lam) {

  x <- force_neg_pi_pi(x)

  # Reparametrize from (-1, 1) to a power for the power batschelet function.
  pwr <- (1 + 0.4052284*lam) / (1 - 0.4052284*lam)

  sign(x) * pi * (abs(x) / pi)^pwr
}




### THESE DERIVATIVES MUST BE CHECKED
# Derivatives of the power function with respect to x.
# tpow_lam_inv_d <- function(x, lam) {
#
#   pwr <- (1 - lam/2) / (1 + lam/2)
#   pi^(1 + pwr) * sign(x) * x * abs(x)^(- 2 - 1 / pwr)
# }
# # Derivatives of the power function
# tpow_lam_inv_d <- function(x, lam) {
#   pwr <- (1 - lam) / (1 + lam)
#   x * pi^(1 - pwr) * pwr * sign(x) * abs(x)^(pwr - 2)
# }


#' Kernel of the von-Mises based symmetric power Batschelet distribution
#'
#' @inheritParams dpowbat
#'
#' @return The unnormalized density value of the power Batschelet distribution.
#'
dpowbatkern <- function(x, mu = 0, kp = 1, lam = 0, log = FALSE) {
  dvm(tpow_lam(x - mu, lam), mu = 0, kp = kp, log = log)
}

# Compute the normalizing constant of the power Batschelet function.
powbat_nc <- function(kp, lam) {

  # Try to use the usual numerical integration, which usually doesn't fail. With
  # some values of lambda, particularly very high, this might fail and the
  # tolerance can be lowered for a sensible answer.
  tryCatch(
    stats::integrate(function(x) dpowbatkern(x, mu = 0, kp, lam, log = FALSE), -pi, pi)$value,

    error = function(e) {
      return(stats::integrate(function(x) dpowbatkern(x, mu = 0, kp, lam, log = FALSE), -pi, pi,
                       abs.tol = 1)$value)
    })
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
#' @return Numeric; Either the probability or log-probability of angle x given the parameters.
#' @export
#'
#' @examples
#' dpowbat(3)
#'
#' # Flat-topped distribution
#' curve(dpowbat(x, lam = -.8), -pi, pi)
#'
#' # Peaked distribution
#' curve(dpowbat(x, lam = .8), -pi, pi)
#'
dpowbat <- function(x, mu = 0, kp = 1, lam = 0, log = FALSE) {
  if (kp < 0 || lam < -1 || lam > 1) return(NA)

  if (log) {
    dpowbatkern(x, mu = mu, kp = kp, lam = lam, log = TRUE)  - log(powbat_nc(kp, lam))
  } else {
    dpowbatkern(x, mu = mu, kp = kp, lam = lam, log = FALSE) / powbat_nc(kp, lam)
  }
}



#' Likelihood of the power Batschelet distribution
#'
#' Two functions to obtain the power Batschelet likelihood of a set parameters,
#' given a data set of angles.
#'
#' @param x A set of angles in radians.
#' @param weights A vector of length \code{length(x)}, which gives importace
#'   weights to be used for \code{x}.
#' @param log If \code{TRUE} (the default), the log-likelihood is used.
#' @param mu A mean direction, in radians.
#' @param kp Numeric, \eqn{> 0,}the concentration parameter.
#' @param lam The shape parameter (peakedness), -1 < \code{lam} < 1.
#'
#' @return \code{likpowbat} returns a value, the likelihood given the data and
#'   parameters. \code{likfunpowbat} returns a function of mu, kp and lam, which
#'   can be evaluated later for a given set of parameters.
#' @export
#'
#' @examples
#' vm_data <- circglmbayes::rvmc(10, 1, 5)
#' llfun <- likfunpowbat(vm_data)
#'
#' # log-likelihood value of true parameters.
#' llfun(mu = 1, kp = 5, lam = 0)
#'
#' # Plot the conditional log-likelihood.
#' kp_conditional_ll <- Vectorize(function(x) llfun(mu = 1, kp = x, lam = 0))
#' curve(kp_conditional_ll, 0, 20)
likpowbat <- function(x, mu, kp, lam, weights = rep(1, length(x)), log = TRUE) {
  if (log) {
    sum(weights * dpowbat(x, mu, kp, lam, log = TRUE))
  } else {
    exp(sum(weights * dpowbat(x, mu, kp, lam, log = TRUE)))
  }
}

#' @describeIn likpowbat Return a likelihood function.
#' @export
#'
likfunpowbat <- function(x, weights = rep(1, length(x)), log = TRUE) {
  if (log) {
    function(mu, kp, lam) sum(weights * dpowbat(x, mu, kp, lam, log = TRUE))
  } else {
    function(mu, kp, lam) exp(sum(weights * dpowbat(x, mu, kp,
                                                    lam, log = TRUE)))
  }
}



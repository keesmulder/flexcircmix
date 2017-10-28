
# The original transformation by Batschelet 1981
# Note that this version does not include peaked distributions.
tau_lam <- function(x, lam) x + lam * sin(x)

# Logistic functions
expit <- function(x) exp(x) / (1 + exp(x)) # R -> [0, 1]
logit <- function(p) log(p) - log(1 - p)   # [0, 1] -> R


# New power transformations
tpow_lam     <- function(x, lam) {

  # Reparametrize from (-1, 1) to a power for the power batschelet function.
  pwr <- (1 + lam) / (1 - lam)

  sign(x) * pi * (abs(x) / pi)^pwr
}

tpow_lam_inv <- function(x, lam) {

  # Reparametrize from (-1, 1) to a power for the power batschelet function.
  pwr <- (1 - lam) / (1 + lam)

  sign(x) * pi * (abs(x) / pi)^pwr
}




### THESE DERIVATIVES MUST BE CHECKED
# Derivatives of the power function with respect to x.
tpow_lam_inv_d <- function(x, lam) {

  pwr <- (1 - lam) / (1 + lam)
  pi^(1 + pwr) * sign(x) * x * abs(x)^(- 2 - 1 / pwr)
}
# Derivatives of the power function
tpow_lam_inv_d <- function(x, lam) {
  pwr <- (1 - lam) / (1 + lam)
  x * pi^(1 - pwr) * pwr * sign(x) * abs(x)^(pwr - 2)
}


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
#' dpowbat(3)
#'
#' # Peaked distribution
#' curve(dpowbat(x, lam = -.8), -pi, pi)
#'
#' # Flat-topped distribution
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



#' Obtain the likelihood of a power Batschelet distribution
#'
#' @param x An set of angles in radians.
#' @param weights A vector of length \code{length(x)}, which gives importace weights to be used for x.
#' @param log If \code{TRUE} (the default), the log-likelihood is used.
#' @param mu A mean direction, in radians.
#' @param kp Numeric, \eqn{> 0,}the concentration parameter.
#' @param lam The shape parameter (peakedness), -1 < \code{lam} < 1.
#'
#' @return \code{likpowbat} returns a value, the likelihood given the data and parameters.
#'   \code{likfunpowbat} returns a function of mu, kp and lam, which can be evaluated later for a
#'   given set of parameters.
#' @export
#'
#' @examples
#'
#'
likfunpowbat <- function(x, weights = rep(1, length(x)), log = TRUE) {
  if (log) {
    function(mu, kp, lam) sum(weights * dpowbat(x, mu, kp, lam, log = TRUE))
  } else {
    function(mu, kp, lam) exp(sum(weights * dpowbat(x, mu, kp, lam, log = FALSE)))
  }
}

#' @describeIn likfunpowbat
likpowbat <- function(x, mu, kp, lam, weights = rep(1, length(x)), log = TRUE) {
  if (log) {
    sum(weights * dpowbat(x, mu, kp, lam, log = TRUE))
  } else {
    exp(sum(weights * dpowbat(x, mu, kp, lam, log = FALSE)))
  }
}



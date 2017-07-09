#' Obtain the likelihood of an inverse Batschelet distribution
#'
#' @param x An set of angles in radians.
#' @param log If \code{TRUE} (the default), the log-likelihood is used.
#' @param mu A mean direction, in radians.
#' @param kp Numeric, \eqn{> 0,}the concentration parameter.
#' @param lam The shape parameter (peakedness), -1 < \code{lam} < 1.
#'
#' @return \code{likinvbat} returns a value, the likelihood given the data and parameters.
#'   \code{likfuninvbat} returns a function of mu, kp and lam, which can be evaluated later for a
#'   given set of parameters.
#' @export
#'
#' @examples
#' x <- rinvbat(5)
#'
#' # Find the likelihood
#' likinvbat(x, mu = 0, kp = 1, lam = 0.1, log = TRUE)
#'
#' # likfuninvbat returns a function.
#' lfib <- likfuninvbat(x)
#' lfib(mu = 0, kp = 1, lam = 0.1, log = TRUE)
#'
likfuninvbat <- function(x, log = TRUE) {
  if (log) {
    function(mu, kp, lam) sum(dinvbat(x, mu, kp, lam, log = TRUE))
  } else {
    exp(function(mu, kp, lam) sum(dinvbat(x, mu, kp, lam, log = TRUE)))
  }
}

#' @describeIn likfuninvbat
likinvbat <- function(x, mu, kp, lam, log = TRUE) {
  if (log) {
    sum(dinvbat(x, mu, kp, lam, log = TRUE))
  } else {
    exp(sum(dinvbat(x, mu, kp, lam, log = TRUE)))
  }
}


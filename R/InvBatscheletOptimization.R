#' Find maximum likelihood estimates for the inverse Batschelet distribution
#'
#' @param x An set of angles in radians.
#' @param fixed_mu If NA, the \code{mu} will be estimated as the mean direction. Else, \code{mu}
#'   will be fixed at \code{fixed_mu}.
#'
#' @return The maximum likelihood estimates for the \code{mu}, \code{kp}, and \code{lam}.
#'
maxlikinvbat <- function(x, fixed_mu = NA) {

  llfib <- likfuninvbat(x, log = TRUE)

  # Maximum likelihood estimate for the mean
  if (is.na(fixed_mu)) {
    mu_hat <- atan2(sum(sin(x)), sum(cos(x)))
  } else {
    mu_hat <- fixed_mu
  }

  # Find maximum likelihood estimates for the other parameters.
  om <- optim(fn = function(params) -llfib(mu = mu_hat, kp = params[1], lam = params[2]),
              method = "Nelder-Mead", control = list(trace = 0), par = c(1, 0))

  c(mu = mu_hat, kp = om$par[1], lam = om$par[2])
}

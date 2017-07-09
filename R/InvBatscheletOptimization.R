#' Find maximum likelihood estimates for the inverse Batschelet distribution
#'
#' @param x An set of angles in radians.
#' @param fixed_mu If NA, the \code{mu} will be estimated as the mean direction. Else, \code{mu}
#'   will be fixed at \code{fixed_mu}.
#'
#' @return The maximum likelihood estimates for the \code{mu}, \code{kp}, and \code{lam}.
#'
maxlikinvbat <- function(x, fixed_mu = NA, fixed_kp = NA, fixed_lam = NA) {

  llfib <- likfuninvbat(x, log = TRUE)

  # Maximum likelihood estimate for the mean
  if (is.na(fixed_mu)) {
    mu_hat <- atan2(sum(sin(x)), sum(cos(x)))
  } else {
    mu_hat <- fixed_mu
  }


  # If kp or lam is fixed, we only need univariate optimization.
  if (!is.na(fixed_kp)) {

    lam_hat <- optimize(f = function(lam) llfib(mu = mu_hat, kp = fixed_kp, lam = lam),
                        lower = -1, upper = 1, maximum = TRUE)$maximum
    return(c(mu = mu_hat, kp = fixed_kp, lam = lam_hat))

  } else if (!is.na(fixed_lam)) {

    kp_hat <- optim(fn = function(kp) -llfib(mu = mu_hat, kp = kp, lam = fixed_lam),
                    par = 1, lower = 0, upper = Inf, method = "L-BFGS-B")$par[1]

    return(c(mu = mu_hat, kp = kp_hat, lam = fixed_lam))

  # Else, maximize for both parameters.
  } else {

    # Find maximum likelihood estimates for the both parameters.
    om <- optim(fn = function(params) -llfib(mu = mu_hat, kp = params[1], lam = params[2]),
                method = "Nelder-Mead", control = list(trace = 0), par = c(1, 0))

    return(c(mu = mu_hat, kp = om$par[1], lam = om$par[2]))
  }
}

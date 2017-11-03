#' Find maximum likelihood estimates for the inverse Batschelet distribution
#'
#' @param x An set of angles in radians.
#' @param fixed_mu If NA, the \code{mu} will be estimated as the mean direction.
#'   Else, \code{mu} will be fixed at \code{fixed_mu}.
#' @param fixed_kp If NA, kp will be estimated. Else, a numeric giving the fixed
#'   value of kp.
#' @param fixed_lam  If NA, lam will be estimated. Else, a numeric giving the
#'   fixed value of lam.
#' @param weights A vector of length \code{length(x)}, which gives importace
#'   weights to be used for x.
#' @param max_its Used only if neither \code{fixed_kp} nor \code{fixed_lam} is
#'   provided. The maximum number of iterations for the Nelder-Mead optimization.
#' @param kp_max Used only if \code{fixed_lam} is provided. The maximum value
#'   for kappa to search for.
#'
#' @return The maximum likelihood estimates for the \code{mu}, \code{kp}, and
#'   \code{lam}.
#'
maxlikbat <- function(x, likfunbat_fun = likfuninvbat, weights, fixed_mu = NA, fixed_kp = NA, fixed_lam = NA,
                         max_its = 20, kp_max = 100) {

  # Default weights if not supplied.
  if (missing(weights)) {
    llfib <- likfunbat_fun(x, log = TRUE)
    weights <- rep(1, length(x))
  } else {
    llfib <- likfunbat_fun(x, weights = weights, log = TRUE)
  }

  # Maximum likelihood estimate for the mean
  if (is.na(fixed_mu)) {
    mu_hat <- atan2(sum(weights*sin(x)), sum(weights*cos(x)))
  } else {
    mu_hat <- fixed_mu
  }


  # If kp or lam is fixed, we only need univariate optimization.
  if (!is.na(fixed_kp)) {

    lam_hat <- optimize(f = function(lam) llfib(mu = mu_hat, kp = fixed_kp, lam = lam),
                        lower = -1, upper = 1, maximum = TRUE)$maximum
    return(c(mu = mu_hat, kp = fixed_kp, lam = lam_hat))

  } else if (!is.na(fixed_lam)) {

    kp_hat <- optimize(f = function(kp) llfib(mu = mu_hat, kp = kp, lam = fixed_lam),
                       lower = 0, upper = kp_max, maximum = TRUE)$maximum

    return(c(mu = mu_hat, kp = kp_hat, lam = fixed_lam))

  # Else, maximize for both parameters.
  } else {

    # Find maximum likelihood estimates for the both parameters.
    om <- optim(fn = function(params) -llfib(mu = mu_hat, kp = params[1], lam = params[2]),
                method = "Nelder-Mead", control = list(trace = 0, maxit = 5*max_its),
                par = c(1, 0))

    return(c(mu = mu_hat, kp = om$par[1], lam = om$par[2]))
  }
}

#' Fit mixtures of inverse Batschelet distributions.
#'
#' @param x A dataset of angles in radians.
#' @param n_comp The number of components to be used in the mixture. This is fixed, so it can not be
#'   estimated.
#' @param init_pmat An \code{n_comp * 4} matrix, the initial values of the parameter matrix. The
#'   parameters are ordered \code{mu}, \code{kp}, \code{lam} and then \code{alph}.
#' @param fixed_pmat An \code{n_comp * 4} matrix, containing a parameter matrix, with \code{NA} for
#'   parameters to be estimated, and a numeric for each parameter that should be kept fixed to a
#'   specific value.
#' @param ll_tol Numeric; the algorithm is stopped after the log-likelihood improves less than \code{ll_tol}.
#' @param max_its The maximum number of iterations.
#' @param verbose Logical; whether to print debug statements.
#'
#' @return A parameter matrix of results.
#' @export
#'
#' @examples
#'
#'
fitinvbatmix <- function(x, n_comp  = 4,
                         init_pmat  = matrix(NA, n_comp, 4),
                         fixed_pmat = matrix(NA, n_comp, 4),
                         ll_tol = 1,
                         max_its = 50,
                         verbose = FALSE) {

  # Force x to be in range -pi, pi.
  x <- force_neg_pi_pi(x)

  n <- length(x)

  # Set initial values if the initial parameter matrix is not given for that parameter (has NAs).
  if (any(is.na(init_pmat[, 1]))) init_pmat[, 1] <- seq(0, 2*pi, length.out = n_comp + 1)[-1]
  if (any(is.na(init_pmat[, 2]))) init_pmat[, 2] <- rep(5, n_comp)
  if (any(is.na(init_pmat[, 3]))) init_pmat[, 3] <- rep(0, n_comp)
  if (any(is.na(init_pmat[, 4]))) init_pmat[, 4] <- rep(1/n_comp, n_comp)

  # The current parameter matrix.
  pmat_cur <- init_pmat
  colnames(pmat_cur) <- c("mu", "kp", "lam", "alph")

  # initialize W matrix.
  W <- matrix(1/n_comp, nrow = n, ncol = n_comp)

  # Accumulate the log likelihoods
  lls    <- numeric(max_its)
  lls[1] <- sum(dinvbatmix_pmat(x, pmat = pmat_cur, log = TRUE))


  if (verbose) cat("Starting log-likelihood: ", lls[1], "\n")

  for (i in 2:max_its) {

    if (verbose) cat("\n Iteration: ", i, ", component: ")

    # E-step
    W <- t(sapply(x, function(xi) {
      sapply(1:n_comp, function(k) {
        pmat_cur[k, 'alph'] * dinvbat(xi,
                                      pmat_cur[k, 'mu'],
                                      pmat_cur[k, 'kp'],
                                      pmat_cur[k, 'lam'])
      })
    }))
    W <- W / rowSums(W)

    # M-step

    # Update alphas, the component weights
    pmat_cur[, 'alph'] <- colSums(W) / n

    # Update the component parameters
    for (ci in 1:n_comp) {

      if (verbose) cat(" (", ci,") ")

      pmat_cur[ci, 1:3] <- maxlikinvbat(x,
                                        weights = pmat_cur[ci, 'alph'] * W[, ci],
                                        fixed_mu  = fixed_pmat[ci, 1],
                                        fixed_kp  = fixed_pmat[ci, 2],
                                        fixed_lam = fixed_pmat[ci, 3])
    }

    lls[i] <- sum(dinvbatmix_pmat(x, pmat = pmat_cur, log = TRUE))

    if (verbose) cat(", log-likelihood: ", lls[i])

    # Finish if the ll is not increasing anymore.
    if (abs(lls[i] - lls[i - 1]) < ll_tol) break
  }

  pmat_cur
}


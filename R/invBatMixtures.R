#' The Inverse Batschelet Mixture distribution
#'
#' @param x An angle in radians.
#' @param mus A vector of component mean directions.
#' @param kps A vector of nonnegative concentration parameters.
#' @param lams A vector of shape parameters.
#' @param alphs A vector of component weights.
#'
#' @return A vector of probabilities.
#' @export
#'
#' @examples
#' x <- rinvbat(10)
#' dinvbatmix(x)
#'
#' curve(dinvbatmix(x), -pi, pi)
#'
dinvbatmix <- function(x, mus = c(-pi/2, 0, pi/2), kps = c(8, 8, 8),
                       lams = c(-.5, 0, .5), alphs = c(.3, .4, .3), log = FALSE) {

  # Check the parameter setting
  if (!all.equal(length(mus), length(kps), length(lams), length(alphs))) {
    stop("Unequal parameter lengths")
  }

  if (log) {
    log(rowSums(sapply(1:length(mus), function(i) alphs[i] * dinvbat(x, mus[i], kps[i], lams[i]))))
  } else {
    rowSums(sapply(1:length(mus), function(i) alphs[i] * dinvbat(x, mus[i], kps[i], lams[i])))
  }
}




fitinvbatmix <- function(x, n_comp  = 4,
                         init_pmat  = matrix(NA, n_comp, 4),
                         fixed_pmat = matrix(NA, n_comp, 4),
                         n_its = 10,
                         verbose = FALSE) {

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
  lls    <- numeric(n_its)
  lls[1] <- sum(dinvbatmix(x, pmat_cur[, 'mu'], pmat_cur[, 'kp'], pmat_cur[, 'lam'], log = TRUE))

  if (verbose) cat("Starting log-likelihood: ", lls[1], "\n")

  for (i in 2:n_its) {

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

    lls[i] <- sum(dinvbatmix(x, pmat_cur[, 'mu'], pmat_cur[, 'kp'], pmat_cur[, 'lam'], log = TRUE))

    if (verbose) cat(", log-likelihood: ", lls[i])

  }

  pmat_cur
}


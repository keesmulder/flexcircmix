#' Fit mixtures of inverse or power Batschelet distributions.
#'
#' @param x A dataset of angles in radians.
#' @param bat_type String; Either "power" or "inverse", denoting the type of
#'   Batschelet distribution to employ.
#' @param n_comp The number of components to be used in the mixture. This is
#'   fixed, so it can not be estimated.
#' @param init_pmat An \code{n_comp * 4} matrix, the initial values of the
#'   parameter matrix. The parameters are ordered \code{mu}, \code{kp},
#'   \code{lam} and then \code{alph}.
#' @param fixed_pmat An \code{n_comp * 4} matrix, containing a parameter matrix,
#'   with \code{NA} for parameters to be estimated, and a numeric for each
#'   parameter that should be kept fixed to a specific value.
#' @param ll_tol Numeric; the algorithm is stopped after the log-likelihood
#'   improves less than \code{ll_tol}.
#' @param max_its The maximum number of iterations.
#' @param verbose Logical; whether to print debug statements.
#' @param optimization_its Integer; The maximum number of iterations to perform
#'   in the M-step (Maximization) part of the algorithm.
#'
#' @return A parameter matrix of results.
#' @export
#'
#' @examples
#' dat <- rinvbatmix(100)
#' batmixEM(dat, n_comp = 3)
#'
batmixEM <- function(x,
                     bat_type = "inverse",
                     n_comp  = 4,
                     init_pmat  = matrix(NA, n_comp, 4),
                     fixed_pmat = matrix(NA, n_comp, 4),
                     ll_tol = .1,
                     max_its = 50,
                     verbose = FALSE,
                     optimization_its = 10) {

  if (bat_type == "inverse") {
    dbat_fun      <- dinvbat
    likfunbat_fun <- likfuninvbat
  } else if (bat_type == "power") {
    dbat_fun      <- dpowbat
    likfunbat_fun <- likfunpowbat
  } else {
    stop("Unknown Batschelet type.")
  }

  # Force x to be in range -pi, pi.
  x <- force_neg_pi_pi(x)

  n <- length(x)

  # Save matrices giving which elements of init_pmat and fixed_pmat were provided.
  na_fixedpmat <- is.na(fixed_pmat)
  na_initpmat  <- is.na(init_pmat)

  # Set initial values if the initial parameter matrix is not given for that parameter (has NAs).
  if (any(na_initpmat[, 1])) init_pmat[, 1] <- seq(0, 2*pi, length.out = n_comp + 1)[-1]
  if (any(na_initpmat[, 2])) init_pmat[, 2] <- rep(5, n_comp)
  if (any(na_initpmat[, 3])) init_pmat[, 3] <- rep(0, n_comp)
  if (any(na_initpmat[, 4])) init_pmat[, 4] <- rep(1/n_comp, n_comp)

  # The current parameter matrix.
  pmat_cur <- init_pmat
  colnames(pmat_cur) <- c("mu", "kp", "lam", "alph")

  if (verbose) cat("Starting log-likelihood: ")

  # initialize W matrix, an n*n_comp matrix that has the probability of each
  # datapoint for each of the components, given the current parameter set.
  W <- matrix(1/n_comp, nrow = n, ncol = n_comp)

  # Accumulate the log likelihoods
  lls    <- numeric(max_its + 1)
  lls[1] <- sum(dbatmix_pmat(x, dbat_fun = dbat_fun, pmat = pmat_cur, log = TRUE))


  if (verbose) cat(lls[1], "\n")

  for (i in 1:max_its) {

    if (verbose) cat(" Iteration: ", sprintf("%4s", i), ", E-step: ", sep = "")

    # E-step
    W <- sapply(1:n_comp, function(k) {
      pmat_cur[k, 'alph'] * dbat_fun(x, pmat_cur[k, 'mu'], pmat_cur[k, 'kp'], pmat_cur[k, 'lam'])
    })


    W <- W / rowSums(W)


    # Update alphas, the component weights, if they are not fixed.
    for (k in 1:n_comp) {
      if (na_fixedpmat[k, 4]) {
        pmat_cur[, 'alph'] <- colSums(W) / n
      }
    }

    if (verbose) cat(" Done. ---  M-step: component: ")

    # M-step

    # Update the component parameters
    for (ci in 1:n_comp) {

      if (verbose) cat(" (", ci,")", sep = "")

      pmat_cur[ci, 1:3] <- maxlikbat(x, likfunbat_fun = likfunbat_fun,
                                     init_kp = pmat_cur[ci, 2], init_lam = pmat_cur[ci, 3],
                                     weights = pmat_cur[ci, 'alph'] * W[, ci],
                                     fixed_mu  = fixed_pmat[ci, 1],
                                     fixed_kp  = fixed_pmat[ci, 2],
                                     fixed_lam = fixed_pmat[ci, 3],
                                     max_its = optimization_its)
    }

    lls[i + 1] <- sum(dbatmix_pmat(x, dbat_fun = dbat_fun, pmat = pmat_cur, log = TRUE))

    if (verbose) cat(", done. Log-likelihood: ", lls[i + 1], ".\n", sep = "")

    # Finish if the ll is not increasing anymore.
    if (abs(lls[i + 1] - lls[i]) < ll_tol) break
  }

  pmat_cur
}


vectorize_pmat <- function(pmat) {
  vec        <- as.vector(pmat)
  names(vec) <- paste(rep(colnames(pmat), each = nrow(pmat)), 1:nrow(pmat), sep = "_")
  vec
}

# Function to compute the circular variance and circular sd and add it to a
# parameter matrix.
add_circ_var_to_pmat <- function(pmat, bat_type = "power") {
  var_mat <- t(apply(pmat, 1, function(row) {
    R_bar <- computeMeanResultantLengthBat(row["kp"], row["lam"], bat_type = bat_type)
    c(circ_var = 1 - R_bar, circ_sd = computeCircSD(R_bar))
  }))

  cbind(pmat, var_mat)
}


#' Fit a mixture of Batschelet distributions
#'
#' This is the main function of the package \code{flexcircmix}, and functions as
#' an interface to fit mixtures of Batschelet-type distributions, using
#' frequentist or Bayesian methods.
#'
#' @param x A dataset of angles in radians.
#' @param method Character; One of "bayes", "EM", or "boot". The method of obtaining a fit.
#' @param bat_type Character; Either "power" or "inverse", denoting the type of
#'   Batschelet distribution to employ.
#' @param n_comp The number of components to be used in the mixture. This is
#'   fixed, so it can not be estimated.
#' @param init_pmat An \code{n_comp * 4} matrix, the initial values of the
#'   parameter matrix. The parameters are ordered \code{mu}, \code{kp},
#'   \code{lam} and then \code{alph}.
#' @param fixed_pmat An \code{n_comp * 4} matrix, containing a parameter matrix,
#'   with \code{NA} for parameters to be estimated, and a numeric for each
#'   parameter that should be kept fixed to a specific value.
#' @param verbose Logical; whether to print debug statements.'
#' @param ... Additional arguments to be passed to the selected \code{method}.
#'
#' @return An object of class 'batmixmod'.
#' @export
#'
#' @examples
#'
#'
fitbatmix <- function(x,
                      method = "bayes",
                      ...) {

  bm_fit <- list(method = method)

  if (method == "bayes") {

    bm_fit$mcmc_sample <- mcmcBatscheletMixture(x, ...)

    # Placeholder, this is obviously not a good idea
    bm_fit$estimates <- colMeans(bm_fit$mcmc_sample)


    bm_fit$est_vector <- vectorize_pmat(bm_fit$estimates)

  } else if (method == "EM") {

    bm_fit$estimates  <- batmixEM(x, ...)
    bm_fit$est_vector <- vectorize_pmat(bm_fit$estimates)

  } else if (method == "boot") {

    bm_fit$estimates <- batmixEM(x, ...)

  } else stop("Method not found.")


  class(bm_fit) <- c("batmixmod", class(bm_fit))

  bm_fit

}





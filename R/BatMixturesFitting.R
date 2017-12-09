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

  if (any(!na_fixedpmat[, 2] & fixed_pmat[,2] < 0))      stop("Invalid fixed kappa value.")
  if (any(!na_fixedpmat[, 3] & abs(fixed_pmat[,3]) > 1)) stop("Invalid fixed lambda value.")

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


# FUNCTION meanDir  ---------------------------------------------------
# Calculates the mean direction from a dataset.
#   th:       A numeric vector containing the angles in the sample, in radians.
#   na.rm:    Whether NA's will be removed.
# Returns:    A scalar, the mean direction in radians, in the range [0, 2pi].
meanDir <- function (th, na.rm = TRUE) {
  C <- sum(cos(th), na.rm = na.rm)
  S <- sum(sin(th), na.rm = na.rm)
  force_neg_pi_pi(atan2(S, C))
}

###
### FUNCTION circularQuantile
###
# A wrapper for quantile() that first reverse-centres the circular data in order
# to prevent the results being influenced by an arbitrary starting point.
#   th:         A numeric vector containing the angles in the sample.
#   ...:        Arguments to be passed to quantile().
# Returns:    A vector containing the angles of the desired quantiles.
circularQuantile <- function (th, ...) {

  # The desired rotation before taking the quantile.
  rotation <- meanDir(th)

  # Centre the data, and move it as far away from 0 radians as possible by
  # th+rot. Then, apply the quantile function, and rotate back.
  force_neg_pi_pi(quantile(force_neg_pi_pi(th - rotation), ...) + rotation)
}

# Create a vector from a parameter matrix for convenience.
vectorize_pmat <- function(pmat) {
  vec        <- as.vector(pmat)
  names(vec) <- paste(rep(colnames(pmat), each = nrow(pmat)), 1:nrow(pmat), sep = "_")
  vec
}

# Create a matrix from a parameter vector for convenience.
matrixize_pvec <- function(pvec) {
  nms    <- names(pvec)
  n_comp <- sum(grepl("alph", nms))
  mat        <- matrix(pvec, nrow = n_comp)
  unique_nms <- nms[3 * (1:(length(nms)/3) -1) + 1]
  colnms     <- substr(unique_nms, 1, nchar(unique_nms) - 2)
  colnames(mat) <- colnms
  mat
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

summarize_one_mu_vector <- function(mu_vec, probs = c(.025, .975)) {
  R_bar <- sqrt(sum(cos(mu_vec))^2 + sum(sin(mu_vec))^2) / length(mu_vec)

  c(mean_dir = meanDir(mu_vec),
    circ_median = circularQuantile(mu_vec, .5),
    circ_se = computeCircSD(R_bar),
    circularQuantile(mu_vec, probs = probs))
}

summarize_one_lin_param <- function(pm_vec, probs = c(.025, .975)) {
  c(mean = mean(pm_vec), median = median(pm_vec), se = sd(pm_vec),
    quantile(pm_vec, probs = probs))
}

# Function to take a sample of parameters and compute a summary of it. Can be an
# mcmc sample or a bootstrap sample, for example.
summarize_batmix_param_sample <- function(bm_sam, probs = c(.025, .975)) {

  # First, treat circular variables
  mu_cols <- grepl("mu_[0-9]",   colnames(bm_sam))
  mu_mat  <- bm_sam[, mu_cols, drop = FALSE]
  mu_summary <- t(apply(mu_mat, 2, summarize_one_mu_vector, probs = probs))

  # Linear parameters
  lp_cols <- !mu_cols
  lp_mat  <- bm_sam[, lp_cols, drop = FALSE]
  lp_summary <- t(apply(lp_mat, 2, summarize_one_lin_param, probs = probs))

  full_summary <- rbind(mu_summary, lp_summary)
  colnames(full_summary) <- colnames(lp_summary)

  full_summary
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
                      bat_type = "power",
                      ...) {

  # Construct fit object.
  bm_fit <- list(method = method, bat_type = bat_type, x = x, call = match.call())

  if (method == "bayes") {

    bm_fit$mcmc_sample <- mcmcBatscheletMixture(x, bat_type = bat_type, ...)

    mcmc_sum <- summarize_batmix_param_sample(bm_fit$mcmc_sample)
    mcmc_sum <- mcmc_sum[!grepl("mean_res_len", rownames(mcmc_sum)), ]

    bm_fit$est_vector <- mcmc_sum[, 2]
    bm_fit$estimates <- matrixize_pvec(bm_fit$est_vector)
    bm_fit$mcmc_summary <- mcmc_sum

  } else if (method == "EM") {

    bm_fit$estimates   <- batmixEM(x, bat_type = bat_type, ...)

    # Augment the results
    bm_fit$estimates   <- add_circ_var_to_pmat(bm_fit$estimates, bat_type = bat_type)
    bm_fit$est_vector <- vectorize_pmat(bm_fit$estimates)

  } else if (method == "boot") {

    bm_fit <- c(bm_fit, bootstrapEMBatMix(x, bat_type = bat_type, ...))
    bm_fit$boot_summary <- summarize_batmix_param_sample(bm_fit$boot_sample)

  } else stop("Method not found.")

  rownames(bm_fit$estimates) <- paste("comp", 1:nrow(bm_fit$estimates), sep = "_")

  class(bm_fit) <- c("batmixmod", class(bm_fit))

  bm_fit

}





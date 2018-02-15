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
                     bat_type = "power",
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

    w_rowsum <- rowSums(W)

    # If some datapoints are numerically equal to zero, assign it equally to
    # each component. Hopefully, this does not happen again in the next
    # iteration.
    W[w_rowsum == 0, ] <- 1/n_comp
    w_rowsum[w_rowsum == 0] <- 1

    W <- W / w_rowsum

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
  unique_nms <- nms[n_comp * (1:(length(nms)/n_comp) - 1) + 1]
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
    circ_median = circularQuantile(mu_vec, .5, na.rm = TRUE),
    circ_se = computeCircSD(R_bar),
    circularQuantile(mu_vec, probs = probs))
}

summarize_one_lin_param <- function(pm_vec, probs = c(.025, .975)) {
  c(mean = mean(pm_vec, na.rm = TRUE),
    median = median(pm_vec, na.rm = TRUE),
    se = sd(pm_vec, na.rm = TRUE),
    quantile(pm_vec, probs = probs, na.rm = TRUE))
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


#' Printing function for \code{batmixmod} objects.
#'
#' @param bm_mod A \code{batmixmod} object.
#' @param ... Further arguments to be passed to print methods.
#'
#' @export
#'
print.batmixmod <- function(bm_mod, ...) {
  cat("Mixture of ", bm_mod$bat_type,
      " Batschelet distributions, using method '",
      bm_mod$method, "'.\n", sep = "")

  print(bm_mod$estimates, ...)
}

summary.batmixmod <- function(bm_mod) {
  if (bm_mod$method == "EM") {
    return(bm_mod$estimates)
  } else if (bm_mod$method == "bayes") {
    return(bm_mod$mcmc_summary)
  } else if (bm_mod$method == "boot") {
    return(bm_mod$boot_summary)
  } else {
    stop("Method not found." )
  }
}

# This function reorganizes a batmixmod object into a summary by component.
bm_summary_by_component <- function(bmm, add_ci = TRUE,
                                    add_circ_sd = TRUE,
                                    add_circ_var = FALSE) {

  # The reordering vector to get by-component results.
  reorder_vector <- as.vector(sapply(1:nrow(bmm$estimates),
                                     function(comp) grep(paste0("_", comp), names(bmm$est_vector))))
  nms <- names(bmm$est_vector[reorder_vector])

  # Gather results
  if (bmm$method == "EM") {
    out <- matrix(bmm$est_vector[reorder_vector])
    rownames(out) <- nms
  } else if(bmm$method == "boot") {
    out <- matrix(bmm$est_vector[reorder_vector])
    rownames(out) <- nms
    if (add_ci) {
      nsumcol <- ncol(bmm$boot_summary)
      out <- cbind(out, bmm$boot_summary[reorder_vector, (nsumcol-1):nsumcol] )
    }
  } else if(bmm$method == "bayes") {
    out <- matrix(bmm$est_vector[reorder_vector])
    rownames(out) <- nms
    if (add_ci) {
      nsumcol <- ncol(bmm$mcmc_summary)
      out <- cbind(out, bmm$mcmc_summary[reorder_vector, (nsumcol-1):nsumcol] )
    }
  }

  # Remove unwanted variances.
  if (!add_circ_sd)  out <- out[!grepl("circ_sd", nms), , drop = FALSE]
  if (!add_circ_var) out <- out[!grepl("circ_var", nms), , drop = FALSE]

  out
}

#' This function takes a list of bat_mix_mods, and provides a table that
#' compares the fits.
#'
#' @param bm_mod_list A list of \code{batmixmod} objects.
#' @param add_ci Logical; Whether to add confidence intervals.
#' @param add_circ_sd Logical; Whether to add circular sd.
#' @param add_circ_var Logical; Whether to add circular variance.
#'
#' @return A data frame with a summary.
#' @export
#'
multisummary.batmixmod <- function(bm_mod_list, add_ci = TRUE,
                                   add_circ_sd = TRUE, add_circ_var = FALSE) {


  # Summarize each model in the list.
  mod_summaries <- sapply(bm_mod_list, bm_summary_by_component,
         add_ci = add_ci, add_circ_sd = add_circ_sd, add_circ_var = add_circ_var)

  # Add NA to the missing variables.
  max_nrow <- max(sapply(mod_summaries, nrow))
  mod_summary_df <- as.data.frame(sapply(mod_summaries, function(x) {
    if (nrow(x) != max_nrow) {
      return(rbind(x, matrix(NA, ncol = ncol(x), nrow = max_nrow - nrow(x))))
    } else {
      return(x)
    }
  }))


  mod_summary_df
}



#' Fit a mixture of Batschelet distributions
#'
#' This is the main function of the package \code{flexcircmix}, and functions as
#' an interface to fit mixtures of Batschelet-type distributions, using
#' frequentist or Bayesian methods.
#'
#' @param x A dataset of angles in radians.
#' @param method Character; One of "bayes", "EM", or "boot". The method of
#'   obtaining a fit.
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
#' @param probs Numeric vector; The probabilities for which to compute quantiles
#'   in summarizing bootstrap or MCMC samples. By default, \code{probs = c(.025,
#'   .975)}, which corresponds to standard 95\% confidence or credible
#'   intervals.
#' @param post_est_median Logical; Only relevant for MCMC. Whether to use the
#'   posterior median as the estimate. If FALSE (the default) we use the mean,
#'   or mean direction for \code{mu}.
#' @param ... Additional arguments to be passed to the selected \code{method}.
#'   In particular, use \code{verbose = TRUE} to print debug statements.
#'
#' @return An object of class 'batmixmod'.
#' @export
#'
#' @examples
#' x <- rinvbatmix(50)
#' fitbatmix(x, method = "EM")
#'
fitbatmix <- function(x,
                      method = "bayes",
                      bat_type = "power",
                      n_comp = 4,
                      init_pmat  = matrix(NA, n_comp, 4),
                      fixed_pmat = matrix(NA, n_comp, 4),
                      probs = c(.025, .975),
                      post_est_median = TRUE,
                      ...) {


  # Check missings
  if (any(is.na(x))) {
    warning("Removing missing values from x.")
    x <- as.numeric(na.omit(x))
  }

  # Construct fit object.
  bm_fit <- list(method = method, bat_type = bat_type, x = x,
                 ic = list(),
                 call = match.call())

  if (method == "bayes") {

    mcmc_result <-  mcmcBatscheletMixture(x, bat_type = bat_type,
                                          n_comp = n_comp,
                                          init_pmat = init_pmat,
                                          fixed_pmat = fixed_pmat,
                                          ...)

    bm_fit$mcmc_sample      <- mcmc_result$mcmc_sample
    bm_fit$acceptance_rates <- mcmc_result$acceptance_rates
    bm_fit$ic               <- mcmc_result$ic
    bm_fit$log_posterior    <- mcmc_result$log_posterior


    mcmc_sum <- summarize_batmix_param_sample(bm_fit$mcmc_sample, probs = c(.025, .975))
    mcmc_sum <- mcmc_sum[!grepl("mean_res_len", rownames(mcmc_sum)), ]

    # If post_est_median == TRUE, we'll use the second column. Otherwise, the first.
    bm_fit$est_vector   <- mcmc_sum[, 1 + post_est_median]
    bm_fit$estimates    <- matrixize_pvec(bm_fit$est_vector)
    bm_fit$mcmc_summary <- mcmc_sum

  } else if (method == "EM") {

    bm_fit$estimates  <- batmixEM(x, bat_type = bat_type,
                                  n_comp = n_comp,
                                  init_pmat = init_pmat,
                                  fixed_pmat = fixed_pmat,
                                  ...)

    # Augment the results
    bm_fit$estimates  <- add_circ_var_to_pmat(bm_fit$estimates, bat_type = bat_type)
    bm_fit$est_vector <- vectorize_pmat(bm_fit$estimates)

  } else if (method == "boot") {

    bm_fit <- c(bm_fit, bootstrapEMBatMix(x, bat_type = bat_type,
                                          n_comp = n_comp,
                                          init_pmat = init_pmat,
                                          fixed_pmat = fixed_pmat,
                                          ...))

    bm_fit$boot_summary <- summarize_batmix_param_sample(bm_fit$boot_sample,
                                                         probs = probs)

  } else stop("Method not found.")


  # Miscellaneous outputs.
  bm_fit$n_components <- nrow(bm_fit$estimates)
  bm_fit$n_parameters <- sum(is.na(fixed_pmat))

  # Fix estimate names.
  rownames(bm_fit$estimates) <- paste("comp", 1:bm_fit$n_components, sep = "_")

  # Add the arguments of the call to the output.
  bm_fit$args <- list(...)

  if (bat_type == "power") {
    ll <- sum(dbatmix_pmat(x, dbat_fun = dpowbat,
                           pmat = bm_fit$estimates, log = TRUE))
  } else if (bat_type == "inverse") {
    ll <- sum(dbatmix_pmat(x, dbat_fun = dinvbat,
                           pmat = bm_fit$estimates, log = TRUE))
  } else {stop("Batschelet type should be 'inverse' or 'power'.")}

  # Parameter-only versions for convenience.
  bm_fit$param_estimates <- bm_fit$estimates[, 1:4]
  bm_fit$param_vector    <- bm_fit$est_vector[1:(bm_fit$n_components * 4)]

  # INFORMATION CRITERIA
  bm_fit$loglik <- ll

  bm_fit$ic <- c(list(loglik = ll,
                      deviance = -2 * ll,
                      n_param = bm_fit$n_parameters,
                      aic = c(p_aic = 2 *  bm_fit$n_parameters,
                              aic = 2 * (bm_fit$n_parameters - ll)),
                      bic = c(
                        p_bic = log(length(x)) * bm_fit$n_parameters,
                        bic = log(length(x)) * bm_fit$n_parameters - 2 * ll)),
                 bm_fit$ic)

  if (method == "bayes") {



    # Log-likelihood at Bayesian estimates.
    D_of_param_bar <- sum(dbatmix_pmat(x,
                                       dbat_fun = ifelse(bat_type == "power",
                                                         dpowbat, dinvbat),
                                       pmat = bm_fit$estimates,
                                       log = TRUE))
    D_bar <- mean(mcmc_result$ll_vec)
    p_d1 <- 2 * (D_of_param_bar - D_bar)
    p_d2 <- 2 * var(mcmc_result$ll_vec)

    bm_fit$ic$dic_1 <- c(p_dic1 = 2 * p_d1, dic1 = -2 * (D_of_param_bar - p_d1))
    bm_fit$ic$dic_2 <- c(p_dic2 = 2 * p_d2, dic2 = -2 * (D_of_param_bar - p_d2))
  }

  # Collect the ICs in a matrix.
  ic_locs       <- vapply(bm_fit$ic, length, 0) > 1
  ic_names      <- names(bm_fit$ic[ic_locs])
  bm_fit$ic_mat <- matrix(unlist(bm_fit$ic[ic_locs]), ncol = 2,
                          dimnames = list(ic_names, c("Penalty", "Value")),
                          byrow = TRUE)



  class(bm_fit) <- c("batmixmod", class(bm_fit))

  bm_fit

}





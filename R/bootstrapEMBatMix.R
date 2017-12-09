one_bat_boot_EM <- function(boot_x, bat_type = "power", n_comp  = 4,
                            init_pmat  = matrix(NA, n_comp, 4), fixed_pmat = matrix(NA, n_comp, 4),
                            ll_tol = .1, max_its = 50, optimization_its = 10) {

  pmat <- batmixEM(x, bat_type = bat_type, n_comp  = n_comp, verbose = FALSE,
                   init_pmat  = init_pmat, fixed_pmat = fixed_pmat,
                   ll_tol = ll_tol, max_its = max_its, optimization_its = optimization_its)

  # Augment output, and return as a vector to fit easily in the sample matrix.
  full_pmat <- add_circ_var_to_pmat(pmat)
  vectorize_pmat(full_pmat)
}



#' Bootstrap the Batschelet mixture parameters
#'
#' @param x Numeric; A set angles in radians.
#' @param B Integer; The number of bootstrap samples.
#' @param parallel Logical; Whether to perform bootstraps in parallel.
#' @param n_comp Integer; Fixed number of components to estimate.
#' @param bat_type Either 'inverse' or 'power', the type of distribution to fit.
#'   The two distributions are similar, but the power Batschelet distribution is
#'   computationally much less demanding.
#' @param init_pmat A numeric matrix with \code{n_comp} rows and four columns,
#'   corresponding to \eqn{\mu, \kappa, \lambda, \alpha}, in that order. Gives
#'   starting values for the parameters. If any element is \code{NA}, it will be
#'   given a default starting value. For \eqn{mu}, the default starting values
#'   are equally spaced on the circle. For \eqn{\kappa}, the default starting
#'   value is 5. For \eqn{\lambda}, the default starting value is 0, which
#'   corresponds to the von Mises distribution. For \eqn{\alpha}, the default
#'   starting value is \code{1/n_comp}.
#' @param fixed_pmat A numeric matrix with \code{n_comp} rows and four columns,
#'   corresponding to \eqn{\mu, \kappa, \lambda, \alpha}, in that order. Any
#'   element that is not \code{NA} in this matrix will be held constant at the
#'   given value and not sampled.
#' @param verbose Logical; Whether to print debug info.
#' @param ll_tol The likelihood tolerance for the EM-algorithm.
#' @param max_its The maximum number of E-M iterations for the EM-algorithm.
#' @param optimization_its The maximum number of maximization iterations within the EM-algorithm.
#'
#' @return A matrix of sampled bootstrap values.
#' @export
#'
#' @examples
#' x <- rinvbatmix(100)
#'
#'  bootstrapEMBatMix(x, B = 5)
#'
bootstrapEMBatMix <- function(x, B = 500, parallel = TRUE, verbose = FALSE,
                              bat_type = "power", n_comp  = 4,
                              init_pmat  = matrix(NA, n_comp, 4), fixed_pmat = matrix(NA, n_comp, 4),
                              ll_tol = .1, max_its = 50, optimization_its = 10) {

  # Original data fit.
  fulldata_fit_pmat <- batmixEM(x, bat_type = bat_type, n_comp  = n_comp, verbose = verbose,
                                init_pmat  = init_pmat, fixed_pmat = fixed_pmat,
                                ll_tol = ll_tol, max_its = max_its, optimization_its = optimization_its)


  if (verbose) cat("\nStarting bootstrap:\n")

  est_vector <- vectorize_pmat(fulldata_fit_pmat)
  out <- list(estimates = fulldata_fit_pmat, est_vector = est_vector)

  if (parallel) {

    # Set up cluster
    no_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(no_cores)
    parallel::clusterEvalQ(cl, library(flexcircmix))
    parallel::clusterExport(cl, c("one_bat_boot_EM",
                                  "x", "bat_type", "n_comp", "fulldata_fit_pmat", "fixed_pmat",
                                  "ll_tol", "max_its", "optimization_its"))

    # Run the bootstrap
    rep_sam <- pbreplicate(B, {
      one_bat_boot_EM(boot_x = sample(x, size = length(x), replace = TRUE),
                      bat_type = bat_type, n_comp  = n_comp,
                      init_pmat  = fulldata_fit_pmat, fixed_pmat = fixed_pmat,
                      ll_tol = ll_tol, max_its = max_its, optimization_its = optimization_its)
    }, cl = cl)

    parallel::stopCluster(cl)
  } else {
    # Run the bootstrap non-parallelized
    rep_sam <- pbreplicate(B, {
      one_bat_boot_EM(boot_x = sample(x, size = length(x), replace = TRUE),
                      bat_type = bat_type, n_comp  = n_comp,
                      init_pmat  = fulldata_fit_pmat, fixed_pmat = fixed_pmat,
                      ll_tol = ll_tol, max_its = max_its, optimization_its = optimization_its)
    })

    parallel::stopCluster(cl)

  }

  if (verbose) cat("\nFinished.\n")

  out$boot_sample <- t(rep_sam)

  out
}


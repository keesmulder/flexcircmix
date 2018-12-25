# Helper function to represent circular variables (such as mean directions) as
# "gapless" numerical representations.
gapless_circ <- function(th) {
  md     <- meanDir(th)
  new_th <- ((th - md + pi) %% (2*pi)) - pi + md
  new_th
}



#' Log Marginal Likelihood via Bridge Sampling
#'
#' Obtain the log of the Marginal Likelihood through Bridge sampling. This is
#' useful to perform any Bayesian hypothesis test between any set of models.
#'
#' More details are provided under
#' \code{\link[bridgesampling:bridge_sampler]{bridge_sampler}}.
#'
#' @param samples A \code{batmixmod} object (from \code{\link{fitbatmix}}),
#'   which contains the MCMC sample required for bridge sampling.
#' @param ... Additional arguments to be passed to \code{bridge_sampler.matrix}.
#'
#' @seealso \code{ \link[brms:bayes_factor]{bayes_factor},
#'   \link[brms:post_prob]{post_prob} }
#'
#' @method bridge_sampler batmixmod
#' @importFrom bridgesampling bridge_sampler
#' @export bridge_sampler
#' @export
#'
bridge_sampler.batmixmod <- function(samples, ...) {

  if (samples$method != "bayes") {
    stop("Bridge sampling is only possible for method =='bayes'.")
  }

  n_comp <- samples$n_components

  sam <- samples$mcmc_sample[, 1:(4*n_comp)]

  lb <- rep(c(-2*pi, 0, -1, 0), each = n_comp)
  ub <- rep(c(2*pi, Inf, 1, 1), each = n_comp)

  names(lb) <- colnames(sam)
  names(ub) <- colnames(sam)

  bridgesampling::bridge_sampler(data = samples$x,
                                 samples = as.matrix(sam),
                                 param_types = rep(c("circular", "real", "real", "simplex"), each = n_comp),
                                 log_posterior = samples$log_posterior,
                                 lb = lb, ub = ub, ...)
}


#' Extract AIC from Batmixmod
#'
#' @param object A \code{batmixmod} object.
#' @param ... Ignored.
#'
#' @return The numerical value of the AIC.
#'
#' @method AIC batmixmod
#' @importFrom stats AIC
#' @export AIC
#' @export
#'
AIC.batmixmod <- function(object, ...) {
  object$ic_mat["aic", 2]
}



#' Extract BIC from Batmixmod
#'
#' @param object A \code{batmixmod} object.
#' @param ... Ignored.
#'
#' @return The numerical value of the BIC.
#' @method BIC batmixmod
#' @importFrom stats BIC
#' @export BIC
#' @export
#'
BIC.batmixmod <- function(object, ...) {
  object$ic_mat["bic", 2]
}
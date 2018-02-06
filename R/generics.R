#' Log Marginal Likelihood via Bridge Sampling
#'
#'
#'   More details are provided under
#'   \code{\link[bridgesampling:bridge_sampler]{bridge_sampler}}.
#'
#' @param bm_mod A \code{batmixmod} object.
#' @param ... Additional arguments to be passed to \code{bridge_sampler.matrix}
#'
#' @seealso \code{
#'   \link[brms:bayes_factor]{bayes_factor},
#'   \link[brms:post_prob]{post_prob}
#' }
#'
#'
#' @method bridge_sampler brmsfit
#' @importFrom bridgesampling bridge_sampler
#' @export bridge_sampler
#' @export
#'
bridge_sampler.batmixmod <- function(bm_mod, ...) {

  if (bm_mod$method != "bayes") {
    stop("Bridge sampling is only possible for method =='bayes'.")
  }

  sam <- bm_mod$mcmc_sample[, 1:bm_mod$n_parameters]

  lb <- rep(c(-pi, 0, -1, 0), each = bm_mod$n_components)
  ub <- rep(c(pi, Inf, 1, 1), each = bm_mod$n_components)

  names(lb) <- colnames(sam)
  names(ub) <- colnames(sam)

  bridgesampling::bridge_sampler(data = x,
                                 samples = as.matrix(sam),
                                 log_posterior = bm_mod$log_posterior,
                                 lb = lb, ub = ub, ...)
}
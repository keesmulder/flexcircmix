# Helper function to represent circular variables (such as mean directions) as
# "gapless" numerical representations.
gapless_circ <- function(th) {
  md     <- flexcircmix:::meanDir(th)
  new_th <- ((th - md + pi) %% (2*pi)) - pi + md
  new_th
}


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
#' @method bridge_sampler batmixmod
#' @importFrom bridgesampling bridge_sampler
#' @export bridge_sampler
#' @export
#'
bridge_sampler.batmixmod <- function(bm_mod, ...) {

  if (bm_mod$method != "bayes") {
    stop("Bridge sampling is only possible for method =='bayes'.")
  }

  n_comp <- bm_mod$n_components

  sam <- bm_mod$mcmc_sample[, 1:(4*n_comp)]

  # Force representation of the posterior samples as "gapless" numbers. This is
  # important because the bridge sampler will approximate the posterior with a
  # normal distribution based on the posterior covariance matrix, will will be
  # too wide if there are gaps, for example a posterior represented in the range
  # 6 - 2*pi as well as 0 - .2.
  sam[, 1:n_comp] <- apply(sam[, 1:n_comp], 2, gapless_circ)

  lb <- rep(c(-pi, 0, -1, 0), each = n_comp)
  ub <- rep(c(pi, Inf, 1, 1), each = n_comp)

  names(lb) <- colnames(sam)
  names(ub) <- colnames(sam)

  bridgesampling::bridge_sampler(data = bm_mod$x,
                                 samples = as.matrix(sam),
                                 log_posterior = bm_mod$log_posterior,
                                 lb = lb, ub = ub, ...)
}
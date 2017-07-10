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
                         fixed_mus  = rep(NA, n_comp),
                         fixed_kps  = rep(NA, n_comp),
                         fixed_lams = rep(NA, n_comp),
                         n_its = 10) {

  n <- length(x)

  cur_mus  <- seq(0, 2*pi, length.out = n_comp + 1)[-1]
  cur_kps  <- rep(5, n_comp)
  cur_lams <- rep(0, n_comp)
  cur_alphs <- rep(1/n_comp, n_comp)

  # initialize W matrix
  W <- matrix(1/n_comp, nrow = n, ncol = n_comp)

  ll_cur <- sum(dinvbatmix(x, cur_mus, cur_kps, cur_lams, cur_alphs, log = TRUE))


  for (i in 1:n_its) {

    print(i)

    # E-step
    W <- t(sapply(x, function(xi) {
      sapply(1:n_comp, function(k) cur_alphs[k] * dinvbat(xi, cur_mus[k], cur_kps[k], cur_lams[k]))
      }))
    W <- W / rowSums(W)

    # M-step

    # Update alphas, the component weights
    cur_alphs <- colSums(W) / n

    # Update the component parameters
    for (ci in 1:n_comp) {

      params <- maxlikinvbat(x, W[, ci],
                             fixed_mu  = fixed_mus[ci],
                             fixed_kp  = fixed_kps[ci],
                             fixed_lam = fixed_lams[ci])

      cur_mus[ci] <- params["mu"]
      cur_kps[ci] <- params["kp"]
      cur_lams[ci] <- params["lam"]
    }

  }

  prs <- cbind(mu = cur_mus, kps = cur_kps, lams = cur_lams, alphs = cur_alphs)

}



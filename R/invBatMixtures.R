#' The Batschelet-type Mixture distribution
#'
#' @param x An angle in radians.
#' @param mus A vector of component mean directions.
#' @param kps A vector of nonnegative concentration parameters.
#' @param lams A vector of shape parameters.
#' @param alphs A vector of component weights.
#' @param bat_fun A function that provides the pdf function of the desired Batschelet density. By default, inverse Batschelet.
#' @param log Logical; indicates if the logarithm of the density must be return.
#'
#' @return A vector of probabilities.
#' @export
#'
#' @examples
#' x <- rinvbat(10)
#' dbatmix(x)
#'
#' curve(dbatmix(x), -pi, pi)
#'
dbatmix <- function(x, dbat_fun = dinvbat,
                       mus = c(-pi/2, 0, pi/2), kps = c(8, 8, 8),
                       lams = c(-.5, 0, .5), alphs = c(.3, .4, .3), log = FALSE) {

  # Check the parameter sizes.
  if (!all.equal(length(mus), length(kps), length(lams), length(alphs))) {
    stop("Unequal parameter lengths")
  }

  # Compute the probability per component for each datapoint.
  pkimat <- sapply(1:length(mus), function(i) alphs[i] * dbat_fun(x, mus[i], kps[i], lams[i]))

  # If length(x) == 1, a vector is return, which must be made into a matrix.
  if (!is.matrix(pkimat)) pkimat <- t(pkimat)

  if (log) {
      log(rowSums(pkimat))
  } else {
      rowSums(pkimat)
  }
}


#' @describeIn dbatmix A version that takes a parameter matrix as input.
dbatmix_pmat <- function(x, dbat_fun = dinvbat,
                            pmat = cbind(mu  = c(-pi/2, 0, pi/2), kp = c(8, 8, 8),
                                         lam = c(-.5, 0, .5),     alph = c(.3, .4, .3)),
                            log = FALSE) {

  # Compute the probability per component for each datapoint.
  pkimat <- sapply(1:nrow(pmat), function(i) pmat[i, 4] * dbat_fun(x, pmat[i, 1], pmat[i, 2], pmat[i, 3]))

  # If length(x) == 1, a vector is return, which must be made into a matrix.
  if (!is.matrix(pkimat)) pkimat <- t(pkimat)

  if (log) {
    log(rowSums(pkimat))
  } else {
    rowSums(pkimat)
  }
}




#' @describeIn dbatmix Random variate generation for Inverse Batschelet mixtures.
rinvbatmix <- function(n, mus = c(-pi/2, 0, pi/2), kps = c(8, 8, 8),
                       lams = c(-.5, 0, .5), alphs = c(.3, .4, .3)) {

  # Check the parameter sizes.
  if (!all.equal(length(mus), length(kps), length(lams), length(alphs))) {
    stop("Unequal parameter lengths")
  }

  n_comp <- length(mus)

  # Each value will be sampled from a specific component.
  chosen_comp <- sample(1:n_comp, size = n, prob = alphs, replace = TRUE)

  sam <- numeric(n)

  # Sample data for each component separately.
  for (i in 1:n_comp) {
    idx <- chosen_comp == i
    sam[idx] <- rinvbat(sum(idx), mus[i], kps[i], lams[i])
  }

  sam
}

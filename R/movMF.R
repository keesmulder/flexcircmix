
#' Compute the polar versions of the movMF results
#'
#' The results of a movMF object are given as a vector called theta. We move to polar coordinates with
#' mu and kp.
#'
#' @param m A movMF object.
#'
#' @return A matrix of mus and kps.
#' @export
#'
coef_polar_movMF <- function(m) {
  t1 <- coef(m)$theta[, 1]
  t2 <- coef(m)$theta[, 2]
  mu <- atan2(t2, t1)
  kp <- sqrt(t1^2 + t2^2)
  cbind(mu = mu, kp = kp)
}

#' Obtain the inverse Batschelet parameter matrix from a movMF object
#'
#' This function allows using the more general inverse Batschelet functions, such as plotting. The
#' results of a movMF object are given as a vector called theta. We move to polar coordinates with
#' mu and kp, set lambda to 0 (the von Mises distribution as a special case of the inverse
#' Batschelet), and obtain the alphas from the movMF object.
#'
#' @param m A movMF object.
#'
#' @return A parameter matrix with inverse Batschelet mixture parameters.
#' @export
#'
invbatmix_pmat_from_movMF <- function(m) {
  cbind(coef_polar_movMF(m), lam = 0, alph = m$alpha)
}


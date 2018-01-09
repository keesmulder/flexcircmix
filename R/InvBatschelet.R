#' Compute a partial normalizing constant of a vM-based inverse Batschelet distribution
#'
#' The value givin is only the additional part of the inverse Batschelet. The von Mises normalizing
#' constant must still be added.
#'
#' @inheritParams dinvbat
#'
#' @return Numeric.
#'
K_kplam <- function(kp, lam) {
  (1 + lam) / (1 - lam) - (2 * lam / (1 - lam)) *
    stats::integrate(function(x) dvmkern(x - 0.5 * (1 - lam) * sin(x), 0, kp), -pi, pi)$value /
    (2 * pi * besselI(nu = 0, x = kp))
}

#' Set an angle to have its bounds in (-pi, pi)
#'
#' This function is used to make sure that an angle is represented as being between -pi and pi.
#'
#' @param x Numeric, representing an angle.
#' @export
#'
#' @return  Numeric, representing an angle, between -pi and pi.
#'
force_neg_pi_pi <- function(x) {
  ((x + pi) %% (2*pi)) - pi
}

#' Inverse Batschelet subfunction
#'
#' @inheritParams t_lam
#'
#' @return Numeric.
#'
s_lam <- function(x, lam) {
  x - 0.5 * (1 + lam) * sin(x)
}


#' Inverse Batschelet inverse subfunction
#'
#' @describeIn s_lam
#'
s_lam_inv <- function(x, lam) {
  # Compute the root to obtain the inverse.
  uniroot(function(y) s_lam(y, lam) - x, lower = -pi, upper = pi)$root
}

#' Vectorized version of the inversion of s_lam
#'
#' @describeIn s_lam
#'
s_lam_inv_vec <-  Vectorize(s_lam_inv)



#' Inverse Batschelet function
#'
#' Function used to generate different shapes for the inverse Batschelet distribution.
#'
#' @param x An angle in radians.
#' @param lam The shape parameter.
#'
#' @return The transformed angle x.
#' @export
#'
#' @examples
#'
#' # Reproduce Figure 5 in Jones & Pewsey
#' curve(t_lam(x, 0), -pi, pi, asp = 1)
#' for (i in seq(-1, 1, by = 1/3)) {
#'   curve(t_lam(x, i), -pi, pi, add = TRUE, lty = "dotted")
#' }
t_lam <- Vectorize(function(x, lam) {

  # Make sure the the input x is in the correct bounds.
  if (x > pi || x < -pi) x <- force_neg_pi_pi(x)

  if (lam != -1) {
    x * (1 - lam) / (1 + lam) + s_lam_inv(x, lam) * 2 * lam / (1 + lam)
  } else {
    x - sin(x)
  }
})



#' Inverse Batschelet inverse function
#'
#' @describeIn t_lam
#'
t_lam_inv <- function(x, lam) t_lam(x, -lam)



#' Kernel of the von-Mises based symmetric inverse Batschelet distribution
#'
#' @inheritParams dinvbat
#'
#' @return The unnormalized density value of the inverse Batschelet distribution.
#'
dinvbatkern <- function(x, mu = 0, kp = 1, lam = 0, log = FALSE) {
  dvm(t_lam(x - mu, lam), mu = 0, kp = kp, log = log)
}



#' The von-Mises based symmetric inverse Batschelet distribution
#'
#' The inverse Batschelet distribution, with mean direction \code{mu}, concentration parameter \code{kp}, and
#' shape (peakedness) parameter \code{lam}. This is the von Mises based version, without a skewness
#' parameter.
#'
#' @param x An angle in radians.
#' @param mu A mean direction, in radians.
#' @param kp Numeric, \eqn{> 0,}the concentration parameter.
#' @param lam The shape parameter (peakedness), -1 < \code{lam} < 1.
#' @param log Logical; whether to return the log of the probability or not.
#'
#' @return Numeric, the probability or log-probability of angle x given the parameters.
#' @export
#'
#' @examples
#' dinvbat(3)
#'
#' # Peaked distribution
#' curve(dinvbat(x, lam = .8), -pi, pi)
#'
#' # Flat-topped distribution
#' curve(dinvbat(x, lam = -.8), -pi, pi)
#'
dinvbat <- function(x, mu = 0, kp = 1, lam = 0, log = FALSE) {
  if (kp < 0) return(NA)
  if (lam <= -1) return(NA)
  if (lam > 1) return(NA)

  if (log) {
    dinvbatkern(x, mu = mu, kp = kp, lam = lam, log = TRUE) - log(K_kplam(kp = kp, lam = lam))
  } else {
    dinvbatkern(x, mu = mu, kp = kp, lam = lam, log = FALSE) / K_kplam(kp = kp, lam = lam)
  }
}


#' Compute the weight function of the inverse Batschelet random generation function
#'
#' @inheritParams dinvbat
#'
#' @return Numeric.
#'
weight_fun_rinvbat <- function(x, lam) {

  wr <- ((1 - 0.5 * (1 + lam) * cos(s_lam_inv(x, -lam))) /
           (1 - 0.5 * (1 - lam) * cos(s_lam_inv(x, -lam))))

  if (lam > 0) {
    ((3 - lam) / (3 + lam)) * wr
  } else {
    ((1 + lam) / (1 - lam)) * wr
  }
}


#' @describeIn dinvbat Random generation
#' @export
rinvbat <- function(n, mu = 0, kp = 1, lam = 0) {

  if (lam > 1 || lam <= -1 || kp < 0) stop("Parameter out of bounds.")

  th_out <- numeric(n)

  for (i in 1:n) {

    # Rejection sampler
    accepted <- FALSE
    while (!accepted) {

      th_can <- force_neg_pi_pi(circglmbayes::rvmc(1, 0, kp))
      u <- runif(1, 0, 1)

      # The weight function alters the von Mises candidate to include peakedness set by the lambda
      w_lam <- weight_fun_rinvbat(th_can, lam)

      if (u <= w_lam) {
        accepted <- TRUE

        th_out[i] <- t_lam(th_can, -lam)
      }
    }

  }

  # Introduce the mean and set sample space to -pi, pi.
  force_neg_pi_pi(th_out + mu)
}




#' Obtain the likelihood of an inverse Batschelet distribution
#'
#' @param x An set of angles in radians.
#' @param weights A vector of length \code{length(x)}, which gives importance weights to be used for x.
#' @param log If \code{TRUE} (the default), the log-likelihood is used.
#' @param mu A mean direction, in radians.
#' @param kp Numeric, \eqn{> 0,}the concentration parameter.
#' @param lam The shape parameter (peakedness), -1 < \code{lam} < 1.
#'
#' @return \code{likinvbat} returns a value, the likelihood given the data and parameters.
#'   \code{likfuninvbat} returns a function of mu, kp and lam, which can be evaluated later for a
#'   given set of parameters.
#' @export
#'
#' @examples
#' x <- rinvbat(5)
#'
#' # Find the likelihood
#' likinvbat(x, mu = 0, kp = 1, lam = 0.1, log = TRUE)
#'
#' # likfuninvbat returns a function.
#' llfib <- likfuninvbat(x)
#' llfib(mu = 0, kp = 1, lam = 0.1)
#'
likfuninvbat <- function(x, weights = rep(1, length(x)), log = TRUE) {
  if (log) {
    function(mu, kp, lam) sum(weights * dinvbat(x, mu, kp, lam, log = TRUE))
  } else {
    function(mu, kp, lam) exp(sum(weights * dinvbat(x, mu, kp, lam, log = TRUE)))
  }
}

#' Likelihood function of inverse Batschelet.
#'
#'  @describeIn likfuninvbat
#'  @export
#'
likinvbat <- function(x, mu, kp, lam, weights = rep(1, length(x)), log = TRUE) {
  if (log) {
    sum(weights * dinvbat(x, mu, kp, lam, log = TRUE))
  } else {
    exp(sum(weights * dinvbat(x, mu, kp, lam, log = TRUE)))
  }
}




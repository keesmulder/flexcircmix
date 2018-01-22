
# The right hand side of the log likelihood of a Batschelet-type distribution
ll_rhs_bat <- function(x, mu, kp, lam, tlam_fun) {
  kp * sum(cos(tlam_fun(x - mu, lam)))
}


# # The left hand side of the log likelihood of a inverse Batschelet distribution
# ll_lhs_invbat <- function(n, kp, lam) {
#   -n * (logBesselI(kp, 0) + log(K_kplam(kp, lam)))
# }


sample_mu_bat <- function(x, mu_cur, kp, lam, tlam_fun, mu_logprior_fun) {

  # Sample a candidate from the distribution of mu when we have a von Mises
  # distribution.
  C_j    <- sum(cos(x))
  S_j    <- sum(sin(x))
  R_j    <- sqrt(C_j^2 + S_j^2)
  mu_can <- circglmbayes::rvmc(1, mu_cur, R_j * kp)

  ll_can <- ll_rhs_bat(x, mu_can, kp, lam, tlam_fun)
  ll_cur <- ll_rhs_bat(x, mu_cur, kp, lam, tlam_fun)

  # The proposal is von Mises and thus symmetric, so the transition
  # probabilities of MH are omitted here.
  mu_lograt <- ll_can + mu_logprior_fun(mu_can) - ll_cur - mu_logprior_fun(mu_cur)

  if (mu_lograt > log(stats::runif(1))) {
    return(mu_can)
  } else {
    return(mu_cur)
  }
}



sample_mu_bat_2 <- function(x, mu_cur, kp, lam, tlam_fun, mu_logprior_fun) {

  # Sample a candidate from the distribution of mu when we have a von Mises
  # distribution.
  C_j    <- sum(cos(x))
  S_j    <- sum(sin(x))
  R_j    <- sqrt(C_j^2 + S_j^2)
  mu_hat <- atan2(S_j, C_j)
  mu_can <- circglmbayes::rvmc(1, mu_hat, R_j * kp)

  ll_can <- ll_rhs_bat(x, mu_can, kp, lam, tlam_fun)
  ll_cur <- ll_rhs_bat(x, mu_cur, kp, lam, tlam_fun)

  logp_mu_can_to_cur <- dvm(mu_cur, mu_hat, kp, log = TRUE)
  logp_mu_cur_to_can <- dvm(mu_can, mu_hat, kp, log = TRUE)

  mu_lograt <- ll_can + mu_logprior_fun(mu_can) + logp_mu_can_to_cur -
    ll_cur - mu_logprior_fun(mu_cur) - logp_mu_cur_to_can

  if (mu_lograt > log(stats::runif(1))) {
    return(mu_can)
  } else {
    return(mu_cur)
  }
}

# Reparametrized gamma proposal to make tuning parameter and mean interpretable.
# This is equal to chi square with 'mean' degrees of freedom if var_tune = 1.
dgammaprop <- function(x, mean = 1, var_tune = 1, log = FALSE) {
  gamma_var   <- 2 * mean * var_tune
  gamma_scale <- gamma_var  / mean
  gamma_shape <- mean / gamma_scale

  stats::dgamma(x, shape = gamma_shape, scale = gamma_scale, log = log)
}
rgammaprop <- function(n = 1, mean = 1, var_tune = 1) {
  gamma_var   <- 2 * mean * var_tune
  gamma_scale <- gamma_var  / mean
  gamma_shape <- mean / gamma_scale

  stats::rgamma(n = n, shape = gamma_shape, scale = gamma_scale)
}



sample_kp_bat <- function(x, mu, kp_cur, lam, llbat, kp_logprior_fun, var_tune = 1) {

  # Sample a candidate
  kp_can <- rgammaprop(1, mean = kp_cur, var_tune = 1)

  ll_can <- llbat(x, mu, kp_can, lam, log = TRUE)
  ll_cur <- llbat(x, mu, kp_cur, lam, log = TRUE)

  logp_kp_can_to_cur <- dgammaprop(kp_cur, kp_can, var_tune, log = TRUE)
  logp_kp_cur_to_can <- dgammaprop(kp_can, kp_cur, var_tune, log = TRUE)

  kp_lograt <- ll_can + kp_logprior_fun(kp_can) + logp_kp_can_to_cur -
    ll_cur - kp_logprior_fun(kp_cur) - logp_kp_cur_to_can

  if (kp_lograt > log(stats::runif(1))) {
    return(kp_can)
  } else {
    return(kp_cur)
  }
}

sample_lam_bat <- function(x, mu, kp, lam_cur, llbat, lam_logprior_fun, lam_bw = .05) {

  # Sample a candidate
  lam_can <- stats::runif(1, max(-1, lam_cur - lam_bw), min(1, lam_cur + lam_bw))

  ll_can <- llbat(x, mu, kp, lam_can)
  ll_cur <- llbat(x, mu, kp, lam_cur)

  logp_lam_can_to_cur <- stats::dunif(lam_cur, max(-1, lam_can - lam_bw),
                                      min(1, lam_can + lam_bw), log = TRUE)
  logp_lam_cur_to_can <- stats::dunif(lam_can, max(-1, lam_cur - lam_bw),
                                      min(1, lam_cur + lam_bw), log = TRUE)

  lam_lograt <- ll_can + lam_logprior_fun(lam_can) + logp_lam_can_to_cur -
    ll_cur - lam_logprior_fun(lam_cur) - logp_lam_cur_to_can

  if (lam_lograt > log(stats::runif(1))) {
    return(lam_can)
  } else {
    return(lam_cur)
  }
}



sample_kp_and_lam_bat <- function(x, mu, kp_cur, lam_cur, llbat, lam_bw = .05,
                                  kp_logprior_fun, lam_logprior_fun, var_tune = 1) {

  # Sample a candidate
  kp_can  <- rgammaprop(1, mean = kp_cur, var_tune = var_tune)
  lam_can <- stats::runif(1,
                          max(-1, lam_cur - lam_bw),
                          min(1, lam_cur + lam_bw))

  logp_kp_can_to_cur <- dgammaprop(kp_cur, kp_can, var_tune, log = TRUE)
  logp_kp_cur_to_can <- dgammaprop(kp_can, kp_cur, var_tune, log = TRUE)
  logp_lam_can_to_cur <- stats::dunif(lam_cur,
                                      max(-1, lam_can - lam_bw),
                                      min(1, lam_can + lam_bw),
                                      log = TRUE)
  logp_lam_cur_to_can <- stats::dunif(lam_can,
                                      max(-1, lam_cur - lam_bw),
                                      min(1, lam_cur + lam_bw),
                                      log = TRUE)

  ll_can <- llbat(x, mu, kp_can, lam_can, log = TRUE)
  ll_cur <- llbat(x, mu, kp_cur, lam_cur, log = TRUE)

  kplam_lograt <- ll_can + kp_logprior_fun(kp_can) + lam_logprior_fun(lam_can) +
    logp_kp_can_to_cur + logp_lam_can_to_cur -
    ll_cur - kp_logprior_fun(kp_cur) - lam_logprior_fun(lam_cur) -
    logp_kp_cur_to_can - logp_lam_cur_to_can

  if (kplam_lograt > log(stats::runif(1))) {
    return(c(kp_can, lam_can))
  } else {
    return(c(kp_cur, lam_cur))
  }
}


#' A rescaled beta prior for lambda
#'
#' This prior is symmetric between -1 and 1, and captures the prior belief that
#' values of \code{lam} on the boundary of the parameter space are a prior
#' unlikely.
#'
#' @param lam Numeric;
#'
#' @return The log-prior probability.
#' @export
#'
lam_beta_log_prior_2_2 <- function(lam) {
  stats::dbeta( (lam + 1) / 2, 2, 2, log = TRUE)
}


#' MCMC sampling for Batschelet-type distributions.
#'
#' @param x A numeric vector of angles, in radians
#' @param Q Integer; The number of iterations to return after taking burn in and
#'   thinning into account.
#' @param burnin Integer; The number of (non-thinned) iterations to discard. No
#'   burn in is performed by default.
#' @param thin Integer; Number of iterations to sample for each saved iteration.
#'   Defaults to 1, which means no thinning.
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
#' @param mu_logprior_fun Function; A function with a single argument, which
#'   returns the log of the prior probability of \eqn{\mu}. Defaults to a
#'   uniform prior function.
#' @param kp_logprior_fun Function; A function with a single argument, which
#'   returns the log of the prior probability of \eqn{\kappa}. Defaults to a
#'   uniform prior function. In contrast to the other parameters, for
#'   \eqn{\kappa} the constant (uniform) prior is improper.
#' @param lam_logprior_fun Function; A function with a single argument, which
#'   returns the log of the prior probability of \eqn{\lambda}. Defaults to a
#'   uniform prior function.
#' @param alph_prior_param Integer vector; The mixture weight parameter vector
#'   \eqn{\alpha} is given its conjugate Dirichlet prior. The default is
#'   \code{rep(1, n_comp)}, which is the noninformative uniform prior over the
#'   \code{n_comp} simplex.
#' @param joint_kp_lam Logical; If \code{TRUE}, the parameters \code{kp} and
#'   \code{lam} are drawn jointly. This can be beneficial if these are strongly
#'   correlated.
#' @param verbose Integer up to 4; Determines the amount of printed debug
#'   information.
#' @param lam_bw Numeric; the maximum distance from the current lambda at which
#'   uniform proposals are drawn.
#' @param kp_bw Numeric; A tuning parameter for kappa proposals. If \code{kp_bw
#'   == 1}, the chi-square distribution is used. Often, this distribution is too
#'   wide, so this parameter can be set to \code{0 < kp_bw < 1} to use a gamma
#'   proposal which has lower variance than the chi-square.
#' @param compute_variance Logical; Whether to add circular variance to the
#'   returned mcmc sample.
#'
#' @return A numeric matrix of sampled parameter values.
#' @export
#'
#' @examples
#' x <- rinvbatmix(100)
#' mcmcBatscheletMixture(x, Q = 10)
mcmcBatscheletMixture <- function(x, Q = 1000,
                                  burnin = 0, thin = 1,
                                  n_comp  = 4,
                                  bat_type = "inverse",
                                  init_pmat  = matrix(NA, n_comp, 4),
                                  fixed_pmat = matrix(NA, n_comp, 4),
                                  joint_kp_lam = FALSE,
                                  kp_bw  = 1,
                                  lam_bw = .05,
                                  mu_logprior_fun   = function(mu)   -log(2*pi),
                                  kp_logprior_fun   = function(kp)   1,
                                  lam_logprior_fun  = function(lam)  -log(2),
                                  alph_prior_param  = rep(1, n_comp),
                                  compute_variance  = TRUE,
                                  verbose = 0) {

  # Select Batschelet type
  if (bat_type == "inverse") {
    dbat_fun <- dinvbat
    llbat    <- likinvbat
    tlam_fun <- t_lam
  } else if (bat_type == "power") {
    dbat_fun <- dpowbat
    llbat    <- likpowbat
    tlam_fun <- tpow_lam
  } else {
    stop("Unknown Batschelet type.")
  }

  # A matrix with logicals for each initial value of the parameter matrix.
  na_fixedpmat <- is.na(fixed_pmat)
  na_initpmat  <- is.na(init_pmat)

  if (any(!na_fixedpmat[, 2] & fixed_pmat[,2] < 0))      stop("Invalid fixed kappa value.")
  if (any(!na_fixedpmat[, 3] & abs(fixed_pmat[,3]) > 1)) stop("Invalid fixed lambda value.")


  # Set initial values if the initial parameter matrix is not given for that parameter (has NAs).
  init_pmat[, 1] <- ifelse(na_initpmat[, 1], seq(0, 2*pi, length.out = n_comp + 1)[-1], init_pmat[, 1])
  init_pmat[, 2] <- ifelse(na_initpmat[, 2], rep(5, n_comp), init_pmat[, 2])
  init_pmat[, 3] <- ifelse(na_initpmat[, 3], rep(0, n_comp), init_pmat[, 3])
  init_pmat[, 4] <- ifelse(na_initpmat[, 4], rep(1/n_comp, n_comp), init_pmat[, 4])

  # Force x to be in range -pi, pi.
  x <- force_neg_pi_pi(x)
  n <- length(x)

  # Initialize parameters
  mu_cur   <- init_pmat[, 1]
  kp_cur   <- init_pmat[, 2]
  lam_cur  <- init_pmat[, 3]
  alph_cur <- init_pmat[, 4]

  # Matrix of acceptance totals for kp and lam, which will be divided by the
  # total later.
  acc_mat <- matrix(0, nrow = n_comp, ncol = 2)
  colnames(acc_mat) <- c("kp", "lam")
  rownames(acc_mat) <- 1:n_comp


  # Initialize latent group labeling
  z_cur <- integer(n)

  output_matrix <- matrix(NA, nrow = Q, ncol = n_comp*4)
  colnames(output_matrix) <- c(paste0("mu_", 1:n_comp), paste0("kp_", 1:n_comp),
                               paste0("lam_", 1:n_comp), paste0("alph_", 1:n_comp))


  if (compute_variance) {
    variance_matrix <-  matrix(NA, nrow = Q, ncol = n_comp*3)
    colnames(variance_matrix) <- c(paste0("mean_res_len_", 1:n_comp),
                                   paste0("circ_var_", 1:n_comp),
                                   paste0("circ_sd_", 1:n_comp))
  }


  Qbythin <- Q * thin + burnin


  if (verbose)     cat("Starting MCMC sampling.\n")
  if (verbose > 1) cat("Iteration:\n")


  for (i in 1:Qbythin) {

    if (verbose > 1) cat(sprintf("%5s, ", i))


    ### Sample group assignments z
    # Probability of each component for each data point.
    W <- sapply(1:n_comp, function(k) {
      alph_cur[k] * dbat_fun(x, mu_cur[k], kp_cur[k], lam_cur[k])
    })

    w_rowsum <- rowSums(W)

    # If some datapoints are numerically equal to zero, assign it equally to
    # each component. Hopefully, this does not happen again in the next
    # iteration.
    W[w_rowsum == 0, ] <- 1/n_comp
    w_rowsum[w_rowsum == 0] <- 1

    W <- W / w_rowsum

    # Randomly sample group assignments from the component probabilities.
    z_cur <- apply(W, 1, function(row_probs) sample(x = 1:n_comp, size = 1, prob = row_probs))

    # Sample weights alph
    dir_asum <- sapply(1:n_comp, function(j) sum(z_cur == j))
    alph_cur <- ifelse(na_fixedpmat[, 4], MCMCpack::rdirichlet(1, dir_asum + alph_prior_param), alph_cur)

    # After sampling the current group assignments, the parameters for each
    # component only can be sampled separately.
    for (j in 1:n_comp) {

      if (verbose > 2) cat(j)

      # Dataset assigned to this component.
      x_j <- x[z_cur == j]

      # Check whether anything is assigned to this components. If not, don't update the parameters.
      if (length(x_j) == 0) {
        if (verbose > 1) cat("---")
        next
      }

      if (verbose > 2) cat("m")

      # Sample mu
      if (na_fixedpmat[j, 1]) {
        mu_cur[j]  <- sample_mu_bat_2(x_j, mu_cur[j], kp_cur[j], lam_cur[j],
                                      tlam_fun, mu_logprior_fun)
      }

      if (joint_kp_lam) {

        if (verbose > 2) cat("kl")

        kplam_curj <- sample_kp_and_lam_bat(x_j,
                                            mu_cur[j], kp_cur[j], lam_cur[j],
                                            llbat, lam_bw = lam_bw,
                                            kp_logprior_fun, lam_logprior_fun,
                                            var_tune = kp_bw)

        # Assign the new values if they are new, and set acceptance count.
        if (kp_cur[j] != kplam_curj[1])  {
          kp_cur[j]     <- kplam_curj[1]
          acc_mat[j, 1] <- acc_mat[j, 1] + 1
        }
        if (lam_cur[j] != kplam_curj[2])  {
          lam_cur[j]    <- kplam_curj[2]
          acc_mat[j, 2] <- acc_mat[j, 2] + 1
        }


      } else {

        if (verbose > 2) cat("k")
        # Sample kp
        if (na_fixedpmat[j, 2]) {
          kp_new <- sample_kp_bat(x_j, mu_cur[j], kp_cur[j], lam_cur[j],
                                      llbat, kp_logprior_fun, var_tune = kp_bw)
          if (kp_cur[j] != kp_new)  {
            kp_cur[j]     <- kp_new
            acc_mat[j, 1] <- acc_mat[j, 1] + 1
          }
        }

        if (verbose > 2) cat("l")
        # Sample lam
        if (na_fixedpmat[j, 3]) {
          lam_new <- sample_lam_bat(x_j, mu_cur[j], kp_cur[j], lam_cur[j],
                                       llbat, lam_logprior_fun, lam_bw = lam_bw)

          if (lam_cur[j] != lam_new)  {
            lam_cur[j]    <- lam_new
            acc_mat[j, 2] <- acc_mat[j, 2] + 1
          }
        }
      }
    }


    # Possibly extremely detailed debugging information.
    if (verbose > 3) {
      cat("\n",
          sprintf("mu:   %8s, ", round(mu_cur, 3)), "\n",
          sprintf("kp:   %8s, ", round(kp_cur, 3)), "\n",
          sprintf("lam:  %8s, ", round(lam_cur, 3)), "\n",
          sprintf("alph: %8s, ", round(alph_cur, 3)), "\n")
    }

    if (i %% 50 == 0 && verbose == 1) cat(i, ", \n")
    if (i %% 5 == 0 && verbose > 1) cat("\n")

    if (i %% thin == 0 && i >= burnin) {
      isav <- (i - burnin) / thin
      output_matrix[isav, ] <- c(mu_cur, kp_cur, lam_cur, alph_cur)

      if (compute_variance) {
        # Compute mean resultant lengths for each component.
        R_bar_cur <- sapply(1:n_comp, function(i) {
          computeMeanResultantLengthBat(kp_cur[i], lam_cur[i], bat_type)
        })

        # Compute variances.
        circ_var_cur <- 1 - R_bar_cur
        circ_sd_cur  <- computeCircSD(R_bar_cur)

        # Put the results in an output matrix.
        variance_matrix[isav, ] <- c(R_bar_cur, circ_var_cur, circ_sd_cur)
      }
    }
  }

  # Compute acceptance ratio.
  acc_mat <- acc_mat / Q

  if (verbose) cat("\nFinished.\n")


  if (compute_variance) {
    return(list(
      mcmc_sample = coda::mcmc(cbind(output_matrix, variance_matrix),
                               start = burnin + 1, end = Qbythin, thin = thin),
      acceptance_rates = acc_mat))
  } else {
    return(list(
      mcmc_sample = coda::mcmc(output_matrix,
                               start = burnin + 1, end = Qbythin, thin = thin),
      acceptance_rates = acc_mat))
  }
}





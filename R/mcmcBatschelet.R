
# The right hand side of the log likelihood of a Batschelet-type distribution
ll_rhs_bat <- function(x, mu, kp, lam, tlam_fun) {
  kp * sum(cos(tlam_fun(x - mu, lam)))
}


# The left hand side of the log likelihood of a Batschelet-type distribution
ll_lhs_bat <- function(x, mu, kp, lam, tlam_fun) {
  kp * sum(cos(tlam_fun(x - mu, lam)))
}


sample_mu_batmix <- function(x, mu_cur, kp, lam, tlam_fun) {

  # Sample a candidate from the distribution of mu when we have a von Mises
  # distribution.
  C_j    <- sum(cos(x))
  S_j    <- sum(sin(x))
  R_j    <- sqrt(C_j^2 + S_j^2)
  mu_hat <- atan2(S_j, C_j)
  mu_can <- suppressWarnings(as.numeric(circular::rvonmises(1, mu_hat, R_j * kp))) # This line should be replaced with rvmc

  ll_can <- ll_rhs_bat(x, mu_can, kp, lam, tlam_fun)
  ll_cur <- ll_rhs_bat(x, mu_cur, kp, lam, tlam_fun)
}





mcmcBatscheletMixture <- function(x, Q = 1000,
                                  n_comp  = 4,
                                  bat_type = "inverse",
                                  init_pmat  = matrix(NA, n_comp, 4),
                                  fixed_pmat = matrix(NA, n_comp, 4)) {

  # Select Batschelet type
  if (bat_type == "inverse") {
    dbat_fun      <- dinvbat
    likfunbat_fun <- likfuninvbat
    tlam_fun      <- t_lam
  } else if (bat_type == "power") {
    dbat_fun      <- dpowbat
    likfunbat_fun <- likfunpowbat
    tlam_fun      <- tpow_lam
  } else {
    stop("Unknown Batschelet type.")
  }

  # A matrix with logicals for each initial value of the parameter matrix.
  na_fixedpmat <- is.na(fixed_pmat)
  na_initpmat  <- is.na(init_pmat)

  # Set initial values if the initial parameter matrix is not given for that parameter (has NAs).
  if (any(na_initpmat[, 1])) init_pmat[, 1] <- seq(0, 2*pi, length.out = n_comp + 1)[-1]
  if (any(na_initpmat[, 2])) init_pmat[, 2] <- rep(5, n_comp)
  if (any(na_initpmat[, 3])) init_pmat[, 3] <- rep(0, n_comp)
  if (any(na_initpmat[, 4])) init_pmat[, 4] <- rep(1/n_comp, n_comp)

  # Force x to be in range -pi, pi.
  x <- force_neg_pi_pi(x)


  n <- length(x)

  # Initialize parameters
  mu_cur   <- init_pmat[, 1]
  kp_cur   <- init_pmat[, 2]
  lam_cur  <- init_pmat[, 3]
  alph_cur <- init_pmat[, 4]

  # Initialize latent group labeling
  z_cur <- integer(n)

  # initialize W matrix.
  W <- matrix(alph_cur, nrow = n, ncol = n_comp, byrow = TRUE)


  for (i in 2:Q) {

    ### Sample group assignments z

    # Probability of each component for each data point.
    W <- t(sapply(x, function(xi) {
      sapply(1:n_comp, function(k) {
        alph_cur[k] * dbat_fun(xi, mu_cur[k], kp_cur[k], lam_cur[k])
      })
    }))
    W <- W / rowSums(W)

    # Randomly sample group assignments from the component probabilities.
    z_cur <- apply(W, 1, function(row_probs) sample(x = 1:n_comp, size = 1, prob = row_probs))

    # Sample weights alph
    dir_asum <- sapply(1:n_comp, function(j) sum(z_cur == j))
    alph_cur <- ifelse(na_fixedpmat[, 4], MCMCpack::rdirichlet(1, dir_asum + 1), alph_cur)

    # After sampling the current group assignments, the parameters for each
    # component only can be sampled separately.
    for (j in 1:n_comp) {

      # Dataset assigned to this component.
      x_j <- x[z_cur == j]

      # Sample mu
      mu_cur    <- sample_mu_bat(x_j, mu_cur[j], kp_cur[j], lam_cur[j], tlam_fun)



      # Sample kp


      # Sample lam
    }

  }
}





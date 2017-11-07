
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
  mu_hat <- atan2(S_j, C_j)
  mu_can <- circular:::RvonmisesRad(1, mu_cur, R_j * kp) # This line should be replaced with rvmc

  ll_can <- ll_rhs_bat(x, mu_can, kp, lam, tlam_fun)
  ll_cur <- ll_rhs_bat(x, mu_cur, kp, lam, tlam_fun)

  # The proposal is von Mises and thus symmetric, so the transition probabilities of MH are omitted here.
  mu_lograt <- ll_can + mu_logprior_fun(mu_can) - ll_cur - mu_logprior_fun(mu)

  if (mu_lograt > log(runif(1))) {
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
  mu_can <- circular:::RvonmisesRad(1, mu_hat, R_j * kp) # This line should be replaced with rvmc

  ll_can <- ll_rhs_bat(x, mu_can, kp, lam, tlam_fun)
  ll_cur <- ll_rhs_bat(x, mu_cur, kp, lam, tlam_fun)

  logp_mu_can_to_cur <- circular:::DvonmisesRad(mu_cur, mu_hat, kp, log = TRUE)
  logp_mu_cur_to_can <- circular:::DvonmisesRad(mu_can, mu_hat, kp, log = TRUE)

  mu_lograt <- ll_can + mu_logprior_fun(mu_can) + logp_mu_can_to_cur - ll_cur - mu_logprior_fun(mu_cur) - logp_mu_cur_to_can

  if (mu_lograt > log(runif(1))) {
    return(mu_can)
  } else {
    return(mu_cur)
  }
}


sample_kp_bat <- function(x, mu, kp_cur, lam, llbat, kp_logprior_fun) {

  # Sample a candidate
  kp_can <- rchisq(1, df = kp_cur)

  ll_can <- llbat(x, mu, kp_can, lam, log = TRUE)
  ll_cur <- llbat(x, mu, kp_cur, lam, log = TRUE)

  logp_kp_can_to_cur <- dchisq(kp_cur, kp_can, log = TRUE)
  logp_kp_cur_to_can <- dchisq(kp_can, kp_cur, log = TRUE)

  kp_lograt <- ll_can + kp_logprior_fun(kp_can) + logp_kp_can_to_cur - ll_cur - kp_logprior_fun(kp_cur) - logp_kp_cur_to_can

  if (kp_lograt > log(runif(1))) {
    return(kp_can)
  } else {
    return(kp_cur)
  }
}




sample_lam_bat <- function(x, mu, kp, lam_cur, llbat, lam_logprior_fun, lam_bw = .05) {

  # Sample a candidate
  lam_can <- runif(1, max(-1, lam_cur - lam_bw), min(1, lam_cur + lam_bw))

  ll_can <- ll_rhs_bat(x, mu, kp, lam_can, tlam_fun)
  ll_cur <- ll_rhs_bat(x, mu, kp, lam_cur, tlam_fun)

  logp_lam_can_to_cur <- dunif(lam_cur, max(-1, lam_can - lam_bw), min(1, lam_can + lam_bw), log = TRUE)
  logp_lam_cur_to_can <- dunif(lam_can, max(-1, lam_cur - lam_bw), min(1, lam_cur + lam_bw), log = TRUE)

  lam_lograt <- ll_can + lam_logprior_fun(lam_can) + logp_lam_can_to_cur -
    ll_cur - lam_logprior_fun(lam_cur) - logp_lam_cur_to_can

  if (lam_lograt > log(runif(1))) {
    return(lam_can)
  } else {
    return(lam_cur)
  }
}

mcmcBatscheletMixture <- function(x, Q = 1000,
                                  burnin = 0, thin = 1,
                                  n_comp  = 4,
                                  bat_type = "inverse",
                                  init_pmat  = matrix(NA, n_comp, 4),
                                  fixed_pmat = matrix(NA, n_comp, 4),
                                  mu_logprior_fun   = function(mu)   -log(2*pi),
                                  kp_logprior_fun   = function(kp)   1,
                                  lam_logprior_fun  = function(lam)  -log(2),
                                  alph_prior_param  = rep(1, n_comp)
                                  ) {

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

  Qbythin <- Q * thin + burnin



  for (i in 2:Qbythin) {

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
    alph_cur <- ifelse(na_fixedpmat[, 4], MCMCpack::rdirichlet(1, dir_asum + alph_prior_param), alph_cur)

    # After sampling the current group assignments, the parameters for each
    # component only can be sampled separately.
    for (j in 1:n_comp) {

      # Dataset assigned to this component.
      x_j <- x[z_cur == j]

      # Sample mu
      mu_cur[j]  <- sample_mu_bat(x_j, mu_cur[j], kp_cur[j], lam_cur[j], tlam_fun, mu_logprior_fun)

      # Sample kp
      kp_cur[j]  <- sample_kp_bat(x_j, mu_cur[j], kp_cur[j], lam_cur[j], llbat, kp_logprior_fun)

      # Sample lam
      lam_cur[j] <- sample_lam_bat(x_j, mu_cur[j], kp_cur[j], lam_cur[j], llbat, lam_logprior_fun)

    }

    if (i %% thin == 0 && i>= burnin) {

    }

  }


}





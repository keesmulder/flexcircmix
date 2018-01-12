# These tests were used to develop and tweak parts of the package. They are
# given here by means of archive.
context("Temporary development tests")



test_that("tpow_lam is comparable to tinv_lam", {
  skip("Initial testing")

  # This code was previously used to tune the tpow_lam function to be as close
  # as possible to the tinv_lam function.

  # Compute the difference between the power and inverse Batschelet functions.
  t_diff <- function(inv_lam, pow_lam, res = 100) {
    tinvsq <- t_lam(seq(-pi, pi, length.out = res), inv_lam)
    tpowsq <- tpow_lam(seq(-pi, pi, length.out = res), pow_lam)

    mean(abs(abs(tpowsq) - abs(tinvsq)))
  }


  # Draw two curves
  t_diff_curve <- function(inv_lam, pow_lam) {
    curve(t_lam(x, inv_lam),    -2*pi, 2*pi, asp = 1, col = "darkolivegreen")
    curve(tpow_lam(x, pow_lam), -2*pi, 2*pi, asp = 1, col = "skyblue", add = TRUE)
  }

  # Find the minimum difference between the two curves, at inv_lam = 1 and inv_lam = -1.
  lam_min     <-  optimize(function(x) t_diff(1, x), interval = c(-1, 1))$minimum
  lam_min_neg <- -optimize(function(x) t_diff(-1, x), interval = c(-1, 1))$minimum

  mean_lam_min <- mean(c(lam_min, lam_min_neg))

  expect_true(t_diff(1, lam_min) < t_diff(1, lam_min + .001))
  expect_true(t_diff(1, lam_min) < t_diff(1, lam_min - .001))

  expect_true(t_diff(-1, -lam_min_neg) < t_diff(-1, -lam_min_neg + .001))
  expect_true(t_diff(-1, -lam_min_neg) < t_diff(-1, -lam_min_neg - .001))

  t_diff_curve(1/3, mean_lam_min/3)
  t_diff_curve(1, mean_lam_min)
})



test_that("dinvbat is comparable to dpowbat", {
  skip("Initial testing")

  # Draw two curves
  dbat_diff_curve <- function(inv_lam, mu = 1, kp = 5, lam = .5) {
    curve(dinvbat(x, mu, kp, lam),  -2*pi, 2*pi, col = "darkolivegreen", n = 500)
    curve(dpowbat(x, mu, kp, lam),  -2*pi, 2*pi, col = "skyblue", add = TRUE, n = 500)
  }

})



test_that("Bootrap paralellization works", {
  skip("Skip all devepment tests.")
  skip("Skip parallel tests because they fail in R CMD Check.")

  # Test with parallelization
  x <- rinvbatmix(200)

  expect_error(fit_boot_2 <- fitbatmix(x, method = "boot", B = 3,
                                       parallel = TRUE), NA)
  expect_error(fit_boot_2, NA)
  expect_error(summary(fit_boot_2), NA)
})



test_that("Can compute sensible functions of parameters", {

  skip("Skip all devepment tests.")

  kp <- 4
  lam <- 0

  computeMeanResultantLengthBat(kp, lam)


  kp <- 4
  lam <- 0.2

  computeMeanResultantLengthBat(kp, lam)

  integrate(function(x) dinvbatkern(x, 0, kp, lam), -pi, pi)
  K_kplam(kp, lam)

  integrate(function(x) dpowbatkern(x, 0, kp, lam), -pi, pi)
  powbat_nc(kp, lam)



  cos_fun <- function(theta) cos(theta) * dbat_fun(theta, 0, kp, lam, log = FALSE)

  nc_batfun(kp, lam) * integrate(cos_fun, -pi, pi)$value


  cos_fun_withnc <- function(theta) {
    nc_batfun(kp, lam) *
      cos(theta) *
      dbatkernfun(theta, 0, kp, lam, log = FALSE)
  }

  integrate(function(x)  dbatkernfun(x, 0, kp, lam)/ nc_batfun(kp, lam), -pi, pi)
  integrate(function(x) dbat_fun(x, 0, kp, lam), -pi, pi)
  dbat_fun(.1, 0, kp, lam)


  integrate(cos_fun_withnc, -pi, pi)$value
})



test_that("Proposals similar to full conditionals, log-likelihood shape", {

  skip("Skip all development tests.")

  mu = 1; kp = 6; lam = .6

  dat   <- rinvbat(100, mu, kp, lam)

  mylik <- likfuninvbat(dat, log = FALSE)

  mull <- Vectorize(function(x) flexcircmix:::ll_rhs_bat(dat, x,  kp, lam, t_lam))
  kpll <- Vectorize(function(x) mylik(mu, x, lam))
  lmll <- Vectorize(function(x) mylik(mu, kp, x))

  curve(mull, -pi, pi)
  curve(kpll, 0, 10, n = 200)
  curve(lmll, -1, 1)



  dat <- rinvbatmix(300, mus = c(-1, 2), kps = c(5, 10), lams = c(-.3, .6), alphs = c(.3, .7))

  hist(dat, breaks = 30)

  curve(ll_rhs_bat(dat, mu_cur, x, lam_cur, t_lam), 0, 20)


})



test_that("Initial MCMC tests", {

  skip("Initial mcmc testing")

  n <- 1000
  mus1 <- replicate(n, sample_mu_bat(x_j, mu_cur[j], kp_cur[j], lam_cur[j], tlam_fun, mu_logprior_fun))
  mus2 <- replicate(n, sample_mu_bat_2(x_j, mu_cur[j], kp_cur[j], lam_cur[j], tlam_fun, mu_logprior_fun))
  plot(density(mus1))
  plot(density(mus2))
  lines(density(mus2))



  # MCMC default arguments for testing
  x <-  rinvbatmix(100)
  Q = 1000
  burnin = 0
  thin = 1
  n_comp  = 4
  bat_type = "inverse"
  init_pmat  = matrix(NA, n_comp, 4)
  fixed_pmat = matrix(NA, n_comp, 4)
  mu_logprior_fun   = function(mu)   -log(2*pi)
  kp_logprior_fun   = function(kp)   1
  lam_logprior_fun  = function(lam)  -log(2)
  alph_prior_param  = rep(1, n_comp)
  i <- 2
  j <- 1
  verbose = TRUE


  skip("Long-chain MCMC testing skipped for because they take too long.")
  x <-  rinvbatmix(100)

  set.seed(1)
  sam_pow <- mcmcBatscheletMixture(x, Q = 1000, thin = 10, bat_type = 'power', verbose = 4)
  plot(sam_pow)

  plot_batmix_sample(x = x, param = sam_pow)


})



test_that("Improvements for kappa proposal", {


  makeKpLik <- function(mu = 1, lam = .4, n = 100, true_kp = 5, log = FALSE) {
    x <- rinvbatmix(n, mu, true_kp, lam, 1)

    llfun <- likfunpowbat(x, log = log)
    Vectorize(function(kp) {llfun(mu, kp, lam)})
  }



  kplik1 <- makeKpLik(1, .6, 100, 5)
  curve(kplik1, 0, 10)

  # with a chi-square proposal, the proposal for a new velue of kp would be:
  curve(dchisq(x, 5), 0, 10)
  # Which is way too wide by comparison.


  circglmbayes::

})

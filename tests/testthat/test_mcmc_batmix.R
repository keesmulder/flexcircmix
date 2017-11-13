library(flexcircmix)

context("Test MCMC functions")




test_that("Proposals similar to full conditionals", {

  skip("Initial testing of log-likelihood shape")


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

test_that("MCMC runs", {

  x <-  rinvbatmix(100)

  time_inv <- system.time(sam_inv <- mcmcBatscheletMixture(x, Q = 10, bat_type = 'inverse'))
  time_pow <- system.time(sam_pow <- mcmcBatscheletMixture(x, Q = 10, bat_type = 'power'))

  time_inv
  time_pow

  plot_batmix_sample(x = x, param = sam_inv)
  plot_batmix_sample(x = x, param = sam_pow)


  skip("Long-chain MCMC testing skipped for because they take too long.")

  sam_pow <- mcmcBatscheletMixture(x, Q = 100, thin = 10, bat_type = 'power', verbose = 0)
  plot(sam_pow)

  plot_batmix_sample(x = x, param = sam_pow)

})





test_that("MCMC runs", {

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


})
library(flexcircmix)

context("Test MCMC functions")




test_that("Proposals similar to full conditionals", {

  skip("Initial testing of log-likelihood shape")


  mu = 1; kp = 6; lam = .6
  dat   <- rinvbat(100, mu, kp, lam)
  mylik <- likfuninvbat(dat, log = FALSE)

  mull <- Vectorize(function(x) ll_rhs_bat(dat, x,  kp, lam, t_lam))
  kpll <- Vectorize(function(x) log(mylik(mu, x, lam)))
  lmll <- Vectorize(function(x) ll_rhs_bat(dat, mu, kp, x, t_lam))

  curve(mull, -pi, pi)
  curve(kpll, 0, 50)
  curve(lmll, -pi, pi)



  dat <- rinvbatmix(300, mus = c(-1, 2), kps = c(5, 10), lams = c(-.3, .6), alphs = c(.3, .7))

  hist(dat, breaks = 30)

  curve(ll_rhs_bat(dat, mu_cur, x, lam_cur, t_lam), 0, 20)


})


test_that("MCMC runs", {


})
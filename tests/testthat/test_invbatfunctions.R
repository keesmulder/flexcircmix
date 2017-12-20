library(flexcircmix)

context("Inverse Batschelet functions")

test_that("Distribution is computed correctly", {
  expect_true(abs(s_lam(s_lam_inv(1, .3), .3) - 1) < .0001 )
  expect_true(abs(t_lam(t_lam_inv(1, .3), .3) - 1) < .0001 )

  expect_equal(dinvbat(3, mu = 2, kp = -1, lam = .2), NA)
  expect_equal(dinvbat(3, mu = 2, kp = 1, lam = 1.2), NA)

  lam <- .4
  kpsq <- seq(.2, 8, by = .2)
  plot(sapply(kpsq, function(x) flexcircmix:::K_kplam(x, lam)),
       sapply(kpsq, function(x) 1 - lam * besselI(x, 1) / besselI(x, 0)),
       xlab = "K_{kappa, lambda}", ylab = " 1 - lam * I(kp, 1) / I(kp, 0)")

})


test_that("Random generation works", {
  expect_true(is.numeric(rinvbat(10)))
  expect_true(length(rinvbat(10)) == 10)

  expect_error(rinvbat(10, 1, 1, -5))
  expect_error(rinvbat(10, 1, 1, -1))
  expect_error(rinvbat(10, 1, 1, 2))

  expect_true(is.numeric(rinvbat(10, 0, 3, -.9)))
  expect_true(is.numeric(rinvbat(10, 0, 3, 1)))

  skip("Initial checks to compare generated samples with pdf")

  mu = 1; kp = 6; lam = .6
  dat <- rinvbat(10000, mu, kp, lam)
  curve(dinvbat(x, mu, kp, lam), -pi, pi, n = 2000)
  hist(dat, breaks = 100, xlim = c(-pi, pi))

})

test_that("Likelihood functions work", {
  dat <- rinvbat(10)

  likfun <- likfuninvbat(dat, log = FALSE)

  loglikfun <- likfuninvbat(dat, log = TRUE)

  expect_equal(log(likfun(0, 1, 0)), loglikfun(0, 1, 0))
})


test_that("Optimization is sensible", {

  set.seed(10)
  x <- rinvbat(20, mu = 2, kp = 2, lam = .3)

  # Inverse Batschelet (default)
  mlpars <- maxlikbat(x)
  mlpars_fixed_mu  <- maxlikbat(x, fixed_mu = 3)
  mlpars_fixed_kp  <- maxlikbat(x, fixed_kp = 3)
  mlpars_fixed_lam <- maxlikbat(x, fixed_lam = -.3)

  expect_true(length(mlpars) == 3)
  expect_true(all(!is.na(mlpars)))

  expect_true(length(mlpars_fixed_mu) == 3)
  expect_true(all(!is.na(mlpars_fixed_mu)))

  expect_true(length(mlpars_fixed_kp) == 3)
  expect_true(all(!is.na(mlpars_fixed_kp)))

  expect_true(length(mlpars_fixed_lam) == 3)
  expect_true(all(!is.na(mlpars_fixed_lam)))

  # Test weights
  mlw <- maxlikbat(x, weights = runif(length(x)))
  expect_true(length(mlw) == 3)
  expect_true(all(!is.na(mlw)))
})

context("Power Batschelet Functions")

test_that("Distribution is computed correctly", {
  expect_true(abs(tpow_lam(tpow_lam_inv(1, .3), .3) - 1) < .0001 )

  expect_equal(dpowbat(3, mu = 2, kp = -1, lam = .2), NA)
  expect_equal(dpowbat(3, mu = 2, kp = 1, lam = 1.2), NA)

})




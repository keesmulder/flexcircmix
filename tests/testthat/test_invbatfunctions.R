library(flexcircmix)

context("Inverse Batschelet functions")

test_that("Distribution is computed correctly", {
  expect_true(abs(s_lam(s_lam_inv(1, .3), .3) - 1) < .0001 )
  expect_true(abs(t_lam(t_lam_inv(1, .3), .3) - 1) < .0001 )

  expect_equal(dinvbat(3, mu = 2, kp = -1, lam = .2), NA)
  expect_equal(dinvbat(3, mu = 2, kp = 1, lam = 1.2), NA)

})


test_that("Random generation works", {
  expect_true(is.numeric(rinvbat(10)))
  expect_true(length(rinvbat(10)) == 10)

  expect_error(rinvbat(10, 1, 1, -5))
  expect_error(rinvbat(10, 1, 1, -1))
  expect_error(rinvbat(10, 1, 1, 2))

  expect_true(is.numeric(rinvbat(10, 0, 3, -.9)))
  expect_true(is.numeric(rinvbat(10, 0, 3, 1)))
})




test_that("Optimization is sensible", {

  set.seed(10)
  x <- rinvbat(20, mu = 2, kp = 2, lam = .3)

  mlpars <- maxlikinvbat(x)
  mlpars_fixed_mu  <- maxlikinvbat(x, fixed_mu = 3)
  mlpars_fixed_kp  <- maxlikinvbat(x, fixed_kp = 3)
  mlpars_fixed_lam <- maxlikinvbat(x, fixed_lam = -.3)

  expect_true(length(mlpars) == 3)
  expect_true(all(!is.na(mlpars)))

  expect_true(length(mlpars_fixed_mu) == 3)
  expect_true(all(!is.na(mlpars_fixed_mu)))

  expect_true(length(mlpars_fixed_kp) == 3)
  expect_true(all(!is.na(mlpars_fixed_kp)))

  expect_true(length(mlpars_fixed_lam) == 3)
  expect_true(all(!is.na(mlpars_fixed_lam)))

  # Test weights
  mlw <- maxlikinvbat(x, weights = runif(length(x)))
  expect_true(length(mlw) == 3)
  expect_true(all(!is.na(mlw)))
})



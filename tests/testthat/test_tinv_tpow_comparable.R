library(flexcircmix)

context("Initial test comparability ")



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


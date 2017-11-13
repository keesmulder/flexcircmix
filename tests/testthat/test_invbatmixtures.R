library(flexcircmix)

context("Inverse Batschelet Mixtures")

test_that("Mixture is computed correctly", {

  # Inverse Batschelet
  dibm <- dbatmix(1, mus = c(-1, 1), kps = c(8, 10), lams = c(-.3, .3), alphs = c(.4, .6))
  expect_true(is.numeric(dibm))
  expect_true(length(dibm) == 1)
  expect_true(abs(dibm - 1.2) < .1)

})

test_that("Random generation works", {
  xibm <- rinvbatmix(10, mus = c(-1, 1), kps = c(8, 10), lams = c(-.3, .3), alphs = c(.4, .6))

  expect_true(length(xibm) == 10)
  expect_true(is.numeric(xibm))


  expect_error(rinvbatmix(10, lams = 1))
  expect_error(rinvbatmix(10, lams = c(-2, -1, 4)))
})



test_that("EM Algorithm works", {
  dat <- rinvbatmix(80, mus = c(-1, 1), kps = c(8, 10), lams = c(-.3, .3), alphs = c(.4, .6))
  pmat <- fitbatmix(dat, n_comp = 2, verbose = FALSE, max_its = 3)

  expect_true(is.matrix(pmat))
  expect_true(all(!is.na(pmat)))

})


context("Power Batschelet Mixtures")



test_that("Mixture is computed correctly", {

  # Inverse Batschelet
  dibm <- dbatmix(1, dbat_fun = dpowbat, mus = c(-1, 1), kps = c(8, 10), lams = c(-.3, .3), alphs = c(.4, .6))
  expect_true(is.numeric(dibm))
  expect_true(length(dibm) == 1)

})


test_that("EM Algorithm works", {


  dat <- rinvbatmix(80, mus = c(-1, 1), kps = c(4, 6), lams = c(-.3, .6), alphs = c(.4, .6))

  expect_error(pmat <- fitbatmix(dat, bat_type = "random_string"))


  # Inverse Batschelet
  pmat_inv <- fitbatmix(dat, bat_type = "inverse", n_comp = 2, verbose = FALSE, max_its = 10, ll_tol = .1)

  expect_true(is.matrix(pmat_inv))
  expect_true(all(!is.na(pmat_inv)))

  # Power Batschelet
  pmat_pow <- fitbatmix(dat, bat_type = "power", n_comp = 2, verbose = FALSE, max_its = 10, ll_tol = .1)

  expect_true(is.matrix(pmat_pow))
  expect_true(all(!is.na(pmat_pow)))

  # Test plotting
  plot_batmixfit(dat, params = pmat_inv, dbat_fun = dinvbat) + ggplot2::ggtitle("Inverse Batschelet fit")
  plot_batmixfit(dat, params = pmat_pow, dbat_fun = dpowbat) + ggplot2::ggtitle("Power Batschelet fit")

})





test_that("EM Algorithm is fast", {
  dat <- rinvbatmix(800, mus = c(-1, 1), kps = c(8, 10), lams = c(-.3, .3), alphs = c(.4, .6))

  cat(system.time(
  pmat <- fitbatmix(dat, n_comp = 2, verbose = TRUE, ll_tol = 1, max_its = 5)
  ))

  expect_true(is.matrix(pmat))
  expect_true(all(!is.na(pmat)))

})
#
#


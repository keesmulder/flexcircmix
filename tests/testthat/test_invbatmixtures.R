library(flexcircmix)

context("Inverse Batschelet Mixtures")



test_that("Mixture is computed correctly", {
  dibm <- dinvbatmix(1, mus = c(-1, 1), kps = c(8, 10), lams = c(-.3, .3), alphs = c(.4, .6))
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

  pmat <- fitinvbatmix(dat, n_comp = 2, verbose = FALSE, max_its = 3)

  expect_true(is.matrix(pmat))
  expect_true(all(!is.na(pmat)))

})



# test_that("EM Algorithm is fast", {
#   dat <- rinvbatmix(800, mus = c(-1, 1), kps = c(8, 10), lams = c(-.3, .3), alphs = c(.4, .6))
#
#   cat(system.time(
#   pmat <- fitinvbatmix(dat, n_comp = 2, verbose = TRUE, ll_tol = 1, max_its = 5)
#   ))
#
#   expect_true(is.matrix(pmat))
#   expect_true(all(!is.na(pmat)))
#
# })
#
#


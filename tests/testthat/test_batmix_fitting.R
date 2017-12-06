library(flexcircmix)

context("General tests for fitting batmix through wrapper fitbatmix")


test_that("EM fitbatmix wrapper works", {

  x <- rinvbatmix(50)


  (system.time(fit_inv <- fitbatmix(x, bat_type = "inverse", method = "EM")))

  (system.time(fit_pow <- fitbatmix(x, bat_type = "power", method = "EM")))


  expect_true(fit_inv$method == "EM")
  expect_true(fit_pow$method == "EM")

})

test_that("Bootstrap works", {

  x <- rinvbatmix(50)
  fit1 <- fitbatmix(x, method = "EM")

  # fit2 <- fitbatmix(x, method = "boot", B = 10)

})

test_that("MCMC works", {

  x <- rinvbatmix(50)
  fit1 <- fitbatmix(x, method = "EM")

  fit2 <- fitbatmix(x, method = "bayes", Q = 10)

})


library(flexcircmix)

context("General tests for fitting batmix through wrapper fitbatmix")


test_that("EM fitbatmix wrapper works", {

  x <- rinvbatmix(50)


  (system.time(fit_inv <- fitbatmix(x, bat_type = "inverse", method = "EM")))

  (system.time(fit_pow <- fitbatmix(x, bat_type = "power", method = "EM")))


  expect_true(fit_inv$method == "EM")
  expect_true(fit_pow$method == "EM")


  # Works with a single component
  expect_error(
    fit_pow <- fitbatmix(x, n_comp = 1, verbose = FALSE,
                         bat_type = "power", method = "EM")
    , NA)



})

test_that("Fixing parameters works", {

  x <- rinvbatmix(50)

  fixed_pmat = matrix(c(NA, NA, 3, NA,
                        1, NA, NA, 3.14,
                        NA, 2, NA, NA),
                      nrow = 3, ncol = 4, byrow = TRUE)

  fit_pow <- fitbatmix(x, n_comp = 3, verbose = FALSE,
                       bat_type = "power", method = "EM")
  fit_inv <- fitbatmix(x, n_comp = 3, verbose = FALSE,
                       bat_type = "inverse", method = "EM")


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


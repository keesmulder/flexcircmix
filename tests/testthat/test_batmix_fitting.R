library(flexcircmix)

context("General tests for fitting batmix through wrapper fitbatmix")


test_that("EM fitbatmix wrapper works", {

  set.seed(15)

  x <- rinvbatmix(100)

  (system.time(fit_inv <- fitbatmix(x, bat_type = "inverse", optimization_its = 2,
                                    max_its = 2, method = "EM")))

  (system.time(fit_pow <- fitbatmix(x, bat_type = "power", optimization_its = 2,
                                    max_its = 2, method = "EM")))



  expect_true(fit_inv$method == "EM")
  expect_true(fit_pow$method == "EM")


  # Works with a single component
  expect_error(
    fit_pow <- fitbatmix(x, n_comp = 1, verbose = FALSE,
                         bat_type = "power", method = "EM")
    , NA)
})

test_that("Plotting works" , {

  x <- rinvbatmix(100)

  fit_pow <- fitbatmix(x, n_comp = 3, verbose = FALSE,
                       bat_type = "power", method = "EM")
  fit_mcmc <- fitbatmix(x, n_comp = 3, verbose = 0, bat_type = "power", Q = 10, method = "bayes")

  plot_batmixfit(x, fit_pow$estimates)

  plot_batmix_sample(x, fit_mcmc$mcmc_sample, plot_n = 2)


  fit_pow <- fitbatmix(x, n_comp = 1, verbose = FALSE,
                       bat_type = "power", method = "EM")
  fit_mcmc <- fitbatmix(x, n_comp = 1, verbose = 0, bat_type = "power", Q = 10, method = "bayes")

  plot_batmixfit(x, fit_pow$estimates)
  plot_batmix_sample(x, fit_mcmc$mcmc_sample, plot_n = 2)

})

test_that("Fixing parameters works", {

  x <- rinvbatmix(50)

  fixed_pmat = matrix(c(NA, NA, .3, NA,
                        1, NA, NA, NA,
                        NA, 2, NA, NA),
                      nrow = 3, ncol = 4, byrow = TRUE)

  fit_pow <- fitbatmix(x, n_comp = 3, verbose = FALSE,
                       fixed_pmat = fixed_pmat,
                       bat_type = "power", method = "EM")


  fit_mcmc <- fitbatmix(x, n_comp = 3, verbose = 0, bat_type = "power",
                       fixed_pmat = fixed_pmat,
                        Q = 10, method = "bayes")



  expect_true(fit_pow$estimates[1, 3] == fixed_pmat[1, 3])
  expect_true(fit_pow$estimates[2, 1] == fixed_pmat[2, 1])
  expect_true(fit_pow$estimates[3, 2] == fixed_pmat[3, 2])
  expect_true(fit_pow$estimates[1, 3] == fixed_pmat[1, 3])
  expect_true(fit_pow$estimates[2, 1] == fixed_pmat[2, 1])
  expect_true(fit_pow$estimates[3, 2] == fixed_pmat[3, 2])


})


test_that("Bootstrap works", {

  x <- rinvbatmix(50)
  fit1 <- fitbatmix(x, method = "EM")

  fit2 <- fitbatmix(x, method = "boot", B = 10)


})

test_that("MCMC works", {

  x <- rinvbatmix(50)

  fit2 <- fitbatmix(x, method = "bayes", Q = 10)

})


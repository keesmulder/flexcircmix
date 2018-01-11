library(flexcircmix)

context("Test MCMC functions")


test_that("MCMC runs", {

  x <-  rinvbatmix(100)

  time_inv <- system.time(sam_inv <- mcmcBatscheletMixture(x, Q = 10, bat_type = 'inverse'))
  time_pow <- system.time(sam_pow <- mcmcBatscheletMixture(x, Q = 10, bat_type = 'power'))

  time_inv
  time_pow

  expect_error(plot_batmix_sample(x = x, param = sam_inv), NA)
  expect_error(plot_batmix_sample(x = x, param = sam_pow), NA)



})




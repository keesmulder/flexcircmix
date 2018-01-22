library(flexcircmix)

context("Test MCMC functions")


test_that("MCMC runs", {

  x <-  rinvbatmix(100)

  time_inv <- system.time(sam_inv <- mcmcBatscheletMixture(x, Q = 10,
                                                           bat_type = 'inverse'))
  time_pow <- system.time(sam_pow <- mcmcBatscheletMixture(x, Q = 10,
                                                           bat_type = 'power'))

  time_inv
  time_pow

  expect_error(plot_batmix_sample(x = x, param = sam_inv), NA)
  expect_error(plot_batmix_sample(x = x, param = sam_pow), NA)


})







test_that("Kappa proposals work", {

  expect_equal(dgammaprop(3, 10), dchisq(3, 10))
  expect_error(rgammaprop(10, 10, 2), NA)


  x <-  rinvbatmix(200)
  sam_pow_1 <- mcmcBatscheletMixture(x, Q = 10000, n_comp = 3, verbose = 1,
                                     kp_bw = 3, bat_type = 'power')
  sam_pow_2 <- mcmcBatscheletMixture(x, Q = 10000, n_comp = 3, verbose = 1,
                                     kp_bw = 1, bat_type = 'power')
  sam_pow_3 <- mcmcBatscheletMixture(x, Q = 10000,  n_comp = 3, verbose = 1,
                                     kp_bw = .2, bat_type = 'power')

  cbind(sam_pow_1$acceptance_rates[, 1],
        sam_pow_2$acceptance_rates[, 1],
        sam_pow_3$acceptance_rates[, 1])

  colMeans(cbind(sam_pow_1$acceptance_rates[, 1],
                 sam_pow_2$acceptance_rates[, 1],
                 sam_pow_3$acceptance_rates[, 1]))

})



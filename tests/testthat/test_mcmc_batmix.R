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

  x <-  rinvbatmix(20)

  # Test if all starting values can be weird.
  weird_init <- matrix(c(NA, 0, NA, NA), ncol = 4, nrow = 3, byrow = TRUE)
  expect_error(
    sam_pow_weird <- mcmcBatscheletMixture(x, Q = 10,
                                         init_pmat = weird_init,
                                         n_comp = 3, verbose = 0,
                                         kp_bw = 3, bat_type = 'power'),
    NA
  )
  # Test if all starting values can be weird.
  weird_init <- matrix(c(NA, 1000, NA, NA), ncol = 4, nrow = 3, byrow = TRUE)
  expect_error(
    sam_pow_weird <- mcmcBatscheletMixture(x, Q = 10,
                                           init_pmat = weird_init,
                                           n_comp = 3, verbose = 0,
                                           kp_bw = 3, bat_type = 'power'),
    NA
  )


  # curve(flexcircmix:::dgammaprop(x, 10, 1), 0, 20)
  # curve(flexcircmix:::dgammaprop(x, 10, 2), 0, 20)
  # curve(flexcircmix:::dgammaprop(x, 10, .1), 0, 20)

  skip("")

  x <-  rinvbatmix(50)
  sam_pow_1 <- mcmcBatscheletMixture(x, Q = 1000, n_comp = 3, verbose = 4,
                                     kp_bw = 1, bat_type = 'power')
  sam_pow_2 <- mcmcBatscheletMixture(x, Q = 1000, n_comp = 3, verbose = 4,
                                     kp_bw = .1, bat_type = 'power')
  sam_pow_3 <- mcmcBatscheletMixture(x, Q = 1000,  n_comp = 3, verbose = 4,
                                     kp_bw = .01, bat_type = 'power')

  cbind(sam_pow_1$acceptance_rates[, 1],
        sam_pow_2$acceptance_rates[, 1],
        sam_pow_3$acceptance_rates[, 1])

  colMeans(cbind(sam_pow_1$acceptance_rates[, 1],
                 sam_pow_2$acceptance_rates[, 1],
                 sam_pow_3$acceptance_rates[, 1]))

  plot.ts(sam_pow_1$mcmc_sample[,4])
  plot.ts(sam_pow_2$mcmc_sample[,4])
  plot.ts(sam_pow_3$mcmc_sample[,4])
})



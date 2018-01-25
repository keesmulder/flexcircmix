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

  parmat <- matrix(c(0, 2, 4,
                     5, 5, 5,
                     -.1, 0, .4,
                     .2, .5, .3), ncol = 4)

  x <-  rinvbatmix(1000, mus = parmat[,1], kps = parmat[,2], lams = parmat[,3], alphs = parmat[,4])
  verb <- 1
  Q <- 2000
  sam_pow_1 <- mcmcBatscheletMixture(x, Q = Q, thin = 1,  n_comp = 3, verbose = verb, init_pmat = parmat,
                                     kp_bw = 1, bat_type = 'power')
  sam_pow_2 <- mcmcBatscheletMixture(x, Q = Q, n_comp = 3, verbose = verb,
                                     kp_bw = .1, bat_type = 'power')
  sam_pow_3 <- mcmcBatscheletMixture(x, Q = Q,  n_comp = 3, verbose = verb,
                                     kp_bw = .001, bat_type = 'power')

  cbind(sam_pow_1$acceptance_rates[, 1],
        sam_pow_2$acceptance_rates[, 1],
        sam_pow_3$acceptance_rates[, 1])

  colMeans(cbind(sam_pow_1$acceptance_rates[, 1],
                 sam_pow_2$acceptance_rates[, 1],
                 sam_pow_3$acceptance_rates[, 1]))

  plot_batmix_sample(x, sam_pow_1$mcmc_sample, plot_n = 20,  dens_darkness = 5)
  plot.ts(sam_pow_1$mcmc_sample[,4:9])


  plot_batmix_sample(x, sam_pow_2$mcmc_sample, plot_n = 20,  dens_darkness = 5)
  plot.ts(sam_pow_2$mcmc_sample[,4:9])

  plot_batmix_sample(x, sam_pow_3$mcmc_sample, plot_n = 20,  dens_darkness = 5)
  plot.ts(sam_pow_3$mcmc_sample[,4:9])
  plot.ts(sam_pow_3$mcmc_sample[,16:21])


})



library(flexcircmix)

context("Test MCMC functions")


test_that("MCMC runs", {

  x <-  rinvbatmix(100)

  time_inv <- system.time(sam_inv <- mcmcBatscheletMixture(x, Q = 10,
                                                           bat_type = 'inverse'))
  time_pow <- system.time(sam_pow <- mcmcBatscheletMixture(x, Q = 10,
                                                           bat_type = 'power'))

  expect_error(plot_batmix_sample(x = x, param = sam_inv$mcmc_sample), NA)
  expect_error(plot_batmix_sample(x = x, param = sam_pow$mcmc_sample), NA)

  expect_error(sam_pow$dic, NA)
  expect_error(sam_inv$dic, NA)
})

test_that("MCMC through fitbatmix", {

  set.seed(30)

  x <-  rinvbatmix(200, kps = c(10, 10, 10))

  bmpow <- fitbatmix(x, n_comp = 3, method = "bayes", Q = 10,
                       bat_type = 'power', compute_waic = TRUE,
                       burnin = 500)

  expect_true(is.function(bmpow$log_posterior))
  expect_warning(bmpow$log_posterior(bmpow$est_vector[1:12]))
  expect_true(is.numeric(bmpow$log_posterior(bmpow$mcmc_sample[1, 1:12])))


  #Test that the environment was correctly simplified
  expect_false(any(grepl("ll_", names(environment(bmpow$log_posterior)))))



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
})


test_that("IC", {

  x <-  rinvbatmix(200, kps = c(10, 10, 10))

  sam_pow <- fitbatmix(x, n_comp = 3, method = "bayes", Q = 10,
                       bat_type = 'power', compute_waic = TRUE,
                       burnin = 500)

  expect_true(is.list(sam_pow$ic))
  expect_true(is.matrix(sam_pow$ic_mat))
  expect_true(all(!is.na(unlist(sam_pow$ic))))
})


test_that("Bridge sampling", {

  set.seed(20)

  x <-  rinvbatmix(200, kps = 2 * c(10, 10, 10))

  bmpow <- fitbatmix(x, n_comp = 3, method = "bayes", Q = 10, burnin = 2,
                       bat_type = 'power', compute_waic = FALSE)

  bmpow$mcmc_sample[, "lam_2"] <- 0
  bmpow$mcmc_sample[, "lam_3"] <- 0

  # bs <- bridge_sampler(bmpow, silent = FALSE)
  #
  # expect_true(class(bs) == "bridge")
})




test_that("Multichain", {

  set.seed(20)


  x <-  rinvbatmix(200, kps = 2 * c(10, 10, 10))

  bmpowpar <- fitbatmix(x, n_comp = 3, method = "bayes", Q = 10, burnin = 2,
                     bat_type = 'power', compute_waic = FALSE, chains = 6,
                     mcmc_parallel = TRUE)

  bmpowseq <- fitbatmix(x, n_comp = 3, method = "bayes", Q = 10, burnin = 2,
                     bat_type = 'power', compute_waic = FALSE, chains = 6,
                     mcmc_parallel = FALSE)


  expect_true(class(bmpowpar$mcmc_list) == "mcmc.list")
  expect_true(class(bmpowseq$mcmc_list) == "mcmc.list")


  # Multichain error handling
  x <-  rep(1, 100)

  bmpowpar <- fitbatmix(x, n_comp = 3, method = "bayes", Q = 10, burnin = 2,
                        bat_type = 'power', compute_waic = FALSE, chains = 6,
                        mcmc_parallel = TRUE)

  bmpowseq <- fitbatmix(x, n_comp = 3, method = "bayes", Q = 10, burnin = 2,
                        bat_type = 'power', compute_waic = FALSE, chains = 6,
                        mcmc_parallel = FALSE)


  expect_true(class(bmpowpar$mcmc_list) == "mcmc.list")
  expect_true(class(bmpowseq$mcmc_list) == "mcmc.list")


  # plot(bmpow$mcmc_list)

  # bs <- bridge_sampler(bmpow, silent = FALSE)
  #
  # expect_true(class(bs) == "bridge")

})







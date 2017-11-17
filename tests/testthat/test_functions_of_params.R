library(flexcircmix)

context("Test functions of parameters")




test_that("Functions of parameters are computed correctly", {


  # Numerical integration should equal the bessel function approach for the von Mises.
  expect_true(is.numeric(computeMeanResultantLengthBat(10, 0.2)))
  expect_true(computeMeanResultantLengthBat(10, 0.2) < 1)
  expect_true(computeMeanResultantLengthBat(10, 0.2) > 0)
  expect_true(abs(computeMeanResultantLengthVM(10) - computeMeanResultantLengthBat(10, 0)) < .0001)


  skip("Initial testing")

  kp <- 4
  lam <- 0

  computeMeanResultantLengthBat(kp, lam)


  kp <- 4
  lam <- 0.2

  computeMeanResultantLengthBat(kp, lam)

  integrate(function(x) dinvbatkern(x, 0, kp, lam), -pi, pi)
  K_kplam(kp, lam)

  integrate(function(x) dpowbatkern(x, 0, kp, lam), -pi, pi)
  powbat_nc(kp, lam)



  cos_fun <- function(theta) cos(theta) * dbat_fun(theta, 0, kp, lam, log = FALSE)

  nc_batfun(kp, lam) * integrate(cos_fun, -pi, pi)$value


  cos_fun_withnc <- function(theta) {
    nc_batfun(kp, lam) *
      cos(theta) *
      dbatkernfun(theta, 0, kp, lam, log = FALSE)
  }

  integrate(function(x)  dbatkernfun(x, 0, kp, lam)/ nc_batfun(kp, lam), -pi, pi)
  integrate(function(x) dbat_fun(x, 0, kp, lam), -pi, pi)
  dbat_fun(.1, 0, kp, lam)


  integrate(cos_fun_withnc, -pi, pi)$value

})











library(flexcircmix)

context("Test functions of parameters")




test_that("Functions of parameters are computed correctly", {


  # Numerical integration should equal the bessel function approach for the von Mises.
  expect_true(is.numeric(computeMeanResultantLengthBat(10, 0.2)))
  expect_true(computeMeanResultantLengthBat(10, 0.2) < 1)
  expect_true(computeMeanResultantLengthBat(10, 0.2) > 0)
  expect_true(abs(computeMeanResultantLengthVM(10) - computeMeanResultantLengthBat(10, 0)) < .0001)


})











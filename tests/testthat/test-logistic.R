# tests/testthat/test-logistic.R
# Comprehensive tests for R/logistic.R: logistic(), L.3(), L.4(), L.5()
# and nested functions: fct, deriv1, derivx, edfct, invfct

# ========================================================================
# Test: logistic() argument validation
# ========================================================================

test_that("logistic() errors on invalid 'names' argument", {
  expect_error(logistic(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(logistic(names = 123), "Not correct 'names' argument")
  expect_error(logistic(names = c("a", "b", "c", "d")), "Not correct 'names' argument")
})

test_that("logistic() errors on invalid 'fixed' argument", {
  expect_error(logistic(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(logistic(fixed = c(NA, NA, NA, NA)), "Not correct 'fixed' argument")
  expect_error(logistic(fixed = c(NA, NA, NA, NA, NA, NA)), "Not correct 'fixed' argument")
})

test_that("logistic() returns object of class 'Boltzmann'", {
  result <- logistic()
  expect_s3_class(result, "Boltzmann")
})

test_that("logistic() return list has correct structure", {
  result <- logistic()
  expect_true(is.function(result$fct))
  expect_true(is.function(result$ssfct))
  expect_true(is.function(result$deriv1))
  expect_null(result$deriv2)
  expect_true(is.function(result$derivx))
  expect_true(is.function(result$edfct))
  expect_true(is.function(result$inversion))
  expect_equal(result$noParm, 5)
  expect_equal(result$fixed, c(NA, NA, NA, NA, NA))
  expect_equal(result$names, c("b", "c", "d", "e", "f"))
})

test_that("logistic() default name and text are correct", {
  result <- logistic()
  expect_equal(result$name, "logistic")
  expect_equal(result$text, "Logistic (ED50 as parameter)")
})

test_that("logistic() custom fctName and fctText override defaults", {
  result <- logistic(fctName = "myModel", fctText = "my description")
  expect_equal(result$name, "myModel")
  expect_equal(result$text, "my description")
})

test_that("logistic() respects custom parameter names", {
  result <- logistic(names = c("slope", "lower", "upper", "mid", "asym"))
  expect_equal(result$names, c("slope", "lower", "upper", "mid", "asym"))
})

test_that("logistic() handles fixed parameters correctly", {
  result <- logistic(fixed = c(NA, 0, NA, NA, 1))
  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("b", "d", "e"))
  expect_equal(result$fixed, c(NA, 0, NA, NA, 1))
})

test_that("logistic() uses ssfct when provided", {
  custom_ssfct <- function(dframe) { c(1, 0, 100, 5, 1) }
  result <- logistic(ssfct = custom_ssfct)
  expect_identical(result$ssfct, custom_ssfct)
})

# ========================================================================
# Test: fct (the nonlinear function)
# ========================================================================

test_that("logistic fct produces expected values with all parameters free", {
  mod <- logistic()
  # Parameters: b, c, d, e, f
  parm <- matrix(c(-1, 0, 100, 5, 1), nrow = 1)
  
  # At dose = e = 5: response = c + (d-c)/(1+exp(b*(dose-e)))^f = 0 + 100/(1+1)^1 = 50
  result <- mod$fct(5, parm)
  expect_equal(result, 50)
  
  # At dose = 0: c + (d-c)/(1+exp(b*(-e)))^f = 0 + 100/(1+exp(5))^1
  result0 <- mod$fct(0, parm)
  expect_equal(result0, 100 / (1 + exp(5)), tolerance = 1e-10)
})

test_that("logistic fct handles fixed parameters", {
  # Fix c=0 and f=1 (like L.3 does)
  mod <- logistic(fixed = c(NA, 0, NA, NA, 1))
  parm <- matrix(c(-1, 100, 5), nrow = 1)  # b, d, e
  
  result <- mod$fct(5, parm)
  expect_equal(result, 50)
})

test_that("logistic fct handles multiple doses", {
  mod <- logistic()
  parm <- matrix(rep(c(-1, 0, 100, 5, 1), each = 3), nrow = 3)
  doses <- c(0, 5, 100)
  
  result <- mod$fct(doses, parm)
  expect_length(result, 3)
  expect_equal(result[2], 50)  # at dose = e
})

# ========================================================================
# Test: deriv1 (first derivatives in parameters)
# ========================================================================

test_that("logistic deriv1 returns correct dimensions with all parameters free", {
  mod <- logistic()
  parm <- matrix(c(-1, 0, 100, 5, 1), nrow = 1)
  
  result <- mod$deriv1(5, parm)
  # Single row: cbind()[,notFixed] drops to vector
  expect_length(result, 5)
})

test_that("logistic deriv1 returns correct dimensions with fixed parameters", {
  mod <- logistic(fixed = c(NA, 0, NA, NA, 1))
  parm <- matrix(c(-1, 100, 5), nrow = 1)
  
  result <- mod$deriv1(5, parm)
  expect_length(result, 3)  # Only free parameters
})

test_that("logistic deriv1 at dose=e has expected structure", {
  mod <- logistic()
  parm <- matrix(c(-1, 0, 100, 5, 1), nrow = 1)
  
  result <- mod$deriv1(5, parm)
  # At dose = e, exp(b*(dose-e)) = 1
  # dc/dc = 1 - 1/(1+1)^1 = 0.5
  # dd/dd = 1/(1+1)^1 = 0.5
  expect_equal(as.numeric(result[2]), 0.5, tolerance = 1e-10)  # dc
  expect_equal(as.numeric(result[3]), 0.5, tolerance = 1e-10)  # dd
})

test_that("logistic deriv1 handles multiple rows", {
  mod <- logistic()
  parm <- matrix(c(-1, 0, 100, 5, 1,
                    -2, 10, 90, 3, 1), nrow = 2, byrow = TRUE)
  doses <- c(5, 3)
  
  result <- mod$deriv1(doses, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 5)
})

# ========================================================================
# Test: derivx (first derivatives in x)
# ========================================================================

test_that("logistic derivx returns correct dimensions", {
  mod <- logistic()
  parm <- matrix(c(-1, 0, 100, 5, 1), nrow = 1)
  
  result <- mod$derivx(5, parm)
  expect_length(result, 1)
})

test_that("logistic derivx at dose=e has expected value", {
  mod <- logistic()
  # b=-1, c=0, d=100, e=5, f=1
  parm <- matrix(c(-1, 0, 100, 5, 1), nrow = 1)
  
  # At dose=e: derivx = (-f*(d-c)*1*b) / (1+1)^(f+1) = (-1*100*1*(-1)) / 4 = 25
  result <- mod$derivx(5, parm)
  expect_equal(as.numeric(result), 25, tolerance = 1e-10)
})

test_that("logistic derivx handles fixed parameters", {
  mod <- logistic(fixed = c(NA, 0, NA, NA, 1))
  parm <- matrix(c(-1, 100, 5), nrow = 1)
  
  result <- mod$derivx(5, parm)
  expect_equal(as.numeric(result), 25, tolerance = 1e-10)
})

test_that("logistic derivx handles multiple rows", {
  mod <- logistic()
  parm <- matrix(c(-1, 0, 100, 5, 1,
                    -2, 10, 90, 3, 1), nrow = 2, byrow = TRUE)
  
  result <- mod$derivx(c(5, 3), parm)
  expect_length(result, 2)
})

# ========================================================================
# Test: edfct (effective dose function)
# ========================================================================

test_that("edfct returns ED50 correctly for symmetric logistic (f=1)", {
  mod <- logistic()
  # b=-1, c=0, d=100, e=5, f=1
  # ED50 = e + log((100/50)^(1/1) - 1) / b = 5 + log(2-1)/(-1) = 5 + 0/(-1) = 5
  result <- mod$edfct(c(-1, 0, 100, 5, 1), 50)
  
  expect_equal(result[[1]], 5, tolerance = 1e-10)
  expect_true(is.numeric(result[[2]]))
  expect_length(result[[2]], 5)  # gradient for all 5 parameters
})

test_that("edfct returns correct ED values for various p", {
  mod <- logistic()
  parms <- c(-1, 0, 100, 5, 1)
  
  # ED10: e + log((100/10)^1 - 1)/(-1) = 5 + log(9)/(-1) = 5 - ln(9) ≈ 2.803
  result10 <- mod$edfct(parms, 10)
  expected10 <- 5 + log((100/10)^1 - 1) / (-1)
  expect_equal(result10[[1]], expected10, tolerance = 1e-10)
  
  # ED90: e + log((100/90)^1 - 1)/(-1) = 5 + log(1/9)/(-1) = 5 + ln(9) ≈ 7.197
  result90 <- mod$edfct(parms, 90)
  expected90 <- 5 + log((100/90)^1 - 1) / (-1)
  expect_equal(result90[[1]], expected90, tolerance = 1e-10)
})

test_that("edfct returns gradients with correct length when parameters are fixed", {
  mod <- logistic(fixed = c(NA, 0, NA, NA, 1))
  # Free: b, d, e (3 params)
  result <- mod$edfct(c(-1, 100, 5), 50)
  
  expect_equal(result[[1]], 5, tolerance = 1e-10)
  expect_length(result[[2]], 3)  # gradient for 3 free parameters
})

test_that("edfct gradient for 'e' is 1", {
  mod <- logistic()
  result <- mod$edfct(c(-1, 0, 100, 5, 1), 50)
  gradient <- result[[2]]
  
  # The 'e' gradient should be 1
  expect_equal(gradient[4], 1, tolerance = 1e-10)
})

test_that("edfct gradient for 'c' and 'd' is 0", {
  mod <- logistic()
  result <- mod$edfct(c(-1, 0, 100, 5, 1), 50)
  gradient <- result[[2]]
  
  # c and d do not appear in ED formula
  expect_equal(gradient[2], 0)
  expect_equal(gradient[3], 0)
})

test_that("edfct works with asymmetric model (f != 1)", {
  mod <- logistic()
  parms <- c(-1, 0, 100, 5, 2)
  
  # ED50 = e + log((100/50)^(1/2) - 1) / b = 5 + log(sqrt(2)-1)/(-1)
  result <- mod$edfct(parms, 50)
  expected <- 5 + log((100/50)^(1/2) - 1) / (-1)
  expect_equal(result[[1]], expected, tolerance = 1e-10)
})

# ========================================================================
# Test: invfct (inverse function)
# ========================================================================

test_that("invfct returns correct inverse for y at midpoint", {
  mod <- logistic()
  # b=-1, c=0, d=100, e=5, f=1
  # invfct(50) = log(((100-0)/(50-0))^(1/1) - 1)/(-1) + 5 = log(2-1)/(-1) + 5 = 5
  result <- mod$inversion(50, c(-1, 0, 100, 5, 1))
  expect_equal(result, 5, tolerance = 1e-10)
})

test_that("invfct is consistent with fct", {
  mod <- logistic()
  parms <- c(-1, 0, 100, 5, 1)
  
  # Compute f(3) then invert
  dose <- 3
  parm_matrix <- matrix(parms, nrow = 1)
  y <- mod$fct(dose, parm_matrix)
  x_recovered <- mod$inversion(y, parms)
  expect_equal(x_recovered, dose, tolerance = 1e-10)
})

test_that("invfct handles fixed parameters", {
  mod <- logistic(fixed = c(NA, 0, NA, NA, 1))
  # Free: b, d, e
  result <- mod$inversion(50, c(-1, 100, 5))
  expect_equal(result, 5, tolerance = 1e-10)
})

# ========================================================================
# Test: L.3()
# ========================================================================

test_that("L.3() returns Boltzmann object with correct structure", {
  result <- L.3()
  expect_s3_class(result, "Boltzmann")
  expect_equal(result$name, "L.3")
  expect_equal(result$text, "Logistic (ED50 as parameter) with lower limit fixed at 0")
  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("b", "d", "e"))
})

test_that("L.3() errors on invalid names argument", {
  expect_error(L.3(names = c("a", "b")), "Not correct names argument")
  expect_error(L.3(names = 123), "Not correct names argument")
})

test_that("L.3() errors on invalid fixed argument", {
  expect_error(L.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
  expect_error(L.3(fixed = c(NA, NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("L.3() fct gives expected response", {
  mod <- L.3()
  # b, d, e (c=0, f=1 are fixed)
  parm <- matrix(c(-1, 100, 5), nrow = 1)
  result <- mod$fct(5, parm)
  expect_equal(result, 50)
})

test_that("L.3() edfct gives correct ED50", {
  mod <- L.3()
  result <- mod$edfct(c(-1, 100, 5), 50)
  expect_equal(result[[1]], 5, tolerance = 1e-10)
  expect_length(result[[2]], 3)
})

# ========================================================================
# Test: L.4()
# ========================================================================

test_that("L.4() returns Boltzmann object with correct structure", {
  result <- L.4()
  expect_s3_class(result, "Boltzmann")
  expect_equal(result$name, "L.4")
  expect_equal(result$text, "Logistic (ED50 as parameter)")
  expect_equal(result$noParm, 4)
  expect_equal(result$names, c("b", "c", "d", "e"))
})

test_that("L.4() errors on invalid names argument", {
  expect_error(L.4(names = c("a", "b")), "Not correct names argument")
  expect_error(L.4(names = 123), "Not correct names argument")
})

test_that("L.4() errors on invalid fixed argument", {
  expect_error(L.4(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
  expect_error(L.4(fixed = c(NA, NA, NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("L.4() fct gives expected response", {
  mod <- L.4()
  parm <- matrix(c(-1, 0, 100, 5), nrow = 1)
  result <- mod$fct(5, parm)
  expect_equal(result, 50)
})

test_that("L.4() edfct gives correct ED50", {
  mod <- L.4()
  result <- mod$edfct(c(-1, 0, 100, 5), 50)
  expect_equal(result[[1]], 5, tolerance = 1e-10)
  expect_length(result[[2]], 4)
})

# ========================================================================
# Test: L.5()
# ========================================================================

test_that("L.5() returns Boltzmann object with correct structure", {
  result <- L.5()
  expect_s3_class(result, "Boltzmann")
  expect_equal(result$name, "L.5")
  expect_equal(result$text, "Generalised logistic (ED50 as parameter)")
  expect_equal(result$noParm, 5)
  expect_equal(result$names, c("b", "c", "d", "e", "f"))
})

test_that("L.5() passes additional arguments to logistic()", {
  custom_ss <- function(dframe) { c(1, 0, 100, 5, 1) }
  result <- L.5(ssfct = custom_ss)
  expect_identical(result$ssfct, custom_ss)
})

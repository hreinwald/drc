# tests/testthat/test-gompertzd.R
# Comprehensive tests for R/gompertzd.R: gompertzd()
# and nested functions: fct, ssfct, deriv1, derivx

# ========================================================================
# Test: gompertzd() argument validation
# ========================================================================

test_that("gompertzd() errors on invalid 'names' argument", {
  expect_error(gompertzd(names = c("a")), "Not correct 'names' argument")
  expect_error(gompertzd(names = 123), "Not correct 'names' argument")
  expect_error(gompertzd(names = c("a", "b", "c")), "Not correct 'names' argument")
})

test_that("gompertzd() errors on invalid 'fixed' argument", {
  expect_error(gompertzd(fixed = c(NA)), "Not correct 'fixed' argument")
  expect_error(gompertzd(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
})

# ========================================================================
# Test: gompertzd() return structure
# ========================================================================

test_that("gompertzd() returns object of class 'gompertzd'", {
  result <- gompertzd()
  expect_s3_class(result, "gompertzd")
})

test_that("gompertzd() return list has correct structure", {
  result <- gompertzd()
  expect_true(is.function(result$fct))
  expect_true(is.function(result$ssfct))
  expect_true(is.function(result$deriv1))
  expect_null(result$deriv2)
  expect_true(is.function(result$derivx))
  expect_null(result$edfct)
  expect_equal(result$noParm, 2)
  expect_equal(result$names, c("a", "b"))
})

test_that("gompertzd() default name and text are correct", {
  result <- gompertzd()
  expect_equal(result$name, "gompertzd")
  expect_equal(result$text, "Gompertz derivative")
})

test_that("gompertzd() respects custom parameter names", {
  result <- gompertzd(names = c("alpha", "beta"))
  expect_equal(result$names, c("alpha", "beta"))
})

test_that("gompertzd() handles fixed parameters correctly", {
  result <- gompertzd(fixed = c(2, NA))
  expect_equal(result$noParm, 1)
  expect_equal(result$names, "b")
})

test_that("gompertzd() handles all parameters fixed", {
  result <- gompertzd(fixed = c(2, 1))
  expect_equal(result$noParm, 0)
  expect_equal(result$names, character(0))
})

# ========================================================================
# Test: fct (the nonlinear function)
# f(x) = a * exp(bx - (a/b)*(exp(bx) - 1))
# ========================================================================

test_that("gompertzd fct produces expected values with all parameters free", {
  mod <- gompertzd()
  # Parameters: a=2, b=0.5
  parm <- matrix(c(2, 0.5), nrow = 1)

  # At dose=0: f(0) = a * exp(0 - (a/b)*(1-1)) = a * exp(0) = a = 2
  result0 <- mod$fct(0, parm)
  expect_equal(as.numeric(result0), 2, tolerance = 1e-10)

  # At dose=1: f(1) = 2*exp(0.5*1 - (2/0.5)*(exp(0.5)-1))
  dose <- 1
  a <- 2; b <- 0.5
  innerT1 <- b * dose
  innerT2 <- (a / b) * (exp(innerT1) - 1)
  expected <- a * exp(innerT1 - innerT2)
  result1 <- mod$fct(dose, parm)
  expect_equal(as.numeric(result1), expected, tolerance = 1e-10)
})

test_that("gompertzd fct handles fixed parameters", {
  # Fix a=2
  mod <- gompertzd(fixed = c(2, NA))
  parm <- matrix(0.5, nrow = 1)  # only b is free

  result <- mod$fct(0, parm)
  expect_equal(as.numeric(result), 2, tolerance = 1e-10)
})

test_that("gompertzd fct handles multiple doses", {
  mod <- gompertzd()
  parm <- matrix(rep(c(2, 0.5), each = 3), nrow = 3)
  doses <- c(0, 1, 5)

  result <- mod$fct(doses, parm)
  expect_length(result, 3)
  # At dose=0: f(0) = a = 2
  expect_equal(as.numeric(result[1]), 2, tolerance = 1e-10)
})

test_that("gompertzd fct is decreasing for a>0 and b>0", {
  mod <- gompertzd()
  parm <- matrix(c(2, 0.5), nrow = 1)

  y0 <- mod$fct(0, parm)
  y5 <- mod$fct(5, parm)
  expect_true(as.numeric(y0) > as.numeric(y5))
})

# ========================================================================
# Test: ssfct (self-starter function)
# ========================================================================

test_that("gompertzd ssfct returns correct initial estimates", {
  mod <- gompertzd()
  doses <- c(0, 1, 2, 3, 5, 8, 10)
  responses <- c(10, 8, 5, 3, 1, 0.5, 0.1)
  dframe <- data.frame(dose = doses, response = responses)

  result <- mod$ssfct(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 2)  # a, b
  # aVal should be max(y) = 10
  expect_equal(result[1], 10)
  # bVal should be 1
  expect_equal(result[2], 1)
})

# ========================================================================
# Test: deriv1 (first derivatives in parameters)
# ========================================================================

test_that("gompertzd deriv1 returns correct dimensions with all parameters free", {
  mod <- gompertzd()
  parm <- matrix(c(2, 0.5), nrow = 1)

  result <- mod$deriv1(1, parm)
  expect_length(result, 2)
})

test_that("gompertzd deriv1 returns correct dimensions with fixed parameters", {
  mod <- gompertzd(fixed = c(2, NA))
  parm <- matrix(0.5, nrow = 1)

  result <- mod$deriv1(1, parm)
  # Only 1 free parameter (b)
  expect_length(result, 1)
})

test_that("gompertzd deriv1 handles multiple rows", {
  mod <- gompertzd()
  parm <- matrix(c(2, 0.5,
                    3, 1.0), nrow = 2, byrow = TRUE)
  doses <- c(1, 2)

  result <- mod$deriv1(doses, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})

test_that("gompertzd deriv1 at dose=0 is correct for parameter a", {
  mod <- gompertzd()
  a <- 2; b <- 0.5
  parm <- matrix(c(a, b), nrow = 1)

  result <- mod$deriv1(0, parm)
  # At dose=0: fct=a, exp(b*0)=1, help3=0, help4=a/b
  # deriva = a*(1/a - 0/b) = 1
  # derivb = a*(0 + (a/b)*0/b + (a/b)*1*0) = 0
  # but wait: derivb = help1*(dose + help4*help3/b + help4*help2*dose)
  # = a*(0 + (a/b)*0/b + (a/b)*1*0) = 0
  # Actually wait, let me re-derive more carefully:
  # At dose=0:
  # help1 = fct(0, parm) = a
  # help2 = exp(b*0) = 1
  # help3 = 1 - 1 = 0
  # help4 = a/b
  # deriva = help1*(1/a - help3/b) = a*(1/a - 0) = 1
  expect_equal(as.numeric(result[1]), 1, tolerance = 1e-10)
})

# ========================================================================
# Test: derivx (first derivatives in x)
# ========================================================================

test_that("gompertzd derivx returns correct dimensions", {
  mod <- gompertzd()
  parm <- matrix(c(2, 0.5), nrow = 1)

  result <- mod$derivx(1, parm)
  expect_length(result, 1)
})

test_that("gompertzd derivx at dose=0 has expected value", {
  mod <- gompertzd()
  a <- 2; b <- 0.5
  parm <- matrix(c(a, b), nrow = 1)

  # derivx = fct(x, parm) * (b - a*exp(b*x))
  # At x=0: fct=a, derivx = a*(b - a*exp(0)) = a*(b - a) = 2*(0.5-2) = -3
  result <- mod$derivx(0, parm)
  expected <- a * (b - a)
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("gompertzd derivx handles fixed parameters", {
  mod <- gompertzd(fixed = c(2, NA))
  parm <- matrix(0.5, nrow = 1)

  result <- mod$derivx(0, parm)
  expected <- 2 * (0.5 - 2)
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("gompertzd derivx handles multiple rows", {
  mod <- gompertzd()
  parm <- matrix(c(2, 0.5,
                    3, 1.0), nrow = 2, byrow = TRUE)

  result <- mod$derivx(c(0, 1), parm)
  expect_length(result, 2)
})

# tests/testthat/test-gompertz-ssf.R
# Comprehensive tests for R/gompertz.ssf.R: gompertz.ssf()
# Internal self-starter function for the Gompertz model

# ========================================================================
# Test: gompertz.ssf() method argument matching
# ========================================================================

test_that("gompertz.ssf() defaults to method '1'", {
  ssf <- gompertz.ssf(fixed = c(NA, NA, NA, NA))
  expect_true(is.function(ssf))
})

test_that("gompertz.ssf() accepts all valid methods", {
  for (m in c("1", "2", "3", "4")) {
    ssf <- gompertz.ssf(method = m, fixed = c(NA, NA, NA, NA))
    expect_true(is.function(ssf), info = paste("Method", m, "should return a function"))
  }
})

test_that("gompertz.ssf() errors on invalid method", {
  expect_error(gompertz.ssf(method = "5", fixed = c(NA, NA, NA, NA)))
  expect_error(gompertz.ssf(method = "invalid", fixed = c(NA, NA, NA, NA)))
})

# ========================================================================
# Test: gompertz.ssf() returned closure functionality
# ========================================================================

test_that("gompertz.ssf() method='1' returns valid initial values", {
  ssf <- gompertz.ssf(method = "1", fixed = c(NA, NA, NA, NA))
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- ssf(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)  # 4 free parameters
  expect_true(all(is.finite(result)))
})

test_that("gompertz.ssf() method='2' returns valid initial values", {
  ssf <- gompertz.ssf(method = "2", fixed = c(NA, NA, NA, NA))
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- ssf(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)
  expect_true(all(is.finite(result)))
})

test_that("gompertz.ssf() method='3' returns valid initial values", {
  ssf <- gompertz.ssf(method = "3", fixed = c(NA, NA, NA, NA))
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- ssf(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)
})

test_that("gompertz.ssf() method='4' returns valid initial values", {
  ssf <- gompertz.ssf(method = "4", fixed = c(NA, NA, NA, NA))
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- ssf(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)
})

# ========================================================================
# Test: gompertz.ssf() with fixed parameters
# ========================================================================

test_that("gompertz.ssf() returns fewer values when some parameters are fixed", {
  # Fix b=0.5 (first parameter)
  ssf <- gompertz.ssf(method = "1", fixed = c(0.5, NA, NA, NA))
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- ssf(dframe)
  expect_length(result, 3)  # Only 3 free parameters (c, d, e)
})

test_that("gompertz.ssf() returns fewer values when two parameters are fixed", {
  # Fix b=0.5, c=10
  ssf <- gompertz.ssf(method = "1", fixed = c(0.5, 10, NA, NA))
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- ssf(dframe)
  expect_length(result, 2)  # Only 2 free parameters (d, e)
})

# ========================================================================
# Test: gompertz.ssf() useFixed parameter
# ========================================================================

test_that("gompertz.ssf() with useFixed=TRUE executes without error", {
  # useFixed=TRUE path is empty {} but should execute without error
  ssf <- gompertz.ssf(method = "1", fixed = c(NA, NA, NA, NA), useFixed = TRUE)
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- ssf(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)
})

# ========================================================================
# Test: gompertz.ssf() with different data patterns
# ========================================================================

test_that("gompertz.ssf() works with decreasing response data", {
  ssf <- gompertz.ssf(method = "1", fixed = c(NA, NA, NA, NA))
  # Decreasing data: high response at low dose, low at high dose
  doses <- c(0, 1, 2, 4, 8, 16, 32)
  responses <- c(98, 95, 85, 60, 20, 5, 2)
  dframe <- data.frame(dose = doses, response = responses)

  result <- ssf(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)
})

test_that("gompertz.ssf() works with increasing response data", {
  ssf <- gompertz.ssf(method = "1", fixed = c(NA, NA, NA, NA))
  # Increasing data
  doses <- c(0, 1, 2, 4, 8, 16, 32)
  responses <- c(2, 5, 20, 60, 85, 95, 98)
  dframe <- data.frame(dose = doses, response = responses)

  result <- ssf(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)
})

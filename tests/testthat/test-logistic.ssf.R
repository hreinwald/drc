# tests/testthat/test-logistic.ssf.R
# Tests for R/logistic.ssf.R: logistic.ssf() self-starter function

# ==============================================================================
# Setup: Create a realistic dose-response data frame
# ==============================================================================

# Standard decreasing logistic data for testing the returned closure
make_logistic_data <- function() {
  data.frame(
    dose = c(0, 1, 2, 3, 5, 7, 10, 15, 20),
    response = c(95, 90, 80, 60, 40, 20, 10, 5, 2)
  )
}

# ==============================================================================
# Test: logistic.ssf() returns a function (closure)
# ==============================================================================

test_that("logistic.ssf returns a closure for each method", {
  for (m in c("1", "2", "3", "4")) {
    result <- drc:::logistic.ssf(method = m, fixed = c(NA, NA, NA, NA, NA))
    expect_type(result, "closure")
    expect_true(is.function(result))
  }
})

# ==============================================================================
# Test: Calling the returned closure with method "1"
# Covers: line 14 (findbe1 + identity lambda), lines 21-34 (closure body),
#         line 8 (ytrans body via findbe1's respTr)
# ==============================================================================

test_that("logistic.ssf method '1' closure returns valid initial values", {
  ssfct <- drc:::logistic.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))
  dframe <- make_logistic_data()

  result <- ssfct(dframe)

  expect_type(result, "double")
  expect_length(result, 5)  # b, c, d, e, f (all free)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# Test: Calling the returned closure with method "2" (Anke)
# Covers: lines 9, 10 (bfct, efct bodies via findbe2)
# ==============================================================================

test_that("logistic.ssf method '2' closure returns valid initial values", {
  ssfct <- drc:::logistic.ssf(method = "2", fixed = c(NA, NA, NA, NA, NA))
  dframe <- make_logistic_data()

  result <- ssfct(dframe)

  expect_type(result, "double")
  expect_length(result, 5)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# Test: Calling the returned closure with method "3"
# ==============================================================================

test_that("logistic.ssf method '3' closure returns valid initial values", {
  ssfct <- drc:::logistic.ssf(method = "3", fixed = c(NA, NA, NA, NA, NA))
  dframe <- make_logistic_data()

  result <- ssfct(dframe)

  expect_type(result, "double")
  expect_length(result, 5)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# Test: Calling the returned closure with method "4" (Normolle)
# Also exercises bfct and efct
# ==============================================================================

test_that("logistic.ssf method '4' closure returns valid initial values", {
  ssfct <- drc:::logistic.ssf(method = "4", fixed = c(NA, NA, NA, NA, NA))
  dframe <- make_logistic_data()

  result <- ssfct(dframe)

  expect_type(result, "double")
  expect_length(result, 5)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# Test: Fixed parameters reduce the returned vector length
# Covers: line 34 subsetting with is.na(fixed)
# ==============================================================================

test_that("logistic.ssf respects fixed parameters (L.3 style: c=0, f=1)", {
  ssfct <- drc:::logistic.ssf(
    method = "1",
    fixed = c(NA, 0, NA, NA, 1)
  )
  dframe <- make_logistic_data()

  result <- ssfct(dframe)

  # Only b, d, e are free (c=0 and f=1 are fixed)
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("logistic.ssf respects fixed parameters (L.4 style: f=1)", {
  ssfct <- drc:::logistic.ssf(
    method = "2",
    fixed = c(NA, NA, NA, NA, 1)
  )
  dframe <- make_logistic_data()

  result <- ssfct(dframe)

  # Only b, c, d, e are free (f=1 is fixed)
  expect_length(result, 4)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# Test: Default method argument
# ==============================================================================

test_that("logistic.ssf defaults to method '1'", {
  ssfct_default <- drc:::logistic.ssf(fixed = c(NA, NA, NA, NA, NA))
  ssfct_explicit <- drc:::logistic.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))
  ssfct_method2 <- drc:::logistic.ssf(method = "2", fixed = c(NA, NA, NA, NA, NA))
  dframe <- make_logistic_data()

  result_default <- ssfct_default(dframe)
  result_explicit <- ssfct_explicit(dframe)
  result_method2 <- ssfct_method2(dframe)

  # Default should match method "1" exactly

  expect_equal(result_default, result_explicit)
  # And differ from method "2" to confirm method "1" is truly the default
  expect_false(isTRUE(all.equal(result_default, result_method2)))
})

# ==============================================================================
# Test: Invalid method argument
# ==============================================================================

test_that("logistic.ssf errors on invalid method", {
  expect_error(
    drc:::logistic.ssf(method = "5", fixed = c(NA, NA, NA, NA, NA)),
    "arg"
  )
})

# ==============================================================================
# Test: f parameter initial value is always 1
# ==============================================================================

test_that("logistic.ssf returns f=1 when f is free", {
  ssfct <- drc:::logistic.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))
  dframe <- make_logistic_data()

  result <- ssfct(dframe)

  # The 5th element is fVal which should be 1
  expect_equal(result[5], 1)
})

# ==============================================================================
# Test: c and d initial values are near min/max of response
# ==============================================================================

test_that("logistic.ssf c and d initial values bracket the response range", {
  ssfct <- drc:::logistic.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))
  dframe <- make_logistic_data()

  result <- ssfct(dframe)

  # result order is: b, c, d, e, f
  c_value <- result[2]  # lower asymptote
  d_value <- result[3]  # upper asymptote
  y <- dframe[, 2]
  expect_true(c_value <= min(y))
  expect_true(d_value >= max(y))
})

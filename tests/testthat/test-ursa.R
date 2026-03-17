# tests/testthat/test-ursa.R
# Comprehensive tests for R/ursa.R: ursa() and all nested functions

# ========================================================================
# Test: ursa() argument validation
# ========================================================================

test_that("ursa() errors on invalid 'names' argument - non-character", {
  expect_error(ursa(names = c(1, 2, 3, 4, 5, 6, 7)), "Not correct 'names' argument")
})

test_that("ursa() errors on invalid 'names' argument - wrong length", {
  expect_error(ursa(names = c("a", "b", "c")), "Not correct 'names' argument")
})

test_that("ursa() errors on invalid 'fixed' argument - wrong length", {
  expect_error(ursa(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
})

# ========================================================================
# Test: ursa() return structure with defaults
# ========================================================================

test_that("ursa() returns object of class 'ursa'", {
  result <- ursa()
  expect_s3_class(result, "ursa")
})

test_that("ursa() return list has correct structure", {
  result <- ursa()
  expect_type(result, "list")
  expect_true(is.function(result$fct))
  expect_true(is.function(result$ssfct))
  expect_null(result$deriv1)
  expect_null(result$deriv2)
  expect_null(result$edfct)
  expect_null(result$sifct)
  expect_equal(result$name, "ursa")
  expect_equal(result$text, "URSA")
  expect_equal(result$noParm, 7)
  expect_equal(result$names, c("b1", "b2", "c", "d", "e1", "e2", "f"))
})

# ========================================================================
# Test: ursa() with fixed parameters
# ========================================================================

test_that("ursa() handles fixed parameters correctly", {
  result <- ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA))
  expect_equal(result$noParm, 6)
  expect_equal(result$names, c("b1", "b2", "d", "e1", "e2", "f"))
})

test_that("ursa() handles multiple fixed parameters", {
  result <- ursa(fixed = c(-2, -2, 0, NA, NA, NA, 0))
  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("d", "e1", "e2"))
})

test_that("ursa() handles all parameters fixed", {
  result <- ursa(fixed = c(-2, -2, 0, 100, 5, 0.5, 0))
  expect_equal(result$noParm, 0)
  expect_length(result$names, 0)
})

# ========================================================================
# Test: ursa() with custom parameter names
# ========================================================================

test_that("ursa() respects custom parameter names", {
  custom <- c("slope1", "slope2", "lower", "upper", "ed1", "ed2", "syn")
  result <- ursa(names = custom)
  expect_equal(result$names, custom)
})

# ========================================================================
# Test: ursa() custom ssfct handling
# ========================================================================

test_that("ursa() uses custom ssfct when provided", {
  custom_ssfct <- function(dframe) {
    c(-2, -2, 0, 100, 5, 0.5, 0)
  }
  result <- ursa(ssfct = custom_ssfct)
  dframe <- data.frame(x = 1:5, y = 5:1, z = 1:5)
  expect_equal(result$ssfct(dframe), c(-2, -2, 0, 100, 5, 0.5, 0))
})

# ========================================================================
# Test: fct (the nonlinear model function)
# ========================================================================

test_that("ursa fct returns d when both doses are zero (both infinite path)", {
  # When both dose components are 0, e1/dose1 and e2/dose2 are Inf
  # The function should return d (the upper asymptote)
  result <- ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA))

  dose <- matrix(c(0, 0), ncol = 2, byrow = TRUE)
  # parm columns: b1, b2, d, e1, e2, f (c fixed at 0)
  parm <- matrix(c(-2, -2, 100, 5, 0.5, 0), ncol = 6, byrow = TRUE)

  res <- result$fct(dose, parm)
  expect_equal(res, 100)
})

test_that("ursa fct computes correct values for normal doses", {
  result <- ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA))

  # Both drugs present
  dose <- matrix(c(10, 1), ncol = 2, byrow = TRUE)
  parm <- matrix(c(-2, -2, 100, 5, 0.5, 0), ncol = 6, byrow = TRUE)

  res <- result$fct(dose, parm)
  expect_type(res, "double")
  expect_length(res, 1)
  expect_true(is.finite(res))
  # Response should be between c (0) and d (100)
  expect_gt(res, 0)
  expect_lt(res, 100)
})

test_that("ursa fct handles multiple observations", {
  result <- ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA))

  dose <- matrix(c(
    0, 0,     # both zero
    10, 0,    # only drug 1
    0, 1,     # only drug 2
    10, 1     # both drugs
  ), ncol = 2, byrow = TRUE)

  parm <- matrix(
    rep(c(-2, -2, 100, 5, 0.5, 0), 4),
    ncol = 6, byrow = TRUE
  )

  res <- result$fct(dose, parm)
  expect_length(res, 4)
  # First observation (0,0) should return d=100
  expect_equal(res[1], 100)
  # All others should be finite
  expect_true(all(is.finite(res)))
})

test_that("ursa fct returns NA when bisection fails (try-error path)", {
  result <- ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA))

  dose <- matrix(c(10, 1), ncol = 2, byrow = TRUE)
  # b1 = 0 causes 1/b1 = Inf, leading to error in bisection
  parm <- matrix(c(0, -2, 100, 5, 0.5, 0), ncol = 6, byrow = TRUE)

  res <- result$fct(dose, parm)
  expect_true(is.na(res))
})

test_that("ursa fct covers both branches of bisec if/else", {
  # The bisec function has: if (fu(fuMiddle) > 0) {fuHigh <- fuMiddle} else {fuLow <- fuMiddle}
  # Normal use of fct exercises both branches during the 25-iteration bisection
  result <- ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA))

  dose <- matrix(c(5, 0.5), ncol = 2, byrow = TRUE)
  parm <- matrix(c(-2, -2, 100, 10, 1, 0), ncol = 6, byrow = TRUE)

  res <- result$fct(dose, parm)
  expect_true(is.finite(res))
  expect_gt(res, 0)
  expect_lt(res, 100)
})

# ========================================================================
# Test: fct with only one dose component zero (one infinite path)
# ========================================================================

test_that("ursa fct handles one drug zero, other nonzero", {
  result <- ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA))

  # Drug 2 is zero → dose[,2] is 0 → parmVec[6] = e2/0 = Inf
  # Only parmVec[5] is finite → else branch taken
  dose_d2zero <- matrix(c(10, 0), ncol = 2, byrow = TRUE)
  parm <- matrix(c(-2, -2, 100, 5, 0.5, 0), ncol = 6, byrow = TRUE)

  res <- result$fct(dose_d2zero, parm)
  expect_true(is.finite(res))
  expect_gt(res, 0)
  expect_lt(res, 100)

  # Drug 1 is zero → dose[,1] is 0 → parmVec[5] = e1/0 = Inf
  dose_d1zero <- matrix(c(0, 1), ncol = 2, byrow = TRUE)

  res2 <- result$fct(dose_d1zero, parm)
  expect_true(is.finite(res2))
  expect_gt(res2, 0)
  expect_lt(res2, 100)
})

# ========================================================================
# Test: default ssfct - both branches
# ========================================================================

test_that("ursa default ssfct works with b >= 0 path", {
  # The Greco example data produces positive b from LL.4
  d1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 2, 5, 10, 20, 50, 2, 2, 2,
    2, 2, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 20, 20, 20, 20,
    20, 50, 50, 50, 50, 50)
  d2 <- c(0, 0, 0, 0.2, 0.5, 1, 2, 5, 0, 0, 0, 0, 0, 0.2,
    0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2,
    0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5)
  effect <- c(106, 99.2, 115, 79.2, 70.1, 49, 21, 3.83, 74.2,
    71.5, 48.1, 30.9, 16.3, 76.3, 48.8, 44.5, 15.5, 3.21,
    56.7, 47.5, 26.8, 16.9, 3.25, 46.7, 35.6, 21.5, 11.1,
    2.94, 24.8, 21.6, 17.3, 7.78, 1.84, 13.6, 11.1, 6.43,
    3.34, 0.89)
  dframe <- data.frame(d1, d2, effect)

  result <- ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA))
  ss <- result$ssfct(dframe)

  expect_type(ss, "double")
  expect_length(ss, 6)  # 6 free parameters (c is fixed)
  expect_true(all(is.finite(ss)))
})

test_that("ursa default ssfct works with b < 0 path", {
  # Inverted data produces negative b from LL.4
  d1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 2, 5, 10, 20, 50, 2, 2, 2,
    2, 2, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 20, 20, 20, 20,
    20, 50, 50, 50, 50, 50)
  d2 <- c(0, 0, 0, 0.2, 0.5, 1, 2, 5, 0, 0, 0, 0, 0, 0.2,
    0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2,
    0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5)
  effect_raw <- c(106, 99.2, 115, 79.2, 70.1, 49, 21, 3.83, 74.2,
    71.5, 48.1, 30.9, 16.3, 76.3, 48.8, 44.5, 15.5, 3.21,
    56.7, 47.5, 26.8, 16.9, 3.25, 46.7, 35.6, 21.5, 11.1,
    2.94, 24.8, 21.6, 17.3, 7.78, 1.84, 13.6, 11.1, 6.43,
    3.34, 0.89)
  # Invert the response to get negative b
  effect <- 120 - effect_raw
  dframe <- data.frame(d1, d2, effect)

  result <- ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA))
  ss <- result$ssfct(dframe)

  expect_type(ss, "double")
  expect_length(ss, 6)  # 6 free parameters (c is fixed)
  expect_true(all(is.finite(ss)))
})

# ========================================================================
# Test: Full model fitting with drm (integration test)
# ========================================================================

test_that("ursa works with drm for the Greco example", {
  d1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 2, 5, 10, 20, 50, 2, 2, 2,
    2, 2, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 20, 20, 20, 20,
    20, 50, 50, 50, 50, 50)
  d2 <- c(0, 0, 0, 0.2, 0.5, 1, 2, 5, 0, 0, 0, 0, 0, 0.2,
    0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2,
    0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5)
  effect <- c(106, 99.2, 115, 79.2, 70.1, 49, 21, 3.83, 74.2,
    71.5, 48.1, 30.9, 16.3, 76.3, 48.8, 44.5, 15.5, 3.21,
    56.7, 47.5, 26.8, 16.9, 3.25, 46.7, 35.6, 21.5, 11.1,
    2.94, 24.8, 21.6, 17.3, 7.78, 1.84, 13.6, 11.1, 6.43,
    3.34, 0.89)
  greco <- data.frame(d1, d2, effect)

  greco_model <- drm(effect ~ d1 + d2, data = greco,
    fct = ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA)))

  expect_s3_class(greco_model, "drc")
  expect_length(coef(greco_model), 6)
})

# ========================================================================
# Test: fct with all parameters free (no fixed)
# ========================================================================

test_that("ursa fct works with all 7 parameters free", {
  result <- ursa()

  dose <- matrix(c(10, 1), ncol = 2, byrow = TRUE)
  # All 7 parameters: b1, b2, c, d, e1, e2, f
  parm <- matrix(c(-2, -2, 0, 100, 5, 0.5, 0), ncol = 7, byrow = TRUE)

  res <- result$fct(dose, parm)
  expect_true(is.finite(res))
})

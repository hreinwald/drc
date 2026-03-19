# Tests for drmEMeventtime() and drmLOFeventtime()
# These are internal functions used for dose-response model fitting
# with an EM algorithm approach for event time data.

# --- Helper setup ---
# A simple CDF-like model: F(dose) = 1 - exp(-a * dose)
# This mimics a cumulative distribution function used in event time analysis.
simple_cdf <- function(dose, parm) {
  1 - exp(-parm[1] * dose)
}

# Test dose as a 2-column matrix (start, end intervals)
# When dose has exactly 2 columns, dose[, -1] returns a vector (not matrix)
test_dose_2col <- matrix(c(0, 1, 2, 3,    # start times (column 1)
                           1, 2, 3, Inf),  # end times (column 2)
                         ncol = 2)

# Test dose as a 3-column matrix (start, extra, end intervals)
# When dose has 3+ columns, dose[, -1] returns a matrix
test_dose_3col <- matrix(c(0, 1, 2, 3,     # column 1 (start times)
                           0.5, 1.5, 2.5, 3.5,  # column 2 (middle)
                           1, 2, 3, Inf),   # column 3 (end times)
                         ncol = 3)

test_resp <- c(5, 10, 8, 2)
test_parm <- c(0.5)


# =============================================================================
# Tests for drmEMeventtime(): structure and return value
# =============================================================================

test_that("drmEMeventtime returns a list with the correct named elements", {
  result <- drmEMeventtime(test_dose_2col, test_resp, simple_cdf)

  expect_true(is.list(result))
  expect_named(result, c("llfct", "opfct", "ssfct", "rvfct", "vcovfct", "parmfct"))
  expect_true(is.function(result$llfct))
  expect_true(is.function(result$opfct))
  expect_null(result$ssfct)
  expect_null(result$rvfct)
  expect_true(is.function(result$vcovfct))
  expect_true(is.function(result$parmfct))
})

test_that("drmEMeventtime works with doseScaling parameter", {
  result <- drmEMeventtime(test_dose_2col, test_resp, simple_cdf, doseScaling = 2)

  expect_true(is.list(result))
  expect_named(result, c("llfct", "opfct", "ssfct", "rvfct", "vcovfct", "parmfct"))
})


# =============================================================================
# Tests for opfct(): objective function (negative log-likelihood)
# =============================================================================

test_that("opfct computes a numeric scalar value with 2-column dose", {
  em <- drmEMeventtime(test_dose_2col, test_resp, simple_cdf)

  val <- em$opfct(test_parm)
  expect_true(is.numeric(val))
  expect_equal(length(val), 1)
  expect_true(is.finite(val))
})

test_that("opfct handles 3-column dose matrix (matrix branch of ifelse)", {
  em <- drmEMeventtime(test_dose_3col, test_resp, simple_cdf)

  val <- em$opfct(test_parm)
  expect_true(is.numeric(val))
  expect_equal(length(val), 1)
  expect_true(is.finite(val))
})

test_that("opfct replaces zero temp values with 1e-9", {
  # Create a model that returns constant values so Fend - Fstart = 0
  constant_model <- function(dose, parm) {
    rep(parm[1], length(dose))
  }

  dose_mat <- matrix(c(0, 1, 1, 2), ncol = 2)
  resp <- c(1, 1)

  em <- drmEMeventtime(dose_mat, resp, constant_model)
  # With a constant model, Fend - Fstart = 0, which gets replaced by 1e-9
  val <- em$opfct(c(0.5))
  expect_true(is.finite(val))
})

test_that("opfct respects doseScaling parameter", {
  em1 <- drmEMeventtime(test_dose_2col, test_resp, simple_cdf, doseScaling = 1)
  em2 <- drmEMeventtime(test_dose_2col * 2, test_resp, simple_cdf, doseScaling = 2)

  val1 <- em1$opfct(test_parm)
  val2 <- em2$opfct(test_parm)
  expect_equal(val1, val2, tolerance = 1e-10)
})

test_that("opfct handles Inf in dose end times (2-column, vector path)", {
  # With 2-column dose, dose[,-1] is a vector
  # Inf values should cause Fend to be set to 1
  dose_with_inf <- matrix(c(0, 1,    # start
                             Inf, 2), # end (first is Inf)
                           ncol = 2)
  resp <- c(3, 5)

  em <- drmEMeventtime(dose_with_inf, resp, simple_cdf)
  val <- em$opfct(test_parm)
  expect_true(is.finite(val))
})

test_that("opfct handles Inf in dose end times (3-column, matrix path)", {
  # With 3-column dose, dose[,-1] is a matrix
  # Inf values should cause Fend to be set to 1
  dose_with_inf <- matrix(c(0, 1,        # column 1
                             0.5, 1.5,    # column 2
                             Inf, 2),     # column 3 (first is Inf)
                           ncol = 3)
  resp <- c(3, 5)

  em <- drmEMeventtime(dose_with_inf, resp, simple_cdf)
  val <- em$opfct(test_parm)
  expect_true(is.finite(val))
})


# =============================================================================
# Tests for llfct(): log-likelihood function
# =============================================================================

test_that("llfct returns a numeric vector of length 2", {
  em <- drmEMeventtime(test_dose_2col, test_resp, simple_cdf)

  mock_object <- list(
    fit = list(value = 10),
    sumList = list(df.residual = 3)
  )

  result <- em$llfct(mock_object)
  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
})

test_that("llfct computes correct values", {
  em <- drmEMeventtime(test_dose_2col, test_resp, simple_cdf)

  mock_object <- list(
    fit = list(value = 10),
    sumList = list(df.residual = 3)
  )

  result <- em$llfct(mock_object)
  # First element: -object$fit$value
  expect_equal(result[1], -10)
  # Second element: df.residual
  expect_equal(result[2], 3)
})


# =============================================================================
# Tests for parmfct(): parameter extraction function
# =============================================================================

test_that("parmfct extracts par from fit object", {
  em <- drmEMeventtime(test_dose_2col, test_resp, simple_cdf)

  mock_fit <- list(par = c(0.5, 1.2, 0.8))
  result <- em$parmfct(mock_fit)
  expect_equal(result, c(0.5, 1.2, 0.8))
})

test_that("parmfct works with fixed argument", {
  em <- drmEMeventtime(test_dose_2col, test_resp, simple_cdf)

  mock_fit <- list(par = c(1, 2))
  result <- em$parmfct(mock_fit, fixed = FALSE)
  expect_equal(result, c(1, 2))
})


# =============================================================================
# Tests for vcovfct(): variance-covariance matrix function
# =============================================================================

test_that("vcovfct returns inverse of hessian", {
  em <- drmEMeventtime(test_dose_2col, test_resp, simple_cdf)

  hessian <- matrix(c(4, 1, 1, 4), nrow = 2)
  mock_object <- list(
    fit = list(hessian = hessian)
  )

  result <- em$vcovfct(mock_object)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)

  expected <- solve(hessian)
  expect_equal(result, expected, tolerance = 1e-10)
})


# =============================================================================
# Tests for drmLOFeventtime()
# =============================================================================

test_that("drmLOFeventtime returns a list with NULL elements", {
  result <- drmLOFeventtime()

  expect_true(is.list(result))
  expect_named(result, c("anovaTest", "gofTest"))
  expect_null(result$anovaTest)
  expect_null(result$gofTest)
})

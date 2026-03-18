# Tests for drmEMstandard() and drmLOFstandard()
# These are internal functions used for dose-response model fitting
# with a standard EM algorithm approach.

# --- Helper setup ---
# A simple linear dose-response model: y = a + b * dose
# This is used as the multCurves function throughout the tests.
simple_model <- function(dose, parm) {
  parm[1] + parm[2] * dose
}

# Test doses with 2 zero-dose observations and 3 non-zero
test_dose <- c(0, 0, 1, 2, 3)
test_resp <- c(10, 11, 8, 5, 3)
test_parm <- c(10, -2)


# =============================================================================
# Tests for drmEMstandard(): structure and return value
# =============================================================================

test_that("drmEMstandard returns a list with the correct named elements", {
  result <- drmEMstandard(test_dose, test_resp, simple_model)

  expect_true(is.list(result))
  expect_named(result, c("llfct", "opfct", "ssfct", "rvfct", "vcovfct", "parmfct"))
  expect_true(is.function(result$llfct))
  expect_true(is.function(result$opfct))
  expect_null(result$ssfct)
  expect_true(is.function(result$rvfct))
  expect_true(is.function(result$vcovfct))
  expect_true(is.function(result$parmfct))
})

test_that("drmEMstandard works with doseScaling parameter", {
  result <- drmEMstandard(test_dose, test_resp, simple_model, doseScaling = 2)

  expect_true(is.list(result))
  expect_named(result, c("llfct", "opfct", "ssfct", "rvfct", "vcovfct", "parmfct"))
})


# =============================================================================
# Tests for opfct(): objective function
# =============================================================================

test_that("opfct computes a numeric scalar value", {
  em <- drmEMstandard(test_dose, test_resp, simple_model)

  val <- em$opfct(test_parm)
  expect_true(is.numeric(val))
  expect_equal(length(val), 1)
  expect_true(is.finite(val))
})

test_that("opfct value is non-negative (weighted sum of squared residuals)", {
  em <- drmEMstandard(test_dose, test_resp, simple_model)

  val <- em$opfct(test_parm)
  expect_true(val >= 0)
})

test_that("opfct respects doseScaling parameter", {
  em1 <- drmEMstandard(test_dose, test_resp, simple_model, doseScaling = 1)
  em2 <- drmEMstandard(test_dose * 2, test_resp, simple_model, doseScaling = 2)

  # With scaled dose and matching doseScaling, multCurves sees the same values
  val1 <- em1$opfct(test_parm)
  val2 <- em2$opfct(test_parm)
  expect_equal(val1, val2, tolerance = 1e-10)
})


# =============================================================================
# Tests for llfct(): log-likelihood function
# =============================================================================

test_that("llfct returns a numeric vector of length 2", {
  em <- drmEMstandard(test_dose, test_resp, simple_model)

  mock_object <- list(
    fit = list(value = 10),
    sumList = list(df.residual = 3)
  )

  result <- em$llfct(mock_object)
  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
})

test_that("llfct computes correct values", {
  em <- drmEMstandard(test_dose, test_resp, simple_model)

  mock_object <- list(
    fit = list(value = 10),
    sumList = list(df.residual = 3)
  )

  result <- em$llfct(mock_object)
  # First element: -object$fit$value + sum(log(gamma(resp+1)))
  expected_ll <- -10 + sum(log(gamma(test_resp + 1)))
  expect_equal(result[1], expected_ll)
  # Second element: df.residual
  expect_equal(result[2], 3)
})


# =============================================================================
# Tests for rvfct(): residual variance function
# =============================================================================

test_that("rvfct returns a numeric scalar", {
  em <- drmEMstandard(test_dose, test_resp, simple_model)

  mock_object <- list(
    fit = list(value = 15),
    df.residual = 3  # used by stats::df.residual.default
  )

  result <- em$rvfct(mock_object)
  expect_true(is.numeric(result))
  expect_equal(length(result), 1)
  expect_equal(result, 15 / 3)
})


# =============================================================================
# Tests for parmfct(): parameter extraction function
# =============================================================================

test_that("parmfct extracts par from fit object", {
  em <- drmEMstandard(test_dose, test_resp, simple_model)

  mock_fit <- list(par = c(10, -2, 0.5))
  result <- em$parmfct(mock_fit)
  expect_equal(result, c(10, -2, 0.5))
})

test_that("parmfct works with fixed argument", {
  em <- drmEMstandard(test_dose, test_resp, simple_model)

  mock_fit <- list(par = c(1, 2))
  result <- em$parmfct(mock_fit, fixed = FALSE)
  expect_equal(result, c(1, 2))
})


# =============================================================================
# Tests for vcovfct(): variance-covariance matrix function
# =============================================================================

test_that("vcovfct returns inverse of scaled hessian when solve succeeds", {
  em <- drmEMstandard(test_dose, test_resp, simple_model)

  # Well-conditioned positive definite hessian
  hessian <- matrix(c(4, 1, 1, 4), nrow = 2)
  mock_object <- list(
    fit = list(value = 10, hessian = hessian),
    df.residual = 3
  )

  result <- em$vcovfct(mock_object)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)

  # Verify it's the inverse of the scaled hessian
  rv <- 10 / 3  # rvfct result
  scaledH <- hessian / (2 * rv)
  expected <- solve(scaledH)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("vcovfct falls back to regularized Cholesky when solve and first chol fail", {
  em <- drmEMstandard(test_dose, test_resp, simple_model)

  # Zero matrix: solve fails (singular), chol fails (not PD)
  # Regularized: 0.99 * 0 + 0.01 * I = 0.01 * I → PD → chol succeeds
  hessian <- matrix(0, nrow = 2, ncol = 2)
  mock_object <- list(
    fit = list(value = 10, hessian = hessian),
    df.residual = 3
  )

  result <- em$vcovfct(mock_object)
  # Regularized = 0.01 * I, so result = chol2inv(chol(0.01 * I)) = 100 * I
  expected <- chol2inv(chol(0.01 * diag(2)))

  expect_true(is.matrix(result))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("vcovfct returns NULL when all decompositions fail", {
  em <- drmEMstandard(test_dose, test_resp, simple_model)

  # Matrix with eigenvalue 0 and a negative eigenvalue:
  # eigenvalues: 0 and -1
  # After regularization: 0.99*0+0.01=0.01 and 0.99*(-1)+0.01=-0.98 → not PD
  hessian <- matrix(c(-0.5, 0.5, 0.5, -0.5), nrow = 2)
  mock_object <- list(
    fit = list(value = 10, hessian = hessian),
    df.residual = 3
  )

  result <- em$vcovfct(mock_object)
  expect_null(result)
})


# =============================================================================
# Tests for drmLOFstandard()
# =============================================================================

test_that("drmLOFstandard returns a list with NULL elements", {
  result <- drmLOFstandard()

  expect_true(is.list(result))
  expect_named(result, c("anovaTest", "gofTest"))
  expect_null(result$anovaTest)
  expect_null(result$gofTest)
})

# Tests for drmEMssd() and drmLOFssd()
# These are internal functions used for species sensitivity distribution (SSD)
# model fitting with an EM algorithm approach.

# --- Helper setup ---
# A simple PDF-like function (always positive) for use as multCurves
simple_pdf <- function(dose, parm) {
  dnorm(dose, mean = parm[1], sd = parm[2])
}

# A simple CDF function for use as multCurves2
simple_cdf <- function(dose, parm) {
  pnorm(dose, mean = parm[1], sd = parm[2])
}

# Test data for uncensored case
test_dose <- c(1, 2, 3, 4, 5)
test_resp <- rep(0, 5)  # resp is unused in drmEMssd
test_parm <- c(3, 1)    # mean=3, sd=1

# Test data for censored case (2-column matrix)
# Some observations are censored (dose1 != dose2), some are exact (dose1 == dose2)
test_dose_cens <- matrix(c(
  1, 1,   # exact (uncensored)
  2, 2,   # exact (uncensored)
  3, 4,   # censored interval [3, 4]
  4, 5,   # censored interval [4, 5]
  5, 5    # exact (uncensored)
), ncol = 2, byrow = TRUE)


# =============================================================================
# Tests for drmEMssd(): structure and return value
# =============================================================================

test_that("drmEMssd returns a list with the correct named elements (uncensored)", {
  result <- drmEMssd(test_dose, test_resp, simple_pdf)

  expect_true(is.list(result))
  expect_named(result, c("llfct", "opfct", "ssfct", "rvfct", "vcovfct", "parmfct"))
  expect_true(is.function(result$llfct))
  expect_true(is.function(result$opfct))
  expect_null(result$ssfct)
  expect_null(result$rvfct)
  expect_true(is.function(result$vcovfct))
  expect_true(is.function(result$parmfct))
})

test_that("drmEMssd returns a list with the correct named elements (censored)", {
  result <- drmEMssd(test_dose_cens, test_resp, simple_pdf,
                     multCurves2 = simple_cdf)

  expect_true(is.list(result))
  expect_named(result, c("llfct", "opfct", "ssfct", "rvfct", "vcovfct", "parmfct"))
  expect_true(is.function(result$llfct))
  expect_true(is.function(result$opfct))
  expect_null(result$ssfct)
  expect_null(result$rvfct)
  expect_true(is.function(result$vcovfct))
  expect_true(is.function(result$parmfct))
})

test_that("drmEMssd works with doseScaling parameter", {
  result <- drmEMssd(test_dose, test_resp, simple_pdf, doseScaling = 2)

  expect_true(is.list(result))
  expect_named(result, c("llfct", "opfct", "ssfct", "rvfct", "vcovfct", "parmfct"))
})


# =============================================================================
# Tests for opfct(): objective function (uncensored path)
# =============================================================================

test_that("opfct computes a numeric scalar value (uncensored)", {
  em <- drmEMssd(test_dose, test_resp, simple_pdf)

  val <- em$opfct(test_parm)
  expect_true(is.numeric(val))
  expect_equal(length(val), 1)
  expect_true(is.finite(val))
})

test_that("opfct computes correct negative log-likelihood (uncensored)", {
  em <- drmEMssd(test_dose, test_resp, simple_pdf)

  val <- em$opfct(test_parm)
  expected <- -sum(log(dnorm(test_dose, mean = 3, sd = 1)))
  expect_equal(val, expected, tolerance = 1e-10)
})

test_that("opfct respects doseScaling parameter (uncensored)", {
  em1 <- drmEMssd(test_dose, test_resp, simple_pdf, doseScaling = 1)
  em2 <- drmEMssd(test_dose * 2, test_resp, simple_pdf, doseScaling = 2)

  val1 <- em1$opfct(test_parm)
  val2 <- em2$opfct(test_parm)
  expect_equal(val1, val2, tolerance = 1e-10)
})


# =============================================================================
# Tests for opfct(): objective function (censored path)
# =============================================================================

test_that("opfct computes a numeric scalar value (censored)", {
  em <- drmEMssd(test_dose_cens, test_resp, simple_pdf,
                 multCurves2 = simple_cdf)

  val <- em$opfct(test_parm)
  expect_true(is.numeric(val))
  expect_equal(length(val), 1)
  expect_true(is.finite(val))
})

test_that("opfct computes correct value (censored)", {
  em <- drmEMssd(test_dose_cens, test_resp, simple_pdf,
                 multCurves2 = simple_cdf)

  val <- em$opfct(test_parm)

  # Manual computation:
  dose1 <- test_dose_cens[, 1]
  dose2 <- test_dose_cens[, 2]
  notCens <- dose1 == dose2  # rows 1, 2, 5 are exact

  # PDF values for uncensored observations
  fValues <- dnorm(dose1[notCens], mean = 3, sd = 1)
  # CDF values for censored observations
  Fvalues1 <- pnorm(dose1[!notCens], mean = 3, sd = 1)
  Fvalues2 <- pnorm(dose2[!notCens], mean = 3, sd = 1)

  expected <- -sum(log(fValues)) + (-sum(log(Fvalues2 - Fvalues1)))
  expect_equal(val, expected, tolerance = 1e-10)
})

test_that("opfct respects doseScaling parameter (censored)", {
  em1 <- drmEMssd(test_dose_cens, test_resp, simple_pdf,
                  doseScaling = 1, multCurves2 = simple_cdf)
  em2 <- drmEMssd(test_dose_cens * 2, test_resp, simple_pdf,
                  doseScaling = 2, multCurves2 = simple_cdf)

  val1 <- em1$opfct(test_parm)
  val2 <- em2$opfct(test_parm)
  expect_equal(val1, val2, tolerance = 1e-10)
})


# =============================================================================
# Tests for llfct(): log-likelihood function
# =============================================================================

test_that("llfct returns a numeric vector of length 2", {
  em <- drmEMssd(test_dose, test_resp, simple_pdf)

  mock_object <- list(
    fit = list(value = 10),
    sumList = list(df.residual = 3)
  )

  result <- em$llfct(mock_object)
  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
})

test_that("llfct computes correct values", {
  em <- drmEMssd(test_dose, test_resp, simple_pdf)

  mock_object <- list(
    fit = list(value = 10),
    sumList = list(df.residual = 3)
  )

  result <- em$llfct(mock_object)
  # First element: -object$fit$value (negated because opfct minimizes negative LL)
  expect_equal(result[1], -10)
  # Second element: df.residual
  expect_equal(result[2], 3)
})


# =============================================================================
# Tests for vcovfct(): variance-covariance matrix function
# =============================================================================

test_that("vcovfct returns inverse of hessian", {
  em <- drmEMssd(test_dose, test_resp, simple_pdf)

  hessian <- matrix(c(4, 1, 1, 4), nrow = 2)
  mock_object <- list(
    fit = list(hessian = hessian)
  )

  result <- em$vcovfct(mock_object)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(result, solve(hessian), tolerance = 1e-10)
})


# =============================================================================
# Tests for parmfct(): parameter extraction function
# =============================================================================

test_that("parmfct extracts par from fit object", {
  em <- drmEMssd(test_dose, test_resp, simple_pdf)

  mock_fit <- list(par = c(3, 1))
  result <- em$parmfct(mock_fit)
  expect_equal(result, c(3, 1))
})

test_that("parmfct works with fixed argument", {
  em <- drmEMssd(test_dose, test_resp, simple_pdf)

  mock_fit <- list(par = c(3, 1))
  result <- em$parmfct(mock_fit, fixed = FALSE)
  expect_equal(result, c(3, 1))
})


# =============================================================================
# Tests for drmLOFssd()
# =============================================================================

test_that("drmLOFssd returns a list with NULL elements", {
  result <- drmLOFssd()

  expect_true(is.list(result))
  expect_named(result, c("anovaTest", "gofTest"))
  expect_null(result$anovaTest)
  expect_null(result$gofTest)
})

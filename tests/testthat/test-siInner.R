# Tests for siInner (R/siInner.R) and fieller (R/EDcomp.R)
# These are internal functions used by EDcomp for selectivity index calculations

# Helper: create a minimal mock sifct function
# Returns a list with $val, $der (and optionally $der1, $der2, $valnum, $valden for fieller)
make_sifct <- function(val = 2.0, der = c(0.1, -0.05, 0.2, -0.1),
                       der1 = NULL, der2 = NULL,
                       valnum = NULL, valden = NULL) {
  function(parm1, parm2, pVec, jInd, kInd, reference, type) {
    list(val = val, der = der,
         der1 = der1, der2 = der2,
         valnum = valnum, valden = valden)
  }
}

# Helper: create basic test inputs for siInner
make_test_inputs <- function(interval = "none", obj_type = "continuous",
                             logBase = NULL, level = 0.95, degfree = 10) {
  npar <- 4
  ncurves <- 2
  # indexMat: maps parameters to positions

  indexMat <- matrix(1:(npar * ncurves), nrow = npar, ncol = ncurves)
  # parmMat: parameter estimates for each curve
  parmMat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = npar, ncol = ncurves)
  # varMat: variance-covariance matrix (must match der length = npar * ncurves)
  ncoef <- npar * ncurves
  set.seed(42)
  A <- matrix(rnorm(ncoef^2), ncoef, ncoef)
  varMat <- A %*% t(A) + diag(ncoef)  # positive definite

  list(
    indPair = c(1, 2),
    pVec = c(50, 50),
    compMatch = NULL,
    object = list(type = obj_type),
    indexMat = indexMat,
    parmMat = parmMat,
    varMat = varMat,
    level = level,
    reference = "control",
    type = "relative",
    sifct = make_sifct(val = 2.0, der = rep(0.1, ncoef)),
    interval = interval,
    degfree = degfree,
    logBase = logBase
  )
}

# -------------------------------------------------------------------
# Tests for siInner
# -------------------------------------------------------------------

test_that("siInner works with interval='none' and continuous data", {
  inputs <- make_test_inputs(interval = "none", obj_type = "continuous")
  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  # Result should be a numeric vector: c(siMatRow[1:4], dSIval)
  expect_true(is.numeric(result))
  ncoef <- nrow(inputs$indexMat) * ncol(inputs$indexMat)
  expect_length(result, 4 + ncoef)

  # First element is SIval = 2.0

  expect_equal(result[1], 2.0)

  # Second element is standard error
  expect_true(result[2] > 0)

  # Third element is t-statistic = (SI - 1) / SE
  expect_equal(result[3], (result[1] - 1) / result[2])

  # Fourth element is p-value (two-sided, t-distribution)
  tstat <- result[3]
  expected_p <- pt(-abs(tstat), inputs$degfree) + (1 - pt(abs(tstat), inputs$degfree))
  expect_equal(result[4], expected_p)
})

test_that("siInner works with interval='none' and non-continuous data (normal dist)", {
  inputs <- make_test_inputs(interval = "none", obj_type = "binomial")
  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  expect_true(is.numeric(result))
  expect_equal(result[1], 2.0)

  # For non-continuous, uses pnorm instead of pt
  tstat <- result[3]
  expected_p <- pnorm(-abs(tstat)) + (1 - pnorm(abs(tstat)))
  expect_equal(result[4], expected_p)
})

test_that("siInner works with interval='delta' and continuous data", {
  inputs <- make_test_inputs(interval = "delta", obj_type = "continuous")
  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  expect_true(is.numeric(result))
  expect_equal(result[1], 2.0)

  # Slots 2 and 3 are lower and upper CI bounds (delta method)
  # Lower bound < SIval < Upper bound
  expect_true(result[2] < result[1])
  expect_true(result[3] > result[1])
  # Slot 4 is NA (not set for delta)
  expect_true(is.na(result[4]))
})

test_that("siInner works with interval='delta' and non-continuous data", {
  inputs <- make_test_inputs(interval = "delta", obj_type = "binomial")
  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  expect_true(is.numeric(result))
  # Uses qnorm for non-continuous
  tquan <- qnorm(1 - (1 - inputs$level) / 2)
  expect_true(result[2] < result[1])
  expect_true(result[3] > result[1])
})

test_that("siInner works with interval='fls' without logBase", {
  inputs <- make_test_inputs(interval = "fls", obj_type = "continuous", logBase = NULL)
  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  expect_true(is.numeric(result))
  expect_equal(result[1], 2.0)
  # Lower and upper CI computed by delta method
  expect_true(result[2] < result[1])
  expect_true(result[3] > result[1])
})

test_that("siInner works with interval='fls' with logBase (from log scale)", {
  inputs <- make_test_inputs(interval = "fls", obj_type = "continuous", logBase = 10)
  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  expect_true(is.numeric(result))
  # With logBase=10, siMatRow values are transformed: logBase^(original)
  # Original SIval = 2.0, so result[1] = 10^2 = 100
  expect_equal(result[1], 10^2)
  # Lower and upper bounds are also transformed
  expect_true(result[2] > 0)
  expect_true(result[3] > 0)
})

test_that("siInner works with interval='tfls' and continuous data", {
  inputs <- make_test_inputs(interval = "tfls", obj_type = "continuous")
  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  expect_true(is.numeric(result))
  expect_equal(result[1], 2.0)
  # Lower and upper are exp(log(SIval) +/- tquan * lsdVal)
  expect_true(result[2] > 0)
  expect_true(result[3] > 0)
  expect_true(result[2] < result[1])
  expect_true(result[3] > result[1])
})

test_that("siInner works with interval='tfls' and non-continuous data", {
  inputs <- make_test_inputs(interval = "tfls", obj_type = "binomial")
  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  expect_true(is.numeric(result))
  expect_equal(result[1], 2.0)
  expect_true(result[2] > 0)
  expect_true(result[3] > 0)
})

test_that("siInner works with interval='fieller' and continuous data", {
  npar <- 4
  ncurves <- 2
  ncoef <- npar * ncurves

  # For fieller, sifct must return der1, der2, valnum, valden
  sifct_fieller <- make_sifct(
    val = 2.0,
    der = rep(0.1, ncoef),
    der1 = rep(0.05, ncoef),
    der2 = rep(0.03, ncoef),
    valnum = 10.0,
    valden = 5.0
  )

  inputs <- make_test_inputs(interval = "fieller", obj_type = "continuous")
  inputs$sifct <- sifct_fieller

  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  expect_true(is.numeric(result))
  expect_equal(result[1], 2.0)
  # Fieller CI: slots 2 and 3
  expect_false(is.na(result[2]))
  expect_false(is.na(result[3]))
})

test_that("siInner works with interval='fieller' and non-continuous data", {
  npar <- 4
  ncurves <- 2
  ncoef <- npar * ncurves

  sifct_fieller <- make_sifct(
    val = 2.0,
    der = rep(0.1, ncoef),
    der1 = rep(0.05, ncoef),
    der2 = rep(0.03, ncoef),
    valnum = 10.0,
    valden = 5.0
  )

  inputs <- make_test_inputs(interval = "fieller", obj_type = "binomial")
  inputs$sifct <- sifct_fieller

  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  expect_true(is.numeric(result))
  expect_equal(result[1], 2.0)
  # Fieller CI: slots 2 and 3
  expect_false(is.na(result[2]))
  expect_false(is.na(result[3]))
})

test_that("siInner return vector has correct structure", {
  inputs <- make_test_inputs(interval = "none", obj_type = "continuous")
  result <- drc:::siInner(
    indPair = inputs$indPair, pVec = inputs$pVec,
    compMatch = inputs$compMatch, object = inputs$object,
    indexMat = inputs$indexMat, parmMat = inputs$parmMat,
    varMat = inputs$varMat, level = inputs$level,
    reference = inputs$reference, type = inputs$type,
    sifct = inputs$sifct, interval = inputs$interval,
    degfree = inputs$degfree, logBase = inputs$logBase
  )
  ncoef <- nrow(inputs$indexMat) * ncol(inputs$indexMat)
  # siMatRow (4 elements) + dSIval (ncoef elements)
  expect_length(result, 4 + ncoef)
})

# -------------------------------------------------------------------
# Tests for fieller (helper in EDcomp.R)
# -------------------------------------------------------------------

test_that("fieller computes standard Fieller CI (finney=FALSE)", {
  mu <- c(10, 5)    # numerator and denominator
  df <- 20
  vcMat <- matrix(c(1, 0.2, 0.2, 0.5), 2, 2)

  result <- drc:::fieller(mu, df, vcMat, level = 0.95)
  expect_true(is.numeric(result))
  expect_length(result, 2)
  # Lower < ratio < Upper
  expect_true(result[1] < mu[1] / mu[2])
  expect_true(result[2] > mu[1] / mu[2])
})

test_that("fieller computes Finney variant (finney=TRUE)", {
  mu <- c(10, 5)
  df <- 20
  vcMat <- matrix(c(1, 0.2, 0.2, 0.5), 2, 2)
  resVar <- 2.0

  result <- drc:::fieller(mu, df, vcMat, level = 0.95, finney = TRUE, resVar = resVar)
  expect_true(is.numeric(result))
  expect_length(result, 2)
  expect_true(result[1] < result[2])
})

test_that("fieller with finney=TRUE throws error when g >= 1", {
  mu <- c(10, 0.5)  # small denominator
  df <- 20
  # Large vcMat[2,2] relative to mu[2]^2 so g >= 1
  vcMat <- matrix(c(1, 0.2, 0.2, 50), 2, 2)
  resVar <- 2.0

  expect_error(
    drc:::fieller(mu, df, vcMat, level = 0.95, finney = TRUE, resVar = resVar),
    "Fieller's theorem not useful"
  )
})

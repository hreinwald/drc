# tests/testthat/test-simFct.R
# Comprehensive tests for simFct() and coverFct() in R/simFct.R

# ==============================================================================
# Helper: access internal functions
# ==============================================================================
simFct  <- drc:::simFct
coverFct <- drc:::coverFct

# ==============================================================================
# Tests for simFct()
# ==============================================================================

# --- Parametric + Binomial + method "p" ---
test_that("simFct parametric binomial method='p' returns correct structure", {
  # Use deguelin dataset for binomial dose-response
  data(deguelin, package = "drc")
  m1 <- drm(r / n ~ dose, weights = n, data = deguelin, fct = LL.2(), type = "binomial")

  # Small simulation (2 sims) for speed
  res <- expect_output(
    simFct(
      noSim   = 2,
      edVal   = c(10, 50),
      type    = "parametric",
      response = "bin",
      fct     = LL.2(),
      coefVec = coef(m1),
      method  = "p",
      doseVec = deguelin$dose,
      nVec    = deguelin$n,
      pfct    = LL.2()
    )
  )

  expect_type(res, "list")
  expect_named(res, c("edArray", "mixVec", "edVal", "aicVec", "spanVec"))
  expect_equal(dim(res$edArray), c(2, 3, 2))
  expect_equal(length(res$mixVec), 2)
  expect_equal(res$edVal, c(10, 50))
  expect_equal(length(res$aicVec), 2)
  expect_equal(length(res$spanVec), 2)
  # For method "p", mixVec should be 0 (purely parametric)
  expect_true(all(res$mixVec == 0, na.rm = TRUE))
})

# --- Parametric + Continuous + method "p" ---
test_that("simFct parametric continuous method='p' returns correct structure", {
  data(ryegrass, package = "drc")
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  res <- expect_output(
    simFct(
      noSim    = 2,
      edVal    = c(50),
      type     = "parametric",
      response = "con",
      fct      = LL.4(),
      coefVec  = coef(m1),
      method   = "p",
      doseVec  = sort(unique(ryegrass$conc)),
      resVar   = summary(m1)$resVar,
      pfct     = LL.4()
    )
  )

  expect_type(res, "list")
  expect_equal(dim(res$edArray)[1], 1)
  expect_equal(dim(res$edArray)[3], 2)
  expect_equal(res$edVal, 50)
})

# --- Non-parametric + Binomial + method "p" ---
test_that("simFct non-parametric binomial method='p' returns correct structure", {
  doseVec <- c(0, 1, 2, 5, 10, 20, 50)
  nVec    <- rep(20, 7)
  pVec    <- c(0.01, 0.05, 0.15, 0.4, 0.7, 0.9, 0.99)

  res <- expect_output(
    simFct(
      noSim    = 2,
      edVal    = c(50),
      type     = "non-parametric",
      response = "bin",
      method   = "p",
      doseVec  = doseVec,
      nVec     = nVec,
      pVec     = pVec,
      pfct     = LL.2()
    )
  )

  expect_type(res, "list")
  expect_equal(dim(res$edArray), c(1, 3, 2))
})

# --- Non-parametric + Continuous + method "p" ---
test_that("simFct non-parametric continuous method='p' returns correct structure", {
  doseVec <- c(0, 1, 2, 5, 10, 20, 50)
  pVec    <- c(7, 6.5, 5.5, 4, 2.5, 1.5, 0.5)

  res <- expect_output(
    simFct(
      noSim    = 2,
      edVal    = c(50),
      type     = "non-parametric",
      response = "con",
      method   = "p",
      doseVec  = doseVec,
      pVec     = pVec,
      resVar   = 0.5,
      pfct     = LL.4()
    )
  )

  expect_type(res, "list")
  expect_equal(dim(res$edArray), c(1, 3, 2))
})

# --- Semi-parametric method "sp" with span provided ---
test_that("simFct parametric continuous method='sp' with span works", {
  data(ryegrass, package = "drc")
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # method "sp" requires semi-parametric (loess + parametric)
  # Using a fixed span to avoid GCV
  res <- expect_output(
    simFct(
      noSim    = 2,
      edVal    = c(50),
      type     = "parametric",
      response = "con",
      fct      = LL.4(),
      coefVec  = coef(m1),
      method   = "sp",
      doseVec  = sort(unique(ryegrass$conc)),
      resVar   = summary(m1)$resVar,
      pfct     = LL.4(),
      span     = 0.75
    )
  )

  expect_type(res, "list")
  expect_equal(dim(res$edArray), c(1, 3, 2))
  # spanVec should all be 0.75
  expect_equal(res$spanVec, rep(0.75, 2))
})

# --- Semi-parametric method "sp" with span=NA triggers GCV ---
test_that("simFct parametric continuous method='sp' with span=NA triggers GCV", {
  data(ryegrass, package = "drc")
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  res <- expect_output(
    simFct(
      noSim    = 2,
      edVal    = c(50),
      type     = "parametric",
      response = "con",
      fct      = LL.4(),
      coefVec  = coef(m1),
      method   = "sp",
      doseVec  = sort(unique(ryegrass$conc)),
      resVar   = summary(m1)$resVar,
      pfct     = LL.4(),
      span     = NA
    )
  )

  expect_type(res, "list")
  # spanVec should NOT be NA (GCV should have found values)
  expect_true(all(!is.na(res$spanVec)))
})

# --- Semi-parametric method "sp" with binomial ---
test_that("simFct parametric binomial method='sp' works", {
  data(deguelin, package = "drc")
  m1 <- drm(r / n ~ dose, weights = n, data = deguelin, fct = LL.2(), type = "binomial")

  res <- expect_output(
    simFct(
      noSim   = 2,
      edVal   = c(50),
      type    = "parametric",
      response = "bin",
      fct     = LL.2(),
      coefVec = coef(m1),
      method  = "sp",
      doseVec = deguelin$dose,
      nVec    = deguelin$n,
      pfct    = LL.2(),
      span    = 0.75
    )
  )

  expect_type(res, "list")
})

# --- Method "p" with drm failure (try-error path) ---
test_that("simFct handles drm failure in method='p' gracefully", {
  # Use dose values and parameters that will cause drm to fail
  # Very few points and extreme values
  doseVec <- c(0, 100)
  pVec    <- c(0.5, 0.5)  # flat response - hard to fit
  nVec    <- rep(2, 2)

  res <- expect_output(
    simFct(
      noSim    = 2,
      edVal    = c(50),
      type     = "non-parametric",
      response = "bin",
      method   = "p",
      doseVec  = doseVec,
      nVec     = nVec,
      pVec     = pVec,
      pfct     = LL.2()
    )
  )

  expect_type(res, "list")
  # Some or all should be NA due to fitting failures
  expect_equal(dim(res$edArray), c(1, 3, 2))
})

# --- Multiple ED values ---
test_that("simFct works with multiple ED values", {
  data(ryegrass, package = "drc")
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  res <- expect_output(
    simFct(
      noSim    = 2,
      edVal    = c(10, 20, 50, 90),
      type     = "parametric",
      response = "con",
      fct      = LL.4(),
      coefVec  = coef(m1),
      method   = "p",
      doseVec  = sort(unique(ryegrass$conc)),
      resVar   = summary(m1)$resVar,
      pfct     = LL.4()
    )
  )

  expect_equal(dim(res$edArray)[1], 4)
  expect_equal(res$edVal, c(10, 20, 50, 90))
})

# --- match.arg validation ---
test_that("simFct validates method argument", {
  expect_error(
    simFct(
      noSim = 1, method = "invalid",
      doseVec = 1:5, type = "parametric",
      response = "con"
    ),
    "arg"
  )
})

test_that("simFct validates response argument", {
  expect_error(
    simFct(
      noSim = 1, method = "p",
      doseVec = 1:5, type = "parametric",
      response = "invalid"
    ),
    "arg"
  )
})

test_that("simFct validates type argument", {
  expect_error(
    simFct(
      noSim = 1, method = "p",
      doseVec = 1:5, type = "invalid",
      response = "con"
    ),
    "arg"
  )
})

# --- mrdrm try-error path (line 133) ---
test_that("simFct handles mrdrm try-error in method='sp'", {
  data(ryegrass, package = "drc")
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  local_mocked_bindings(
    mrdrm = function(...) structure("mock mrdrm error", class = "try-error"),
    .package = "drc"
  )

  res <- expect_output(
    simFct(
      noSim    = 1,
      edVal    = c(50),
      type     = "parametric",
      response = "con",
      fct      = LL.4(),
      coefVec  = coef(m1),
      method   = "sp",
      doseVec  = sort(unique(ryegrass$conc)),
      resVar   = summary(m1)$resVar,
      pfct     = LL.4(),
      span     = 0.75
    )
  )

  expect_true(all(is.na(res$edArray)))
})

# --- ED try-error path in method "p" (lines 158-159) ---
test_that("simFct handles ED try-error in method='p'", {
  data(ryegrass, package = "drc")
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  local_mocked_bindings(
    ED.drc = function(...) stop("mock ED failure"),
    .package = "drc"
  )

  res <- expect_output(
    simFct(
      noSim    = 1,
      edVal    = c(50),
      type     = "parametric",
      response = "con",
      fct      = LL.4(),
      coefVec  = coef(m1),
      method   = "p",
      doseVec  = sort(unique(ryegrass$conc)),
      resVar   = summary(m1)$resVar,
      pfct     = LL.4()
    )
  )

  expect_true(all(is.na(res$edArray)))
  expect_true(all(is.na(res$mixVec)))
})

# ==============================================================================
# Tests for coverFct()
# ==============================================================================

test_that("coverFct with provided edVec returns correct structure", {
  # Create a mock simulation result
  # 2 ED levels, 3 columns (estimate, lower, upper), 10 sims
  edArray <- array(NA, c(2, 3, 10))
  # Fill with plausible values
  set.seed(42)
  for (i in 1:10) {
    edArray[1, 1, i] <- rnorm(1, 5, 0.5)    # ED10 estimate ~5
    edArray[1, 2, i] <- edArray[1, 1, i] - 1  # lower bound
    edArray[1, 3, i] <- edArray[1, 1, i] + 1  # upper bound
    edArray[2, 1, i] <- rnorm(1, 10, 1)       # ED50 estimate ~10
    edArray[2, 2, i] <- edArray[2, 1, i] - 2
    edArray[2, 3, i] <- edArray[2, 1, i] + 2
  }

  simres <- list(
    edArray = edArray,
    mixVec  = rep(0, 10),
    edVal   = c(10, 50),
    aicVec  = rep(-10, 10),
    spanVec = rep(0.75, 10)
  )

  res <- coverFct(mfit = NULL, simres = simres, edVec = c(5, 10))

  expect_type(res, "list")
  expect_named(res, c("coverage", "covLow", "covUp", "true", "mean", "width", "notNAs", "NAs", "mixingAverage"))
  expect_equal(length(res$coverage), 2)
  expect_equal(res$true, c(5, 10))
  expect_true(is.numeric(res$coverage))
  expect_true(all(res$coverage >= 0 & res$coverage <= 1))
  expect_equal(res$mixingAverage, 0)
  expect_equal(length(res$notNAs), 2)
  expect_equal(length(res$NAs), 2)
})

test_that("coverFct with NULL edVec computes ED from model", {
  data(ryegrass, package = "drc")
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Create simulation result matching m1
  edArray <- array(NA, c(1, 3, 5))
  trueED <- ED(m1, 50, display = FALSE)[, 1]
  for (i in 1:5) {
    edArray[1, 1, i] <- trueED + rnorm(1, 0, 0.1)
    edArray[1, 2, i] <- edArray[1, 1, i] - 0.5
    edArray[1, 3, i] <- edArray[1, 1, i] + 0.5
  }

  simres <- list(
    edArray = edArray,
    mixVec  = rep(0, 5),
    edVal   = c(50),
    aicVec  = rep(-10, 5),
    spanVec = rep(NA, 5)
  )

  res <- coverFct(mfit = m1, simres = simres, edVec = NULL)

  expect_type(res, "list")
  expect_equal(length(res$coverage), 1)
  expect_true(!is.na(res$true))
})

test_that("coverFct handles NA values in edArray", {
  # Create simulation result with some NAs
  edArray <- array(NA, c(1, 3, 10))
  for (i in 1:7) {
    edArray[1, 1, i] <- rnorm(1, 5, 0.5)
    edArray[1, 2, i] <- edArray[1, 1, i] - 1
    edArray[1, 3, i] <- edArray[1, 1, i] + 1
  }
  # Leave i=8,9,10 as NA

  simres <- list(
    edArray = edArray,
    mixVec  = c(rep(0, 7), rep(NA, 3)),
    edVal   = c(50),
    aicVec  = rep(NA, 10),
    spanVec = rep(NA, 10)
  )

  res <- coverFct(mfit = NULL, simres = simres, edVec = c(5))

  expect_type(res, "list")
  # notNAs should be 7 (the ones that have both lower and upper)
  expect_equal(res$notNAs, 7)
  expect_equal(res$NAs, 3)
})

test_that("coverFct handles partial NAs (only lower or only upper is NA)", {
  # Create simulation result where some have only lower NA, some only upper NA
  edArray <- array(NA, c(1, 3, 6))
  # Normal cases
  for (i in 1:2) {
    edArray[1, 1, i] <- 5
    edArray[1, 2, i] <- 4
    edArray[1, 3, i] <- 6
  }
  # Lower NA, upper present and > true value
  edArray[1, 1, 3] <- 5
  edArray[1, 2, 3] <- NA
  edArray[1, 3, 3] <- 6  # upper > 5

  # Lower NA, upper present but < true value
  edArray[1, 1, 4] <- 5
  edArray[1, 2, 4] <- NA
  edArray[1, 3, 4] <- 4  # upper < 5

  # Upper NA, lower present and < true value
  edArray[1, 1, 5] <- 5
  edArray[1, 2, 5] <- 4  # lower < 5
  edArray[1, 3, 5] <- NA

  # Upper NA, lower present but > true value
  edArray[1, 1, 6] <- 5
  edArray[1, 2, 6] <- 6  # lower > 5
  edArray[1, 3, 6] <- NA

  simres <- list(
    edArray = edArray,
    mixVec  = rep(0, 6),
    edVal   = c(50),
    aicVec  = rep(NA, 6),
    spanVec = rep(NA, 6)
  )

  res <- coverFct(mfit = NULL, simres = simres, edVec = c(5))

  expect_type(res, "list")
  # Only observations with both lower and upper non-NA count for notNAs
  expect_equal(res$notNAs, 2)
  expect_equal(res$NAs, 4)
  # covLow counts: lower NA AND upper > true. Obs 3 qualifies (NA lower, upper=6>5)
  expect_equal(res$covLow, 1)
  # covUp counts: upper NA AND lower < true. Obs 5 qualifies (NA upper, lower=4<5)
  expect_equal(res$covUp, 1)
})

# ==============================================================================
# Integration test: simFct + coverFct together
# ==============================================================================

test_that("simFct and coverFct work together end-to-end", {
  data(deguelin, package = "drc")
  m1 <- drm(r / n ~ dose, weights = n, data = deguelin, fct = LL.2(), type = "binomial")

  simres <- expect_output(
    simFct(
      noSim   = 3,
      edVal   = c(50),
      type    = "parametric",
      response = "bin",
      fct     = LL.2(),
      coefVec = coef(m1),
      method  = "p",
      doseVec = deguelin$dose,
      nVec    = deguelin$n,
      pfct    = LL.2()
    )
  )

  covRes <- coverFct(mfit = m1, simres = simres)

  expect_type(covRes, "list")
  expect_equal(length(covRes$coverage), 1)
  expect_true(is.numeric(covRes$coverage))
})

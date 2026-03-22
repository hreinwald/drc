# Test mselect() function

# Create test dataset (ryegrass)
ryegrass <- data.frame(
  rootl = c(
    7.58, 8.00, 8.33, 7.25, 7.17, 7.00, 7.17, 7.83, 7.92, 7.58,
    6.17, 5.75, 5.83, 6.00, 5.83, 4.92, 4.50, 4.17, 4.42, 4.00,
    2.67, 2.08, 2.42, 2.50, 2.25, 1.17, 0.75, 0.92, 1.00, 0.58
  ),
  conc = c(
    rep(0, 5), rep(0.94, 5), rep(1.88, 5),
    rep(3.75, 5), rep(7.50, 5), rep(15, 5)
  )
)

test_that("mselect returns non-zero Lack of fit p-values when nested = FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fctList <- list(LL.3(), W1.4(), W2.4())

  result <- mselect(m1, fctList = fctList, nested = FALSE)

  # "Lack of fit" column should contain meaningful p-values, not all zeros
  lof_col <- result[, "Lack of fit"]
  expect_true(any(lof_col > 0, na.rm = TRUE),
    info = "Lack of fit p-values should not all be zero when nested = FALSE")
})

test_that("mselect Lack of fit values match between nested = TRUE and nested = FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fctList <- list(LL.3(), W1.4(), W2.4())

  result_nested <- mselect(m1, fctList = fctList, nested = TRUE)
  result_plain <- mselect(m1, fctList = fctList, nested = FALSE)

  # The "Lack of fit" column should have the same values regardless of 'nested'
  expect_equal(
    result_nested[, "Lack of fit"],
    result_plain[, "Lack of fit"],
    info = "Lack of fit p-values should be the same whether nested is TRUE or FALSE"
  )
})

test_that("mselect with nested = TRUE does not produce NaN in Nested F test column", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fctList <- list(LL.3(), W1.4(), W2.4())

  # Suppress benign warnings from the optimization process (log of negative values
  # during parameter search), which are unrelated to the F-test computation
  result <- suppressWarnings(mselect(m1, fctList = fctList, nested = TRUE))

  # The "Nested F test" column should have valid p-values (or NA for first model),
  # but never NaN
  ftest_col <- result[, "Nested F test"]
  expect_true(all(is.na(ftest_col) | is.finite(ftest_col)),
    info = "Nested F test column should not contain NaN values")
})

test_that("anova.drc returns p-value of 1 when F statistic is negative", {
  # Comparing non-nested models with same df can produce negative F statistics
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m2 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())

  # Should not produce NaN warnings from pf() and p-value should be valid
  result <- anova(m1, m2, details = FALSE)
  pval <- result[2, 5]
  fstat <- result[2, 4]
  expect_true(!is.nan(pval), info = "p-value should not be NaN")
  # If the F statistic was negative, p-value should be 1
  if (!is.na(fstat) && fstat < 0) {
    expect_equal(pval, 1, info = "Negative F statistic should give p-value of 1")
  }
  # If the F statistic is a valid positive number, p-value should be between 0 and 1
  if (!is.na(fstat) && is.finite(fstat) && fstat >= 0) {
    expect_true(pval >= 0 && pval <= 1,
      info = "Positive F statistic should give p-value between 0 and 1")
  }
})

test_that("anova.drc returns valid p-values for truly nested models", {
  # LL.3 is nested within LL.4 (LL.3 fixes one parameter)
  m_full <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m_reduced <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())

  result <- anova(m_reduced, m_full, details = FALSE)
  pval <- result[2, 5]
  fstat <- result[2, 4]

  # For truly nested models, F statistic should be non-negative and p-value valid
  expect_true(!is.na(pval), info = "p-value should not be NA for nested models")
  expect_true(!is.nan(pval), info = "p-value should not be NaN for nested models")
  expect_true(pval >= 0 && pval <= 1,
    info = "p-value should be between 0 and 1 for nested models")
})

test_that("mselect with nested = FALSE returns correct column names", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fctList <- list(LL.3())

  result <- mselect(m1, fctList = fctList, nested = FALSE)
  expect_true("Lack of fit" %in% colnames(result))
  expect_true("logLik" %in% colnames(result))
  expect_true("IC" %in% colnames(result))
  expect_true("Res var" %in% colnames(result))
  expect_false("Nested F test" %in% colnames(result))
})

test_that("mselect with nested = TRUE returns Nested F test column", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fctList <- list(LL.3())

  result <- mselect(m1, fctList = fctList, nested = TRUE)
  expect_true("Nested F test" %in% colnames(result))
})

test_that("mselect Res var column contains numeric values, not try-error objects", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fctList <- list(LL.3(), W1.4())

  result <- mselect(m1, fctList = fctList, nested = FALSE)

  # The "Res var" column should contain numeric values, never try-error objects
  # This tests the fix for inherits("tryRV", "try-error") bug where the string
  # literal was checked instead of the variable
  resvar_col <- result[, "Res var"]
  expect_true(is.numeric(resvar_col),
    info = "Res var column should be numeric")
  expect_true(all(resvar_col > 0, na.rm = TRUE),
    info = "Res var values should be positive for continuous data")
})

test_that("modelFit does not produce NaN p-values", {
  # Test that modelFit returns valid p-values even in edge cases
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  result <- modelFit(m1)
  pval <- result[2, 5]

  # p-value should not be NaN
  expect_true(!is.nan(pval),
    info = "modelFit p-value should not be NaN")
  # p-value should be between 0 and 1 (or NA)
  if (!is.na(pval)) {
    expect_true(pval >= 0 && pval <= 1,
      info = "modelFit p-value should be between 0 and 1")
  }
})

# --- Tests for uncovered code paths ---

test_that("mselect errors when nested is not logical", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(
    mselect(m1, nested = "yes"),
    "'nested' argument takes only the values: FALSE, TRUE"
  )
  expect_error(
    mselect(m1, nested = 1),
    "'nested' argument takes only the values: FALSE, TRUE"
  )
})

test_that("mselect with linreg = TRUE includes polynomial fits", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  result <- mselect(m1, fctList = list(LL.3()), linreg = TRUE)

  # Should have rows for Lin, Quad, Cubic in addition to DRC models
  expect_true("Lin" %in% rownames(result))
  expect_true("Quad" %in% rownames(result))
  expect_true("Cubic" %in% rownames(result))

  # Result should be a matrix with correct columns
  expect_true(is.matrix(result))
  expect_true("logLik" %in% colnames(result))
  expect_true("IC" %in% colnames(result))
  expect_true("Lack of fit" %in% colnames(result))
  expect_true("Res var" %in% colnames(result))

  # Polynomial fit rows should have NA for "Lack of fit"
  expect_true(is.na(result["Lin", "Lack of fit"]))
  expect_true(is.na(result["Quad", "Lack of fit"]))
  expect_true(is.na(result["Cubic", "Lack of fit"]))

  # logLik, IC, and Res var should be numeric for polynomial fits
  expect_true(is.finite(result["Lin", "logLik"]))
  expect_true(is.finite(result["Quad", "IC"]))
  expect_true(is.finite(result["Cubic", "Res var"]))
})

test_that("mselect with linreg = TRUE and nested = TRUE drops Nested F test column", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  result <- mselect(m1, fctList = list(LL.3()), linreg = TRUE, nested = TRUE)

  # When linreg = TRUE, the Nested F test column should be removed
  expect_false("Nested F test" %in% colnames(result))
  # Should still have the first 4 columns
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("logLik", "IC", "Lack of fit", "Res var"))
})

test_that("mselect handles model fitting failure gracefully (line 109)", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Create a dose-response function whose starting values always error,
  # guaranteeing update() returns a try-error
  bad_fct <- LL.3()
  bad_fct$name <- "BadModel"
  bad_fct$ssfct <- function(...) stop("Cannot compute starting values")

  result <- suppressWarnings(
    mselect(m1, fctList = list(bad_fct), nested = FALSE)
  )

  expect_true(is.matrix(result))
  # The first model (LL.4) should always have valid values
  expect_true(is.finite(result["LL.4", "logLik"]))
  # The bad model should have all NA values
  expect_true(all(is.na(result["BadModel", ])))
})

test_that("mselect handles summary failure for initial model resVar (line 71)", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Corrupt the initial model so summary.drc fails on it:
  # Setting robust to "Tukey's biweight" triggers a code path in summary.drc
  # that calls solve(hessian); a zero hessian causes solve() to error
  m1$robust <- "Tukey's biweight"
  m1$fit$hessian <- matrix(0, 4, 4)

  result <- mselect(m1, nested = FALSE, sorted = "no")

  expect_true(is.matrix(result))
  # The initial model's Res var should be NA because summary failed
  expect_true(is.na(result[1, "Res var"]))
})

test_that("mselect handles summary failure for updated model resVar (line 100)", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Use namespace-level mock to make summary.drc fail only on the second call
  # (first call is for the initial model, second is for the updated model)
  orig_summary <- drc:::summary.drc
  call_count <- 0L
  mock_summary <- function(object, ...) {
    call_count <<- call_count + 1L
    if (call_count > 1L) {
      stop("Simulated summary failure for updated model")
    }
    orig_summary(object, ...)
  }

  assignInNamespace("summary.drc", mock_summary, ns = "drc")
  on.exit(assignInNamespace("summary.drc", orig_summary, ns = "drc"), add = TRUE)

  result <- mselect(m1, fctList = list(LL.3()), nested = FALSE, sorted = "no")

  expect_true(is.matrix(result))
  # The initial model's Res var should be valid (first summary call succeeds)
  expect_true(is.finite(result["LL.4", "Res var"]))
  # The updated model's Res var should be NA (second summary call fails)
  expect_true(is.na(result["LL.3", "Res var"]))
})

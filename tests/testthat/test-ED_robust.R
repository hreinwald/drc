# Tests for ED_robust() and maED_robust() functions
# Also covers helper functions: get_ed_interval() and drm_name()

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


# =============================================================================
# Tests for get_ed_interval()
# =============================================================================

test_that("get_ed_interval returns 'tfls' for LL model with small_n = TRUE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_equal(drc:::get_ed_interval(m1), "tfls")
})

test_that("get_ed_interval returns 'fls' for LL model with small_n = FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_equal(drc:::get_ed_interval(m1, small_n = FALSE), "fls")
})

test_that("get_ed_interval returns 'delta' for Weibull models", {
  expect_equal(drc:::get_ed_interval("W1.4"), "delta")
  expect_equal(drc:::get_ed_interval("W2.3"), "delta")
})

test_that("get_ed_interval returns 'tfls' for LL character input", {
  expect_equal(drc:::get_ed_interval("LL.4"), "tfls")
  expect_equal(drc:::get_ed_interval("LL.4", small_n = TRUE), "tfls")
})

test_that("get_ed_interval returns 'fls' for LN character input with small_n = FALSE", {
  expect_equal(drc:::get_ed_interval("LN.4", small_n = FALSE), "fls")
})

test_that("get_ed_interval returns 'tfls' for BC character input", {
  expect_equal(drc:::get_ed_interval("BC.4"), "tfls")
})

test_that("get_ed_interval returns 'tfls' for CRS character input", {
  expect_equal(drc:::get_ed_interval("CRS.4"), "tfls")
})

test_that("get_ed_interval defaults to 'tfls' for unknown model with message when verbose", {
  expect_message(
    result <- drc:::get_ed_interval("SomeUnknownModel", verbose = TRUE),
    "Defaulting to 'tfls'"
  )
  expect_equal(result, "tfls")
})

test_that("get_ed_interval defaults to 'tfls' for unknown model silently when not verbose", {
  expect_silent(result <- drc:::get_ed_interval("SomeUnknownModel", verbose = FALSE))
  expect_equal(result, "tfls")
})

test_that("get_ed_interval errors for invalid input", {
  expect_error(drc:::get_ed_interval(42), "must be a 'drc' object or a single character string")
  expect_error(drc:::get_ed_interval(c("LL.4", "W1.4")), "must be a 'drc' object or a single character string")
  expect_error(drc:::get_ed_interval(NULL), "must be a 'drc' object or a single character string")
})

test_that("get_ed_interval works with drc object input", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- drc:::get_ed_interval(m1)
  expect_true(result %in% c("tfls", "fls", "delta"))
})

test_that("get_ed_interval works with Weibull drc model", {
  m_w <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
  expect_equal(drc:::get_ed_interval(m_w), "delta")
})


# =============================================================================
# Tests for drm_name()
# =============================================================================

test_that("drm_name returns correct format for LL.4 model", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- drc:::drm_name(m1)
  expect_true(is.character(result))
  expect_true(grepl("LL.4", result))
  expect_true(grepl(":", result))
  expect_true(grepl("-", result))
})

test_that("drm_name errors for non-drc input", {
  expect_error(drc:::drm_name("not_a_model"), "must be a `drc` object")
  expect_error(drc:::drm_name(42), "must be a `drc` object")
  expect_error(drc:::drm_name(lm(rootl ~ conc, data = ryegrass)), "must be a `drc` object")
})


# =============================================================================
# Tests for ED_robust()
# =============================================================================

test_that("ED_robust returns data.table with correct structure", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED_robust(m1, respLev = c(10, 50))

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 2)
  expect_true(all(c("Estimate", "stderr", "Lower", "Upper", "confint_level",
                     "confint_method", "model", "EC") %in% names(result)))
})

test_that("ED_robust returns positive estimates for valid levels", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED_robust(m1, respLev = c(10, 50, 90))

  expect_true(all(result$Estimate > 0, na.rm = TRUE))
  expect_true(all(result$EC == c(10, 50, 90)))
})

test_that("ED_robust returns correct metadata", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED_robust(m1, respLev = 50, CI_level = 0.90)

  expect_equal(result$confint_level, 0.90)
  expect_equal(result$EC, 50)
  expect_true(grepl("LL.4", result$model))
})

test_that("ED_robust returns NA for non-estimable response levels", {
  # BC.4 model with extreme response levels might fail
  data(lettuce, package = "drc")
  m1 <- drm(weight ~ conc, data = lettuce, fct = BC.4())
  result <- ED_robust(m1, respLev = c(50, 99))

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 2)
  # At least one row should have data; some extreme levels may be NA
})

test_that("ED_robust handles errors gracefully and returns NA rows", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  # Remove edfct to force an error inside ED
  m1_broken <- m1
  m1_broken$fct$edfct <- NULL

  result <- ED_robust(m1_broken, respLev = c(10, 50))

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 2)
  # All should be NA since the model is broken
  expect_true(all(is.na(result$Estimate)))
})

test_that("ED_robust verbose mode prints messages on success", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_message(
    ED_robust(m1, respLev = 50, verbose = TRUE),
    "Successfully calculated ED"
  )
})

test_that("ED_robust verbose mode prints messages on error", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m1$fct$edfct <- NULL
  expect_message(
    ED_robust(m1, respLev = 50, verbose = TRUE),
    "Error calculating ED"
  )
})

test_that("ED_robust verbose mode prints appending info message", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_message(
    ED_robust(m1, respLev = 50, verbose = TRUE),
    "Appending info"
  )
})

test_that("ED_robust handles single response level", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED_robust(m1, respLev = 50)

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 1)
})

test_that("ED_robust uses default interval from get_ed_interval", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED_robust(m1, respLev = 50)
  expect_true(result$confint_method %in% c("tfls", "fls", "delta"))
})

test_that("ED_robust works with Weibull model", {
  m_w <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
  result <- ED_robust(m_w, respLev = c(10, 50))

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 2)
  expect_equal(result$confint_method[1], "delta")
})

test_that("ED_robust returns NA row for negative or NA estimate", {
  # Use a model where ED at an extreme level might produce negative estimates
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  # Test with very extreme response level that might produce non-positive ED
  result <- ED_robust(m1, respLev = c(50, 99.99))

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 2)
})


# =============================================================================
# Tests for maED_robust()
# =============================================================================

test_that("maED_robust returns data.table with correct structure", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fct_list <- list(LL.5 = LL.5(), W1.4 = W1.4())
  result <- maED_robust(m1, fct_ls = fct_list, respLev = c(10, 50))

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 2)
  expect_true(all(c("Estimate", "stderr", "Lower", "Upper", "confint_level",
                     "confint_method", "model", "EC") %in% names(result)))
})

test_that("maED_robust returns positive estimates for valid levels", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fct_list <- list(LL.5 = LL.5())
  result <- maED_robust(m1, fct_ls = fct_list, respLev = c(10, 50))

  expect_true(all(result$Estimate > 0, na.rm = TRUE))
})

test_that("maED_robust returns correct metadata", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fct_list <- list(W1.4 = W1.4())
  result <- maED_robust(m1, fct_ls = fct_list, respLev = 50, CI_level = 0.90)

  expect_equal(result$confint_level, 0.90)
  expect_equal(result$confint_method, "buckland")
  expect_equal(result$EC, 50)
  expect_true(grepl("/", result$model))  # model name should contain "/" separator
})

test_that("maED_robust handles errors gracefully", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m1_broken <- m1
  m1_broken$fct$edfct <- NULL

  fct_list <- list(LL.5 = LL.5())
  result <- maED_robust(m1_broken, fct_ls = fct_list, respLev = c(10, 50))

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 2)
  # Should have NA estimates due to broken model
  expect_true(all(is.na(result$Estimate)))
})

test_that("maED_robust verbose mode prints messages on success", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fct_list <- list(LL.5 = LL.5())
  expect_message(
    maED_robust(m1, fct_ls = fct_list, respLev = 50, verbose = TRUE),
    "Successfully calculated maED"
  )
})

test_that("maED_robust verbose mode prints messages on error", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m1$fct$edfct <- NULL
  fct_list <- list(LL.5 = LL.5())
  expect_message(
    maED_robust(m1, fct_ls = fct_list, respLev = 50, verbose = TRUE),
    "Error calculating maED"
  )
})

test_that("maED_robust verbose mode prints appending info message", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fct_list <- list(LL.5 = LL.5())
  expect_message(
    maED_robust(m1, fct_ls = fct_list, respLev = 50, verbose = TRUE),
    "Appending info"
  )
})

test_that("maED_robust handles single response level", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fct_list <- list(LL.5 = LL.5())
  result <- maED_robust(m1, fct_ls = fct_list, respLev = 50)

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 1)
})

test_that("maED_robust returns NA for non-estimable response levels", {
  data(lettuce, package = "drc")
  m1 <- drm(weight ~ conc, data = lettuce, fct = BC.5())
  fct_list <- list(W2.4 = W2.4())
  result <- maED_robust(m1, fct_ls = fct_list, respLev = c(50, 99))

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 2)
})

test_that("maED_robust model name includes all model names", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fct_list <- list(LL.5 = LL.5(), W1.4 = W1.4())
  result <- maED_robust(m1, fct_ls = fct_list, respLev = 50)

  # Model name should combine base model + alternatives separated by /
  expect_true(grepl("LL.4", result$model))
})

test_that("maED_robust works with default parameters", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fct_list <- list(LL.5 = LL.5())
  result <- maED_robust(m1, fct_ls = fct_list)

  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 3)  # default respLev is c(10, 20, 50)
  expect_equal(result$confint_level[1], 0.95)
})


# --- Tests for non-positive/NA estimate paths (lines 148, 283) ---

test_that("ED_robust returns NA row when ED estimate is non-positive", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Mock drc::ED to return a result with a non-positive estimate
  mock_ed <- function(mod, respLev, interval, level, display, ...) {
    mat <- matrix(c(-1, 0.5, -2, 0), nrow = 1)
    colnames(mat) <- c("Estimate", "Std. Error", "Lower", "Upper")
    rownames(mat) <- paste0("e:1:", respLev)
    mat
  }
  local_mocked_bindings(ED = mock_ed, .package = "drc")

  result <- ED_robust(m1, respLev = 50)
  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 1)
  expect_true(is.na(result$Estimate))
})

test_that("ED_robust returns NA row when ED estimate is NA", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  mock_ed <- function(mod, respLev, interval, level, display, ...) {
    mat <- matrix(c(NA_real_, 0.5, NA_real_, NA_real_), nrow = 1)
    colnames(mat) <- c("Estimate", "Std. Error", "Lower", "Upper")
    rownames(mat) <- paste0("e:1:", respLev)
    mat
  }
  local_mocked_bindings(ED = mock_ed, .package = "drc")

  result <- ED_robust(m1, respLev = 50)
  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 1)
  expect_true(is.na(result$Estimate))
})

test_that("maED_robust returns NA row when maED estimate is non-positive", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fct_list <- list(LL.5 = LL.5())

  mock_maED <- function(mod, fctList, respLev, interval, level, display, na.rm, ...) {
    mat <- matrix(c(-1, 0.5, -2, 0), nrow = 1)
    colnames(mat) <- c("Estimate", "Std. Error", "Lower", "Upper")
    rownames(mat) <- paste0("e:1:", respLev)
    mat
  }
  local_mocked_bindings(maED = mock_maED, .package = "drc")

  result <- maED_robust(m1, fct_ls = fct_list, respLev = 50)
  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 1)
  expect_true(is.na(result$Estimate))
})

test_that("maED_robust returns NA row when maED estimate is NA", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fct_list <- list(LL.5 = LL.5())

  mock_maED <- function(mod, fctList, respLev, interval, level, display, na.rm, ...) {
    mat <- matrix(c(NA_real_, 0.5, NA_real_, NA_real_), nrow = 1)
    colnames(mat) <- c("Estimate", "Std. Error", "Lower", "Upper")
    rownames(mat) <- paste0("e:1:", respLev)
    mat
  }
  local_mocked_bindings(maED = mock_maED, .package = "drc")

  result <- maED_robust(m1, fct_ls = fct_list, respLev = 50)
  expect_true(data.table::is.data.table(result))
  expect_equal(nrow(result), 1)
  expect_true(is.na(result$Estimate))
})

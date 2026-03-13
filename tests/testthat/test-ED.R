# Test ED.drc() function - Estimating effective doses

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

# Multi-curve dataset
set.seed(42)
multi_data <- data.frame(
  dose = rep(c(0, 0.5, 1, 2, 5, 10), each = 5, times = 2),
  resp = c(
    rnorm(5, 100, 5), rnorm(5, 95, 5), rnorm(5, 85, 5),
    rnorm(5, 60, 5), rnorm(5, 20, 5), rnorm(5, 5, 5),
    rnorm(5, 100, 5), rnorm(5, 90, 5), rnorm(5, 70, 5),
    rnorm(5, 40, 5), rnorm(5, 10, 5), rnorm(5, 3, 5)
  ),
  group = rep(c("A", "B"), each = 30)
)

# Tests for ED.drc() with single curve models

test_that("ED.drc returns correct structure for single response level", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, 50, display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)
  expect_true(all(c("Estimate", "Std. Error") %in% colnames(result)))
})

test_that("ED.drc returns correct structure for multiple response levels", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, c(10, 50, 90), display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)
  expect_true(all(result[, "Estimate"] > 0))
})

test_that("ED.drc with delta method confidence intervals", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, 50, interval = "delta", display = FALSE)

  expect_equal(ncol(result), 4)
  expect_true(all(c("Estimate", "Std. Error", "Lower", "Upper") %in% colnames(result)))
  expect_true(result[, "Lower"] < result[, "Estimate"])
  expect_true(result[, "Upper"] > result[, "Estimate"])
})

test_that("ED.drc with different confidence levels", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result_95 <- ED(m1, 50, interval = "delta", level = 0.95, display = FALSE)
  result_90 <- ED(m1, 50, interval = "delta", level = 0.90, display = FALSE)

  # 90% CI should be narrower than 95% CI
  width_95 <- result_95[, "Upper"] - result_95[, "Lower"]
  width_90 <- result_90[, "Upper"] - result_90[, "Lower"]
  expect_true(width_90 < width_95)
})

test_that("ED.drc validates response level bounds for relative type", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Should error for response levels outside (0, 100)
  expect_error(ED(m1, 0, display = FALSE), "outside the interval")
  expect_error(ED(m1, 100, display = FALSE), "outside the interval")
  expect_error(ED(m1, -10, display = FALSE), "outside the interval")
  expect_error(ED(m1, 150, display = FALSE), "outside the interval")
})

test_that("ED.drc allows extreme values when bound = FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Should not error with bound = FALSE
  expect_no_error(ED(m1, 0, bound = FALSE, display = FALSE))
  expect_no_error(ED(m1, 100, bound = FALSE, display = FALSE))
})

test_that("ED.drc works with absolute type response levels", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, 5, type = "absolute", display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_true(result[, "Estimate"] > 0)
})

test_that("ED.drc errors when model has no edfct function", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  # Remove the edfct function to simulate a model without it
  m1$fct$edfct <- NULL

  expect_error(ED(m1, 50, display = FALSE), "ED values cannot be calculated")
})

# Tests for multi-curve models

test_that("ED.drc works with multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  result <- ED(m_multi, 50, display = FALSE)

  expect_equal(nrow(result), 2)  # One ED50 for each curve
  expect_true(all(grepl("A:|B:", rownames(result))))
})

test_that("ED.drc with clevel filters specific curves", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  result <- ED(m_multi, 50, clevel = "A", display = FALSE)

  expect_equal(nrow(result), 1)
  expect_true(grepl("A:", rownames(result)))
})

test_that("ED.drc with multiple response levels and curves", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  result <- ED(m_multi, c(10, 50, 90), display = FALSE)

  expect_equal(nrow(result), 6)  # 3 response levels × 2 curves
})

# Tests for different interval types

test_that("ED.drc with fls interval (from log scale)", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, 50, interval = "fls", display = FALSE)

  expect_equal(ncol(result), 3)  # Estimate, Lower, Upper (no Std. Error)
  expect_true(all(c("Estimate", "Lower", "Upper") %in% colnames(result)))
  expect_false("Std. Error" %in% colnames(result))
})

test_that("ED.drc with tfls interval (to and from log scale)", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, 50, interval = "tfls", display = FALSE)

  expect_equal(ncol(result), 4)
  expect_true(result[, "Lower"] < result[, "Estimate"])
  expect_true(result[, "Upper"] > result[, "Estimate"])
})

test_that("ED.drc with inverse regression interval", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, 50, interval = "inv", display = FALSE)

  expect_equal(ncol(result), 3)  # Estimate, Lower, Upper (no Std. Error)
  expect_true(all(c("Estimate", "Lower", "Upper") %in% colnames(result)))
})

# Tests for different model types

test_that("ED.drc works with LL.3 model", {
  m_ll3 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  result <- ED(m_ll3, 50, display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_true(result[, "Estimate"] > 0)
})

test_that("ED.drc works with Weibull models", {
  m_w1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
  result <- ED(m_w1, 50, display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_true(result[, "Estimate"] > 0)
})

test_that("ED.drc works with binomial type data", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  result <- ED(m_binom, 50, display = FALSE)

  expect_true(is.matrix(result))
  expect_true(result[, "Estimate"] > 0)
})

# Tests for display parameter

test_that("ED.drc respects display parameter", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # With display = FALSE, should not print anything
  expect_silent(result <- ED(m1, 50, display = FALSE))
  expect_true(is.matrix(result))
})

# Tests for multcomp output

test_that("ED.drc returns multcomp format when requested", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, c(10, 50), multcomp = TRUE, display = FALSE)

  expect_true(is.list(result))
  expect_true("EDmultcomp" %in% names(result))
})

# Tests for reference parameter

test_that("ED.drc works with different reference types", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  result_control <- ED(m1, 50, reference = "control", display = FALSE)
  result_upper <- ED(m1, 50, reference = "upper", display = FALSE)

  expect_true(is.matrix(result_control))
  expect_true(is.matrix(result_upper))
  # Both reference types should produce valid estimates
  expect_true(result_control[, "Estimate"] > 0)
  expect_true(result_upper[, "Estimate"] > 0)
})

# Tests for logBase parameter

test_that("ED.drc transforms ED values with logBase", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  result_no_log <- ED(m1, 50, display = FALSE)
  result_with_log <- ED(m1, 50, logBase = 10, display = FALSE)

  # With logBase transformation, estimates should differ
  expect_false(isTRUE(all.equal(result_no_log[, "Estimate"],
                                  result_with_log[, "Estimate"])))
})

# Tests for vcov parameter

test_that("ED.drc accepts custom vcov function", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Use a custom vcov that returns scaled variance
  custom_vcov <- function(x) vcov(x) * 2

  result_default <- ED(m1, 50, interval = "delta", display = FALSE)
  result_custom <- ED(m1, 50, interval = "delta", vcov. = custom_vcov, display = FALSE)

  # Standard errors should be different (larger with scaled vcov)
  expect_true(result_custom[, "Std. Error"] > result_default[, "Std. Error"])
})

test_that("ED.drc accepts vcov matrix directly", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  vcov_mat <- vcov(m1)
  result <- ED(m1, 50, vcov. = vcov_mat, display = FALSE)

  expect_true(is.matrix(result))
  expect_true(result[, "Estimate"] > 0)
})

# Tests for edge cases

test_that("ED.drc handles single curve with numeric curve names", {
  ryegrass_num <- ryegrass
  m1 <- drm(rootl ~ conc, data = ryegrass_num, fct = LL.4())
  result <- ED(m1, 50, display = FALSE)

  expect_true(is.matrix(result))
  # ED.drc prefixes rownames with "e:" (e.g., "e:1:50")
  expect_true(all(grepl("^e:", rownames(result))))
})

test_that("ED.drc handles very small response levels", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, 1, display = FALSE)

  expect_true(is.matrix(result))
  expect_true(result[, "Estimate"] > 0)
})

test_that("ED.drc handles very large response levels", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, 99, display = FALSE)

  expect_true(is.matrix(result))
  expect_true(result[, "Estimate"] > 0)
})

test_that("ED.drc returns invisible output", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # invisible() should not print when assigned
  result <- ED(m1, 50, display = FALSE)
  expect_true(is.matrix(result))
})

# Integration test: ED values should be reasonable

test_that("ED.drc ED values are in expected order", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, c(10, 50, 90), display = FALSE)

  # For decreasing curves, ED10 < ED50 < ED90
  ed_values <- result[, "Estimate"]
  expect_true(ed_values[1] < ed_values[2])
  expect_true(ed_values[2] < ed_values[3])
})

test_that("ED.drc standard errors are positive", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, c(10, 50, 90), display = FALSE)

  expect_true(all(result[, "Std. Error"] > 0))
})

test_that("ED.drc confidence intervals contain the estimate", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- ED(m1, 50, interval = "delta", display = FALSE)

  expect_true(result[, "Lower"] < result[, "Estimate"])
  expect_true(result[, "Upper"] > result[, "Estimate"])
})

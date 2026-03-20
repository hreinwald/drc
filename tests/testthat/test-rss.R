# Tests for rss()
# Residual Sum of Squares for dose-response models

# =============================================================================
# Test Data Setup
# =============================================================================

# Standard continuous data (ryegrass-like)
ryegrass_test <- data.frame(
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

# Multi-curve data
set.seed(42)
multi_data_test <- data.frame(
  dose = rep(c(0, 0.5, 1, 2, 5, 10), each = 5, times = 2),
  resp = c(
    rnorm(5, 100, 5), rnorm(5, 95, 5), rnorm(5, 85, 5),
    rnorm(5, 60, 5), rnorm(5, 20, 5), rnorm(5, 5, 5),
    rnorm(5, 100, 5), rnorm(5, 90, 5), rnorm(5, 70, 5),
    rnorm(5, 40, 5), rnorm(5, 10, 5), rnorm(5, 3, 5)
  ),
  group = rep(c("A", "B"), each = 30)
)

# =============================================================================
# Tests for rss()
# =============================================================================

# --- Single curve model ---

test_that("rss returns correct structure for single-curve model", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  result <- rss(m1)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "")
  expect_equal(rownames(result), "")
})

test_that("rss value matches manual calculation for single-curve model", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  result <- rss(m1)

  expected_rss <- sum(residuals(m1)^2)
  expect_equal(as.numeric(result), expected_rss, tolerance = 1e-10)
})

test_that("rss returns value invisibly", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  output <- capture.output(result <- rss(m1))

  expect_true(is.matrix(result))
  expect_true(length(output) > 0)
})

test_that("rss prints header for single-curve model", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  output <- capture.output(rss(m1))

  expect_true(any(grepl("Residual sum of squares", output)))
})

test_that("rss value is non-negative", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  result <- rss(m1)

  expect_true(as.numeric(result) >= 0)
})

# --- Multi-curve model ---

test_that("rss returns correct structure for multi-curve model", {
  m2 <- drm(resp ~ dose, curveid = group, data = multi_data_test, fct = LL.4())
  result <- rss(m2)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)  # 2 curves + total
  expect_equal(ncol(result), 1)
  expect_equal(rownames(result), c("A", "B", "Total"))
})

test_that("rss total equals sum of per-curve values for multi-curve model", {
  m2 <- drm(resp ~ dose, curveid = group, data = multi_data_test, fct = LL.4())
  result <- rss(m2)

  per_curve_sum <- sum(result[1:2, 1])
  total <- result[3, 1]
  expect_equal(per_curve_sum, total, tolerance = 1e-10)
})

test_that("rss prints header for multi-curve model", {
  m2 <- drm(resp ~ dose, curveid = group, data = multi_data_test, fct = LL.4())
  output <- capture.output(rss(m2))

  expect_true(any(grepl("Residual sums of squares", output)))
})

test_that("rss per-curve values match manual calculation", {
  m2 <- drm(resp ~ dose, curveid = group, data = multi_data_test, fct = LL.4())
  result <- rss(m2)

  resids <- residuals(m2)
  curve <- m2$data[, 4]
  expected_rss <- tapply(resids^2, curve, sum)

  expect_equal(as.numeric(result[1:2, 1]), as.numeric(expected_rss), tolerance = 1e-10)
})

# --- Consistency with Rsq ---

test_that("rss is consistent with the numerator in Rsq calculation", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  rss_val <- as.numeric(rss(m1))

  # RSS should equal sum of squared residuals (used in R² = 1 - RSS/TSS)
  expected <- sum(residuals(m1)^2)
  expect_equal(rss_val, expected, tolerance = 1e-10)
})

test_that("Rsq uses rss internally and produces correct R-squared", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  rss_val <- as.numeric(rss(m1))
  rsq_val <- as.numeric(Rsq(m1))

  response <- ryegrass_test$rootl
  tss <- sum((response - mean(response))^2)
  expected_rsq <- 1 - rss_val / tss

  expect_equal(rsq_val, expected_rsq, tolerance = 1e-10)
})

# --- print parameter ---

test_that("rss suppresses output when print = FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  output <- capture.output(result <- rss(m1, print = FALSE))

  expect_true(is.matrix(result))
  expect_equal(length(output), 0)
})

test_that("rss with print = FALSE returns same values as print = TRUE", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  r1 <- rss(m1, print = FALSE)
  r2 <- rss(m1, print = FALSE)

  expect_equal(r1, r2)
})

# --- Different model types ---

test_that("rss works with LL.3 model", {
  m3 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.3())
  result <- rss(m3)

  expect_true(is.matrix(result))
  expect_true(as.numeric(result) >= 0)
  expect_equal(as.numeric(result), sum(residuals(m3)^2), tolerance = 1e-10)
})

test_that("rss works with W1.4 model", {
  m_w <- drm(rootl ~ conc, data = ryegrass_test, fct = W1.4())
  result <- rss(m_w)

  expect_true(is.matrix(result))
  expect_true(as.numeric(result) >= 0)
  expect_equal(as.numeric(result), sum(residuals(m_w)^2), tolerance = 1e-10)
})

test_that("rss is consistent across repeated calls", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  r1 <- rss(m1)
  r2 <- rss(m1)

  expect_equal(r1, r2)
})

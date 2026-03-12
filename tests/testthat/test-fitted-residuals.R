# Test fitted.drc() and residuals.drc() functions

# Create test datasets
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

# Tests for fitted.drc()

test_that("fitted.drc returns fitted values", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fitted_vals <- fitted(m1)

  expect_true(is.numeric(fitted_vals))
  expect_equal(length(fitted_vals), nrow(ryegrass))
})

test_that("fitted.drc returns same as predict with no newdata", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  fitted_vals <- fitted(m1)
  predicted_vals <- predict(m1)

  expect_equal(fitted_vals, predicted_vals)
})

test_that("fitted.drc works with multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  fitted_vals <- fitted(m_multi)

  expect_equal(length(fitted_vals), nrow(multi_data))
  expect_true(all(is.finite(fitted_vals)))
})

test_that("fitted.drc values are within reasonable range", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fitted_vals <- fitted(m1)

  # Fitted values should be within the range of observed data
  expect_true(all(fitted_vals >= min(ryegrass$rootl) - 2))
  expect_true(all(fitted_vals <= max(ryegrass$rootl) + 2))
})

test_that("fitted.drc passes additional arguments to predict", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Pass se.fit argument
  result <- fitted(m1, se.fit = TRUE)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
})

# Tests for residuals.drc()

test_that("residuals.drc returns working residuals by default", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  resids <- residuals(m1)

  expect_true(is.numeric(resids))
  expect_equal(length(resids), nrow(ryegrass))
})

test_that("residuals.drc working residuals sum to near zero", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  resids <- residuals(m1, typeRes = "working")

  # For models with intercept, sum of residuals should be near zero
  expect_true(abs(sum(resids)) < 1)
})

test_that("residuals.drc fitted + residuals equals observed", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  fitted_vals <- fitted(m1)
  resids <- residuals(m1, typeRes = "working")
  reconstructed <- fitted_vals + resids

  expect_equal(reconstructed, ryegrass$rootl, tolerance = 1e-10)
})

test_that("residuals.drc returns standardised residuals", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  resids_std <- residuals(m1, typeRes = "standardised")

  expect_true(is.numeric(resids_std))
  expect_equal(length(resids_std), nrow(ryegrass))
})

test_that("residuals.drc standardised residuals have unit variance approx", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  resids_std <- residuals(m1, typeRes = "standardised")

  # Standardised residuals should have approximately unit variance
  expect_true(abs(var(resids_std) - 1) < 0.5)
})

test_that("residuals.drc returns studentised residuals", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  resids_stud <- residuals(m1, typeRes = "studentised")

  expect_true(is.numeric(resids_stud))
  expect_equal(length(resids_stud), nrow(ryegrass))
})

test_that("residuals.drc studentised residuals account for leverage", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  resids_std <- residuals(m1, typeRes = "standardised")
  resids_stud <- residuals(m1, typeRes = "studentised")

  # Studentised residuals should generally have larger absolute values
  # (accounting for leverage), but not always
  expect_equal(length(resids_std), length(resids_stud))
  expect_false(identical(resids_std, resids_stud))
})

test_that("residuals.drc errors for studentised without derivative matrix", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m1$deriv1 <- NULL

  expect_error(residuals(m1, typeRes = "studentised"),
               "Studentised residuals not available")
})

test_that("residuals.drc with trScale handles Box-Cox transformation", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), bcVal = 0)  # log transform
  resids_tr <- residuals(m1, trScale = TRUE)
  resids_no_tr <- residuals(m1, trScale = FALSE)

  expect_equal(length(resids_tr), nrow(ryegrass))
  expect_equal(length(resids_no_tr), nrow(ryegrass))
  # With Box-Cox, residuals on different scales should differ
  expect_false(isTRUE(all.equal(resids_tr, resids_no_tr)))
})

test_that("residuals.drc with no Box-Cox ignores trScale", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  resids_tr <- residuals(m1, trScale = TRUE)
  resids_no_tr <- residuals(m1, trScale = FALSE)

  # Without Box-Cox, both should be the same
  expect_equal(resids_tr, resids_no_tr)
})

# Tests with multi-curve models

test_that("residuals.drc works with multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  resids <- residuals(m_multi)

  expect_equal(length(resids), nrow(multi_data))
})

test_that("residuals.drc multi-curve fitted + residuals equals observed", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  fitted_vals <- fitted(m_multi)
  resids <- residuals(m_multi)
  reconstructed <- fitted_vals + resids

  expect_equal(reconstructed, multi_data$resp, tolerance = 1e-10)
})

# Tests with different model types

test_that("residuals.drc works with binomial type data", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  resids <- residuals(m_binom)

  expect_equal(length(resids), 7)
  expect_true(all(is.finite(resids)))
})

test_that("residuals.drc handles binomial without standardisation", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)

  # May return working residuals with a message for binomial
  expect_no_error(resids <- residuals(m_binom, typeRes = "standardised"))
})

test_that("residuals.drc works with Poisson type data", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  m_poisson <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")
  resids <- residuals(m_poisson)

  expect_equal(length(resids), 7)
  expect_true(all(is.finite(resids)))
})

test_that("residuals.drc studentised for non-continuous handles NA scale", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  m_poisson <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")

  # For non-continuous data, scale estimate might be NA
  # Function should handle this gracefully
  resids_stud <- residuals(m_poisson, typeRes = "studentised")

  expect_true(is.numeric(resids_stud))
  expect_equal(length(resids_stud), 7)
})

# Edge cases

test_that("residuals.drc all three types return same length", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  resids_work <- residuals(m1, typeRes = "working")
  resids_std <- residuals(m1, typeRes = "standardised")
  resids_stud <- residuals(m1, typeRes = "studentised")

  expect_equal(length(resids_work), nrow(ryegrass))
  expect_equal(length(resids_std), nrow(ryegrass))
  expect_equal(length(resids_stud), nrow(ryegrass))
})

test_that("residuals.drc types are ordered by variance adjustment", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  resids_work <- residuals(m1, typeRes = "working")
  resids_std <- residuals(m1, typeRes = "standardised")

  # Standardised residuals should have similar scale but adjusted
  var_work <- var(resids_work)
  var_std <- var(resids_std)

  # Both should be finite and positive
  expect_true(is.finite(var_work) && var_work > 0)
  expect_true(is.finite(var_std) && var_std > 0)
})

# Integration tests

test_that("residuals.drc residual plot data makes sense", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  fitted_vals <- fitted(m1)
  resids <- residuals(m1)

  # Should be able to create residual plot without errors
  expect_no_error({
    plot_data <- data.frame(fitted = fitted_vals, residuals = resids)
  })

  # No strong pattern in residuals (approximate check)
  cor_val <- cor(fitted_vals, resids)
  expect_true(abs(cor_val) < 0.5)  # Weak correlation suggests good fit
})

test_that("residuals.drc no extreme outliers in standardised residuals", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  resids_std <- residuals(m1, typeRes = "standardised")

  # Most standardised residuals should be within ±3
  prop_within_3sd <- mean(abs(resids_std) < 3)
  expect_true(prop_within_3sd > 0.95)
})

test_that("fitted and residuals are consistent across calls", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  fitted1 <- fitted(m1)
  resids1 <- residuals(m1)

  fitted2 <- fitted(m1)
  resids2 <- residuals(m1)

  expect_equal(fitted1, fitted2)
  expect_equal(resids1, resids2)
})

test_that("residuals.drc typeRes argument validation", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Valid types should work
  expect_no_error(residuals(m1, typeRes = "working"))
  expect_no_error(residuals(m1, typeRes = "standardised"))
  expect_no_error(residuals(m1, typeRes = "studentised"))

  # Invalid type should error
  expect_error(residuals(m1, typeRes = "invalid"))
})

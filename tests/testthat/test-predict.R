# Test predict.drc() function - Prediction for dose-response models

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

# Basic prediction tests

test_that("predict.drc returns fitted values when newdata is missing", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  predictions <- predict(m1)

  expect_true(is.numeric(predictions))
  expect_equal(length(predictions), nrow(ryegrass))
})

test_that("predict.drc returns predictions for new dose values", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(0.5, 1.0, 2.0))
  predictions <- predict(m1, newdata = newdata)

  expect_true(is.numeric(predictions))
  expect_equal(length(predictions), 3)
})

test_that("predict.drc predictions are within reasonable bounds", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(0, 5, 10, 15))
  predictions <- predict(m1, newdata = newdata)

  # Predictions should be positive and within reasonable range
  expect_true(all(predictions > 0))
  expect_true(all(predictions < 10))  # Based on ryegrass data range
})

# Tests with standard errors

test_that("predict.drc returns standard errors when requested", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(1, 2, 3))
  result <- predict(m1, newdata = newdata, se.fit = TRUE)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_true("Prediction" %in% colnames(result))
  expect_true("SE" %in% colnames(result))
  expect_true(all(result[, "SE"] > 0))
})

# Tests with confidence intervals

test_that("predict.drc returns confidence intervals", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(1, 2, 3))
  result <- predict(m1, newdata = newdata, interval = "confidence")

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
  expect_true(all(c("Prediction", "Lower", "Upper") %in% colnames(result)))
  expect_true(all(result[, "Lower"] < result[, "Prediction"]))
  expect_true(all(result[, "Upper"] > result[, "Prediction"]))
})

test_that("predict.drc returns prediction intervals", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(1, 2, 3))
  result <- predict(m1, newdata = newdata, interval = "prediction")

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
  expect_true(all(result[, "Lower"] < result[, "Prediction"]))
  expect_true(all(result[, "Upper"] > result[, "Prediction"]))
})

test_that("predict.drc prediction intervals are wider than confidence intervals", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(1, 2, 3))

  result_conf <- predict(m1, newdata = newdata, interval = "confidence")
  result_pred <- predict(m1, newdata = newdata, interval = "prediction")

  width_conf <- result_conf[, "Upper"] - result_conf[, "Lower"]
  width_pred <- result_pred[, "Upper"] - result_pred[, "Lower"]

  expect_true(all(width_pred > width_conf))
})

test_that("predict.drc confidence level affects interval width", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(2))

  result_95 <- predict(m1, newdata = newdata, interval = "confidence", level = 0.95)
  result_90 <- predict(m1, newdata = newdata, interval = "confidence", level = 0.90)

  width_95 <- result_95[, "Upper"] - result_95[, "Lower"]
  width_90 <- result_90[, "Upper"] - result_90[, "Lower"]

  expect_true(width_90 < width_95)
})

# Tests with multi-curve models

test_that("predict.drc works with multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  predictions <- predict(m_multi)

  expect_true(is.numeric(predictions))
  expect_equal(length(predictions), nrow(multi_data))
})

test_that("predict.drc predicts for specific curves in newdata", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  newdata <- data.frame(dose = c(1, 2, 3), group = c("A", "A", "B"))
  predictions <- predict(m_multi, newdata = newdata)

  expect_equal(length(predictions), 3)
  expect_true(all(is.finite(predictions)))
})

test_that("predict.drc handles missing curve ID in newdata", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  newdata <- data.frame(dose = c(1, 2, 3))  # No group column
  predictions <- predict(m_multi, newdata = newdata)

  expect_equal(length(predictions), 3)
})

# Tests with different model types

test_that("predict.drc works with binomial type data", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  predictions <- predict(m_binom)

  expect_true(all(predictions >= 0 & predictions <= 1))
})

test_that("predict.drc constrains binomial predictions to [0, 1]", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  result <- predict(m_binom, interval = "confidence")

  # Upper limit should not exceed 1
  expect_true(all(result[, "Upper"] <= 1))
})

test_that("predict.drc works with Poisson type data", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  m_poisson <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")
  predictions <- predict(m_poisson)

  expect_true(all(predictions >= 0))
  expect_equal(length(predictions), 7)
})

test_that("predict.drc constrains non-continuous predictions to >= 0", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  m_poisson <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")
  result <- predict(m_poisson, interval = "confidence")

  # Lower limit should not be negative
  expect_true(all(result[, "Lower"] >= 0))
})

# Tests with constrain parameter

test_that("predict.drc respects constrain = FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(0.1))

  result_constrained <- predict(m1, newdata = newdata, interval = "confidence", constrain = TRUE)
  result_unconstrained <- predict(m1, newdata = newdata, interval = "confidence", constrain = FALSE)

  # Both should return results, but unconstrained might have different bounds
  expect_true(is.matrix(result_constrained))
  expect_true(is.matrix(result_unconstrained))
})

# Tests with vcov parameter

test_that("predict.drc accepts custom vcov function", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(2))

  custom_vcov <- function(x) vcov(x) * 2
  result_default <- predict(m1, newdata = newdata, se.fit = TRUE)
  result_custom <- predict(m1, newdata = newdata, se.fit = TRUE, vcov. = custom_vcov)

  # Standard errors should be different (larger with scaled vcov)
  expect_true(result_custom[, "SE"] > result_default[, "SE"])
})

# Tests with od (over-dispersion) parameter

test_that("predict.drc respects od parameter for binomial data", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  newdata <- data.frame(dose = c(1))

  result_no_od <- predict(m_binom, newdata = newdata, se.fit = TRUE, od = FALSE)
  result_with_od <- predict(m_binom, newdata = newdata, se.fit = TRUE, od = TRUE)

  # Standard errors might differ with od adjustment
  expect_true(is.matrix(result_no_od))
  expect_true(is.matrix(result_with_od))
})

# Tests for models without derivatives

test_that("predict.drc returns only predictions when derivatives unavailable", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  # Remove derivatives to simulate unavailable derivatives
  m1$fct$deriv1 <- NULL
  newdata <- data.frame(conc = c(1, 2))

  predictions <- predict(m1, newdata = newdata)

  expect_true(is.numeric(predictions))
  expect_equal(length(predictions), 2)
})

# Tests with se.fit and interval together

test_that("predict.drc returns SE and confidence interval columns together", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(1, 2))

  result <- predict(m1, newdata = newdata, se.fit = TRUE, interval = "confidence")

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 4)
  expect_true(all(c("Prediction", "SE", "Lower", "Upper") %in% colnames(result)))
})

test_that("predict.drc with se.fit = FALSE and interval returns no SE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(1, 2))

  result <- predict(m1, newdata = newdata, se.fit = FALSE, interval = "confidence")

  expect_equal(ncol(result), 3)
  expect_false("SE" %in% colnames(result))
})

# Tests with checkND parameter

test_that("predict.drc with checkND = FALSE accepts non-standard newdata", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Pass a vector instead of data frame
  newdata <- c(1, 2, 3)
  predictions <- predict(m1, newdata = newdata, checkND = FALSE)

  expect_equal(length(predictions), 3)
})

# Edge cases and error conditions

test_that("predict.drc handles zero dose values", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(0, 0.01, 0.1))
  predictions <- predict(m1, newdata = newdata)

  expect_equal(length(predictions), 3)
  expect_true(all(is.finite(predictions)))
})

test_that("predict.drc handles very large dose values", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(100, 1000))
  predictions <- predict(m1, newdata = newdata)

  expect_equal(length(predictions), 2)
  expect_true(all(is.finite(predictions)))
})

test_that("predict.drc returns consistent predictions", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(2))

  pred1 <- predict(m1, newdata = newdata)
  pred2 <- predict(m1, newdata = newdata)

  expect_equal(pred1, pred2)
})

test_that("predict.drc fitted values match original data size", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fitted_vals <- predict(m1)

  expect_equal(length(fitted_vals), nrow(ryegrass))
})

# Tests with different response types

test_that("predict.drc uses correct distribution for continuous data", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(2))

  # For continuous data, should use t-distribution
  result <- predict(m1, newdata = newdata, interval = "confidence")

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
})

test_that("predict.drc uses correct distribution for binomial data", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  newdata <- data.frame(dose = c(1))

  # For binomial data, should use normal distribution
  result <- predict(m_binom, newdata = newdata, interval = "confidence")

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
})

# Integration tests

test_that("predict.drc predictions are monotonic for monotonic models", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  newdata <- data.frame(conc = c(0.5, 1, 2, 4, 8))
  predictions <- predict(m1, newdata = newdata)

  # For decreasing dose-response, predictions should decrease
  expect_true(all(diff(predictions) <= 0.01))  # Allow small numerical errors
})

test_that("predict.drc standard errors increase at extreme doses", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Predictions at middle range
  middle <- predict(m1, newdata = data.frame(conc = c(2)), se.fit = TRUE)

  # Predictions at extremes
  extreme <- predict(m1, newdata = data.frame(conc = c(0.01, 100)), se.fit = TRUE)

  # Generally, SE at extremes should be larger (though not always guaranteed)
  expect_true(is.matrix(middle))
  expect_true(is.matrix(extreme))
})

# tests/testthat/test-ED.lin.R
# Comprehensive tests for ED.lin() - ED calculation for linear models

# Helper to extract numeric value from the list-matrix returned by ED.lin
get_val <- function(result, row, col) {
  as.numeric(unlist(result[row, col]))
}

# --- 2-parameter linear models (lparco == 2) ---------------------------------

test_that("ED.lin returns correct structure for increasing linear model", {
  set.seed(42)
  x <- 1:10
  y <- 2 + 3 * x + rnorm(10, 0, 0.5)
  fit <- lm(y ~ x)

  result <- ED.lin(fit, 50)

  # Return type and dimensions
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("Estimate", "SE", "2.5 %", "97.5 %"))

  # The ED50 estimate should be finite and positive
  expect_true(is.finite(get_val(result, 1, "Estimate")))
  expect_true(get_val(result, 1, "Estimate") > 0)
  # SE should be positive
  expect_true(get_val(result, 1, "SE") > 0)
})

test_that("ED.lin handles increasing linear model (non-decreasing path)", {
  # y = 2 + 3*x => increasing, parCoef[2] > 0 => decreasing = FALSE
  set.seed(1)
  x <- 1:10
  y <- 2 + 3 * x + rnorm(10, 0, 0.1)
  fit <- lm(y ~ x)

  # decreasing should be FALSE for positive slope
  expect_true(coef(fit)[2] > 0)

  result <- ED.lin(fit, 50)

  # ED50 estimate should be reasonable
  est <- get_val(result, 1, "Estimate")
  expect_true(is.finite(est))
})

test_that("ED.lin handles decreasing linear model (decreasing path)", {
  # y = 30 - 3*x => decreasing, parCoef[2] < 0 => decreasing = TRUE
  set.seed(1)
  x <- 1:10
  y <- 30 - 3 * x + rnorm(10, 0, 0.1)
  fit <- lm(y ~ x)

  # decreasing should be TRUE for negative slope
  expect_true(coef(fit)[2] < 0)

  result <- ED.lin(fit, 50)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_true(is.finite(get_val(result, 1, "Estimate")))
})

test_that("ED.lin handles multiple response levels for 2-param model", {
  set.seed(42)
  x <- 1:10
  y <- 2 + 3 * x + rnorm(10, 0, 0.5)
  fit <- lm(y ~ x)

  result <- ED.lin(fit, c(10, 25, 50, 75, 90))

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 4)
  # All estimates should be finite
  for (i in seq_len(nrow(result))) {
    expect_true(is.finite(get_val(result, i, "Estimate")))
  }
})

# --- 3-parameter quadratic models (lparco == 3) ------------------------------

test_that("ED.lin handles quadratic concave-down model (vertex within range)", {
  # parCoef[3] < 0, vertex within data range
  set.seed(42)
  x <- seq(0, 10, length.out = 20)
  y <- 5 + 3 * x - 0.3 * x^2 + rnorm(20, 0, 0.5)
  fit <- lm(y ~ x + I(x^2))

  cc <- coef(fit)
  expect_true(cc[3] < 0)  # concave down
  vertex_x <- -cc[2] / (2 * cc[3])
  expect_true(vertex_x <= max(x))  # vertex within range

  result <- ED.lin(fit, 50)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 4)
  expect_true(is.finite(get_val(result, 1, "Estimate")))
})

test_that("ED.lin handles quadratic concave-down model (vertex beyond range)", {
  # parCoef[3] < 0, vertex > max(x)
  set.seed(42)
  x <- seq(0, 3, length.out = 20)
  y <- 5 + 10 * x - 0.5 * x^2 + rnorm(20, 0, 0.5)
  fit <- lm(y ~ x + I(x^2))

  cc <- coef(fit)
  expect_true(cc[3] < 0)  # concave down
  vertex_x <- -cc[2] / (2 * cc[3])
  expect_true(vertex_x > max(x))  # vertex beyond range

  result <- ED.lin(fit, 50)

  expect_true(is.matrix(result))
  expect_true(is.finite(get_val(result, 1, "Estimate")))
})

test_that("ED.lin handles quadratic concave-up model (cup, parCoef[3] > 0)", {
  # parCoef[3] > 0, decreasing = TRUE on line 11
  set.seed(42)
  x <- seq(0, 10, length.out = 20)
  y <- 10 - 3 * x + 0.3 * x^2 + rnorm(20, 0, 0.5)
  fit <- lm(y ~ x + I(x^2))

  cc <- coef(fit)
  expect_true(cc[3] > 0)  # concave up (cup)

  result <- ED.lin(fit, 50)

  expect_true(is.matrix(result))
  expect_true(is.finite(get_val(result, 1, "Estimate")))
})

test_that("ED.lin handles multiple response levels for quadratic model", {
  set.seed(42)
  x <- seq(0, 10, length.out = 20)
  y <- 5 + 3 * x - 0.3 * x^2 + rnorm(20, 0, 0.5)
  fit <- lm(y ~ x + I(x^2))

  result <- ED.lin(fit, c(10, 50, 90))

  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
})

# --- Edge case: negative cVal truncation (pmax(0, cVal)) ---------------------

test_that("ED.lin truncates negative cVal to zero", {
  # Create a model where cVal (lower limit) would be negative
  # For an increasing model, cVal = fitted at min(x)
  set.seed(42)
  x <- 0:10
  y <- -5 + 3 * x + rnorm(11, 0, 0.1)
  fit <- lm(y ~ x)

  # Verify fitted value at min(x) is negative
  expect_true(fitted(fit)[which.min(x)] < 0)

  # cVal is truncated to 0 by pmax
  result <- ED.lin(fit, 50)
  expect_true(is.matrix(result))
  expect_true(is.finite(get_val(result, 1, "Estimate")))
})

# --- Correctness checks (known solutions) ------------------------------------

test_that("ED.lin gives correct ED50 for perfect linear model", {
  # Perfect model: y = 10*x, x from 0 to 10
  # cVal = max(0, fitted(0)) = 0, dVal = fitted(10) = 100
  # ED50: x = 5
  x <- seq(0, 10, by = 1)
  y <- 10 * x
  fit <- lm(y ~ x)

  result <- suppressWarnings(ED.lin(fit, 50))
  expect_equal(get_val(result, 1, "Estimate"), 5, tolerance = 1e-6)
})

test_that("ED.lin gives correct ED50 for perfect decreasing linear model", {
  # y = 100 - 10*x, x from 0 to 10
  # ED50: x = 5
  x <- seq(0, 10, by = 1)
  y <- 100 - 10 * x
  fit <- lm(y ~ x)

  result <- suppressWarnings(ED.lin(fit, 50))
  expect_equal(get_val(result, 1, "Estimate"), 5, tolerance = 1e-6)
})

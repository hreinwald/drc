# Tests for yieldLoss() and genFixedFct() from R/yieldLoss.R

# --- Setup: Fit a Michaelis-Menten model using the methionine dataset ---
# This mirrors the example from the function's documentation
met.mm.m1 <- drm(gain ~ dose, product, data = methionine, fct = MM.3(),
                 pmodels = list(~1, ~factor(product), ~factor(product)))

# =============================================================================
# Tests for yieldLoss()
# =============================================================================

# --- Happy Path: interval = "none" (default) ---
test_that("yieldLoss returns correct structure with interval='none'", {
  result <- yieldLoss(met.mm.m1, display = FALSE)

  expect_type(result, "list")
  expect_named(result, c("A", "I"))

  # A and I should be matrices

  expect_true(is.matrix(result$A))
  expect_true(is.matrix(result$I))

  # With interval="none", should have 2 columns: Estimate and Std. Error
  expect_equal(ncol(result$A), 2)
  expect_equal(ncol(result$I), 2)
  expect_equal(colnames(result$A), c("Estimate", "Std. Error"))
  expect_equal(colnames(result$I), c("Estimate", "Std. Error"))

  # Row names should match curve names from the model
  expect_equal(rownames(result$A), colnames(met.mm.m1$parmMat))
  expect_equal(rownames(result$I), colnames(met.mm.m1$parmMat))

  # Estimates and standard errors should be finite numeric values
  expect_true(all(is.finite(result$A)))
  expect_true(all(is.finite(result$I)))

  # Standard errors should be positive
  expect_true(all(result$A[, "Std. Error"] > 0))
  expect_true(all(result$I[, "Std. Error"] > 0))
})

# --- Happy Path: interval = "as" with continuous model ---
test_that("yieldLoss returns confidence intervals with interval='as'", {
  result <- yieldLoss(met.mm.m1, interval = "as", display = FALSE)

  expect_type(result, "list")
  expect_named(result, c("A", "I"))

  # With interval="as", should have 4 columns
  expect_equal(ncol(result$A), 4)
  expect_equal(ncol(result$I), 4)
  expect_equal(colnames(result$A), c("Estimate", "Std. Error", "Lower", "Upper"))
  expect_equal(colnames(result$I), c("Estimate", "Std. Error", "Lower", "Upper"))

  # Lower should be less than Estimate, Upper should be greater
  expect_true(all(result$A[, "Lower"] < result$A[, "Estimate"]))
  expect_true(all(result$A[, "Upper"] > result$A[, "Estimate"]))
  expect_true(all(result$I[, "Lower"] < result$I[, "Estimate"]))
  expect_true(all(result$I[, "Upper"] > result$I[, "Estimate"]))
})

# --- Happy Path: different confidence levels ---
test_that("yieldLoss respects the level parameter", {
  result_95 <- yieldLoss(met.mm.m1, interval = "as", level = 0.95, display = FALSE)
  result_99 <- yieldLoss(met.mm.m1, interval = "as", level = 0.99, display = FALSE)

  # 99% CIs should be wider than 95% CIs
  width_95_A <- result_95$A[, "Upper"] - result_95$A[, "Lower"]
  width_99_A <- result_99$A[, "Upper"] - result_99$A[, "Lower"]
  expect_true(all(width_99_A > width_95_A))
})

# --- Display behavior ---
test_that("yieldLoss prints output when display=TRUE", {
  output <- capture.output(result <- yieldLoss(met.mm.m1, display = TRUE))
  expect_true(length(output) > 0)
  expect_true(any(grepl("Estimated A parameters", output)))
  expect_true(any(grepl("Estimated I parameters", output)))
})

test_that("yieldLoss suppresses output when display=FALSE", {
  output <- capture.output(result <- yieldLoss(met.mm.m1, display = FALSE))
  expect_equal(length(output), 0)
})

# --- Display with interval = "as" ---
test_that("yieldLoss displays correctly with interval='as'", {
  output <- capture.output(result <- yieldLoss(met.mm.m1, interval = "as", display = TRUE))
  expect_true(length(output) > 0)
  expect_true(any(grepl("Estimated A parameters", output)))
})

# --- Error handling: invalid interval argument ---
test_that("yieldLoss errors on invalid interval argument", {
  expect_error(yieldLoss(met.mm.m1, interval = "invalid"))
})

# --- Non-continuous model type path (covers qnorm branch in ciFct) ---
test_that("yieldLoss uses qnorm for non-continuous model types with interval='as'", {
  # Create a modified copy of the model with non-continuous type
  mock_model <- met.mm.m1
  mock_model$type <- "binomial"

  result <- yieldLoss(mock_model, interval = "as", display = FALSE)

  expect_type(result, "list")
  expect_named(result, c("A", "I"))
  expect_equal(ncol(result$A), 4)
  expect_equal(ncol(result$I), 4)

  # Verify the CIs differ from the continuous case (qt vs qnorm)
  result_cont <- yieldLoss(met.mm.m1, interval = "as", display = FALSE)
  # CIs should differ because qt and qnorm give different quantiles
  expect_false(identical(result$A[, "Lower"], result_cont$A[, "Lower"]))
})

# =============================================================================
# Tests for genFixedFct() (internal helper)
# =============================================================================

test_that("genFixedFct works with allComp=TRUE", {
  # fixed has some NAs (free) and some fixed values
  # MM.3 has fixed = c(-1, NA, NA, NA, 1) for b, c, d, e, f
  fixed <- c(-1, NA, NA, NA, 1)
  fct <- drc:::genFixedFct(fixed)

  # Provide values for the 3 free parameters (c, d, e)
  result <- fct(c(10, 20, 30), allComp = TRUE)
  expect_equal(result, c(-1, 10, 20, 30, 1))
})

test_that("genFixedFct works with allComp=FALSE", {
  fixed <- c(-1, NA, NA, NA, 1)
  fct <- drc:::genFixedFct(fixed)

  # With allComp=FALSE, it should filter the parm vector to free positions
  parm <- c(100, 200, 300, 400, 500)
  result <- fct(parm, allComp = FALSE)
  # notFixed positions are 2, 3, 4 (the NA positions)
  expect_equal(result, c(200, 300, 400))
})

test_that("genFixedFct handles all-fixed parameters", {
  fixed <- c(1, 2, 3)
  fct <- drc:::genFixedFct(fixed)

  # allComp=TRUE: all positions are fixed, parm should be ignored
  result <- fct(numeric(0), allComp = TRUE)
  expect_equal(result, c(1, 2, 3))

  # allComp=FALSE: no free parameters, return empty
  result2 <- fct(c(10, 20, 30), allComp = FALSE)
  expect_length(result2, 0)
})

test_that("genFixedFct handles all-free parameters", {
  fixed <- c(NA, NA, NA)
  fct <- drc:::genFixedFct(fixed)

  result <- fct(c(10, 20, 30), allComp = TRUE)
  expect_equal(result, c(10, 20, 30))

  result2 <- fct(c(10, 20, 30), allComp = FALSE)
  expect_equal(result2, c(10, 20, 30))
})

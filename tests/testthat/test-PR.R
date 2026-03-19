# tests/testthat/test-PR.R
# Comprehensive tests for PR() function (R/pr.R)

# --- Setup: Fit models used across multiple tests ---

# Single-curve model using ryegrass dataset
ryegrass_model <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

# Multi-curve model using S.alba dataset (two herbicide curves)
salba_model <- drm(DryMatter ~ Dose, curveid = Herbicide,
                   data = S.alba, fct = LL.4())

# -------------------------------------------------------------------
# Section 1: Correctness tests ("Happy Path") — single-curve models
# -------------------------------------------------------------------

test_that("PR returns named numeric vector for single-curve model", {
  result <- drc:::PR(ryegrass_model, c(1, 5, 10))

  # Should be a numeric vector, not a matrix

  expect_false(is.matrix(result))
  expect_true(is.numeric(result))
  expect_length(result, 3)

  # Names should correspond to xVec values

  expect_equal(names(result), c("1", "5", "10"))
})

test_that("PR predictions match predict.drc for single-curve model", {
  xvals <- c(0.5, 2, 8)
  pr_result <- drc:::PR(ryegrass_model, xvals)
  pred_result <- predict(ryegrass_model, data.frame(xvals))

  expect_equal(unname(pr_result), as.numeric(pred_result))
})

test_that("PR handles single dose value for single-curve model", {
  result <- drc:::PR(ryegrass_model, 5)

  expect_false(is.matrix(result))
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_equal(names(result), "5")
})

# -------------------------------------------------------------------
# Section 2: Single-curve model with se.fit = TRUE (matrix return)
# Covers lines 34-37: the is.matrix(retMat) == TRUE branch
# -------------------------------------------------------------------

test_that("PR returns matrix with rownames when se.fit = TRUE (single-curve)", {
  result <- drc:::PR(ryegrass_model, c(2, 10), se.fit = TRUE)

  # predict.drc with se.fit=TRUE returns a matrix
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_true(ncol(result) >= 2)  # At least Prediction + SE

  # Rownames should be the dose values as strings
  expect_equal(rownames(result), c("2", "10"))
})

# -------------------------------------------------------------------
# Section 3: Correctness tests — multi-curve models
# Covers lines 29-32: the lenCI > 1 branch
# -------------------------------------------------------------------

test_that("PR returns matrix with rownames for multi-curve model", {
  xvals <- c(1, 10)
  result <- drc:::PR(salba_model, xvals)

  # Multi-curve always returns a matrix (se.fit = TRUE is hardcoded)
  expect_true(is.matrix(result))

  # Curve IDs are stored as numeric codes internally
  curve_ids <- as.character(unique(salba_model$data[, 3]))
  expected_nrow <- length(xvals) * length(curve_ids)
  expect_equal(nrow(result), expected_nrow)

  # Rownames should be curveId:dose format
  expected_names <- paste(rep(curve_ids, each = length(xvals)),
                          rep(as.character(xvals), length(curve_ids)),
                          sep = ":")
  expect_equal(rownames(result), expected_names)
})

test_that("PR handles single dose value for multi-curve model", {
  result <- drc:::PR(salba_model, 5)

  expect_true(is.matrix(result))

  curve_ids <- as.character(unique(salba_model$data[, 3]))
  expect_equal(nrow(result), length(curve_ids))
})

# -------------------------------------------------------------------
# Section 4: Edge cases
# -------------------------------------------------------------------

test_that("PR works with dose value of zero", {
  result <- drc:::PR(ryegrass_model, 0)
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_equal(names(result), "0")
})

test_that("PR works with very large dose values", {
  result <- drc:::PR(ryegrass_model, c(1e6))
  expect_true(is.numeric(result))
  expect_length(result, 1)
})

test_that("PR works with many dose values", {
  xvals <- seq(0, 30, by = 0.5)
  result <- drc:::PR(ryegrass_model, xvals)
  expect_length(result, length(xvals))
  expect_equal(names(result), as.character(xvals))
})

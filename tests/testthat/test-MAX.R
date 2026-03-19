# tests/testthat/test-MAX.R
# Comprehensive tests for drc:::MAX to achieve 100% code coverage.

# ──────────────────────────────────────────────────────────────────────
# Helper: fit a simple hormesis model for the happy-path tests
# ──────────────────────────────────────────────────────────────────────
make_crs_model <- function() {
  data(lettuce, package = "drc", envir = environment())
  drm(weight ~ conc, data = lettuce, fct = CRS.4c())
}

make_bc_model <- function() {
  data(lettuce, package = "drc", envir = environment())
  drm(weight ~ conc, data = lettuce, fct = BC.4())
}

# ====================================================================
# 1. Input validation – object class (lines 81-83)
# ====================================================================
test_that("MAX errors when object is not class 'drc'", {
  expect_error(
    drc:::MAX("not_a_model"),
    "'object' must be of class 'drc'"
  )
  expect_error(
    drc:::MAX(list(a = 1)),
    "'object' must be of class 'drc'"
  )
  expect_error(
    drc:::MAX(42),
    "'object' must be of class 'drc'"
  )
})

# ====================================================================
# 2. Input validation – no maxfct method (lines 88-94)
# ====================================================================
test_that("MAX errors when model has no 'maxfct' method", {
  # LL.4 is a log-logistic model with no hormesis => maxfct is NULL
  data(ryegrass, package = "drc", envir = environment())
  m_ll4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(
    drc:::MAX(m_ll4),
    "No 'maxfct' method available"
  )
})

# ====================================================================
# 3. Input validation – lower bound (lines 99-101)
# ====================================================================
test_that("MAX errors for invalid 'lower' argument", {
  m <- make_crs_model()
  # Non-numeric

expect_error(drc:::MAX(m, lower = "a"),
               "'lower' must be a single finite numeric value")
  # Length > 1
  expect_error(drc:::MAX(m, lower = c(1, 2)),
               "'lower' must be a single finite numeric value")
  # NA
  expect_error(drc:::MAX(m, lower = NA_real_),
               "'lower' must be a single finite numeric value")
  # Inf
  expect_error(drc:::MAX(m, lower = Inf),
               "'lower' must be a single finite numeric value")
  # -Inf
  expect_error(drc:::MAX(m, lower = -Inf),
               "'lower' must be a single finite numeric value")
  # NaN
  expect_error(drc:::MAX(m, lower = NaN),
               "'lower' must be a single finite numeric value")
})

# ====================================================================
# 4. Input validation – upper bound (lines 102-104)
# ====================================================================
test_that("MAX errors for invalid 'upper' argument", {
  m <- make_crs_model()
  # Non-numeric
  expect_error(drc:::MAX(m, upper = "b"),
               "'upper' must be a single finite numeric value")
  # Length > 1
  expect_error(drc:::MAX(m, upper = c(100, 200)),
               "'upper' must be a single finite numeric value")
  # NA
  expect_error(drc:::MAX(m, upper = NA_real_),
               "'upper' must be a single finite numeric value")
  # Inf
  expect_error(drc:::MAX(m, upper = Inf),
               "'upper' must be a single finite numeric value")
})

# ====================================================================
# 5. Input validation – lower >= upper (lines 105-109)
# ====================================================================
test_that("MAX errors when lower >= upper", {
  m <- make_crs_model()
  expect_error(drc:::MAX(m, lower = 100, upper = 100),
               "'lower'.*must be strictly less than 'upper'")
  expect_error(drc:::MAX(m, lower = 200, upper = 100),
               "'lower'.*must be strictly less than 'upper'")
})

# ====================================================================
# 6. Happy path – CRS.4c model (cedergreen class)
# ====================================================================
test_that("MAX returns correct structure for CRS.4c model", {
  m <- make_crs_model()
  result <- drc:::MAX(m)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("Dose", "Response"))
  expect_true(nrow(result) >= 1)
  # Values should be finite positive numbers for this dataset
  expect_true(all(is.finite(result)))
  expect_true(all(result[, "Dose"] > 0))
  expect_true(all(result[, "Response"] > 0))
})

# ====================================================================
# 7. Happy path – BC.4 model (braincousens class)
# ====================================================================
test_that("MAX returns correct structure for BC.4 model", {
  m <- make_bc_model()
  result <- drc:::MAX(m)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("Dose", "Response"))
  expect_true(nrow(result) >= 1)
  expect_true(all(is.finite(result)))
})

# ====================================================================
# 8. Custom search interval
# ====================================================================
test_that("MAX works with custom lower/upper bounds", {
  m <- make_crs_model()
  result <- drc:::MAX(m, lower = 1e-5, upper = 500)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_true(all(is.finite(result)))
})

# ====================================================================
# 9. pool = FALSE path (lines 119-125, alternate vcov call)
# ====================================================================
test_that("MAX works with pool = FALSE", {
  m <- make_crs_model()
  result <- drc:::MAX(m, pool = FALSE)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_true(all(is.finite(result)))
})

# ====================================================================
# 10. vcov computation failure – warning path (lines 119-125)
# ====================================================================
test_that("MAX warns when vcov computation fails", {
  m <- make_crs_model()
  # Corrupt the model's fit component so vcov.drc errors internally
  m$fit <- NULL

  expect_warning(
    result <- drc:::MAX(m),
    "Could not compute variance-covariance matrix"
  )

  # Should still return valid results despite vcov failure
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
})

# ====================================================================
# 11. Boundary warning – dose at lower bound (lines 151-158)
# ====================================================================
test_that("MAX warns when maximum dose is at the lower boundary", {
  m <- make_crs_model()
  # Use a very high lower bound so the result must be at the boundary
  # The lettuce model max dose is around 0.02, so lower = 50 forces boundary
  expect_warning(
    result <- drc:::MAX(m, lower = 50, upper = 100),
    "is at the boundary"
  )
})

# ====================================================================
# 12. Boundary warning – dose at upper bound (lines 151-158)
# ====================================================================
test_that("MAX warns when maximum dose is at the upper boundary", {
  m <- make_crs_model()
  # Use a very small upper bound so the result sits at the boundary
  expect_warning(
    result <- drc:::MAX(m, lower = 1e-8, upper = 1e-7),
    "is at the boundary"
  )
})

# ====================================================================
# 13. maxfct error → warning + NA (lines 161-167)
# ====================================================================
test_that("MAX returns NA and warns when maxfct computation fails for a curve", {
  m <- make_crs_model()

  # Replace the maxfct function with one that errors
  m$fct$maxfct <- function(parm, lower, upper, ...) {
    stop("simulated maxfct failure")
  }

  expect_warning(
    result <- drc:::MAX(m),
    "MAX computation failed for curve"
  )
  # All values should be NA
  expect_true(all(is.na(result)))
})

# ====================================================================
# 14. Curve name fallback when strParm is NULL (line 141)
# ====================================================================
test_that("MAX falls back to 'Curve_i' when strParm is NULL", {
  m <- make_crs_model()
  # Remove column names from parmMat to trigger fallback
  colnames(m$parmMat) <- NULL

  result <- drc:::MAX(m)
  expect_true(is.matrix(result))
  expect_true(all(grepl("^Curve_", rownames(result))))
})

# ====================================================================
# 15. Curve name fallback when strParm is NA (line 141)
# ====================================================================
test_that("MAX falls back to 'Curve_i' when strParm contains NA", {
  m <- make_crs_model()
  # Set column names to NA to trigger the fallback
  colnames(m$parmMat) <- NA_character_

  result <- drc:::MAX(m)
  expect_true(is.matrix(result))
  expect_true(all(grepl("^Curve_", rownames(result))))
})

# ====================================================================
# 16. Return value is invisible (line 176)
# ====================================================================
test_that("MAX returns invisibly", {
  m <- make_crs_model()
  # expect_invisible checks that the return value is invisible
  expect_invisible(drc:::MAX(m))
})

# ====================================================================
# 17. Multiple curves
# ====================================================================
test_that("MAX works with multiple curves (multi-curve data)", {
  data(lettuce, package = "drc", envir = environment())
  # Create a fake curveid column to have two "curves"
  lettuce2 <- rbind(
    transform(lettuce, curveid = "A"),
    transform(lettuce, curveid = "B", weight = weight * 1.1)
  )
  m <- drm(weight ~ conc, curveid = curveid, data = lettuce2, fct = CRS.4c())
  result <- drc:::MAX(m)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("Dose", "Response"))
})

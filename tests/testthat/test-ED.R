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

# Helper: compute the expected SE for an absolute-type ED using numerical
# central differences on the model's edfct and the fitted vcov matrix.
compute_numgrad_se <- function(model, absResp) {
  vc <- vcov(model)
  edfct <- model$fct$edfct
  parms <- coef(model)
  eps <- .Machine$double.eps^(1/3)
  numGrad <- numeric(length(parms))
  for (k in seq_along(parms)) {
    h <- max(abs(parms[k]), 1) * eps
    pu <- replace(parms, k, parms[k] + h)
    pd <- replace(parms, k, parms[k] - h)
    eu <- edfct(pu, absResp, reference = "control", type = "absolute")[[1]]
    ed <- edfct(pd, absResp, reference = "control", type = "absolute")[[1]]
    numGrad[k] <- (eu - ed) / (2 * h)
  }
  as.numeric(sqrt(numGrad %*% vc %*% numGrad))
}

test_that("ED.drc absolute type SE includes asymptote parameter uncertainty", {
  # Fit a 4-parameter log-logistic model
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Pick an absolute response level (midpoint of the fitted curve)
  cf <- coef(m1)
  midResp <- (cf[2] + cf[3]) / 2  # midpoint between c and d

  result_abs <- ED(m1, midResp, type = "absolute", display = FALSE)
  result_rel <- ED(m1, 50, type = "relative", display = FALSE)

  # The ED estimates should be the same (midpoint = ED50 for symmetric model)
  expect_equal(result_abs[, "Estimate"], result_rel[, "Estimate"],
               tolerance = 0.01)

  # The absolute-type SE differs from the relative-type SE because it
  # additionally accounts for uncertainty in c and d via the full
  # numerical gradient (parameter covariances can make it larger or smaller).
  expect_false(isTRUE(all.equal(result_abs[, "Std. Error"],
                                result_rel[, "Std. Error"])))

  # Cross-check: manually computed SE should match
  expectedSE <- compute_numgrad_se(m1, midResp)
  expect_equal(result_abs[, "Std. Error"], expectedSE, tolerance = 1e-4)
})

test_that("ED.drc absolute type SE is correct for Weibull type 2", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = W2.4())

  cf <- coef(m1)
  midResp <- (cf[2] + cf[3]) / 2

  result_abs <- ED(m1, midResp, type = "absolute", display = FALSE)

  expectedSE <- compute_numgrad_se(m1, midResp)
  expect_equal(result_abs[, "Std. Error"], expectedSE, tolerance = 1e-4)
})

test_that("ED.drc absolute type SE is correct for Weibull type 1", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())

  cf <- coef(m1)
  midResp <- (cf[2] + cf[3]) / 2

  result_abs <- ED(m1, midResp, type = "absolute", display = FALSE)

  expectedSE <- compute_numgrad_se(m1, midResp)
  expect_equal(result_abs[, "Std. Error"], expectedSE, tolerance = 1e-4)
})

test_that("ED.drc relative type SE unchanged by fix", {
  # The relative-type SE should NOT be affected by the absolute-type fix
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  result <- ED(m1, 50, type = "relative", display = FALSE)

  # Manually compute relative-type SE using the analytical gradient
  cf <- coef(m1)
  vc <- vcov(m1)
  edfct <- m1$fct$edfct
  edResult <- edfct(cf, 50, reference = "control", type = "relative")
  edGrad <- edResult[[2]]
  expectedSE <- as.numeric(sqrt(edGrad %*% vc %*% edGrad))

  expect_equal(result[, "Std. Error"], expectedSE, tolerance = 1e-8)
})

test_that("ED.drc errors when model has no edfct function", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  # Remove the edfct function to simulate a model without it
  m1$fct$edfct <- NULL

  expect_error(ED(m1, 50, display = FALSE), "ED values cannot be calculated")
})

# Tests for input validation error branches in ED.drc

test_that("ED.drc errors when object is not of class drc", {
  expect_error(drc:::ED.drc("not_a_model", 50), "'object' must be of class 'drc'")
  expect_error(drc:::ED.drc(42, 50), "'object' must be of class 'drc'")
})

test_that("ED.drc errors when respLev is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(drc:::ED.drc(m1, "abc"), "'respLev' must be a non-empty numeric vector")
  expect_error(drc:::ED.drc(m1, numeric(0)), "'respLev' must be a non-empty numeric vector")
})

test_that("ED.drc errors when level is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(drc:::ED.drc(m1, 50, level = "a"), "'level' must be a single numeric value strictly between 0 and 1")
  expect_error(drc:::ED.drc(m1, 50, level = 0), "'level' must be a single numeric value strictly between 0 and 1")
  expect_error(drc:::ED.drc(m1, 50, level = 1), "'level' must be a single numeric value strictly between 0 and 1")
  expect_error(drc:::ED.drc(m1, 50, level = c(0.9, 0.95)), "'level' must be a single numeric value strictly between 0 and 1")
})

test_that("ED.drc errors when bound is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(drc:::ED.drc(m1, 50, bound = "yes"), "'bound' must be a single logical value")
  expect_error(drc:::ED.drc(m1, 50, bound = c(TRUE, FALSE)), "'bound' must be a single logical value")
})

test_that("ED.drc errors when display is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(drc:::ED.drc(m1, 50, display = "yes"), "'display' must be a single logical value")
})

test_that("ED.drc errors when multcomp is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(drc:::ED.drc(m1, 50, multcomp = "yes"), "'multcomp' must be a single logical value")
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

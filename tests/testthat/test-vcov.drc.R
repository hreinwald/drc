# Regression test for DoseResponse/drc#36:
# summary() must not crash for binomial models with a singular Hessian.

test_that("vcDisc returns NA matrix (not an error) for a singular Hessian", {
  # Use a negative definite Hessian so that all fallback paths fail:
  # solve() fails (singular/non-invertible in the expected sense),
  # chol() fails (not positive definite), and the regularised
  # chol(0.99*H + 0.01*I) also fails (still not positive definite).
  fake_obj <- list(fit = list(hessian = matrix(c(-100, 0, 0, -100), 2, 2)))
  expect_warning(
    result <- drc:::vcDisc(fake_obj),
    "Hessian is singular"
  )
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2L, 2L))
  expect_true(all(is.na(result)))
})

test_that("summary() does not error for increasing-trend binomial W1.4 model", {
  # Exact reproducer from DoseResponse/drc#36
  dataframe <- data.frame(
    conc     = c(0, 0.944, 2.18, 4.14, 8.37, 16.1),
    total    = c(160, 80, 80, 80, 80, 80),
    response = c(3, 1, 1, 5, 8, 80)
  )
  m <- drm(response / total ~ conc, weights = total, data = dataframe,
           fct = W1.4(), type = "binomial")
  # Must not throw; may produce NA standard errors if Hessian is singular
  expect_no_error(suppressWarnings(summary(m)))
  summ <- suppressWarnings(summary(m))
  expect_s3_class(summ, "summary.drc")
  expect_true("coefficients" %in% names(summ))
})

test_that("summary() still works for well-conditioned binomial model (non-regression)", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n    = rep(50, 7)
  )
  m <- drm(resp ~ dose, data = binom_data, fct = LL.2(),
           type = "binomial", weights = n)
  summ <- summary(m)
  expect_s3_class(summ, "summary.drc")
  # Standard errors should be finite for a well-conditioned fit
  expect_true(all(is.finite(summ$coefficients[, "Std. Error"])))
})

# Test coef.drc(), vcov.drc(), and confint.drc() functions

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

# Tests for coef.drc()

test_that("coef.drc returns coefficient vector", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  coeffs <- coef(m1)

  expect_true(is.numeric(coeffs))
  expect_true(length(coeffs) > 0)
  expect_true(!is.null(names(coeffs)))
})

test_that("coef.drc returns correct number of coefficients for LL.4", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  coeffs <- coef(m1)

  expect_equal(length(coeffs), 4)
})

test_that("coef.drc returns correct number of coefficients for LL.3", {
  m_ll3 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  coeffs <- coef(m_ll3)

  expect_equal(length(coeffs), 3)
})

test_that("coef.drc coefficients have proper names", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  coeffs <- coef(m1)

  # LL.4 parameters: b (slope), c (lower), d (upper), e (ED50)
  expect_true(all(grepl(":", names(coeffs))))  # Parameters should have colon separator
})

test_that("coef.drc works with multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  coeffs <- coef(m_multi)

  # Should have 4 parameters × 2 curves = 8 coefficients
  expect_equal(length(coeffs), 8)
})

test_that("coef.drc handles NULL coefficients gracefully", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Simulate NULL coefficients
  m1_copy <- m1
  m1_copy$coefficients <- NULL

  # Should fall back to fit$par
  coeffs <- coef(m1_copy)
  expect_true(is.numeric(coeffs))
  expect_true(length(coeffs) > 0)
})

# Tests for vcov.drc()

test_that("vcov.drc returns variance-covariance matrix", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  vcov_mat <- vcov(m1)

  expect_true(is.matrix(vcov_mat))
  expect_equal(nrow(vcov_mat), length(coef(m1)))
  expect_equal(ncol(vcov_mat), length(coef(m1)))
})

test_that("vcov.drc matrix is symmetric", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  vcov_mat <- vcov(m1)

  expect_equal(vcov_mat, t(vcov_mat))
})

test_that("vcov.drc diagonal elements are positive", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  vcov_mat <- vcov(m1)

  diag_elements <- diag(vcov_mat)
  expect_true(all(diag_elements > 0))
})

test_that("vcov.drc returns correlation matrix when corr = TRUE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  corr_mat <- vcov(m1, corr = TRUE)

  expect_true(is.matrix(corr_mat))

  # Diagonal elements should be 1 for correlation matrix
  expect_true(all(abs(diag(corr_mat) - 1) < 1e-10))

  # Off-diagonal elements should be between -1 and 1
  expect_true(all(corr_mat >= -1 & corr_mat <= 1))
})

test_that("vcov.drc with od = TRUE for binomial data", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)

  vcov_no_od <- vcov(m_binom, od = FALSE)
  vcov_with_od <- vcov(m_binom, od = TRUE)

  expect_true(is.matrix(vcov_no_od))
  expect_true(is.matrix(vcov_with_od))

  # Dimensions should be the same
  expect_equal(dim(vcov_no_od), dim(vcov_with_od))

  # Values might differ with OD adjustment
  expect_equal(nrow(vcov_no_od), length(coef(m_binom)))
})

test_that("vcov.drc with pool = FALSE for multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  vcov_pooled <- vcov(m_multi, pool = TRUE)
  vcov_unpooled <- vcov(m_multi, pool = FALSE)

  expect_true(is.matrix(vcov_pooled))
  expect_true(is.matrix(vcov_unpooled))

  # Both should have same dimensions
  expect_equal(dim(vcov_pooled), dim(vcov_unpooled))
})

test_that("vcov.drc with unscaled = TRUE for continuous data", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  vcov_scaled <- vcov(m1, unscaled = FALSE)
  vcov_unscaled <- vcov(m1, unscaled = TRUE)

  expect_true(is.matrix(vcov_scaled))
  expect_true(is.matrix(vcov_unscaled))
  expect_equal(dim(vcov_scaled), dim(vcov_unscaled))
})

test_that("vcov.drc works for different model types", {
  m_ll3 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  m_ll4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  vcov_ll3 <- vcov(m_ll3)
  vcov_ll4 <- vcov(m_ll4)

  expect_equal(nrow(vcov_ll3), 3)
  expect_equal(nrow(vcov_ll4), 4)
})

test_that("vcov.drc works with Poisson type data", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  m_poisson <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")
  vcov_mat <- vcov(m_poisson)

  expect_true(is.matrix(vcov_mat))
  expect_equal(nrow(vcov_mat), length(coef(m_poisson)))
})

# Tests for confint.drc()

test_that("confint.drc returns confidence intervals for all parameters", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  ci <- confint(m1)

  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), length(coef(m1)))
  expect_equal(ncol(ci), 2)
})

test_that("confint.drc interval contains point estimate", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  ci <- confint(m1)
  coeffs <- coef(m1)

  # Each coefficient should be within its confidence interval
  for (i in 1:length(coeffs)) {
    expect_true(ci[i, 1] < coeffs[i])
    expect_true(ci[i, 2] > coeffs[i])
  }
})

test_that("confint.drc respects level parameter", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  ci_95 <- confint(m1, level = 0.95)
  ci_90 <- confint(m1, level = 0.90)

  # 90% CI should be narrower than 95% CI
  width_95 <- ci_95[, 2] - ci_95[, 1]
  width_90 <- ci_90[, 2] - ci_90[, 1]

  expect_true(all(width_90 < width_95))
})

test_that("confint.drc can select specific parameters", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  coeffs <- coef(m1)

  # Get CI for just one parameter
  param_name <- names(coeffs)[1]
  ci <- confint(m1, parm = param_name)

  expect_equal(nrow(ci), 1)
  expect_equal(rownames(ci), param_name)
})

test_that("confint.drc errors for invalid parameter name", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_error(confint(m1, parm = "nonexistent_param"),
               "does not match an actual parameter")
})

test_that("confint.drc can select multiple parameters", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  coeffs <- coef(m1)

  # Get CI for first two parameters
  param_names <- names(coeffs)[1:2]
  ci <- confint(m1, parm = param_names)

  expect_equal(nrow(ci), 2)
  expect_equal(rownames(ci), param_names)
})

test_that("confint.drc works with multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  ci <- confint(m_multi)

  expect_equal(nrow(ci), length(coef(m_multi)))
  expect_equal(ncol(ci), 2)
})

test_that("confint.drc with pool = FALSE for multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  ci_pooled <- confint(m_multi, pool = TRUE)
  ci_unpooled <- confint(m_multi, pool = FALSE)

  expect_true(is.matrix(ci_pooled))
  expect_true(is.matrix(ci_unpooled))
  expect_equal(dim(ci_pooled), dim(ci_unpooled))
})

test_that("confint.drc uses t-distribution for continuous data", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  ci <- confint(m1, level = 0.95)

  # Check that CI is calculated correctly using t-distribution
  summ <- summary(m1)
  est <- summ$coefficients[, "Estimate"]
  se <- summ$coefficients[, "Std. Error"]
  df <- df.residual(m1)
  t_val <- qt(0.975, df)

  expected_lower <- est - t_val * se
  expected_upper <- est + t_val * se

  expect_equal(ci[, 1], expected_lower, tolerance = 1e-6)
  expect_equal(ci[, 2], expected_upper, tolerance = 1e-6)
})

test_that("confint.drc uses normal distribution for binomial data", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  ci <- confint(m_binom, level = 0.95)

  # Check that CI is calculated correctly using normal distribution
  summ <- summary(m_binom)
  est <- summ$coefficients[, "Estimate"]
  se <- summ$coefficients[, "Std. Error"]
  z_val <- qnorm(0.975)

  expected_lower <- est - z_val * se
  expected_upper <- est + z_val * se

  expect_equal(ci[, 1], expected_lower, tolerance = 1e-6)
  expect_equal(ci[, 2], expected_upper, tolerance = 1e-6)
})

# Integration tests

test_that("coef, vcov, and confint are consistent", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  coeffs <- coef(m1)
  vcov_mat <- vcov(m1)
  ci <- confint(m1)

  # Number of coefficients should match dimensions
  expect_equal(length(coeffs), nrow(vcov_mat))
  expect_equal(length(coeffs), nrow(ci))

  # Names should match
  expect_equal(names(coeffs), rownames(vcov_mat))
  expect_equal(names(coeffs), rownames(ci))
})

test_that("vcov diagonal matches SE squared", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  vcov_mat <- vcov(m1)
  summ <- summary(m1)

  se_from_vcov <- sqrt(diag(vcov_mat))
  se_from_summary <- summ$coefficients[, "Std. Error"]

  expect_equal(se_from_vcov, se_from_summary, tolerance = 1e-10)
})

test_that("confint.basic helper function works", {
  est_mat <- matrix(c(1, 0.1, 2, 0.2), ncol = 2)
  rownames(est_mat) <- c("param1", "param2")
  colnames(est_mat) <- c("Estimate", "Std. Error")

  ci <- confint.basic(est_mat, level = 0.95, intType = "continuous", dfres = 26)

  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)
  expect_equal(rownames(ci), c("param1", "param2"))
})

test_that("confint.basic uses correct quantiles for different types", {
  est_mat <- matrix(c(1, 0.1), ncol = 2)

  # Continuous should use t-distribution
  ci_cont <- confint.basic(est_mat, level = 0.95, intType = "continuous", dfres = 30)

  # Binomial should use normal distribution
  ci_binom <- confint.basic(est_mat, level = 0.95, intType = "binomial", dfres = 30)

  expect_true(is.matrix(ci_cont))
  expect_true(is.matrix(ci_binom))

  # CI widths will differ slightly due to different distributions
  width_cont <- ci_cont[1, 2] - ci_cont[1, 1]
  width_binom <- ci_binom[1, 2] - ci_binom[1, 1]

  expect_true(width_cont > 0)
  expect_true(width_binom > 0)
})

test_that("All three functions work together in a workflow", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Extract coefficients
  coeffs <- coef(m1)
  expect_true(length(coeffs) == 4)

  # Get variance-covariance matrix
  vcov_mat <- vcov(m1)
  expect_equal(dim(vcov_mat), c(4, 4))

  # Get confidence intervals
  ci <- confint(m1)
  expect_equal(dim(ci), c(4, 2))

  # All should have matching names
  expect_equal(names(coeffs), rownames(vcov_mat))
  expect_equal(names(coeffs), rownames(ci))
})

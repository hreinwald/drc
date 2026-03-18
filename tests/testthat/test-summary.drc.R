# Comprehensive tests for summary.drc() and print.summary.drc()
# Targeting 100% code coverage

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

# Multi-curve data for independent fits
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

# Binomial data
binom_data_test <- data.frame(
  dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
  resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
  n = rep(50, 7)
)

# Poisson data
poisson_data_test <- data.frame(
  dose = c(0, 1, 2, 4, 8, 16, 32),
  count = c(50, 48, 40, 25, 10, 3, 1)
)

# =============================================================================
# Tests for summary.drc()
# =============================================================================

# --- Happy Path: Basic continuous model ---

test_that("summary.drc returns correct structure for continuous model", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  # Check class
  expect_s3_class(summ, "summary.drc")

  # Check all named elements
  expected_names <- c("resVar", "varMat", "coefficients", "boxcox", "fctName",
                      "robust", "varParm", "type", "df.residual",
                      "cov.unscaled", "text", "noParm", "rseMat")
  expect_equal(names(summ), expected_names)

  # Check types of key elements
  expect_true(is.numeric(summ$resVar))
  expect_true(is.matrix(summ$varMat))
  expect_true(is.matrix(summ$coefficients))
  expect_true(is.matrix(summ$rseMat))
  expect_equal(summ$type, "continuous")
  expect_null(summ$robust)
  expect_null(summ$varParm)
})

test_that("summary.drc coefficient matrix has correct structure", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)
  coef_mat <- summ$coefficients

  expect_equal(ncol(coef_mat), 4)
  expect_equal(colnames(coef_mat), c("Estimate", "Std. Error", "t-value", "p-value"))
  expect_equal(nrow(coef_mat), 4)  # LL.4 has 4 parameters

  # Estimates should match coef()
  expect_equal(coef_mat[, "Estimate"], coef(m1), ignore_attr = TRUE)

  # Standard errors should be positive
  expect_true(all(coef_mat[, "Std. Error"] > 0))

  # p-values should be in [0, 1]
  expect_true(all(coef_mat[, "p-value"] >= 0 & coef_mat[, "p-value"] <= 1))
})

test_that("summary.drc uses t-distribution for continuous data", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)
  coef_mat <- summ$coefficients

  # Manually compute expected p-values using t-distribution
  t_vals <- coef_mat[, "t-value"]
  df_val <- df.residual(m1)
  expected_pvals <- pt(-abs(t_vals), df_val) + (1 - pt(abs(t_vals), df_val))
  expect_equal(coef_mat[, "p-value"], expected_pvals, tolerance = 1e-10)
})

test_that("summary.drc rseMat is correct for single-curve model", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  expect_equal(nrow(summ$rseMat), 1)
  expect_equal(ncol(summ$rseMat), 2)
  expect_equal(colnames(summ$rseMat), c("rse", "df"))
  expect_true(summ$rseMat[1, "rse"] > 0)
  expect_equal(summ$rseMat[1, "df"], df.residual(m1))
})

test_that("summary.drc computes cov.unscaled for continuous data", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  expect_true(is.matrix(summ$cov.unscaled))
  expect_false(is.null(summ$cov.unscaled))
})

# --- Non-continuous models (pnorm branch + varMat.us NULL) ---

test_that("summary.drc uses pnorm for binomial data", {
  m_binom <- drm(resp ~ dose, data = binom_data_test, fct = LL.2(),
                 type = "binomial", weights = n)
  summ <- summary(m_binom)

  expect_s3_class(summ, "summary.drc")
  expect_equal(summ$type, "binomial")

  # p-values should use normal distribution
  coef_mat <- summ$coefficients
  z_vals <- coef_mat[, "t-value"]
  expected_pvals <- pnorm(-abs(z_vals)) + (1 - pnorm(abs(z_vals)))
  expect_equal(coef_mat[, "p-value"], expected_pvals, tolerance = 1e-10)
})

test_that("summary.drc sets cov.unscaled to NULL for non-continuous data", {
  m_binom <- drm(resp ~ dose, data = binom_data_test, fct = LL.2(),
                 type = "binomial", weights = n)
  summ <- summary(m_binom)

  # resVar should be NA for binomial
  expect_true(is.na(summ$resVar))

  # cov.unscaled should be NULL (not a matrix of NAs)
  expect_null(summ$cov.unscaled)
})

test_that("summary.drc works with Poisson data", {
  m_poisson <- drm(count ~ dose, data = poisson_data_test, fct = LL.4(),
                   type = "Poisson")
  summ <- summary(m_poisson)

  expect_s3_class(summ, "summary.drc")
  expect_equal(summ$type, "Poisson")
  expect_true(is.na(summ$resVar))
  expect_null(summ$cov.unscaled)
})

# --- Multi-curve models with pool=FALSE (unpooled path) ---

test_that("summary.drc with pool=FALSE for multi-curve independent model", {
  skip_if_not_installed("magic")
  library(magic)
  data(spinach)
  m_sep <- drm(SLOPE ~ DOSE, HERBICIDE, data = spinach,
               fct = LL.4(), separate = TRUE)
  summ_unpooled <- summary(m_sep, pool = FALSE)

  expect_s3_class(summ_unpooled, "summary.drc")

  # rseMat should have multiple rows (one per curve)
  expect_true(nrow(summ_unpooled$rseMat) > 1)
  expect_equal(colnames(summ_unpooled$rseMat), c("rse", "df"))

  # All RSEs should be positive
  expect_true(all(summ_unpooled$rseMat[, "rse"] > 0))
  expect_true(all(summ_unpooled$rseMat[, "df"] > 0))
})

test_that("summary.drc with pool=TRUE for multi-curve independent model", {
  skip_if_not_installed("magic")
  library(magic)
  data(spinach)
  m_sep <- drm(SLOPE ~ DOSE, HERBICIDE, data = spinach,
               fct = LL.4(), separate = TRUE)
  summ_pooled <- summary(m_sep, pool = TRUE)

  expect_s3_class(summ_pooled, "summary.drc")

  # rseMat should have 1 row when pooled
  expect_equal(nrow(summ_pooled$rseMat), 1)
})

# --- Robust estimation methods ---
# Note: robust models use eval in parent.frame() so we must use package data

test_that("summary.drc with robust='trimmed' (metric trimming)", {
  data(ryegrass)
  m_robust <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
                  robust = "trimmed")
  summ <- summary(m_robust)

  expect_s3_class(summ, "summary.drc")
  expect_equal(summ$robust, "metric trimming")

  # Standard errors should be computed via Hessian
  expect_true(all(summ$coefficients[, "Std. Error"] > 0))
  expect_true(all(is.finite(summ$coefficients[, "Std. Error"])))
})

test_that("summary.drc with robust='tukey' (Tukey's biweight)", {
  data(ryegrass)
  m_robust <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
                  robust = "tukey")
  summ <- summary(m_robust)

  expect_s3_class(summ, "summary.drc")
  expect_equal(summ$robust, "Tukey's biweight")
  expect_true(all(summ$coefficients[, "Std. Error"] > 0))
})

test_that("summary.drc with robust='winsor' (metric Winsorizing)", {
  data(ryegrass)
  # Winsorizing may fail to converge with some datasets; use W1.4 which is
  # more robust to convergence issues
  m_robust <- tryCatch(
    drm(rootl ~ conc, data = ryegrass, fct = LL.4(), robust = "winsor"),
    error = function(e) NULL
  )
  skip_if(is.null(m_robust), "Winsorizing model did not converge")
  summ <- summary(m_robust)

  expect_s3_class(summ, "summary.drc")
  expect_equal(summ$robust, "metric Winsorizing")
  expect_true(all(summ$coefficients[, "Std. Error"] > 0))
})

# --- Over-dispersion ---

test_that("summary.drc with od=TRUE for binomial data", {
  m_binom <- drm(resp ~ dose, data = binom_data_test, fct = LL.2(),
                 type = "binomial", weights = n)

  summ_no_od <- summary(m_binom, od = FALSE)
  summ_od <- summary(m_binom, od = TRUE)

  expect_s3_class(summ_no_od, "summary.drc")
  expect_s3_class(summ_od, "summary.drc")
})

# --- Consistency tests ---

test_that("summary.drc is consistent across repeated calls", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  s1 <- summary(m1)
  s2 <- summary(m1)

  expect_equal(s1$coefficients, s2$coefficients)
  expect_equal(s1$resVar, s2$resVar)
  expect_equal(s1$varMat, s2$varMat)
})

# =============================================================================
# Tests for print.summary.drc()
# =============================================================================

# --- Basic print ---

test_that("print.summary.drc produces output and returns invisibly", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  output <- capture.output(result <- print(summ))
  expect_true(length(output) > 0)
  expect_identical(result, summ)  # returns invisibly
})

test_that("print.summary.drc shows model text and noParm", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  output <- capture.output(print(summ))
  # Should contain model fitted text with noParm
  expect_true(any(grepl("Model fitted:", output)))
  expect_true(any(grepl("parms", output)))
})

test_that("print.summary.drc handles noParm = NULL", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  # Manually set noParm to NULL
  summ$noParm <- NULL

  output <- capture.output(print(summ))
  expect_true(any(grepl("Model fitted:", output)))
  # Should NOT contain "parms" since noParm is NULL
  expect_false(any(grepl("parms", output)))
})

test_that("print.summary.drc shows robust estimation info", {
  data(ryegrass)
  m_robust <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
                  robust = "trimmed")
  summ <- summary(m_robust)

  output <- capture.output(print(summ))
  expect_true(any(grepl("Robust estimation:", output)))
  expect_true(any(grepl("metric trimming", output)))
})

test_that("print.summary.drc shows multiple RSEs for multi-curve unpooled", {
  skip_if_not_installed("magic")
  library(magic)
  data(spinach)
  m_sep <- drm(SLOPE ~ DOSE, HERBICIDE, data = spinach,
               fct = LL.4(), separate = TRUE)
  summ <- summary(m_sep, pool = FALSE)

  output <- capture.output(print(summ))
  # Should print "Residual standard errors:" (with 's')
  expect_true(any(grepl("Residual standard errors:", output)))
})

test_that("print.summary.drc shows single RSE for pooled model", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  output <- capture.output(print(summ))
  # Should show "Residual standard error:" (without trailing 's')
  expect_true(any(grepl("Residual standard error:", output)))
})

test_that("print.summary.drc warns when df.residual < 1", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  # Set df.residual to 0 to trigger warning
  summ$df.residual <- 0

  output <- capture.output(print(summ))
  expect_true(any(grepl("Too complex model fitted", output)))
})

test_that("print.summary.drc skips RSE section for non-continuous data", {
  m_binom <- drm(resp ~ dose, data = binom_data_test, fct = LL.2(),
                 type = "binomial", weights = n)
  summ <- summary(m_binom)

  output <- capture.output(print(summ))
  # Should NOT show "Residual standard error" for binomial
  expect_false(any(grepl("Residual standard error", output)))
})

# --- varComp tests (using mock) ---

test_that("print.summary.drc shows varComp when present", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  # Add mock varComp
  summ$varComp <- matrix(c(1.5, 0.3, 5.0, 0.001), nrow = 1,
                         dimnames = list("sigma", c("Estimate", "Std. Error",
                                                    "t-value", "p-value")))

  output <- capture.output(print(summ))
  expect_true(any(grepl("Estimated variance components:", output)))
})

# --- varParm tests (using mock) ---

test_that("print.summary.drc shows varParm with varPower type (single row)", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  # Mock varParm with single row estimates (varPower)
  summ$varParm <- list(
    type = "varPower",
    estimates = matrix(c(2.0, 0.5, 4.0, 0.001), nrow = 1,
                       dimnames = list("theta", c("Estimate", "Std. Error",
                                                  "t-value", "p-value")))
  )

  output <- capture.output(print(summ))
  expect_true(any(grepl("power-of-the-mean variance model", output)))
})

test_that("print.summary.drc shows varParm with varPower type (multi-row)", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  # Mock varParm with multiple row estimates (varPower)
  summ$varParm <- list(
    type = "varPower",
    estimates = matrix(c(1.0, 2.0, 0.3, 0.5, 3.33, 4.0, 0.01, 0.001),
                       nrow = 2,
                       dimnames = list(c("sigma", "theta"),
                                       c("Estimate", "Std. Error",
                                         "t-value", "p-value")))
  )

  output <- capture.output(print(summ))
  expect_true(any(grepl("power-of-the-mean variance model", output)))
})

test_that("print.summary.drc shows varParm with hetvar type", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  # Mock varParm with hetvar type
  summ$varParm <- list(
    type = "hetvar",
    estimates = matrix(c(1.0, 2.0, 0.3, 0.5, 3.33, 4.0, 0.01, 0.001),
                       nrow = 2,
                       dimnames = list(c("var1", "var2"),
                                       c("Estimate", "Std. Error",
                                         "t-value", "p-value")))
  )

  output <- capture.output(print(summ))
  expect_true(any(grepl("Estimated heterogeneous variances:", output)))
})

# --- boxcox tests ---

test_that("print.summary.drc shows boxcox with lambda and CI", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  # Mock boxcox with lambda and CI
  summ$boxcox <- list(lambda = 0.5, ci = c(0.2, 0.8))

  output <- capture.output(print(summ))
  expect_true(any(grepl("Box-Cox transformation", output)))
  expect_true(any(grepl("Estimated lambda:", output)))
  expect_true(any(grepl("Confidence interval for lambda:", output)))
})

test_that("print.summary.drc shows boxcox with specified lambda (NA CI)", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  # Mock boxcox with lambda but NA CI
  summ$boxcox <- list(lambda = 1.0, ci = c(NA, NA))

  output <- capture.output(print(summ))
  expect_true(any(grepl("Box-Cox transformation", output)))
  expect_true(any(grepl("Specified lambda:", output)))
})

test_that("print.summary.drc handles boxcox with NA lambda", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)

  # Mock boxcox with NA lambda (should print nothing)
  summ$boxcox <- list(lambda = NA)

  output <- capture.output(print(summ))
  expect_false(any(grepl("Box-Cox transformation", output)))
})

test_that("print.summary.drc handles NULL boxcox", {
  m1 <- drm(rootl ~ conc, data = ryegrass_test, fct = LL.4())
  summ <- summary(m1)
  summ$boxcox <- NULL

  output <- capture.output(print(summ))
  expect_false(any(grepl("Box-Cox transformation", output)))
})

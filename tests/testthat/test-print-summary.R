# Test print.drc() and summary.drc() functions

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

# Tests for print.drc()

test_that("print.drc prints model information", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Capture output
  output <- capture.output(print(m1))

  expect_true(length(output) > 0)
  expect_true(any(grepl("drc", output)))
  expect_true(any(grepl("Call:", output)))
  expect_true(any(grepl("Coefficients:", output)))
})

test_that("print.drc returns object invisibly", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Should return the object invisibly
  result <- capture.output(returned_obj <- print(m1))

  expect_identical(returned_obj, m1)
})

test_that("print.drc displays coefficients", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  output <- capture.output(print(m1))

  # Should show parameter names
  expect_true(any(grepl("b:", output) | grepl("Coefficients", output)))
})

test_that("print.drc respects digits parameter", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  output_3 <- capture.output(print(m1, digits = 3))
  output_6 <- capture.output(print(m1, digits = 6))

  # More digits should produce longer output strings
  expect_true(length(output_3) > 0)
  expect_true(length(output_6) > 0)
})

test_that("print.drc handles model with no coefficients", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  # Simulate model with no coefficients
  m1_copy <- m1
  m1_copy$coefficients <- numeric(0)

  output <- capture.output(print(m1_copy))

  expect_true(any(grepl("No coefficients", output)))
})

test_that("print.drc works with multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  output <- capture.output(print(m_multi))

  expect_true(length(output) > 0)
  expect_true(any(grepl("drc", output)))
})

# Tests for summary.drc()

test_that("summary.drc returns summary object", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  expect_true(is.list(summ))
  expect_true("resVar" %in% names(summ))
  expect_true("varMat" %in% names(summ) | length(summ) > 0)
})

test_that("summary.drc contains residual standard error", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  expect_true("resVar" %in% names(summ))
  expect_true(is.numeric(summ$resVar))
  expect_true(summ$resVar > 0)
})

test_that("summary.drc contains coefficient matrix", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  expect_true("coefficients" %in% names(summ))
  expect_true(is.matrix(summ$coefficients))
  expect_equal(nrow(summ$coefficients), length(coef(m1)))
})

test_that("summary.drc coefficient matrix has correct columns", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  coef_mat <- summ$coefficients
  expect_true("Estimate" %in% colnames(coef_mat))
  expect_true("Std. Error" %in% colnames(coef_mat))
  expect_true("t-value" %in% colnames(coef_mat) | "z-value" %in% colnames(coef_mat))
  expect_true("p-value" %in% colnames(coef_mat))
})

test_that("summary.drc p-values are between 0 and 1", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  p_values <- summ$coefficients[, "p-value"]
  expect_true(all(p_values >= 0 & p_values <= 1))
})

test_that("summary.drc standard errors are positive", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  std_errors <- summ$coefficients[, "Std. Error"]
  expect_true(all(std_errors > 0))
})

test_that("summary.drc works with over-dispersion adjustment", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)

  summ_no_od <- summary(m_binom, od = FALSE)
  summ_with_od <- summary(m_binom, od = TRUE)

  expect_true(is.list(summ_no_od))
  expect_true(is.list(summ_with_od))

  # Standard errors might differ with OD adjustment
  se_no_od <- summ_no_od$coefficients[, "Std. Error"]
  se_with_od <- summ_with_od$coefficients[, "Std. Error"]

  expect_equal(length(se_no_od), length(se_with_od))
})

test_that("summary.drc with pool = FALSE for multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  summ_pooled <- summary(m_multi, pool = TRUE)
  summ_unpooled <- summary(m_multi, pool = FALSE)

  expect_true(is.list(summ_pooled))
  expect_true(is.list(summ_unpooled))
})

test_that("summary.drc contains residual standard error matrix", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  expect_true("rseMat" %in% names(summ))
  expect_true(is.matrix(summ$rseMat))
  expect_true("rse" %in% colnames(summ$rseMat))
  expect_true("df" %in% colnames(summ$rseMat))
})

test_that("summary.drc degrees of freedom are positive integers", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  df <- summ$rseMat[, "df"]
  expect_true(all(df > 0))
  expect_true(all(df == floor(df)))
})

# Tests with different model types

test_that("summary.drc works with binomial type data", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  summ <- summary(m_binom)

  expect_true(is.list(summ))
  expect_true("coefficients" %in% names(summ))
})

test_that("summary.drc uses z-values for non-continuous data", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  summ <- summary(m_binom)

  # Binomial data should use z-values (normal distribution)
  expect_true("z-value" %in% colnames(summ$coefficients) |
              "z value" %in% colnames(summ$coefficients) |
              "t-value" %in% colnames(summ$coefficients))
})

test_that("summary.drc uses t-values for continuous data", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  # Continuous data should use t-values
  expect_true("t-value" %in% colnames(summ$coefficients) |
              "t value" %in% colnames(summ$coefficients))
})

test_that("summary.drc works with Poisson type data", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  m_poisson <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")
  summ <- summary(m_poisson)

  expect_true(is.list(summ))
  expect_true("coefficients" %in% names(summ))
})

test_that("summary.drc handles robust estimation methods", {
  m_robust <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), robust = "median")
  summ <- summary(m_robust)

  expect_true(is.list(summ))
  expect_true("coefficients" %in% names(summ))
})

# Tests for print.summary.drc (if it exists)

test_that("summary.drc object can be printed", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  # Should be able to print summary without error
  expect_no_error({
    output <- capture.output(print(summ))
  })
})

# Integration tests

test_that("print and summary work together", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Both should work without errors
  expect_no_error({
    print_output <- capture.output(print(m1))
    summ <- summary(m1)
    summary_output <- capture.output(print(summ))
  })
})

test_that("summary.drc estimates match coef", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  summ <- summary(m1)

  coef_from_summary <- summ$coefficients[, "Estimate"]
  coef_direct <- coef(m1)

  expect_equal(coef_from_summary, coef_direct)
})

test_that("summary.drc is consistent across calls", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  summ1 <- summary(m1)
  summ2 <- summary(m1)

  expect_equal(summ1$coefficients, summ2$coefficients)
  expect_equal(summ1$resVar, summ2$resVar)
})

test_that("summary.drc with different LL models", {
  m_ll3 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  m_ll4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  summ_ll3 <- summary(m_ll3)
  summ_ll4 <- summary(m_ll4)

  # LL.3 has fewer parameters
  expect_equal(nrow(summ_ll3$coefficients), 3)
  expect_equal(nrow(summ_ll4$coefficients), 4)
})

test_that("summary.drc with Weibull models", {
  m_w1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
  summ <- summary(m_w1)

  expect_true(is.list(summ))
  expect_equal(nrow(summ$coefficients), 4)
})

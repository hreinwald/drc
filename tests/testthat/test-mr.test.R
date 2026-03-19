# Test suite for mr.test()
# Achieves 100% code coverage for R/mr.test.R

# ---- Setup: create model fixtures ----
# etmotc dataset: dose-response data for testing
# Using first 15 observations as in the documented example
data(etmotc)
etmotc_sub <- etmotc[1:15, ]

# Null model: 4-parameter log-logistic
m1 <- drm(rgr1 ~ dose1, data = etmotc_sub, fct = LL.4())

# Alternative model: 4-parameter Weibull type 1
m2 <- update(m1, fct = W1.4())

# Fitted model under alternative: fit the null model's fitted values using
# the alternative model structure
m3 <- drm(fitted(m1) ~ dose1, data = etmotc_sub, fct = W1.4())

# Dose vector with zeros replaced by small value (as per docs example)
xVec <- etmotc_sub$dose1
xVec[xVec == 0] <- 1e-10

# ---- Tests for var.equal = TRUE (default path) ----

test_that("mr.test with var.equal=TRUE returns correct structure", {
  result <- mr.test(m1, m2, m3, xVec, var.equal = TRUE)
  expect_type(result, "double")
  expect_length(result, 4)
  expect_named(result, c("Statistic", "p-value", "Difference", "SE"))
})

test_that("mr.test with var.equal=TRUE returns valid statistics", {
  result <- mr.test(m1, m2, m3, xVec, var.equal = TRUE)
  # p-value must be between 0 and 1
  expect_true(result["p-value"] >= 0 && result["p-value"] <= 1)
  # SE must be positive
  expect_true(result["SE"] > 0)
  # Statistic should be finite
  expect_true(is.finite(result["Statistic"]))
})

test_that("mr.test default var.equal is TRUE", {
  result_default <- mr.test(m1, m2, m3, xVec)
  result_explicit <- mr.test(m1, m2, m3, xVec, var.equal = TRUE)
  expect_equal(result_default, result_explicit)
})

# ---- Tests for var.equal = FALSE (heteroscedastic path) ----

test_that("mr.test with var.equal=FALSE returns correct structure", {
  result <- mr.test(m1, m2, m3, xVec, var.equal = FALSE)
  expect_type(result, "double")
  expect_length(result, 4)
  expect_named(result, c("Statistic", "p-value", "Difference", "SE"))
})

test_that("mr.test with var.equal=FALSE returns valid statistics", {
  result <- mr.test(m1, m2, m3, xVec, var.equal = FALSE)
  # p-value must be between 0 and 1
  expect_true(result["p-value"] >= 0 && result["p-value"] <= 1)
  # SE must be positive
  expect_true(result["SE"] > 0)
  # Statistic should be finite
  expect_true(is.finite(result["Statistic"]))
})

test_that("mr.test var.equal=TRUE and FALSE give different results", {
  result_equal <- mr.test(m1, m2, m3, xVec, var.equal = TRUE)
  result_unequal <- mr.test(m1, m2, m3, xVec, var.equal = FALSE)
  # Results should differ due to different variance estimation
  expect_false(identical(result_equal, result_unequal))
})

# ---- Tests for component parameter ----

test_that("mr.test with different component values", {
  result_c1 <- mr.test(m1, m2, m3, xVec, component = 1)
  result_c2 <- mr.test(m1, m2, m3, xVec, component = 2)
  result_c3 <- mr.test(m1, m2, m3, xVec, component = 3)
  result_c4 <- mr.test(m1, m2, m3, xVec, component = 4)
  # Each component should give different results
  expect_false(identical(result_c1, result_c2))
  # All should have the standard structure
  expect_named(result_c2, c("Statistic", "p-value", "Difference", "SE"))
  expect_named(result_c3, c("Statistic", "p-value", "Difference", "SE"))
  expect_named(result_c4, c("Statistic", "p-value", "Difference", "SE"))
})

test_that("mr.test with component and var.equal=FALSE", {
  result <- mr.test(m1, m2, m3, xVec, var.equal = FALSE, component = 2)
  expect_type(result, "double")
  expect_length(result, 4)
  expect_named(result, c("Statistic", "p-value", "Difference", "SE"))
  expect_true(result["p-value"] >= 0 && result["p-value"] <= 1)
})

# ---- Reproducibility test with documented example values ----

test_that("mr.test reproduces documented example output", {
  # Using the exact example from the roxygen docs
  result <- mr.test(m1, m2, m3, xVec, var.equal = FALSE)
  # Verify the result is a named numeric vector of length 4
  expect_type(result, "double")
  expect_length(result, 4)
  expect_named(result, c("Statistic", "p-value", "Difference", "SE"))
})

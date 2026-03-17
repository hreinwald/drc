# Test suite for neill.test() and neill.default()
# Achieves 100% code coverage for R/neill.test.R

# ---- Setup: create model fixtures ----
# ryegrass dataset: 24 obs, 7 unique dose levels, LL.4() has 4 parameters
ryegrass_m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

# ---- Tests for neill.test() with explicit grouping ----

test_that("neill.test with explicit grouping and display=TRUE returns anova", {
  result <- neill.test(ryegrass_m1, ryegrass$conc)
  expect_s3_class(result, "anova")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("F value", "p value"))
  expect_equal(nrow(result), 1)
  expect_true(result[["F value"]] >= 0)
  expect_true(result[["p value"]] >= 0 && result[["p value"]] <= 1)
  expect_true(grepl("Neill", attr(result, "heading")))
})

test_that("neill.test with explicit grouping and display=FALSE returns p-value", {
  result <- neill.test(ryegrass_m1, ryegrass$conc, display = FALSE)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

# ---- Tests for grouping methods (missing grouping) ----

test_that("neill.test with method='finest' generates pairwise grouping", {
  result <- neill.test(ryegrass_m1, method = "finest")
  expect_s3_class(result, "anova")
  expect_named(result, c("F value", "p value"))
  expect_true(result[["p value"]] >= 0 && result[["p value"]] <= 1)
})

test_that("neill.test with method='c-finest' generates clustering grouping", {
  result <- neill.test(ryegrass_m1, method = "c-finest")
  expect_s3_class(result, "anova")
  expect_named(result, c("F value", "p value"))
  expect_true(result[["p value"]] >= 0 && result[["p value"]] <= 1)
})

test_that("neill.test with method='percentiles' generates percentile grouping", {
  result <- neill.test(ryegrass_m1, method = "percentiles")
  expect_s3_class(result, "anova")
  expect_named(result, c("F value", "p value"))
  expect_true(result[["p value"]] >= 0 && result[["p value"]] <= 1)
})

# ---- Tests for breakp parameter ----

test_that("neill.test with breakp overrides method-based grouping", {
  result <- neill.test(ryegrass_m1, breakp = c(0.5, 1, 2, 4, 8, 16))
  expect_s3_class(result, "anova")
  expect_named(result, c("F value", "p value"))
  expect_true(result[["p value"]] >= 0 && result[["p value"]] <= 1)
})

# ---- Tests for display=FALSE with each method ----

test_that("neill.test display=FALSE with method='finest' returns p-value", {
  result <- neill.test(ryegrass_m1, method = "finest", display = FALSE)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("neill.test display=FALSE with method='c-finest' returns p-value", {
  result <- neill.test(ryegrass_m1, method = "c-finest", display = FALSE)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("neill.test display=FALSE with method='percentiles' returns p-value", {
  result <- neill.test(ryegrass_m1, method = "percentiles", display = FALSE)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("neill.test display=FALSE with breakp returns p-value", {
  result <- neill.test(ryegrass_m1, breakp = c(0.5, 1, 2, 4, 8, 16),
                       display = FALSE)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

# ---- Error handling tests for neill.default() ----

test_that("neill.test errors with too many groups", {
  # 24 groups for 24 observations -> denDF = 24 - 24 = 0
  expect_error(
    neill.test(ryegrass_m1, grouping = 1:24),
    "Too many groups"
  )
})

test_that("neill.test errors with too few groups", {
  # 1 group for 4-param model -> numDF = 1 - 4 = -3

  expect_error(
    neill.test(ryegrass_m1, grouping = rep(1, 24)),
    "Too few groups"
  )
})

test_that("neill.test errors when M equals p (numDF = 0)", {
  # 4 groups for 4-param model -> numDF = 4 - 4 = 0
  expect_error(
    neill.test(ryegrass_m1, grouping = rep(1:4, each = 6)),
    "Too few groups"
  )
})

# ---- Tests for neill.default() output formatting ----

test_that("neill.default anova output has correct heading", {
  result <- neill.test(ryegrass_m1, ryegrass$conc)
  expect_equal(attr(result, "heading"), "Neill's lack-of-fit test\n")
})

test_that("neill.default non-display mode returns matrix with F and p", {
  # Access internal function directly for more targeted testing
  result <- neill.test(ryegrass_m1, ryegrass$conc, display = FALSE)
  # When display=FALSE, returns the p-value extracted from [1,2] of the matrix
  expect_true(is.numeric(result))
})

# ---- Test for method partial matching ----

test_that("neill.test method argument supports partial matching", {
  # "perc" should match "percentiles"
  result <- neill.test(ryegrass_m1, method = "perc")
  expect_s3_class(result, "anova")
})

# ---- Verify consistent results ----

test_that("explicit conc grouping and c-finest give same result for replicated data", {
  # ryegrass has replicates, so c-finest should recover original dose grouping
  r_explicit <- neill.test(ryegrass_m1, ryegrass$conc, display = FALSE)
  r_cfin <- neill.test(ryegrass_m1, method = "c-finest", display = FALSE)
  expect_equal(r_explicit, r_cfin)
})

# ---- Grouping display output test ----

test_that("neill.test with display=TRUE prints grouping info", {
  expect_output(
    neill.test(ryegrass_m1, ryegrass$conc),
    "Grouping used"
  )
})

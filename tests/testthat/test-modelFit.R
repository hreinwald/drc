# Tests for modelFit.R and all related functions
# Functions: modelFit, lofTest, gofTest, returnFct

# ==============================================================================
# returnFct tests
# ==============================================================================

test_that("returnFct with default arguments returns 'No test available' anova", {
  result <- drc:::returnFct()
  expect_s3_class(result, "anova")
  expect_s3_class(result, "data.frame")
  expect_equal(attr(result, "heading"), "No test available\n")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 5)
  expect_true(all(is.na(result)))
})

test_that("returnFct with provided arguments creates proper anova table", {
  result <- drc:::returnFct(
    dfModel = c(10, 15),
    loglik = c(100, 120),
    dfDiff = c(NA, 5),
    testStat = c(NA, 3.5),
    pVal = c(NA, 0.02),
    headName = "Test heading\n",
    colNames = c("A", "B", "C", "D", "E"),
    rowNames = c("Row1", "Row2")
  )
  expect_s3_class(result, "anova")
  expect_equal(attr(result, "heading"), "Test heading\n")
  expect_equal(rownames(result), c("Row1", "Row2"))
  expect_equal(colnames(result), c("A", "B", "C", "D", "E"))
  expect_equal(result[1, 1], 10)
  expect_equal(result[2, 1], 15)
  expect_equal(result[2, 4], 3.5)
  expect_equal(result[2, 5], 0.02)
})

# ==============================================================================
# modelFit tests: continuous data (F-test path via lofTest)
# ==============================================================================

test_that("modelFit with continuous data performs lack-of-fit F-test", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- modelFit(m1)

  expect_s3_class(result, "anova")
  expect_s3_class(result, "data.frame")
  expect_equal(attr(result, "heading"), "Lack-of-fit test\n")
  expect_equal(nrow(result), 2)
  expect_equal(rownames(result), c("ANOVA", "DRC model"))
  expect_equal(colnames(result), c("ModelDf", "RSS", "Df", "F value", "p value"))

  # ANOVA row should have ModelDf and RSS but NA for test stats
  expect_false(is.na(result[1, "ModelDf"]))
  expect_false(is.na(result[1, "RSS"]))
  expect_true(is.na(result[1, "Df"]))
  expect_true(is.na(result[1, "F value"]))
  expect_true(is.na(result[1, "p value"]))

  # DRC model row should have all values
  expect_false(is.na(result[2, "ModelDf"]))
  expect_false(is.na(result[2, "RSS"]))
  expect_false(is.na(result[2, "Df"]))
  expect_false(is.na(result[2, "F value"]))
  expect_false(is.na(result[2, "p value"]))

  # p-value should be valid
  pval <- result[2, "p value"]
  expect_true(pval >= 0 && pval <= 1)
})

test_that("modelFit with continuous data and Box-Cox uses bcAdd", {
  m_bc <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), bcVal = 0.5, bcAdd = 1)
  expect_false(is.null(m_bc$boxcox))
  expect_equal(m_bc$boxcox$bcAdd, 1)

  result <- modelFit(m_bc)
  expect_s3_class(result, "anova")
  expect_equal(attr(result, "heading"), "Lack-of-fit test\n")
  expect_false(is.na(result[2, "p value"]))
})

test_that("modelFit returns 'No test available' when ANOVA has 0 residual df", {
  # Create data with all unique doses (no replicates) so ANOVA has 0 df
  set.seed(42)
  dose_unique <- 1:20
  resp_unique <- 5 + 10 / (1 + (dose_unique / 5)^(-2)) + rnorm(20, 0, 0.5)
  df_unique <- data.frame(dose = dose_unique, resp = resp_unique)
  m_unique <- drm(resp ~ dose, data = df_unique, fct = LL.4())

  result <- modelFit(m_unique)
  expect_s3_class(result, "anova")
  expect_equal(attr(result, "heading"), "No test available\n")
})

test_that("modelFit method argument is validated", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  # Valid methods
  result_gof <- modelFit(m1, method = "gof")
  expect_s3_class(result_gof, "anova")
  result_cum <- modelFit(m1, method = "cum")
  expect_s3_class(result_cum, "anova")
  # Invalid method should error
  expect_error(modelFit(m1, method = "invalid"))
})

# ==============================================================================
# modelFit tests: binomial data (gofTest path)
# ==============================================================================

test_that("modelFit with binomial data performs goodness-of-fit test", {
  m_binom <- drm(r / n ~ dose, weights = n, data = deguelin, fct = LL.2(),
                 type = "binomial")
  result <- modelFit(m_binom)

  expect_s3_class(result, "anova")
  expect_equal(attr(result, "heading"), "Goodness-of-fit test\n")
  expect_equal(nrow(result), 2)

  # Check Chisq value and p value
  chisq_val <- result[2, "Chisq value"]
  pval <- result[2, "p value"]
  expect_true(!is.na(chisq_val))
  expect_true(!is.na(pval))
  expect_true(pval >= 0 && pval <= 1)
})

# ==============================================================================
# modelFit tests: Poisson data
# ==============================================================================

test_that("modelFit with Poisson type returns NULL", {
  set.seed(42)
  dose_p <- rep(c(0, 1, 2, 5, 10, 20), each = 3)
  count_p <- rpois(18, lambda = 50 * exp(-0.1 * dose_p))
  df_p <- data.frame(dose = dose_p, count = count_p)
  m_p <- drm(count ~ dose, data = df_p, fct = LL.4(), type = "Poisson")

  result <- modelFit(m_p)
  expect_null(result)
})

# ==============================================================================
# lofTest: direct tests for edge cases
# ==============================================================================

test_that("lofTest with NULL anovaTest returns 'No test available'", {
  result <- drc:::lofTest(object = NULL, anovaTest = NULL)
  expect_s3_class(result, "anova")
  expect_equal(attr(result, "heading"), "No test available\n")
})

test_that("lofTest with F-test and NaN test statistic returns NA p-value", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Create a mock anovaTest that returns a fit with 0 deviance and 0 residual df
  # This causes testStat = (nlsSS - 0) / dfDiff / (0 / 0) = NaN
  mock_anovaTest_nan <- function(formula, ds) {
    fit <- lm(formula, data = ds)
    # Create a mock fit where deviance = 0 and df.residual = 0
    mock_fit <- list(
      residuals = rep(0, nrow(ds)),
      rank = ncol(model.matrix(formula, data = ds)),
      fitted.values = fit$fitted.values,
      assign = fit$assign,
      df.residual = 0L  # Force 0 residual df
    )
    class(mock_fit) <- "lm"
    list(test = "F", anovaFit = mock_fit)
  }

  result <- drc:::lofTest(m1, mock_anovaTest_nan)
  expect_s3_class(result, "anova")
  # p-value should be NA when test stat is NaN
  expect_true(is.na(result[2, "p value"]))
})

test_that("lofTest with F-test and infinite test statistic returns NA p-value", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Create a mock anovaTest that returns a fit with 0 deviance but positive df
  # This causes testStat = (nlsSS - 0) / dfDiff / (0 / anovaDF) = Inf
  mock_anovaTest_inf <- function(formula, ds) {
    fit <- lm(formula, data = ds)
    mock_fit <- list(
      residuals = rep(0, nrow(ds)),
      rank = ncol(model.matrix(formula, data = ds)),
      fitted.values = fit$fitted.values,
      assign = fit$assign,
      df.residual = fit$df.residual  # Keep positive df
    )
    class(mock_fit) <- "lm"
    list(test = "F", anovaFit = mock_fit)
  }

  result <- drc:::lofTest(m1, mock_anovaTest_inf)
  expect_s3_class(result, "anova")
  # p-value should be NA when test stat is infinite
  expect_true(is.na(result[2, "p value"]))
})

test_that("lofTest with F-test and negative test statistic returns p-value of 1", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Create a mock anovaTest where anovaSS > nlsSS, producing negative F-stat
  # testStat = (nlsSS - anovaSS) / dfDiff / (anovaSS / anovaDF)
  # If anovaSS is very large, testStat will be negative
  mock_anovaTest_neg <- function(formula, ds) {
    fit <- lm(formula, data = ds)
    # Create a mock fit with very large residuals (large deviance)
    mock_fit <- fit
    mock_fit$residuals <- fit$residuals * 1000  # Makes deviance much larger than nlsSS
    list(test = "F", anovaFit = mock_fit)
  }

  result <- drc:::lofTest(m1, mock_anovaTest_neg)
  expect_s3_class(result, "anova")
  # p-value should be 1 for negative F-stat
  expect_equal(result[2, "p value"], 1)
})

test_that("lofTest lr test path works correctly", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Create a mock anovaTest that returns test = "lr"
  mock_anovaTest_lr <- function(formula, ds) {
    fit <- glm(formula, data = ds)
    list(test = "lr", anovaFit = fit)
  }

  result <- drc:::lofTest(m1, mock_anovaTest_lr)
  expect_s3_class(result, "anova")
  expect_equal(attr(result, "heading"), "Goodness-of-fit test\n")
  expect_equal(rownames(result), c("ANOVA", "DRC model"))
  expect_equal(colnames(result), c("ModelDf", "Log lik", "Df", "Chisq value", "p value"))
})

# ==============================================================================
# gofTest: direct tests
# ==============================================================================

test_that("gofTest with NULL gofTest result returns 'No test available'", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Create a mock gofTest function that returns NULL
  mock_gof <- function(resp, weights, fitted, dfres) NULL

  result <- drc:::gofTest(m1, mock_gof)
  expect_s3_class(result, "anova")
  expect_equal(attr(result, "heading"), "No test available\n")
})

test_that("gofTest with valid result returns goodness-of-fit table", {
  m_binom <- drm(r / n ~ dose, weights = n, data = deguelin, fct = LL.2(),
                 type = "binomial")

  # Use the actual gofTest from drmLOFbinomial
  binom_fns <- drc:::drmLOFbinomial()
  result <- drc:::gofTest(m_binom, binom_fns$gofTest)

  expect_s3_class(result, "anova")
  expect_equal(attr(result, "heading"), "Goodness-of-fit test\n")
})

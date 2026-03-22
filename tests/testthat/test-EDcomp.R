# Tests for EDcomp (R/EDcomp.R) and related functions:
#   - EDcomp: main function for comparing relative potencies
#   - fieller: Fieller's confidence interval (also tested in test-siInner.R)
#   - splitInd: split index vectors into shared/unique components
#   - createsifct: factory for selectivity index functions

# =============================================================================
# Tests for EDcomp
# =============================================================================

test_that("EDcomp basic call with numeric curve names returns correct structure", {
  # spinach has numeric curve names: 1,2,3,4,5
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(50, 50), display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("Estimate", "Std. Error", "t-value", "p-value"))
  # With 5 curves and 2 percentages: C(5,2) * C(2,2) = 10 * 1 = 10 comparisons
  expect_equal(nrow(result), 10)
})

test_that("EDcomp with non-numeric curve names triggers alphabetical ordering", {
  # Use HERBICIDE factor for non-numeric curve names
  spinach_herb <- drm(SLOPE ~ DOSE, HERBICIDE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach_herb, c(50, 50), display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)  # 2 curves -> 1 comparison
  expect_true(grepl("bentazon", rownames(result)[1]))
  expect_true(grepl("diuron", rownames(result)[1]))
})

test_that("EDcomp errors when interval='fls' and logBase is NULL", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  expect_error(
    EDcomp(spinach.LL.4, c(50, 50), interval = "fls"),
    "Argument 'logBase' not specified"
  )
})

test_that("EDcomp errors when relative percentages are outside (0, 100)", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  expect_error(
    EDcomp(spinach.LL.4, c(0, 50)),
    "Percentages outside the interval"
  )
  expect_error(
    EDcomp(spinach.LL.4, c(100, 50)),
    "Percentages outside the interval"
  )
  expect_error(
    EDcomp(spinach.LL.4, c(-5, 50)),
    "Percentages outside the interval"
  )
})

test_that("EDcomp with compMatch filters comparisons", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(50, 50), compMatch = c("1", "2"), display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_true(grepl("1/2", rownames(result)[1]))
})

test_that("EDcomp with percMat restricts percentage comparisons", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(10, 50, 90), percMat = matrix(c(1, 2), ncol = 2),
                   display = FALSE)

  expect_true(is.matrix(result))
  # Only one percentage comparison (10 vs 50), but 10 curve pairs
  expect_equal(nrow(result), 10)
})

test_that("EDcomp with reverse=TRUE reverses comparison order", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result_fwd <- EDcomp(spinach.LL.4, c(10, 50), display = FALSE)
  result_rev <- EDcomp(spinach.LL.4, c(10, 50), reverse = TRUE, display = FALSE)

  expect_true(is.matrix(result_rev))
  # Reversed order: ratio should be reciprocal
  # Row names should be reversed
  expect_true(grepl("2/1", rownames(result_rev)[1]))
})

test_that("EDcomp with interval='delta' returns CI columns", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(50, 50), interval = "delta", display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
  expect_equal(colnames(result), c("Estimate", "Lower", "Upper"))
  # Lower < Estimate < Upper for most comparisons
  expect_true(all(result[, "Lower"] < result[, "Estimate"]))
  expect_true(all(result[, "Upper"] > result[, "Estimate"]))
})

test_that("EDcomp with interval='fieller' returns CI columns", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(50, 50), interval = "fieller", display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
  expect_equal(colnames(result), c("Estimate", "Lower", "Upper"))
})

test_that("EDcomp with interval='fls' and logBase returns CI columns", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(50, 50), interval = "fls", logBase = 10,
                   display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
  expect_equal(colnames(result), c("Estimate", "Lower", "Upper"))
})

test_that("EDcomp with display=TRUE prints output", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  expect_output(
    EDcomp(spinach.LL.4, c(50, 50), display = TRUE),
    "Estimated ratios of effect doses"
  )
})

test_that("EDcomp with multcomp=TRUE returns parm object", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(50, 50), multcomp = TRUE, display = FALSE)

  expect_true(is.list(result))
  expect_true("multcomp" %in% names(result))
  expect_s3_class(result$multcomp, "parm")
})

test_that("EDcomp with multcomp=FALSE returns matrix invisibly", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(50, 50), multcomp = FALSE, display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 4)
})

test_that("EDcomp with logBase (no fls) applies logBase transformation", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(50, 50), logBase = 10, display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 4)
  # All estimates should be positive (10^x is always positive)
  expect_true(all(result[, "Estimate"] > 0))
})

test_that("EDcomp with 3 percentages generates correct number of comparisons", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(10, 50, 90), display = FALSE)

  # 5 curves -> C(5,2) = 10 curve pairs
  # 3 percentages -> C(3,2) = 3 percentage pairs
  # Total: 10 * 3 = 30 comparisons
  expect_equal(nrow(result), 30)
})

test_that("EDcomp switch statement covers 'delta' ciLabel path", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  # Delta interval triggers "Delta method" ciLabel
  expect_output(
    EDcomp(spinach.LL.4, c(50, 50), interval = "delta", display = TRUE),
    "Estimated ratios of effect doses"
  )
})

test_that("EDcomp switch statement covers 'fieller' ciLabel path", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  expect_output(
    EDcomp(spinach.LL.4, c(50, 50), interval = "fieller", display = TRUE),
    "Estimated ratios of effect doses"
  )
})

test_that("EDcomp switch statement covers 'fls' ciLabel path", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  expect_output(
    EDcomp(spinach.LL.4, c(50, 50), interval = "fls", logBase = 10, display = TRUE),
    "Estimated ratios of effect doses"
  )
})

test_that("EDcomp with compMatch that matches no curves returns empty matrix", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(50, 50), compMatch = c("nonexistent1", "nonexistent2"),
                   display = FALSE)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 0)
})

test_that("EDcomp with two percentages and two curves in compMatch", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  result <- EDcomp(spinach.LL.4, c(10, 50), compMatch = c("1", "2"),
                   display = FALSE)

  expect_true(is.matrix(result))
  # Only 1 curve pair, 1 percentage pair -> 1 row
  expect_equal(nrow(result), 1)
})

# =============================================================================
# Tests for splitInd
# =============================================================================

test_that("splitInd correctly identifies common and unique elements", {
  result <- drc:::splitInd(c(1, 2, 3, 4), c(3, 4, 5, 6))

  # only1: elements in ind1 but not ind2
  expect_equal(result[[1]], matrix(c(1, 2, 1, 2), ncol = 2))
  # only2: elements in ind2 but not ind1
  expect_equal(result[[2]], matrix(c(3, 4, 5, 6), ncol = 2))
  # inCommon: shared elements with positions in both vectors
  expect_true(is.matrix(result[[3]]))
  expect_equal(nrow(result[[3]]), 2)  # elements 3 and 4 are common
  expect_equal(result[[3]][, 3], c(3, 4))  # the common values
})

test_that("splitInd with no overlap returns NULL for inCommon", {
  result <- drc:::splitInd(c(1, 2), c(3, 4))

  expect_equal(result[[1]], matrix(c(1, 2, 1, 2), ncol = 2))
  expect_equal(result[[2]], matrix(c(1, 2, 3, 4), ncol = 2))
  expect_null(result[[3]])
})

test_that("splitInd with full overlap returns empty only1 and only2", {
  result <- drc:::splitInd(c(1, 2, 3), c(1, 2, 3))

  # only1 and only2 should have 0 rows
  expect_equal(nrow(result[[1]]), 0)
  expect_equal(nrow(result[[2]]), 0)
  # inCommon should have 3 rows
  expect_equal(nrow(result[[3]]), 3)
  expect_equal(result[[3]][, 3], c(1, 2, 3))
})

test_that("splitInd with single elements", {
  result <- drc:::splitInd(c(5), c(5))

  expect_equal(nrow(result[[1]]), 0)
  expect_equal(nrow(result[[2]]), 0)
  expect_equal(nrow(result[[3]]), 1)
  expect_equal(result[[3]][, 3], 5)
})

test_that("splitInd with single elements no overlap", {
  result <- drc:::splitInd(c(1), c(2))

  expect_equal(result[[1]], matrix(c(1, 1), ncol = 2))
  expect_equal(result[[2]], matrix(c(1, 2), ncol = 2))
  expect_null(result[[3]])
})

# =============================================================================
# Tests for createsifct
# =============================================================================

test_that("createsifct errors when edfct is NULL", {
  indexMat <- matrix(1:4, nrow = 2, ncol = 2)
  expect_error(
    drc:::createsifct(NULL, NULL, FALSE, indexMat, 4),
    "SI values cannot be calculated"
  )
})

test_that("createsifct returns function when fls=FALSE and logBase=NULL", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  sifct <- drc:::createsifct(
    spinach.LL.4$fct$edfct, logBase = NULL, fls = FALSE,
    indexMat = spinach.LL.4$indexMat, lenCoef = length(coef(spinach.LL.4))
  )
  expect_true(is.function(sifct))

  # Test the returned function
  parm1 <- spinach.LL.4$parmMat[, 1]
  parm2 <- spinach.LL.4$parmMat[, 2]
  result <- sifct(parm1, parm2, c(50, 50), 1, 2, "control", "relative")
  expect_true(is.list(result))
  expect_true("val" %in% names(result))
  expect_true("der" %in% names(result))
  expect_true("der1" %in% names(result))
  expect_true("der2" %in% names(result))
  expect_true("valnum" %in% names(result))
  expect_true("valden" %in% names(result))
  # val should be ratio of ED values
  expect_true(is.numeric(result$val))
})

test_that("createsifct returns function when fls=FALSE and logBase is provided", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  sifct <- drc:::createsifct(
    spinach.LL.4$fct$edfct, logBase = 10, fls = FALSE,
    indexMat = spinach.LL.4$indexMat, lenCoef = length(coef(spinach.LL.4))
  )
  expect_true(is.function(sifct))

  parm1 <- spinach.LL.4$parmMat[, 1]
  parm2 <- spinach.LL.4$parmMat[, 2]
  result <- sifct(parm1, parm2, c(50, 50), 1, 2, "control", "relative")
  expect_true(is.list(result))
  # With logBase, SIpair = logBase^(ED1v - ED2v)
  expect_true(result$val > 0)  # logBase^x is always positive
})

test_that("createsifct returns function when fls=TRUE", {
  spinach.LL.4 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4())
  sifct <- drc:::createsifct(
    spinach.LL.4$fct$edfct, logBase = 10, fls = TRUE,
    indexMat = spinach.LL.4$indexMat, lenCoef = length(coef(spinach.LL.4))
  )
  expect_true(is.function(sifct))

  parm1 <- spinach.LL.4$parmMat[, 1]
  parm2 <- spinach.LL.4$parmMat[, 2]
  result <- sifct(parm1, parm2, c(50, 50), 1, 2, "control", "relative")
  expect_true(is.list(result))
  # With fls=TRUE, SIpair = ED1v - ED2v (difference, not ratio)
  expect_true(is.numeric(result$val))
})

# =============================================================================
# Tests for fieller (additional coverage beyond test-siInner.R)
# =============================================================================

test_that("fieller standard (finney=FALSE) returns two numeric values", {
  mu <- c(10, 5)
  df <- 20
  vcMat <- matrix(c(1, 0.2, 0.2, 0.5), 2, 2)

  result <- drc:::fieller(mu, df, vcMat, level = 0.95)
  expect_true(is.numeric(result))
  expect_length(result, 2)
  expect_true(result[1] < mu[1] / mu[2])
  expect_true(result[2] > mu[1] / mu[2])
})

test_that("fieller Finney variant (finney=TRUE) returns two numeric values", {
  mu <- c(10, 5)
  df <- 20
  vcMat <- matrix(c(1, 0.2, 0.2, 0.5), 2, 2)
  resVar <- 2.0

  result <- drc:::fieller(mu, df, vcMat, level = 0.95, finney = TRUE, resVar = resVar)
  expect_true(is.numeric(result))
  expect_length(result, 2)
  expect_true(result[1] < result[2])
})

test_that("fieller finney=TRUE errors when g >= 1", {
  mu <- c(10, 0.5)
  df <- 20
  vcMat <- matrix(c(1, 0.2, 0.2, 50), 2, 2)
  resVar <- 2.0

  expect_error(
    drc:::fieller(mu, df, vcMat, level = 0.95, finney = TRUE, resVar = resVar),
    "Fieller's theorem not useful"
  )
})

test_that("fieller with different confidence levels", {
  mu <- c(10, 5)
  df <- 20
  vcMat <- matrix(c(1, 0.2, 0.2, 0.5), 2, 2)

  result_90 <- drc:::fieller(mu, df, vcMat, level = 0.90)
  result_99 <- drc:::fieller(mu, df, vcMat, level = 0.99)

  # 99% CI should be wider than 90% CI
  width_90 <- result_90[2] - result_90[1]
  width_99 <- result_99[2] - result_99[1]
  expect_true(width_99 > width_90)
})

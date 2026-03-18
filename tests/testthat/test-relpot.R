# tests/testthat/test-relpot.R
# Comprehensive tests for relpot() and commatFct()

# ============================================================
# Setup: create model fixtures used across multiple tests
# ============================================================

# Decreasing 2-curve model (monoton < 0): spinach data
data(spinach, package = "drc")
m_dec <- drm(SLOPE ~ DOSE, HERBICIDE, data = spinach, fct = LL.4())

# Increasing 2-curve model (monoton > 0): synthetic data
set.seed(42)
inc_data <- rbind(
  data.frame(
    dose = rep(c(0, 0.5, 1, 2, 5, 10), each = 3),
    resp = c(rnorm(3, 0.5, 0.1), rnorm(3, 1, 0.1), rnorm(3, 2, 0.1),
             rnorm(3, 3, 0.1), rnorm(3, 4.5, 0.1), rnorm(3, 5, 0.1)),
    group = "A"
  ),
  data.frame(
    dose = rep(c(0, 0.5, 1, 2, 5, 10), each = 3),
    resp = c(rnorm(3, 0.3, 0.1), rnorm(3, 0.8, 0.1), rnorm(3, 1.5, 0.1),
             rnorm(3, 2.5, 0.1), rnorm(3, 3.5, 0.1), rnorm(3, 4.8, 0.1)),
    group = "B"
  )
)
m_inc <- drm(resp ~ dose, group, data = inc_data, fct = LL.4())

# ============================================================
# commatFct tests
# ============================================================

test_that("commatFct returns full parmMat when compMatch is NULL", {
  result <- drc:::commatFct(m_dec, NULL)
  expect_equal(result, m_dec$parmMat)
})

test_that("commatFct filters columns when compMatch is specified", {
  result <- drc:::commatFct(m_dec, c("bentazon", "diuron"))
  expect_equal(ncol(result), 2)
  expect_true(all(c("bentazon", "diuron") %in% colnames(result)))
})

# ============================================================
# relpot: return structure (happy path)
# ============================================================

test_that("relpot returns a list with x, y, percVec components (plotit=FALSE)", {
  r <- relpot(m_dec, plotit = FALSE)
  expect_type(r, "list")
  expect_named(r, c("x", "y", "percVec"))
  expect_true(is.numeric(r$x))
  expect_true(is.numeric(r$y))
  expect_true(is.numeric(r$percVec))
  # Default type="relative", scale="original" gives 99 x-values
  expect_length(r$x, 99)
  expect_length(r$y, 99)
  # No NAs in the output
  expect_false(any(is.na(r$y)))
})

# ============================================================
# relpot: type="relative" with all three scale options
# ============================================================

test_that("relpot with type='relative', scale='original' (default)", {
  r <- relpot(m_dec, plotit = FALSE, type = "relative", scale = "original")
  expect_length(r$x, 99)
  expect_length(r$y, 99)
  # percVec should have 99 elements (101 minus endpoints)
  expect_length(r$percVec, 99)
})

test_that("relpot with type='relative', scale='percent'", {
  r <- relpot(m_dec, plotit = FALSE, type = "relative", scale = "percent")
  expect_length(r$x, 99)
  # xVec should equal percVec when scale is "percent"
  expect_equal(r$x, r$percVec)
})

test_that("relpot with type='relative', scale='unconstrained'", {
  r <- relpot(m_dec, plotit = FALSE, type = "relative", scale = "unconstrained")
  # unconstrained overrides percVec to 1:99
  expect_equal(r$percVec, 1:99)
  expect_length(r$y, 99)
  # xVec should equal percVec for unconstrained
  expect_equal(r$x, r$percVec)
})

# ============================================================
# relpot: type="absolute"
# ============================================================

test_that("relpot with type='absolute' auto-generates percVec", {
  r <- relpot(m_dec, plotit = FALSE, type = "absolute")
  expect_length(r$percVec, 100)
  expect_length(r$y, 100)
  # xVec should equal percVec for absolute type
  expect_equal(r$x, r$percVec)
})

test_that("relpot with type='absolute' and custom percVec", {
  # Use percVec values within the model's response range to avoid NaN
  parmMat <- m_dec$parmMat
  low <- max(apply(parmMat, 2, m_dec$fct$lowerAs))
  up  <- min(apply(parmMat, 2, m_dec$fct$upperAs))
  custom_perc <- seq(low * 1.1, up * 0.9, length.out = 20)
  r <- relpot(m_dec, plotit = FALSE, type = "absolute", percVec = custom_perc)
  expect_equal(r$percVec, custom_perc)
  expect_length(r$y, 20)
  expect_equal(r$x, custom_perc)
})

# ============================================================
# relpot: monoton < 0 vs monoton > 0 (decreasing vs increasing)
# ============================================================

test_that("relpot handles decreasing curves (monoton < 0)", {
  # m_dec has monoton < 0
  expect_true(m_dec$fct$monoton(m_dec$parmMat[, 1]) < 0)
  r <- relpot(m_dec, plotit = FALSE)
  expect_false(any(is.na(r$y)))
})

test_that("relpot handles increasing curves (monoton > 0)", {
  # m_inc has monoton > 0
  expect_true(m_inc$fct$monoton(m_inc$parmMat[, 1]) > 0)
  r <- relpot(m_inc, plotit = FALSE)
  expect_false(any(is.na(r$y)))
})

# ============================================================
# relpot: interval != "none" (confidence interval path)
# ============================================================

test_that("relpot with interval='delta' computes confidence bands", {
  custom_perc <- seq(40, 60, by = 10)
  r <- relpot(m_dec, plotit = FALSE, interval = "delta", percVec = custom_perc)
  expect_type(r, "list")
  expect_length(r$y, 3)
  expect_false(any(is.na(r$y)))
})

# ============================================================
# relpot: custom percVec with type="relative"
# ============================================================

test_that("relpot with custom percVec skips auto-determination", {
  custom_perc <- c(30, 40, 50, 60, 70)
  r <- relpot(m_dec, plotit = FALSE, percVec = custom_perc)
  expect_equal(r$percVec, custom_perc)
  expect_length(r$y, 5)
})

# ============================================================
# relpot: plotit=TRUE tests (all plot code paths)
# ============================================================

test_that("relpot plots with type='relative', scale='original', interval='none'", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  r <- relpot(m_dec, plotit = TRUE, type = "relative", scale = "original")
  expect_type(r, "list")
})

test_that("relpot plots with type='relative', scale='percent'", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  r <- relpot(m_dec, plotit = TRUE, type = "relative", scale = "percent")
  expect_type(r, "list")
})

test_that("relpot plots with type='relative', scale='unconstrained'", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  r <- relpot(m_dec, plotit = TRUE, type = "relative", scale = "unconstrained")
  expect_type(r, "list")
})

test_that("relpot plots with type='absolute', interval='none'", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  r <- relpot(m_dec, plotit = TRUE, type = "absolute")
  expect_type(r, "list")
})

test_that("relpot plots with interval='delta' (confidence bands)", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  # Use auto-generated percVec for consistency
  r <- relpot(m_dec, plotit = TRUE, interval = "delta")
  expect_type(r, "list")
})

test_that("relpot plots with interval='delta', type='absolute' (no reference line)", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  r <- relpot(m_dec, plotit = TRUE, interval = "delta", type = "absolute")
  expect_type(r, "list")
})

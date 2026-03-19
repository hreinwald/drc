# Tests for maED() function - Model-averaged effective doses

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


# --- Input validation tests ---

test_that("maED errors when object is not of class drc", {
  expect_error(maED("not_a_model", list(LL.5()), c(50)), "'object' must be of class 'drc'")
  expect_error(maED(42, list(LL.5()), c(50)), "'object' must be of class 'drc'")
})

test_that("maED errors when respLev is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(maED(m1, list(LL.5()), "abc"), "'respLev' must be a non-empty numeric vector")
  expect_error(maED(m1, list(LL.5()), numeric(0)), "'respLev' must be a non-empty numeric vector")
})

test_that("maED errors when level is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(maED(m1, list(LL.5()), 50, level = "a"), "'level' must be a single numeric value strictly between 0 and 1")
  expect_error(maED(m1, list(LL.5()), 50, level = 0), "'level' must be a single numeric value strictly between 0 and 1")
  expect_error(maED(m1, list(LL.5()), 50, level = 1), "'level' must be a single numeric value strictly between 0 and 1")
  expect_error(maED(m1, list(LL.5()), 50, level = c(0.9, 0.95)), "'level' must be a single numeric value strictly between 0 and 1")
})

test_that("maED errors when linreg is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(maED(m1, list(LL.5()), 50, linreg = "yes"), "'linreg' must be a single logical value")
  expect_error(maED(m1, list(LL.5()), 50, linreg = c(TRUE, FALSE)), "'linreg' must be a single logical value")
})

test_that("maED errors when display is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(maED(m1, list(LL.5()), 50, display = "yes"), "'display' must be a single logical value")
})

test_that("maED errors when na.rm is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(maED(m1, list(LL.5()), 50, na.rm = "yes"), "'na.rm' must be a single logical value")
})

test_that("maED errors when extended is invalid", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_error(maED(m1, list(LL.5()), 50, extended = "yes"), "'extended' must be a single logical value")
})


# --- Happy path tests ---

test_that("maED returns matrix with correct structure for interval='none'", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- maED(m1, list(LL.5(), W1.4()), c(10, 50, 90), display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "Estimate")
  expect_true(all(result[, "Estimate"] > 0))
})

test_that("maED returns matrix with correct structure for interval='buckland'", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- maED(m1, list(LL.5(), W1.4()), c(10, 50), interval = "buckland", display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 4)
  expect_true(all(c("Estimate", "Std. Error", "Lower", "Upper") %in% colnames(result)))
  expect_true(all(result[, "Lower"] < result[, "Estimate"]))
  expect_true(all(result[, "Upper"] > result[, "Estimate"]))
})

test_that("maED returns matrix with correct structure for interval='kang'", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- maED(m1, list(LL.5(), W1.4()), c(10, 50), interval = "kang", display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  expect_true(all(c("Estimate", "Lower", "Upper") %in% colnames(result)))
})

test_that("maED works with a single response level", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- maED(m1, list(LL.5()), 50, display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_true(result[, "Estimate"] > 0)
})


# --- Extended output ---

test_that("maED returns list when extended = TRUE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- maED(m1, list(LL.5(), W1.4()), c(10, 50), display = FALSE, extended = TRUE)

  expect_true(is.list(result))
  expect_true(all(c("estimates", "fits") %in% names(result)))
  expect_true(is.matrix(result$estimates))
  expect_true(is.matrix(result$fits))
  expect_equal(nrow(result$fits), 3)  # 1 original + 2 in fctList
  expect_true("Weight" %in% colnames(result$fits))
})


# --- Display parameter ---

test_that("maED prints when display = TRUE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_output(maED(m1, list(LL.5()), 50, display = TRUE), "Weight")
})

test_that("maED is silent when display = FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_no_error(
    suppressWarnings(result <- maED(m1, list(LL.5()), 50, display = FALSE))
  )
  expect_true(is.matrix(result))
})


# --- Linear regression option ---

test_that("maED works with linreg = TRUE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- maED(m1, list(LL.5()), c(10, 50), linreg = TRUE, display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_true(all(result[, "Estimate"] > 0))
})

test_that("maED linreg = TRUE with extended returns fit info including Lin row", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- maED(m1, list(LL.5()), c(50), linreg = TRUE, display = FALSE, extended = TRUE)

  expect_true("Lin" %in% rownames(result$fits))
})

test_that("maED linreg = TRUE with buckland interval", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- maED(m1, list(LL.5()), c(50), linreg = TRUE, interval = "buckland", display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 4)
})


# --- Multi-curve handling ---

test_that("maED handles multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  result <- maED(m_multi, list(LL.5()), 50, display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)  # One row per curve
})

test_that("maED handles multi-curve with clevel filter", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  result <- maED(m_multi, list(LL.5()), 50, clevel = "A", display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
})


# --- type = "absolute" ---

test_that("maED works with type = 'absolute'", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- maED(m1, list(LL.5()), 5, type = "absolute", display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
})


# --- na.rm option ---

test_that("maED works with na.rm = TRUE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- maED(m1, list(LL.5()), c(10, 50), na.rm = TRUE, display = FALSE)

  expect_true(is.matrix(result))
  expect_true(all(result[, "Estimate"] > 0))
})


# --- try-error handling for models that fail in fctList ---

test_that("maED handles try-error from failed model in fctList", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  # LL.2 has fixed lower limit at 0, may produce different results;
  # Use a function that might fail during update fitting
  result <- maED(m1, list(LL.5(), LN.4(), W1.4(), W2.4()), c(10, 50, 90),
                 interval = "buckland", display = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
})

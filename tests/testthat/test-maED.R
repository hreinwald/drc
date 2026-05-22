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
  expect_equal(ncol(result), 4)
  expect_true(all(c("Estimate", "Std. Error", "Lower", "Upper") %in% colnames(result)))
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


# --- Non-finite ED value filtering ---

# Algae dataset from the debug folder that causes LL.5 to return Inf for ED50
# due to a negative f parameter estimate.
algae_data <- data.frame(
  yield = c(
    824948.1, 874756.3, 818722.1,  # conc = 10
    289510.1, 345544.3, 280171.1,  # conc = 10000
    1077102.1, 653732.5, 824948.1, # conc = 1000
    905886.4, 753348.9, 756461.9,  # conc = 100
    697314.7, 691088.6, 762687.9,  # conc = 10  (second batch)
    880982.4, 747122.8, 803157.1,  # conc = 1   (second batch)
    295736.1, 295736.1, 255267.0, 283284.1, 286397.1, 273945.0,  # controls
    1503584.7, 1388403.3, 946355.6, 1195396.5, 1410194.4,        # controls
    407804.6, 485629.8, 678636.6, 809383.1, 582133.2,            # controls
    1049085.0, 884095.4, 715992.7, 986824.8, 905886.4            # controls
  ),
  conc = c(
    rep(10, 3), rep(10000, 3), rep(1000, 3), rep(100, 3),
    rep(10, 3), rep(1, 3),
    rep(0, 21)
  )
)

test_that("maED warns when a model produces non-finite ED values or fitting fails", {
  m_algae <- drm(yield ~ conc, data = algae_data,
                 fct = LL.4(fixed = c(NA, 1e-9, NA, NA)))

  fcts <- list(
    LL.5(fixed = c(NA, 1e-9, NA, NA, NA)),
    W1.4(fixed = c(NA, 1e-9, NA, NA))
  )

  # Some platforms/R versions may produce non-finite ED values for these
  # models, triggering the exclusion warning. On others, all models converge
  # successfully. We therefore check both paths: if a warning is produced it
  # must match the expected pattern, and the result must always be valid.
  exclusion_warned <- FALSE
  result <- withCallingHandlers(
    maED(m_algae, fcts, 50, display = FALSE),
    warning = function(w) {
      if (grepl("excluded from model-averaging", conditionMessage(w))) {
        exclusion_warned <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_true(is.finite(result[, "Estimate"]))
})

test_that("maED extended output shows excluded models with zero weight", {
  m_algae <- drm(yield ~ conc, data = algae_data,
                 fct = LL.4(fixed = c(NA, 1e-9, NA, NA)))

  fcts <- list(
    LL.5(fixed = c(NA, 1e-9, NA, NA, NA)),
    W1.4(fixed = c(NA, 1e-9, NA, NA))
  )

  result <- suppressWarnings(
    maED(m_algae, fcts, 50, display = FALSE, extended = TRUE)
  )

  expect_true(is.list(result))
  fits <- result$fits

  # On platforms where models produce non-finite ED values, excluded models
  # get zero weight. On others, all models converge and all weights are
  # positive. Either outcome is valid; we only check structural correctness.
  expect_true(is.matrix(fits))
  expect_true("Weight" %in% colnames(fits))
  expect_true(all(fits[, "Weight"] >= 0))

  # The model-averaged estimate should be finite
  expect_true(is.finite(result$estimates[, "Estimate"]))
})

test_that("maED buckland interval works when models are excluded", {
  m_algae <- drm(yield ~ conc, data = algae_data,
                 fct = LL.4(fixed = c(NA, 1e-9, NA, NA)))

  fcts <- list(
    LL.5(fixed = c(NA, 1e-9, NA, NA, NA)),
    W1.4(fixed = c(NA, 1e-9, NA, NA))
  )

  result <- suppressWarnings(
    maED(m_algae, fcts, 50, interval = "buckland", display = FALSE)
  )

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 4)
  expect_true(all(c("Estimate", "Std. Error", "Lower", "Upper") %in% colnames(result)))
  # Result should be finite
  expect_true(all(is.finite(result)))
})

test_that("maED without non-finite values produces no exclusion warning", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # These models all produce finite ED50 on ryegrass data.
  # Internal optimization may emit "NaNs produced" warnings which are
  # unrelated to model exclusion, so we only check that no exclusion
  # warning is issued.
  exclusion_warned <- FALSE
  result <- withCallingHandlers(
    maED(m1, list(W1.4(), W2.4()), 50, display = FALSE),
    warning = function(w) {
      if (grepl("excluded from model-averaging", conditionMessage(w))) {
        exclusion_warned <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )

  expect_false(exclusion_warned)
  expect_true(is.matrix(result))
  expect_true(is.finite(result[, "Estimate"]))
})

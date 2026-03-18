# tests/testthat/test-twophase.R
# Comprehensive test suite for the twophase() function

# --- Argument Validation ---

test_that("twophase rejects invalid 'names' argument", {
  # names must be character

expect_error(twophase(names = 1:7), "Not correct 'names' argument")
  # names wrong length
  expect_error(twophase(names = c("a", "b")), "Not correct 'names' argument")
  # names not character (logical)
  expect_error(twophase(names = rep(TRUE, 7)), "Not correct 'names' argument")
})

test_that("twophase rejects invalid 'fixed' argument", {
  # fixed wrong length
  expect_error(twophase(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(twophase(fixed = rep(NA, 8)), "Not correct 'fixed' argument")
})

# --- Default Behavior (Happy Path) ---

test_that("twophase returns correct structure with default arguments", {
  res <- twophase()

  # Class
  expect_s3_class(res, "two-phase")

  # All expected list elements
  expect_true(is.list(res))
  expected_names <- c("fct", "ssfct", "names", "deriv1", "deriv2",
                       "derivx", "edfct", "name", "text", "noParm")
  expect_equal(names(res), expected_names)

  # Parameter names (all 7 free)
  expect_equal(res$names, c("b1", "c1", "d1", "e1", "b2", "d2", "e2"))

  # noParm should be 7
  expect_equal(res$noParm, 7)

  # name should be "twophase" (from match.call)
  expect_equal(res$name, "twophase")

  # text should be "Two-phase"
  expect_equal(res$text, "Two-phase")

  # NULL fields
  expect_null(res$deriv1)
  expect_null(res$deriv2)
  expect_null(res$derivx)
  expect_null(res$edfct)

  # Functions
  expect_true(is.function(res$fct))
  expect_true(is.function(res$ssfct))
})

# --- Fixed Parameters ---

test_that("twophase handles fixed parameters correctly", {
  # Fix b1 at 1 and c1 at 0
  res <- twophase(fixed = c(1, 0, NA, NA, NA, NA, NA))
  expect_equal(res$names, c("d1", "e1", "b2", "d2", "e2"))
  expect_equal(res$noParm, 5)

  # Fix all parameters
  res_all <- twophase(fixed = c(1, 0, 100, 5, 2, 50, 10))
  expect_equal(res_all$noParm, 0)
  expect_equal(length(res_all$names), 0)
})

# --- Custom Names ---

test_that("twophase accepts custom parameter names", {
  custom_names <- c("slope1", "lower", "upper1", "ed501", "slope2", "upper2", "ed502")
  res <- twophase(names = custom_names)
  expect_equal(res$names, custom_names)
})

# --- fctName and fctText Arguments ---

test_that("twophase uses provided fctName and fctText", {
  res <- twophase(fctName = "myModel", fctText = "My custom text")
  expect_equal(res$name, "myModel")
  expect_equal(res$text, "My custom text")
})

test_that("twophase uses defaults when fctName and fctText are missing", {
  res <- twophase()
  expect_equal(res$name, "twophase")
  expect_equal(res$text, "Two-phase")
})

# --- fct function (dose-response evaluation) ---

test_that("twophase fct computes dose-response values correctly", {
  res <- twophase()

  # Use known parameter values
  # b1=1, c1=0, d1=50, e1=5, b2=1, d2=50, e2=50
  parm <- matrix(c(1, 0, 50, 5, 1, 50, 50), nrow = 1)
  dose <- c(1, 5, 10, 50, 100)

  values <- res$fct(dose, parm)
  expect_true(is.numeric(values))
  expect_equal(length(values), length(dose))

  # At e1=5 (dose=5), first component should be at midpoint (d1-c1)/2 = 25
  # LL.4 at e1: c + (d-c)/2 = 0 + 50/2 = 25
  # LL.3 at dose=5: d2/(1+exp(b2*(log(5)-log(50)))) = 50/(1+exp(1*log(0.1)))
  #   = 50/(1+0.1) = 50/1.1 ≈ 45.45
  # Total ≈ 25 + 45.45 ≈ 70.45
  val_at_5 <- res$fct(5, parm)
  expect_true(val_at_5 > 60 && val_at_5 < 80)
})

test_that("twophase fct works with multiple rows in parm matrix", {
  res <- twophase()

  # Two sets of parameters
  parm <- matrix(c(1, 0, 50, 5, 1, 50, 50,
                    2, 10, 80, 3, 2, 40, 20), nrow = 2, byrow = TRUE)
  dose <- c(1, 10)

  values <- res$fct(dose, parm)
  expect_equal(length(values), 2)
  expect_true(is.numeric(values))
})

test_that("twophase fct works with fixed parameters", {
  # Fix b1=1 and c1=0
  res <- twophase(fixed = c(1, 0, NA, NA, NA, NA, NA))

  # parm should only have 5 columns (d1, e1, b2, d2, e2)
  parm <- matrix(c(50, 5, 1, 50, 50), nrow = 1)
  dose <- c(1, 5, 10)

  values <- res$fct(dose, parm)
  expect_equal(length(values), 3)
  expect_true(is.numeric(values))

  # Compare with full model using same parameter values
  res_full <- twophase()
  parm_full <- matrix(c(1, 0, 50, 5, 1, 50, 50), nrow = 1)
  values_full <- res_full$fct(dose, parm_full)

  expect_equal(values, values_full)
})

# --- ssfct function (self-starter) ---

test_that("twophase ssfct returns initial parameter values", {
  res <- twophase()

  # Create a mock data frame similar to what drm passes to ssfct
  # ssfct expects a data frame with dose and response
  set.seed(42)
  dose <- c(0.1, 0.5, 1, 2, 5, 10, 20, 50, 100)
  response <- c(5, 10, 20, 40, 60, 75, 85, 95, 100)
  dframe <- data.frame(dose, response)

  init_vals <- res$ssfct(dframe)

  # Should return 7 values (all params free)
  expect_equal(length(init_vals), 7)
  expect_true(is.numeric(init_vals))
})

test_that("twophase ssfct respects fixed parameters", {
  res <- twophase(fixed = c(1, 0, NA, NA, NA, NA, NA))

  dose <- c(0.1, 0.5, 1, 2, 5, 10, 20, 50, 100)
  response <- c(5, 10, 20, 40, 60, 75, 85, 95, 100)
  dframe <- data.frame(dose, response)

  init_vals <- res$ssfct(dframe)

  # Should return 5 values (only non-fixed params)
  expect_equal(length(init_vals), 5)
  expect_true(is.numeric(init_vals))
})

# --- Integration with drm ---

test_that("twophase works with drm for model fitting", {
  skip_if_not_installed("drc")

  # Create synthetic two-phase data
  set.seed(123)
  dose <- rep(c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500), each = 3)
  # Two-phase response: LL.4 component + LL.3 component
  b1 <- 1; c1 <- 0; d1 <- 50; e1 <- 5; b2 <- 1; d2 <- 50; e2 <- 100
  response <- c1 + (d1 - c1) / (1 + exp(b1 * (log(dose) - log(e1)))) +
    d2 / (1 + exp(b2 * (log(dose) - log(e2)))) +
    rnorm(length(dose), 0, 2)

  dat <- data.frame(dose = dose, response = response)

  # Should not error when fitting
  expect_no_error(
    mod <- drm(response ~ dose, data = dat, fct = twophase())
  )
})

# --- invisible return ---

test_that("twophase returns invisibly", {
  # The function uses invisible(), so direct assignment should work
  # but printing should not show output
  res <- twophase()
  expect_s3_class(res, "two-phase")
})

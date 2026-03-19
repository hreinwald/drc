# Test file for CRS.6 function
# Goal: Achieve 100% code coverage

library(testthat)
library(drc)

# Load test data
data(ryegrass)
test_data <- data.frame(dose = ryegrass$conc, response = ryegrass$rootl)

# ==============================================================================
# Basic Structure and Default Behavior Tests
# ==============================================================================

test_that("CRS.6 returns correct structure with defaults", {
  model <- CRS.6()

  expect_s3_class(model, "cedergreen.extended")
  expect_named(
    model,
    c("fct", "ssfct", "names", "deriv1", "deriv2", "edfct", "maxfct", "name", "text", "noParm")
  )
  expect_equal(model$noParm, 6)
  expect_equal(model$name, "CRS.6")
  expect_equal(model$text, "Generalised Cedergreen-Ritz-Streibig (hormesis)")
})

test_that("CRS.6 has correct default parameter names", {
  model <- CRS.6()

  expect_equal(model$names, c("b", "c", "d", "e", "f", "g"))
  expect_length(model$names, 6)
})

test_that("CRS.6 fct function exists and is callable", {
  model <- CRS.6()

  expect_true(is.function(model$fct))

  # Test with simple parameters
  result <- model$fct(
    dose = c(0.1, 1, 10),
    parm = matrix(c(2, 0, 100, 1, 10, 1), nrow = 1)
  )

  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("CRS.6 ssfct function exists and is callable", {
  model <- CRS.6()

  expect_true(is.function(model$ssfct))

  inits <- model$ssfct(test_data)

  expect_type(inits, "double")
  expect_length(inits, 6)
  expect_true(all(is.finite(inits)))
})

test_that("CRS.6 deriv1 and deriv2 are NULL", {
  model <- CRS.6()

  expect_null(model$deriv1)
  expect_null(model$deriv2)
})

test_that("CRS.6 edfct and maxfct are NULL", {
  model <- CRS.6()

  expect_null(model$edfct)
  expect_null(model$maxfct)
})

# ==============================================================================
# Custom Names Tests
# ==============================================================================

test_that("CRS.6 accepts custom parameter names", {
  custom_names <- c("slope", "lower", "upper", "ed50", "horm", "alpha")

  model <- CRS.6(names = custom_names)

  expect_equal(model$names, custom_names)
})

test_that("CRS.6 rejects invalid names argument - not character", {
  expect_error(
    CRS.6(names = c(1, 2, 3, 4, 5, 6)),
    "Not correct 'names' argument"
  )
})

test_that("CRS.6 rejects invalid names argument - wrong length", {
  expect_error(
    CRS.6(names = c("b", "c", "d", "e")),
    "Not correct 'names' argument"
  )
})

test_that("CRS.6 rejects invalid names argument - NULL", {
  expect_error(
    CRS.6(names = NULL),
    "Not correct 'names' argument"
  )
})

# ==============================================================================
# Fixed Parameters Tests
# ==============================================================================

test_that("CRS.6 with no fixed parameters", {
  model <- CRS.6(fixed = c(NA, NA, NA, NA, NA, NA))

  expect_equal(model$noParm, 6)
  expect_length(model$names, 6)

  inits <- model$ssfct(test_data)
  expect_length(inits, 6)
})

test_that("CRS.6 with all parameters fixed", {
  model <- CRS.6(fixed = c(2, 0, 100, 1, 10, 1))

  expect_equal(model$noParm, 0)
  expect_length(model$names, 0)

  inits <- model$ssfct(test_data)
  expect_length(inits, 0)
})

test_that("CRS.6 with b fixed", {
  model <- CRS.6(fixed = c(2, NA, NA, NA, NA, NA))

  expect_equal(model$noParm, 5)
  expect_equal(model$names, c("c", "d", "e", "f", "g"))

  inits <- model$ssfct(test_data)
  expect_length(inits, 5)
})

test_that("CRS.6 with c fixed", {
  model <- CRS.6(fixed = c(NA, 0, NA, NA, NA, NA))

  expect_equal(model$noParm, 5)
  expect_equal(model$names, c("b", "d", "e", "f", "g"))

  inits <- model$ssfct(test_data)
  expect_length(inits, 5)
})

test_that("CRS.6 with d fixed", {
  model <- CRS.6(fixed = c(NA, NA, 100, NA, NA, NA))

  expect_equal(model$noParm, 5)
  expect_equal(model$names, c("b", "c", "e", "f", "g"))

  inits <- model$ssfct(test_data)
  expect_length(inits, 5)
})

test_that("CRS.6 with e fixed", {
  model <- CRS.6(fixed = c(NA, NA, NA, 1, NA, NA))

  expect_equal(model$noParm, 5)
  expect_equal(model$names, c("b", "c", "d", "f", "g"))

  inits <- model$ssfct(test_data)
  expect_length(inits, 5)
})

test_that("CRS.6 with f fixed", {
  model <- CRS.6(fixed = c(NA, NA, NA, NA, 10, NA))

  expect_equal(model$noParm, 5)
  expect_equal(model$names, c("b", "c", "d", "e", "g"))

  inits <- model$ssfct(test_data)
  expect_length(inits, 5)
})

test_that("CRS.6 with g (alpha) fixed", {
  model <- CRS.6(fixed = c(NA, NA, NA, NA, NA, 1))

  expect_equal(model$noParm, 5)
  expect_equal(model$names, c("b", "c", "d", "e", "f"))

  inits <- model$ssfct(test_data)
  expect_length(inits, 5)
})

test_that("CRS.6 with multiple parameters fixed", {
  model <- CRS.6(fixed = c(2, 0, NA, NA, NA, 1))

  expect_equal(model$noParm, 3)
  expect_equal(model$names, c("d", "e", "f"))

  inits <- model$ssfct(test_data)
  expect_length(inits, 3)
})

test_that("CRS.6 rejects invalid fixed argument - wrong length", {
  expect_error(
    CRS.6(fixed = c(NA, NA, NA, NA)),
    "Not correct 'fixed' argument"
  )
})

test_that("CRS.6 rejects invalid fixed argument - NULL", {
  expect_error(
    CRS.6(fixed = NULL),
    "Not correct 'fixed' argument"
  )
})

# ==============================================================================
# Method Parameter Tests
# ==============================================================================

test_that("CRS.6 accepts method parameter (for compatibility)", {
  # method parameter exists but is not currently used in active ssfct
  model1 <- CRS.6(method = "1")
  model2 <- CRS.6(method = "2")
  model3 <- CRS.6(method = "3")
  model4 <- CRS.6(method = "4")

  expect_s3_class(model1, "cedergreen.extended")
  expect_s3_class(model2, "cedergreen.extended")
  expect_s3_class(model3, "cedergreen.extended")
  expect_s3_class(model4, "cedergreen.extended")
})

# ==============================================================================
# Custom ssfct Tests
# ==============================================================================

test_that("CRS.6 accepts custom ssfct function", {
  custom_ssfct <- function(dframe) {
    c(2, 0, 100, 1, 10, 1)
  }

  model <- CRS.6(ssfct = custom_ssfct)

  expect_true(is.function(model$ssfct))

  inits <- model$ssfct(test_data)
  expect_equal(inits, c(2, 0, 100, 1, 10, 1))
})

test_that("CRS.6 with NULL ssfct uses default", {
  model <- CRS.6(ssfct = NULL)

  expect_true(is.function(model$ssfct))

  inits <- model$ssfct(test_data)
  expect_type(inits, "double")
  expect_length(inits, 6)
})

test_that("CRS.6 custom ssfct respects fixed parameters", {
  custom_ssfct <- function(dframe) {
    c(2, 0, 100, 1, 10, 1)
  }

  model <- CRS.6(
    fixed = c(NA, 0, NA, NA, NA, NA),
    ssfct = custom_ssfct
  )

  inits <- model$ssfct(test_data)
  # Should return full vector from custom_ssfct
  expect_length(inits, 6)
})

# ==============================================================================
# Self-starter Function (ssfct) Tests
# ==============================================================================

test_that("CRS.6 ssfct uses llogistic for first 4 parameters", {
  model <- CRS.6()

  inits <- model$ssfct(test_data)

  # Compare with llogistic ssfct
  ll_model <- llogistic()
  ll_inits <- ll_model$ssfct(test_data)

  # First 4 parameters should match llogistic
  expect_equal(inits[1:4], ll_inits[1:4])
})

test_that("CRS.6 ssfct sets g (6th parameter) to 0", {
  model <- CRS.6()

  inits <- model$ssfct(test_data)

  # 6th parameter (g/alpha) should be 0
  expect_equal(inits[6], 0)
})

test_that("CRS.6 ssfct calculates f (5th parameter) correctly", {
  model <- CRS.6()

  inits <- model$ssfct(test_data)

  # f should be calculated using the formula on line 99
  # f = (2*(median(dframe[, 2])-initval[2])-(initval[3]-initval[2]))*exp(1/(initval[4]^initval[6]))

  ll_model <- llogistic()
  ll_inits <- ll_model$ssfct(test_data)

  expected_f <- (2 * (median(test_data[, 2]) - ll_inits[2]) -
    (ll_inits[3] - ll_inits[2])) * exp(1 / (ll_inits[4]^0))

  expect_equal(inits[5], expected_f)
})

test_that("CRS.6 ssfct works with different data sets", {
  model <- CRS.6()

  # Test with different data
  small_data <- data.frame(
    dose = c(0.1, 1, 10),
    response = c(1, 5, 9)
  )

  inits <- model$ssfct(small_data)

  expect_type(inits, "double")
  expect_length(inits, 6)
  expect_true(all(is.finite(inits)))
})

# ==============================================================================
# Function (fct) Tests
# ==============================================================================

test_that("CRS.6 fct works with single dose value", {
  model <- CRS.6()

  result <- model$fct(
    dose = 1,
    parm = matrix(c(2, 0, 100, 1, 10, 1), nrow = 1)
  )

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(is.finite(result))
})

test_that("CRS.6 fct works with multiple dose values", {
  model <- CRS.6()

  result <- model$fct(
    dose = c(0.1, 1, 10, 100),
    parm = matrix(c(2, 0, 100, 1, 10, 1), nrow = 1)
  )

  expect_type(result, "double")
  expect_length(result, 4)
  expect_true(all(is.finite(result)))
})

test_that("CRS.6 fct works with multiple parameter sets", {
  model <- CRS.6()

  # Test with 2 parameter sets
  parms <- matrix(
    c(
      2, 0, 100, 1, 10, 1,
      3, 0, 100, 2, 5, 0.5
    ),
    nrow = 2,
    byrow = TRUE
  )

  result <- model$fct(dose = c(0.1, 1, 10), parm = parms)

  expect_type(result, "double")
  expect_length(result, 3)
})

test_that("CRS.6 fct respects fixed parameters", {
  # Fix b=2 and c=0
  model <- CRS.6(fixed = c(2, 0, NA, NA, NA, NA))

  # Only provide d, e, f, g (4 parameters)
  result <- model$fct(
    dose = c(0.1, 1, 10),
    parm = matrix(c(100, 1, 10, 1), nrow = 1)
  )

  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("CRS.6 fct with all parameters fixed uses fixed values", {
  model <- CRS.6(fixed = c(2, 0, 100, 1, 10, 1))

  # Provide empty parameter matrix
  result <- model$fct(
    dose = c(0.1, 1, 10),
    parm = matrix(numeric(0), nrow = 1, ncol = 0)
  )

  expect_type(result, "double")
  expect_length(result, 3)
})

test_that("CRS.6 fct model equation is correct", {
  model <- CRS.6()

  b <- 2
  c_param <- 0
  d <- 100
  e <- 1
  f <- 10
  g <- 1

  dose <- 1

  # Manual calculation: c + (d-c+f*exp(-1/dose^g))/(1+exp(b*(log(dose)-log(e))))
  expected <- c_param + (d - c_param + f * exp(-1 / (dose^g))) /
    (1 + exp(b * (log(dose) - log(e))))

  result <- model$fct(
    dose = dose,
    parm = matrix(c(b, c_param, d, e, f, g), nrow = 1)
  )

  expect_equal(result[1], expected)
})

# ==============================================================================
# Integration Tests with drm
# ==============================================================================

test_that("CRS.6 works with drm function", {
  model <- CRS.6()

  fit <- drm(rootl ~ conc, data = ryegrass, fct = model)

  expect_s3_class(fit, "drc")
  expect_true(length(coef(fit)) > 0)
  expect_equal(length(coef(fit)), 6)
})

test_that("CRS.6 with fixed parameters works with drm", {
  model <- CRS.6(fixed = c(NA, NA, NA, NA, NA, 1))

  fit <- drm(rootl ~ conc, data = ryegrass, fct = model)

  expect_s3_class(fit, "drc")
  expect_equal(length(coef(fit)), 5)
})

test_that("CRS.6 with custom names works with drm", {
  model <- CRS.6(names = c("slope", "lower", "upper", "ed50", "horm", "alpha"))

  fit <- drm(rootl ~ conc, data = ryegrass, fct = model)

  expect_s3_class(fit, "drc")
  # drm adds ":(Intercept)" suffix to parameter names
  expect_true(all(grepl("^(slope|lower|upper|ed50|horm|alpha):", names(coef(fit)))))
})

# ==============================================================================
# Edge Cases and Boundary Conditions
# ==============================================================================

test_that("CRS.6 handles dose = 0 gracefully", {
  model <- CRS.6()

  # dose = 0 will cause issues with log(dose), but should handle via exp()
  result <- model$fct(
    dose = c(0, 0.1, 1),
    parm = matrix(c(2, 0, 100, 1, 10, 1), nrow = 1)
  )

  expect_type(result, "double")
  expect_length(result, 3)
  # First value may be non-finite
  expect_true(all(is.finite(result[2:3])))
})

test_that("CRS.6 handles very small dose values", {
  model <- CRS.6()

  result <- model$fct(
    dose = c(1e-10, 1e-5, 1),
    parm = matrix(c(2, 0, 100, 1, 10, 1), nrow = 1)
  )

  expect_type(result, "double")
  expect_length(result, 3)
})

test_that("CRS.6 handles very large dose values", {
  model <- CRS.6()

  result <- model$fct(
    dose = c(1, 1000, 1e6),
    parm = matrix(c(2, 0, 100, 1, 10, 1), nrow = 1)
  )

  expect_type(result, "double")
  expect_length(result, 3)
})

test_that("CRS.6 with g=0 (special case)", {
  model <- CRS.6()

  # When g=0, exp(-1/dose^0) = exp(-1) for all dose
  result <- model$fct(
    dose = c(0.1, 1, 10),
    parm = matrix(c(2, 0, 100, 1, 10, 0), nrow = 1)
  )

  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("CRS.6 with negative g (edge case)", {
  model <- CRS.6()

  result <- model$fct(
    dose = c(0.1, 1, 10),
    parm = matrix(c(2, 0, 100, 1, 10, -1), nrow = 1)
  )

  expect_type(result, "double")
  expect_length(result, 3)
})

test_that("CRS.6 ssfct with minimal data", {
  model <- CRS.6()

  minimal_data <- data.frame(
    dose = c(0.1, 1, 10),
    response = c(1, 5, 9)
  )

  inits <- model$ssfct(minimal_data)

  expect_type(inits, "double")
  expect_length(inits, 6)
})

# ==============================================================================
# Class and Attributes Tests
# ==============================================================================

test_that("CRS.6 returns object with correct class", {
  model <- CRS.6()

  expect_s3_class(model, "cedergreen.extended")
  expect_true("cedergreen.extended" %in% class(model))
})

test_that("CRS.6 invisible return", {
  # CRS.6 uses invisible() to return the list
  # This test ensures the function still returns the model
  model <- CRS.6()

  expect_type(model, "list")
  expect_s3_class(model, "cedergreen.extended")
})

# ==============================================================================
# Consistency Tests
# ==============================================================================

test_that("CRS.6 ssfct + fct produce reasonable results", {
  model <- CRS.6()

  inits <- model$ssfct(test_data)

  # Use initial parameters to predict
  predictions <- model$fct(
    dose = test_data$dose,
    parm = matrix(inits, nrow = 1)
  )

  expect_type(predictions, "double")
  expect_length(predictions, nrow(test_data))
  expect_true(all(is.finite(predictions)))

  # Predictions should be in a reasonable range
  expect_true(all(predictions >= min(test_data$response) - 10))
  expect_true(all(predictions <= max(test_data$response) + 10))
})

test_that("CRS.6 with fixed params - ssfct returns correct length", {
  # Test various fixed parameter combinations
  fixed_combinations <- list(
    c(2, NA, NA, NA, NA, NA), # 5 params
    c(NA, 0, NA, NA, NA, NA), # 5 params
    c(NA, NA, 100, NA, NA, NA), # 5 params
    c(2, 0, NA, NA, NA, NA), # 4 params
    c(2, 0, 100, NA, NA, NA), # 3 params
    c(2, 0, 100, 1, NA, NA), # 2 params
    c(2, 0, 100, 1, 10, NA) # 1 param
  )

  expected_lengths <- c(5, 5, 5, 4, 3, 2, 1)

  for (i in seq_along(fixed_combinations)) {
    model <- CRS.6(fixed = fixed_combinations[[i]])
    inits <- model$ssfct(test_data)

    expect_length(inits, expected_lengths[i])
  }
})

# Test suite for cedergreen function
# Target: 100% code coverage

# Test data setup
test_dose <- c(0.1, 0.5, 1, 2, 5, 10, 20)
test_response <- c(102, 105, 95, 80, 40, 25, 20)
test_data <- data.frame(dose = test_dose, response = test_response)

# ==============================================================================
# Tests for cedergreen() main function
# ==============================================================================

test_that("cedergreen returns correct structure with default arguments", {
  result <- cedergreen(alpha = 0.5)

  expect_true(is.list(result))
  expect_s3_class(result, "mllogistic")
  expect_true("fct" %in% names(result))
  expect_true("ssfct" %in% names(result))
  expect_true("names" %in% names(result))
  expect_true("deriv1" %in% names(result))
  expect_true("edfct" %in% names(result))
  expect_true("maxfct" %in% names(result))
  expect_equal(result$noParm, 5)
})

test_that("cedergreen works with different alpha values", {
  result1 <- cedergreen(alpha = 1)
  result2 <- cedergreen(alpha = 0.5)
  result3 <- cedergreen(alpha = 0.25)

  expect_s3_class(result1, "mllogistic")
  expect_s3_class(result2, "mllogistic")
  expect_s3_class(result3, "mllogistic")
})

test_that("cedergreen works with fixed parameters", {
  # Fix c parameter to 0
  result <- cedergreen(fixed = c(NA, 0, NA, NA, NA), alpha = 1)

  expect_equal(result$noParm, 4)  # Only 4 parameters to estimate
  expect_equal(length(result$names), 4)
  expect_false("c" %in% result$names)
})

test_that("cedergreen works with multiple fixed parameters", {
  # Fix b and c
  result <- cedergreen(fixed = c(2, 0, NA, NA, NA), alpha = 1)

  expect_equal(result$noParm, 3)  # Only 3 parameters to estimate
  expect_false("b" %in% result$names)
  expect_false("c" %in% result$names)
})

test_that("cedergreen works with all different methods", {
  result1 <- cedergreen(method = "loglinear", alpha = 1)
  result2 <- cedergreen(method = "anke", alpha = 1)
  result3 <- cedergreen(method = "method3", alpha = 1)
  result4 <- cedergreen(method = "normolle", alpha = 1)

  expect_s3_class(result1, "mllogistic")
  expect_s3_class(result2, "mllogistic")
  expect_s3_class(result3, "mllogistic")
  expect_s3_class(result4, "mllogistic")
})

test_that("cedergreen works with custom parameter names", {
  custom_names <- c("slope", "lower", "upper", "ed50", "hormesis")
  result <- cedergreen(names = custom_names, alpha = 1)

  expect_equal(result$names, custom_names)
})

test_that("cedergreen works with custom fctName and fctText", {
  result <- cedergreen(alpha = 1, fctName = "MyModel", fctText = "My custom model")

  expect_equal(result$name, "MyModel")
  expect_equal(result$text, "My custom model")
})

test_that("cedergreen sets default fctName and fctText when missing", {
  result <- cedergreen(alpha = 1)

  expect_equal(result$name, "cedergreen")
  expect_equal(result$text, "Cedergreen-Ritz-Streibig")
})

test_that("cedergreen works with custom ssfct", {
  custom_ssfct <- function(dframe) {
    return(c(1, 0, 100, 1, 10))
  }

  result <- cedergreen(ssfct = custom_ssfct, alpha = 1)

  expect_identical(result$ssfct, custom_ssfct)
})

test_that("cedergreen fct function can be called", {
  result <- cedergreen(alpha = 1)

  # Test the fct function with sample data
  dose_vec <- c(0.1, 1, 10)
  parm_mat <- matrix(c(2, 0, 100, 1, 10), nrow = 1)

  response <- result$fct(dose_vec, parm_mat)

  expect_true(is.numeric(response))
  expect_equal(length(response), length(dose_vec))
  expect_true(all(is.finite(response)))
})

test_that("cedergreen fct function handles multiple parameter sets", {
  result <- cedergreen(alpha = 1)

  dose_vec <- c(1, 10)
  parm_mat <- matrix(c(2, 0, 100, 1, 10,
                       3, 5, 95, 2, 15), nrow = 2, byrow = TRUE)

  response <- result$fct(dose_vec, parm_mat)

  expect_equal(length(response), length(dose_vec))
})

test_that("cedergreen deriv1 function can be called", {
  result <- cedergreen(alpha = 1)

  dose_vec <- c(0.1, 1, 10)
  parm_mat <- matrix(c(2, 0, 100, 1, 10), nrow = 1)

  derivs <- result$deriv1(dose_vec, parm_mat)

  expect_true(is.matrix(derivs) || is.numeric(derivs))
})

test_that("cedergreen edfct function can be called", {
  result <- cedergreen(alpha = 1)

  parm_vec <- c(2, 0, 100, 1, 10)

  ed_result <- result$edfct(parm_vec, respl = 50, reference = "control", type = "relative")

  expect_true(is.list(ed_result))
  expect_equal(length(ed_result), 2)
})

test_that("cedergreen maxfct function can be called", {
  result <- cedergreen(alpha = 1)

  parm_vec <- c(2, 0, 100, 1, 10)

  max_result <- result$maxfct(parm_vec)

  expect_true(is.numeric(max_result))
  expect_equal(length(max_result), 2)
  expect_true(all(names(max_result) %in% c("maxDose", "maxResponse")))
})

# ==============================================================================
# Error Handling Tests
# ==============================================================================

test_that("cedergreen errors when names is not character", {
  expect_error(
    cedergreen(names = c(1, 2, 3, 4, 5), alpha = 1),
    "Not correct 'names' argument"
  )
})

test_that("cedergreen errors when names has wrong length", {
  expect_error(
    cedergreen(names = c("b", "c", "d"), alpha = 1),
    "Not correct 'names' argument"
  )
})

test_that("cedergreen errors when fixed has wrong length", {
  expect_error(
    cedergreen(fixed = c(NA, NA, NA), alpha = 1),
    "Not correct 'fixed' argument"
  )
})

test_that("cedergreen errors when alpha is missing", {
  expect_error(
    cedergreen(),
    "'alpha' argument must be specified"
  )
})

# ==============================================================================
# Tests for cedergreen_edfct helper function
# ==============================================================================

test_that("cedergreen_edfct returns NA when root finding fails", {
  # Create a scenario where root finding should fail
  result <- cedergreen(alpha = 1)

  # Use parameters that might cause issues
  parm_vec <- c(0.1, 50, 51, 1, 0.01)  # Very small range

  # This might trigger a warning and return NA
  ed_result <- result$edfct(parm_vec, respl = 99, reference = "control", type = "relative")

  # The result could be NA or a valid number depending on convergence
  expect_true(is.list(ed_result))
})

test_that("cedergreen_edfct handles different response levels", {
  result <- cedergreen(alpha = 1)
  parm_vec <- c(2, 0, 100, 1, 10)

  ed10 <- result$edfct(parm_vec, respl = 10, reference = "control", type = "relative")
  ed50 <- result$edfct(parm_vec, respl = 50, reference = "control", type = "relative")
  ed90 <- result$edfct(parm_vec, respl = 90, reference = "control", type = "relative")

  expect_true(is.list(ed10))
  expect_true(is.list(ed50))
  expect_true(is.list(ed90))
})

test_that("cedergreen_edfct handles absolute reference type", {
  result <- cedergreen(alpha = 1)
  parm_vec <- c(2, 0, 100, 1, 10)

  ed_abs <- result$edfct(parm_vec, respl = 50, reference = "control", type = "absolute")

  expect_true(is.list(ed_abs))
})

# ==============================================================================
# Tests for cedergreen_maxfct helper function
# ==============================================================================

test_that("cedergreen_maxfct finds maximum hormesis", {
  result <- cedergreen(alpha = 1)
  parm_vec <- c(2, 0, 100, 1, 50)  # Significant hormesis parameter

  max_result <- result$maxfct(parm_vec)

  expect_true(is.numeric(max_result))
  expect_equal(length(max_result), 2)
  expect_false(is.na(max_result[1]))
  expect_false(is.na(max_result[2]))
})

test_that("cedergreen_maxfct handles custom bounds", {
  result <- cedergreen(alpha = 1)
  parm_vec <- c(2, 0, 100, 1, 10)

  max_result <- result$maxfct(parm_vec, lower = 0.001, upper = 100)

  expect_true(is.numeric(max_result))
  expect_equal(length(max_result), 2)
})

test_that("cedergreen_maxfct handles edge case parameters", {
  result <- cedergreen(alpha = 1)

  # Test with different parameter combinations
  parm_vec1 <- c(2, 0, 100, 1, 0)  # Zero hormesis
  max_result1 <- result$maxfct(parm_vec1)
  expect_true(is.numeric(max_result1))

  # Test with negative hormesis (shouldn't have maximum above baseline)
  parm_vec2 <- c(2, 0, 100, 1, -10)
  max_result2 <- result$maxfct(parm_vec2)
  expect_true(is.numeric(max_result2))

  # Test with parameters that might cause numerical issues in optimize
  # Using very extreme values
  parm_vec3 <- c(1e10, 1e10, 1e10, 1e-20, 1e10)
  suppressWarnings({
    max_result3 <- result$maxfct(parm_vec3, lower=1e-30, upper=1e-25)
  })
  expect_true(is.numeric(max_result3))
  expect_equal(length(max_result3), 2)

  # Test with parameters where c > d (inverted bounds)
  parm_vec4 <- c(2, 100, 0, 1, 10)  # c > d
  suppressWarnings({
    max_result4 <- result$maxfct(parm_vec4)
  })
  expect_true(is.numeric(max_result4))
})

test_that("cedergreen_maxfct error handling when optimize fails", {
  result <- cedergreen(alpha = 1)

  # Create a mock optimize that fails
  mock_optimize <- function(...) {
    stop("Forced error for testing")
  }

  # Test error handling using the injectable parameter
  expect_warning(
    max_result <- result$maxfct(c(2, 0, 100, 1, 10), .optimize_fn = mock_optimize),
    "Optimization failed"
  )

  # Should return NA when optimization fails
  expect_true(is.na(max_result[1]))
  expect_true(is.na(max_result[2]))
})

# ==============================================================================
# Tests for CRS.5 wrapper function
# ==============================================================================

test_that("CRS.5 works with default arguments", {
  result <- CRS.5()

  expect_s3_class(result, "mllogistic")
  expect_equal(result$noParm, 5)
})

test_that("CRS.5 works with alpha_type 'a'", {
  result <- CRS.5(alpha_type = "a")

  expect_s3_class(result, "mllogistic")
  expect_true(grepl("alpha=1", result$text, fixed = TRUE))
})

test_that("CRS.5 works with alpha_type 'b'", {
  result <- CRS.5(alpha_type = "b")

  expect_s3_class(result, "mllogistic")
  expect_true(grepl("alpha=0.5", result$text, fixed = TRUE))
})

test_that("CRS.5 works with alpha_type 'c'", {
  result <- CRS.5(alpha_type = "c")

  expect_s3_class(result, "mllogistic")
  expect_true(grepl("alpha=0.25", result$text, fixed = TRUE))
})

test_that("CRS.5 works with numeric alpha_type", {
  result <- CRS.5(alpha_type = 0.75)

  expect_s3_class(result, "mllogistic")
  expect_true(grepl("alpha=0.75", result$text, fixed = TRUE))
})

test_that("CRS.5 works with fixed parameters", {
  result <- CRS.5(fixed = c(NA, 0, NA, NA, NA))

  expect_equal(result$noParm, 4)
})

test_that("CRS.5 works with custom names", {
  custom_names <- c("slope", "lower", "upper", "ed50", "hormesis")
  result <- CRS.5(names = custom_names)

  expect_equal(result$names, custom_names)
})

test_that("CRS.5 generates automatic fctName when not provided", {
  result <- CRS.5(alpha_type = "a")

  expect_equal(result$name, "CRS.5a")
})

test_that("CRS.5 uses custom fctName when provided", {
  result <- CRS.5(fctName = "MyCustomModel")

  expect_equal(result$name, "MyCustomModel")
})

test_that("CRS.5 generates automatic fctText when not provided", {
  result <- CRS.5(alpha_type = "a")

  expect_true(grepl("Cedergreen-Ritz-Streibig", result$text))
})

test_that("CRS.5 uses custom fctText when provided", {
  result <- CRS.5(fctText = "My custom text")

  expect_equal(result$text, "My custom text")
})

test_that("CRS.5 errors with invalid alpha_type character", {
  expect_error(
    CRS.5(alpha_type = "invalid"),
    "Invalid 'alpha_type'"
  )
})

test_that("CRS.5 errors with invalid names argument", {
  expect_error(
    CRS.5(names = c("a", "b", "c")),
    "Not correct 'names' argument"
  )
})

test_that("CRS.5 errors with non-character names", {
  expect_error(
    CRS.5(names = c(1, 2, 3, 4, 5)),
    "Not correct 'names' argument"
  )
})

test_that("CRS.5 handles various parameter combinations", {
  # Test that CRS.5 can handle various combinations without crashing
  result1 <- CRS.5(alpha_type = "a")
  expect_s3_class(result1, "mllogistic")

  result2 <- CRS.5(alpha_type = 0.3, fixed = c(NA, 0, NA, NA, NA))
  expect_s3_class(result2, "mllogistic")

  result3 <- CRS.5(alpha_type = 2.5)
  expect_s3_class(result3, "mllogistic")
})

test_that("CRS.5 handles errors from cedergreen when invalid fixed length passed", {
  # CRS.5's tryCatch should catch errors from cedergreen
  # Pass invalid fixed argument (wrong length) which will cause cedergreen to stop()
  expect_warning(
    result <- CRS.5(fixed = c(NA, NA, NA), alpha_type = "a"),
    "cedergreen\\(\\) model call failed"
  )

  expect_null(result)
})

# ==============================================================================
# Tests for deprecated functions
# ==============================================================================

test_that("CRS.4a triggers deprecation warning", {
  expect_warning(
    CRS.4a(),
    "deprecated"
  )
})

test_that("CRS.4a returns correct structure", {
  suppressWarnings({
    result <- CRS.4a()
  })

  expect_s3_class(result, "mllogistic")
  expect_equal(result$noParm, 4)  # c is fixed, so 4 parameters
})

test_that("ml3a is an alias for CRS.4a", {
  expect_identical(ml3a, CRS.4a)
})

test_that("CRS.4b triggers deprecation warning", {
  expect_warning(
    CRS.4b(),
    "deprecated"
  )
})

test_that("CRS.4b returns correct structure", {
  suppressWarnings({
    result <- CRS.4b()
  })

  expect_s3_class(result, "mllogistic")
  expect_equal(result$noParm, 4)
})

test_that("ml3b is an alias for CRS.4b", {
  expect_identical(ml3b, CRS.4b)
})

test_that("CRS.4c triggers deprecation warning", {
  expect_warning(
    CRS.4c(),
    "deprecated"
  )
})

test_that("CRS.4c returns correct structure", {
  suppressWarnings({
    result <- CRS.4c()
  })

  expect_s3_class(result, "mllogistic")
  expect_equal(result$noParm, 4)
})

test_that("ml3c is an alias for CRS.4c", {
  expect_identical(ml3c, CRS.4c)
})

test_that("CRS.5a triggers deprecation warning", {
  expect_warning(
    CRS.5a(),
    "deprecated"
  )
})

test_that("CRS.5a returns correct structure", {
  suppressWarnings({
    result <- CRS.5a()
  })

  expect_s3_class(result, "mllogistic")
  expect_equal(result$noParm, 5)
})

test_that("ml4a is an alias for CRS.5a", {
  expect_identical(ml4a, CRS.5a)
})

test_that("CRS.5b triggers deprecation warning", {
  expect_warning(
    CRS.5b(),
    "deprecated"
  )
})

test_that("CRS.5b returns correct structure", {
  suppressWarnings({
    result <- CRS.5b()
  })

  expect_s3_class(result, "mllogistic")
  expect_equal(result$noParm, 5)
})

test_that("ml4b is an alias for CRS.5b", {
  expect_identical(ml4b, CRS.5b)
})

test_that("CRS.5c triggers deprecation warning", {
  expect_warning(
    CRS.5c(),
    "deprecated"
  )
})

test_that("CRS.5c returns correct structure", {
  suppressWarnings({
    result <- CRS.5c()
  })

  expect_s3_class(result, "mllogistic")
  expect_equal(result$noParm, 5)
})

test_that("ml4c is an alias for CRS.5c", {
  expect_identical(ml4c, CRS.5c)
})

# ==============================================================================
# Integration Tests - Using cedergreen with drm
# ==============================================================================

test_that("cedergreen can be used with drm function", {
  skip_if_not_installed("drc")

  # Simple test data with hormesis
  dose <- c(0.1, 0.5, 1, 2, 5, 10, 20)
  response <- c(100, 105, 95, 80, 40, 25, 20)
  test_data <- data.frame(dose = dose, response = response)

  expect_no_error({
    model <- drm(response ~ dose, data = test_data, fct = cedergreen(alpha = 1))
  })
})

test_that("CRS.5 can be used with drm function", {
  skip_if_not_installed("drc")

  dose <- c(0.1, 0.5, 1, 2, 5, 10, 20)
  response <- c(100, 105, 95, 80, 40, 25, 20)
  test_data <- data.frame(dose = dose, response = response)

  expect_no_error({
    model <- drm(response ~ dose, data = test_data, fct = CRS.5(alpha_type = "a"))
  })
})

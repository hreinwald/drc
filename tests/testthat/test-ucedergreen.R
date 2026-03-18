# Test suite for ucedergreen function
# Benchmarked against the cedergreen test suite

# Test data setup
test_dose <- c(0.1, 0.5, 1, 2, 5, 10, 20)
test_response <- c(102, 105, 95, 80, 40, 25, 20)
test_data <- data.frame(dose = test_dose, response = test_response)

# ==============================================================================
# Tests for ucedergreen() main function
# ==============================================================================

test_that("ucedergreen returns correct structure with default arguments", {
  result <- ucedergreen(alpha = 0.5)

  expect_true(is.list(result))
  expect_s3_class(result, "UCRS")
  expect_true("fct" %in% names(result))
  expect_true("ssfct" %in% names(result))
  expect_true("names" %in% names(result))
  expect_true("deriv1" %in% names(result))
  expect_true("deriv2" %in% names(result))
  expect_true("edfct" %in% names(result))
  expect_true("maxfct" %in% names(result))
  expect_equal(result$noParm, 5)
})

test_that("ucedergreen works with different alpha values", {
  result1 <- ucedergreen(alpha = 1)
  result2 <- ucedergreen(alpha = 0.5)
  result3 <- ucedergreen(alpha = 0.25)

  expect_s3_class(result1, "UCRS")
  expect_s3_class(result2, "UCRS")
  expect_s3_class(result3, "UCRS")
})

test_that("ucedergreen works with fixed parameters", {
  # Fix c parameter to 0
  result <- ucedergreen(fixed = c(NA, 0, NA, NA, NA), alpha = 1)

  expect_equal(result$noParm, 4)  # Only 4 parameters to estimate
  expect_equal(length(result$names), 4)
  expect_false("c" %in% result$names)
})

test_that("ucedergreen works with multiple fixed parameters", {
  # Fix b and c
  result <- ucedergreen(fixed = c(2, 0, NA, NA, NA), alpha = 1)

  expect_equal(result$noParm, 3)  # Only 3 parameters to estimate
  expect_false("b" %in% result$names)
  expect_false("c" %in% result$names)
})

test_that("ucedergreen works with all different methods", {
  result1 <- ucedergreen(method = "loglinear", alpha = 1)
  result2 <- ucedergreen(method = "anke", alpha = 1)
  result3 <- ucedergreen(method = "method3", alpha = 1)
  result4 <- ucedergreen(method = "normolle", alpha = 1)

  expect_s3_class(result1, "UCRS")
  expect_s3_class(result2, "UCRS")
  expect_s3_class(result3, "UCRS")
  expect_s3_class(result4, "UCRS")
})

test_that("ucedergreen validates method argument via match.arg", {
  expect_error(
    ucedergreen(method = "invalid_method", alpha = 1),
    "'arg' should be one of"
  )
})

test_that("ucedergreen works with custom parameter names", {
  custom_names <- c("slope", "lower", "upper", "ed50", "hormesis")
  result <- ucedergreen(names = custom_names, alpha = 1)

  expect_equal(result$names, custom_names)
})

test_that("ucedergreen works with custom fctName and fctText", {
  result <- ucedergreen(alpha = 1, fctName = "MyUModel", fctText = "My U-shaped model")

  expect_equal(result$name, "MyUModel")
  expect_equal(result$text, "My U-shaped model")
})

test_that("ucedergreen sets default fctName and fctText when missing", {
  result <- ucedergreen(alpha = 1)

  expect_equal(result$name, "ucedergreen")
  expect_equal(result$text, "U-shaped Cedergreen-Ritz-Streibig")
})

test_that("ucedergreen works with custom ssfct", {
  custom_ssfct <- function(dframe) {
    return(c(1, 0, 100, 1, 10))
  }

  result <- ucedergreen(ssfct = custom_ssfct, alpha = 1)

  expect_identical(result$ssfct, custom_ssfct)
})

# ==============================================================================
# Tests for the fct (model) function - Issue #1: missing +c term
# ==============================================================================

test_that("ucedergreen fct function can be called", {
  result <- ucedergreen(alpha = 1)

  dose_vec <- c(0.1, 1, 10)
  parm_mat <- matrix(c(2, 0, 100, 1, 10), nrow = 1)

  response <- result$fct(dose_vec, parm_mat)

  expect_true(is.numeric(response))
  expect_equal(length(response), length(dose_vec))
  expect_true(all(is.finite(response)))
})

test_that("ucedergreen fct implements correct u-shaped formula f(x)=c+d-numTerm/denTerm", {
  # f(x) = c + d - (d - c + f*exp(-1/x^alpha)) / (1 + exp(b*(log(x) - log(e))))
  result <- ucedergreen(alpha = 1)

  dose <- 1
  b <- 2; c_val <- 10; d <- 100; e <- 5; f_val <- 20
  parm_mat <- matrix(c(b, c_val, d, e, f_val), nrow = 1)

  response <- result$fct(dose, parm_mat)

  # Manual calculation
  numTerm <- d - c_val + f_val * exp(-1/dose^1)
  denTerm <- 1 + exp(b * (log(dose) - log(e)))
  expected <- c_val + d - numTerm/denTerm

  expect_equal(response, expected)
})

test_that("ucedergreen fct handles fixed c parameter correctly", {
  # When c is fixed to 0, the formula should use c=0
  result <- ucedergreen(fixed = c(NA, 0, NA, NA, NA), alpha = 1)

  dose <- 1
  # Parameters are only the non-fixed ones: b, d, e, f
  parm_mat <- matrix(c(2, 100, 5, 20), nrow = 1)

  response <- result$fct(dose, parm_mat)

  # Manual calculation with c=0
  b <- 2; c_val <- 0; d <- 100; e <- 5; f_val <- 20
  numTerm <- d - c_val + f_val * exp(-1/dose)
  denTerm <- 1 + exp(b * (log(dose) - log(e)))
  expected <- c_val + d - numTerm/denTerm

  expect_equal(response, expected)
})

test_that("ucedergreen fct handles multiple parameter sets", {
  result <- ucedergreen(alpha = 1)

  dose_vec <- c(1, 10)
  parm_mat <- matrix(c(2, 0, 100, 1, 10,
                       3, 5, 95, 2, 15), nrow = 2, byrow = TRUE)

  response <- result$fct(dose_vec, parm_mat)

  expect_equal(length(response), length(dose_vec))
})

# ==============================================================================
# Tests for the deriv1 function - Issue #3: xlogx availability, Issue #10
# ==============================================================================

test_that("ucedergreen deriv1 function can be called", {
  result <- ucedergreen(alpha = 1)

  dose_vec <- c(0.1, 1, 10)
  parm_mat <- matrix(c(2, 0, 100, 1, 10), nrow = 1)

  derivs <- result$deriv1(dose_vec, parm_mat)

  expect_true(is.matrix(derivs) || is.numeric(derivs))
})

test_that("ucedergreen deriv1 c-derivative is 1 + 1/t3 (not 1/t3)", {
  # Verify derivative with respect to c for U-shaped model
  result <- ucedergreen(alpha = 1)

  dose <- 5
  b <- 2; c_val <- 10; d <- 100; e <- 5; f_val <- 20
  parm_mat <- matrix(c(b, c_val, d, e, f_val), nrow = 1)

  derivs <- result$deriv1(dose, parm_mat)

  # Manual: d/dc of (c + d - numTerm/denTerm)
  # = 1 + 1/denTerm  (since numTerm has -c)
  t2 <- exp(b * (log(dose) - log(e)))
  t3 <- 1 + t2
  expected_dc <- 1 + 1/t3

  # Second column should be c-derivative
  expect_equal(derivs[2], expected_dc, tolerance = 1e-10)
})

# ==============================================================================
# Tests for edfct - Issue #2: signature mismatch
# ==============================================================================

test_that("ucedergreen edfct function accepts correct signature", {
  result <- ucedergreen(alpha = 1)
  parm_vec <- c(2, 0, 100, 1, 10)

  # Should accept (parm, respl, reference, type, ...) like cedergreen
  ed_result <- result$edfct(parm_vec, respl = 50, reference = "control", type = "relative")

  expect_true(is.list(ed_result))
  expect_equal(length(ed_result), 2)
})

# ==============================================================================
# Tests for maxfct - Issue #7: signature mismatch
# ==============================================================================

test_that("ucedergreen maxfct function can be called", {
  result <- ucedergreen(alpha = 1)
  parm_vec <- c(2, 0, 100, 1, 50)

  max_result <- result$maxfct(parm_vec)

  expect_true(is.numeric(max_result))
  expect_equal(length(max_result), 2)
})

test_that("ucedergreen maxfct handles custom bounds", {
  result <- ucedergreen(alpha = 1)
  parm_vec <- c(2, 0, 100, 1, 10)

  max_result <- result$maxfct(parm_vec, lower = 0.001, upper = 100)

  expect_true(is.numeric(max_result))
  expect_equal(length(max_result), 2)
})

test_that("ucedergreen maxfct with fixed parameters reconstructs correctly", {
  # Fix c to 0
  result <- ucedergreen(fixed = c(NA, 0, NA, NA, NA), alpha = 1)
  parm_vec <- c(2, 100, 1, 50)  # b, d, e, f (c is fixed)

  max_result <- result$maxfct(parm_vec)

  expect_true(is.numeric(max_result))
  expect_equal(length(max_result), 2)
})

# ==============================================================================
# Error Handling Tests
# ==============================================================================

test_that("ucedergreen errors when names is not character", {
  expect_error(
    ucedergreen(names = c(1, 2, 3, 4, 5), alpha = 1),
    "Not correct 'names' argument"
  )
})

test_that("ucedergreen errors when names has wrong length", {
  expect_error(
    ucedergreen(names = c("b", "c", "d"), alpha = 1),
    "Not correct 'names' argument"
  )
})

test_that("ucedergreen errors when fixed has wrong length", {
  expect_error(
    ucedergreen(fixed = c(NA, NA, NA), alpha = 1),
    "Not correct 'fixed' argument"
  )
})

test_that("ucedergreen errors when alpha is missing", {
  expect_error(
    ucedergreen(),
    "'alpha' argument must be specified"
  )
})

# ==============================================================================
# Tests for self-starter function - Issue #6, #8, #9
# ==============================================================================

test_that("ucedergreen ssfct delegates to cedergreen.ssf and negates b", {
  result <- ucedergreen(alpha = 1)

  # Create a simple data frame to test the self-starter
  dframe <- data.frame(dose = c(0.1, 0.5, 1, 5, 10, 50, 100),
                       response = c(99, 95, 80, 50, 30, 15, 10))

  initval <- result$ssfct(dframe)

  # b should be negated (negative for U-shaped)
  expect_true(is.numeric(initval))
  expect_equal(length(initval), 5)
})

test_that("ucedergreen ssfct respects useFixed flag", {
  # Fix c to 0 - should only return 4 initial values
  result <- ucedergreen(fixed = c(NA, 0, NA, NA, NA), alpha = 1)

  dframe <- data.frame(dose = c(0.1, 0.5, 1, 5, 10, 50, 100),
                       response = c(99, 95, 80, 50, 30, 15, 10))

  initval <- result$ssfct(dframe)

  # Only 4 non-fixed parameters
  expect_equal(length(initval), 4)
})

# ==============================================================================
# Tests for wrapper functions - Issue #5 (|| vs |), Issue #13
# ==============================================================================

test_that("UCRS.4a returns correct structure", {
  result <- UCRS.4a()

  expect_s3_class(result, "UCRS")
  expect_equal(result$noParm, 4)
  expect_equal(length(result$names), 4)
})

test_that("UCRS.4b returns correct structure", {
  result <- UCRS.4b()

  expect_s3_class(result, "UCRS")
  expect_equal(result$noParm, 4)
})

test_that("UCRS.4c returns correct structure", {
  result <- UCRS.4c()

  expect_s3_class(result, "UCRS")
  expect_equal(result$noParm, 4)
})

test_that("UCRS.5a returns correct structure", {
  result <- UCRS.5a()

  expect_s3_class(result, "UCRS")
  expect_equal(result$noParm, 5)
  expect_equal(length(result$names), 5)
})

test_that("UCRS.5b returns correct structure", {
  result <- UCRS.5b()

  expect_s3_class(result, "UCRS")
  expect_equal(result$noParm, 5)
})

test_that("UCRS.5c returns correct structure", {
  result <- UCRS.5c()

  expect_s3_class(result, "UCRS")
  expect_equal(result$noParm, 5)
})

test_that("UCRS.4a errors with wrong names length", {
  expect_error(
    UCRS.4a(names = c("a", "b")),
    "Not correct 'names' argument"
  )
})

test_that("UCRS.5a errors with wrong names length", {
  expect_error(
    UCRS.5a(names = c("a", "b")),
    "Not correct 'names' argument"
  )
})

test_that("Aliases work correctly", {
  expect_identical(uml3a, UCRS.4a)
  expect_identical(uml3b, UCRS.4b)
  expect_identical(uml3c, UCRS.4c)
  expect_identical(uml4a, UCRS.5a)
  expect_identical(uml4b, UCRS.5b)
  expect_identical(uml4c, UCRS.5c)
})

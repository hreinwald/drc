# Test file for gaussian function

# Test basic functionality and correctness

test_that("gaussian returns correct structure with default arguments", {
  result <- gaussian()

  expect_s3_class(result, "gaussian")
  expect_type(result, "list")
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2", "derivx", "edfct",
                     "name", "text", "noParm", "lowerAs", "upperAs", "monoton", "fixed") %in% names(result)))
  expect_equal(result$noParm, 5)
  expect_equal(result$names, c("b", "c", "d", "e", "f"))
  expect_equal(result$name, "gaussian")
  expect_equal(result$text, "Gaussian")
  expect_null(result$deriv2)
  expect_true(is.na(result$monoton))
})

test_that("gaussian works with custom fctName and fctText", {
  result <- gaussian(fctName = "custom_name", fctText = "custom text")

  expect_equal(result$name, "custom_name")
  expect_equal(result$text, "custom text")
})

test_that("gaussian works with fixed parameters", {
  result <- gaussian(fixed = c(1, NA, NA, NA, 1))

  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("c", "d", "e"))
})

test_that("gaussian works with custom parameter names", {
  custom_names <- c("scale", "lower", "upper", "loc", "shape")
  result <- gaussian(names = custom_names)

  expect_equal(result$names, custom_names)
})

test_that("gaussian works with all parameters fixed", {
  result <- gaussian(fixed = c(1, 0, 1, 5, 2))

  expect_equal(result$noParm, 0)
  expect_length(result$names, 0)
})

test_that("gaussian works with only one parameter free", {
  result <- gaussian(fixed = c(1, 0, 1, 5, NA))

  expect_equal(result$noParm, 1)
  expect_equal(result$names, "f")
})

test_that("gaussian works with partial fixed parameters in different positions", {
  # Fix first and last
  result1 <- gaussian(fixed = c(1, NA, NA, NA, 1))
  expect_equal(result1$noParm, 3)

  # Fix middle parameters
  result2 <- gaussian(fixed = c(NA, 0, 1, NA, NA))
  expect_equal(result2$noParm, 3)

  # Fix alternating parameters
  result3 <- gaussian(fixed = c(NA, 0, NA, 5, NA))
  expect_equal(result3$noParm, 3)
})

test_that("gaussian method argument accepts different values", {
  result1 <- gaussian(method = "1")
  result2 <- gaussian(method = "2")
  result3 <- gaussian(method = "3")
  result4 <- gaussian(method = "4")

  expect_s3_class(result1, "gaussian")
  expect_s3_class(result2, "gaussian")
  expect_s3_class(result3, "gaussian")
  expect_s3_class(result4, "gaussian")
})

test_that("gaussian works with custom ssfct", {
  custom_ssfct <- function(dframe) {
    return(c(1, 0, 1, 5, 1))
  }

  result <- gaussian(ssfct = custom_ssfct)

  dframe <- data.frame(dose = 1:5, response = 5:1)
  init_vals <- result$ssfct(dframe)

  expect_equal(init_vals, c(1, 0, 1, 5, 1))
})

# Test fct (dose-response function)

test_that("gaussian fct evaluates correctly", {
  result <- gaussian()

  # Test with known parameters: c + (d-c) * exp(-0.5 * (sqrt(((dose-e)/b)^2))^f)
  dose <- c(0, 5, 10)
  # b=2, c=0, d=1, e=5, f=2
  parm <- matrix(c(2, 0, 1, 5, 2), nrow = 1)

  output <- result$fct(dose, parm)

  expect_type(output, "double")
  expect_length(output, length(dose))
  expect_true(all(is.finite(output)))

  # At dose=e=5, exponent is 0, so output should be d=1
  expect_equal(as.numeric(output[2]), 1, tolerance = 1e-10)
})

test_that("gaussian fct handles multiple dose values with single parameter set", {
  result <- gaussian()

  dose <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  parm <- matrix(c(2, 0, 1, 5, 2), nrow = 1)

  output <- result$fct(dose, parm)

  expect_length(output, 10)
  expect_true(all(is.finite(output)))
})

test_that("gaussian fct works with matrix parm (multiple rows)", {
  result <- gaussian()

  dose <- c(1, 5, 10)
  parm <- matrix(c(2, 0, 1, 5, 2,
                   2, 0, 1, 5, 2,
                   2, 0, 1, 5, 2),
                 nrow = 3, byrow = TRUE)

  output <- result$fct(dose, parm)

  expect_type(output, "double")
  expect_length(output, length(dose))
})

test_that("gaussian fct works with fixed parameters", {
  # Fix b=2 and f=2
  result <- gaussian(fixed = c(2, NA, NA, NA, 2))

  dose <- c(0, 5, 10)
  # Only 3 free parameters: c, d, e
  parm <- matrix(c(0, 1, 5), nrow = 1)

  output <- result$fct(dose, parm)

  expect_type(output, "double")
  expect_length(output, length(dose))
  expect_true(all(is.finite(output)))
})

test_that("gaussian fct handles edge case doses", {
  result <- gaussian()
  parm <- matrix(c(2, 0, 1, 5, 2), nrow = 1)

  # Very small dose
  output_small <- result$fct(1e-10, parm)
  expect_true(is.finite(output_small))

  # Very large dose
  output_large <- result$fct(1e10, parm)
  expect_true(is.finite(output_large))
})

# Test deriv1 (derivatives w.r.t. parameters)

test_that("gaussian deriv1 evaluates correctly", {
  result <- gaussian()

  # Avoid dose=e exactly (causes numerical singularity in gradient)
  dose <- c(1, 4, 10)
  parm <- matrix(c(2, 0, 1, 5, 2), nrow = 1)

  derivs <- result$deriv1(dose, parm)

  expect_true(is.matrix(derivs))
  expect_equal(nrow(derivs), length(dose))
  expect_equal(ncol(derivs), 5)
  expect_true(all(is.finite(derivs)))
})

test_that("gaussian deriv1 works with fixed parameters", {
  result <- gaussian(fixed = c(2, NA, NA, NA, 2))

  dose <- c(1, 5, 10)
  parm <- matrix(c(0, 1, 5), nrow = 1)

  derivs <- result$deriv1(dose, parm)

  expect_true(is.matrix(derivs))
  expect_equal(nrow(derivs), length(dose))
  expect_equal(ncol(derivs), 3)  # Only 3 free parameters
})

test_that("gaussian deriv1 handles edge case doses", {
  result <- gaussian()
  parm <- matrix(c(2, 0, 1, 5, 2), nrow = 1)

  derivs_small <- result$deriv1(1e-10, parm)
  expect_true(all(is.finite(derivs_small)))

  derivs_large <- result$deriv1(1e10, parm)
  expect_true(all(is.finite(derivs_large)))
})

# Test derivx (derivative w.r.t. dose)

test_that("gaussian derivx evaluates correctly", {
  result <- gaussian()

  # Avoid dose=e exactly (causes numerical singularity in gradient)
  dose <- c(1, 4, 10)
  parm <- matrix(c(2, 0, 1, 5, 2), nrow = 1)

  derivx_result <- result$derivx(dose, parm)

  expect_true(is.matrix(derivx_result))
  expect_equal(nrow(derivx_result), length(dose))
  expect_equal(ncol(derivx_result), 1)
  expect_true(all(is.finite(derivx_result)))
})

test_that("gaussian derivx works with fixed parameters", {
  result <- gaussian(fixed = c(2, NA, NA, NA, 2))

  dose <- c(1, 5, 10)
  parm <- matrix(c(0, 1, 5), nrow = 1)

  derivx_result <- result$derivx(dose, parm)

  expect_true(is.matrix(derivx_result))
  expect_equal(nrow(derivx_result), length(dose))
})

test_that("gaussian derivx near dose=e returns near-zero derivative", {
  result <- gaussian()

  # Near dose=e, the Gaussian is near its peak, derivative is near 0
  # Avoid exact dose=e due to numerical singularity with f>1
  dose <- c(5.001)
  parm <- matrix(c(2, 0, 1, 5, 2), nrow = 1)

  derivx_result <- result$derivx(dose, parm)
  expect_equal(as.numeric(derivx_result), 0, tolerance = 1e-3)
})

# Test edfct (effective dose function)

test_that("gaussian edfct works with relative type", {
  result <- gaussian()

  # Parameters: b=2, c=0, d=1, e=5, f=2
  parm <- c(2, 0, 1, 5, 2)

  ed_result <- result$edfct(parm, respl = 50, reference = "control", type = "relative")

  expect_type(ed_result, "list")
  expect_length(ed_result, 2)
  expect_type(ed_result[[1]], "double")
  expect_type(ed_result[[2]], "double")
  expect_length(ed_result[[2]], 5)
})

test_that("gaussian edfct works with absolute type", {
  result <- gaussian()

  parm <- c(2, 0, 1, 5, 2)

  ed_result <- result$edfct(parm, respl = 0.5, reference = "control", type = "absolute")

  expect_type(ed_result, "list")
  expect_length(ed_result, 2)
})

test_that("gaussian edfct works with relative type and negative b and control reference", {
  result <- gaussian()

  # Parameters with negative b to trigger the control reference path
  parm <- c(-2, 0, 1, 5, 2)

  ed_result <- result$edfct(parm, respl = 50, reference = "control", type = "relative")

  expect_type(ed_result, "list")
  expect_length(ed_result, 2)
})

test_that("gaussian edfct works with relative type and positive b", {
  result <- gaussian()

  parm <- c(2, 0, 1, 5, 2)

  ed_result <- result$edfct(parm, respl = 50, reference = "upper", type = "relative")

  expect_type(ed_result, "list")
  expect_length(ed_result, 2)
})

test_that("gaussian edfct works with fixed parameters", {
  result <- gaussian(fixed = c(2, NA, NA, NA, 2))

  parm <- c(0, 1, 5)

  ed_result <- result$edfct(parm, respl = 50, reference = "control", type = "relative")

  expect_type(ed_result, "list")
  expect_length(ed_result, 2)
  expect_length(ed_result[[2]], 3)  # Only 3 free parameters
})

# Test ssfct (self-starter)

test_that("gaussian ssfct works with valid data", {
  result <- gaussian()

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- result$ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
  expect_true(all(is.finite(init_vals)))
})

# Test lowerAs and upperAs

test_that("gaussian lowerAs and upperAs return correct values", {
  result <- gaussian()

  parm <- c(2, 0.5, 1.5, 5, 2)
  lower_val <- result$lowerAs(parm)
  upper_val <- result$upperAs(parm)

  expect_equal(lower_val, 0.5)
  expect_equal(upper_val, 1.5)
})

test_that("gaussian lowerAs and upperAs with fixed parameters", {
  # Fix c=0.2 and d=0.9
  result <- gaussian(fixed = c(NA, 0.2, 0.9, NA, NA))

  # Only free: b, e, f
  parm <- c(2, 5, 2)
  lower_val <- result$lowerAs(parm)
  upper_val <- result$upperAs(parm)

  expect_equal(lower_val, 0.2)
  expect_equal(upper_val, 0.9)
})

# Test error handling

test_that("gaussian errors with incorrect names argument - not character", {
  expect_error(
    gaussian(names = c(1, 2, 3, 4, 5)),
    "Not correct 'names' argument"
  )
})

test_that("gaussian errors with incorrect names argument - wrong length", {
  expect_error(
    gaussian(names = c("b", "c", "d")),
    "Not correct 'names' argument"
  )
})

test_that("gaussian errors with incorrect fixed argument - wrong length", {
  expect_error(
    gaussian(fixed = c(NA, NA, NA)),
    "Not correct 'fixed' argument"
  )
})

# Test fixed field in return value

test_that("gaussian returns fixed argument in result", {
  fixed_vec <- c(1, NA, NA, NA, 2)
  result <- gaussian(fixed = fixed_vec)

  expect_equal(result$fixed, fixed_vec)
})

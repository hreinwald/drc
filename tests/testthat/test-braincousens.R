# Test file for braincousens function

# Test basic functionality and correctness

test_that("braincousens returns correct structure with default arguments", {
  result <- braincousens()

  expect_s3_class(result, "braincousens")
  expect_type(result, "list")
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2", "edfct", "maxfct", "name", "text", "noParm") %in% names(result)))
  expect_equal(result$noParm, 5)
  expect_equal(result$names, c("b", "c", "d", "e", "f"))
  expect_equal(result$name, "braincousens")
  expect_equal(result$text, "Brain-Cousens (hormesis)")
  expect_null(result$deriv2)
})

test_that("braincousens works with custom fctName and fctText", {
  result <- braincousens(fctName = "custom_name", fctText = "custom text")

  expect_equal(result$name, "custom_name")
  expect_equal(result$text, "custom text")
})

test_that("braincousens works with fixed parameters", {
  result <- braincousens(fixed = c(1, NA, NA, NA, 0))

  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("c", "d", "e"))
})

test_that("braincousens works with custom parameter names", {
  custom_names <- c("slope", "lower", "upper", "ed50", "horm")
  result <- braincousens(names = custom_names)

  expect_equal(result$names, custom_names)
})

test_that("braincousens fct evaluates correctly", {
  result <- braincousens()

  # Test with simple parameters
  dose <- c(0.1, 1, 10)
  parm <- matrix(c(1, 0, 1, 1, 0.1), nrow = 1)

  output <- result$fct(dose, parm)

  expect_type(output, "double")
  expect_length(output, length(dose))
  expect_true(all(is.finite(output)))
})

test_that("braincousens ssfct works with valid data", {
  result <- braincousens()

  # Create test data frame
  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- result$ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
  expect_true(all(is.finite(init_vals)))
})

test_that("braincousens ssfct works with custom ssfct", {
  custom_ssfct <- function(dframe) {
    return(c(1, 0, 1, 1, 0))
  }

  result <- braincousens(ssfct = custom_ssfct)

  dframe <- data.frame(dose = 1:5, response = 5:1)
  init_vals <- result$ssfct(dframe)

  expect_equal(init_vals, c(1, 0, 1, 1, 0))
})

test_that("braincousens deriv1 evaluates correctly", {
  result <- braincousens()

  dose <- c(0.1, 1, 10)
  parm <- matrix(c(1, 0, 1, 1, 0.1), nrow = 1)

  derivs <- result$deriv1(dose, parm)

  expect_true(is.matrix(derivs))
  expect_equal(nrow(derivs), length(dose))
  expect_equal(ncol(derivs), 5)
  expect_true(all(is.finite(derivs)))
})

test_that("braincousens edfct works correctly", {
  result <- braincousens()

  # Use parameters that work well with the Brain-Cousens model
  parm <- c(2, 0, 1, 1, 0.5)

  ed_result <- result$edfct(parm, respl = 50, reference = "control", type = "relative")

  expect_type(ed_result, "list")
  expect_length(ed_result, 2)
  expect_type(ed_result[[1]], "double")
  expect_type(ed_result[[2]], "double")
  expect_length(ed_result[[2]], 5)
})

test_that("braincousens maxfct works correctly", {
  result <- braincousens()

  parm <- c(2, 0, 1, 1, 0.5)

  max_result <- result$maxfct(parm)

  expect_type(max_result, "double")
  expect_length(max_result, 2)
  expect_true(all(is.finite(max_result)))
})

# Test error handling

test_that("braincousens errors with incorrect names argument - not character", {
  expect_error(
    braincousens(names = c(1, 2, 3, 4, 5)),
    "Not correct 'names' argument"
  )
})

test_that("braincousens errors with incorrect names argument - wrong length", {
  expect_error(
    braincousens(names = c("b", "c", "d")),
    "Not correct 'names' argument"
  )
})

test_that("braincousens errors with incorrect fixed argument - wrong length", {
  expect_error(
    braincousens(fixed = c(NA, NA, NA)),
    "Not correct 'fixed' argument"
  )
})

test_that("braincousens maxfct errors when b < 1", {
  result <- braincousens()

  parm <- c(0.5, 0, 1, 1, 0.5)

  expect_error(
    result$maxfct(parm),
    "Brain-Cousens model with b<1 not meaningful"
  )
})

test_that("braincousens maxfct errors when f < 0", {
  result <- braincousens()

  parm <- c(2, 0, 1, 1, -0.5)

  expect_error(
    result$maxfct(parm),
    "Brain-Cousens model with f<0 not meaningful"
  )
})

# Test edge cases

test_that("braincousens works with all parameters fixed", {
  result <- braincousens(fixed = c(1, 0, 1, 1, 0))

  expect_equal(result$noParm, 0)
  expect_length(result$names, 0)
})

test_that("braincousens works with only one parameter free", {
  result <- braincousens(fixed = c(1, 0, 1, 1, NA))

  expect_equal(result$noParm, 1)
  expect_equal(result$names, "f")
})

test_that("braincousens method argument accepts different values", {
  result1 <- braincousens(method = "1")
  result2 <- braincousens(method = "2")
  result3 <- braincousens(method = "3")
  result4 <- braincousens(method = "4")

  expect_s3_class(result1, "braincousens")
  expect_s3_class(result2, "braincousens")
  expect_s3_class(result3, "braincousens")
  expect_s3_class(result4, "braincousens")
})

test_that("braincousens fct handles edge case doses", {
  result <- braincousens()
  parm <- matrix(c(1, 0, 1, 1, 0.1), nrow = 1)

  # Very small dose
  output_small <- result$fct(1e-10, parm)
  expect_true(is.finite(output_small))

  # Very large dose
  output_large <- result$fct(1e10, parm)
  expect_true(is.finite(output_large))
})

test_that("braincousens deriv1 handles edge case doses", {
  result <- braincousens()
  parm <- matrix(c(1, 0, 1, 1, 0.1), nrow = 1)

  # Very small dose
  derivs_small <- result$deriv1(1e-10, parm)
  expect_true(all(is.finite(derivs_small)))

  # Very large dose
  derivs_large <- result$deriv1(1e10, parm)
  expect_true(all(is.finite(derivs_large)))
})

test_that("braincousens edfct works with different bounds", {
  result <- braincousens()
  parm <- c(2, 0, 1, 1, 0.5)

  ed_result <- result$edfct(parm, respl = 50, reference = "control", type = "relative",
                            lower = 1e-6, upper = 10000)

  expect_type(ed_result, "list")
  expect_length(ed_result, 2)
})

test_that("braincousens maxfct works with different bounds", {
  result <- braincousens()
  parm <- c(2, 0, 1, 1, 0.5)

  max_result <- result$maxfct(parm, lower = 1e-6, upper = 10000)

  expect_type(max_result, "double")
  expect_length(max_result, 2)
})

test_that("braincousens fct works with matrix parm", {
  result <- braincousens()

  dose <- c(0.1, 1, 10)
  # The function expects the number of rows in parm to match the number of doses
  parm <- matrix(c(1, 0, 1, 1, 0.1,
                   1, 0, 1, 1, 0.1,
                   1, 0, 1, 1, 0.1),
                 nrow = 3, byrow = TRUE)

  output <- result$fct(dose, parm)

  expect_type(output, "double")
  expect_length(output, length(dose))
})

test_that("braincousens works with partial fixed parameters in different positions", {
  # Fix first and last
  result1 <- braincousens(fixed = c(1, NA, NA, NA, 0))
  expect_equal(result1$noParm, 3)

  # Fix middle parameters
  result2 <- braincousens(fixed = c(NA, 0, 1, NA, NA))
  expect_equal(result2$noParm, 3)

  # Fix alternating parameters
  result3 <- braincousens(fixed = c(NA, 0, NA, 1, NA))
  expect_equal(result3$noParm, 3)
})

# Test BC.4 wrapper function

test_that("BC.4 returns correct structure", {
  result <- BC.4()

  expect_s3_class(result, "braincousens")
  expect_equal(result$noParm, 4)
  expect_equal(result$names, c("b", "d", "e", "f"))
  expect_equal(result$name, "BC.4")
  expect_equal(result$text, "Brain-Cousens (hormesis) with lower limit fixed at 0")
})

test_that("BC.4 works with custom names", {
  custom_names <- c("slope", "upper", "ed50", "horm")
  result <- BC.4(names = custom_names)

  expect_equal(result$names, custom_names)
})

test_that("BC.4 works with fixed parameters", {
  result <- BC.4(fixed = c(1, NA, NA, 0))

  expect_equal(result$noParm, 2)
  expect_equal(result$names, c("d", "e"))
})

test_that("BC.4 errors with incorrect names length", {
  expect_error(
    BC.4(names = c("b", "d", "e")),
    "Not correct 'names' argument"
  )
})

test_that("BC.4 errors with non-character names", {
  expect_error(
    BC.4(names = c(1, 2, 3, 4)),
    "Not correct 'names' argument"
  )
})

# Test BC.5 wrapper function

test_that("BC.5 returns correct structure", {
  result <- BC.5()

  expect_s3_class(result, "braincousens")
  expect_equal(result$noParm, 5)
  expect_equal(result$names, c("b", "c", "d", "e", "f"))
  expect_equal(result$name, "BC.5")
})

test_that("BC.5 works with custom names", {
  custom_names <- c("slope", "lower", "upper", "ed50", "horm")
  result <- BC.5(names = custom_names)

  expect_equal(result$names, custom_names)
})

test_that("BC.5 works with fixed parameters", {
  result <- BC.5(fixed = c(1, 0, NA, NA, 0))

  expect_equal(result$noParm, 2)
  expect_equal(result$names, c("d", "e"))
})

test_that("BC.5 errors with incorrect names length", {
  expect_error(
    BC.5(names = c("b", "c", "d")),
    "Not correct 'names' argument"
  )
})

test_that("BC.5 errors with non-character names", {
  expect_error(
    BC.5(names = c(1, 2, 3, 4, 5)),
    "Not correct 'names' argument"
  )
})

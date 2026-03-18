# Test file for threephase function

# ---- Correctness: Default arguments ----

test_that("threephase returns correct structure with default arguments", {
  result <- threephase()

  expect_s3_class(result, "three-phase")
  expect_type(result, "list")
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2",
                     "derivx", "edfct", "name", "text", "noParm") %in% names(result)))
  expect_equal(result$noParm, 10)
  expect_equal(result$names, c("b1", "c1", "d1", "e1", "b2", "d2", "e2", "b3", "d3", "e3"))
  expect_equal(result$name, "threephase")
  expect_equal(result$text, "Three-phase")
  expect_null(result$deriv1)
  expect_null(result$deriv2)
  expect_null(result$derivx)
  expect_null(result$edfct)
})

# ---- Correctness: Custom fctName and fctText ----

test_that("threephase works with custom fctName and fctText", {
  result <- threephase(fctName = "custom_name", fctText = "custom text")

  expect_equal(result$name, "custom_name")
  expect_equal(result$text, "custom text")
})

# ---- Correctness: Fixed parameters ----

test_that("threephase works with some fixed parameters", {
  result <- threephase(fixed = c(1, NA, NA, NA, NA, NA, NA, NA, NA, NA))

  expect_equal(result$noParm, 9)
  expect_equal(result$names, c("c1", "d1", "e1", "b2", "d2", "e2", "b3", "d3", "e3"))
})

test_that("threephase works with multiple fixed parameters", {
  result <- threephase(fixed = c(1, 0, NA, NA, 2, NA, NA, 3, NA, NA))

  expect_equal(result$noParm, 6)
  expect_equal(result$names, c("d1", "e1", "d2", "e2", "d3", "e3"))
})

# ---- Correctness: Custom parameter names ----

test_that("threephase works with custom parameter names", {
  custom_names <- c("a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10")
  result <- threephase(names = custom_names)

  expect_equal(result$names, custom_names)
})

# ---- Correctness: fct evaluates correctly ----

test_that("threephase fct evaluates correctly with single row parm", {
  result <- threephase()

  dose <- c(0.1, 1, 10)
  # 10 parameters: b1, c1, d1, e1, b2, d2, e2, b3, d3, e3
  parm <- matrix(c(1, 0, 1, 1, 1, 1, 1, 1, 1, 1), nrow = 1)

  output <- result$fct(dose, parm)

  expect_type(output, "double")
  expect_length(output, length(dose))
  expect_true(all(is.finite(output)))
})

test_that("threephase fct evaluates correctly with multiple rows", {
  result <- threephase()

  dose <- c(0.1, 1, 10)
  parm <- matrix(rep(c(1, 0, 1, 1, 1, 1, 1, 1, 1, 1), 3), nrow = 3, byrow = TRUE)

  output <- result$fct(dose, parm)

  expect_type(output, "double")
  expect_length(output, 3)
  expect_true(all(is.finite(output)))
})

test_that("threephase fct works with fixed parameters", {
  # Fix b1=1 and b2=2
  result <- threephase(fixed = c(1, NA, NA, NA, 2, NA, NA, NA, NA, NA))

  dose <- c(0.1, 1, 10)
  # parm should have 8 values (the 8 free parameters)
  parm <- matrix(c(0, 1, 1, 1, 1, 1, 1, 1), nrow = 1)

  output <- result$fct(dose, parm)

  expect_type(output, "double")
  expect_length(output, length(dose))
  expect_true(all(is.finite(output)))
})

# ---- Correctness: ssfct self-starter function ----

test_that("threephase ssfct works with valid data", {
  result <- threephase()

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- result$ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 10)
  expect_true(all(is.finite(init_vals)))
})

test_that("threephase ssfct works with fixed parameters", {
  result <- threephase(fixed = c(1, NA, NA, NA, NA, NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- result$ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 9)
  expect_true(all(is.finite(init_vals)))
})

# ---- Error Handling: Invalid names argument ----

test_that("threephase errors with non-character names", {
  expect_error(threephase(names = 1:10), "Not correct 'names' argument")
})

test_that("threephase errors with wrong length names", {
  expect_error(threephase(names = c("a", "b")), "Not correct 'names' argument")
})

# ---- Error Handling: Invalid fixed argument ----

test_that("threephase errors with wrong length fixed", {
  expect_error(threephase(fixed = c(NA, NA)), "Not correct 'fixed' argument")
})

test_that("threephase errors with too long fixed", {
  expect_error(threephase(fixed = rep(NA, 11)), "Not correct 'fixed' argument")
})

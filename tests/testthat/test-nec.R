# ============================================================================
# Tests for NEC (No Effect Concentration) dose-response model functions
# Functions: NEC, NEC.2, NEC.3, NEC.4 from R/nec.R
# ============================================================================

# --- NEC.4 (Full four-parameter model) ---

test_that("NEC.4() returns object of class 'NEC'", {
  result <- NEC.4()
  expect_s3_class(result, "NEC")
})

test_that("NEC.4() returns correct list structure", {
  result <- NEC.4()
  expected_names <- c("fct", "ssfct", "names", "deriv1", "deriv2",
                      "derivx", "edfct", "name", "text", "noParm")
  expect_named(result, expected_names)
  expect_true(is.function(result$fct))
  expect_true(is.function(result$ssfct))
  expect_equal(result$names, c("b", "c", "d", "e"))
  expect_null(result$deriv1)
  expect_null(result$deriv2)
  expect_null(result$derivx)
  expect_null(result$edfct)
  expect_equal(result$name, "NEC.4")
  expect_equal(result$text, "NEC")
  expect_equal(result$noParm, 4)
})

test_that("NEC.4() fct computes correct response below threshold", {
  result <- NEC.4()
  # b=1, c=0, d=100, e=5 -> dose <= 5 should give d=100
  parm <- matrix(c(1, 0, 100, 5), nrow = 1)
  dose <- c(0, 2, 5)
  response <- result$fct(dose, parm)
  expect_equal(response, c(100, 100, 100))
})

test_that("NEC.4() fct computes correct response above threshold", {
  result <- NEC.4()
  # b=1, c=0, d=100, e=5 -> dose > 5: c + (d-c)*exp(-b*(dose-e))
  parm <- matrix(c(1, 0, 100, 5), nrow = 1)
  dose <- c(7, 10)
  response <- result$fct(dose, parm)
  expected <- 0 + (100 - 0) * exp(-1 * (dose - 5))
  expect_equal(response, expected)
})

test_that("NEC.4() fct handles multiple rows in parm", {
  result <- NEC.4()
  # Two different parameter sets
  parm <- matrix(c(1, 2, 0, 10, 100, 50, 5, 3), nrow = 2)
  dose <- c(7, 1)  # first above threshold, second below
  response <- result$fct(dose, parm)
  # Row 1: b=1, c=0, d=100, e=5 -> dose=7 > 5: 0 + 100*exp(-1*2) = 100*exp(-2)
  # Row 2: b=2, c=10, d=50, e=3 -> dose=1 < 3: 10 + (50-10)*exp(0) = 50
  expected <- c(100 * exp(-2), 50)
  expect_equal(response, expected)
})

test_that("NEC.4() with fixed parameters works correctly", {
  # Fix c=0
  result <- NEC.4(fixed = c(NA, 0, NA, NA))
  expect_equal(result$names, c("b", "d", "e"))
  expect_equal(result$noParm, 3)

  # Test fct with 3 free parameters (b, d, e)
  parm <- matrix(c(1, 100, 5), nrow = 1)
  dose <- c(0, 7)
  response <- result$fct(dose, parm)
  # dose=0 <= 5: c + (d-c)*1 = 0 + 100 = 100
  # dose=7 > 5: 0 + 100*exp(-1*2) = 100*exp(-2)
  expected <- c(100, 100 * exp(-2))
  expect_equal(response, expected)
})

test_that("NEC.4() with all parameters fixed works", {
  result <- NEC.4(fixed = c(1, 0, 100, 5))
  expect_equal(result$names, character(0))
  expect_equal(result$noParm, 0)
})

test_that("NEC.4() with custom names works", {
  result <- NEC.4(names = c("slope", "lower", "upper", "threshold"))
  expect_equal(result$names, c("slope", "lower", "upper", "threshold"))
})

test_that("NEC.4() errors on invalid names", {
  expect_error(NEC.4(names = c("a", "b")), "Not correct names argument")
  expect_error(NEC.4(names = 1:4), "Not correct names argument")
})

test_that("NEC.4() errors on invalid fixed", {
  expect_error(NEC.4(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
  expect_error(NEC.4(fixed = c(NA, NA, NA, NA, NA)), "Not correct length of 'fixed' argument")
})


# --- NEC.3 (Three-parameter, lower limit fixed at 0) ---

test_that("NEC.3() returns object of class 'NEC'", {
  result <- NEC.3()
  expect_s3_class(result, "NEC")
})

test_that("NEC.3() returns correct structure", {
  result <- NEC.3()
  expect_equal(result$names, c("b", "d", "e"))
  expect_equal(result$noParm, 3)
  expect_equal(result$name, "NEC.3")
  expect_match(result$text, "NEC")
  expect_match(result$text, "lower limit at 0")
})

test_that("NEC.3() fct computes correctly with lower limit fixed at 0", {
  result <- NEC.3()
  # b=1, d=100, e=5 (c is fixed at 0)
  parm <- matrix(c(1, 100, 5), nrow = 1)
  dose <- c(0, 5, 7, 10)
  response <- result$fct(dose, parm)
  # dose <= 5: 0 + (100-0)*1 = 100
  # dose=7: 0 + 100*exp(-1*2) = 100*exp(-2)
  # dose=10: 0 + 100*exp(-1*5) = 100*exp(-5)
  expected <- c(100, 100, 100 * exp(-2), 100 * exp(-5))
  expect_equal(response, expected)
})

test_that("NEC.3() with fixed parameters works", {
  result <- NEC.3(fixed = c(1, NA, NA))
  expect_equal(result$names, c("d", "e"))
  expect_equal(result$noParm, 2)
})

test_that("NEC.3() errors on invalid names", {
  expect_error(NEC.3(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(NEC.3(names = 1:3), "Not correct 'names' argument")
})

test_that("NEC.3() errors on invalid fixed length", {
  expect_error(NEC.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
  expect_error(NEC.3(fixed = c(NA, NA, NA, NA)), "Not correct length of 'fixed' argument")
})


# --- NEC.2 (Two-parameter, lower=0, upper fixed) ---

test_that("NEC.2() returns object of class 'NEC'", {
  result <- NEC.2()
  expect_s3_class(result, "NEC")
})

test_that("NEC.2() returns correct structure", {
  result <- NEC.2()
  expect_equal(result$names, c("b", "e"))
  expect_equal(result$noParm, 2)
  expect_equal(result$name, "NEC.2")
  expect_match(result$text, "NEC")
  expect_match(result$text, "lower limit at 0")
  expect_match(result$text, "upper limit at 1")
})

test_that("NEC.2() fct computes correctly with defaults (upper=1)", {
  result <- NEC.2()
  # b=1, e=0.5 (c=0, d=1)
  parm <- matrix(c(1, 0.5), nrow = 1)
  dose <- c(0, 0.3, 0.5, 0.8, 1)
  response <- result$fct(dose, parm)
  # dose <= 0.5: 0 + (1-0)*1 = 1
  # dose=0.8: 0 + 1*exp(-1*0.3) = exp(-0.3)
  # dose=1: 0 + 1*exp(-1*0.5) = exp(-0.5)
  expected <- c(1, 1, 1, exp(-0.3), exp(-0.5))
  expect_equal(response, expected)
})

test_that("NEC.2() with custom upper limit works", {
  result <- NEC.2(upper = 50)
  parm <- matrix(c(1, 5), nrow = 1)
  dose <- c(0, 7)
  response <- result$fct(dose, parm)
  # dose=0 <= 5: 0 + (50-0)*1 = 50
  # dose=7 > 5: 0 + 50*exp(-1*2) = 50*exp(-2)
  expected <- c(50, 50 * exp(-2))
  expect_equal(response, expected)
  expect_match(result$text, "upper limit at 50")
})

test_that("NEC.2() with fixed parameters works", {
  result <- NEC.2(fixed = c(1, NA))
  expect_equal(result$names, c("e"))
  expect_equal(result$noParm, 1)
})

test_that("NEC.2() errors on invalid names", {
  expect_error(NEC.2(names = c("a")), "Not correct 'names' argument")
  expect_error(NEC.2(names = 1:2), "Not correct 'names' argument")
})

test_that("NEC.2() errors on invalid fixed length", {
  expect_error(NEC.2(fixed = c(NA)), "Not correct length of 'fixed' argument")
  expect_error(NEC.2(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
})


# --- Base NEC function (called through NEC.4/NEC.3/NEC.2) ---

test_that("NEC base function errors on invalid names", {
  expect_error(NEC.4(names = c("a", "b", "c")), "Not correct names argument")
})

test_that("NEC base function errors on invalid fixed length", {
  expect_error(NEC.4(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("NEC base function uses fctName when provided (via NEC.4)", {
  result <- NEC.4()
  expect_equal(result$name, "NEC.4")
})

test_that("NEC base function uses fctText when provided (via NEC.3)", {
  result <- NEC.3()
  expect_match(result$text, "NEC with lower limit at 0")
})

test_that("NEC base function errors on invalid names directly", {
  necFn <- drc:::NEC
  expect_error(necFn(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(necFn(names = 1:4), "Not correct 'names' argument")
})

test_that("NEC base function errors on invalid fixed directly", {
  necFn <- drc:::NEC
  expect_error(necFn(fixed = c(NA, NA)), "Not correct 'fixed' argument")
})

test_that("NEC base function defaults name/text when missing", {
  # Call NEC directly (not through convenience functions)
  # NEC is not exported, so assign to a local variable to get a clean match.call()
  necFn <- drc:::NEC
  result <- necFn()
  expect_equal(result$name, "necFn")
  expect_equal(result$text, "NEC")
})


# --- ssfct (self-starter function) tests ---

test_that("NEC.4() ssfct returns initial parameter estimates", {
  result <- NEC.4()

  # Create a simple data frame mimicking dose-response data
  # Need columns: dose, response (at minimum)
  dose <- c(0, 0.1, 0.5, 1, 2, 5, 10, 20)
  response <- c(100, 100, 98, 90, 60, 20, 5, 1)
  dframe <- data.frame(dose = dose, response = response)

  initVals <- result$ssfct(dframe)
  expect_type(initVals, "double")
  expect_length(initVals, 4)  # 4 free parameters
  expect_true(all(is.finite(initVals)))
})

test_that("NEC.3() ssfct returns 3 initial values", {
  result <- NEC.3()
  dose <- c(0, 0.1, 0.5, 1, 2, 5, 10, 20)
  response <- c(100, 100, 98, 90, 60, 20, 5, 1)
  dframe <- data.frame(dose = dose, response = response)

  initVals <- result$ssfct(dframe)
  expect_type(initVals, "double")
  expect_length(initVals, 3)
})

test_that("NEC.2() ssfct returns 2 initial values", {
  result <- NEC.2(upper = 100)
  dose <- c(0, 0.1, 0.5, 1, 2, 5, 10, 20)
  response <- c(100, 100, 98, 90, 60, 20, 5, 1)
  dframe <- data.frame(dose = dose, response = response)

  initVals <- result$ssfct(dframe)
  expect_type(initVals, "double")
  expect_length(initVals, 2)
})


# --- Integration test with drm ---

test_that("NEC.4() works with drm() on ryegrass data", {
  data(ryegrass, package = "drc")
  model <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4())
  expect_s3_class(model, "drc")

  coefs <- coef(model)
  expect_length(coefs, 4)
  expect_true(all(is.finite(coefs)))
})

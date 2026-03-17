# tests/testthat/test-gompertz.R
# Comprehensive tests for R/gompertz.R: gompertz(), G.2(), G.3(), G.3u(), G.4()
# and nested functions: fct, deriv1, derivx, edfct

# ========================================================================
# Test: gompertz() argument validation
# ========================================================================

test_that("gompertz() errors on invalid 'names' argument", {
  expect_error(gompertz(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(gompertz(names = 123), "Not correct 'names' argument")
  expect_error(gompertz(names = c("a", "b", "c")), "Not correct 'names' argument")
})

test_that("gompertz() errors on invalid 'fixed' argument", {
  expect_error(gompertz(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(gompertz(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
  expect_error(gompertz(fixed = c(NA, NA, NA, NA, NA)), "Not correct 'fixed' argument")
})

# ========================================================================
# Test: gompertz() return structure
# ========================================================================

test_that("gompertz() returns object of class 'gompertz'", {
  result <- gompertz()
  expect_s3_class(result, "gompertz")
})

test_that("gompertz() return list has correct structure", {
  result <- gompertz()
  expect_true(is.function(result$fct))
  expect_true(is.function(result$ssfct))
  expect_true(is.function(result$deriv1))
  expect_null(result$deriv2)
  expect_true(is.function(result$derivx))
  expect_true(is.function(result$edfct))
  expect_equal(result$noParm, 4)
  expect_equal(result$names, c("b", "c", "d", "e"))
})

test_that("gompertz() default name and text are correct", {
  result <- gompertz()
  expect_equal(result$name, "gompertz")
  expect_equal(result$text, "Gompertz")
})

test_that("gompertz() custom fctName and fctText override defaults", {
  result <- gompertz(fctName = "myModel", fctText = "my description")
  expect_equal(result$name, "myModel")
  expect_equal(result$text, "my description")
})

test_that("gompertz() respects custom parameter names", {
  result <- gompertz(names = c("slope", "lower", "upper", "mid"))
  expect_equal(result$names, c("slope", "lower", "upper", "mid"))
})

test_that("gompertz() handles fixed parameters correctly", {
  result <- gompertz(fixed = c(NA, 0, NA, NA))
  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("b", "d", "e"))
})

test_that("gompertz() uses ssfct when provided", {
  custom_ssfct <- function(dframe) { c(1, 0, 100, 5) }
  result <- gompertz(ssfct = custom_ssfct)
  expect_identical(result$ssfct, custom_ssfct)
})

test_that("gompertz() uses default ssfct when not provided", {
  result <- gompertz()
  expect_true(is.function(result$ssfct))
})

test_that("gompertz() lowerAs, upperAs, monoton are functions", {
  result <- gompertz()
  expect_true(is.function(result$lowerAs))
  expect_true(is.function(result$upperAs))
  expect_true(is.function(result$monoton))
})

# ========================================================================
# Test: fct (the nonlinear function)
# Gompertz: f(x) = c + (d-c)*exp(-exp(b*(x-e)))
# ========================================================================

test_that("gompertz fct produces expected values with all parameters free", {
  mod <- gompertz()
  # Parameters: b, c, d, e
  # b=1, c=0, d=100, e=5
  parm <- matrix(c(1, 0, 100, 5), nrow = 1)

  # At dose = e = 5: f(5) = 0 + (100-0)*exp(-exp(1*(5-5))) = 100*exp(-1) = 36.7879...
  result <- mod$fct(5, parm)
  expected <- 100 * exp(-exp(0))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)

  # At dose = 0: f(0) = 100*exp(-exp(1*(0-5))) = 100*exp(-exp(-5))
  result0 <- mod$fct(0, parm)
  expected0 <- 100 * exp(-exp(-5))
  expect_equal(as.numeric(result0), expected0, tolerance = 1e-10)
})

test_that("gompertz fct handles fixed parameters", {
  # Fix c=0 (like G.3 does)
  mod <- gompertz(fixed = c(NA, 0, NA, NA))
  parm <- matrix(c(1, 100, 5), nrow = 1)  # b, d, e

  result <- mod$fct(5, parm)
  expected <- 100 * exp(-exp(0))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("gompertz fct handles multiple doses", {
  mod <- gompertz()
  parm <- matrix(rep(c(1, 0, 100, 5), each = 3), nrow = 3)
  doses <- c(0, 5, 100)

  result <- mod$fct(doses, parm)
  expect_length(result, 3)
  # At dose = 5 (inflection point): 100*exp(-1)
  expect_equal(as.numeric(result[2]), 100 * exp(-1), tolerance = 1e-10)
})

test_that("gompertz fct returns gradient attribute", {
  mod <- gompertz()
  parm <- matrix(c(1, 0, 100, 5), nrow = 1)

  result <- mod$fct(5, parm)
  expect_true(!is.null(attr(result, "gradient")))
  grad <- attr(result, "gradient")
  expect_equal(ncol(grad), 4)
})

test_that("gompertz fct decreasing curve when b > 0", {
  mod <- gompertz()
  parm <- matrix(c(1, 0, 100, 5), nrow = 1)

  # Higher dose should give lower value for b > 0 (decreasing)
  y_low <- mod$fct(0, parm)
  y_high <- mod$fct(100, parm)
  expect_true(as.numeric(y_low) > as.numeric(y_high))
})

test_that("gompertz fct increasing curve when b < 0", {
  mod <- gompertz()
  parm <- matrix(c(-1, 0, 100, 5), nrow = 1)

  # Lower dose should give lower value for b < 0 (increasing)
  y_low <- mod$fct(0, parm)
  y_high <- mod$fct(100, parm)
  expect_true(as.numeric(y_low) < as.numeric(y_high))
})

# ========================================================================
# Test: deriv1 (first derivatives in parameters)
# ========================================================================

test_that("gompertz deriv1 returns correct dimensions with all parameters free", {
  mod <- gompertz()
  parm <- matrix(c(1, 0, 100, 5), nrow = 1)

  result <- mod$deriv1(5, parm)
  # Single row: gradient drops to vector
  expect_length(result, 4)
})

test_that("gompertz deriv1 returns correct dimensions with fixed parameters", {
  mod <- gompertz(fixed = c(NA, 0, NA, NA))
  parm <- matrix(c(1, 100, 5), nrow = 1)

  result <- mod$deriv1(5, parm)
  expect_length(result, 3)  # Only free parameters
})

test_that("gompertz deriv1 handles multiple rows", {
  mod <- gompertz()
  parm <- matrix(c(1, 0, 100, 5,
                    2, 10, 90, 3), nrow = 2, byrow = TRUE)
  doses <- c(5, 3)

  result <- mod$deriv1(doses, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 4)
})

test_that("gompertz deriv1 gradient at dose=e is correct", {
  mod <- gompertz()
  # b=1, c=0, d=100, e=5
  parm <- matrix(c(1, 0, 100, 5), nrow = 1)

  result <- mod$deriv1(5, parm)
  # At dose=e=5: exp(b*(dose-e)) = exp(0) = 1, exp(-1) = 1/e
  # dc/dc = 1 - exp(-1) 
  # dd/dd = exp(-1)
  expect_equal(as.numeric(result[2]), 1 - exp(-1), tolerance = 1e-10)  # dc
  expect_equal(as.numeric(result[3]), exp(-1), tolerance = 1e-10)  # dd
})

# ========================================================================
# Test: derivx (first derivatives in x)
# ========================================================================

test_that("gompertz derivx returns correct dimensions", {
  mod <- gompertz()
  parm <- matrix(c(1, 0, 100, 5), nrow = 1)

  result <- mod$derivx(5, parm)
  expect_length(result, 1)
})

test_that("gompertz derivx at dose=e has expected value", {
  mod <- gompertz()
  # b=1, c=0, d=100, e=5
  parm <- matrix(c(1, 0, 100, 5), nrow = 1)

  # df/dx = -(d-c)*exp(-exp(b*(x-e)))*exp(b*(x-e))*b
  # At x=e: -(100)*exp(-1)*1*1 = -100/e
  result <- mod$derivx(5, parm)
  expected <- -100 * exp(-1) * 1
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("gompertz derivx handles fixed parameters", {
  mod <- gompertz(fixed = c(NA, 0, NA, NA))
  parm <- matrix(c(1, 100, 5), nrow = 1)

  result <- mod$derivx(5, parm)
  expected <- -100 * exp(-1) * 1
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("gompertz derivx handles multiple rows", {
  mod <- gompertz()
  parm <- matrix(c(1, 0, 100, 5,
                    2, 10, 90, 3), nrow = 2, byrow = TRUE)

  result <- mod$derivx(c(5, 3), parm)
  expect_length(result, 2)
})

# ========================================================================
# Test: edfct (effective dose function)
# ========================================================================

test_that("edfct returns ED correctly for relative type", {
  mod <- gompertz()
  # b=1, c=0, d=100, e=5
  # For Gompertz: ED = e + log(-log((100-p)/100)) / b
  parms <- c(1, 0, 100, 5)

  # ED50: e + log(-log(0.5)) / b = 5 + log(-log(0.5)) / 1
  result50 <- mod$edfct(parms, 50, reference = "control", type = "relative")
  expected50 <- 5 + log(-log(0.5)) / 1
  expect_equal(result50[[1]], expected50, tolerance = 1e-10)
  expect_true(is.numeric(result50[[2]]))
  expect_length(result50[[2]], 4)
})

test_that("edfct returns correct ED values for various p", {
  mod <- gompertz()
  parms <- c(1, 0, 100, 5)

  # ED10: e + log(-log(0.9)) / b
  result10 <- mod$edfct(parms, 10, reference = "control", type = "relative")
  expected10 <- 5 + log(-log(0.9)) / 1
  expect_equal(result10[[1]], expected10, tolerance = 1e-10)

  # ED90: e + log(-log(0.1)) / b
  result90 <- mod$edfct(parms, 90, reference = "control", type = "relative")
  expected90 <- 5 + log(-log(0.1)) / 1
  expect_equal(result90[[1]], expected90, tolerance = 1e-10)
})

test_that("edfct returns gradients with correct length when parameters are fixed", {
  mod <- gompertz(fixed = c(NA, 0, NA, NA))
  # Free: b, d, e (3 params)
  result <- mod$edfct(c(1, 100, 5), 50, reference = "control", type = "relative")

  expect_true(is.numeric(result[[1]]))
  expect_length(result[[2]], 3)  # gradient for 3 free parameters
})

test_that("edfct gradient for 'e' equals EDp (current implementation)", {
  mod <- gompertz()
  result <- mod$edfct(c(1, 0, 100, 5), 50, reference = "control", type = "relative")
  gradient <- result[[2]]
  EDp <- result[[1]]

  # In the current implementation, EDder = EDp * c(-tempVal/b^2, 0, 0, 1)
  # So gradient[4] = EDp * 1 = EDp
  expect_equal(gradient[4], EDp, tolerance = 1e-10)
})

test_that("edfct gradient for 'c' and 'd' is 0", {
  mod <- gompertz()
  result <- mod$edfct(c(1, 0, 100, 5), 50, reference = "control", type = "relative")
  gradient <- result[[2]]

  # c and d do not appear in ED formula directly
  expect_equal(gradient[2], 0)
  expect_equal(gradient[3], 0)
})

test_that("edfct works with negative b (increasing curve)", {
  mod <- gompertz()
  parms <- c(-1, 0, 100, 5)

  # For increasing curve (b < 0), EDhelper swaps p → 100-p
  result <- mod$edfct(parms, 50, reference = "control", type = "relative")
  expect_true(is.numeric(result[[1]]))
  expect_length(result[[2]], 4)
})

# ========================================================================
# Test: lowerAs, upperAs, monoton functions
# ========================================================================

test_that("lowerAs returns correct parameter", {
  mod <- gompertz()
  # lowerAs picks parameter 2 (c)
  result <- mod$lowerAs(c(1, 0, 100, 5))
  expect_equal(result, 0)
})

test_that("upperAs returns correct parameter", {
  mod <- gompertz()
  # upperAs picks parameter 3 (d)
  result <- mod$upperAs(c(1, 0, 100, 5))
  expect_equal(result, 100)
})

test_that("monoton returns -1 * b", {
  mod <- gompertz()
  # monoton picks parameter 1 (b) with sign -1
  result <- mod$monoton(c(1, 0, 100, 5))
  expect_equal(result, -1)
})

# ========================================================================
# Test: G.2() (two-parameter Gompertz, c=0, d=upper fixed)
# ========================================================================

test_that("G.2() returns gompertz object with correct structure", {
  result <- G.2()
  expect_s3_class(result, "gompertz")
  expect_equal(result$name, "G.2")
  expect_equal(result$text, "Gompertz with lower limit at 0 and upper limit at 1")
  expect_equal(result$noParm, 2)
  expect_equal(result$names, c("b", "e"))
})

test_that("G.2() errors on invalid names argument", {
  expect_error(G.2(names = c("a")), "Not correct 'names' argument")
  expect_error(G.2(names = 123), "Not correct 'names' argument")
})

test_that("G.2() errors on invalid fixed argument", {
  expect_error(G.2(fixed = c(NA)), "Not correct length of 'fixed' argument")
  expect_error(G.2(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("G.2() custom upper limit", {
  result <- G.2(upper = 50)
  expect_equal(result$text, "Gompertz with lower limit at 0 and upper limit at 50")
  # Check fct produces expected value
  parm <- matrix(c(1, 5), nrow = 1)  # b, e
  res <- result$fct(5, parm)
  expected <- 50 * exp(-exp(0))  # c=0, d=50
  expect_equal(as.numeric(res), expected, tolerance = 1e-10)
})

test_that("G.2() fct gives expected response", {
  mod <- G.2()
  # b, e (c=0, d=1 are fixed)
  parm <- matrix(c(1, 5), nrow = 1)
  result <- mod$fct(5, parm)
  expected <- 1 * exp(-exp(0))  # c=0, d=1
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

# ========================================================================
# Test: G.3() (three-parameter Gompertz, c=0 fixed)
# ========================================================================

test_that("G.3() returns gompertz object with correct structure", {
  result <- G.3()
  expect_s3_class(result, "gompertz")
  expect_equal(result$name, "G.3")
  expect_equal(result$text, "Gompertz with lower limit at 0")
  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("b", "d", "e"))
})

test_that("G.3() errors on invalid names argument", {
  expect_error(G.3(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(G.3(names = 123), "Not correct 'names' argument")
})

test_that("G.3() errors on invalid fixed argument", {
  expect_error(G.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
  expect_error(G.3(fixed = c(NA, NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("G.3() fct gives expected response", {
  mod <- G.3()
  # b, d, e (c=0 fixed)
  parm <- matrix(c(1, 100, 5), nrow = 1)
  result <- mod$fct(5, parm)
  expected <- 100 * exp(-exp(0))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("G.3() edfct gives correct ED", {
  mod <- G.3()
  result <- mod$edfct(c(1, 100, 5), 50, reference = "control", type = "relative")
  expect_true(is.numeric(result[[1]]))
  expect_length(result[[2]], 3)
})

# ========================================================================
# Test: G.3u() (three-parameter Gompertz, d=upper fixed)
# ========================================================================

test_that("G.3u() returns gompertz object with correct structure", {
  result <- G.3u()
  expect_s3_class(result, "gompertz")
  expect_equal(result$name, "G.3u")
  expect_equal(result$text, "Gompertz with upper limit at 1")
  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("b", "c", "e"))
})

test_that("G.3u() errors on invalid names argument", {
  expect_error(G.3u(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(G.3u(names = 123), "Not correct 'names' argument")
})

test_that("G.3u() errors on invalid fixed argument", {
  expect_error(G.3u(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
  expect_error(G.3u(fixed = c(NA, NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("G.3u() custom upper limit", {
  result <- G.3u(upper = 200)
  expect_equal(result$text, "Gompertz with upper limit at 200")
})

test_that("G.3u() fct gives expected response", {
  mod <- G.3u()
  # b, c, e (d=1 fixed)
  parm <- matrix(c(1, 0, 5), nrow = 1)
  result <- mod$fct(5, parm)
  expected <- 0 + (1 - 0) * exp(-exp(1 * (5 - 5)))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

# ========================================================================
# Test: G.4() (four-parameter Gompertz)
# ========================================================================

test_that("G.4() returns gompertz object with correct structure", {
  result <- G.4()
  expect_s3_class(result, "gompertz")
  expect_equal(result$name, "G.4")
  expect_equal(result$noParm, 4)
  expect_equal(result$names, c("b", "c", "d", "e"))
})

test_that("G.4() errors on invalid names argument", {
  expect_error(G.4(names = c("a", "b")), "Not correct names argument")
  expect_error(G.4(names = 123), "Not correct names argument")
})

test_that("G.4() errors on invalid fixed argument", {
  expect_error(G.4(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
  expect_error(G.4(fixed = c(NA, NA, NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("G.4() fct gives expected response", {
  mod <- G.4()
  parm <- matrix(c(1, 0, 100, 5), nrow = 1)
  result <- mod$fct(5, parm)
  expected <- 100 * exp(-exp(0))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("G.4() passes additional arguments to gompertz()", {
  custom_ss <- function(dframe) { c(1, 0, 100, 5) }
  result <- G.4(ssfct = custom_ss)
  expect_identical(result$ssfct, custom_ss)
})

# ========================================================================
# Test: ssfct (self-starter function) via gompertz()
# ========================================================================

test_that("gompertz ssfct works with simple data frame", {
  mod <- gompertz()
  # Create a simple dose-response like data frame
  set.seed(42)
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  # Gompertz with b=0.5, c=10, d=100, e=5
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- mod$ssfct(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)  # 4 free parameters
  expect_true(all(is.finite(result)))
})

test_that("gompertz ssfct with method='2' works", {
  mod <- gompertz(method = "2")
  set.seed(42)
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- mod$ssfct(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)
})

test_that("gompertz ssfct with method='3' works", {
  mod <- gompertz(method = "3")
  set.seed(42)
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- mod$ssfct(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)
})

test_that("gompertz ssfct with method='4' works", {
  mod <- gompertz(method = "4")
  set.seed(42)
  doses <- c(0, 1, 2, 3, 5, 8, 10, 15, 20)
  responses <- 10 + (100 - 10) * exp(-exp(0.5 * (doses - 5)))
  dframe <- data.frame(dose = doses, response = responses)

  result <- mod$ssfct(dframe)
  expect_true(is.numeric(result))
  expect_length(result, 4)
})

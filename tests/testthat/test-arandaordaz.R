# tests/testthat/test-arandaordaz.R
# Comprehensive tests for the arandaordaz() function in R/arandaordaz.R.
# Note: The AR.2() and AR.3() convenience wrappers were removed from
# arandaordaz.R because they were dead code (overridden by R/weibull2.R).

# ============================================================================
# arandaordaz() — Input Validation (Error Handling)
# ============================================================================

test_that("arandaordaz errors on invalid 'names' argument", {
  # Non-character names
  expect_error(arandaordaz(names = c(1, 2, 3)), "Not correct 'names' argument")
  # Wrong length names
  expect_error(arandaordaz(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(arandaordaz(names = c("a", "b", "c", "d")), "Not correct 'names' argument")
})

test_that("arandaordaz errors on invalid 'fixed' argument", {
  expect_error(arandaordaz(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(arandaordaz(fixed = c(NA, NA, NA, NA)), "Not correct 'fixed' argument")
})

# ============================================================================
# arandaordaz() — Return Structure (Happy Path)
# ============================================================================

test_that("arandaordaz returns drcMean object with correct structure", {
  result <- arandaordaz()
  expect_s3_class(result, "drcMean")

  expected_names <- c("fct", "ssfct", "names", "deriv1", "deriv2", "derivx",
                       "edfct", "inversion", "name", "text", "noParm")
  expect_named(result, expected_names)

  expect_true(is.function(result$fct))
  expect_true(is.function(result$ssfct))
  expect_true(is.function(result$edfct))
  expect_true(is.function(result$inversion))
  expect_null(result$deriv1)
  expect_null(result$deriv2)
  expect_null(result$derivx)
})

test_that("arandaordaz default call uses correct name, text, and noParm", {
  result <- arandaordaz()
  expect_equal(result$name, "arandaordaz")
  expect_equal(result$text, "Asymptotic regression")
  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("a", "b", "c"))
})

test_that("arandaordaz uses fctName and fctText when provided", {
  result <- arandaordaz(fctName = "myName", fctText = "myText")
  expect_equal(result$name, "myName")
  expect_equal(result$text, "myText")
})

test_that("arandaordaz handles partially fixed parameters", {
  # Fix the first parameter 'a' at 0
  result <- arandaordaz(fixed = c(0, NA, NA))
  expect_equal(result$names, c("b", "c"))
  expect_equal(result$noParm, 2)

  # Fix two parameters
  result2 <- arandaordaz(fixed = c(0, 10, NA))
  expect_equal(result2$names, "c")
  expect_equal(result2$noParm, 1)

  # Fix all parameters (noParm = 0)
  result3 <- arandaordaz(fixed = c(0, 10, 0.5))
  expect_equal(result3$names, character(0))
  expect_equal(result3$noParm, 0)
})

# ============================================================================
# arandaordaz()$fct — Mean Function
# ============================================================================

test_that("fct computes correct mean values (all params free)", {
  result <- arandaordaz()
  # f(x) = a + (b - a) * (1 - exp(-c * x))
  parm <- matrix(c(0, 10, 0.5), nrow = 1)
  dose <- 1
  val <- result$fct(dose, parm)
  expected <- 0 + (10 - 0) * (1 - exp(-0.5 * 1))
  expect_equal(val, expected)
})

test_that("fct computes correct mean with fixed parameter", {
  # Fix a=0
  result <- arandaordaz(fixed = c(0, NA, NA))
  parm <- matrix(c(10, 0.5), nrow = 1)
  dose <- 1
  val <- result$fct(dose, parm)
  expected <- 0 + (10 - 0) * (1 - exp(-0.5 * 1))
  expect_equal(val, expected)
})

test_that("fct handles multiple doses and observations", {
  result <- arandaordaz()
  parm <- matrix(c(0, 10, 0.5, 1, 20, 1), nrow = 2, byrow = TRUE)
  dose <- c(1, 2)
  vals <- result$fct(dose, parm)
  expected1 <- 0 + (10 - 0) * (1 - exp(-0.5 * 1))
  expected2 <- 1 + (20 - 1) * (1 - exp(-1 * 2))
  expect_equal(vals, c(expected1, expected2))
})

test_that("fct returns lower limit when dose is 0", {
  result <- arandaordaz()
  parm <- matrix(c(2, 10, 0.5), nrow = 1)
  dose <- 0
  val <- result$fct(dose, parm)
  expect_equal(val, 2)
})

# ============================================================================
# arandaordaz()$ssfct — Self-Starter Function
# ============================================================================

test_that("ssfct returns starting values for all-free parameters", {
  result <- arandaordaz()
  x <- c(0.5, 1, 2, 3, 5, 10)
  y <- 10 * (1 - exp(-0.5 * x))
  dataf <- data.frame(x = x, y = y)
  starts <- result$ssfct(dataf)
  expect_length(starts, 3)
  expect_true(is.numeric(starts))
})

test_that("ssfct returns correct number of starting values with fixed params", {
  result <- arandaordaz(fixed = c(0, NA, NA))
  x <- c(0.5, 1, 2, 3, 5, 10)
  y <- 10 * (1 - exp(-0.5 * x))
  dataf <- data.frame(x = x, y = y)
  starts <- result$ssfct(dataf)
  expect_length(starts, 2)
})

# ============================================================================
# arandaordaz()$edfct — ED Function (all 4 branch combinations)
# ============================================================================

test_that("edfct with type='relative' and reference='upper'", {
  result <- arandaordaz()
  parm <- c(0, 10, 0.5)
  ed_result <- result$edfct(parm, respl = 50, reference = "upper", type = "relative")
  # p = 50, pProp = 0.5, EDp = -log(0.5)/0.5
  expect_length(ed_result, 2)
  expect_equal(ed_result[[1]], -log(0.5) / 0.5)
  expect_length(ed_result[[2]], 3)
})

test_that("edfct with type='absolute' and reference='upper'", {
  result <- arandaordaz()
  parm <- c(0, 10, 0.5)
  # type="absolute": p = 100*((b - respl)/(b - a)) = 100*(10-5)/(10-0) = 50
  ed_result <- result$edfct(parm, respl = 5, reference = "upper", type = "absolute")
  expect_equal(ed_result[[1]], -log(0.5) / 0.5)
})

test_that("edfct with type='relative' and reference='control'", {
  result <- arandaordaz()
  parm <- c(0, 10, 0.5)
  # p = 50, then p = 100 - 50 = 50 (control reference)
  ed_result <- result$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_length(ed_result, 2)
  expect_true(is.numeric(ed_result[[1]]))
})

test_that("edfct with type='absolute' and reference='control'", {
  result <- arandaordaz()
  parm <- c(0, 10, 0.5)
  ed_result <- result$edfct(parm, respl = 5, reference = "control", type = "absolute")
  expect_length(ed_result, 2)
  expect_true(is.numeric(ed_result[[1]]))
})

test_that("edfct returns correct derivative length with fixed params", {
  result <- arandaordaz(fixed = c(0, NA, NA))
  parm <- c(10, 0.5)
  ed_result <- result$edfct(parm, respl = 50, reference = "upper", type = "relative")
  # Only b and c are free, so derivative vector should have 2 elements
  expect_length(ed_result[[2]], 2)
})

# ============================================================================
# arandaordaz()$inversion — Inverse Function
# ============================================================================

test_that("inversion function correctly inverts the mean function", {
  result <- arandaordaz()
  parm <- c(0, 10, 0.5)
  # inv(y) = -log((y - a) / (b - a)) / c
  inv_val <- result$inversion(5, parm)
  expected <- -log((5 - 0) / (10 - 0)) / 0.5
  expect_equal(inv_val, expected)
})

test_that("inversion round-trips with fct", {
  result <- arandaordaz()
  parm_vec <- c(2, 20, 0.3)
  dose <- 3
  parm_mat <- matrix(parm_vec, nrow = 1)
  y <- result$fct(dose, parm_mat)
  x_back <- result$inversion(y, parm_vec)
  expect_equal(x_back, dose, tolerance = 1e-10)
})

test_that("inversion works with fixed parameters", {
  result <- arandaordaz(fixed = c(0, NA, NA))
  parm <- c(10, 0.5)
  inv_val <- result$inversion(5, parm)
  expected <- -log((5 - 0) / (10 - 0)) / 0.5
  expect_equal(inv_val, expected)
})

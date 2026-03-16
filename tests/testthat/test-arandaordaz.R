# Tests for arandaordaz.R: arandaordaz() function

# --- Input Validation: Error Handling ---

test_that("arandaordaz errors when fixed is not numeric", {
  expect_error(
    arandaordaz(fixed = c("a", "b", "c")),
    "'fixed' must be a numeric vector"
  )
  expect_error(
    arandaordaz(fixed = list(NA, NA, NA)),
    "'fixed' must be a numeric vector"
  )
})

test_that("arandaordaz errors when fixed has wrong length", {
  expect_error(
    arandaordaz(fixed = c(NA, NA)),
    "'fixed' must have length 3"
  )
  expect_error(
    arandaordaz(fixed = c(NA, NA, NA, NA)),
    "'fixed' must have length 3"
  )
  expect_error(
    arandaordaz(fixed = numeric(0)),
    "'fixed' must have length 3"
  )
})

test_that("arandaordaz errors when names is not character", {
  expect_error(
    arandaordaz(names = c(1, 2, 3)),
    "'names' must be a character vector of length 3"
  )
  expect_error(
    arandaordaz(names = list("a", "b", "c")),
    "'names' must be a character vector of length 3"
  )
})

test_that("arandaordaz errors when names has wrong length", {
  expect_error(
    arandaordaz(names = c("a", "b")),
    "'names' must be a character vector of length 3"
  )
  expect_error(
    arandaordaz(names = c("a", "b", "c", "d")),
    "'names' must be a character vector of length 3"
  )
  expect_error(
    arandaordaz(names = character(0)),
    "'names' must be a character vector of length 3"
  )
})

# --- Correctness: Happy Path ---

test_that("arandaordaz returns correct class and structure with defaults", {
  ar <- arandaordaz()
  expect_s3_class(ar, "drcMean")
  expect_true(is.list(ar))
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2", "derivx",
                     "edfct", "inversion", "name", "text", "noParm") %in% names(ar)))
})

test_that("arandaordaz has correct default parameter names", {
  ar <- arandaordaz()
  expect_equal(ar$names, c("a", "b", "c"))
})

test_that("arandaordaz has correct noParm with no fixed parameters", {
  ar <- arandaordaz()
  expect_equal(ar$noParm, 3)
})

test_that("arandaordaz uses default name when fctName not provided", {
  ar <- arandaordaz()
  expect_equal(ar$name, "arandaordaz")
})

test_that("arandaordaz uses default text when fctText not provided", {
  ar <- arandaordaz()
  expect_equal(ar$text, "Asymptotic regression")
})

test_that("arandaordaz uses provided fctName", {
  ar <- arandaordaz(fctName = "CustomName")
  expect_equal(ar$name, "CustomName")
})

test_that("arandaordaz uses provided fctText", {
  ar <- arandaordaz(fctText = "Custom text description")
  expect_equal(ar$text, "Custom text description")
})

test_that("arandaordaz returns invisible result", {
  result <- withVisible(arandaordaz())
  expect_false(result$visible)
})

# --- Fixed Parameters ---

test_that("arandaordaz handles fixed parameters correctly", {
  ar <- arandaordaz(fixed = c(1, NA, NA))
  expect_equal(ar$noParm, 2)
  expect_equal(ar$names, c("b", "c"))
})

test_that("arandaordaz handles multiple fixed parameters", {
  ar <- arandaordaz(fixed = c(1, 10, NA))
  expect_equal(ar$noParm, 1)
  expect_equal(ar$names, c("c"))
})

test_that("arandaordaz handles all parameters fixed", {
  ar <- arandaordaz(fixed = c(1, 10, 0.5))
  expect_equal(ar$noParm, 0)
  expect_equal(length(ar$names), 0)
})

test_that("arandaordaz custom names are preserved for non-fixed parameters", {
  ar <- arandaordaz(fixed = c(NA, 5, NA), names = c("alpha", "beta", "gamma"))
  expect_equal(ar$names, c("alpha", "gamma"))
})

# --- Mean Function (fct) ---

test_that("arandaordaz fct computes correct values", {
  ar <- arandaordaz()
  # Parameters: a=0, b=10, c=1
  dose <- c(0, 1, 5, 10)
  parm <- matrix(c(0, 10, 1), nrow = 4, ncol = 3, byrow = TRUE)
  result <- ar$fct(dose, parm)
  # f(x) = a + (b-a)(1-exp(-c*x))
  expected <- 0 + (10 - 0) * (1 - exp(-1 * dose))
  expect_equal(result, expected)
})

test_that("arandaordaz fct works with fixed parameters", {
  ar <- arandaordaz(fixed = c(0, NA, NA))
  dose <- c(0, 1, 5)
  parm <- matrix(c(10, 1), nrow = 3, ncol = 2, byrow = TRUE)
  result <- ar$fct(dose, parm)
  expected <- 0 + (10 - 0) * (1 - exp(-1 * dose))
  expect_equal(result, expected)
})

test_that("arandaordaz fct handles zero dose", {
  ar <- arandaordaz()
  dose <- 0
  parm <- matrix(c(1, 5, 0.5), nrow = 1, ncol = 3)
  result <- ar$fct(dose, parm)
  # At dose=0: f(0) = a + (b-a)(1-exp(0)) = a + (b-a)*0 = a
  expect_equal(result, 1)
})

# --- Self-Starter Function (ssfct) ---

test_that("arandaordaz ssfct returns initial parameter estimates", {
  ar <- arandaordaz()
  # Create simple test data
  dataf <- data.frame(
    x = c(0, 1, 2, 5, 10),
    y = c(1, 2, 3, 5, 7)
  )
  init_params <- ar$ssfct(dataf)
  expect_equal(length(init_params), 3)
  expect_true(all(is.finite(init_params)))
})

test_that("arandaordaz ssfct works with fixed parameters", {
  ar <- arandaordaz(fixed = c(NA, NA, 0.5))
  dataf <- data.frame(
    x = c(0, 1, 2, 5, 10),
    y = c(1, 2, 3, 5, 7)
  )
  init_params <- ar$ssfct(dataf)
  expect_equal(length(init_params), 2)
  expect_true(all(is.finite(init_params)))
})

test_that("arandaordaz ssfct warns on invalid log argument", {
  ar <- arandaordaz()
  # Create data that will trigger warning
  # Need y values that cause innerVal <= 0
  # innerVal = -((y - aPar) / (bPar - aPar) - 1)
  # aPar = min(y) * 0.95, bPar = max(y) * 1.05
  # For y = min(y), innerVal = -(0.95/1 - 1) = -(-0.05) = 0.05
  # For y = max(y), innerVal = -(1/1.05 - 1) = -(-0.0476) ≈ 0.0476
  # Need y values outside this range to trigger warning
  # If y > bPar, then (y - aPar)/(bPar - aPar) > 1, so innerVal < 0
  dataf <- data.frame(
    x = c(0, 1, 2, 3, 4),
    y = c(1, 2, 3, 4, 10)  # Last value much higher
  )
  expect_warning(
    ar$ssfct(dataf),
    "Self-starter encountered invalid log argument"
  )
})

# --- Effective Dose Function (edfct) ---

test_that("arandaordaz edfct works with relative type and control reference", {
  ar <- arandaordaz()
  parm <- c(0, 10, 1)
  result <- ar$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_true(is.numeric(result[[1]]))
  expect_true(is.numeric(result[[2]]))
})

test_that("arandaordaz edfct works with absolute type and control reference", {
  ar <- arandaordaz()
  parm <- c(0, 10, 1)
  result <- ar$edfct(parm, respl = 5, reference = "control", type = "absolute")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("arandaordaz edfct works with relative type and non-control reference", {
  ar <- arandaordaz()
  parm <- c(0, 10, 1)
  result <- ar$edfct(parm, respl = 50, reference = "upper", type = "relative")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("arandaordaz edfct works with absolute type and non-control reference", {
  ar <- arandaordaz()
  parm <- c(0, 10, 1)
  result <- ar$edfct(parm, respl = 5, reference = "upper", type = "absolute")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("arandaordaz edfct works with fixed parameters", {
  ar <- arandaordaz(fixed = c(0, NA, NA))
  parm <- c(10, 1)  # Only non-fixed parameters
  result <- ar$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_equal(length(result[[2]]), 2)  # Derivatives only for non-fixed
})

# --- Inverse Function (inversion) ---

test_that("arandaordaz inversion computes correct values", {
  ar <- arandaordaz()
  parm <- c(0, 10, 1)
  y <- c(2, 5, 8)
  result <- ar$inversion(y, parm)
  expect_true(all(is.finite(result)))
  expect_equal(length(result), 3)
})

test_that("arandaordaz inversion works with fixed parameters", {
  ar <- arandaordaz(fixed = c(0, NA, NA))
  parm <- c(10, 1)
  y <- c(2, 5, 8)
  result <- ar$inversion(y, parm)
  expect_true(all(is.finite(result)))
})

test_that("arandaordaz inversion is inverse of fct", {
  ar <- arandaordaz()
  parm <- c(1, 10, 0.5)
  dose <- c(1, 5, 10)

  # Calculate response from dose
  parm_mat <- matrix(parm, nrow = 3, ncol = 3, byrow = TRUE)
  y <- ar$fct(dose, parm_mat)

  # Back-calculate dose from response
  dose_back <- ar$inversion(y, parm)

  expect_equal(dose, dose_back, tolerance = 1e-10)
})

# --- Derivative Slots ---

test_that("arandaordaz deriv1 is NULL", {
  ar <- arandaordaz()
  expect_null(ar$deriv1)
})

test_that("arandaordaz deriv2 is NULL", {
  ar <- arandaordaz()
  expect_null(ar$deriv2)
})

test_that("arandaordaz derivx is NULL", {
  ar <- arandaordaz()
  expect_null(ar$derivx)
})

# --- Edge Cases ---

test_that("arandaordaz handles NA in fixed at different positions", {
  ar1 <- arandaordaz(fixed = c(1, NA, NA))
  expect_equal(ar1$names, c("b", "c"))

  ar2 <- arandaordaz(fixed = c(NA, 5, NA))
  expect_equal(ar2$names, c("a", "c"))

  ar3 <- arandaordaz(fixed = c(NA, NA, 0.5))
  expect_equal(ar3$names, c("a", "b"))
})

test_that("arandaordaz handles numeric fixed values including zero", {
  ar <- arandaordaz(fixed = c(0, NA, NA))
  expect_equal(ar$noParm, 2)
  expect_true(all(c("fct", "ssfct") %in% names(ar)))
})

test_that("arandaordaz handles negative fixed values", {
  ar <- arandaordaz(fixed = c(-5, NA, NA))
  expect_equal(ar$noParm, 2)
  expect_s3_class(ar, "drcMean")
})

test_that("arandaordaz ssfct handles single data point", {
  ar <- arandaordaz()
  dataf <- data.frame(x = 1, y = 5)
  init_params <- ar$ssfct(dataf)
  expect_equal(length(init_params), 3)
})

test_that("arandaordaz ssfct handles identical y values", {
  ar <- arandaordaz()
  dataf <- data.frame(
    x = c(0, 1, 2, 3),
    y = c(5, 5, 5, 5)
  )
  # This will produce aPar = 4.75, bPar = 5.25
  # All y values will be the same, potentially causing issues
  init_params <- ar$ssfct(dataf)
  expect_equal(length(init_params), 3)
})

# tests/testthat/test-multi2.R
# Comprehensive tests for R/multi2.R: multi2() and all nested functions
# (fct, fd, dFct, deriv1, derivx, edfct, ssfct)

# ========================================================================
# Test: multi2() argument validation
# ========================================================================

test_that("multi2() errors on invalid 'names' argument", {
  expect_error(multi2(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(multi2(names = 123), "Not correct 'names' argument")
  expect_error(multi2(names = c("a", "b", "c", "d")), "Not correct 'names' argument")
})

test_that("multi2() errors on invalid 'fixed' argument", {
  expect_error(multi2(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(multi2(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
  expect_error(multi2(fixed = c(NA, NA, NA, NA, NA, NA)), "Not correct 'fixed' argument")
})

# ========================================================================
# Test: multi2() return structure
# ========================================================================

test_that("multi2() returns object of class 'multistage'", {
  result <- multi2()
  expect_s3_class(result, "multistage")
})

test_that("multi2() return list has correct structure", {
  result <- multi2()
  expect_true(is.function(result$fct))
  expect_true(is.function(result$ssfct))
  expect_true(is.function(result$deriv1))
  expect_null(result$deriv2)
  expect_true(is.function(result$derivx))
  expect_true(is.function(result$edfct))
  expect_equal(result$noParm, 5)
  expect_equal(result$names, c("b1", "b2", "b3", "c", "d"))
})

test_that("multi2() default name and text are correct", {
  result <- multi2()
  expect_equal(result$name, "multi2")
  expect_equal(result$text, "Multistage")
})

test_that("multi2() custom fctName and fctText override defaults", {
  result <- multi2(fctName = "myModel", fctText = "my description")
  expect_equal(result$name, "myModel")
  expect_equal(result$text, "my description")
})

test_that("multi2() custom names are applied", {
  result <- multi2(names = c("p1", "p2", "p3", "p4", "p5"))
  expect_equal(result$names, c("p1", "p2", "p3", "p4", "p5"))
})

# ========================================================================
# Test: multi2() with fixed parameters
# ========================================================================

test_that("multi2() with some fixed parameters updates noParm and names", {
  result <- multi2(fixed = c(0, NA, NA, NA, NA))
  expect_equal(result$noParm, 4)
  expect_equal(result$names, c("b2", "b3", "c", "d"))

  result2 <- multi2(fixed = c(0, NA, 0, NA, 1))
  expect_equal(result2$noParm, 2)
  expect_equal(result2$names, c("b2", "c"))
})

# ========================================================================
# Test: fct (mean function) calculations
# ========================================================================

test_that("multi2 fct computes correct model values", {
  m <- multi2()
  # f(x) = c + (d-c)*(1 - exp(-b1 - b2*x - b3*x^2))
  # With b1=0, b2=1, b3=0, c=0, d=100:
  #   f(0) = 0 + 100*(1 - exp(0 - 0 - 0)) = 100*(1-1) = 0
  #   f(1) = 0 + 100*(1 - exp(0 - 1 - 0)) = 100*(1 - exp(-1)) ≈ 63.21
  dose <- c(0, 1, 2)
  parm <- matrix(c(0, 1, 0, 0, 100), nrow = 1)
  result <- m$fct(dose, parm)
  expected <- 0 + (100 - 0) * (1 - exp(-0 - 1 * dose - 0 * dose^2))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("multi2 fct works with multiple observations (rows)", {
  m <- multi2()
  dose <- c(0, 1, 2, 3)
  parm <- matrix(c(0, 1, 0, 0, 100), nrow = 4, ncol = 5, byrow = TRUE)
  result <- m$fct(dose, parm)
  expected <- 0 + 100 * (1 - exp(-0 - 1 * dose - 0 * dose^2))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("multi2 fct handles fixed parameters correctly", {
  # Fix b1=0 and b3=0
  m <- multi2(fixed = c(0, NA, 0, NA, NA))
  dose <- c(0, 1, 2)
  # Only free parameters: b2, c, d
  parm <- matrix(c(1, 0, 100), nrow = 1)
  result <- m$fct(dose, parm)
  expected <- 0 + 100 * (1 - exp(-0 - 1 * dose - 0 * dose^2))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("multi2 fct with quadratic term produces correct results", {
  m <- multi2()
  dose <- c(0, 0.5, 1, 2)
  # b1=0.1, b2=0.5, b3=0.2, c=10, d=90
  parm <- matrix(c(0.1, 0.5, 0.2, 10, 90), nrow = 1)
  result <- m$fct(dose, parm)
  expected <- 10 + (90 - 10) * (1 - exp(-0.1 - 0.5 * dose - 0.2 * dose^2))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

# ========================================================================
# Test: deriv1 (parameter derivatives)
# ========================================================================

test_that("multi2 deriv1 returns gradient matrix with correct dimensions", {
  m <- multi2()
  dose <- c(0, 1, 2)
  parm <- matrix(c(0, 1, 0, 0, 100), nrow = 1)
  grad <- m$deriv1(dose, parm)
  # Should have nrow = length(dose), ncol = 5 (all free parameters)
  expect_equal(nrow(grad), length(dose))
  expect_equal(ncol(grad), 5)
})

test_that("multi2 deriv1 with fixed parameters returns reduced gradient", {
  m <- multi2(fixed = c(0, NA, 0, NA, NA))
  dose <- c(0, 1, 2)
  # free: b2, c, d
  parm <- matrix(c(1, 0, 100), nrow = 1)
  grad <- m$deriv1(dose, parm)
  expect_equal(nrow(grad), length(dose))
  expect_equal(ncol(grad), 3)
})

test_that("multi2 deriv1 computes correct gradient values (numerical check)", {
  m <- multi2()
  dose <- c(1)
  b1 <- 0; b2 <- 1; b3 <- 0.5; cc <- 0; dd <- 100
  parm <- matrix(c(b1, b2, b3, cc, dd), nrow = 1)
  grad <- m$deriv1(dose, parm)

  # Numerical derivatives via finite differences
  eps <- 1e-7
  params <- c(b1, b2, b3, cc, dd)
  numgrad <- numeric(5)
  for (i in 1:5) {
    params_up <- params
    params_up[i] <- params_up[i] + eps
    f_up <- cc + (dd - cc) * (1 - exp(-params_up[1] - params_up[2] * dose - params_up[3] * dose^2))
    # Need to recompute with modified c or d properly
    if (i == 4) {
      f_up <- params_up[4] + (dd - params_up[4]) * (1 - exp(-b1 - b2 * dose - b3 * dose^2))
    }
    if (i == 5) {
      f_up <- cc + (params_up[5] - cc) * (1 - exp(-b1 - b2 * dose - b3 * dose^2))
    }
    f_base <- cc + (dd - cc) * (1 - exp(-b1 - b2 * dose - b3 * dose^2))
    numgrad[i] <- (f_up - f_base) / eps
  }
  expect_equal(as.numeric(grad), numgrad, tolerance = 1e-5)
})

# ========================================================================
# Test: derivx (dose derivative)
# ========================================================================

test_that("multi2 derivx returns gradient matrix with correct dimensions", {
  m <- multi2()
  dose <- c(0, 1, 2)
  parm <- matrix(c(0, 1, 0, 0, 100), nrow = 1)
  dfdx <- m$derivx(dose, parm)
  expect_equal(nrow(dfdx), length(dose))
  expect_equal(ncol(dfdx), 1)
})

test_that("multi2 derivx computes correct dose derivative (numerical check)", {
  m <- multi2()
  dose <- c(1)
  b1 <- 0; b2 <- 1; b3 <- 0.5; cc <- 0; dd <- 100
  parm <- matrix(c(b1, b2, b3, cc, dd), nrow = 1)
  dfdx <- m$derivx(dose, parm)

  # Numerical derivative w.r.t. dose
  eps <- 1e-7
  f <- function(x) cc + (dd - cc) * (1 - exp(-b1 - b2 * x - b3 * x^2))
  numderiv <- (f(dose + eps) - f(dose)) / eps
  expect_equal(as.numeric(dfdx), numderiv, tolerance = 1e-5)
})

test_that("multi2 derivx with fixed parameters", {
  m <- multi2(fixed = c(0, NA, 0, NA, NA))
  dose <- c(1)
  parm <- matrix(c(1, 0, 100), nrow = 1)
  dfdx <- m$derivx(dose, parm)
  expect_equal(nrow(dfdx), 1)
  expect_equal(ncol(dfdx), 1)
})

# ========================================================================
# Test: edfct (effective dose)
# ========================================================================

test_that("multi2 edfct returns list with ED value and derivatives", {
  m <- multi2()
  # b1=0, b2=1, b3=0.1, c=0, d=100
  parm <- c(0, 1, 0.1, 0, 100)
  result <- m$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_true(is.numeric(result[[1]]))
  expect_true(is.numeric(result[[2]]))
})

test_that("multi2 edfct with type='absolute'", {
  m <- multi2()
  # b1=0, b2=1, b3=0.1, c=0, d=100
  parm <- c(0, 1, 0.1, 0, 100)
  result <- m$edfct(parm, respl = 50, reference = "control", type = "absolute")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("multi2 edfct with type='relative' and negative b1 and reference='control'", {
  m <- multi2()
  # b1=-0.5 (negative), b2=1, b3=0.1, c=0, d=100
  parm <- c(-0.5, 1, 0.1, 0, 100)
  result <- m$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("multi2 edfct with type='relative' and positive b1", {
  m <- multi2()
  # b1=0.5 (positive), b2=1, b3=0.1, c=0, d=100
  # This path does NOT trigger the `100 - p` reversal
  parm <- c(0.5, 1, 0.1, 0, 100)
  result <- m$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("multi2 edfct with type='relative' and reference != 'control'", {
  m <- multi2()
  parm <- c(-0.5, 1, 0.1, 0, 100)
  # reference is NOT "control", so second reversal should not apply
  result <- m$edfct(parm, respl = 50, reference = "upper", type = "relative")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("multi2 edfct with fixed parameters", {
  m <- multi2(fixed = c(0, NA, 0, NA, NA))
  # free params: b2, c, d
  parm <- c(1, 0, 100)
  result <- m$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result[[2]]), 3)  # 3 free parameters
})

# ========================================================================
# Test: ssfct (self-starter function)
# ========================================================================

test_that("multi2 default ssfct returns initial parameter estimates", {
  m <- multi2()
  # Create a simple dose-response data frame
  dose <- c(0, 0.5, 1, 2, 5, 10, 20)
  response <- c(0, 5, 15, 40, 70, 90, 98)
  dframe <- data.frame(dose, response)
  result <- m$ssfct(dframe)
  expect_true(is.numeric(result))
  expect_equal(length(result), 5)  # 5 free parameters
})

test_that("multi2 custom ssfct is used when provided", {
  custom_ss <- function(dframe) {
    rep(1, 5)
  }
  m <- multi2(ssfct = custom_ss)
  dframe <- data.frame(dose = c(0, 1, 10), response = c(0, 50, 100))
  result <- m$ssfct(dframe)
  expect_equal(result, rep(1, 5))
})

test_that("multi2 default ssfct with fixed parameters returns correct length", {
  m <- multi2(fixed = c(0, NA, 0, NA, NA))
  dose <- c(0, 0.5, 1, 2, 5, 10, 20)
  response <- c(0, 5, 15, 40, 70, 90, 98)
  dframe <- data.frame(dose, response)
  result <- m$ssfct(dframe)
  expect_true(is.numeric(result))
  expect_equal(length(result), 3)  # 3 free parameters
})

# ========================================================================
# Test: Model mathematical properties
# ========================================================================

test_that("multi2 model: at dose=0 with b1=0, response equals c (lower asymptote)", {
  m <- multi2()
  dose <- 0
  # b1=0, b2=1, b3=0.5, c=10, d=90
  parm <- matrix(c(0, 1, 0.5, 10, 90), nrow = 1)
  result <- m$fct(dose, parm)
  # f(0) = c + (d-c)*(1 - exp(-b1)) = 10 + 80*(1 - 1) = 10
  expect_equal(as.numeric(result), 10, tolerance = 1e-10)
})

test_that("multi2 model: at large dose, response approaches d (upper asymptote)", {
  m <- multi2()
  dose <- 1000
  # b1=0, b2=1, b3=0.5, c=10, d=90
  parm <- matrix(c(0, 1, 0.5, 10, 90), nrow = 1)
  result <- m$fct(dose, parm)
  # At large dose, exp(-b1-b2*x-b3*x^2) → 0, so f → c + (d-c)*1 = d = 90
  expect_equal(as.numeric(result), 90, tolerance = 1e-10)
})

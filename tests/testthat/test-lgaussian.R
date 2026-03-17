# Tests for lgaussian.R: lgaussian() function and its internal components

# --- lgaussian() main function: structure and class ---

test_that("lgaussian returns correct class and structure", {
  lg <- lgaussian()
  expect_s3_class(lg, "lgaussian")
  expect_true(is.list(lg))
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2", "derivx",
                     "edfct", "name", "text", "noParm", "fixed") %in% names(lg)))
})

test_that("lgaussian default names are b, c, d, e, f", {
  lg <- lgaussian()
  expect_equal(lg$names, c("b", "c", "d", "e", "f"))
})

test_that("lgaussian noParm reflects number of NA in fixed", {
  lg_full <- lgaussian()
  expect_equal(lg_full$noParm, 5)

  lg_partial <- lgaussian(fixed = c(1, NA, NA, NA, NA))
  expect_equal(lg_partial$noParm, 4)
  expect_equal(lg_partial$names, c("c", "d", "e", "f"))
})

test_that("lgaussian uses default text when fctText not provided", {
  lg <- lgaussian()
  expect_equal(lg$text, "Log-Gaussian")
})

test_that("lgaussian uses provided fctText", {
  lg <- lgaussian(fctText = "Custom text")
  expect_equal(lg$text, "Custom text")
})

test_that("lgaussian uses provided fctName", {
  lg <- lgaussian(fctName = "myFunc")
  expect_equal(lg$name, "myFunc")
})

test_that("lgaussian default name when fctName not provided", {
  lg <- lgaussian()
  expect_equal(lg$name, "lgaussian")
})

# --- Error handling ---

test_that("lgaussian errors on invalid names argument", {
  expect_error(lgaussian(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(lgaussian(names = 123), "Not correct 'names' argument")
})

test_that("lgaussian errors on invalid fixed argument", {
  expect_error(lgaussian(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(lgaussian(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
})

# --- Self-starter function (ssfct) ---

test_that("lgaussian uses provided ssfct when not NULL", {
  custom_ss <- function(dframe) { c(1, 0, 10, 5, 1) }
  lg <- lgaussian(ssfct = custom_ss)
  expect_identical(lg$ssfct, custom_ss)
})

test_that("lgaussian uses default ssfct when ssfct is NULL", {
  lg <- lgaussian(ssfct = NULL)
  expect_true(is.function(lg$ssfct))
})

test_that("lgaussian method argument works for self-starter", {
  for (m in c("1", "2", "3", "4")) {
    lg <- lgaussian(method = m)
    expect_true(is.function(lg$ssfct))
  }
})

# --- loge parameter (present for API compatibility) ---

test_that("lgaussian with loge=FALSE (default)", {
  lg <- lgaussian(loge = FALSE)
  expect_s3_class(lg, "lgaussian")
  expect_true(is.function(lg$fct))
})

test_that("lgaussian with loge=TRUE", {
  lg <- lgaussian(loge = TRUE)
  expect_s3_class(lg, "lgaussian")
  expect_true(is.function(lg$fct))
})

# --- fct (internal nonlinear function) ---

test_that("lgaussian fct computes correct values", {
  lg <- lgaussian()
  # Parameters: b=1, c=0, d=100, e=5, f=2
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5, 2), nrow = 3, ncol = 5, byrow = TRUE)
  result <- lg$fct(dose, parm)
  # Model: c + (d-c)*exp(-0.5*(sqrt(((log(dose)-log(e))/b)^2))^f)
  b <- 1; c_val <- 0; d_val <- 100; e_val <- 5; f_val <- 2
  expected <- c_val + (d_val - c_val) * exp(-0.5 * (sqrt(((log(dose) - log(e_val)) / b)^2))^f_val)
  expect_equal(as.numeric(result), expected)
})

test_that("lgaussian fct at peak dose equals d parameter", {
  lg <- lgaussian()
  # At dose=e, log(dose)-log(e)=0, so result = c + (d-c)*exp(0) = d
  dose <- c(5)
  parm <- matrix(c(1, 0, 100, 5, 2), nrow = 1, ncol = 5, byrow = TRUE)
  result <- lg$fct(dose, parm)
  expect_equal(as.numeric(result), 100)
})

test_that("lgaussian fct works with fixed parameters", {
  lg <- lgaussian(fixed = c(1, 0, NA, NA, 2))
  dose <- c(1, 5, 10)
  # Only d and e are free
  parm <- matrix(c(100, 5), nrow = 3, ncol = 2, byrow = TRUE)
  result <- lg$fct(dose, parm)
  b <- 1; c_val <- 0; d_val <- 100; e_val <- 5; f_val <- 2
  expected <- c_val + (d_val - c_val) * exp(-0.5 * (sqrt(((log(dose) - log(e_val)) / b)^2))^f_val)
  expect_equal(as.numeric(result), expected)
})

test_that("lgaussian fct has gradient attribute", {
  lg <- lgaussian()
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5, 2), nrow = 3, ncol = 5, byrow = TRUE)
  result <- lg$fct(dose, parm)
  expect_true(!is.null(attr(result, "gradient")))
})

# --- deriv1 (parameter derivatives) ---

test_that("lgaussian deriv1 returns gradient matrix with correct dimensions", {
  lg <- lgaussian()
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5, 2), nrow = 3, ncol = 5, byrow = TRUE)
  result <- lg$deriv1(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 5)  # all 5 params free
})

test_that("lgaussian deriv1 with fixed params reduces columns", {
  lg <- lgaussian(fixed = c(1, 0, NA, NA, 2))
  dose <- c(1, 5, 10)
  parm <- matrix(c(100, 5), nrow = 3, ncol = 2, byrow = TRUE)
  result <- lg$deriv1(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)  # only d and e are free
})

# --- derivx (dose derivatives) ---

test_that("lgaussian derivx returns gradient with correct dimensions", {
  lg <- lgaussian()
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5, 2), nrow = 3, ncol = 5, byrow = TRUE)
  result <- lg$derivx(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 1)
})

test_that("lgaussian derivx with fixed parameters", {
  lg <- lgaussian(fixed = c(1, 0, NA, NA, 2))
  dose <- c(1, 5, 10)
  parm <- matrix(c(100, 5), nrow = 3, ncol = 2, byrow = TRUE)
  result <- lg$derivx(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 1)
})

# --- edfct (effective dose function) ---

test_that("lgaussian edfct returns list with ED and gradient (relative)", {
  lg <- lgaussian()
  parm <- c(1, 0, 100, 5, 2)
  result <- lg$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_true(is.numeric(result[[1]]))
  # Gradient should have 5 elements (all params free)
  expect_equal(length(result[[2]]), 5)
})

test_that("lgaussian edfct works with absolute type", {
  lg <- lgaussian()
  parm <- c(1, 0, 100, 5, 2)
  result <- lg$edfct(parm, respl = 50, reference = "control", type = "absolute")
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_true(is.numeric(result[[1]]))
})

test_that("lgaussian edfct works with negative b and control reference (relative)", {
  lg <- lgaussian()
  parm <- c(-1, 0, 100, 5, 2)
  result <- lg$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_true(is.numeric(result[[1]]))
})

test_that("lgaussian edfct works with positive b and control reference (relative)", {
  lg <- lgaussian()
  parm <- c(1, 0, 100, 5, 2)
  result <- lg$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("lgaussian edfct with fixed parameters", {
  lg <- lgaussian(fixed = c(1, 0, NA, NA, 2))
  parm <- c(100, 5)
  result <- lg$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result[[2]]), 2)  # only free params
})

test_that("lgaussian edfct at different response levels", {
  lg <- lgaussian()
  parm <- c(1, 0, 100, 5, 2)
  result10 <- lg$edfct(parm, respl = 10, reference = "control", type = "relative")
  result90 <- lg$edfct(parm, respl = 90, reference = "control", type = "relative")
  expect_true(is.numeric(result10[[1]]))
  expect_true(is.numeric(result90[[1]]))
})

# --- lowerAs, upperAs, monoton ---

test_that("lgaussian lowerAs and upperAs return correct values", {
  lg <- lgaussian()
  # lowerAs extracts parameter 2 (c), upperAs extracts parameter 3 (d)
  parm <- c(1, 0, 100, 5, 2)
  expect_equal(lg$lowerAs(parm), 0)
  expect_equal(lg$upperAs(parm), 100)
})

test_that("lgaussian monoton is NA", {
  lg <- lgaussian()
  expect_true(is.na(lg$monoton))
})

# --- Fixed parameter edge cases ---

test_that("lgaussian with all parameters fixed", {
  lg <- lgaussian(fixed = c(1, 0, 100, 5, 2))
  expect_equal(lg$noParm, 0)
  expect_equal(length(lg$names), 0)
})

test_that("lgaussian deriv2 is NULL", {
  lg <- lgaussian()
  expect_null(lg$deriv2)
})

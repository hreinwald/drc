# Tests for llogistic2.R: llogistic2(), LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()
# and helper functions: lowFixed(), upFixed(), lowupFixed()

# ==============================================================================
# Test: llogistic2() argument validation
# ==============================================================================

test_that("llogistic2 errors on invalid 'names' argument", {
  expect_error(llogistic2(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(llogistic2(names = 123), "Not correct 'names' argument")
})

test_that("llogistic2 errors on invalid 'fixed' argument", {
  expect_error(llogistic2(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(llogistic2(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
})

# ==============================================================================
# Test: llogistic2() return structure
# ==============================================================================

test_that("llogistic2 returns object of class 'llogistic'", {
  result <- llogistic2()
  expect_s3_class(result, "llogistic")
  expect_true(is.list(result))
})

test_that("llogistic2 return list has correct structure", {
  result <- llogistic2()
  expect_true(is.function(result$fct))
  expect_true(is.function(result$ssfct))
  expect_true(is.function(result$deriv1))
  expect_null(result$deriv2)
  expect_true(is.function(result$derivx))
  expect_true(is.function(result$edfct))
  expect_true(is.function(result$inversion))
  expect_true(is.function(result$bfct))
  expect_equal(result$noParm, 5)
  expect_equal(result$names, c("b", "c", "d", "e", "f"))
})

test_that("llogistic2 default name and text are correct", {
  result <- llogistic2()
  expect_equal(result$name, "llogistic2")
  expect_equal(result$text, "Log-logistic (log(ED50) as parameter)")
})

test_that("llogistic2 custom fctName and fctText override defaults", {
  result <- llogistic2(fctName = "myModel", fctText = "my description")
  expect_equal(result$name, "myModel")
  expect_equal(result$text, "my description")
})

test_that("llogistic2 handles fixed parameters correctly", {
  result <- llogistic2(fixed = c(NA, 0, NA, NA, 1))
  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("b", "d", "e"))
})

test_that("llogistic2 uses custom ssfct when provided", {
  custom_ss <- function(dframe) { c(1, 0, 100, 5, 1) }
  result <- llogistic2(ssfct = custom_ss)
  expect_identical(result$ssfct, custom_ss)
})

# ==============================================================================
# Test: bfct (basic nonlinear function)
# ==============================================================================

test_that("llogistic2 bfct produces expected values", {
  mod <- llogistic2()
  # parm: b, c, d, e, f
  # bfct(x, parm) = parm[2] + (parm[3]-parm[2])/((1+(x/exp(parm[4]))^parm[1]))^parm[5]
  # With b=1, c=0, d=100, e=log(5), f=1:
  # At x=5: 0 + 100/((1+(5/5)^1))^1 = 100/2 = 50
  result <- mod$bfct(5, c(1, 0, 100, log(5), 1))
  expect_equal(result, 50, tolerance = 1e-10)
})

# ==============================================================================
# Test: fct (nonlinear function)
# ==============================================================================

test_that("llogistic2 fct produces expected values with all parameters free", {
  mod <- llogistic2()
  # Parameters: b, c, d, e(log scale), f
  # fct = c + (d-c)/((1+exp(b*(log(dose)-e)))^f)
  # With b=1, c=0, d=100, e=log(5), f=1:
  # At dose=5: 0 + 100/((1+exp(1*(log(5)-log(5))))^1) = 100/(1+1) = 50
  parm <- matrix(c(1, 0, 100, log(5), 1), nrow = 1)
  result <- mod$fct(5, parm)
  expect_equal(result, 50, tolerance = 1e-10)
})

test_that("llogistic2 fct handles fixed parameters", {
  mod <- llogistic2(fixed = c(NA, 0, NA, NA, 1))
  parm <- matrix(c(1, 100, log(5)), nrow = 1)
  result <- mod$fct(5, parm)
  expect_equal(result, 50, tolerance = 1e-10)
})

test_that("llogistic2 fct handles multiple doses", {
  mod <- llogistic2()
  parm <- matrix(rep(c(1, 0, 100, log(5), 1), each = 3), nrow = 3)
  doses <- c(1, 5, 25)
  result <- mod$fct(doses, parm)
  expect_length(result, 3)
  expect_equal(result[2], 50, tolerance = 1e-10)
})

# ==============================================================================
# Test: deriv1 (first derivatives in parameters)
# ==============================================================================

test_that("llogistic2 deriv1 returns correct dimensions", {
  mod <- llogistic2()
  parm <- matrix(c(1, 0, 100, log(5), 1), nrow = 1)
  result <- mod$deriv1(5, parm)
  expect_length(result, 5)
})

test_that("llogistic2 deriv1 works with fixed parameters", {
  mod <- llogistic2(fixed = c(NA, 0, NA, NA, 1))
  parm <- matrix(c(1, 100, log(5)), nrow = 1)
  result <- mod$deriv1(5, parm)
  expect_length(result, 3)
})

test_that("llogistic2 deriv1 handles multiple rows", {
  mod <- llogistic2()
  parm <- matrix(c(1, 0, 100, log(5), 1,
                    2, 10, 90, log(3), 1), nrow = 2, byrow = TRUE)
  doses <- c(5, 3)
  result <- mod$deriv1(doses, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 5)
})

# ==============================================================================
# Test: derivx (first derivative in dose)
# ==============================================================================

test_that("llogistic2 derivx returns correct structure", {
  mod <- llogistic2()
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, log(5), 1), nrow = 3, ncol = 5, byrow = TRUE)
  result <- mod$derivx(dose, parm)
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("llogistic2 derivx works with fixed parameters", {
  mod <- llogistic2(fixed = c(NA, 0, NA, NA, 1))
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 100, log(5)), nrow = 3, ncol = 3, byrow = TRUE)
  result <- mod$derivx(dose, parm)
  expect_length(result, 3)
})

# ==============================================================================
# Test: edfct (effective dose function)
# ==============================================================================

test_that("llogistic2 edfct works with relative type", {
  mod <- llogistic2()
  parm <- c(1, 0, 100, log(5), 1)
  result <- mod$edfct(parm, 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_length(result, 2)
  expect_true(is.numeric(result[[1]]))
  # ED50: lEDp = e + log(100/(100-50)^(1/f) - 1)/b = log(5) + log(2^1 - 1)/1 = log(5) + 0
  expect_equal(result[[1]], log(5), tolerance = 1e-10)
})

test_that("llogistic2 edfct returns gradient of correct length with fixed params", {
  mod <- llogistic2(fixed = c(NA, 0, NA, NA, 1))
  parm <- c(1, 100, log(5))
  result <- mod$edfct(parm, 50, reference = "control", type = "relative")
  expect_length(result[[2]], 3)
})

# ==============================================================================
# Test: invfct (inverse function)
# ==============================================================================

test_that("llogistic2 invfct returns correct value", {
  mod <- llogistic2()
  # b=1, c=0, d=100, e=log(5), f=1
  # invfct(y) = exp(log(((d-c)/(y-c))^(1/f) - 1)/b + e)
  # invfct(50) = exp(log(((100-0)/(50-0))^1 - 1)/1 + log(5)) = exp(0 + log(5)) = 5
  result <- mod$inversion(50, c(1, 0, 100, log(5), 1))
  expect_equal(result, 5, tolerance = 1e-10)
})

test_that("llogistic2 invfct is consistent with fct", {
  mod <- llogistic2()
  parms <- c(1, 0, 100, log(5), 1)
  dose <- 3
  parm_matrix <- matrix(parms, nrow = 1)
  y <- mod$fct(dose, parm_matrix)
  x_recovered <- mod$inversion(y, parms)
  expect_equal(x_recovered, dose, tolerance = 1e-10)
})

test_that("llogistic2 invfct handles fixed parameters", {
  mod <- llogistic2(fixed = c(NA, 0, NA, NA, 1))
  result <- mod$inversion(50, c(1, 100, log(5)))
  expect_equal(result, 5, tolerance = 1e-10)
})

# ==============================================================================
# Test: Self-starter ss="1" (version 1, default)
# ==============================================================================

test_that("llogistic2 ss=1 ssfct computes start values for normal data", {
  mod <- llogistic2(ss = "1")
  # Create typical dose-response data
  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10, 50, 100),
    resp = c(95, 85, 75, 55, 25, 10, 3, 1)
  )
  result <- mod$ssfct(dframe)
  expect_length(result, 5)
  expect_true(all(is.finite(result)))
})

test_that("llogistic2 ss=1 ssfct returns early for single unique dose", {
  mod <- llogistic2(ss = "1")
  # All same dose
  dframe <- data.frame(dose = rep(5, 5), resp = c(45, 50, 55, 48, 52))
  result <- mod$ssfct(dframe)
  # Should return c(NA, NA, max(y)+0.001, NA, NA)[notFixed] which is all 5 values
  # since all 5 params are free
  expect_length(result, 5)
  # Only d should be non-NA
  expect_true(is.na(result[1]))  # b is NA
  expect_true(is.na(result[2]))  # c is NA
  expect_true(!is.na(result[3]))  # d is max(y) + 0.001
  expect_true(is.na(result[4]))  # e is NA
  expect_true(is.na(result[5]))  # f is NA
})

# ==============================================================================
# Test: Self-starter ss="2" (version 2)
# ==============================================================================

test_that("llogistic2 ss=2 ssfct computes start values for normal data", {
  mod <- llogistic2(ss = "2")
  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10, 50, 100),
    resp = c(95, 85, 75, 55, 25, 10, 3, 1)
  )
  result <- mod$ssfct(dframe)
  expect_length(result, 5)
  expect_true(all(is.finite(result)))
})

test_that("llogistic2 ss=2 ssfct returns early for single unique dose", {
  mod <- llogistic2(ss = "2")
  dframe <- data.frame(dose = rep(5, 5), resp = c(45, 50, 55, 48, 52))
  result <- mod$ssfct(dframe)
  expect_length(result, 5)
  expect_true(is.na(result[1]))   # b
  expect_true(is.na(result[2]))   # c
  expect_true(!is.na(result[3]))  # d
  expect_true(is.na(result[4]))   # e
  expect_true(is.na(result[5]))   # f
})

test_that("llogistic2 ss=2 ssfct uses fixed c and d when specified", {
  mod <- llogistic2(ss = "2", fixed = c(NA, 0, 100, NA, NA))
  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10, 50, 100),
    resp = c(95, 85, 75, 55, 25, 10, 3, 1)
  )
  result <- mod$ssfct(dframe)
  # With c=0 and d=100 fixed, only 3 free params: b, e, f
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# Test: Self-starter ss="3" (version 3)
# ==============================================================================

test_that("llogistic2 ss=3 ssfct computes start values for normal data", {
  mod <- llogistic2(ss = "3")
  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 10, 20, 50, 100),
    resp = c(95, 85, 75, 55, 25, 10, 3, 1)
  )
  result <- mod$ssfct(dframe)
  expect_length(result, 5)
  expect_true(all(is.finite(result)))
})

test_that("llogistic2 ss=3 ssfct returns early for single unique dose", {
  mod <- llogistic2(ss = "3")
  dframe <- data.frame(dose = rep(5, 5), resp = c(45, 50, 55, 48, 52))
  result <- mod$ssfct(dframe)
  expect_length(result, 5)
  expect_true(is.na(result[1]))
  expect_true(is.na(result[2]))
  expect_true(!is.na(result[3]))
  expect_true(is.na(result[4]))
  expect_true(is.na(result[5]))
})

test_that("llogistic2 ss=3 ssfct uses fixed c and d when specified", {
  mod <- llogistic2(ss = "3", fixed = c(NA, 0, 100, NA, NA))
  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 10, 20, 50, 100),
    resp = c(95, 85, 75, 55, 25, 10, 3, 1)
  )
  result <- mod$ssfct(dframe)
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# Test: lowerAs, upperAs, monoton helper functions
# ==============================================================================

test_that("llogistic2 lowerAs and upperAs functions work", {
  mod <- llogistic2()
  parms <- c(1, 0, 100, log(5), 1)
  expect_equal(mod$lowerAs(parms), 0)
  expect_equal(mod$upperAs(parms), 100)
})

test_that("llogistic2 monoton function works", {
  mod <- llogistic2()
  parms <- c(1, 0, 100, log(5), 1)
  result <- mod$monoton(parms)
  # monoParm with signVal=-1 and parmNo=1: -1 * parmVec[1] = -1 * 1 = -1
  expect_equal(result, -1)
})

# ==============================================================================
# Test: LL2.2() wrapper
# ==============================================================================

test_that("LL2.2 returns correct class and structure", {
  result <- LL2.2()
  expect_s3_class(result, "llogistic")
  expect_equal(result$noParm, 2)
  expect_equal(result$names, c("b", "e"))
})

test_that("LL2.2 with custom upper limit", {
  result <- LL2.2(upper = 100)
  expect_true(grepl("upper limit at 100", result$text))
})

test_that("LL2.2 errors on invalid names", {
  expect_error(LL2.2(names = c("a")), "Not correct 'names' argument")
  expect_error(LL2.2(names = 42), "Not correct 'names' argument")
})

test_that("LL2.2 errors on invalid fixed", {
  expect_error(LL2.2(fixed = c(NA)), "Not correct length of 'fixed' argument")
  expect_error(LL2.2(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("LL2.2 fct computes correctly", {
  mod <- LL2.2(upper = 1)
  dose <- c(0.5, 1, 2)
  # free: b, e; fixed: c=0, d=1, f=1
  parm <- matrix(c(1, log(1)), nrow = 3, ncol = 2, byrow = TRUE)
  result <- mod$fct(dose, parm)
  # f(x) = 0 + (1-0)/((1+exp(1*(log(x)-log(1))))^1) = 1/(1+x)
  expected <- 1 / (1 + dose)
  expect_equal(result, expected, tolerance = 1e-10)
})

# ==============================================================================
# Test: LL2.3() wrapper
# ==============================================================================

test_that("LL2.3 returns correct class and structure", {
  result <- LL2.3()
  expect_s3_class(result, "llogistic")
  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("b", "d", "e"))
})

test_that("LL2.3 text indicates lower limit fixed at 0", {
  result <- LL2.3()
  expect_true(grepl("lower limit at 0", result$text))
})

test_that("LL2.3 errors on invalid names", {
  expect_error(LL2.3(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(LL2.3(names = 123), "Not correct 'names' argument")
})

test_that("LL2.3 errors on invalid fixed", {
  expect_error(LL2.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

# ==============================================================================
# Test: LL2.3u() wrapper
# ==============================================================================

test_that("LL2.3u returns correct class and structure", {
  result <- LL2.3u()
  expect_s3_class(result, "llogistic")
  expect_equal(result$noParm, 3)
  expect_equal(result$names, c("b", "c", "e"))
})

test_that("LL2.3u text indicates upper limit fixed", {
  result <- LL2.3u(upper = 1)
  expect_true(grepl("upper limit at 1", result$text))
})

test_that("LL2.3u errors on invalid names", {
  expect_error(LL2.3u(names = c("x")), "Not correct 'names' argument")
  expect_error(LL2.3u(names = 99), "Not correct 'names' argument")
})

test_that("LL2.3u errors on invalid fixed", {
  expect_error(LL2.3u(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

# ==============================================================================
# Test: LL2.4() wrapper
# ==============================================================================

test_that("LL2.4 returns correct class and structure", {
  result <- LL2.4()
  expect_s3_class(result, "llogistic")
  expect_equal(result$noParm, 4)
  expect_equal(result$names, c("b", "c", "d", "e"))
})

test_that("LL2.4 errors on invalid names", {
  expect_error(LL2.4(names = c("a", "b")), "Not correct names argument")
  expect_error(LL2.4(names = 123), "Not correct names argument")
})

test_that("LL2.4 errors on invalid fixed", {
  expect_error(LL2.4(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
})

# ==============================================================================
# Test: LL2.5() wrapper
# ==============================================================================

test_that("LL2.5 returns correct class and structure", {
  result <- LL2.5()
  expect_s3_class(result, "llogistic")
  expect_equal(result$noParm, 5)
  expect_equal(result$names, c("b", "c", "d", "e", "f"))
  expect_equal(result$text, "Generalised log-logistic (log(ED50) as parameter)")
})

test_that("LL2.5 passes additional arguments to llogistic2", {
  custom_ss <- function(dframe) { c(1, 0, 100, 5, 1) }
  result <- LL2.5(ssfct = custom_ss)
  expect_identical(result$ssfct, custom_ss)
})

# ==============================================================================
# Test: lowupFixed, lowFixed, upFixed helpers
# ==============================================================================

test_that("lowupFixed returns correct string", {
  result <- lowupFixed("Model A", 100)
  expect_equal(result, "Model A with lower limit at 0 and upper limit at 100")
})

test_that("lowFixed returns correct string", {
  result <- lowFixed("Model A")
  expect_equal(result, "Model A with lower limit at 0")
})

test_that("upFixed returns correct string", {
  result <- upFixed("Model A", 50)
  expect_equal(result, "Model A with upper limit at 50")
})

# ==============================================================================
# Integration test: drm model fitting
# ==============================================================================

test_that("LL2.4 works in drm model fitting", {
  ryegrass <- data.frame(
    rootl = c(
      7.58, 8.00, 8.33, 7.25, 7.17, 7.00, 7.17, 7.83, 7.92, 7.58,
      6.17, 5.75, 5.83, 6.00, 5.83, 4.92, 4.50, 4.17, 4.42, 4.00,
      2.67, 2.08, 2.42, 2.50, 2.25, 1.17, 0.75, 0.92, 1.00, 0.58
    ),
    conc = c(
      rep(0, 5), rep(0.94, 5), rep(1.88, 5),
      rep(3.75, 5), rep(7.50, 5), rep(15, 5)
    )
  )
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.4())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 4)
  preds <- predict(m1)
  expect_true(all(is.finite(preds)))
})

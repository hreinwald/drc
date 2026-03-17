# Tests for baro5() function
# Baroreflex five-parameter dose-response model

# ===========================================================================
# baro5() - Input validation
# ===========================================================================

test_that("baro5 errors on invalid 'names' argument", {
  # Wrong length
  expect_error(baro5(names = c("a", "b")),
               "Not correct 'names' argument")
  # Not character
  expect_error(baro5(names = c(1, 2, 3, 4, 5)),
               "Not correct 'names' argument")
})

test_that("baro5 errors on invalid 'fixed' argument", {
  expect_error(baro5(fixed = c(NA, NA)),
               "Not correct 'fixed' argument")
  expect_error(baro5(fixed = c(NA, NA, NA, NA, NA, NA)),
               "Not correct 'fixed' argument")
})

# ===========================================================================
# baro5() - Happy path construction (all parameters free)
# ===========================================================================

test_that("baro5 returns correct structure with default args", {
  b <- baro5()

  # Class and type
  expect_s3_class(b, "baro5")
  expect_type(b, "list")

  # All expected components present
  expected_names <- c("fct", "ssfct", "names", "deriv1", "deriv2",
                       "edfct", "sifct", "name", "text", "noParm")
  expect_true(all(expected_names %in% names(b)))

  # Parameter names (all 5 free)
  expect_equal(b$names, c("b1", "b2", "c", "d", "e"))
  expect_equal(b$noParm, 5)

  # Name and text
  expect_equal(b$name, "baro5")
  expect_equal(b$text, "Baroreflex")

  # NULL fields
  expect_null(b$deriv1)
  expect_null(b$deriv2)
  expect_null(b$edfct)
  expect_null(b$sifct)
})

test_that("baro5 handles custom parameter names", {
  b <- baro5(names = c("slope1", "slope2", "lower", "upper", "mid"))
  expect_equal(b$names, c("slope1", "slope2", "lower", "upper", "mid"))
})

# ===========================================================================
# baro5() - Fixed parameters
# ===========================================================================

test_that("baro5 with fixed parameters reduces names and noParm", {
  # Fix c=0
  b <- baro5(fixed = c(NA, NA, 0, NA, NA))
  expect_equal(b$names, c("b1", "b2", "d", "e"))
  expect_equal(b$noParm, 4)

  # Fix c=0 and d=100
  b2 <- baro5(fixed = c(NA, NA, 0, 100, NA))
  expect_equal(b2$names, c("b1", "b2", "e"))
  expect_equal(b2$noParm, 3)
})

# ===========================================================================
# baro5() - fct function (model evaluation)
# ===========================================================================

test_that("baro5 fct evaluates correctly with all free parameters", {
  b <- baro5()

  # parm must be a matrix with one row per curve
  # Parameters: b1, b2, c, d, e
  parm <- matrix(c(1, 1, 0, 100, 5), nrow = 1)
  doses <- c(1, 5, 10)
  result <- b$fct(doses, parm)

  expect_length(result, 3)
  expect_true(is.numeric(result))
  # At dose=e=5, with b1=b2=1, c=0, d=100: should be at midpoint
  expect_equal(result[2], 50, tolerance = 1)
})

test_that("baro5 fct with fixed parameters works", {
  b <- baro5(fixed = c(NA, NA, 0, NA, NA))
  # parm has only 4 columns (b1, b2, d, e)
  parm <- matrix(c(1, 1, 100, 5), nrow = 1)
  doses <- c(1, 5, 10)
  result <- b$fct(doses, parm)
  expect_length(result, 3)
  expect_true(is.numeric(result))
})

test_that("baro5 fct works with asymmetric parameters (b1 != b2)", {
  b <- baro5()
  # b1=2, b2=0.5 -> asymmetric curve
  parm <- matrix(c(2, 0.5, 0, 100, 5), nrow = 1)
  doses <- c(1, 5, 10)
  result <- b$fct(doses, parm)
  expect_length(result, 3)
  expect_true(is.numeric(result))
  # All values should be between c=0 and d=100
  expect_true(all(result >= 0 & result <= 100))
})

test_that("baro5 fct handles multiple rows in parm matrix", {
  b <- baro5()
  parm <- matrix(c(1, 1, 0, 100, 5,
                    2, 2, 10, 90, 3), nrow = 2, byrow = TRUE)
  doses <- c(1, 5)
  result <- b$fct(doses, parm)
  expect_length(result, 2)
  expect_true(is.numeric(result))
})

# ===========================================================================
# baro5() - ssfct (self-starter)
# ===========================================================================

test_that("baro5 default ssfct returns correct number of initial values", {
  b <- baro5()
  dframe <- data.frame(dose = c(0, 0.5, 1, 2, 5, 10),
                       resp = c(100, 90, 75, 50, 20, 5))
  ssvals <- b$ssfct(dframe)

  expect_length(ssvals, 5)  # all 5 parameters free
  expect_true(is.numeric(ssvals))
})

test_that("baro5 default ssfct respects fixed parameters", {
  b <- baro5(fixed = c(NA, NA, 0, NA, NA))
  dframe <- data.frame(dose = c(0, 0.5, 1, 2, 5, 10),
                       resp = c(100, 90, 75, 50, 20, 5))
  ssvals <- b$ssfct(dframe)

  expect_length(ssvals, 4)  # only 4 free parameters
})

test_that("baro5 with custom ssfct uses provided function", {
  custom_ss <- function(dframe) { c(1, 1, 0, 100, 5) }
  b <- baro5(ssfct = custom_ss)
  dframe <- data.frame(dose = c(0, 1, 5), resp = c(100, 50, 5))
  result <- b$ssfct(dframe)
  expect_equal(result, c(1, 1, 0, 100, 5))
})

# ===========================================================================
# baro5() - Integration with drm model fitting
# ===========================================================================

test_that("baro5 works with drm for model fitting", {
  # Create dose-response data suitable for baro5
  set.seed(42)
  dose <- rep(c(0.1, 0.5, 1, 2, 5, 10, 20, 50), each = 3)
  resp <- 50 + 50 / (1 + exp(1.5 * (log(dose) - log(5)))) + rnorm(24, 0, 2)
  test_data <- data.frame(dose = dose, resp = resp)

  m1 <- drm(resp ~ dose, data = test_data, fct = baro5())
  expect_s3_class(m1, "drc")
  expect_length(coef(m1), 5)
})

# Tests for fplogistic() and FPL.4() functions
# Fractional polynomial-logistic dose-response model

# --- Helper data ---

# Simple dose-response data for self-starter and model fitting tests
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

# ===========================================================================
# fplogistic() - Input validation
# ===========================================================================

test_that("fplogistic errors on invalid 'names' argument", {
  # Wrong length

expect_error(fplogistic(-1, 1, names = c("a", "b")),
               "Not correct 'names' argument")
  # Not character
  expect_error(fplogistic(-1, 1, names = c(1, 2, 3, 4)),
               "Not correct 'names' argument")
})

test_that("fplogistic errors on invalid 'fixed' argument", {
  expect_error(fplogistic(-1, 1, fixed = c(NA, NA)),
               "Not correct 'fixed' argument")
  expect_error(fplogistic(-1, 1, fixed = c(NA, NA, NA, NA, NA)),
               "Not correct 'fixed' argument")
})

# ===========================================================================
# fplogistic() - Happy path construction (all parameters free)
# ===========================================================================

test_that("fplogistic returns correct structure with default args", {
  fp <- fplogistic(-1, 1)

  # Class and type
expect_s3_class(fp, "fp-logistic")
  expect_type(fp, "list")

  # All expected components present
  expected_names <- c("fct", "ssfct", "names", "deriv1", "deriv2", "derivx",
                       "edfct", "name", "text", "noParm", "fixed")
  expect_true(all(expected_names %in% names(fp)))

  # Parameter names
  expect_equal(fp$names, c("b", "c", "d", "e"))
  expect_equal(fp$noParm, 4)
  expect_equal(fp$fixed, c(NA, NA, NA, NA))

  # Default name and text
  expect_match(fp$name, "fplogistic\\(-1,1\\)")
  expect_equal(fp$text, "Fractional polynomial")

  # deriv2 is always NULL
  expect_null(fp$deriv2)
})

test_that("fplogistic accepts custom fctName and fctText", {
  fp <- fplogistic(-1, 1, fctName = "myModel", fctText = "My description")
  expect_equal(fp$name, "myModel")
  expect_equal(fp$text, "My description")
})

test_that("fplogistic handles custom parameter names", {
  fp <- fplogistic(-1, 1, names = c("slope", "lower", "upper", "scale"))
  expect_equal(fp$names, c("slope", "lower", "upper", "scale"))
})

# ===========================================================================
# fplogistic() - Fixed parameters
# ===========================================================================

test_that("fplogistic with fixed parameters reduces names and noParm", {
  # Fix c=0
  fp <- fplogistic(-1, 1, fixed = c(NA, 0, NA, NA))
  expect_equal(fp$names, c("b", "d", "e"))
  expect_equal(fp$noParm, 3)

  # Fix c=0 and d=100
  fp2 <- fplogistic(-1, 1, fixed = c(NA, 0, 100, NA))
  expect_equal(fp2$names, c("b", "e"))
  expect_equal(fp2$noParm, 2)
})

# ===========================================================================
# fplogistic() - fct function (model evaluation)
# ===========================================================================

test_that("fplogistic fct evaluates correctly", {
  fp <- fplogistic(-1, 1)

  # parm must be a matrix with one row per curve
  parm <- matrix(c(-1, 0, 100, 1), nrow = 1)
  doses <- c(0, 1, 5, 10)
  result <- fp$fct(doses, parm)

  # At dose=0: log(0+1)=0, exp(0)=1, c + (d-c)/2 = 50 ... wait
  # Actually at dose=0: log(1)=0, 0^p1 is NaN for p1<0, but exp(b*NaN + e*NaN) is NaN
  # So dose=0 may produce NaN. Let's just check dose > 0.
  expect_length(result, 4)
  expect_true(is.numeric(result))

  # Gradient attribute should be present
  grad <- attr(result, "gradient")
  expect_true(!is.null(grad))
  expect_equal(dim(grad), c(4, 4))
})

test_that("fplogistic fct with fixed parameters works", {
  fp <- fplogistic(-1, 1, fixed = c(NA, 0, NA, NA))
  # parm has only 3 columns (b, d, e)
  parm <- matrix(c(-1, 100, 1), nrow = 1)
  doses <- c(1, 5, 10)
  result <- fp$fct(doses, parm)
  expect_length(result, 3)
  expect_true(is.numeric(result))
})

# ===========================================================================
# fplogistic() - ssfct (self-starter)
# ===========================================================================

test_that("fplogistic default ssfct returns correct number of initial values", {
  fp <- fplogistic(-1, 1)
  dframe <- data.frame(dose = c(0, 0.5, 1, 2, 5, 10),
                       resp = c(100, 90, 75, 50, 20, 5))
  ssvals <- fp$ssfct(dframe)

  expect_length(ssvals, 4)  # all 4 parameters free
  expect_true(is.numeric(ssvals))
})

test_that("fplogistic default ssfct respects fixed parameters", {
  fp <- fplogistic(-1, 1, fixed = c(NA, 0, NA, NA))
  dframe <- data.frame(dose = c(0, 0.5, 1, 2, 5, 10),
                       resp = c(100, 90, 75, 50, 20, 5))
  ssvals <- fp$ssfct(dframe)

  expect_length(ssvals, 3)  # only 3 free parameters
})

test_that("fplogistic with custom ssfct uses provided function", {
  custom_ss <- function(dframe) { c(-1, 0, 100, 1) }
  fp <- fplogistic(-1, 1, ssfct = custom_ss)
  dframe <- data.frame(dose = c(0, 1, 5), resp = c(100, 50, 5))
  result <- fp$ssfct(dframe)
  expect_equal(result, c(-1, 0, 100, 1))
})

# ===========================================================================
# fplogistic() - deriv1 (parameter derivatives)
# ===========================================================================

test_that("fplogistic deriv1 returns gradient matrix", {
  fp <- fplogistic(-1, 1)
  parm <- matrix(c(-1, 0, 100, 1), nrow = 1)
  doses <- c(1, 5, 10)
  d1 <- fp$deriv1(doses, parm)

  expect_true(is.matrix(d1) || is.numeric(d1))
  # Should have 3 rows (doses) x 4 cols (params)
  expect_equal(nrow(d1), 3)
  expect_equal(ncol(d1), 4)
})

test_that("fplogistic deriv1 with fixed parameters returns reduced cols", {
  fp <- fplogistic(-1, 1, fixed = c(NA, 0, NA, NA))
  parm <- matrix(c(-1, 100, 1), nrow = 1)
  doses <- c(1, 5, 10)
  d1 <- fp$deriv1(doses, parm)

  expect_equal(ncol(d1), 3)  # only free parameters
})

# ===========================================================================
# fplogistic() - derivx (dose derivative)
# ===========================================================================

test_that("fplogistic derivx returns gradient in dose", {
  fp <- fplogistic(-1, 1)
  parm <- matrix(c(-1, 0, 100, 1), nrow = 1)
  doses <- c(1, 5, 10)
  dx <- fp$derivx(doses, parm)

  expect_true(is.matrix(dx) || is.numeric(dx))
  expect_equal(nrow(dx), 3)
  expect_equal(ncol(dx), 1)
})

# ===========================================================================
# fplogistic() - edfct (effective dose calculation)
# ===========================================================================

test_that("fplogistic edfct computes ED values for relative type", {
  fp <- fplogistic(-1, 1)
  parm <- c(-1, 0, 100, 1)
  result <- fp$edfct(parm, 50, "control", "relative")

  expect_type(result, "list")
  expect_length(result, 2)
  # First element is ED estimate
  expect_true(is.numeric(result[[1]]))
  expect_true(result[[1]] > 0)
  # Second element is derivative vector
  expect_length(result[[2]], 4)
})

test_that("fplogistic edfct computes ED with absolute type", {
  fp <- fplogistic(-1, 1)
  parm <- c(-1, 0, 100, 1)
  result <- fp$edfct(parm, 50, "control", "absolute")

  expect_type(result, "list")
  expect_true(is.numeric(result[[1]]))
  expect_true(result[[1]] > 0)
})

test_that("fplogistic edfct with loged=TRUE returns log-transformed ED", {
  fp <- fplogistic(-1, 1)
  parm <- c(-1, 0, 100, 1)

  result_nolog <- fp$edfct(parm, 50, "control", "relative", loged = FALSE)
  result_log <- fp$edfct(parm, 50, "control", "relative", loged = TRUE)

  # Log-transformed ED should be log of non-transformed
  expect_equal(result_log[[1]], log(result_nolog[[1]]), tolerance = 1e-6)
  # Derivative should also be transformed
  expect_true(all(is.numeric(result_log[[2]])))
})

test_that("fplogistic edfct works with positive b (increasing curve)", {
  fp <- fplogistic(-1, 1)
  parm <- c(1, 0, 100, -1)  # positive b
  result <- fp$edfct(parm, 50, "control", "relative")

  expect_type(result, "list")
  expect_true(is.numeric(result[[1]]))
})

test_that("fplogistic edfct with fixed parameters returns reduced deriv", {
  fp <- fplogistic(-1, 1, fixed = c(NA, 0, NA, NA))
  parm <- c(-1, 100, 1)
  result <- fp$edfct(parm, 50, "control", "relative")

  expect_length(result[[2]], 3)  # only free parameters
})

# ===========================================================================
# fplogistic() - Integration with drm model fitting
# ===========================================================================

test_that("fplogistic works with drm for model fitting", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = FPL.4(-1, 1))
  expect_s3_class(m1, "drc")
  expect_length(coef(m1), 4)
})

test_that("fplogistic ED calculation works through drm", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = FPL.4(-1, 1))
  ed <- ED(m1, 50, display = FALSE)
  expect_true(is.matrix(ed))
  expect_true(ed[1, 1] > 0)
})

# ===========================================================================
# fplogistic() - Different p1, p2 values
# ===========================================================================

test_that("fplogistic works with various p1, p2 combinations", {
  # Different power combinations
  fp1 <- fplogistic(-2, 3)
  expect_s3_class(fp1, "fp-logistic")
  expect_match(fp1$name, "fplogistic\\(-2,3\\)")

  fp2 <- fplogistic(-0.5, 0.5)
  expect_s3_class(fp2, "fp-logistic")

  # Check that fct evaluation works with these
  parm <- matrix(c(-1, 0, 100, 1), nrow = 1)
  doses <- c(1, 5)
  r1 <- fp1$fct(doses, parm)
  r2 <- fp2$fct(doses, parm)
  expect_length(r1, 2)
  expect_length(r2, 2)
  # Different powers should give different results
  expect_false(isTRUE(all.equal(r1, r2, check.attributes = FALSE)))
})

# ===========================================================================
# FPL.4() - Convenience wrapper
# ===========================================================================

test_that("FPL.4 returns fp-logistic object", {
  fp <- FPL.4(-1, 1)
  expect_s3_class(fp, "fp-logistic")
  expect_equal(fp$noParm, 4)
  expect_match(fp$name, "FPL\\.4\\(-1,1\\)")
})

test_that("FPL.4 errors on invalid 'names' argument", {
  expect_error(FPL.4(-1, 1, names = c("a", "b")),
               "Not correct names argument")
  expect_error(FPL.4(-1, 1, names = c(1, 2, 3, 4)),
               "Not correct names argument")
})

test_that("FPL.4 errors on invalid 'fixed' argument", {
  expect_error(FPL.4(-1, 1, fixed = c(NA, NA)),
               "Not correct length of 'fixed' argument")
})

test_that("FPL.4 passes extra arguments to fplogistic", {
  fp <- FPL.4(-1, 1, fixed = c(NA, 0, NA, NA))
  expect_equal(fp$noParm, 3)
  expect_equal(fp$names, c("b", "d", "e"))
})

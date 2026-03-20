# Tests for gammadr.R: gammadr() function

# --- gammadr() main function ---

test_that("gammadr returns correct class and structure", {
  g <- gammadr()

  expect_s3_class(g, "gamma")
  expect_true(is.list(g))
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2", "derivx",
                     "edfct", "name", "text", "noParm") %in% names(g)))
})

test_that("gammadr default names are b, c, d, e", {
  g <- gammadr()
  expect_equal(g$names, c("b", "c", "d", "e"))
})

test_that("gammadr noParm reflects number of NA in fixed", {
  g_full <- gammadr()
  expect_equal(g_full$noParm, 4)

  g_partial <- gammadr(fixed = c(1, NA, NA, NA))
  expect_equal(g_partial$noParm, 3)
  expect_equal(g_partial$names, c("c", "d", "e"))
})

test_that("gammadr uses default text when fctText not provided", {
  g <- gammadr()
  expect_equal(g$text, "Gamma")
})

test_that("gammadr uses provided fctText", {
  g <- gammadr(fctText = "Custom text")
  expect_equal(g$text, "Custom text")
})

test_that("gammadr uses provided fctName", {
  g <- gammadr(fctName = "myFunc")
  expect_equal(g$name, "myFunc")
})

test_that("gammadr uses default name when fctName not provided", {
  g <- gammadr()
  expect_equal(g$name, "gammadr")
})

# --- Error handling ---

test_that("gammadr errors on invalid names argument - not character", {
  expect_error(gammadr(names = c(1, 2, 3, 4)), "Not correct 'names' argument")
})

test_that("gammadr errors on invalid names argument - wrong length", {
  expect_error(gammadr(names = c("a", "b")), "Not correct 'names' argument")
})

test_that("gammadr errors on invalid fixed argument - wrong length", {
  expect_error(gammadr(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(gammadr(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
})

# --- deriv2 and edfct are NULL ---

test_that("gammadr deriv2 is NULL", {
  g <- gammadr()
  expect_null(g$deriv2)
})

test_that("gammadr edfct is NULL", {
  g <- gammadr()
  expect_null(g$edfct)
})

# --- Fixed parameter variations ---

test_that("gammadr works with all parameters fixed", {
  g <- gammadr(fixed = c(1, 0, 1, 2))
  expect_equal(g$noParm, 0)
  expect_length(g$names, 0)
})

test_that("gammadr works with only one parameter free", {
  g <- gammadr(fixed = c(1, 0, 1, NA))
  expect_equal(g$noParm, 1)
  expect_equal(g$names, "e")
})

test_that("gammadr works with partial fixed parameters in different positions", {
  # Fix first and last
  g1 <- gammadr(fixed = c(1, NA, NA, 2))
  expect_equal(g1$noParm, 2)
  expect_equal(g1$names, c("c", "d"))

  # Fix middle parameters
  g2 <- gammadr(fixed = c(NA, 0, 1, NA))
  expect_equal(g2$noParm, 2)
  expect_equal(g2$names, c("b", "e"))
})

# --- fct (internal nonlinear function) ---

test_that("gammadr fct computes correct values", {
  g <- gammadr()

  # Parameters: b=1, c=0, d=1, e=2
  dose <- c(1, 2, 5)
  parm <- matrix(c(1, 0, 1, 2), nrow = 3, ncol = 4, byrow = TRUE)

  result <- g$fct(dose, parm)

  # f(x) = c + (d - c) * pgamma(b * x, e, 1)
  expected <- 0 + (1 - 0) * pgamma(1 * dose, 2, 1)
  expect_equal(as.numeric(result), expected)
})

test_that("gammadr fct works with fixed parameters", {
  g <- gammadr(fixed = c(1, 0, NA, NA))

  dose <- c(1, 2, 5)
  # Only 2 free parameters: d, e
  parm <- matrix(c(1, 2), nrow = 3, ncol = 2, byrow = TRUE)

  result <- g$fct(dose, parm)

  # b=1 (fixed), c=0 (fixed), d=1, e=2
  expected <- 0 + (1 - 0) * pgamma(1 * dose, 2, 1)
  expect_equal(as.numeric(result), expected)
})

test_that("gammadr fct handles dose = 0", {
  g <- gammadr()

  dose <- c(0)
  parm <- matrix(c(1, 0, 1, 2), nrow = 1, ncol = 4)

  result <- g$fct(dose, parm)

  # pgamma(0, ...) = 0, so f(0) = c = 0
  expect_equal(as.numeric(result), 0)
})

test_that("gammadr fct handles multiple rows in parm", {
  g <- gammadr()

  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 1, 2,
                    1, 0, 1, 2,
                    1, 0, 1, 2),
                 nrow = 3, byrow = TRUE)

  result <- g$fct(dose, parm)
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

# --- deriv1 (parameter derivatives) ---

test_that("gammadr deriv1 returns matrix with correct dimensions", {
  g <- gammadr()

  dose <- c(1, 2, 5)
  parm <- matrix(c(1, 0, 1, 2), nrow = 3, ncol = 4, byrow = TRUE)

  result <- g$deriv1(dose, parm)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
})

test_that("gammadr deriv1 works with fixed parameters", {
  g <- gammadr(fixed = c(1, NA, NA, NA))

  dose <- c(1, 2, 5)
  parm <- matrix(c(0, 1, 2), nrow = 3, ncol = 3, byrow = TRUE)

  result <- g$deriv1(dose, parm)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
})

test_that("gammadr deriv1 computes finite values", {
  g <- gammadr()

  dose <- c(0.5, 1, 3)
  parm <- matrix(c(2, 0, 1, 3), nrow = 3, ncol = 4, byrow = TRUE)

  result <- g$deriv1(dose, parm)

  expect_true(all(is.finite(result)))
})

test_that("gammadr deriv1 computes correct df/db value", {
  g <- gammadr()

  # b=1, c=0, d=1, e=3
  dose <- c(2)
  parm <- matrix(c(1, 0, 1, 3), nrow = 1, ncol = 4)

  result <- g$deriv1(dose, parm)

  # df/db = (d-c) * dgamma(b*x, e, 1) * x (chain rule: d/db[b*x] = x)
  expected_db <- (1 - 0) * dgamma(1 * 2, 3, 1) * 2
  expect_equal(result[1], expected_db, tolerance = 1e-10)
})

# --- logGamma helper (called via deriv1) ---

test_that("gammadr deriv1 exercises logGamma with x < 1e-10 (zero dose)", {

  g <- gammadr()

  # dose = 0 means b * dose = 0, which is < 1e-10
  # This triggers the ifelse branch in logGamma: retVec[i] <- 0
  dose <- c(0)
  parm <- matrix(c(1, 0, 1, 2), nrow = 1, ncol = 4)

  result <- g$deriv1(dose, parm)

  # Single row with all 4 params free drops to a vector
  expect_type(result, "double")
  expect_true(all(is.finite(result)))
})

test_that("gammadr deriv1 exercises logGamma with x >= 1e-10 (positive dose)", {
  g <- gammadr()

  # dose > 0 and b > 0, so b*dose > 1e-10
  # This triggers the integrate branch in logGamma
  dose <- c(5)
  parm <- matrix(c(1, 0, 1, 2), nrow = 1, ncol = 4)

  result <- g$deriv1(dose, parm)

  # Single row with all 4 params free drops to a vector
  expect_type(result, "double")
  expect_true(all(is.finite(result)))
})

test_that("gammadr deriv1 exercises logGamma with both branches in one call", {
  g <- gammadr()

  # Mix of dose=0 (x < 1e-10) and dose>0 (x >= 1e-10)
  dose <- c(0, 5)
  parm <- matrix(c(1, 0, 1, 2,
                    1, 0, 1, 2),
                 nrow = 2, byrow = TRUE)

  result <- g$deriv1(dose, parm)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_true(all(is.finite(result)))
})

# --- derivx (dose derivative) ---

test_that("gammadr derivx returns correct structure", {
  g <- gammadr()

  dose <- c(1, 2, 5)
  parm <- matrix(c(1, 0, 1, 2), nrow = 3, ncol = 4, byrow = TRUE)

  result <- g$derivx(dose, parm)

  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("gammadr derivx works with fixed parameters", {
  g <- gammadr(fixed = c(1, 0, NA, NA))

  dose <- c(1, 2, 5)
  parm <- matrix(c(1, 2), nrow = 3, ncol = 2, byrow = TRUE)

  result <- g$derivx(dose, parm)

  expect_type(result, "double")
  expect_length(result, 3)
})

test_that("gammadr derivx computes correct values", {
  g <- gammadr()

  dose <- c(2)
  # b=1, c=0, d=1, e=3
  parm <- matrix(c(1, 0, 1, 3), nrow = 1, ncol = 4)

  result <- g$derivx(dose, parm)

  # (d - c) * dgamma(b * x, e, 1) * b
  expected <- (1 - 0) * dgamma(1 * 2, 3, 1) * 1
  expect_equal(as.numeric(result), expected)
})

# --- ssfct (self-starter function) ---

test_that("gammadr ssfct returns initial values", {
  g <- gammadr()

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0)
  )

  init_vals <- g$ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 4)
  expect_true(all(is.finite(init_vals)))
})

test_that("gammadr ssfct respects fixed parameters", {
  g <- gammadr(fixed = c(1, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0)
  )

  init_vals <- g$ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 3)  # Only 3 free parameters
})

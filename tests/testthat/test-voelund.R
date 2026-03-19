# ==============================================================================
# Tests for voelund() function from R/voelund.R
# ==============================================================================

# --- Correctness: Default construction ----------------------------------------

test_that("voelund returns correct structure with default parameters", {
  v <- voelund()
  expect_s3_class(v, "Voelund")
  expect_equal(v$name, "voelund")
  expect_equal(v$text, "Voelund mixture")
  expect_equal(v$noParm, 7)
  expect_equal(v$names, c("b", "c", "d", "e", "f", "g", "h"))
  expect_true(is.function(v$fct))
  expect_true(is.function(v$ssfct))
  expect_true(is.function(v$scaleFct))
  expect_null(v$deriv1)
  expect_null(v$deriv2)
  expect_null(v$edfct)
  expect_null(v$sifct)
})

# --- Correctness: Fixed parameters --------------------------------------------

test_that("voelund returns correct structure with fixed parameters", {
  v <- voelund(fixed = c(NA, 0, 1, NA, NA, NA, NA))
  expect_equal(v$noParm, 5)
  expect_equal(v$names, c("b", "e", "f", "g", "h"))
})

# --- Error handling: invalid names argument -----------------------------------

test_that("voelund errors on incorrect names argument", {
  # Wrong length
  expect_error(voelund(names = c("a")), "Not correct 'names' argument")
  # Not character type
  expect_error(voelund(names = c(1, 2, 3, 4, 5, 6, 7)), "Not correct 'names' argument")
})

# --- Error handling: invalid fixed argument -----------------------------------

test_that("voelund errors on incorrect fixed argument", {
  expect_error(voelund(fixed = c(NA, NA)), "Not correct 'fixed' argument")
})

# --- Correctness: fct computes correct values for normal doses ----------------

test_that("voelund fct computes correct values for normal doses", {
  v <- voelund()
  dose <- c(0.1, 1, 10, 100)
  parm <- matrix(c(-2, 0, 1, 10, 10, 1, 1), nrow = 4, ncol = 7, byrow = TRUE)
  result <- v$fct(dose, parm)
  expect_length(result, 4)
  expect_true(all(is.finite(result)))
})

# --- Edge case: zero dose handled by eps threshold ----------------------------

test_that("voelund fct handles zero dose correctly", {
  v <- voelund()
  dose <- c(0, 1, 10)
  parm <- matrix(c(-2, 0, 1, 10, 10, 1, 1), nrow = 3, ncol = 7, byrow = TRUE)
  result <- v$fct(dose, parm)
  # When dose < eps (zero), result should be d parameter (column 3)
  expect_equal(result[1], 1)
})

# --- Edge case: infinite e parameter -----------------------------------------

test_that("voelund fct handles infinite e parameter", {
  v <- voelund()
  dose <- c(0.1, 1, 10)
  parm <- matrix(c(-2, 0, 1, 10, 10, 1, 1), nrow = 3, ncol = 7, byrow = TRUE)
  parm[2, 4] <- Inf  # e parameter = Inf for row 2
  result <- v$fct(dose, parm)
  expect_length(result, 3)
  # When e is Inf, loge should use log(f) instead
  expect_true(is.finite(result[1]))
  expect_true(is.finite(result[3]))
})

# --- Edge case: infinite f parameter -----------------------------------------

test_that("voelund fct handles infinite f parameter", {
  v <- voelund()
  dose <- c(0.1, 1, 10)
  parm <- matrix(c(-2, 0, 1, 10, 10, 1, 1), nrow = 3, ncol = 7, byrow = TRUE)
  parm[2, 5] <- Inf  # f parameter = Inf for row 2
  result <- v$fct(dose, parm)
  expect_length(result, 3)
  # When f is Inf, loge should use log(e) instead
  expect_true(is.finite(result[1]))
  expect_true(is.finite(result[3]))
})

# --- Correctness: default ssfct returns valid starting values -----------------

test_that("voelund default ssfct returns valid starting values", {
  v <- voelund()
  df <- data.frame(dose = c(0.01, 0.1, 1, 10, 100),
                   resp = c(1, 0.95, 0.5, 0.1, 0.01))
  ss <- v$ssfct(df)
  expect_length(ss, 7)
  expect_true(all(is.finite(ss)))
})

# --- Correctness: custom ssfct is used when provided --------------------------

test_that("voelund with custom ssfct", {
  custom_ss <- function(dframe) rep(1, 7)
  v <- voelund(ssfct = custom_ss)
  df <- data.frame(dose = 1:5, resp = 5:1)
  result <- v$ssfct(df)
  expect_equal(result, rep(1, 7))
})

# --- Correctness: scaleFct returns correct scaling ----------------------------

test_that("voelund scaleFct returns correct scaling", {
  v <- voelund()
  sf <- v$scaleFct(10, 100)
  expect_equal(sf, c(1, 100, 100, 10, 10, 1, 1))
})

# --- Correctness: scaleFct respects fixed parameters --------------------------

test_that("voelund scaleFct respects fixed parameters", {
  v <- voelund(fixed = c(NA, 0, NA, NA, NA, NA, NA))
  sf <- v$scaleFct(10, 100)
  expect_equal(sf, c(1, 100, 10, 10, 1, 1))
})

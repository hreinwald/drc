# Tests for weibull2x.R: weibull2x(), W2x.3(), W2x.4()

# --- weibull2x() main function ---

test_that("weibull2x returns correct class and structure", {
  w2x <- weibull2x()
  expect_s3_class(w2x, "Weibull-2")
  expect_true(is.list(w2x))
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2", "derivx",
                     "edfct", "name", "text", "noParm", "fixed") %in% names(w2x)))
})

test_that("weibull2x default names are b, c, d, e, t0", {
  w2x <- weibull2x()
  expect_equal(w2x$names, c("b", "c", "d", "e", "t0"))
})

test_that("weibull2x noParm reflects number of NA in fixed", {
  w2x_full <- weibull2x()
  expect_equal(w2x_full$noParm, 5)

  w2x_partial <- weibull2x(fixed = c(1, NA, NA, NA, NA))
  expect_equal(w2x_partial$noParm, 4)
  expect_equal(w2x_partial$names, c("c", "d", "e", "t0"))
})

test_that("weibull2x uses default text when fctText not provided", {
  w2x <- weibull2x()
  expect_equal(w2x$text, "Weibull (type 2)")
})

test_that("weibull2x uses provided fctText", {
  w2x <- weibull2x(fctText = "Custom text")
  expect_equal(w2x$text, "Custom text")
})

test_that("weibull2x uses provided fctName", {
  w2x <- weibull2x(fctName = "myFunc")
  expect_equal(w2x$name, "myFunc")
})

test_that("weibull2x uses default fctName from match.call", {
  w2x <- weibull2x()
  expect_equal(w2x$name, "weibull2x")
})

test_that("weibull2x derivatives are NULL", {
  w2x <- weibull2x()
  expect_null(w2x$deriv1)
  expect_null(w2x$deriv2)
  expect_null(w2x$derivx)
})

test_that("weibull2x fixed field is preserved", {
  fixed_vals <- c(1, NA, NA, NA, NA)
  w2x <- weibull2x(fixed = fixed_vals)
  expect_equal(w2x$fixed, fixed_vals)
})

# --- Error handling for weibull2x ---

test_that("weibull2x errors on invalid names argument", {
  expect_error(weibull2x(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(weibull2x(names = 123), "Not correct 'names' argument")
})

test_that("weibull2x errors on invalid fixed argument", {
  expect_error(weibull2x(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(weibull2x(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
})

test_that("weibull2x errors when t0 (fixed[5]) is not NA", {
  expect_error(weibull2x(fixed = c(NA, NA, NA, NA, 0)), "The lag time cannot be fixed")
  expect_error(weibull2x(fixed = c(1, 0, 100, 5, 10)), "The lag time cannot be fixed")
})

# --- fct (internal nonlinear function) ---

test_that("weibull2x fct computes correct values for dose > t0", {
  w2x <- weibull2x()
  # Parameters: b=1, c=0, d=100, e=5, t0=1
  dose <- c(3, 6, 11)
  parm <- matrix(c(1, 0, 100, 5, 1), nrow = 3, ncol = 5, byrow = TRUE)
  result <- w2x$fct(dose, parm)
  expected <- 0 + (100 - 0) * (1 - exp(-exp(1 * (log(dose - 1) - log(5)))))
  expect_equal(result, expected)
})

test_that("weibull2x fct returns c when dose <= t0", {
  w2x <- weibull2x()
  # Parameters: b=1, c=5, d=100, e=10, t0=3
  dose <- c(1, 2, 3)  # all <= t0=3
  parm <- matrix(c(1, 5, 100, 10, 3), nrow = 3, ncol = 5, byrow = TRUE)
  result <- w2x$fct(dose, parm)
  expect_equal(result, c(5, 5, 5))
})

test_that("weibull2x fct handles mix of dose > t0 and dose <= t0", {
  w2x <- weibull2x()
  # Parameters: b=1, c=0, d=100, e=5, t0=5
  dose <- c(2, 5, 10)  # first two <= t0, last > t0
  parm <- matrix(c(1, 0, 100, 5, 5), nrow = 3, ncol = 5, byrow = TRUE)
  # NaN warning is expected from log(dose - t0) when dose <= t0
  result <- suppressWarnings(w2x$fct(dose, parm))
  # dose=2 <= t0=5: result = c = 0
  # dose=5 <= t0=5: result = c = 0 (not strictly > t0)
  # dose=10 > t0=5: Weibull formula
  expected_10 <- 0 + (100 - 0) * (1 - exp(-exp(1 * (log(10 - 5) - log(5)))))
  expect_equal(result[1], 0)
  expect_equal(result[2], 0)
  expect_equal(result[3], expected_10)
})

test_that("weibull2x fct works with fixed parameters", {
  w2x <- weibull2x(fixed = c(1, 0, NA, NA, NA))
  dose <- c(6, 11)
  # Only free params: d, e, t0
  parm <- matrix(c(100, 5, 1), nrow = 2, ncol = 3, byrow = TRUE)
  result <- w2x$fct(dose, parm)
  expected <- 0 + (100 - 0) * (1 - exp(-exp(1 * (log(dose - 1) - log(5)))))
  expect_equal(result, expected)
})

# --- ssfct (self-starter function) ---

test_that("weibull2x uses provided ssfct when not NULL", {
  custom_ss <- function(dframe) { c(1, 0, 10, 5, 0) }
  w2x <- weibull2x(ssfct = custom_ss)
  expect_identical(w2x$ssfct, custom_ss)
})

test_that("weibull2x uses default ssfct when ssfct is NULL", {
  w2x <- weibull2x(ssfct = NULL)
  expect_true(is.function(w2x$ssfct))
})

test_that("weibull2x default ssfct returns correct number of values", {
  w2x <- weibull2x()
  # Create a simple test data frame
  dframe <- data.frame(x = c(1, 2, 5, 10, 20), y = c(0, 10, 40, 80, 95))
  result <- w2x$ssfct(dframe)
  # Should return 5 values for the 5 free parameters
  expect_length(result, 5)
  expect_true(all(is.finite(result)))
})

test_that("weibull2x method argument works for self-starter", {
  for (m in c("1", "2", "3", "4")) {
    w2x <- weibull2x(method = m)
    expect_true(is.function(w2x$ssfct))
  }
})

# --- edfct (ED function) ---

test_that("weibull2x edfct works with b > 0 and control reference", {
  w2x <- weibull2x()
  # b=1 (>0), reference="control" -> p = 100-p, reference set to "upper"
  parm <- c(1, 0, 100, 5, 1)
  result <- w2x$edfct(parm, 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_true(length(result) == 2)
  expect_true(is.numeric(result[[1]]))
})

test_that("weibull2x edfct works with b < 0 and control reference", {
  w2x <- weibull2x()
  # b=-1 (<0), reference="control" -> p = 100-p (no reference change)
  parm <- c(-1, 0, 100, 5, 1)
  result <- w2x$edfct(parm, 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("weibull2x edfct works with non-control reference", {
  w2x <- weibull2x()
  parm <- c(1, 0, 100, 5, 1)
  result <- w2x$edfct(parm, 50, reference = "upper", type = "relative")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("weibull2x edfct adds t0 to the result", {
  # The ED result should be offset by the t0 value
  w2x_t0_1 <- weibull2x()
  parm_t0_1 <- c(1, 0, 100, 5, 1)
  result_1 <- w2x_t0_1$edfct(parm_t0_1, 50, reference = "upper", type = "relative")

  w2x_t0_2 <- weibull2x()
  parm_t0_2 <- c(1, 0, 100, 5, 2)
  result_2 <- w2x_t0_2$edfct(parm_t0_2, 50, reference = "upper", type = "relative")

  # Result with t0=2 should be 1 more than result with t0=1 (same underlying Weibull)
  expect_equal(result_2[[1]] - result_1[[1]], 1)
})

# --- W2x.3 wrapper ---

test_that("W2x.3 returns correct class and structure", {
  w2x3 <- W2x.3()
  expect_s3_class(w2x3, "Weibull-2")
  expect_equal(w2x3$noParm, 3)
  expect_equal(w2x3$names, c("d", "e", "t0"))
})

test_that("W2x.3 has b=1 and c=0 fixed", {
  w2x3 <- W2x.3()
  expect_equal(w2x3$fixed, c(1, 0, NA, NA, NA))
})

test_that("W2x.3 text indicates lower limit fixed", {
  w2x3 <- W2x.3()
  expect_true(grepl("lower limit at 0", w2x3$text))
})

test_that("W2x.3 errors on invalid names", {
  expect_error(W2x.3(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(W2x.3(names = 123), "Not correct 'names' argument")
})

test_that("W2x.3 errors on invalid fixed length", {
  expect_error(W2x.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
  expect_error(W2x.3(fixed = c(NA)), "Not correct length of 'fixed' argument")
})

test_that("W2x.3 with partial fixed values", {
  w2x3 <- W2x.3(fixed = c(100, NA, NA))
  expect_equal(w2x3$noParm, 2)
  expect_equal(w2x3$names, c("e", "t0"))
})

# --- W2x.4 wrapper ---

test_that("W2x.4 returns correct class and structure", {
  w2x4 <- W2x.4()
  expect_s3_class(w2x4, "Weibull-2")
  expect_equal(w2x4$noParm, 4)
  expect_equal(w2x4$names, c("c", "d", "e", "t0"))
})

test_that("W2x.4 has b=1 fixed", {
  w2x4 <- W2x.4()
  expect_equal(w2x4$fixed, c(1, NA, NA, NA, NA))
})

test_that("W2x.4 text indicates lower limit fixed", {
  w2x4 <- W2x.4()
  expect_true(grepl("lower limit at 0", w2x4$text))
})

test_that("W2x.4 errors on invalid names", {
  expect_error(W2x.4(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(W2x.4(names = 42), "Not correct 'names' argument")
})

test_that("W2x.4 errors on invalid fixed length", {
  expect_error(W2x.4(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
  expect_error(W2x.4(fixed = c(NA)), "Not correct length of 'fixed' argument")
})

test_that("W2x.4 with partial fixed values", {
  w2x4 <- W2x.4(fixed = c(0, 100, NA, NA))
  expect_equal(w2x4$noParm, 2)
  expect_equal(w2x4$names, c("e", "t0"))
})

# --- Integration tests using drm ---

test_that("W2x.4 works in drm model fitting", {
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
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = W2x.4())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 4)
  preds <- predict(m1)
  expect_true(all(is.finite(preds)))
})

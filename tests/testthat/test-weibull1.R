# Tests for weibull1.R: weibull1(), W1.2(), W1.3(), W1.3u(), W1.4()

# Create test dataset used throughout
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

# --- weibull1() main function ---

test_that("weibull1 returns correct class and structure", {
  w1 <- weibull1()
  expect_s3_class(w1, "Weibull-1")
  expect_true(is.list(w1))
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2", "derivx",
                     "edfct", "name", "text", "noParm") %in% names(w1)))
})

test_that("weibull1 default names are b, c, d, e", {
  w1 <- weibull1()
  expect_equal(w1$names, c("b", "c", "d", "e"))
})

test_that("weibull1 noParm reflects number of NA in fixed", {
  w1_full <- weibull1()
  expect_equal(w1_full$noParm, 4)

  w1_partial <- weibull1(fixed = c(1, NA, NA, NA))
  expect_equal(w1_partial$noParm, 3)
  expect_equal(w1_partial$names, c("c", "d", "e"))
})

test_that("weibull1 errors on invalid names argument", {
  expect_error(weibull1(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(weibull1(names = 123), "Not correct 'names' argument")
})

test_that("weibull1 errors on invalid fixed argument", {
  expect_error(weibull1(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(weibull1(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
})

test_that("weibull1 uses provided ssfct when not NULL", {
  custom_ss <- function(dframe) { c(1, 0, 10, 5) }
  w1 <- weibull1(ssfct = custom_ss)
  expect_identical(w1$ssfct, custom_ss)
})

# --- fct function tests ---

test_that("weibull1 fct computes correct values", {
  w1 <- weibull1()

  # f(x) = c + (d - c) * exp(-exp(b*(log(x) - log(e))))
  # b=1, c=0, d=1, e=2: f(2) = 0 + (1-0)*exp(-exp(1*(log(2)-log(2)))) = exp(-exp(0)) = exp(-1)
  dose <- 2
  parm <- matrix(c(1, 0, 1, 2), nrow = 1, ncol = 4)
  result <- w1$fct(dose, parm)
  expect_equal(as.numeric(result), exp(-1), tolerance = 1e-10)
})

test_that("weibull1 fct handles multiple doses", {
  w1 <- weibull1()

  dose <- c(1, 2, 4)
  parm <- matrix(c(1, 0, 1, 2), nrow = 3, ncol = 4, byrow = TRUE)
  result <- w1$fct(dose, parm)

  # Manual calculations
  expected <- exp(-exp(1 * (log(dose) - log(2))))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("weibull1 fct works with fixed parameters", {
  # Fix b=1 and c=0
  w1 <- weibull1(fixed = c(1, 0, NA, NA))
  dose <- c(2)
  parm <- matrix(c(1, 2), nrow = 1, ncol = 2)  # only d and e free
  result <- w1$fct(dose, parm)
  expected <- 0 + (1 - 0) * exp(-exp(1 * (log(2) - log(2))))
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

# --- deriv1 tests ---

test_that("weibull1 deriv1 returns matrix with correct dimensions", {
  w1 <- weibull1()
  dose <- c(1, 2, 5)
  parm <- matrix(c(1, 0, 1, 2), nrow = 3, ncol = 4, byrow = TRUE)
  result <- w1$deriv1(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
})

test_that("weibull1 deriv1 computes finite values", {
  w1 <- weibull1()
  dose <- c(0.5, 1, 3)
  parm <- matrix(c(2, 0, 1, 3), nrow = 3, ncol = 4, byrow = TRUE)
  result <- w1$deriv1(dose, parm)
  expect_true(all(is.finite(result)))
})

# --- derivx tests ---

test_that("weibull1 derivx returns correct structure", {
  w1 <- weibull1()
  dose <- c(1, 2, 5)
  parm <- matrix(c(1, 0, 1, 2), nrow = 3, ncol = 4, byrow = TRUE)
  result <- w1$derivx(dose, parm)
  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

# --- W1.2 wrapper ---

test_that("W1.2 returns correct class and structure", {
  w12 <- W1.2()
  expect_s3_class(w12, "Weibull-1")
  expect_equal(w12$noParm, 2)
  expect_equal(w12$names, c("b", "e"))
})

test_that("W1.2 text indicates fixed limits", {
  w12 <- W1.2(upper = 1)
  expect_true(grepl("lower limit at 0", w12$text))
  expect_true(grepl("upper limit at 1", w12$text))
})

test_that("W1.2 errors on invalid names", {
  expect_error(W1.2(names = c("x")), "Not correct 'names' argument")
  expect_error(W1.2(names = 99), "Not correct 'names' argument")
})

test_that("W1.2 errors on invalid fixed", {
  expect_error(W1.2(fixed = c(NA)), "Not correct length of 'fixed' argument")
})

# --- W1.3 wrapper ---

test_that("W1.3 returns correct class and structure", {
  w13 <- W1.3()
  expect_s3_class(w13, "Weibull-1")
  expect_equal(w13$noParm, 3)
  expect_equal(w13$names, c("b", "d", "e"))
})

test_that("W1.3 text indicates lower limit fixed", {
  w13 <- W1.3()
  expect_true(grepl("lower limit at 0", w13$text))
})

test_that("W1.3 errors on invalid names", {
  expect_error(W1.3(names = c("x")), "Not correct 'names' argument")
  expect_error(W1.3(names = 99), "Not correct 'names' argument")
})

test_that("W1.3 errors on invalid fixed", {
  expect_error(W1.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

# --- W1.3u wrapper ---

test_that("W1.3u returns correct class and structure", {
  w13u <- W1.3u()
  expect_s3_class(w13u, "Weibull-1")
  expect_equal(w13u$noParm, 3)
  expect_equal(w13u$names, c("b", "c", "e"))
})

test_that("W1.3u text indicates upper limit fixed", {
  w13u <- W1.3u(upper = 1)
  expect_true(grepl("upper limit at 1", w13u$text))
})

test_that("W1.3u errors on invalid names", {
  expect_error(W1.3u(names = c("x")), "Not correct 'names' argument")
  expect_error(W1.3u(names = 99), "Not correct 'names' argument")
})

test_that("W1.3u errors on invalid fixed", {
  expect_error(W1.3u(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

# --- W1.4 wrapper ---

test_that("W1.4 returns correct class and structure", {
  w14 <- W1.4()
  expect_s3_class(w14, "Weibull-1")
  expect_equal(w14$noParm, 4)
  expect_equal(w14$names, c("b", "c", "d", "e"))
})

test_that("W1.4 text is standard", {
  w14 <- W1.4()
  expect_equal(w14$text, "Weibull (type 1)")
})

test_that("W1.4 errors on invalid fixed", {
  expect_error(W1.4(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("W1.4 errors on invalid names", {
  expect_error(W1.4(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(W1.4(names = 123), "Not correct 'names' argument")
})

# --- Aliases ---

test_that("w2, w3, w4 are aliases for W1.2, W1.3, W1.4", {
  expect_identical(w2, W1.2)
  expect_identical(w3, W1.3)
  expect_identical(w4, W1.4)
})

# --- deriv2 is NULL ---

test_that("weibull1 deriv2 is NULL", {
  w1 <- weibull1()
  expect_null(w1$deriv2)
})

# --- Integration tests using drm ---

test_that("W1.4 works in drm model fitting", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 4)
  preds <- predict(m1)
  expect_true(all(is.finite(preds)))
})

test_that("W1.3 works in drm model fitting", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.3())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 3)
})

# --- EXD.2 wrapper ---

test_that("EXD.2 returns correct class and structure", {
  exd2 <- EXD.2()
  expect_s3_class(exd2, "Weibull-1")
  expect_equal(exd2$noParm, 2)
  expect_equal(exd2$names, c("d", "e"))
})

test_that("EXD.2 text indicates lower limit fixed", {
  exd2 <- EXD.2()
  expect_true(grepl("lower limit at 0", exd2$text))
})

test_that("EXD.2 errors on invalid names", {
  expect_error(EXD.2(names = c("x")), "Not correct 'names' argument")
  expect_error(EXD.2(names = 99), "Not correct 'names' argument")
})

test_that("EXD.2 errors on invalid fixed", {
  expect_error(EXD.2(fixed = c(NA)), "Not correct length of 'fixed' argument")
})

# --- EXD.3 wrapper ---

test_that("EXD.3 returns correct class and structure", {
  exd3 <- EXD.3()
  expect_s3_class(exd3, "Weibull-1")
  expect_equal(exd3$noParm, 3)
  expect_equal(exd3$names, c("c", "d", "e"))
})

test_that("EXD.3 text is correct", {
  exd3 <- EXD.3()
  expect_equal(exd3$text, "Shifted exponential decay")
})

test_that("EXD.3 errors on invalid fixed", {
  expect_error(EXD.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("EXD.3 errors on invalid names", {
  expect_error(EXD.3(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(EXD.3(names = 123), "Not correct 'names' argument")
})

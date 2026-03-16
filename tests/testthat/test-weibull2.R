# Tests for weibull2.R: weibull2(), W2.2(), W2.3(), W2.3u(), W2.4(), AR.2(), AR.3()

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

# --- weibull2() main function ---

test_that("weibull2 returns correct class and structure", {
  w2 <- weibull2()
  expect_s3_class(w2, "Weibull-2")
  expect_true(is.list(w2))
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2", "derivx",
                     "edfct", "name", "text", "noParm", "fixed") %in% names(w2)))
})

test_that("weibull2 default names are b, c, d, e", {
  w2 <- weibull2()
  expect_equal(w2$names, c("b", "c", "d", "e"))
})

test_that("weibull2 noParm reflects number of NA in fixed", {
  w2_full <- weibull2()
  expect_equal(w2_full$noParm, 4)

  w2_partial <- weibull2(fixed = c(1, NA, NA, NA))
  expect_equal(w2_partial$noParm, 3)
  expect_equal(w2_partial$names, c("c", "d", "e"))
})

test_that("weibull2 uses default text when fctText not provided", {
  w2 <- weibull2()
  expect_equal(w2$text, "Weibull (type 2)")
})

test_that("weibull2 uses provided fctText", {
  w2 <- weibull2(fctText = "Custom text")
  expect_equal(w2$text, "Custom text")
})

test_that("weibull2 uses provided fctName", {
  w2 <- weibull2(fctName = "myFunc")
  expect_equal(w2$name, "myFunc")
})

test_that("weibull2 errors on invalid names argument", {
  expect_error(weibull2(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(weibull2(names = 123), "Not correct 'names' argument")
})

test_that("weibull2 errors on invalid fixed argument", {
  expect_error(weibull2(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(weibull2(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
})

test_that("weibull2 uses provided ssfct when not NULL", {
  custom_ss <- function(dframe) { c(1, 0, 10, 5) }
  w2 <- weibull2(ssfct = custom_ss)
  expect_identical(w2$ssfct, custom_ss)
})

test_that("weibull2 uses default ssfct when ssfct is NULL", {
  w2 <- weibull2(ssfct = NULL)
  expect_true(is.function(w2$ssfct))
})

test_that("weibull2 method argument works for self-starter", {
  for (m in c("1", "2", "3", "4")) {
    w2 <- weibull2(method = m)
    expect_true(is.function(w2$ssfct))
  }
})

# --- fct (internal nonlinear function) ---

test_that("weibull2 fct computes correct values", {
  w2 <- weibull2()
  # Parameters: b=1, c=0, d=100, e=5
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5), nrow = 3, ncol = 4, byrow = TRUE)
  result <- w2$fct(dose, parm)
  # f(x) = c + (d - c)(1 - exp(-exp(b(log(x) - log(e)))))
  expected <- 0 + (100 - 0) * (1 - exp(-exp(1 * (log(dose) - log(5)))))
  expect_equal(result, expected)
})

test_that("weibull2 fct works with fixed parameters", {
  w2 <- weibull2(fixed = c(1, 0, NA, NA))
  dose <- c(1, 5, 10)
  parm <- matrix(c(100, 5), nrow = 3, ncol = 2, byrow = TRUE)
  result <- w2$fct(dose, parm)
  expected <- 0 + (100 - 0) * (1 - exp(-exp(1 * (log(dose) - log(5)))))
  expect_equal(result, expected)
})

# --- deriv1 (parameter derivatives) ---

test_that("weibull2 deriv1 returns matrix with correct dimensions", {
  w2 <- weibull2()
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5), nrow = 3, ncol = 4, byrow = TRUE)
  result <- w2$deriv1(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
})

test_that("weibull2 deriv1 works with fixed parameters", {
  w2 <- weibull2(fixed = c(1, NA, NA, NA))
  dose <- c(1, 5, 10)
  parm <- matrix(c(0, 100, 5), nrow = 3, ncol = 3, byrow = TRUE)
  result <- w2$deriv1(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
})

# --- derivx (dose derivative) ---

test_that("weibull2 derivx returns correct structure", {
  w2 <- weibull2()
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5), nrow = 3, ncol = 4, byrow = TRUE)
  result <- w2$derivx(dose, parm)
  expect_true(is.matrix(result) || is.array(result))
  expect_equal(nrow(result), 3)
})

test_that("weibull2 derivx works with fixed parameters", {
  w2 <- weibull2(fixed = c(1, 0, NA, NA))
  dose <- c(1, 5, 10)
  parm <- matrix(c(100, 5), nrow = 3, ncol = 2, byrow = TRUE)
  result <- w2$derivx(dose, parm)
  expect_true(is.matrix(result) || is.array(result))
})

# --- edfct (ED function) ---

test_that("weibull2 edfct works with relative type", {
  w2 <- weibull2()
  parm <- c(1, 0, 100, 5)
  result <- w2$edfct(parm, 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_true(length(result) == 2)
  expect_true(is.numeric(result[[1]]))
})

test_that("weibull2 edfct works with absolute type and b>0 and control", {
  w2 <- weibull2()
  # b=1 (>0), reference = "control", type = "absolute"
  parm <- c(1, 0, 100, 5)
  result <- w2$edfct(parm, 50, reference = "control", type = "absolute")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("weibull2 edfct with absolute type, b<=0, control", {
  w2 <- weibull2()
  # b=-1 (<0), reference = "control", type = "absolute"
  parm <- c(-1, 0, 100, 5)
  result <- w2$edfct(parm, 50, reference = "control", type = "absolute")
  expect_true(is.list(result))
})

test_that("weibull2 edfct with absolute type and non-control reference", {
  w2 <- weibull2()
  parm <- c(1, 0, 100, 5)
  result <- w2$edfct(parm, 50, reference = "upper", type = "absolute")
  expect_true(is.list(result))
})

# --- W2.2 wrapper ---

test_that("W2.2 returns correct class and structure", {
  w22 <- W2.2()
  expect_s3_class(w22, "Weibull-2")
  expect_equal(w22$noParm, 2)
  expect_equal(w22$names, c("b", "e"))
})

test_that("W2.2 with custom upper limit", {
  w22 <- W2.2(upper = 100)
  expect_s3_class(w22, "Weibull-2")
  expect_true(grepl("upper limit at 100", w22$text))
})

test_that("W2.2 errors on invalid names", {
  expect_error(W2.2(names = c("a")), "Not correct 'names' argument")
  expect_error(W2.2(names = 42), "Not correct 'names' argument")
})

test_that("W2.2 errors on invalid fixed", {
  expect_error(W2.2(fixed = c(NA)), "Not correct length of 'fixed' argument")
  expect_error(W2.2(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("W2.2 fct computes correctly", {
  w22 <- W2.2(upper = 1)
  dose <- c(0.5, 1, 2)
  parm <- matrix(c(1, 5), nrow = 3, ncol = 2, byrow = TRUE)
  result <- w22$fct(dose, parm)
  expected <- 0 + (1 - 0) * (1 - exp(-exp(1 * (log(dose) - log(5)))))
  expect_equal(result, expected)
})

# --- W2.3 wrapper ---

test_that("W2.3 returns correct class and structure", {
  w23 <- W2.3()
  expect_s3_class(w23, "Weibull-2")
  expect_equal(w23$noParm, 3)
  expect_equal(w23$names, c("b", "d", "e"))
})

test_that("W2.3 text indicates lower limit fixed at 0", {
  w23 <- W2.3()
  expect_true(grepl("lower limit at 0", w23$text))
})

test_that("W2.3 errors on invalid names", {
  expect_error(W2.3(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(W2.3(names = 123), "Not correct 'names' argument")
})

test_that("W2.3 errors on invalid fixed", {
  expect_error(W2.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

# --- W2.3u wrapper ---

test_that("W2.3u returns correct class and structure", {
  w23u <- W2.3u()
  expect_s3_class(w23u, "Weibull-2")
  expect_equal(w23u$noParm, 3)
  expect_equal(w23u$names, c("b", "c", "e"))
})

test_that("W2.3u text indicates upper limit fixed", {
  w23u <- W2.3u(upper = 1)
  expect_true(grepl("upper limit at 1", w23u$text))
})

test_that("W2.3u errors on invalid names", {
  expect_error(W2.3u(names = c("x")), "Not correct 'names' argument")
  expect_error(W2.3u(names = 99), "Not correct 'names' argument")
})

test_that("W2.3u errors on invalid fixed", {
  expect_error(W2.3u(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

# --- W2.4 wrapper ---

test_that("W2.4 returns correct class and structure", {
  w24 <- W2.4()
  expect_s3_class(w24, "Weibull-2")
  expect_equal(w24$noParm, 4)
  expect_equal(w24$names, c("b", "c", "d", "e"))
})

test_that("W2.4 text is standard", {
  w24 <- W2.4()
  expect_equal(w24$text, "Weibull (type 2)")
})

test_that("W2.4 errors on invalid fixed", {
  expect_error(W2.4(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("W2.4 errors on invalid names", {
  expect_error(W2.4(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(W2.4(names = 123), "Not correct 'names' argument")
})

# --- AR.2 wrapper ---

test_that("AR.2 returns correct class and structure", {
  ar2 <- AR.2()
  expect_s3_class(ar2, "Weibull-2")
  expect_equal(ar2$noParm, 2)
  expect_equal(ar2$names, c("d", "e"))
})

test_that("AR.2 text indicates asymptotic regression with lower fixed", {
  ar2 <- AR.2()
  expect_true(grepl("Asymptotic regression", ar2$text))
  expect_true(grepl("lower limit at 0", ar2$text))
})

test_that("AR.2 errors on invalid names", {
  expect_error(AR.2(names = c("a")), "Not correct 'names' argument")
  expect_error(AR.2(names = 42), "Not correct 'names' argument")
})

test_that("AR.2 errors on invalid fixed", {
  expect_error(AR.2(fixed = c(NA)), "Not correct length of 'fixed' argument")
})

# --- AR.3 wrapper ---

test_that("AR.3 returns correct class and structure", {
  ar3 <- AR.3()
  expect_s3_class(ar3, "Weibull-2")
  expect_equal(ar3$noParm, 3)
  expect_equal(ar3$names, c("c", "d", "e"))
})

test_that("AR.3 text is shifted asymptotic regression", {
  ar3 <- AR.3()
  expect_equal(ar3$text, "Shifted asymptotic regression")
})

test_that("AR.3 errors on invalid names", {
  expect_error(AR.3(names = c("a")), "Not correct 'names' argument")
  expect_error(AR.3(names = 42), "Not correct 'names' argument")
})

test_that("AR.3 errors on invalid fixed", {
  expect_error(AR.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

# --- Integration tests using drm ---

test_that("W2.4 works in drm model fitting", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = W2.4())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 4)
  preds <- predict(m1)
  expect_true(all(is.finite(preds)))
})

test_that("W2.3 works in drm model fitting", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = W2.3())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 3)
})

test_that("AR.2 works in drm model fitting", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = AR.2())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 2)
})

test_that("AR.3 works in drm model fitting", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = AR.3())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 3)
})

# --- deriv2 is NULL ---

test_that("weibull2 deriv2 is NULL", {
  w2 <- weibull2()
  expect_null(w2$deriv2)
})

# --- fixed field preserved ---

test_that("weibull2 fixed field is preserved", {
  fixed_vals <- c(1, NA, NA, NA)
  w2 <- weibull2(fixed = fixed_vals)
  expect_equal(w2$fixed, fixed_vals)
})

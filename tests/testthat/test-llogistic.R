# Tests for llogistic.R: llogistic(), LL.2(), LL.3(), LL.3u(), LL.4(), LL.5()

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

# --- llogistic() main function ---

test_that("llogistic returns correct class and structure", {
  ll <- llogistic()
  expect_s3_class(ll, "log-logistic")
  expect_true(is.list(ll))
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2",
                     "edfct", "name", "text", "noParm") %in% names(ll)))
})

test_that("llogistic default names are b, c, d, e, f", {
  ll <- llogistic()
  expect_equal(ll$names, c("b", "c", "d", "e", "f"))
})

test_that("llogistic noParm reflects number of NA in fixed", {
  ll_full <- llogistic()
  expect_equal(ll_full$noParm, 5)

  ll_partial <- llogistic(fixed = c(1, NA, NA, NA, NA))
  expect_equal(ll_partial$noParm, 4)
  expect_equal(ll_partial$names, c("c", "d", "e", "f"))
})

test_that("llogistic errors on invalid names argument", {
  expect_error(llogistic(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(llogistic(names = 123), "Not correct 'names' argument")
})

test_that("llogistic errors on invalid fixed argument", {
  expect_error(llogistic(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(llogistic(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
})

# --- fct function tests ---

test_that("llogistic fct computes correct values", {
  ll <- llogistic()

  # f(x) = c + (d-c) / (1 + exp(b*log(x/e)))^f
  # b=1, c=0, d=1, e=2, f=1: f(2) = 0 + (1-0) / (1+exp(1*log(2/2)))^1 = 1 / (1+1)^1 = 0.5
  dose <- 2
  parm <- matrix(c(1, 0, 1, 2, 1), nrow = 1, ncol = 5)
  result <- ll$fct(dose, parm)
  expect_equal(as.numeric(result), 0.5, tolerance = 1e-10)
})

test_that("llogistic fct handles multiple doses", {
  ll <- llogistic()

  dose <- c(1, 2, 4)
  parm <- matrix(c(1, 0, 1, 2, 1), nrow = 3, ncol = 5, byrow = TRUE)
  result <- ll$fct(dose, parm)

  # Manual calculations: f(x) = 1/(1 + exp(log(x/2)))^1 = 1/(1 + x/2)
  expected <- 1 / (1 + dose / 2)
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("llogistic fct at ED50 gives midpoint (f=1)", {
  ll <- llogistic()

  # At dose = e, response should be (c + d)/2 for f=1
  dose <- 5
  parm <- matrix(c(2, 0, 10, 5, 1), nrow = 1, ncol = 5)
  result <- ll$fct(dose, parm)
  expected <- (0 + 10) / 2
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
})

test_that("llogistic fct works with f != 1 (asymmetric)", {
  ll <- llogistic()

  # b=1, c=0, d=1, e=1, f=2
  dose <- 1
  parm <- matrix(c(1, 0, 1, 1, 2), nrow = 1, ncol = 5)
  result <- ll$fct(dose, parm)
  # f(1) = 0 + (1-0)/(1 + exp(0))^2 = 1/4
  expect_equal(as.numeric(result), 0.25, tolerance = 1e-10)
})

# --- LL.2 wrapper ---

test_that("LL.2 returns correct class and structure", {
  ll2 <- LL.2()
  expect_s3_class(ll2, "log-logistic")
  expect_equal(ll2$noParm, 2)
  expect_equal(ll2$names, c("b", "e"))
})

test_that("LL.2 text indicates fixed limits", {
  ll2 <- LL.2(upper = 1)
  expect_true(grepl("lower limit at 0", ll2$text))
  expect_true(grepl("upper limit at 1", ll2$text))
})

test_that("LL.2 errors on invalid names", {
  expect_error(LL.2(names = c("x")), "Not correct 'names' argument")
  expect_error(LL.2(names = 99), "Not correct 'names' argument")
})

test_that("LL.2 errors on invalid fixed", {
  expect_error(LL.2(fixed = c(NA)), "Not correct length of 'fixed' argument")
})

# --- LL.3 wrapper ---

test_that("LL.3 returns correct class and structure", {
  ll3 <- LL.3()
  expect_s3_class(ll3, "log-logistic")
  expect_equal(ll3$noParm, 3)
  expect_equal(ll3$names, c("b", "d", "e"))
})

test_that("LL.3 text indicates lower limit fixed", {
  ll3 <- LL.3()
  expect_true(grepl("lower limit at 0", ll3$text))
})

test_that("LL.3 errors on invalid names", {
  expect_error(LL.3(names = c("x")), "Not correct 'names' argument")
  expect_error(LL.3(names = 99), "Not correct 'names' argument")
})

test_that("LL.3 errors on invalid fixed", {
  expect_error(LL.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

# --- LL.3u wrapper ---

test_that("LL.3u returns correct class and structure", {
  ll3u <- LL.3u()
  expect_s3_class(ll3u, "log-logistic")
  expect_equal(ll3u$noParm, 3)
  expect_equal(ll3u$names, c("b", "c", "e"))
})

test_that("LL.3u text indicates upper limit fixed", {
  ll3u <- LL.3u(upper = 1)
  expect_true(grepl("upper limit at 1", ll3u$text))
})

test_that("LL.3u errors on invalid names", {
  expect_error(LL.3u(names = c("x")), "Not correct 'names' argument")
  expect_error(LL.3u(names = 99), "Not correct 'names' argument")
})

test_that("LL.3u errors on invalid fixed", {
  expect_error(LL.3u(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

# --- LL.4 wrapper ---

test_that("LL.4 returns correct class and structure", {
  ll4 <- LL.4()
  expect_s3_class(ll4, "log-logistic")
  expect_equal(ll4$noParm, 4)
  expect_equal(ll4$names, c("b", "c", "d", "e"))
})

test_that("LL.4 errors on invalid fixed", {
  expect_error(LL.4(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("LL.4 errors on invalid names", {
  expect_error(LL.4(names = c("a", "b")), "Not correct names argument")
  expect_error(LL.4(names = 123), "Not correct names argument")
})

# --- LL.5 wrapper ---

test_that("LL.5 returns correct class and structure", {
  ll5 <- LL.5()
  expect_s3_class(ll5, "log-logistic")
  expect_equal(ll5$noParm, 5)
  expect_equal(ll5$names, c("b", "c", "d", "e", "f"))
})

test_that("LL.5 text indicates generalized model", {
  ll5 <- LL.5()
  expect_true(grepl("Generalized", ll5$text))
})

# --- Aliases ---

test_that("l2, l3, l3u, l4, l5 are aliases for LL.2, LL.3, LL.3u, LL.4, LL.5", {
  expect_identical(l2, LL.2)
  expect_identical(l3, LL.3)
  expect_identical(l3u, LL.3u)
  expect_identical(l4, LL.4)
  expect_identical(l5, LL.5)
})

# --- deriv1 tests ---

test_that("llogistic deriv1 returns matrix with correct dimensions", {
  ll <- llogistic()
  dose <- c(1, 2, 5)
  parm <- matrix(c(1, 0, 1, 2, 1), nrow = 3, ncol = 5, byrow = TRUE)
  result <- ll$deriv1(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 5)
})

test_that("llogistic deriv1 computes finite values", {
  ll <- llogistic()
  dose <- c(0.5, 1, 3)
  parm <- matrix(c(2, 0, 1, 3, 1), nrow = 3, ncol = 5, byrow = TRUE)
  result <- ll$deriv1(dose, parm)
  expect_true(all(is.finite(result)))
})

# --- Integration tests using drm ---

test_that("LL.4 works in drm model fitting", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 4)
  preds <- predict(m1)
  expect_true(all(is.finite(preds)))
})

test_that("LL.3 works in drm model fitting", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 3)
})

test_that("LL.5 works in drm model fitting", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.5())
  expect_s3_class(m1, "drc")
  expect_equal(length(coef(m1)), 5)
})

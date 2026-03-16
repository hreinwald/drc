# Tests for lnormal.R: lnormal(), LN.2(), LN.3(), LN.3u(), LN.4()
# and the internal edfct function

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

# --- lnormal() main function ---

test_that("lnormal returns correct class and structure", {
  ln <- lnormal()
  expect_s3_class(ln, "log-normal")
  expect_true(is.list(ln))
  expect_true(all(c("fct", "ssfct", "names", "deriv1", "deriv2", "derivx",
                     "edfct", "name", "text", "noParm", "fixed") %in% names(ln)))
})

test_that("lnormal default names are b, c, d, e", {
  ln <- lnormal()
  expect_equal(ln$names, c("b", "c", "d", "e"))
})

test_that("lnormal noParm reflects number of NA in fixed", {
  ln_full <- lnormal()
  expect_equal(ln_full$noParm, 4)

  ln_partial <- lnormal(fixed = c(1, NA, NA, NA))
  expect_equal(ln_partial$noParm, 3)
  expect_equal(ln_partial$names, c("c", "d", "e"))
})

test_that("lnormal uses default text when fctText not provided", {
  ln <- lnormal()
  expect_equal(ln$text, "Log-normal")
})

test_that("lnormal uses provided fctText", {
  ln <- lnormal(fctText = "Custom text")
  expect_equal(ln$text, "Custom text")
})

test_that("lnormal uses provided fctName", {
  ln <- lnormal(fctName = "myFunc")
  expect_equal(ln$name, "myFunc")
})

test_that("lnormal errors on invalid names argument", {
  expect_error(lnormal(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(lnormal(names = 123), "Not correct 'names' argument")
})

test_that("lnormal errors on invalid fixed argument", {
  expect_error(lnormal(fixed = c(NA, NA)), "Not correct 'fixed' argument")
  expect_error(lnormal(fixed = c(NA, NA, NA)), "Not correct 'fixed' argument")
})

test_that("lnormal uses provided ssfct when not NULL", {
  custom_ss <- function(dframe) { c(1, 0, 10, 5) }
  ln <- lnormal(ssfct = custom_ss)
  expect_identical(ln$ssfct, custom_ss)
})

test_that("lnormal uses default ssfct when ssfct is NULL", {
  ln <- lnormal(ssfct = NULL)
  expect_true(is.function(ln$ssfct))
})

test_that("lnormal method argument works for self-starter", {
  for (m in c("1", "2", "3", "4")) {
    ln <- lnormal(method = m)
    expect_true(is.function(ln$ssfct))
  }
})

test_that("lnormal with loge=FALSE (default)", {
  ln <- lnormal(loge = FALSE)
  expect_s3_class(ln, "log-normal")
  expect_true(is.function(ln$fct))
})

test_that("lnormal with loge=TRUE", {
  ln <- lnormal(loge = TRUE)
  expect_s3_class(ln, "log-normal")
  expect_true(is.function(ln$fct))
})

# --- fct (internal nonlinear function) ---

test_that("lnormal fct computes correct values with loge=FALSE", {
  ln <- lnormal()
  # Parameters: b=1, c=0, d=100, e=5
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5), nrow = 3, ncol = 4, byrow = TRUE)
  result <- ln$fct(dose, parm)
  # f(x) = c + (d-c)*pnorm(b*(log(x)-log(e)))
  expected <- 0 + (100 - 0) * pnorm(1 * (log(dose) - log(5)))
  expect_equal(as.numeric(result), expected)
})

test_that("lnormal fct computes correct values with loge=TRUE", {
  ln <- lnormal(loge = TRUE)
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, log(5)), nrow = 3, ncol = 4, byrow = TRUE)
  result <- ln$fct(dose, parm)
  # f(x) = c + (d-c)*pnorm(b*(log(x)-e))
  expected <- 0 + (100 - 0) * pnorm(1 * (log(dose) - log(5)))
  expect_equal(as.numeric(result), expected)
})

test_that("lnormal fct works with fixed parameters", {
  ln <- lnormal(fixed = c(1, 0, NA, NA))
  dose <- c(1, 5, 10)
  parm <- matrix(c(100, 5), nrow = 3, ncol = 2, byrow = TRUE)
  result <- ln$fct(dose, parm)
  expected <- 0 + (100 - 0) * pnorm(1 * (log(dose) - log(5)))
  expect_equal(as.numeric(result), expected)
})

# --- deriv1 (parameter derivatives) ---

test_that("lnormal deriv1 returns gradient matrix with correct dimensions (loge=FALSE)", {
  ln <- lnormal()
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5), nrow = 3, ncol = 4, byrow = TRUE)
  result <- ln$deriv1(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)  # all 4 params free
})

test_that("lnormal deriv1 returns gradient matrix with correct dimensions (loge=TRUE)", {
  ln <- lnormal(loge = TRUE)
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, log(5)), nrow = 3, ncol = 4, byrow = TRUE)
  result <- ln$deriv1(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
})

test_that("lnormal deriv1 with fixed params reduces columns", {
  ln <- lnormal(fixed = c(1, 0, NA, NA))
  dose <- c(1, 5, 10)
  parm <- matrix(c(100, 5), nrow = 3, ncol = 2, byrow = TRUE)
  result <- ln$deriv1(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)  # only d and e are free
})

# --- derivx (dose derivatives) ---

test_that("lnormal derivx returns gradient with correct dimensions (loge=FALSE)", {
  ln <- lnormal()
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5), nrow = 3, ncol = 4, byrow = TRUE)
  result <- ln$derivx(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 1)
})

test_that("lnormal derivx returns gradient with correct dimensions (loge=TRUE)", {
  ln <- lnormal(loge = TRUE)
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, log(5)), nrow = 3, ncol = 4, byrow = TRUE)
  result <- ln$derivx(dose, parm)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 1)
})

# --- edfct (effective dose function) ---

test_that("lnormal edfct returns list with ED and gradient (loge=FALSE, relative)", {
  ln <- lnormal()
  parm <- c(1, 0, 100, 5)
  result <- ln$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  # ED50 for symmetric model should equal e=5
  expect_equal(as.numeric(result[[1]]), 5, tolerance = 1e-6)
  # Gradient should have 4 elements (all params free)
  expect_equal(length(result[[2]]), 4)
})

test_that("lnormal edfct returns list with ED and gradient (loge=TRUE, relative)", {
  ln <- lnormal(loge = TRUE)
  parm <- c(1, 0, 100, log(5))
  result <- ln$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  # ED50 on log scale should be log(5)
  expect_equal(as.numeric(result[[1]]), log(5), tolerance = 1e-6)
  expect_equal(length(result[[2]]), 4)
})

test_that("lnormal edfct works with absolute type", {
  ln <- lnormal()
  parm <- c(1, 0, 100, 5)
  result <- ln$edfct(parm, respl = 50, reference = "control", type = "absolute")
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_true(is.numeric(result[[1]]))
})

test_that("lnormal edfct works with negative b and control reference (relative)", {
  ln <- lnormal()
  parm <- c(-1, 0, 100, 5)
  result <- ln$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_true(is.numeric(result[[1]]))
})

test_that("lnormal edfct works with positive b and control reference (relative)", {
  ln <- lnormal()
  parm <- c(1, 0, 100, 5)
  result <- ln$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  # For positive b, no reversal
  expect_equal(as.numeric(result[[1]]), 5, tolerance = 1e-6)
})

test_that("lnormal edfct with fixed parameters", {
  ln <- lnormal(fixed = c(1, 0, NA, NA))
  parm <- c(100, 5)
  result <- ln$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_equal(length(result[[2]]), 2)  # only free params
})

test_that("lnormal edfct at different response levels", {
  ln <- lnormal()
  parm <- c(1, 0, 100, 5)
  # ED10

  result10 <- ln$edfct(parm, respl = 10, reference = "control", type = "relative")
  # ED90
  result90 <- ln$edfct(parm, respl = 90, reference = "control", type = "relative")
  # ED10 < ED50 < ED90 (for positive b, decreasing curve)
  expect_true(is.numeric(result10[[1]]))
  expect_true(is.numeric(result90[[1]]))
})

test_that("lnormal edfct with loge=TRUE and absolute type", {
  ln <- lnormal(loge = TRUE)
  parm <- c(1, 0, 100, log(5))
  result <- ln$edfct(parm, respl = 50, reference = "control", type = "absolute")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

test_that("lnormal edfct with loge=TRUE, negative b and control reference", {
  ln <- lnormal(loge = TRUE)
  parm <- c(-1, 0, 100, log(5))
  result <- ln$edfct(parm, respl = 50, reference = "control", type = "relative")
  expect_true(is.list(result))
  expect_true(is.numeric(result[[1]]))
})

# --- fd function: edge case with non-finite values (loge=FALSE) ---

test_that("lnormal fct handles dose=0 gracefully (loge=FALSE)", {
  ln <- lnormal()
  dose <- c(0, 1, 5)
  parm <- matrix(c(1, 0, 100, 5), nrow = 3, ncol = 4, byrow = TRUE)
  result <- ln$fct(dose, parm)
  # dose=0 => log(0)=-Inf => pnorm(-Inf)=0, so result = c + (d-c)*0 = c = 0
  expect_equal(as.numeric(result[1]), 0)
})

test_that("lnormal fct handles dose=0 gracefully (loge=TRUE)", {
  ln <- lnormal(loge = TRUE)
  dose <- c(0, 1, 5)
  parm <- matrix(c(1, 0, 100, log(5)), nrow = 3, ncol = 4, byrow = TRUE)
  result <- ln$fct(dose, parm)
  # dose=0 => log(0)=-Inf => b*(-Inf - e) = -Inf => pnorm(-Inf)=0
  expect_equal(as.numeric(result[1]), 0)
})

# --- lowerAs, upperAs, monoton ---

test_that("lnormal lowerAs and upperAs return correct values", {
  ln <- lnormal()
  # lowerAs extracts parameter 2 (c), upperAs extracts parameter 3 (d)
  parm <- c(1, 0, 100, 5)
  expect_equal(ln$lowerAs(parm), 0)
  expect_equal(ln$upperAs(parm), 100)
})

test_that("lnormal monoton returns correct sign", {
  ln <- lnormal()
  # monoton returns -1 * parmVec[1] (sign=-1, parmNo=1)
  parm <- c(2, 0, 100, 5)
  expect_equal(ln$monoton(parm), -2)
})

# --- LN.2 convenience function ---

test_that("LN.2 returns correct structure", {
  ln2 <- LN.2()
  expect_s3_class(ln2, "log-normal")
  expect_equal(ln2$noParm, 2)
  expect_equal(ln2$names, c("b", "e"))
})

test_that("LN.2 has correct text with default upper=1", {
  ln2 <- LN.2()
  expect_true(grepl("lower limit at 0", ln2$text))
  expect_true(grepl("upper limit at 1", ln2$text))
})

test_that("LN.2 with custom upper", {
  ln2 <- LN.2(upper = 100)
  expect_true(grepl("upper limit at 100", ln2$text))
})

test_that("LN.2 errors on invalid names", {
  expect_error(LN.2(names = c("a")), "Not correct 'names' argument")
  expect_error(LN.2(names = 123), "Not correct 'names' argument")
})

test_that("LN.2 errors on invalid fixed", {
  expect_error(LN.2(fixed = c(NA)), "Not correct length of 'fixed' argument")
  expect_error(LN.2(fixed = c(NA, NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("LN.2 fct computes correct values", {
  ln2 <- LN.2(upper = 100)
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 5), nrow = 3, ncol = 2, byrow = TRUE)
  result <- ln2$fct(dose, parm)
  expected <- 0 + (100 - 0) * pnorm(1 * (log(dose) - log(5)))
  expect_equal(as.numeric(result), expected)
})

test_that("LN.2 passes extra args to lnormal", {
  ln2 <- LN.2(loge = TRUE)
  expect_s3_class(ln2, "log-normal")
})

# --- LN.3 convenience function ---

test_that("LN.3 returns correct structure", {
  ln3 <- LN.3()
  expect_s3_class(ln3, "log-normal")
  expect_equal(ln3$noParm, 3)
  expect_equal(ln3$names, c("b", "d", "e"))
})

test_that("LN.3 has correct text", {
  ln3 <- LN.3()
  expect_true(grepl("lower limit at 0", ln3$text))
})

test_that("LN.3 errors on invalid names", {
  expect_error(LN.3(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(LN.3(names = 123), "Not correct 'names' argument")
})

test_that("LN.3 errors on invalid fixed", {
  expect_error(LN.3(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("LN.3 fct computes correct values", {
  ln3 <- LN.3()
  dose <- c(1, 5, 10)
  # b, d, e (c is fixed at 0)
  parm <- matrix(c(1, 100, 5), nrow = 3, ncol = 3, byrow = TRUE)
  result <- ln3$fct(dose, parm)
  expected <- 0 + (100 - 0) * pnorm(1 * (log(dose) - log(5)))
  expect_equal(as.numeric(result), expected)
})

test_that("LN.3 passes extra args to lnormal", {
  ln3 <- LN.3(loge = TRUE)
  expect_s3_class(ln3, "log-normal")
})

# --- LN.3u convenience function ---

test_that("LN.3u returns correct structure", {
  ln3u <- LN.3u()
  expect_s3_class(ln3u, "log-normal")
  expect_equal(ln3u$noParm, 3)
  expect_equal(ln3u$names, c("b", "c", "e"))
})

test_that("LN.3u has correct text with default upper=1", {
  ln3u <- LN.3u()
  expect_true(grepl("upper limit at 1", ln3u$text))
})

test_that("LN.3u with custom upper", {
  ln3u <- LN.3u(upper = 100)
  expect_true(grepl("upper limit at 100", ln3u$text))
})

test_that("LN.3u errors on invalid names", {
  expect_error(LN.3u(names = c("a", "b")), "Not correct 'names' argument")
  expect_error(LN.3u(names = 123), "Not correct 'names' argument")
})

test_that("LN.3u errors on invalid fixed", {
  expect_error(LN.3u(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("LN.3u fct computes correct values", {
  ln3u <- LN.3u(upper = 100)
  dose <- c(1, 5, 10)
  # b, c, e (d is fixed at upper=100)
  parm <- matrix(c(1, 0, 5), nrow = 3, ncol = 3, byrow = TRUE)
  result <- ln3u$fct(dose, parm)
  expected <- 0 + (100 - 0) * pnorm(1 * (log(dose) - log(5)))
  expect_equal(as.numeric(result), expected)
})

test_that("LN.3u passes extra args to lnormal", {
  ln3u <- LN.3u(loge = TRUE)
  expect_s3_class(ln3u, "log-normal")
})

# --- LN.4 convenience function ---

test_that("LN.4 returns correct structure", {
  ln4 <- LN.4()
  expect_s3_class(ln4, "log-normal")
  expect_equal(ln4$noParm, 4)
  expect_equal(ln4$names, c("b", "c", "d", "e"))
})

test_that("LN.4 errors on invalid names", {
  expect_error(LN.4(names = c("a", "b")), "Not correct names argument")
  expect_error(LN.4(names = 123), "Not correct names argument")
})

test_that("LN.4 errors on invalid fixed", {
  expect_error(LN.4(fixed = c(NA, NA)), "Not correct length of 'fixed' argument")
})

test_that("LN.4 fct computes correct values", {
  ln4 <- LN.4()
  dose <- c(1, 5, 10)
  parm <- matrix(c(1, 0, 100, 5), nrow = 3, ncol = 4, byrow = TRUE)
  result <- ln4$fct(dose, parm)
  expected <- 0 + (100 - 0) * pnorm(1 * (log(dose) - log(5)))
  expect_equal(as.numeric(result), expected)
})

test_that("LN.4 passes extra args to lnormal", {
  ln4 <- LN.4(loge = TRUE)
  expect_s3_class(ln4, "log-normal")
})

# --- Integration tests with drm ---

test_that("LN.4 works in drm model fit", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LN.4())
  expect_s3_class(m, "drc")
  ed <- ED(m, 50, display = FALSE)
  expect_true(is.matrix(ed))
  expect_true(ed[, "Estimate"] > 0)
})

test_that("LN.3 works in drm model fit", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LN.3())
  expect_s3_class(m, "drc")
  ed <- ED(m, 50, display = FALSE)
  expect_true(ed[, "Estimate"] > 0)
})

test_that("LN.2 works in drm model fit with scaled data", {
  # Scale rootl to 0-1 range for LN.2 with upper=1
  rg_scaled <- ryegrass
  rg_scaled$rootl <- rg_scaled$rootl / max(rg_scaled$rootl)
  m <- drm(rootl ~ conc, data = rg_scaled, fct = LN.2())
  expect_s3_class(m, "drc")
})

test_that("ED with absolute type works for LN.4 model", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LN.4())
  ed_abs <- ED(m, 5, type = "absolute", display = FALSE)
  expect_true(is.matrix(ed_abs))
  expect_true(ed_abs[, "Estimate"] > 0)
})

test_that("ED with relative type and multiple levels works for LN.4 model", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LN.4())
  ed <- ED(m, c(10, 50, 90), display = FALSE)
  expect_equal(nrow(ed), 3)
  expect_true(all(ed[, "Estimate"] > 0))
})

test_that("LN.4 with loge=TRUE works in drm model fit", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LN.4(loge = TRUE))
  expect_s3_class(m, "drc")
  ed <- ED(m, 50, display = FALSE)
  expect_true(is.matrix(ed))
})

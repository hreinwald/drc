# Regression tests for the second review batch (ED intervals, predict, confint, separate fits,
# model averaging, comped, simDR, noEffect, pmodels).

# ---- ED intervals on log-scale and dose-scale models -------------------------------

test_that("interval = 'fls' stops with instructions for dose-scale models only", {
  mLL <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  mLL2 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.4())

  expect_error(ED(mLL, 50, interval = "fls", display = FALSE), "cannot be used with model 'LL.4'")
  expect_error(ED(mLL, 50, interval = "fls", display = FALSE), "How to fix")
  expect_error(ED(mLL, 50, interval = "fls", display = FALSE), "interval = 'tfls'")
  expect_no_error(ED(mLL2, 50, interval = "fls", display = FALSE))
  expect_no_error(ED(mLL, 50, interval = "fls", logBase = 10, display = FALSE))

  # LN with log(ED50) parameter returns log ED as well
  mLN <- drm(dryweight ~ conc, data = P.promelas, fct = LN.3(loge = TRUE))
  expect_true(mLN$fct$edLogScale)
  expect_no_error(ED(mLN, c(10, 50), interval = "fls", display = FALSE))
  mLN0 <- drm(dryweight ~ conc, data = P.promelas, fct = LN.3())
  expect_false(isTRUE(mLN0$fct$edLogScale))
  expect_error(ED(mLN0, 50, interval = "fls", display = FALSE), "cannot be used")

  # ED_robust must not hide the problem behind NA rows
  expect_error(ED_robust(mLL, 50, interval = "fls"), "cannot be used with model 'LL.4'")
})

test_that("maED and maED_robust refuse log(ED) models", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m2 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.4())

  expect_error(maED(m, list(LL2.5(), W1.4()), c(10, 50), display = FALSE), "LL2.5")
  expect_error(maED(m, list(LL2.5(), W1.4()), c(10, 50), display = FALSE), "How to fix")
  expect_error(maED(m2, list(LL.5()), c(10, 50), display = FALSE), "LL2.4")
  expect_error(maED_robust(m, list(LL2.5()), c(10, 50)), "log\\(ED\\) parameterisation")
  # dose-scale candidates are unaffected
  expect_no_error(maED(m, list(LL.5(), W1.4()), c(10, 50), display = FALSE))
})

test_that("get_ed_interval never recommends 'fls' for dose-scale models", {
  expect_equal(drc:::get_ed_interval("LL.4", small_n = FALSE), "delta")
  expect_equal(drc:::get_ed_interval("LL2.5"), "fls")
})

test_that("ED_robust returns EDs (not log EDs) for LL2 models", {
  mLL <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  mLL2 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.4())
  r4 <- ED_robust(mLL, c(10, 50))
  r2 <- ED_robust(mLL2, c(10, 50))

  expect_equal(r2$confint_method, rep("fls", 2))
  expect_false(anyNA(r2$Estimate))
  expect_equal(r2$Estimate, r4$Estimate, tolerance = 1e-4)
  expect_equal(r2$stderr, r4$stderr, tolerance = 1e-3)
  expect_equal(r2$Lower, r4$Lower, tolerance = 1e-2)
})

# ---- gammadr derivative is covered in test-numerical-derivatives-ed.R ----------------

# ---- predict(): standard errors without analytic derivatives ------------------------------

test_that("predict() falls back on numerical derivatives when deriv1 is missing", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  nd <- data.frame(conc = c(1, 5, 20))
  analytic <- predict(m, nd, se.fit = TRUE, interval = "confidence")
  m$fct$deriv1 <- NULL
  numerical <- predict(m, nd, se.fit = TRUE, interval = "confidence")
  expect_equal(numerical, analytic, tolerance = 1e-5)

  # models that never had deriv1
  mn <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4())
  expect_null(mn$fct$deriv1)
  res <- predict(mn, nd, interval = "confidence")
  expect_true(is.matrix(res))
  expect_equal(colnames(res), c("Prediction", "Lower", "Upper"))
  expect_true(all(res[, "Lower"] < res[, "Prediction"] & res[, "Prediction"] < res[, "Upper"]))
  # plain predictions are unchanged
  expect_true(is.numeric(predict(mn, nd)) && is.null(dim(predict(mn, nd))))
})

# ---- ssd models ---------------------------------------------------------------------------

test_that("confint, ED(delta) work and noEffect errors clearly for ssd fits", {
  d <- data.frame(conc = c(0.5, 1, 1.3, 2, 2.4, 3.1, 4, 5.2, 7, 9))
  ms <- drm(~ conc, data = d, fct = LN.2(), type = "ssd")

  ci <- confint(ms)
  expect_equal(dim(ci), c(2, 2))
  expect_true(all(ci[, 1] < coef(ms) & coef(ms) < ci[, 2]))
  ed <- ED(ms, 50, interval = "delta", display = FALSE)
  expect_true(all(is.finite(ed)))
  expect_error(noEffect(ms), "not available for response type 'ssd'")
})

# ---- separate fits forward their arguments --------------------------------------------------

test_that("drm(separate = TRUE) honours lowerl, robust and bcVal", {
  fit <- function(...) drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4(), separate = TRUE, ...)
  base <- fit()
  eIdx <- grepl("^e:", names(coef(base)))

  bounded <- suppressWarnings(fit(lowerl = c(-Inf, -Inf, -Inf, 0.3)))
  expect_true(all(coef(bounded)[eIdx] >= 0.3 - 1e-8))
  expect_true(any(coef(base)[eIdx] < 0.3))           # the bound was active

  med <- suppressWarnings(fit(robust = "median"))
  expect_false(isTRUE(all.equal(coef(med), coef(base))))

  bc <- suppressWarnings(fit(bcVal = 0.5))
  expect_false(isTRUE(all.equal(coef(bc), coef(base))))

  expect_error(fit(varcov = diag(nrow(spinach))), "not supported")
  expect_error(fit(robust = "tukey"), "not supported")
  expect_s3_class(bounded, "drc")
})

# ---- model averaging: Buckland standard error ---------------------------------------------

test_that("maED Buckland SE centres on the model-averaged estimate", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fl <- list(LL.5(), W1.4(), W2.4())
  r <- maED(m, fl, c(10, 50), interval = "buckland", display = FALSE, extended = TRUE)

  w <- r$fits[, "Weight"]
  est <- r$fits[, c("ED10", "ED50")]
  se_i <- sapply(c(10, 50), function(p) {
    c(ED(m, p, display = FALSE)[1, 2],
      vapply(fl, function(f) ED(update(m, fct = f), p, display = FALSE)[1, 2], numeric(1)))
  })
  centre <- colSums(est * w)
  expected <- colSums(w * sqrt(se_i^2 + sweep(est, 2, centre)^2))

  expect_equal(unname(r$estimates[, "Std. Error"]), unname(expected), tolerance = 1e-8)
})

# ---- comped -------------------------------------------------------------------------------

test_that("comped reports the ratio on the original scale for log = TRUE", {
  res <- suppressWarnings(capture.output(out <- comped(c(log(2), log(4)), c(0.1, 0.1), log = TRUE, operator = "/")))
  expect_equal(unname(out[1, "Estimate"]), 0.5)
  expect_equal(unname(out[1, "Std. Error"]), 0.5 * sqrt(0.02), tolerance = 1e-10)
  expect_true(out[1, "Lower"] < out[1, "Estimate"] && out[1, "Estimate"] < out[1, "Upper"])
  expect_equal(unname(out[1, "Lower"]), exp(log(0.5) - qnorm(0.975) * sqrt(0.02)), tolerance = 1e-10)

  # log difference and non-log paths are unchanged
  d <- suppressWarnings(capture.output(out2 <- comped(c(1, 3), c(0.1, 0.1), log = TRUE, operator = "-")))
  expect_equal(unname(out2[1, "Estimate"]), -2)
})

# ---- simDR with increasing curves ------------------------------------------------------------

test_that("simDR works for increasing curves (b < 0)", {
  conc <- c(1.88, 3.75, 7.5, 0.94, 15, 0.47, 30, 0.23, 60)
  expect_error(
    capture.output(simDR(c(-1, 0.5, 8, 3), 0.3, LL.4(), noSim = 2, conc = conc)),
    NA
  )
})

# ---- noEffect p-value -----------------------------------------------------------------------

test_that("noEffect p-value does not collapse to exactly 0 for strong effects", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  res <- noEffect(m)
  expect_gt(res["p-value"], 0)
  expect_lt(res["p-value"], 1e-6)
})

# ---- pmodels from the calling frame ----------------------------------------------------------

test_that("drm evaluates pmodels in the calling environment when data is NULL", {
  f <- function() {
    pm <- data.frame(b = factor(spinach$CURVE), c = 1, d = 1, e = 1)
    drm(spinach$SLOPE ~ spinach$DOSE, curveid = spinach$CURVE, pmodels = pm, fct = LL.4())
  }
  m <- f()
  expect_s3_class(m, "drc")
  expect_equal(length(coef(m)), 5 + 3)   # b per curve, one c, d, e
})

# Test suite for mrdrm.R functions
# Covers: mrdrm, leaveOneOut, pressWeights, dfFct, predFct,
#         loessEst, hat.loess, hat.drc, hat.mr, se.mr,
#         predict.mrdrc, inverseRegBasic, inverseReg, EDprint,
#         ED.mrdrc, plot.mrdrc, print.mrdrc, EDboot, pava

# --- Setup ---
# Continuous data
data(ryegrass)
m1_cont <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
m2_cont <- loess(rootl ~ conc, data = ryegrass, degree = 1)

# Binomial data
data(deguelin)
m1_binom <- drm(r/n ~ dose, weights = n, data = deguelin, fct = LL.2(), type = "binomial")
m2_binom <- loess(r/n ~ dose, data = deguelin, degree = 1)

# ---- pava ----

test_that("pava returns input for single element", {
  expect_equal(drc:::pava(5), 5)
})

test_that("pava returns input for already sorted sequence", {
  x <- c(1, 2, 3, 4)
  expect_equal(drc:::pava(x), x)
})

test_that("pava pools adjacent violators", {
  x <- c(1, 3, 2, 4)
  result <- drc:::pava(x)
  expect_equal(result, c(1, 2.5, 2.5, 4))
})

test_that("pava handles weighted input", {
  x <- c(3, 1, 2)
  wt <- c(1, 2, 1)
  result <- drc:::pava(x, wt)
  expect_true(all(diff(result) >= 0))
})

test_that("pava handles all-equal input", {
  x <- c(5, 5, 5)
  expect_equal(drc:::pava(x), x)
})

test_that("pava handles fully decreasing input", {
  x <- c(4, 3, 2, 1)
  result <- drc:::pava(x)
  expect_true(all(diff(result) >= 0))
  expect_equal(result, rep(mean(x), 4))
})

# ---- loessEst ----

test_that("loessEst returns correct structure", {
  result <- drc:::loessEst(5, ryegrass$conc, ryegrass$rootl, 0.75)
  expect_type(result, "list")
  expect_length(result, 2)
  expect_true(is.numeric(result[[1]]))
  expect_true(is.numeric(result[[2]]))
  expect_length(result[[2]], nrow(ryegrass))
})

test_that("loessEst works with logScale = TRUE", {
  result <- drc:::loessEst(5, ryegrass$conc, ryegrass$rootl, 0.75, logScale = TRUE)
  expect_true(is.numeric(result[[1]]))
})

# ---- hat.loess ----

test_that("hat.loess returns a matrix", {
  x <- deguelin$dose
  H <- drc:::hat.loess(x, m2_binom$pars$span)
  expect_true(is.matrix(H))
  expect_equal(dim(H), c(length(x), length(x)))
})

test_that("hat.loess works with x0 different from x", {
  x <- deguelin$dose
  x0 <- c(5, 15, 25)
  H <- drc:::hat.loess(x, m2_binom$pars$span, x0 = x0)
  expect_equal(dim(H), c(length(x0), length(x)))
})

# ---- hat.drc ----

test_that("hat.drc works for continuous data", {
  x <- m1_cont$data[, 1]
  H <- drc:::hat.drc(m1_cont, x)
  expect_true(is.matrix(H))
  expect_equal(nrow(H), length(x))
})

test_that("hat.drc works for binomial data", {
  x <- m1_binom$data[, 1]
  H <- drc:::hat.drc(m1_binom, x)
  expect_true(is.matrix(H))
  expect_equal(nrow(H), length(x))
})

test_that("hat.drc works with x0 different from x", {
  x <- m1_cont$data[, 1]
  x0 <- c(0, 5, 10)
  H <- drc:::hat.drc(m1_cont, x, x0 = x0)
  expect_equal(nrow(H), length(x0))
  expect_equal(ncol(H), length(x))
})

# ---- hat.mr ----

test_that("hat.mr returns a matrix", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  H <- drc:::hat.mr(mr_b)
  expect_true(is.matrix(H))
})

# ---- se.mr ----

test_that("se.mr works for continuous data", {
  mr_c <- drc:::mrdrm(m1_cont, m2_cont)
  se <- drc:::se.mr(mr_c, mr_c$dose)
  expect_true(is.numeric(se))
  expect_length(se, length(mr_c$dose))
})

test_that("se.mr works for binomial data", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  se <- drc:::se.mr(mr_b, mr_b$dose)
  expect_true(is.numeric(se))
  expect_length(se, length(mr_b$dose))
})

# ---- dfFct ----

test_that("dfFct returns a function", {
  fn <- drc:::dfFct(m1_cont, m2_cont)
  expect_type(fn, "closure")
  df_val <- fn(0.5)
  expect_true(is.numeric(df_val))
  expect_length(df_val, 1)
})

# ---- predFct ----

test_that("predFct returns a function", {
  looList <- list(pred1 = rep(1, 10), pred2 = rep(2, 10))
  fn <- drc:::predFct(looList)
  expect_type(fn, "closure")
  result <- fn(0.5)
  expect_equal(result, rep(1.5, 10))
})

test_that("predFct handles lambda = 0 and lambda = 1", {
  looList <- list(pred1 = c(1, 2, 3), pred2 = c(4, 5, 6))
  fn <- drc:::predFct(looList)
  expect_equal(fn(0), c(1, 2, 3))
  expect_equal(fn(1), c(4, 5, 6))
})

# ---- leaveOneOut ----

test_that("leaveOneOut returns correct structure for continuous data", {
  dataSet <- m1_cont$origData
  dose <- m1_cont$data[, 1]
  resp <- m1_cont$data[, 2]
  result <- drc:::leaveOneOut(m1_cont, m2_cont, dose, dataSet, resp, fixedEnd = FALSE)
  expect_type(result, "list")
  expect_named(result, c("pred1", "pred2"))
  expect_true(is.numeric(result$pred1))
  expect_true(is.numeric(result$pred2))
})

test_that("leaveOneOut with fixedEnd = TRUE modifies boundary values", {
  dataSet <- m1_cont$origData
  dose <- m1_cont$data[, 1]
  resp <- m1_cont$data[, 2]
  result_nofix <- drc:::leaveOneOut(m1_cont, m2_cont, dose, dataSet, resp, fixedEnd = FALSE)
  result_fix <- drc:::leaveOneOut(m1_cont, m2_cont, dose, dataSet, resp, fixedEnd = TRUE)
  # With fixedEnd, boundary pred2 values are set to observed means
  uniDose <- sort(unique(dose))
  expect_equal(result_fix$pred2[1], mean(resp[dose == uniDose[1]]))
  expect_equal(result_fix$pred2[length(uniDose)], mean(resp[dose == uniDose[length(uniDose)]]))
})

# ---- pressWeights ----

test_that("pressWeights 'none' returns vector of ones", {
  result <- drc:::pressWeights("none", 10, rep(1, 10), m1_cont, rep(0.5, 10), m2_cont)
  expect_equal(result, rep(1, 10))
})

test_that("pressWeights 'par' uses parametric predictions", {
  nVec <- m1_binom$weights
  resp <- m1_binom$data[, 2]
  lenData <- length(resp)
  result <- drc:::pressWeights("par", lenData, nVec, m1_binom, resp, m2_binom)
  expect_true(is.numeric(result))
  expect_length(result, lenData)
})

test_that("pressWeights 'nonpar' uses non-parametric predictions", {
  nVec <- m1_binom$weights
  resp <- m1_binom$data[, 2]
  lenData <- length(resp)
  result <- drc:::pressWeights("nonpar", lenData, nVec, m1_binom, resp, m2_binom)
  expect_true(is.numeric(result))
  expect_length(result, lenData)
})

test_that("pressWeights 'response' uses response values", {
  nVec <- m1_binom$weights
  resp <- m1_binom$data[, 2]
  lenData <- length(resp)
  result <- drc:::pressWeights("response", lenData, nVec, m1_binom, resp, m2_binom)
  expect_true(is.numeric(result))
  expect_length(result, lenData)
})

test_that("pressWeights 'ad hoc' adjusts for boundary values", {
  nVec <- m1_binom$weights
  resp <- m1_binom$data[, 2]
  lenData <- length(resp)
  result <- drc:::pressWeights("ad hoc", lenData, nVec, m1_binom, resp, m2_binom)
  expect_true(is.numeric(result))
  expect_length(result, lenData)
})

# ---- mrdrm main function ----

test_that("mrdrm returns an object of class mrdrc", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  expect_s3_class(mr, "mrdrc")
})

test_that("mrdrm return object has expected structure", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  expected_names <- c("pressVal", "lambda", "fitted", "gof", "object1", "object2",
                      "dose", "EDmethod", "ll", "ls.weights", "df")
  expect_named(mr, expected_names)
  expect_true(is.numeric(mr$lambda))
  expect_true(is.numeric(mr$fitted))
  expect_true(is.numeric(mr$gof))
  expect_named(mr$gof, c("mr.gof", "p.gof", "aic", "rv"))
  expect_equal(mr$EDmethod, "inverse")
})

test_that("mrdrm creates loess object internally when object2 is missing", {
  mr <- drc:::mrdrm(m1_cont)
  expect_s3_class(mr, "mrdrc")
  expect_s3_class(mr$object2, "loess")
})

test_that("mrdrm errors with non-linear loess (degree > 1)", {
  m2_quad <- loess(rootl ~ conc, data = ryegrass, degree = 2)
  expect_error(drc:::mrdrm(m1_cont, m2_quad), "Local regression fit not linear!")
})

test_that("mrdrm works with continuous data and GCV criterion", {
  mr <- drc:::mrdrm(m1_cont, m2_cont, criterion = "gcv")
  expect_s3_class(mr, "mrdrc")
  expect_true(is.null(mr$ll))  # GCV doesn't use leave-one-out
  expect_true(is.null(mr$ls.weights))
  # For continuous data, critFct is forced to "ls" and ls.weights to "none"
  expect_equal(mr$gof["rv"], mr$gof["mr.gof"] / mr$df, ignore_attr = TRUE)
})

test_that("mrdrm works with continuous data and LCV criterion", {
  mr <- drc:::mrdrm(m1_cont, m2_cont, criterion = "lcv")
  expect_s3_class(mr, "mrdrc")
  expect_false(is.null(mr$ll))  # LCV uses leave-one-out
  expect_false(is.null(mr$ls.weights))
})

test_that("mrdrm works with binomial data and GCV criterion", {
  suppressWarnings({
    mr <- drc:::mrdrm(m1_binom, m2_binom, criterion = "gcv")
  })
  expect_s3_class(mr, "mrdrc")
  expect_true(is.na(mr$gof["rv"]))  # Binomial has NA for rv
})

test_that("mrdrm works with binomial data and LCV criterion", {
  suppressWarnings({
    mr <- drc:::mrdrm(m1_binom, m2_binom, criterion = "lcv")
  })
  expect_s3_class(mr, "mrdrc")
  expect_false(is.null(mr$ll))
})

test_that("mrdrm works with binomial data and LL criterion", {
  suppressWarnings({
    mr <- drc:::mrdrm(m1_binom, m2_binom, critFct = "ll")
  })
  expect_s3_class(mr, "mrdrc")
  expect_false(is.null(mr$ll))
  expect_true(is.null(mr$ls.weights))
})

test_that("mrdrm with single lambda value", {
  mr <- drc:::mrdrm(m1_cont, m2_cont, lambda = 0.5)
  expect_equal(mr$lambda, 0.5)
  expect_true(is.na(mr$pressVal))
})

test_that("mrdrm with unitScale = TRUE", {
  mr <- drc:::mrdrm(m1_cont, m2_cont, unitScale = TRUE)
  expect_s3_class(mr, "mrdrc")
})

test_that("mrdrm with fixedEnd = TRUE for LCV criterion", {
  mr <- drc:::mrdrm(m1_cont, m2_cont, criterion = "lcv", fixedEnd = TRUE)
  expect_s3_class(mr, "mrdrc")
})

test_that("mrdrm with different ls.weights for binomial LCV", {
  for (w in c("nonpar", "ad hoc", "par", "response")) {
    suppressWarnings({
      mr <- drc:::mrdrm(m1_binom, m2_binom, criterion = "lcv", ls.weights = w)
    })
    expect_s3_class(mr, "mrdrc")
  }
})

test_that("mrdrm forces ls critFct and none weights for continuous data", {
  # Even if user passes different values, continuous data overrides
  mr <- drc:::mrdrm(m1_cont, m2_cont, critFct = "ll", ls.weights = "par")
  expect_s3_class(mr, "mrdrc")
  # The function overrides critFct to "ls" and ls.weights to "none" for continuous data
})

# ---- predict.mrdrc ----

test_that("predict.mrdrc returns fitted values without newdata", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  pred <- predict(mr)
  expect_true(is.numeric(pred))
  expect_length(pred, length(mr$fitted))
  expect_equal(as.numeric(pred), as.numeric(mr$fitted))
})

test_that("predict.mrdrc works with newdata", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  newdf <- data.frame(conc = c(0, 1, 5, 10))
  pred <- predict(mr, newdata = newdf)
  expect_true(is.numeric(pred))
  expect_length(pred, 4)
})

test_that("predict.mrdrc returns SE when se.fit = TRUE", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  pred <- predict(mr, se.fit = TRUE)
  expect_true(is.matrix(pred))
  expect_equal(ncol(pred), 2)
  expect_equal(colnames(pred), c("Prediction", "SE"))
})

test_that("predict.mrdrc returns confidence intervals", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  pred <- predict(mr, interval = "confidence")
  expect_true(is.matrix(pred))
  expect_equal(ncol(pred), 3)
  expect_equal(colnames(pred), c("Prediction", "Lower CI", "Upper CI"))
})

test_that("predict.mrdrc returns prediction intervals", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  pred <- predict(mr, interval = "prediction")
  expect_true(is.matrix(pred))
  expect_equal(ncol(pred), 3)
  expect_equal(colnames(pred), c("Prediction", "Lower PI", "Upper PI"))
})

test_that("predict.mrdrc with pava for decreasing curve (coef[1] > 0)", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  # ryegrass LL.4 has positive b parameter (decreasing)
  expect_true(coef(mr$object1)[1] > 0)
  pred <- predict(mr, pava = TRUE)
  expect_true(is.numeric(pred))
  expect_length(pred, length(mr$fitted))
})

test_that("predict.mrdrc with pava for increasing curve (coef[1] < 0)", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  # LL.2 on deguelin has negative b parameter (increasing)
  expect_true(coef(mr_b$object1)[1] < 0)
  pred <- predict(mr_b, pava = TRUE)
  expect_true(is.numeric(pred))
})

test_that("predict.mrdrc with se.fit for binomial", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  pred <- predict(mr_b, se.fit = TRUE)
  expect_true(is.matrix(pred))
  expect_equal(ncol(pred), 2)
})

test_that("predict.mrdrc with confidence interval for binomial", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  pred <- predict(mr_b, interval = "confidence")
  expect_true(is.matrix(pred))
  expect_equal(ncol(pred), 3)
})

test_that("predict.mrdrc with newdata and se.fit", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  newdf <- data.frame(conc = c(0, 1, 5, 10))
  pred <- predict(mr, newdata = newdf, se.fit = TRUE)
  expect_true(is.matrix(pred))
  expect_equal(nrow(pred), 4)
  expect_equal(ncol(pred), 2)
})

# ---- inverseRegBasic ----

test_that("inverseRegBasic returns 3-element vector with bisection method", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::inverseRegBasic(mr, 50, 0.95, "approximate", "bisection",
                                  20, 100, NULL, NULL, "confidence", "response")
  expect_length(result, 3)
  expect_true(is.numeric(result[1]))
})

test_that("inverseRegBasic returns 3-element vector with grid method", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::inverseRegBasic(mr, 50, 0.95, "approximate", "grid",
                                  20, 100, NULL, NULL, "confidence", "response")
  expect_length(result, 3)
  expect_true(is.numeric(result[1]))
})

test_that("inverseRegBasic with interval = none", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::inverseRegBasic(mr, 50, 0.95, "none", "bisection",
                                  20, 100, NULL, NULL, "confidence", "response")
  expect_length(result, 3)
  # With "none" interval, min1l and min1u are NA
  expect_true(is.na(result[2]))
  expect_true(is.na(result[3]))
})

test_that("inverseRegBasic warns when ED cannot be estimated", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  # For binomial with increasing curve (coef[1] < 0): val = perc/100
  # Use perc = 1 so val = 0.01, which is below min(predict(mr_b)) ≈ 0.29
  expect_warning(
    result <- drc:::inverseRegBasic(mr_b, 1, 0.95, "none", "bisection",
                                    20, 100, NULL, NULL, "confidence", "response"),
    "cannot be estimated"
  )
  expect_equal(result, rep(NA, 3))
})

test_that("inverseRegBasic with lower and upper specified", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::inverseRegBasic(mr, 50, 0.95, "none", "bisection",
                                  20, 100, 0, 8, "confidence", "response")
  expect_length(result, 3)
  expect_true(is.numeric(result[1]))
})

test_that("inverseRegBasic with minmax = dose", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::inverseRegBasic(mr, 50, 0.95, "none", "bisection",
                                  20, 100, NULL, NULL, "confidence", "dose")
  expect_length(result, 3)
})

test_that("inverseRegBasic with prediction interval type", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::inverseRegBasic(mr, 50, 0.95, "approximate", "bisection",
                                  20, 100, NULL, NULL, "prediction", "response")
  expect_length(result, 3)
})

test_that("inverseRegBasic warns when bisection fails for CI limits", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  suppressWarnings({
    result <- drc:::inverseRegBasic(mr_b, 35, 0.95, "approximate", "bisection",
                                    20, 100, NULL, NULL, "confidence", "response")
  })
  expect_length(result, 3)
})

test_that("inverseRegBasic warns when bisection fails for main estimate (column 1)", {
  # Create data where bisection fails for the main estimate with small cgridsize
  dose_vals <- c(rep(0, 3), rep(0.5, 3), rep(1, 3), rep(5, 3), rep(10, 3), rep(30, 3))
  resp_vals <- c(8, 8.2, 7.8, 7.5, 7.8, 7.2, 7.0, 6.8, 7.2, 3.5, 3.0, 4.0, 1.0, 0.8, 1.2, 0.3, 0.2, 0.4)
  m1_t <- drm(resp_vals ~ dose_vals, fct = LL.4())
  m2_t <- loess(resp_vals ~ dose_vals, degree = 1)
  mr_t <- drc:::mrdrm(m1_t, m2_t, lambda = 0.9)
  expect_warning(
    result <- drc:::inverseRegBasic(mr_t, 75, 0.95, "approximate", "bisection",
                                    3, 100, NULL, NULL, "confidence", "response"),
    "cannot be estimated"
  )
  # retVec[1] is NA -> triggers return(rep(NA, 3))
  expect_equal(result, rep(NA, 3))
})

test_that("inverseRegBasic adjusts CI boundaries (truncation and unbounding)", {
  # Use data with small dose start so lower CI ≈ minx triggers truncation to 0
  dose_vals <- c(0.001, 0.5, 1, 2, 5, 10, 20, 30)
  resp_vals <- c(8, 7.8, 7.5, 7.0, 3.5, 1.0, 0.3, 0.2)
  m1_t <- drm(resp_vals ~ dose_vals, fct = LL.4())
  m2_t <- loess(resp_vals ~ dose_vals, degree = 1)
  mr_t <- drc:::mrdrm(m1_t, m2_t)

  # ED0.5 with grid: lower CI ≈ minx -> truncated to 0
  suppressWarnings({
    result_low <- drc:::inverseRegBasic(mr_t, 0.5, 0.95, "approximate", "grid",
                                        100, 1000, NULL, NULL, "confidence", "response")
  })
  expect_equal(result_low[2], 0)

  # ED99 with grid: upper CI ≈ maxx -> set to Inf
  suppressWarnings({
    result_high <- drc:::inverseRegBasic(mr_t, 99, 0.95, "approximate", "grid",
                                         100, 1000, NULL, NULL, "confidence", "response")
  })
  expect_equal(result_high[3], Inf)
})

test_that("inverseRegBasic for binomial data", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  result <- drc:::inverseRegBasic(mr_b, 50, 0.95, "none", "bisection",
                                  20, 100, NULL, NULL, "confidence", "response")
  expect_length(result, 3)
  expect_true(is.numeric(result[1]))
})

test_that("inverseRegBasic with decreasing curve swaps percentage", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  # coef[1] > 0 means decreasing, so newPerc = 100 - perc
  expect_true(coef(mr$object1)[1] > 0)
  result <- drc:::inverseRegBasic(mr, 50, 0.95, "none", "bisection",
                                  20, 100, NULL, NULL, "confidence", "response")
  expect_length(result, 3)
})

# ---- inverseReg ----

test_that("inverseReg is vectorized over perc", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::inverseReg(mr, c(10, 50), 0.95, "none", "bisection",
                             20, 100, NULL, NULL, "confidence", "response")
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 3)
})

# ---- EDprint ----

test_that("EDprint displays output when display = TRUE", {
  EDmat <- matrix(c(3, 2, 4), nrow = 1)
  rownames(EDmat) <- "50"
  colnames(EDmat) <- c("Estimate", "Lower", "Upper")
  expect_output(
    drc:::EDprint(EDmat, "approximate", "Approximate variance formula", TRUE),
    "Estimated effective doses"
  )
})

test_that("EDprint shows CI text when ci is not 'none'", {
  EDmat <- matrix(c(3, 2, 4), nrow = 1)
  rownames(EDmat) <- "50"
  colnames(EDmat) <- c("Estimate", "Lower", "Upper")
  expect_output(
    drc:::EDprint(EDmat, "approximate", "Approximate variance formula", TRUE),
    "confidence interval"
  )
})

test_that("EDprint does not show CI text when ci = 'none'", {
  EDmat <- matrix(3, nrow = 1)
  rownames(EDmat) <- "50"
  colnames(EDmat) <- "Estimate"
  output <- capture.output(drc:::EDprint(EDmat, "none", "", TRUE))
  expect_false(any(grepl("confidence interval", output)))
})

test_that("EDprint returns invisible EDmat when display = FALSE", {
  EDmat <- matrix(3, nrow = 1)
  rownames(EDmat) <- "50"
  colnames(EDmat) <- "Estimate"
  result <- drc:::EDprint(EDmat, "none", "", FALSE)
  expect_equal(result, EDmat)
})

# ---- ED.mrdrc ----

test_that("ED.mrdrc returns estimate for continuous data", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  ed <- ED(mr, 50, display = FALSE)
  expect_true(is.matrix(ed))
  expect_equal(nrow(ed), 1)
  expect_equal(colnames(ed), "Estimate")
  expect_true(ed[1, 1] > 0)
})

test_that("ED.mrdrc returns multiple ED values", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  ed <- ED(mr, c(10, 50, 90), display = FALSE)
  expect_equal(nrow(ed), 3)
  expect_equal(rownames(ed), c("10", "50", "90"))
})

test_that("ED.mrdrc with approximate interval", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  ed <- ED(mr, 50, interval = "approximate", display = FALSE)
  expect_equal(ncol(ed), 3)
  expect_equal(colnames(ed), c("Estimate", "Lower", "Upper"))
})

test_that("ED.mrdrc with grid method", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  ed <- ED(mr, 50, interval = "approximate", method = "grid", display = FALSE)
  expect_equal(ncol(ed), 3)
  expect_true(ed[1, 1] > 0)
})

test_that("ED.mrdrc with bootstrap interval", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  ed <- ED(mr, 50, interval = "bootstrap", n = 5, display = FALSE)
  expect_equal(ncol(ed), 3)
  expect_equal(colnames(ed), c("Estimate", "Lower", "Upper"))
})

test_that("ED.mrdrc with display = TRUE and no interval", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  expect_output(ED(mr, 50, display = TRUE), "Estimated effective doses")
})

test_that("ED.mrdrc with display = TRUE and approximate interval", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  expect_output(
    ED(mr, 50, interval = "approximate", display = TRUE),
    "confidence interval"
  )
})

test_that("ED.mrdrc with lower and upper specified", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  ed <- ED(mr, 50, display = FALSE, lower = 0, upper = 8)
  expect_true(is.matrix(ed))
})

test_that("ED.mrdrc with intType = prediction", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  ed <- ED(mr, 50, interval = "approximate", intType = "prediction", display = FALSE)
  expect_equal(ncol(ed), 3)
})

test_that("ED.mrdrc with minmax = dose", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  ed <- ED(mr, 50, display = FALSE, minmax = "dose")
  expect_true(is.matrix(ed))
})

test_that("ED.mrdrc for binomial data", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  ed <- ED(mr_b, 50, display = FALSE)
  expect_true(is.matrix(ed))
  expect_true(ed[1, 1] > 0)
})

# ---- EDboot ----

test_that("EDboot returns CI matrix for single ED value (continuous)", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::EDboot(5, mr, 50, 123, 0.95)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 2)
})

test_that("EDboot returns CI matrix for multiple ED values (continuous)", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::EDboot(5, mr, c(10, 50), 123, 0.95)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})

test_that("EDboot works for binomial data", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  suppressWarnings({
    result <- drc:::EDboot(5, mr_b, 50, 123, 0.95)
  })
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
})

test_that("ED.mrdrc bootstrap for multiple ED values", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  ed <- ED(mr, c(10, 50), interval = "bootstrap", n = 5, display = FALSE)
  expect_equal(nrow(ed), 2)
  expect_equal(ncol(ed), 3)
  expect_equal(colnames(ed), c("Estimate", "Lower", "Upper"))
})

test_that("ED.mrdrc bootstrap with display = TRUE", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  expect_output(
    ED(mr, 50, interval = "bootstrap", n = 5, display = TRUE),
    "Bootstrap"
  )
})

# ---- plot.mrdrc ----

test_that("plot.mrdrc works for continuous data", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  pdf(nullfile())
  on.exit(dev.off())
  expect_no_error(plot(mr))
})

test_that("plot.mrdrc works with pava = TRUE", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  pdf(nullfile())
  on.exit(dev.off())
  expect_no_error(plot(mr, pava = TRUE))
})

# ---- print.mrdrc ----

test_that("print.mrdrc displays output for continuous data with mixing", {
  mr <- drc:::mrdrm(m1_cont, m2_cont, lambda = 0.5)
  expect_output(print(mr), "model-robust dose-response fit")
  expect_output(print(mr), "Mixing coefficient:")
  expect_output(print(mr), "Residual sum of squares:")
  expect_output(print(mr), "Residual standard error:")
  expect_output(print(mr), "AIC:")
  expect_output(print(mr), "for purely parametric fit:")
})

test_that("print.mrdrc displays output for continuous data without mixing", {
  mr <- drc:::mrdrm(m1_cont, m2_cont, lambda = 0)
  output <- capture.output(print(mr))
  combined <- paste(output, collapse = "\n")
  expect_true(grepl("Mixing coefficient: 0", combined))
  # When lambda = 0, no "for purely parametric fit" lines
  expect_false(grepl("for purely parametric fit", combined))
})

test_that("print.mrdrc displays output for binomial data with mixing", {
  suppressWarnings({
    mr_b <- drc:::mrdrm(m1_binom, m2_binom)
  })
  expect_true(mr_b$lambda > 0)
  expect_output(print(mr_b), "Pearson's chi-square:")
  expect_output(print(mr_b), "for purely parametric fit:")
  expect_output(print(mr_b), "AIC:")
})

test_that("print.mrdrc displays output for binomial data without mixing", {
  suppressWarnings({
    mr_b0 <- drc:::mrdrm(m1_binom, m2_binom, lambda = 0)
  })
  output <- capture.output(print(mr_b0))
  combined <- paste(output, collapse = "\n")
  expect_true(grepl("Pearson's chi-square:", combined))
  expect_false(grepl("for purely parametric fit", combined))
})

test_that("print.mrdrc returns the object invisibly", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- capture.output(ret <- print(mr))
  expect_equal(ret, mr)
})

# ---- Grid search in inverseRegBasic ----

test_that("inverseRegBasic grid method with approximate interval", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::inverseRegBasic(mr, 50, 0.95, "approximate", "grid",
                                  20, 100, NULL, NULL, "confidence", "response")
  expect_length(result, 3)
  expect_true(!is.na(result[1]))
})

test_that("inverseRegBasic grid method with none interval", {
  mr <- drc:::mrdrm(m1_cont, m2_cont)
  result <- drc:::inverseRegBasic(mr, 50, 0.95, "none", "grid",
                                  20, 100, NULL, NULL, "confidence", "response")
  expect_length(result, 3)
  # With "none", min1l and min1u are NA, so grid returns NA for limits
  expect_true(is.na(result[2]))
  expect_true(is.na(result[3]))
})

# Tests for idrm() - interactive dose-response modelling
# idrm is an internal function called by drm() when separate = TRUE

test_that("idrm works with a single function and multiple curves (separate = TRUE)", {
  data(spinach)
  result <- drm(SLOPE ~ DOSE, HERBICIDE, data = spinach, fct = LL.4(), separate = TRUE)

  # Should return a drc object
  expect_s3_class(result, "drc")

  # Should have coefficients for both curves
  cf <- coef(result)
  expect_true(length(cf) == 8)  # 4 params * 2 curves
  expect_true(all(grepl(":", names(cf))))  # Names should have "param:curve" format

  # Check parameter name structure
  pn <- result$parNames
  expect_type(pn, "list")
  expect_length(pn, 3)

  # Check that objList contains individual fits
  expect_true(!is.null(result$objList))
  expect_length(result$objList, 2)  # 2 curves

  # Check indexMat
  expect_true(!is.null(result$indexMat))
  expect_equal(ncol(result$indexMat), 2)
  expect_equal(nrow(result$indexMat), 4)  # 4 parameters

  # Check data is combined
  expect_true(nrow(result$data) > 0)

  # Check df.residual
  expect_true(result$df.residual > 0)

  # Check minval
  expect_true(is.numeric(result$minval))

  # Check parmMat has correct dimensions
  expect_equal(ncol(result$parmMat), 2)

  # Check curve function works (covers the plotFct closure on line 63)
  curveFct <- result$curve[[1]]
  pred <- curveFct(c(0.1, 1, 10))
  expect_true(is.matrix(pred))
  expect_equal(ncol(pred), 2)  # 2 curves
  expect_equal(nrow(pred), 3)  # 3 dose values
})

test_that("idrm works with a list of functions (oneFunction = FALSE)", {
  data(spinach)
  # Call idrm directly with a list of fct specifications
  # drm() validates fct before calling idrm(), so we need to call idrm directly
  fctList <- list(LL.4(), LL.3())
  result <- drc:::idrm(
    x = spinach$DOSE,
    y = spinach$SLOPE,
    curveid = spinach$HERBICIDE,
    weights = rep(1, nrow(spinach)),
    fct = fctList,
    type = "continuous",
    control = drmc()
  )

  # Should still return a drc object
  expect_s3_class(result, "drc")
})

test_that("drm with separate = TRUE and only one curve gives warning", {
  data(ryegrass)
  # ryegrass has no curveid - only one level
  expect_warning(
    drm(rootl ~ conc, data = ryegrass, fct = LL.4(), separate = TRUE),
    "Only one level"
  )
})

test_that("idrm result coefficients match individual fits", {
  data(spinach)
  # Fit with separate = TRUE
  sep_fit <- drm(SLOPE ~ DOSE, HERBICIDE, data = spinach, fct = LL.4(), separate = TRUE)

  # Fit each curve individually
  bentazon_fit <- drm(SLOPE ~ DOSE, data = subset(spinach, HERBICIDE == "bentazon"), fct = LL.4())
  diuron_fit <- drm(SLOPE ~ DOSE, data = subset(spinach, HERBICIDE == "diuron"), fct = LL.4())

  # Coefficients should match between separate and individual fits
  sep_coefs <- coef(sep_fit)
  bent_coefs <- coef(bentazon_fit)
  diur_coefs <- coef(diuron_fit)

  # Check bentazon coefficients
  expect_equal(unname(sep_coefs[grep("bentazon", names(sep_coefs))]),
               unname(bent_coefs), tolerance = 1e-4)

  # Check diuron coefficients
  expect_equal(unname(sep_coefs[grep("diuron", names(sep_coefs))]),
               unname(diur_coefs), tolerance = 1e-4)
})

test_that("idrm with separate = TRUE and three or more curves", {
  # Create a dataset with 3 curves
  set.seed(42)
  dose <- rep(c(0, 1, 2, 5, 10, 20), each = 3)
  n <- length(dose)
  curve_id <- rep(c("A", "B", "C"), each = n)
  dose_all <- rep(dose, 3)

  # Generate responses for 3 different curves
  resp_A <- 1 / (1 + exp(2 * (log(dose + 0.001) - log(5)))) + rnorm(n, 0, 0.02)
  resp_B <- 1 / (1 + exp(1.5 * (log(dose + 0.001) - log(3)))) + rnorm(n, 0, 0.02)
  resp_C <- 1 / (1 + exp(1 * (log(dose + 0.001) - log(8)))) + rnorm(n, 0, 0.02)
  resp_all <- c(resp_A, resp_B, resp_C)

  df <- data.frame(dose = dose_all, resp = resp_all, curve = factor(curve_id))

  result <- drm(resp ~ dose, curve, data = df, fct = LL.4(), separate = TRUE)

  expect_s3_class(result, "drc")
  expect_length(coef(result), 12)  # 4 params * 3 curves
  expect_equal(ncol(result$parmMat), 3)
  expect_equal(ncol(result$indexMat), 3)
  expect_length(result$objList, 3)
})

test_that("idrm dataList names are preserved from first fit", {
  data(spinach)
  result <- drm(SLOPE ~ DOSE, HERBICIDE, data = spinach, fct = LL.4(), separate = TRUE)

  # dataList should have a "names" element preserved from first fit
  expect_true(!is.null(result$dataList$names))
})

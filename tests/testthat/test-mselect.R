# Test mselect() function

# Create test dataset (ryegrass)
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

test_that("mselect returns non-zero Lack of fit p-values when nested = FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fctList <- list(LL.3(), W1.4(), W2.4())

  result <- mselect(m1, fctList = fctList, nested = FALSE)

  # "Lack of fit" column should contain meaningful p-values, not all zeros
  lof_col <- result[, "Lack of fit"]
  expect_true(any(lof_col > 0, na.rm = TRUE),
    info = "Lack of fit p-values should not all be zero when nested = FALSE")
})

test_that("mselect Lack of fit values match between nested = TRUE and nested = FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fctList <- list(LL.3(), W1.4(), W2.4())

  result_nested <- mselect(m1, fctList = fctList, nested = TRUE)
  result_plain <- mselect(m1, fctList = fctList, nested = FALSE)

  # The "Lack of fit" column should have the same values regardless of 'nested'
  expect_equal(
    result_nested[, "Lack of fit"],
    result_plain[, "Lack of fit"],
    info = "Lack of fit p-values should be the same whether nested is TRUE or FALSE"
  )
})

test_that("mselect with nested = FALSE returns correct column names", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fctList <- list(LL.3())

  result <- mselect(m1, fctList = fctList, nested = FALSE)
  expect_true("Lack of fit" %in% colnames(result))
  expect_true("logLik" %in% colnames(result))
  expect_true("IC" %in% colnames(result))
  expect_true("Res var" %in% colnames(result))
  expect_false("Nested F test" %in% colnames(result))
})

test_that("mselect with nested = TRUE returns Nested F test column", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fctList <- list(LL.3())

  result <- mselect(m1, fctList = fctList, nested = TRUE)
  expect_true("Nested F test" %in% colnames(result))
})

# Test backfit() function - Calculation of backfit values from fitted dose-response model

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

test_that("backfit returns correct structure", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- suppressWarnings(backfit(m1))

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("dose", "Estimate"))
  expect_null(rownames(result))
})

test_that("backfit returns correct number of rows matching unique dose levels", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- suppressWarnings(backfit(m1))

  unique_doses <- sort(unique(ryegrass$conc))
  expect_equal(nrow(result), length(unique_doses))
  expect_equal(result[, "dose"], unique_doses)
})

test_that("backfit values are numeric", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- suppressWarnings(backfit(m1))

  expect_true(is.numeric(result[, "Estimate"]))
})

test_that("backfit produces reasonable values for mid-range doses", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- suppressWarnings(backfit(m1))

  # For mid-range doses (within the dynamic range of the curve),

  # backfit values should approximate the original dose within a tolerance
  mid_idx <- which(result[, "dose"] == 3.75)
  expect_true(abs(result[mid_idx, "Estimate"] - 3.75) < 2)
})

test_that("backfit works with different dose-response models", {
  m_w1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
  result <- suppressWarnings(backfit(m_w1))

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("dose", "Estimate"))
})

# Test file for gaussian.ssf function
# Note: gaussian.ssf is an internal function, accessed via :::

# Test basic functionality and correctness

test_that("gaussian.ssf returns a function", {
  result <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  expect_type(result, "closure")
  expect_true(is.function(result))
})

test_that("gaussian.ssf with method 1 returns valid initial values", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
  expect_true(all(is.finite(init_vals)))
})

test_that("gaussian.ssf with different methods returns valid initial values", {
  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  for (m in c("1", "2", "3", "4")) {
    ssfct <- drc:::gaussian.ssf(method = m, fixed = c(NA, NA, NA, NA, NA))
    init_vals <- ssfct(dframe)
    expect_type(init_vals, "double")
    expect_length(init_vals, 5)
  }
})

test_that("gaussian.ssf f parameter is always 1", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  # The 5th parameter (f) should be 1
  expect_equal(init_vals[5], 1)
})

test_that("gaussian.ssf e parameter is x at max y", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  # The 4th parameter (e) should be x[which.max(y)] = 5
  expect_equal(init_vals[4], 5)
})

test_that("gaussian.ssf c and d parameters use findcd", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  x <- c(1, 2, 3, 5, 7, 10)
  y <- c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  dframe <- data.frame(dose = x, response = y)

  init_vals <- ssfct(dframe)

  # findcd(x, y) returns c(min(y) - 0.001*diff(range(y)), max(y) + 0.001*diff(range(y)))
  expected_cd <- drc:::findcd(x, y)
  expect_equal(init_vals[2], expected_cd[1])
  expect_equal(init_vals[3], expected_cd[2])
})

# Test logg parameter

test_that("gaussian.ssf with logg=FALSE uses sd(x[...])", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA), logg = FALSE)

  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  # b should be 0.75 * sd(x[y > quantile(y, .75)])
  x <- dframe[, 1]
  y <- dframe[, 2]
  expected_b <- 0.75 * sd(x[y > quantile(y, .75)])
  expect_equal(init_vals[1], expected_b)
})

test_that("gaussian.ssf with logg=TRUE uses sd(log(x[...]))", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA), logg = TRUE)

  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  # b should be 0.75 * sd(log(x[y > quantile(y, .75)]))
  x <- dframe[, 1]
  y <- dframe[, 2]
  expected_b <- 0.75 * sd(log(x[y > quantile(y, .75)]))
  expect_equal(init_vals[1], expected_b)
})

# Test useFixed parameter

test_that("gaussian.ssf with useFixed=TRUE does not error", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA), useFixed = TRUE)

  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

test_that("gaussian.ssf with useFixed=FALSE works normally", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA), useFixed = FALSE)

  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

# Test fixed parameters

test_that("gaussian.ssf works with fixed parameters", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(1, NA, NA, NA, 1))

  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  # Should only return values for non-fixed parameters (c, d, e)
  expect_length(init_vals, 3)
})

test_that("gaussian.ssf works with all parameters fixed", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(1, 0, 1, 5, 1))

  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  expect_length(init_vals, 0)
})

test_that("gaussian.ssf works with different fixed parameter combinations", {
  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  # Fix only b
  ssfct1 <- drc:::gaussian.ssf(method = "1", fixed = c(1, NA, NA, NA, NA))
  init1 <- ssfct1(dframe)
  expect_length(init1, 4)

  # Fix c and d
  ssfct2 <- drc:::gaussian.ssf(method = "1", fixed = c(NA, 0, 1, NA, NA))
  init2 <- ssfct2(dframe)
  expect_length(init2, 3)

  # Fix e
  ssfct3 <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, 5, NA))
  init3 <- ssfct3(dframe)
  expect_length(init3, 4)
})

# Test error handling

test_that("gaussian.ssf errors with invalid method", {
  expect_error(
    drc:::gaussian.ssf(method = "5", fixed = c(NA, NA, NA, NA, NA)),
    "'arg' should be one of"
  )
})

test_that("gaussian.ssf errors with non-character method", {
  expect_error(
    drc:::gaussian.ssf(method = 1, fixed = c(NA, NA, NA, NA, NA))
  )
})

# Test edge cases

test_that("gaussian.ssf works with minimal data", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(1, 5, 10),
    response = c(0.1, 1.0, 0.1)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

test_that("gaussian.ssf returns correct parameter order", {
  ssfct <- drc:::gaussian.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(1, 2, 3, 5, 7, 10),
    response = c(0.1, 0.5, 0.9, 1.0, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  # Order should be: b, c, d, e, f
  # f should be 1
  expect_equal(init_vals[5], 1)
  # e should be x[which.max(y)] = 5
  expect_equal(init_vals[4], 5)
})

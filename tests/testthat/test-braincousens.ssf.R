# Test file for braincousens.ssf function
# Note: braincousens.ssf is an internal function, accessed via :::

# Test basic functionality and correctness

test_that("braincousens.ssf returns a function", {
  result <- drc:::braincousens.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  expect_type(result, "closure")
  expect_true(is.function(result))
})

test_that("braincousens.ssf with method 1 returns valid initial values", {
  ssfct <- drc:::braincousens.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  # Create test data frame
  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

test_that("braincousens.ssf with method 2 returns valid initial values", {
  ssfct <- drc:::braincousens.ssf(method = "2", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

test_that("braincousens.ssf with method 3 returns valid initial values", {
  ssfct <- drc:::braincousens.ssf(method = "3", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

test_that("braincousens.ssf with method 4 returns valid initial values", {
  ssfct <- drc:::braincousens.ssf(method = "4", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

test_that("braincousens.ssf works with fixed parameters", {
  # Fix first and last parameters
  ssfct <- drc:::braincousens.ssf(method = "1", fixed = c(1, NA, NA, NA, 0))

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- ssfct(dframe)

  # Should only return values for non-fixed parameters (c, d, e)
  expect_length(init_vals, 3)
})

test_that("braincousens.ssf works with all parameters fixed", {
  ssfct <- drc:::braincousens.ssf(method = "1", fixed = c(1, 0, 1, 1, 0))

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- ssfct(dframe)

  # Should return empty vector since all parameters are fixed
  expect_length(init_vals, 0)
})

test_that("braincousens.ssf f parameter is always 0", {
  ssfct <- drc:::braincousens.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- ssfct(dframe)

  # The 5th parameter (f) should be 0
  expect_equal(init_vals[5], 0)
})

test_that("braincousens.ssf works with useFixed parameter", {
  # Test with useFixed = TRUE (although not implemented, should not error)
  ssfct <- drc:::braincousens.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA), useFixed = TRUE)

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

test_that("braincousens.ssf works with useFixed = FALSE", {
  ssfct <- drc:::braincousens.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA), useFixed = FALSE)

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

# Test error handling

test_that("braincousens.ssf errors with invalid method", {
  expect_error(
    drc:::braincousens.ssf(method = "5", fixed = c(NA, NA, NA, NA, NA)),
    "'arg' should be one of"
  )
})

test_that("braincousens.ssf errors with non-character method", {
  expect_error(
    drc:::braincousens.ssf(method = 1, fixed = c(NA, NA, NA, NA, NA))
  )
})

# Test edge cases

test_that("braincousens.ssf works with minimal data", {
  ssfct <- drc:::braincousens.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  # Minimal data frame
  dframe <- data.frame(
    dose = c(0.1, 1, 10),
    response = c(0.9, 0.5, 0.1)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

test_that("braincousens.ssf works with larger dataset", {
  ssfct <- drc:::braincousens.ssf(method = "2", fixed = c(NA, NA, NA, NA, NA))

  # Larger data frame
  dframe <- data.frame(
    dose = rep(c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20), each = 3),
    response = c(
      rep(0.95, 3), rep(0.9, 3), rep(0.7, 3), rep(0.5, 3),
      rep(0.3, 3), rep(0.15, 3), rep(0.05, 3), rep(0.02, 3)
    )
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

test_that("braincousens.ssf handles data with varying responses", {
  ssfct <- drc:::braincousens.ssf(method = "3", fixed = c(NA, NA, NA, NA, NA))

  # Data with some variation
  set.seed(123)
  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05) + rnorm(6, 0, 0.05)
  )

  init_vals <- ssfct(dframe)

  expect_type(init_vals, "double")
  expect_length(init_vals, 5)
})

test_that("braincousens.ssf method 1 uses findbe1", {
  # Test that method 1 produces expected behavior
  ssfct <- drc:::braincousens.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 1, 10, 100),
    response = c(0.95, 0.75, 0.25, 0.05)
  )

  init_vals <- ssfct(dframe)

  # All values should be finite
  expect_true(all(is.finite(init_vals)))
})

test_that("braincousens.ssf method 2 uses findbe2 with Anke", {
  ssfct <- drc:::braincousens.ssf(method = "2", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 1, 10, 100),
    response = c(0.95, 0.75, 0.25, 0.05)
  )

  init_vals <- ssfct(dframe)

  expect_true(all(is.finite(init_vals)))
})

test_that("braincousens.ssf method 3 uses findbe3", {
  ssfct <- drc:::braincousens.ssf(method = "3", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 1, 10, 100),
    response = c(0.95, 0.75, 0.25, 0.05)
  )

  init_vals <- ssfct(dframe)

  expect_true(all(is.finite(init_vals)))
})

test_that("braincousens.ssf method 4 uses findbe2 with Normolle", {
  ssfct <- drc:::braincousens.ssf(method = "4", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 1, 10, 100),
    response = c(0.95, 0.75, 0.25, 0.05)
  )

  init_vals <- ssfct(dframe)

  expect_true(all(is.finite(init_vals)))
})

test_that("braincousens.ssf works with different fixed parameter combinations", {
  # Fix only b parameter
  ssfct1 <- drc:::braincousens.ssf(method = "1", fixed = c(1, NA, NA, NA, NA))

  # Fix c and d parameters
  ssfct2 <- drc:::braincousens.ssf(method = "1", fixed = c(NA, 0, 1, NA, NA))

  # Fix e parameter
  ssfct3 <- drc:::braincousens.ssf(method = "1", fixed = c(NA, NA, NA, 1, NA))

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init1 <- ssfct1(dframe)
  init2 <- ssfct2(dframe)
  init3 <- ssfct3(dframe)

  expect_length(init1, 4)
  expect_length(init2, 3)
  expect_length(init3, 4)
})

test_that("braincousens.ssf helper functions are defined", {
  # The helper functions should be created inside braincousens.ssf
  ssfct <- drc:::braincousens.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  # Should work without error
  dframe <- data.frame(
    dose = c(0.1, 1, 10),
    response = c(0.9, 0.5, 0.1)
  )

  expect_no_error(ssfct(dframe))
})

test_that("braincousens.ssf returns correct parameter order", {
  ssfct <- drc:::braincousens.ssf(method = "1", fixed = c(NA, NA, NA, NA, NA))

  dframe <- data.frame(
    dose = c(0.1, 0.5, 1, 2, 5, 10),
    response = c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05)
  )

  init_vals <- ssfct(dframe)

  # Should return c(b, c, d, e, f) in that order
  # We know f should be 0
  expect_equal(init_vals[5], 0)
})

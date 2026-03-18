# Test file for iband function
# Note: iband is an internal function, accessed via :::

test_that("iband returns NULL", {
  result <- drc:::iband(object = NULL)

  expect_null(result)
})

test_that("iband returns NULL for various input types", {
  expect_null(drc:::iband(object = 1))
  expect_null(drc:::iband(object = "test"))
  expect_null(drc:::iband(object = list()))
})

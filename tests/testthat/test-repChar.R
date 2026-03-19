# Tests for repChar (R/repChar.R) and its nested helper functions
# repChar is an internal function used in 'mixdrc' for replacing characters
# in strings and building function/formula strings.

# -------------------------------------------------------------------
# Tests for repChar
# -------------------------------------------------------------------

test_that("repChar works with basic inputs and no fixed values (NULL fixed)", {

  # When fixed is NULL, it should be replaced with rep(NA, length(names))
  # This tests the is.null(fixed) == TRUE branch
  result <- drc:::repChar(
    str = "b + c * DOSE",
    names = c("b", "c"),
    fixed = NULL,
    keep = c("DOSE")
  )

  expect_true(is.list(result))
  expect_length(result, 2)

  # First element is the function string
  expect_true(is.character(result[[1]]))
  # Second element is the formula string
  expect_true(is.character(result[[2]]))

  # Since fixed is NULL (all NA), all names should appear in argNames
  expect_true(grepl("b", result[[1]]))
  expect_true(grepl("c", result[[1]]))
  expect_true(grepl("DOSE", result[[1]]))

  # The formula string should contain all unfixed argument names
  expect_true(grepl("b,c", result[[2]]))
})

test_that("repChar works with some fixed values (fixed is not NULL)", {
  # This tests the is.null(fixed) == FALSE branch
  # and the !is.na(fixed[i]) == TRUE branch in the inner loop

  result <- drc:::repChar(
    str = "b + c * DOSE",
    names = c("b", "c"),
    fixed = c(5, NA),
    keep = c("DOSE")
  )

  expect_true(is.list(result))
  expect_length(result, 2)

  # "b" should be replaced by its fixed value "5" in the body string
  expect_true(grepl("5", result[[1]]))
  # "c" is unfixed (NA), so it should remain in the function string
  expect_true(grepl("c", result[[1]]))
  # DOSE should be preserved (kept)
  expect_true(grepl("DOSE", result[[1]]))

  # The formula string should only contain unfixed argument names
  # Only "c" is unfixed
  expect_true(grepl("c", result[[2]]))
})

test_that("repChar works when all values are fixed", {
  # All names are fixed, no unfixed names should appear in argNames
  result <- drc:::repChar(
    str = "b + c * DOSE",
    names = c("b", "c"),
    fixed = c(5, 10),
    keep = c("DOSE")
  )

  expect_true(is.list(result))
  expect_length(result, 2)

  # Both b and c should be substituted with their fixed values
  expect_true(grepl("5", result[[1]]))
  expect_true(grepl("10", result[[1]]))
  # DOSE should be preserved
  expect_true(grepl("DOSE", result[[1]]))
})

test_that("repChar handles multiple keep patterns", {
  # Tests that the for-loop over keep works with multiple elements
  # The sep parameter defaults to c(",", ";") so two keep patterns use "," and ";"
  result <- drc:::repChar(
    str = "b * DOSE + c * CONC",
    names = c("b", "c"),
    fixed = c(NA, 3),
    keep = c("DOSE", "CONC")
  )

  expect_true(is.list(result))
  expect_length(result, 2)

  # "c" should be replaced by "3"
  expect_true(grepl("3", result[[1]]))
  # "b" should remain (unfixed)
  expect_true(grepl("b", result[[1]]))
  # Both DOSE and CONC should be preserved
  expect_true(grepl("DOSE", result[[1]]))
  expect_true(grepl("CONC", result[[1]]))
})

test_that("repChar builds correct function header string", {
  result <- drc:::repChar(
    str = "b + c * DOSE",
    names = c("b", "c"),
    fixed = c(NA, NA),
    keep = c("DOSE")
  )

  # The function string should have the form:
  # "function(DOSE, b,c){( b + c * DOSE ^lambda - 1)/lambda}"
  expect_true(grepl("function\\(DOSE,", result[[1]]))
  expect_true(grepl("lambda", result[[1]]))
})

test_that("repChar builds correct formula string", {
  result <- drc:::repChar(
    str = "b + c * DOSE",
    names = c("b", "c"),
    fixed = c(NA, NA),
    keep = c("DOSE")
  )

  # The formula string should have the form:
  # "formula(respVar ~ opfct(doseVar, b,c))"
  expect_true(grepl("formula\\(respVar", result[[2]]))
  expect_true(grepl("opfct\\(doseVar,", result[[2]]))
  expect_true(grepl("b,c", result[[2]]))
})

test_that("repChar with single name and single keep", {
  result <- drc:::repChar(
    str = "a * DOSE",
    names = c("a"),
    fixed = c(NA),
    keep = c("DOSE")
  )

  expect_true(is.list(result))
  expect_length(result, 2)
  expect_true(grepl("a", result[[1]]))
  expect_true(grepl("DOSE", result[[1]]))
  expect_true(grepl("a", result[[2]]))
})

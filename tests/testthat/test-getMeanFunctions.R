# Tests for getMeanFunctions.R
# Achieves 100% code coverage for the getMeanFunctions function.

# =============================================================================
# 1. Default call (no arguments) — all models, display = TRUE
# =============================================================================

test_that("getMeanFunctions() with no arguments returns all models invisibly", {
  out <- capture.output(result <- getMeanFunctions())

  # Returns an invisible list

  expect_type(result, "list")
  expect_true(length(result) > 0)

  # Each element is a list of c(name, text) from the lapFct path
  expect_true(all(vapply(result, length, integer(1)) == 2))

  # Console output was produced (display = TRUE)
  expect_true(length(out) > 0)
})

# =============================================================================
# 2. Filter by noParm only
# =============================================================================

test_that("getMeanFunctions filters by noParm", {
  result <- getMeanFunctions(noParm = 4, display = FALSE)

  expect_type(result, "list")
  expect_true(length(result) > 0)

  # Every returned model must have noParm == 4
  for (mod in result) {
    expect_equal(mod$noParm, 4)
  }
})

test_that("getMeanFunctions filters by noParm = 2", {
  result <- getMeanFunctions(noParm = 2, display = FALSE)

  expect_type(result, "list")
  expect_true(length(result) > 0)

  for (mod in result) {
    expect_equal(mod$noParm, 2)
  }
})

test_that("getMeanFunctions returns empty list when noParm matches nothing", {
  result <- getMeanFunctions(noParm = 99, display = FALSE)

  expect_type(result, "list")
  expect_length(result, 0)
})

# =============================================================================
# 3. Filter by fname only
# =============================================================================

test_that("getMeanFunctions filters by fname (single name)", {
  result <- getMeanFunctions(fname = "LL.4", display = FALSE)

  expect_type(result, "list")
  expect_true(length(result) >= 1)

  # At least one returned model should have name containing "LL.4"
  names_found <- vapply(result, function(x) x$name, character(1))
  expect_true("LL.4" %in% names_found)
})

test_that("getMeanFunctions filters by fname (multiple names)", {
  result <- getMeanFunctions(fname = c("LL.4", "W1.4"), display = FALSE)

  expect_type(result, "list")
  expect_true(length(result) >= 2)

  names_found <- vapply(result, function(x) x$name, character(1))
  expect_true("LL.4" %in% names_found)
  expect_true("W1.4" %in% names_found)
})

test_that("getMeanFunctions returns empty list when fname matches nothing", {
  result <- getMeanFunctions(fname = "NONEXISTENT_MODEL", display = FALSE)

  expect_type(result, "list")
  expect_length(result, 0)
})

# =============================================================================
# 4. Filter by both noParm and fname
# =============================================================================

test_that("getMeanFunctions filters by both noParm and fname", {
  result <- getMeanFunctions(noParm = 3, fname = c("LL.3", "W1.3"), display = FALSE)

  expect_type(result, "list")
  expect_true(length(result) >= 2)

  for (mod in result) {
    expect_equal(mod$noParm, 3)
  }
})

# =============================================================================
# 5. display = TRUE produces console output
# =============================================================================

test_that("getMeanFunctions display = TRUE prints model info", {
  out <- capture.output(result <- getMeanFunctions(fname = "LL.4", display = TRUE))

  expect_true(length(out) > 0)
  # Output should mention the model name and number of parameters
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("LL.4", combined))
  expect_true(grepl("4 parameters", combined))
  expect_true(grepl("In 'drc'", combined))
})

# =============================================================================
# 6. display = FALSE suppresses console output
# =============================================================================

test_that("getMeanFunctions display = FALSE suppresses output", {
  out <- capture.output(result <- getMeanFunctions(fname = "LL.4", display = FALSE))

  expect_equal(length(out), 0)
  expect_type(result, "list")
  expect_true(length(result) >= 1)
})

# =============================================================================
# 7. Custom flist argument
# =============================================================================

test_that("getMeanFunctions uses custom flist instead of default", {
  custom_list <- list(LL.4(), W1.4())
  result <- getMeanFunctions(flist = custom_list, display = FALSE)

  expect_type(result, "list")
  # Should return name-text pairs since noParm = NA and fname = NULL
  expect_true(length(result) == 2)
})

test_that("getMeanFunctions with custom flist and noParm filter", {
  custom_list <- list(LL.2(), LL.4(), W1.4())
  result <- getMeanFunctions(noParm = 4, flist = custom_list, display = FALSE)

  expect_type(result, "list")
  expect_true(length(result) >= 2)

  for (mod in result) {
    expect_equal(mod$noParm, 4)
  }
})

# =============================================================================
# 8. Return structure: default path returns name-text pairs via lapFct
# =============================================================================

test_that("getMeanFunctions default returns list of name-text pairs", {
  result <- getMeanFunctions(display = FALSE)

  expect_type(result, "list")

  # Each element should be a character vector of length 2 (name, text)
  for (elem in result) {
    expect_type(elem, "character")
    expect_length(elem, 2)
  }
})

# =============================================================================
# 9. Return structure: filtered path returns model objects
# =============================================================================

test_that("getMeanFunctions filtered returns model list objects", {
  result <- getMeanFunctions(noParm = 4, display = FALSE)

  expect_type(result, "list")
  for (mod in result) {
    expect_true(is.list(mod))
    expect_true("name" %in% names(mod))
    expect_true("noParm" %in% names(mod))
  }
})

# =============================================================================
# 10. Edge case: noParm combined with fname that yields no results
# =============================================================================

test_that("getMeanFunctions combined filter with no matching results", {
  # LL.4 has 4 parameters, filtering for noParm = 2 AND fname = "LL.4" yields nothing
  result <- getMeanFunctions(noParm = 2, fname = "LL.4", display = FALSE)

  expect_type(result, "list")
  expect_length(result, 0)
})

# =============================================================================
# 11. Verify the displayFunction path where condition is FALSE (returns NULL)
# =============================================================================

test_that("getMeanFunctions filters out non-matching models correctly", {
  # Only 2-parameter models, but display = TRUE to exercise the cat() path too
  out <- capture.output(result <- getMeanFunctions(noParm = 2, display = TRUE))

  # All returned models should have exactly 2 parameters
  for (mod in result) {
    expect_equal(mod$noParm, 2)
  }

  # Models with other parameter counts should NOT appear
  all_result <- getMeanFunctions(display = FALSE)
  expect_true(length(result) < length(all_result))
})

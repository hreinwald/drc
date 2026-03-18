# Tests for fct2list.R: vec2mat, nParm, fParm, fct2list

# =============================================================================
# Helper: define simple test functions
# =============================================================================

# A simple 2-parameter function using vector indexing
simple_2p <- function(x, b) {
  b[1] + b[2] * x
}

# A 3-parameter function using vector indexing
simple_3p <- function(x, b) {
  b[1] + b[2] * x + b[3] * x^2
}

# A 4-parameter log-logistic function using vector indexing
ll4_raw <- function(dose, b) {
  b[2] + (b[3] - b[2]) / (1 + exp(b[1] * (log(dose) - log(b[4]))))
}

# A single-line function (no braces in body)
single_line_fct <- function(x, b) b[1] + b[2] * x

# =============================================================================
# Tests for vec2mat()
# =============================================================================

test_that("vec2mat converts vector indexing to matrix indexing for a 2-param function", {
  result <- drc:::vec2mat(simple_2p, 2)

  # Returns a list of length 3

  expect_type(result, "list")
  expect_length(result, 3)

  # First element is a function

  expect_true(is.function(result[[1]]))

  # Second element is a character string (body text)
  expect_type(result[[2]], "character")

  # Third element is the parameter name
  expect_equal(result[[3]], "b")

  # The body string should contain matrix-style indexing [,
  expect_true(grepl("\\[,", result[[2]]))
})

test_that("vec2mat works with a 3-param function", {
  result <- drc:::vec2mat(simple_3p, 2)

  expect_type(result, "list")
  expect_length(result, 3)
  expect_true(is.function(result[[1]]))
  expect_equal(result[[3]], "b")
  expect_true(grepl("\\[,", result[[2]]))
})

test_that("vec2mat works with a 4-param log-logistic function", {
  result <- drc:::vec2mat(ll4_raw, 2)

  expect_type(result, "list")
  expect_length(result, 3)
  expect_true(is.function(result[[1]]))
  expect_equal(result[[3]], "b")
})

test_that("vec2mat works with a single-line (no braces) function body", {
  result <- drc:::vec2mat(single_line_fct, 2)

  expect_type(result, "list")
  expect_length(result, 3)
  expect_true(is.function(result[[1]]))
  expect_equal(result[[3]], "b")
})

test_that("vec2mat errors when argument number does not exist", {
  expect_error(
    drc:::vec2mat(simple_2p, 5),
    "Argument number does not exist"
  )
})

test_that("vec2mat works with first argument index", {
  f <- function(a, b) a[1] + a[2] * b
  result <- drc:::vec2mat(f, 1)

  expect_type(result, "list")
  expect_equal(result[[3]], "a")
  expect_true(grepl("\\[,", result[[2]]))
})

# =============================================================================
# Tests for nParm()
# =============================================================================

test_that("nParm counts 2 unique parameters correctly", {
  # vec2mat result for 2-param function
  v2m <- drc:::vec2mat(simple_2p, 2)
  result <- drc:::nParm(v2m[[2]])
  expect_equal(result, 2)
})

test_that("nParm counts 3 unique parameters correctly", {
  v2m <- drc:::vec2mat(simple_3p, 2)
  result <- drc:::nParm(v2m[[2]])
  expect_equal(result, 3)
})

test_that("nParm counts 4 unique parameters correctly", {
  v2m <- drc:::vec2mat(ll4_raw, 2)
  result <- drc:::nParm(v2m[[2]])
  expect_equal(result, 4)
})

# =============================================================================
# Tests for fParm()
# =============================================================================

test_that("fParm returns the vec2mat function when all fixed are NA", {
  # When all elements of fixed are NA, fParm returns v2m[[1]] early
  fixed_all_na <- c(NA, NA)
  result <- drc:::fParm(simple_2p, 2, fixed_all_na)

  expect_true(is.function(result))
})

test_that("fParm fixes one parameter and frees another", {
  # Fix b[1] = 5, let b[2] be free
  fixed <- c(5, NA)
  result <- drc:::fParm(simple_2p, 2, fixed)

  expect_true(is.function(result))
})

test_that("fParm fixes all parameters", {
  # Fix both parameters
  fixed <- c(5, 3)
  result <- drc:::fParm(simple_2p, 2, fixed)

  expect_true(is.function(result))
})

test_that("fParm works with 3-param function fixing one parameter", {
  # Fix b[2] = 10, leave b[1] and b[3] free
  fixed <- c(NA, 10, NA)
  result <- drc:::fParm(simple_3p, 2, fixed)

  expect_true(is.function(result))
})

test_that("fParm works with 4-param function fixing multiple parameters", {
  # Fix b[1] = 1, b[3] = 100, free b[2] and b[4]
  fixed <- c(1, NA, 100, NA)
  result <- drc:::fParm(ll4_raw, 2, fixed)

  expect_true(is.function(result))
})

test_that("fParm handles renumbering when parameter count changes digit length", {
  # Create a function with many parameters (>= 10) to trigger the digit-count branch
  # This tests lines 106-112 where nchar(inStr3) < nchar(as.character(numVec[i]))
  many_params_fct <- function(x, b) {
    b[1] + b[2] * x + b[3] * x^2 + b[4] * x^3 + b[5] * x^4 +
      b[6] * x^5 + b[7] * x^6 + b[8] * x^7 + b[9] * x^8 + b[10] * x^9
  }
  # Fix parameters 1-9, leave parameter 10 free
  # This means parameter 10 gets renumbered to 1 (nchar("1") < nchar("10"))
  fixed <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, NA)
  result <- drc:::fParm(many_params_fct, 2, fixed)
  expect_true(is.function(result))
})

# =============================================================================
# Tests for fct2list()
# =============================================================================

test_that("fct2list returns a proper list for a 2-param function", {
  result <- drc:::fct2list(simple_2p, 2)

  expect_type(result, "list")
  expect_length(result, 3)

  # First element is a function
  expect_true(is.function(result[[1]]))

  # Second element is NULL
  expect_null(result[[2]])

  # Third element is parameter names as letters
  expect_equal(result[[3]], c("a", "b"))
})

test_that("fct2list returns correct parameter names for a 3-param function", {
  result <- drc:::fct2list(simple_3p, 2)

  expect_equal(result[[3]], c("a", "b", "c"))
  expect_null(result[[2]])
  expect_true(is.function(result[[1]]))
})

test_that("fct2list returns correct parameter names for a 4-param function", {
  result <- drc:::fct2list(ll4_raw, 2)

  expect_equal(result[[3]], c("a", "b", "c", "d"))
  expect_null(result[[2]])
  expect_true(is.function(result[[1]]))
})

test_that("fct2list works with single-line function", {
  result <- drc:::fct2list(single_line_fct, 2)

  expect_type(result, "list")
  expect_length(result, 3)
  expect_true(is.function(result[[1]]))
  expect_null(result[[2]])
  expect_equal(result[[3]], c("a", "b"))
})

test_that("fct2list errors when argument number is invalid", {
  expect_error(
    drc:::fct2list(simple_2p, 10),
    "Argument number does not exist"
  )
})

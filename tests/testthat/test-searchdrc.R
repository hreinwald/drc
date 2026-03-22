# Tests for searchdrc function

# ---------- Helper: create a valid drc model for reuse ----------
local_model <- function() {
  drm(rootl ~ conc, data = ryegrass, fct = LL.4())
}

# ================================================================
# 1. Happy-path / Correctness tests
# ================================================================

test_that("searchdrc returns a drc object when convergence is achieved", {
  m1 <- local_model()
  result <- searchdrc(m1, which = "b", range = c(0.1, 10))
  expect_s3_class(result, "drc")
})

test_that("searchdrc respects len parameter", {
  m1 <- local_model()
  result <- searchdrc(m1, which = "e", range = c(1, 10), len = 5)
  expect_s3_class(result, "drc")
})

test_that("searchdrc works with different model specifications", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
  result <- searchdrc(m1, which = "e", range = c(1, 10))
  expect_s3_class(result, "drc")
})

# ================================================================
# 2. Verbose output tests
# ================================================================

test_that("searchdrc prints attempt messages when verbose = TRUE", {
  m1 <- local_model()
  expect_message(
    searchdrc(m1, which = "b", range = c(0.1, 10), len = 3, verbose = TRUE),
    "\\[searchdrc\\] Attempt"
  )
})

test_that("searchdrc prints convergence message when verbose = TRUE", {
  m1 <- local_model()
  expect_message(
    searchdrc(m1, which = "b", range = c(0.1, 10), len = 3, verbose = TRUE),
    "Convergence achieved"
  )
})

# ================================================================
# 3. Input validation error tests
# ================================================================

test_that("searchdrc errors when object is not class 'drc'", {
  expect_error(
    searchdrc(object = list(a = 1), which = "b", range = c(1, 10)),
    "'object' must be of class 'drc'"
  )
  expect_error(
    searchdrc(object = "not_a_model", which = "b", range = c(1, 10)),
    "'object' must be of class 'drc'"
  )
})

test_that("searchdrc errors when object has no $start", {
  m1 <- local_model()
  m1$start <- NULL
  expect_error(
    searchdrc(m1, which = "b", range = c(1, 10)),
    "\\$start.*\\$parNames"
  )
})

test_that("searchdrc errors when object has no $parNames", {
  m1 <- local_model()
  m1$parNames <- NULL
  expect_error(
    searchdrc(m1, which = "b", range = c(1, 10)),
    "\\$start.*\\$parNames"
  )
})

test_that("searchdrc errors when 'which' is not a single non-empty string", {
  m1 <- local_model()
  # Not a character
  expect_error(
    searchdrc(m1, which = 42, range = c(1, 10)),
    "'which' must be a single non-empty character string"
  )
  # Length > 1
  expect_error(
    searchdrc(m1, which = c("b", "c"), range = c(1, 10)),
    "'which' must be a single non-empty character string"
  )
  # Empty string
  expect_error(
    searchdrc(m1, which = "", range = c(1, 10)),
    "'which' must be a single non-empty character string"
  )
  # Whitespace-only
  expect_error(
    searchdrc(m1, which = "  ", range = c(1, 10)),
    "'which' must be a single non-empty character string"
  )
})

test_that("searchdrc errors when 'range' is not numeric or not length 2", {
  m1 <- local_model()
  # Not numeric
  expect_error(
    searchdrc(m1, which = "b", range = c("a", "b")),
    "'range' must be a numeric vector of exactly length 2"
  )
  # Length 1
  expect_error(
    searchdrc(m1, which = "b", range = 5),
    "'range' must be a numeric vector of exactly length 2"
  )
  # Length 3
  expect_error(
    searchdrc(m1, which = "b", range = c(1, 2, 3)),
    "'range' must be a numeric vector of exactly length 2"
  )
})

test_that("searchdrc errors when range endpoints are equal", {
  m1 <- local_model()
  expect_error(
    searchdrc(m1, which = "b", range = c(5, 5)),
    "two endpoints of 'range' must be different"
  )
})

test_that("searchdrc errors when 'len' is invalid", {
  m1 <- local_model()
  # Not numeric
  expect_error(
    searchdrc(m1, which = "b", range = c(1, 10), len = "ten"),
    "'len' must be a single numeric value of at least 2"
  )
  # Less than 2
  expect_error(
    searchdrc(m1, which = "b", range = c(1, 10), len = 1),
    "'len' must be a single numeric value of at least 2"
  )
  # Length > 1
  expect_error(
    searchdrc(m1, which = "b", range = c(1, 10), len = c(5, 10)),
    "'len' must be a single numeric value of at least 2"
  )
})

test_that("searchdrc errors when 'verbose' is invalid", {
  m1 <- local_model()
  # Not logical
  expect_error(
    searchdrc(m1, which = "b", range = c(1, 10), verbose = "yes"),
    "'verbose' must be a single logical value"
  )
  # Length > 1
  expect_error(
    searchdrc(m1, which = "b", range = c(1, 10), verbose = c(TRUE, FALSE)),
    "'verbose' must be a single logical value"
  )
})

# ================================================================
# 4. Parameter matching tests
# ================================================================

test_that("searchdrc errors when parameter name is not found", {
  m1 <- local_model()
  expect_error(
    searchdrc(m1, which = "nonexistent", range = c(0.1, 10)),
    "No such parameter"
  )
})

test_that("searchdrc warns when multiple parameters match", {
  m1 <- local_model()
  # Duplicate the "b" parameter entry so two indices match "^b(:|$)"
  m1$parNames[[2]] <- c("b", "b", "c", "d")
  m1$start <- rep(m1$start[1], 4)  # match length
  expect_warning(
    searchdrc(m1, which = "b", range = c(0.1, 10), len = 3),
    "Multiple parameters matched"
  )
})

# ================================================================
# 5. Convergence failure test
# ================================================================

test_that("searchdrc warns when convergence fails across all attempts", {
  m1 <- local_model()
  # Force convergence failure by corrupting the model call
  m1$call$fct <- quote(LL.2())

  expect_warning(
    result <- searchdrc(m1, which = "b", range = c(1e10, 1e11), len = 2),
    "Convergence failed"
  )
  expect_null(result)
})

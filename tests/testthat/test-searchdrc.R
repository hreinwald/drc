# Tests for searchdrc function

test_that("searchdrc returns a drc object when convergence is achieved", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- searchdrc(m1, which = "b", range = c(0.1, 10))

  expect_s3_class(result, "drc")
})

test_that("searchdrc errors on invalid parameter name", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_error(searchdrc(m1, which = "nonexistent", range = c(0.1, 10)),
               "No parameter matching")
})

test_that("searchdrc warns when convergence fails", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_error(
    searchdrc(m1, which = "b", range = c(1e10, 1e11), len = 2),
    "Convergence failed"
  )
})

test_that("searchdrc respects len parameter", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- searchdrc(m1, which = "e", range = c(1, 10), len = 5)

  expect_s3_class(result, "drc")
})

test_that("searchdrc works with different model specifications", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
  result <- searchdrc(m1, which = "e", range = c(1, 10))

  expect_s3_class(result, "drc")
})

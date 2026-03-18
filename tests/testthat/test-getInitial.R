# Tests for getInitial function

test_that("getInitial returns named vector of starting values for LL.4 model", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- getInitial(m1)

  expect_true(is.numeric(result))
  expect_true(!is.null(names(result)))
  expect_length(result, 4)
  expect_equal(names(result), c("b", "c", "d", "e"))
})

test_that("getInitial returns named vector for LL.3 model", {
  m2 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  result <- getInitial(m2)

  expect_true(is.numeric(result))
  expect_length(result, 3)
  expect_equal(names(result), c("b", "d", "e"))
})

test_that("getInitial returns named vector for LL.2 model with binomial type", {
  m3 <- drm(r/n ~ dose, weights = n, data = deguelin, fct = LL.2(), type = "binomial")
  result <- getInitial(m3)

  expect_true(is.numeric(result))
  expect_length(result, 2)
  expect_equal(names(result), c("b", "e"))
})

test_that("getInitial values match object$start", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- getInitial(m1)

  expect_equal(as.numeric(result), m1$start)
})

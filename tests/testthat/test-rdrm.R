# tests/testthat/test-rdrm.R
# Comprehensive tests for the rdrm() function

# ---- Happy Path: Continuous (rnorm) with numeric xerror, onlyY = FALSE ----
test_that("rdrm returns x and y with numeric xerror and rnorm (default)", {
  set.seed(42)
  doses <- c(0, 1, 2, 5, 10)
  mpar <- c(1, 1, 0, 1)  # LL.4 parameters: b, c, d, e
  result <- rdrm(3, LL.4(), mpar, xerror = doses)


  expect_type(result, "list")
  expect_named(result, c("x", "y"))
  expect_true(is.matrix(result$x))
  expect_true(is.matrix(result$y))
  expect_equal(nrow(result$x), 3)
  expect_equal(ncol(result$x), length(doses))
  expect_equal(nrow(result$y), 3)
  expect_equal(ncol(result$y), length(doses))
  # x values should be sorted doses repeated

  expect_equal(as.numeric(result$x[1, ]), sort(doses))
  expect_equal(as.numeric(result$x[2, ]), sort(doses))
})

# ---- Happy Path: Continuous (rnorm) with numeric xerror, onlyY = TRUE ----
test_that("rdrm returns only y when onlyY = TRUE with rnorm", {
  set.seed(42)
  doses <- c(0, 1, 5, 10)
  mpar <- c(1, 1, 0, 1)
  result <- rdrm(2, LL.4(), mpar, xerror = doses, onlyY = TRUE)

  expect_type(result, "list")
  expect_named(result, "y")
  expect_true(is.matrix(result$y))
  expect_equal(nrow(result$y), 2)
  expect_equal(ncol(result$y), length(doses))
})

# ---- Happy Path: xerror as character string (function name) ----
test_that("rdrm works with xerror as character string (runif)", {
  set.seed(42)
  mpar <- c(1, 1, 0, 1)
  result <- rdrm(2, LL.4(), mpar, xerror = "runif", xpar = c(5, 0, 10))

  expect_type(result, "list")
  expect_named(result, c("x", "y"))
  expect_true(is.matrix(result$x))
  expect_true(is.matrix(result$y))
  expect_equal(nrow(result$x), 2)
  expect_equal(ncol(result$x), 5)  # 5 dose values from runif(5, 0, 10)
})

# ---- Happy Path: Binomial (rbinom) with ypar length > 1, onlyY = FALSE ----
test_that("rdrm works with rbinom and ypar length > 1", {
  set.seed(42)
  doses <- c(0, 1, 5, 10)
  mpar <- c(1, 0.5)  # LL.2 parameters: b, e
  weights <- c(20, 20, 20, 20)
  result <- rdrm(2, LL.2(), mpar, xerror = doses,
                 yerror = "rbinom", ypar = weights)

  expect_type(result, "list")
  expect_named(result, c("x", "w", "y"))
  expect_true(is.matrix(result$x))
  expect_true(is.matrix(result$w))
  expect_true(is.matrix(result$y))
  expect_equal(nrow(result$x), 2)
  expect_equal(ncol(result$x), length(doses))
  expect_equal(nrow(result$w), 2)
  expect_equal(ncol(result$w), length(doses))
  expect_equal(nrow(result$y), 2)
  expect_equal(ncol(result$y), length(doses))
})

# ---- Happy Path: Binomial (rbinom) with ypar length == 1 ----
test_that("rdrm works with rbinom and ypar length == 1", {
  set.seed(42)
  doses <- c(0, 1, 5, 10)
  mpar <- c(1, 0.5)  # LL.2 parameters
  result <- rdrm(2, LL.2(), mpar, xerror = doses,
                 yerror = "rbinom", ypar = 20)

  expect_type(result, "list")
  expect_named(result, c("x", "w", "y"))
  expect_true(is.matrix(result$x))
  expect_true(is.matrix(result$w))
  expect_true(is.matrix(result$y))
  # All weights should be 20
  expect_true(all(result$w == 20))
})

# ---- Binomial with onlyY = TRUE ----
test_that("rdrm with rbinom and onlyY = TRUE returns only y", {
  set.seed(42)
  doses <- c(0, 1, 5, 10)
  mpar <- c(1, 0.5)
  result <- rdrm(2, LL.2(), mpar, xerror = doses,
                 yerror = "rbinom", ypar = c(20, 20, 20, 20), onlyY = TRUE)

  expect_type(result, "list")
  expect_named(result, "y")
  expect_true(is.matrix(result$y))
  expect_equal(nrow(result$y), 2)
  expect_equal(ncol(result$y), length(doses))
})

# ---- Binomial with onlyY = TRUE and ypar length == 1 ----
test_that("rdrm with rbinom, ypar length 1, and onlyY = TRUE", {
  set.seed(42)
  doses <- c(0, 1, 5)
  mpar <- c(1, 0.5)
  result <- rdrm(2, LL.2(), mpar, xerror = doses,
                 yerror = "rbinom", ypar = 10, onlyY = TRUE)

  expect_type(result, "list")
  expect_named(result, "y")
  expect_true(is.matrix(result$y))
})

# ---- Edge case: Single dose value ----
test_that("rdrm works with a single dose value", {
  set.seed(42)
  result <- rdrm(2, LL.4(), c(1, 0, 1, 5), xerror = 5)

  expect_type(result, "list")
  expect_equal(ncol(result$x), 1)
  expect_equal(ncol(result$y), 1)
})

# ---- Edge case: nosim = 1 ----
test_that("rdrm works with nosim = 1", {
  set.seed(42)
  doses <- c(0, 1, 5, 10)
  result <- rdrm(1, LL.4(), c(1, 0, 1, 5), xerror = doses)

  expect_type(result, "list")
  expect_equal(nrow(result$x), 1)
  expect_equal(nrow(result$y), 1)
})

# ---- Reproducibility with set.seed ----
test_that("rdrm produces reproducible results with set.seed", {
  doses <- c(0, 1, 5, 10)
  mpar <- c(1, 0, 1, 5)

  set.seed(123)
  r1 <- rdrm(3, LL.4(), mpar, xerror = doses)

  set.seed(123)
  r2 <- rdrm(3, LL.4(), mpar, xerror = doses)

  expect_equal(r1$x, r2$x)
  expect_equal(r1$y, r2$y)
})

# ---- xerror as character with xpar default ----
test_that("rdrm works with xerror as character and xpar = 1 (default)", {
  set.seed(42)
  mpar <- c(1, 0, 1, 5)
  result <- rdrm(2, LL.4(), mpar, xerror = "rnorm", xpar = c(5, 0, 1))

  expect_type(result, "list")
  expect_named(result, c("x", "y"))
  expect_equal(nrow(result$x), 2)
  expect_equal(ncol(result$x), 5)  # 5 values from rnorm(5, 0, 1)
})

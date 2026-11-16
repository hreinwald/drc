# Tests for the robust distance functions in R/drmRobust.R

test_that("drm() rejects unknown 'robust' values with a clear error", {
  expect_error(drm(rootl ~ conc, data = ryegrass, fct = LL.4(), robust = "huber"),
               "Invalid 'robust' argument")
})

test_that("all documented robust options are accepted", {
  for (r in c("mean", "median", "trimmed", "tukey", "winsor", "lms", "lts")) {
    f <- suppressWarnings(try(drm(rootl ~ conc, data = ryegrass, fct = LL.4(), robust = r), silent = TRUE))
    expect_false(inherits(f, "try-error") && grepl("Invalid 'robust'", f), info = r)
  }
})

test_that("winsorised metric handles missing values without erroring", {
  call <- quote(drm(rootl ~ conc, data = ryegrass, fct = LL.4(), robust = "winsor"))
  fct <- drmRobust("winsor", call, lenData = nrow(ryegrass), lenPar = 4)
  x <- c(0.1, NA, 5, -20, 0.3)
  res <- fct(x)
  expect_length(res, length(x))
  expect_true(is.na(res[2]))
  expect_true(all(is.finite(res[-2])))
})

test_that("trimmed and Tukey metrics propagate missing values", {
  call <- quote(drm(rootl ~ conc, data = ryegrass, fct = LL.4(), robust = "trimmed"))
  expect_true(is.na(drmRobust("trimmed", call, nrow(ryegrass), 4)(c(1, NA))[2]))
  call$robust <- "tukey"
  expect_true(is.na(drmRobust("tukey", call, nrow(ryegrass), 4)(c(1, NA))[2]))
})

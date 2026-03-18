# tests/testthat/test-simDR.R
# Comprehensive tests for simDR() in R/simDR.R

# ==============================================================================
# Helper: fit a model and extract parameters for tests
# ==============================================================================
get_test_params <- function() {
  data(ryegrass, package = "drc")
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  list(
    mpar  = coef(m1),
    sigma = sqrt(summary(m1)$resVar),
    conc  = c(1.88, 3.75, 7.50, 0.94, 15, 0.47, 30, 0.23, 60)
  )
}

# ==============================================================================
# Happy path tests
# ==============================================================================

test_that("simDR returns correct structure with default edVec", {
  p <- get_test_params()

  result <- expect_output(
    simDR(p$mpar, p$sigma, LL.4(), noSim = 2, conc = p$conc, seedVal = 12345),
    "Concentrations used"
  )

  # Return value is a list with element "se"
  expect_type(result, "list")
  expect_named(result, "se")

  # se is an array: (length(conc)-4) x 6 x length(edVec)
  # Default edVec = c(10, 50) => 2 ED values
  expect_equal(dim(result$se), c(5, 6, 2))
  expect_true(is.numeric(result$se))
})

test_that("simDR handles single ED value in edVec", {
  p <- get_test_params()

  result <- expect_output(
    simDR(p$mpar, p$sigma, LL.4(), noSim = 2, conc = p$conc,
          edVec = c(50), seedVal = 12345),
    "ED value considered: 50"
  )

  # Single ED value => third dimension is 1
  expect_equal(dim(result$se), c(5, 6, 1))
})

test_that("simDR handles three ED values in edVec", {
  p <- get_test_params()

  result <- expect_output(
    simDR(p$mpar, p$sigma, LL.4(), noSim = 2, conc = p$conc,
          edVec = c(10, 50, 90), seedVal = 12345),
    "ED value considered: 90"
  )

  expect_equal(dim(result$se), c(5, 6, 3))
})

# ==============================================================================
# Output tests
# ==============================================================================

test_that("simDR output includes concentrations and ED information", {
  p <- get_test_params()

  output <- capture.output(
    result <- simDR(p$mpar, p$sigma, LL.4(), noSim = 2, conc = p$conc,
                    edVec = c(10, 50), seedVal = 12345)
  )

  expect_true(any(grepl("Concentrations used:", output)))
  expect_true(any(grepl("ED value considered: 10", output)))
  expect_true(any(grepl("ED value considered: 50", output)))
  expect_true(any(grepl("Conc. no.\\Replicates:", output, fixed = TRUE)))
})

test_that("simDR returns result invisibly", {
  p <- get_test_params()

  expect_invisible(
    simDR(p$mpar, p$sigma, LL.4(), noSim = 2, conc = p$conc, seedVal = 12345)
  )
})

# ==============================================================================
# Row/column naming tests (bug fix: dynamic rownames)
# ==============================================================================

test_that("simDR correctly labels rows for different conc lengths", {
  p <- get_test_params()

  # 9 concentrations: rows should be named 5:9
  output9 <- capture.output(
    res9 <- simDR(p$mpar, p$sigma, LL.4(), noSim = 2, conc = p$conc, seedVal = 12345)
  )
  mat9 <- res9$se[, , 1]
  expect_equal(nrow(mat9), 5)

  # Verify row names are correctly set in the output by checking for "5" and "9"
  expect_true(any(grepl("^5 ", output9)))
  expect_true(any(grepl("^9 ", output9)))
})

# ==============================================================================
# Seed reproducibility test
# ==============================================================================

test_that("simDR is reproducible with same seed", {
  p <- get_test_params()

  output1 <- capture.output(
    res1 <- simDR(p$mpar, p$sigma, LL.4(), noSim = 2, conc = p$conc, seedVal = 999)
  )

  output2 <- capture.output(
    res2 <- simDR(p$mpar, p$sigma, LL.4(), noSim = 2, conc = p$conc, seedVal = 999)
  )

  expect_equal(res1$se, res2$se)
})

# ==============================================================================
# try-error branch test
# ==============================================================================

test_that("simDR handles drm fitting failures gracefully", {
  p <- get_test_params()

  # Use extremely large sigma to increase chance of fitting failures.
  # Even if all fits succeed, the function should still run without error.
  # Use very small noSim to keep test fast.
  output <- capture.output(
    result <- simDR(p$mpar, sigma = 1e6, LL.4(), noSim = 2,
                    conc = p$conc, seedVal = 42)
  )

  expect_type(result, "list")
  expect_named(result, "se")
  expect_equal(dim(result$se), c(5, 6, 2))
})

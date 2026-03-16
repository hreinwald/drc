# Test file for findbe.R functions (findbe1, findbe2, findbe3)
# These are internal helper functions for finding initial b and e parameter estimates

# ==============================================================================
# Common helper functions (mirroring those used in self-starter functions)
# ==============================================================================

ytrans <- function(y, cVal, dVal) { log((dVal - y) / (y - cVal)) }
bfct_helper <- function(x, y, cVal, dVal, eVal) { ytrans(y, cVal, dVal) / log(x / eVal) }
efct_helper <- function(x, y, bVal, cVal, dVal) { x * exp(-ytrans(y, cVal, dVal) / bVal) }
doseTr_log <- function(x) { rVec <- log(x); rVec[!x > 0] <- NA; rVec }

# ==============================================================================
# Tests for findbe1
# ==============================================================================

test_that("findbe1 returns a closure", {
  fn <- drc:::findbe1(doseTr_log, ytrans)
  expect_type(fn, "closure")
  expect_true(is.function(fn))
})

test_that("findbe1 returns correct estimates with standard decreasing data", {
  fn <- drc:::findbe1(doseTr_log, ytrans)

  x <- c(0.1, 0.5, 1, 2, 5, 10, 50, 100)
  y <- c(0.95, 0.85, 0.75, 0.55, 0.25, 0.1, 0.03, 0.01)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(all(is.finite(result)))
  # b should be positive for decreasing response
  expect_true(result[1] > 0)
  # e should be positive
  expect_true(result[2] > 0)
})

test_that("findbe1 works with sgnb = -1", {
  fn <- drc:::findbe1(doseTr_log, ytrans, sgnb = -1)

  x <- c(0.1, 0.5, 1, 2, 5, 10, 50, 100)
  y <- c(0.95, 0.85, 0.75, 0.55, 0.25, 0.1, 0.03, 0.01)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
})

test_that("findbe1 works with custom back function", {
  fn <- drc:::findbe1(doseTr_log, ytrans, back = function(x) 10^x)

  x <- c(0.1, 0.5, 1, 2, 5, 10, 50, 100)
  y <- c(0.95, 0.85, 0.75, 0.55, 0.25, 0.1, 0.03, 0.01)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
})

# ==============================================================================
# Tests for findbe2 - Anke method
# ==============================================================================

test_that("findbe2 Anke returns a closure", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Anke")
  expect_type(fn, "closure")
  expect_true(is.function(fn))
})

test_that("findbe2 Anke returns valid estimates with standard data", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Anke")

  x <- c(0.1, 0.5, 1, 2, 5, 10, 50, 100)
  y <- c(0.95, 0.85, 0.7, 0.55, 0.25, 0.1, 0.03, 0.01)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
})

test_that("findbe2 Anke: mixed responses at dose triggers aboveVec trim (line 51)", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Anke")

  # At dose=2, one response above midResp (0.6) and one below (0.3)
  # This triggers: length(aboveVec) < sum(x %in% uniAbove) on line 49
  x <- c(0.1, 0.5, 1, 2, 2, 5, 10, 50)
  y <- c(0.95, 0.85, 0.7, 0.6, 0.3, 0.2, 0.05, 0.01)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
})

test_that("findbe2 Anke: mixed responses triggers belowVec trim (line 62)", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Anke")

  # At dose=5, one response above midResp (0.6) and one below (0.3)
  # This triggers: length(belowVec) < sum(x %in% uniBelow) on line 59
  x <- c(0.1, 0.5, 1, 5, 5, 10, 50, 100)
  y <- c(0.95, 0.85, 0.7, 0.6, 0.3, 0.1, 0.03, 0.01)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
})

test_that("findbe2 Anke: NaN eVal when no doses between max and min (line 69-71)", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Anke")

  # Adjacent dose levels with all responses above/below midResp
  # subsetInd = (x > maxDose) & (x < minDose) is all FALSE
  x <- c(1, 2, 3, 4)
  y <- c(0.9, 0.7, 0.3, 0.1)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(all(is.finite(result)))
})

test_that("findbe2 Anke: eVal < sort1 triggers correction (line 74-76)", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Anke")

  # maxDose=1, minDose=3 → eVal=(3+1)/2=2 but sort1=3
  # Actually need (NaN fallback) eVal < sort1
  x <- c(1, 3, 5, 7)
  y <- c(0.6, 0.3, 0.1, 0.01)
  cVal <- 0
  dVal <- 1
  # midResp = 0.5, aboveVec = c(1), maxDose = 1
  # belowVec = c(3, 5, 7), minDose = 3
  # subsetInd: no x between 1 and 3 → NaN
  # eVal = (3+1)/2 = 2, sort1 = 3 → eVal < sort1 → eVal = 3

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
})

test_that("findbe2 Anke: eVal > sort2 triggers correction (line 79-81)", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Anke")

  # maxDose=5, minDose=7 → eVal=(7+5)/2=6 but sort2=5
  x <- c(1, 3, 5, 7)
  y <- c(0.95, 0.7, 0.6, 0.1)
  cVal <- 0
  dVal <- 1
  # midResp = 0.5, aboveVec = c(1,3,5), maxDose = 5
  # belowVec = c(7), minDose = 7
  # subsetInd: no x between 5 and 7 → NaN
  # eVal = (7+5)/2 = 6, sort2 = 5 → eVal > sort2 → eVal = 5

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
})

test_that("findbe2 Anke: sign correction triggers with sgnb=-1 (line 92-98)", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Anke", sgnb = -1)

  # Standard decreasing data: regSlope < 0, bVal > 0 initially
  # sgnb*regSlope/bVal = (-1)*(neg)/(pos) = pos > 0 → triggers correction
  x <- c(0.1, 0.5, 1, 2, 5, 10, 50, 100)
  y <- c(0.95, 0.85, 0.7, 0.55, 0.25, 0.1, 0.03, 0.01)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
})

test_that("findbe2 Anke: NA bVal triggers fallback (line 99-102)", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Anke")

  # y values outside [cVal, dVal] make all ytrans = NaN
  # so bFct returns NaN (treated as NA)
  x <- c(1, 2, 3, 4, 5)
  y <- c(1.5, 1.3, 1.1, -0.1, -0.3)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
  # bVal should be from the fallback: sgnb * (-regSlope)
  expect_true(is.finite(result[1]))
})

# ==============================================================================
# Tests for findbe2 - Normolle method
# ==============================================================================

test_that("findbe2 Normolle returns a closure", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Normolle")
  expect_type(fn, "closure")
  expect_true(is.function(fn))
})

test_that("findbe2 Normolle returns valid estimates with standard data", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Normolle")

  x <- c(0.1, 0.5, 1, 2, 5, 10, 50, 100)
  y <- c(0.95, 0.85, 0.7, 0.55, 0.25, 0.1, 0.03, 0.01)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
})

test_that("findbe2 Normolle works with sgnb = -1", {
  fn <- drc:::findbe2(bfct_helper, efct_helper, "Normolle", sgnb = -1)

  x <- c(0.1, 0.5, 1, 2, 5, 10, 50, 100)
  y <- c(0.95, 0.85, 0.7, 0.55, 0.25, 0.1, 0.03, 0.01)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
})

# ==============================================================================
# Tests for findbe3
# ==============================================================================

test_that("findbe3 returns a closure", {
  fn <- drc:::findbe3()
  expect_type(fn, "closure")
  expect_true(is.function(fn))
})

test_that("findbe3 returns valid estimates with decreasing response (crit2 path)", {
  fn <- drc:::findbe3()

  # Crossing at i=4: uniy[4]=0.3 < 0.5, uniy[3]=0.6 > 0.5 → crit2 TRUE
  x <- c(1, 2, 3, 4, 5)
  y <- c(0.95, 0.8, 0.6, 0.3, 0.1)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(all(is.finite(result)))
  # For decreasing: bVal should be positive with sgnb=1
  expect_equal(result[1], 1)
  # eVal should be (unix[4]+unix[3])/2 = 3.5
  expect_equal(result[2], 3.5)
})

test_that("findbe3 returns valid estimates with increasing response (crit1 path)", {
  fn <- drc:::findbe3()

  # Crossing at i=3: uniy[3]=0.6 > 0.5, uniy[2]=0.3 < 0.5 → crit1 TRUE
  x <- c(1, 2, 3, 4, 5)
  y <- c(0.1, 0.3, 0.6, 0.8, 0.95)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(all(is.finite(result)))
  # For increasing: bVal should be negative with sgnb=1 (sign of -(uniy[j]-uniy[j-1]))
  expect_equal(result[1], -1)
  # eVal should be (unix[3]+unix[2])/2 = 2.5
  expect_equal(result[2], 2.5)
})

test_that("findbe3 works with sgnb = -1", {
  fn <- drc:::findbe3(sgnb = -1)

  x <- c(1, 2, 3, 4, 5)
  y <- c(0.95, 0.8, 0.6, 0.3, 0.1)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
  # With sgnb=-1, bVal sign flips
  expect_equal(result[1], -1)
})

test_that("findbe3 works with replicated doses", {
  fn <- drc:::findbe3()

  x <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
  y <- c(0.95, 0.9, 0.85, 0.75, 0.6, 0.55, 0.3, 0.25, 0.1, 0.05)
  cVal <- 0
  dVal <- 1

  result <- fn(x, y, cVal, dVal)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(all(is.finite(result)))
})

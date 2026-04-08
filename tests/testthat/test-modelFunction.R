# Tests for modelFunction (R/modelFunction.R)
# An internal function that creates a model evaluation closure (multCurves)

# --- Helper setup ---

# Simple parm2mat: converts a parameter vector into a matrix
# For n observations and p parameters, returns n x p matrix
make_parm2mat <- function(n, p) {
  function(parm) {
    matrix(parm, nrow = n, ncol = p, byrow = TRUE)
  }
}

# Simple drcFct: linear dose-response f(dose, parm) = parm[,1] * dose + parm[,2]
simple_drcFct <- function(dose, parm) {
  parm[, 1] * dose + parm[, 2]
}

# ============================================================
# Test Block 1: Basic functionality with cm = NULL, retFct = NULL, pshifts = NULL
# ============================================================
test_that("modelFunction returns a function when cm=NULL, retFct=NULL, pshifts=NULL", {
  n <- 5
  dose <- c(0, 1, 2, 3, 4)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)

  result <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = NULL, assayNoOld = rep(1, n), upperPos = 2,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = NULL
  )

  expect_true(is.function(result))
})

test_that("modelFunction (cm=NULL) evaluates correctly with simple linear function", {
  n <- 5
  dose <- c(0, 1, 2, 3, 4)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = NULL, assayNoOld = rep(1, n), upperPos = 2,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = NULL
  )

  # parm = c(2, 3) -> each row is [2, 3]
  # f(dose, parm) = 2*dose + 3
  out <- multCurves(dose, c(2, 3))
  expected <- 2 * dose + 3
  expect_equal(out, expected)
})

# ============================================================
# Test Block 2: retFct is not NULL
# ============================================================
test_that("modelFunction uses retFct to replace drcFct when retFct is not NULL", {
  n <- 3
  dose <- c(1, 2, 3)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)

  # retFct returns a function that scales dose response by doseScaling * respScaling
  my_retFct <- function(doseScaling, respScaling) {
    scale <- doseScaling * respScaling
    function(dose, parm) {
      scale * (parm[, 1] * dose + parm[, 2])
    }
  }

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = NULL, assayNoOld = rep(1, n), upperPos = 2,
    retFct = my_retFct, doseScaling = 2, respScaling = 3,
    isFinite = isFinite, pshifts = NULL
  )

  out <- multCurves(dose, c(1, 0))
  # scale = 2*3 = 6; f = 6 * (1*dose + 0) = 6*dose
  expect_equal(out, 6 * dose)
})

# ============================================================
# Test Block 3: pshifts is not NULL (cm = NULL path)
# ============================================================
test_that("modelFunction applies pshifts when dimensions match (cm=NULL path)", {
  n <- 3
  dose <- c(1, 2, 3)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)

  # pshifts adds [1, 1] to each row of parmVal
  pshifts <- matrix(1, nrow = n, ncol = 2)

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = NULL, assayNoOld = rep(1, n), upperPos = 2,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = pshifts
  )

  # parm = c(2, 3) -> parmVal = [[2,3],[2,3],[2,3]]
  # pshifts = [[1,1],[1,1],[1,1]]
  # parmVal + pshifts = [[3,4],[3,4],[3,4]]
  # f(dose, parm) = 3*dose + 4
  out <- multCurves(dose, c(2, 3))
  expected <- 3 * dose + 4
  expect_equal(out, expected)
})

test_that("modelFunction does NOT apply pshifts when dimensions do not match (cm=NULL)", {
  n <- 3
  dose <- c(1, 2, 3)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)

  # Wrong dimensions: 4 x 2 instead of 3 x 2
  pshifts <- matrix(1, nrow = 4, ncol = 2)

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = NULL, assayNoOld = rep(1, n), upperPos = 2,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = pshifts
  )

  # pshifts not applied because dim mismatch
  out <- multCurves(dose, c(2, 3))
  expected <- 2 * dose + 3
  expect_equal(out, expected)
})

# ============================================================
# Test Block 4: cm is NOT NULL (control measurement path)
# ============================================================
test_that("modelFunction with cm != NULL separates control and non-control observations", {
  n <- 4
  dose <- c(0, 1, 2, 3)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)
  assayNoOld <- c(1, 1, 2, 2)
  cm <- 1  # assay 1 is control
  upperPos <- 2  # second column is the upper asymptote

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = cm, assayNoOld = assayNoOld, upperPos = upperPos,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = NULL
  )

  parm <- c(2, 5)
  out <- multCurves(dose, parm)

  # iv = isFinite & (assayNoOld == cm) -> TRUE, TRUE, FALSE, FALSE
  # niv = FALSE, FALSE, TRUE, TRUE
  # parmVal (each row [2, 5]):
  # For iv (obs 1,2): fctEval = parmVal[, upperPos] = 5
  # For niv (obs 3,4): fctEval = simple_drcFct(dose[niv], parmVal[niv,]) = 2*dose + 5
  expected <- c(5, 5, 2 * 2 + 5, 2 * 3 + 5)
  expect_equal(out, expected)
})

# ============================================================
# Test Block 5: cm != NULL with pshifts
# ============================================================
test_that("modelFunction with cm != NULL applies pshifts when dimensions match", {
  n <- 4
  dose <- c(0, 1, 2, 3)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)
  assayNoOld <- c(1, 1, 2, 2)
  cm <- 1
  upperPos <- 2

  pshifts <- matrix(c(0.5, 1), nrow = n, ncol = 2, byrow = TRUE)

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = cm, assayNoOld = assayNoOld, upperPos = upperPos,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = pshifts
  )

  parm <- c(2, 5)
  out <- multCurves(dose, parm)

  # parmVal (each row [2,5]) + pshifts (each row [0.5,1]) = each row [2.5, 6]
  # iv (obs 1,2): fctEval = parmVal[, 2] = 6
  # niv (obs 3,4): fctEval = 2.5 * dose + 6
  expected <- c(6, 6, 2.5 * 2 + 6, 2.5 * 3 + 6)
  expect_equal(out, expected)
})

test_that("modelFunction with cm != NULL does NOT apply pshifts when dims don't match", {
  n <- 4
  dose <- c(0, 1, 2, 3)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)
  assayNoOld <- c(1, 1, 2, 2)
  cm <- 1
  upperPos <- 2

  # Wrong dims
  pshifts <- matrix(1, nrow = 5, ncol = 2)

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = cm, assayNoOld = assayNoOld, upperPos = upperPos,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = pshifts
  )

  parm <- c(2, 5)
  out <- multCurves(dose, parm)

  # pshifts NOT applied
  expected <- c(5, 5, 2 * 2 + 5, 2 * 3 + 5)
  expect_equal(out, expected)
})

# ============================================================
# Test Block 6: retFct with cm != NULL
# ============================================================
test_that("modelFunction with retFct and cm != NULL uses retFct for drcFct", {
  n <- 4
  dose <- c(0, 1, 2, 3)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)
  assayNoOld <- c(1, 1, 2, 2)
  cm <- 1
  upperPos <- 2

  # retFct: returns a scaled drcFct
  my_retFct <- function(doseScaling, respScaling) {
    scale <- doseScaling * respScaling
    function(dose, parm) {
      scale * (parm[, 1] * dose + parm[, 2])
    }
  }

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = cm, assayNoOld = assayNoOld, upperPos = upperPos,
    retFct = my_retFct, doseScaling = 2, respScaling = 3,
    isFinite = isFinite, pshifts = NULL
  )

  parm <- c(1, 10)
  out <- multCurves(dose, parm)

  # Note: in the cm != NULL path, drcFct is replaced by retFct result
  # BUT the cm path uses drcFct directly (not drcFct1), so retFct replacement affects both paths
  # iv: fctEval = parmVal[, upperPos] = 10 (not through drcFct)
  # niv: drcFct(dose[niv], parmVal[niv,]) = 6 * (1*dose + 10)
  expected <- c(10, 10, 6 * (1 * 2 + 10), 6 * (1 * 3 + 10))
  expect_equal(out, expected)
})

# ============================================================
# Test Block 7: isFinite with some FALSE values (cm=NULL path)
# ============================================================
test_that("modelFunction with isFinite subset (cm=NULL)", {
  n <- 4
  dose <- c(0, 1, 2, 3)
  isFinite <- rep(TRUE, n)

  # parm2mat needs to return n x p matrix
  parm2mat <- make_parm2mat(n, 2)

  # drcFct receives only the rows where isFinite is TRUE
  drcFct_sub <- function(dose, parm) {
    parm[, 1] * dose + parm[, 2]
  }

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = drcFct_sub,
    cm = NULL, assayNoOld = rep(1, n), upperPos = 2,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = NULL
  )

  out <- multCurves(dose, c(2, 3))
  # All isFinite = TRUE, so all rows selected
  # drcFct(dose, parmVal) = 2*dose + 3
  expected <- 2 * dose + 3
  expect_equal(out, expected)
  expect_equal(length(out), n)
})

# ============================================================
# Test Block 8: Edge case - single observation
# ============================================================
test_that("modelFunction works with single observation (cm=NULL)", {
  n <- 1
  dose <- 5
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- TRUE

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = NULL, assayNoOld = 1, upperPos = 2,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = NULL
  )

  out <- multCurves(dose, c(2, 3))
  expect_equal(out, 2 * 5 + 3)
})

test_that("modelFunction works with single observation (cm != NULL, control)", {
  n <- 1
  dose <- 5
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- TRUE
  assayNoOld <- 1
  cm <- 1  # this obs is control
  upperPos <- 2

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = cm, assayNoOld = assayNoOld, upperPos = upperPos,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = NULL
  )

  out <- multCurves(dose, c(2, 5))
  # iv = TRUE (obs matches cm), niv = FALSE
  # fctEval[iv] = parmVal[, upperPos] = 5
  expect_equal(out, 5)
})

# ============================================================
# Test Block 9: Both retFct and pshifts (cm=NULL)
# ============================================================
test_that("modelFunction with both retFct and pshifts (cm=NULL)", {
  n <- 3
  dose <- c(1, 2, 3)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)
  pshifts <- matrix(c(0.1, 0.2), nrow = n, ncol = 2, byrow = TRUE)

  my_retFct <- function(doseScaling, respScaling) {
    function(dose, parm) {
      parm[, 1] * dose + parm[, 2]
    }
  }

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = NULL, assayNoOld = rep(1, n), upperPos = 2,
    retFct = my_retFct, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = pshifts
  )

  out <- multCurves(dose, c(1, 2))
  # parmVal = [[1,2],[1,2],[1,2]] + [[0.1,0.2],[0.1,0.2],[0.1,0.2]] = [[1.1, 2.2],...]
  # f(dose, parm) = 1.1*dose + 2.2
  expected <- 1.1 * dose + 2.2
  expect_equal(out, expected)
})

# ============================================================
# Test Block 10: cm != NULL, all obs are non-control
# ============================================================
test_that("modelFunction cm != NULL where no obs match control", {
  n <- 3
  dose <- c(1, 2, 3)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)
  assayNoOld <- c(2, 2, 2)
  cm <- 1  # no obs have assayNoOld==1

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = cm, assayNoOld = assayNoOld, upperPos = 2,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = NULL
  )

  out <- multCurves(dose, c(2, 5))
  # all are niv -> drcFct applied to all
  expected <- 2 * dose + 5
  expect_equal(out, expected)
})

# ============================================================
# Test Block 11: cm != NULL, all obs are control
# ============================================================
test_that("modelFunction cm != NULL where ALL obs match control", {
  n <- 3
  dose <- c(1, 2, 3)
  parm2mat <- make_parm2mat(n, 2)
  isFinite <- rep(TRUE, n)
  assayNoOld <- c(1, 1, 1)
  cm <- 1  # all obs are control

  multCurves <- drc:::modelFunction(
    dose = dose, parm2mat = parm2mat, drcFct = simple_drcFct,
    cm = cm, assayNoOld = assayNoOld, upperPos = 2,
    retFct = NULL, doseScaling = 1, respScaling = 1,
    isFinite = isFinite, pshifts = NULL
  )

  out <- multCurves(dose, c(2, 5))
  # all are iv -> fctEval = parmVal[, upperPos] = 5
  expected <- c(5, 5, 5)
  expect_equal(out, expected)
})

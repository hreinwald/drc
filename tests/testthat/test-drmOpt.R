# tests/testthat/test-drmOpt.R
# Comprehensive tests for drmOpt() — internal optim() wrapper

# --- Helper: simple quadratic objective (min at x=3, y=1) ----
quad_obj <- function(par) (par[1] - 3)^2 + (par[2] - 1)^2
quad_grad <- function(par) c(2 * (par[1] - 3), 2 * (par[2] - 1))
quad_hess <- function(par) matrix(c(2, 0, 0, 2), 2, 2)

# A function that always errors, to make optim() fail via try()
bad_obj <- function(par) stop("forced failure")
bad_grad <- function(par) stop("forced gradient failure")

# Pre-define a matchCall value to avoid match.call() errors inside test wrappers
fake_call <- call("drmOpt")

# ============================================================
# 1. No derivatives, unconstrained, success  (happy path)
# ============================================================
test_that("no derivatives, unconstrained, success", {
  res <- drc:::drmOpt(
    opfct       = quad_obj,
    opdfct1     = NULL,
    startVec    = c(0, 0),
    optMethod   = "Nelder-Mead",
    constrained = FALSE,
    warnVal     = 0,
    upperLimits = NULL,
    lowerLimits = NULL,
    errorMessage = TRUE,
    maxIt       = 500,
    relTol      = 1e-8,
    opdfct2     = NULL,
    parmVec     = c("a", "b"),
    traceVal    = 0,
    silentVal   = TRUE,
    matchCall   = fake_call
  )
  expect_true(res$convergence)
  expect_true(!is.null(res$hessian))
  expect_equal(res$par, c(3, 1), tolerance = 1e-3)
  # ovalue is the raw optim value; value is recomputed via opfct
  expect_equal(res$value, quad_obj(res$par))
})

# ============================================================
# 2. No derivatives, constrained, success
# ============================================================
test_that("no derivatives, constrained, success", {
  res <- drc:::drmOpt(
    opfct       = quad_obj,
    opdfct1     = NULL,
    startVec    = c(0, 0),
    optMethod   = "L-BFGS-B",
    constrained = TRUE,
    warnVal     = 0,
    upperLimits = c(10, 10),
    lowerLimits = c(-10, -10),
    errorMessage = TRUE,
    maxIt       = 500,
    relTol      = 1e-8,
    opdfct2     = NULL,
    parmVec     = c("a", "b"),
    traceVal    = 0,
    silentVal   = TRUE,
    matchCall   = fake_call
  )
  expect_true(res$convergence)
  expect_equal(res$par, c(3, 1), tolerance = 1e-3)
})

# ============================================================
# 3. No derivatives, convergence failure + errorMessage=TRUE → stop()
# ============================================================
test_that("no derivatives, failure with errorMessage=TRUE throws stop()", {
  expect_error(
    drc:::drmOpt(
      opfct       = bad_obj,
      opdfct1     = NULL,
      startVec    = c(0, 0),
      optMethod   = "Nelder-Mead",
      constrained = FALSE,
      warnVal     = 0,
      upperLimits = NULL,
      lowerLimits = NULL,
      errorMessage = TRUE,
      maxIt       = 10,
      relTol      = 1e-8,
      opdfct2     = NULL,
      parmVec     = c("a", "b"),
      traceVal    = 0,
      silentVal   = TRUE,
      matchCall   = fake_call
    ),
    "Convergence failed"
  )
})

# ============================================================
# 4. No derivatives, convergence failure + errorMessage=FALSE → warning
# ============================================================
test_that("no derivatives, failure with errorMessage=FALSE gives warning", {
  expect_warning(
    {
      res <- drc:::drmOpt(
        opfct       = bad_obj,
        opdfct1     = NULL,
        startVec    = c(0, 0),
        optMethod   = "Nelder-Mead",
        constrained = FALSE,
        warnVal     = 0,
        upperLimits = NULL,
        lowerLimits = NULL,
        errorMessage = FALSE,
        maxIt       = 10,
        relTol      = 1e-8,
        opdfct2     = NULL,
        parmVec     = c("a", "b"),
        traceVal    = 0,
        silentVal   = TRUE,
        matchCall   = fake_call
      )
    },
    "Convergence failed"
  )
  expect_false(res$convergence)
  expect_equal(res$startVal, c(0, 0))
  expect_equal(res$parNames, c("a", "b"))
})

# ============================================================
# 5. No derivatives, constrained, convergence failure + errorMessage=FALSE
# ============================================================
test_that("no derivatives, constrained, failure with errorMessage=FALSE gives warning", {
  expect_warning(
    {
      res <- drc:::drmOpt(
        opfct       = bad_obj,
        opdfct1     = NULL,
        startVec    = c(0, 0),
        optMethod   = "L-BFGS-B",
        constrained = TRUE,
        warnVal     = 0,
        upperLimits = c(10, 10),
        lowerLimits = c(-10, -10),
        errorMessage = FALSE,
        maxIt       = 10,
        relTol      = 1e-8,
        opdfct2     = NULL,
        parmVec     = c("a", "b"),
        traceVal    = 0,
        silentVal   = TRUE,
        matchCall   = fake_call
      )
    },
    "Convergence failed"
  )
  expect_false(res$convergence)
})

# ============================================================
# 6. With derivatives, unconstrained, success (opdfct2=NULL → hes=TRUE)
# ============================================================
test_that("with derivatives, unconstrained, success, hessian from optim", {
  res <- drc:::drmOpt(
    opfct       = quad_obj,
    opdfct1     = quad_grad,
    startVec    = c(0, 0),
    optMethod   = "BFGS",
    constrained = FALSE,
    warnVal     = 0,
    upperLimits = NULL,
    lowerLimits = NULL,
    errorMessage = TRUE,
    maxIt       = 500,
    relTol      = 1e-8,
    opdfct2     = NULL,
    parmVec     = c("a", "b"),
    traceVal    = 0,
    silentVal   = TRUE,
    matchCall   = fake_call
  )
  expect_true(res$convergence)
  expect_equal(res$par, c(3, 1), tolerance = 1e-5)
  expect_true(!is.null(res$hessian))
  expect_equal(res$value, quad_obj(res$par))
})

# ============================================================
# 7. With derivatives, constrained, success
# ============================================================
test_that("with derivatives, constrained, success", {
  res <- drc:::drmOpt(
    opfct       = quad_obj,
    opdfct1     = quad_grad,
    startVec    = c(0, 0),
    optMethod   = "L-BFGS-B",
    constrained = TRUE,
    warnVal     = 0,
    upperLimits = c(10, 10),
    lowerLimits = c(-10, -10),
    errorMessage = TRUE,
    maxIt       = 500,
    relTol      = 1e-8,
    opdfct2     = NULL,
    parmVec     = c("a", "b"),
    traceVal    = 0,
    silentVal   = TRUE,
    matchCall   = fake_call
  )
  expect_true(res$convergence)
  expect_equal(res$par, c(3, 1), tolerance = 1e-5)
})

# ============================================================
# 8. With derivatives, unconstrained, success, opdfct2 provided (hes=FALSE)
# ============================================================
test_that("with derivatives, opdfct2 provided, hessian computed externally", {
  res <- drc:::drmOpt(
    opfct       = quad_obj,
    opdfct1     = quad_grad,
    startVec    = c(0, 0),
    optMethod   = "BFGS",
    constrained = FALSE,
    warnVal     = 0,
    upperLimits = NULL,
    lowerLimits = NULL,
    errorMessage = TRUE,
    maxIt       = 500,
    relTol      = 1e-8,
    opdfct2     = quad_hess,
    parmVec     = c("a", "b"),
    traceVal    = 0,
    silentVal   = TRUE,
    matchCall   = fake_call
  )
  expect_true(res$convergence)
  expect_equal(res$par, c(3, 1), tolerance = 1e-5)
  # hessian should be from quad_hess, i.e. diag(2,2)
  expect_equal(res$hessian, matrix(c(2, 0, 0, 2), 2, 2))
})

# ============================================================
# 9. With derivatives, convergence failure → warning + return
# ============================================================
test_that("with derivatives, convergence failure gives warning", {
  expect_warning(
    {
      res <- drc:::drmOpt(
        opfct       = bad_obj,
        opdfct1     = bad_grad,
        startVec    = c(0, 0),
        optMethod   = "BFGS",
        constrained = FALSE,
        warnVal     = 0,
        upperLimits = NULL,
        lowerLimits = NULL,
        errorMessage = TRUE,
        maxIt       = 10,
        relTol      = 1e-8,
        opdfct2     = NULL,
        parmVec     = c("a", "b"),
        traceVal    = 0,
        silentVal   = TRUE,
        matchCall   = fake_call
      )
    },
    "Convergence failed"
  )
  expect_false(res$convergence)
  expect_equal(res$startVal, c(0, 0))
  expect_equal(res$parNames, c("a", "b"))
})

# ============================================================
# 10. With derivatives, constrained, convergence failure
# ============================================================
test_that("with derivatives, constrained, convergence failure gives warning", {
  expect_warning(
    {
      res <- drc:::drmOpt(
        opfct       = bad_obj,
        opdfct1     = bad_grad,
        startVec    = c(0, 0),
        optMethod   = "L-BFGS-B",
        constrained = TRUE,
        warnVal     = 0,
        upperLimits = c(10, 10),
        lowerLimits = c(-10, -10),
        errorMessage = TRUE,
        maxIt       = 10,
        relTol      = 1e-8,
        opdfct2     = NULL,
        parmVec     = c("a", "b"),
        traceVal    = 0,
        silentVal   = TRUE,
        matchCall   = fake_call
      )
    },
    "Convergence failed"
  )
  expect_false(res$convergence)
})

# ============================================================
# 11. Edge case: startVec with very small values (psVec clamping)
# ============================================================
test_that("small startVec values are clamped in parscale", {
  # Start with very small values — triggers psVec[psVec < 1e-4] <- 1
  res <- drc:::drmOpt(
    opfct       = quad_obj,
    opdfct1     = NULL,
    startVec    = c(1e-6, 1e-6),
    optMethod   = "Nelder-Mead",
    constrained = FALSE,
    warnVal     = 0,
    upperLimits = NULL,
    lowerLimits = NULL,
    errorMessage = TRUE,
    maxIt       = 500,
    relTol      = 1e-8,
    opdfct2     = NULL,
    parmVec     = c("a", "b"),
    traceVal    = 0,
    silentVal   = TRUE,
    matchCall   = fake_call
  )
  expect_true(res$convergence)
  expect_equal(res$par, c(3, 1), tolerance = 1e-2)
})

# ============================================================
# 12. With derivatives, constrained, opdfct2 provided (hes=FALSE)
# ============================================================
test_that("with derivatives, constrained, opdfct2 provided", {
  res <- drc:::drmOpt(
    opfct       = quad_obj,
    opdfct1     = quad_grad,
    startVec    = c(0, 0),
    optMethod   = "L-BFGS-B",
    constrained = TRUE,
    warnVal     = 0,
    upperLimits = c(10, 10),
    lowerLimits = c(-10, -10),
    errorMessage = TRUE,
    maxIt       = 500,
    relTol      = 1e-8,
    opdfct2     = quad_hess,
    parmVec     = c("a", "b"),
    traceVal    = 0,
    silentVal   = TRUE,
    matchCall   = fake_call
  )
  expect_true(res$convergence)
  expect_equal(res$par, c(3, 1), tolerance = 1e-5)
  # Hessian should come from opdfct2
  expect_equal(res$hessian, matrix(c(2, 0, 0, 2), 2, 2))
})

# ============================================================
# 13. No derivatives, constrained, convergence failure + errorMessage=TRUE → stop
# ============================================================
test_that("no derivatives, constrained, failure with errorMessage=TRUE throws stop()", {
  expect_error(
    drc:::drmOpt(
      opfct       = bad_obj,
      opdfct1     = NULL,
      startVec    = c(0, 0),
      optMethod   = "L-BFGS-B",
      constrained = TRUE,
      warnVal     = 0,
      upperLimits = c(10, 10),
      lowerLimits = c(-10, -10),
      errorMessage = TRUE,
      maxIt       = 10,
      relTol      = 1e-8,
      opdfct2     = NULL,
      parmVec     = c("a", "b"),
      traceVal    = 0,
      silentVal   = TRUE,
      matchCall   = fake_call
    ),
    "Convergence failed"
  )
})

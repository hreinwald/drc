# Tests that the sandwich generics dispatch to the drc methods (estfun/bread
# must be registered as S3 methods, not just exported as plain functions).

# Helper: sandwich estimator rebuilt from its parts (Zeileis 2006)
manual_sandwich <- function(m) {
  ef <- sandwich::estfun(m)
  br <- sandwich::bread(m)
  n <- NROW(ef)
  br %*% (crossprod(ef) / n) %*% br / n
}

expect_sandwich_ok <- function(m) {
  np <- length(coef(m))
  s <- sandwich::sandwich(m)
  expect_true(is.matrix(s))
  expect_equal(dim(s), c(np, np))
  expect_true(all(is.finite(s)))
  expect_equal(unname(s), unname(t(s)), tolerance = 1e-8)
  expect_equal(unname(s), unname(manual_sandwich(m)), tolerance = 1e-6)
  expect_equal(rownames(s), names(coef(m)))
  expect_equal(colnames(s), names(coef(m)))
  expect_equal(colnames(sandwich::estfun(m)), names(coef(m)))
  expect_equal(rownames(sandwich::bread(m)), names(coef(m)))
}

test_that("estfun and bread are registered S3 methods for class 'drc'", {
  expect_true(is.function(getS3method("estfun", "drc", optional = TRUE)))
  expect_true(is.function(getS3method("bread", "drc", optional = TRUE)))
  # the direct-call exports are kept for backward compatibility
  expect_true(all(c("estfun.drc", "bread.drc") %in% getNamespaceExports("drc")))
})

test_that("sandwich generics work without attaching sandwich (continuous)", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_sandwich_ok(m)
  expect_equal(unname(sandwich::estfun(m)), unname(estfun.drc(m)))
  expect_equal(unname(sandwich::bread(m)), unname(bread.drc(m)))
})

test_that("sandwich generics work with fixed parameters", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  expect_sandwich_ok(m)
})

test_that("sandwich generics work for multiple curves", {
  mc <- data.frame(
    resp = c(ryegrass$rootl, ryegrass$rootl * 0.8),
    dose = rep(ryegrass$conc, 2),
    g = factor(rep(c("A", "B"), each = nrow(ryegrass)))
  )
  expect_sandwich_ok(drm(resp ~ dose, curveid = g, data = mc, fct = LL.4()))
  expect_sandwich_ok(drm(resp ~ dose, curveid = g, pmodels = data.frame(g, 1, 1, 1),
                         data = mc, fct = LL.4()))
})

test_that("sandwich generics work for binomial and Poisson models", {
  mb <- drm(dead / total ~ conc, weights = total, data = selenium,
            fct = LL.2(), type = "binomial")
  expect_sandwich_ok(mb)

  set.seed(42)
  pd <- data.frame(dose = rep(c(0, 0.1, 0.5, 1, 2, 5, 10), each = 3),
                   count = rpois(21, rep(c(20, 18, 15, 10, 5, 2, 1), each = 3)))
  expect_sandwich_ok(drm(count ~ dose, data = pd, fct = LL.4(), type = "Poisson"))
})

test_that("sandwich generics work for weighted fits and hormesis models", {
  w <- rep(c(1, 2), length.out = nrow(ryegrass))
  expect_sandwich_ok(drm(rootl ~ conc, data = ryegrass, fct = LL.4(), weights = w))
  expect_sandwich_ok(drm(weight ~ conc, data = lettuce, fct = BC.4()))
})

test_that("sandwich generics work when only one parameter is free (issue #50 setting)", {
  d <- data.frame(
    conc = rep(c(0, 0.1, 1, 10, 100, 1000), each = 4),
    yield = c(rep(100, 4), 96, 98, 94, 99, 85, 88, 80, 90, 55, 60, 50, 58, 20, 24, 15, 22, 5, 3, 8, 6)
  )
  m <- drm(yield ~ conc, data = d, fct = EXD.3(fixed = c(0, 100, NA)))
  expect_equal(length(coef(m)), 1)
  expect_sandwich_ok(m)
  ed <- ED(m, c(10, 50), vcov. = sandwich::sandwich, display = FALSE)
  expect_true(all(is.finite(ed[, "Std. Error"])))
})

test_that("sandwich generics work for event-time data", {
  m <- drm(count ~ start + end, data = chickweed, fct = LL.3(), type = "event")
  ef <- sandwich::estfun(m)
  expect_equal(ncol(ef), length(coef(m)))
  expect_equal(unname(ef), unname(estfun.drc(m)))
})

test_that("ED() accepts sandwich as vcov. and gives robust standard errors", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  ed_model <- ED(m, c(10, 50), display = FALSE)
  ed_robust <- ED(m, c(10, 50), vcov. = sandwich::sandwich, display = FALSE)

  expect_equal(ed_robust[, "Estimate"], ed_model[, "Estimate"])
  expect_true(all(is.finite(ed_robust[, "Std. Error"])))
  expect_false(isTRUE(all.equal(ed_robust[, "Std. Error"], ed_model[, "Std. Error"])))

  # a precomputed matrix gives the same answer as the function
  ed_mat <- ED(m, c(10, 50), vcov. = sandwich::sandwich(m), display = FALSE)
  expect_equal(ed_mat, ed_robust)
})

test_that("vcovCL() works for drc models", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  v <- sandwich::vcovCL(m, cluster = seq_len(nrow(ryegrass)))
  expect_equal(dim(v), rep(length(coef(m)), 2))
  expect_true(all(is.finite(v)))
})

test_that("every <generic>.drc function in the namespace is a registered S3 method", {
  ns <- asNamespace("drc")
  reg <- getNamespaceInfo(ns, "S3methods")
  registered <- paste(reg[, 1], reg[, 2], sep = ".")
  candidates <- grep("\\.(drc|drclist|mrdrc)$", ls(ns, all.names = TRUE), value = TRUE)
  # internal helpers that are deliberately not S3 methods: no object has class "drclist"
  # (anova.drclist is called directly by anova.drc) and hat.drc is called directly by mrdrm
  not_methods <- c("hat.drc", "anova.drclist")
  missing <- setdiff(candidates, c(registered, not_methods))
  expect_equal(missing, character(0))
})

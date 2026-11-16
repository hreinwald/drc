# Numerical cross-checks of analytic derivatives, ED values and ED standard errors.
# Guards against wrong gradients / ED conventions in the model functions
# (logistic, Gompertz, LL2, W2x, Cedergreen-Ritz-Streibig, ...).

# ---- helpers ---------------------------------------------------------------

# central finite differences of a vector-valued function f(p); returns n x length(p)
num_jac <- function(f, p) {
  cols <- lapply(seq_along(p), function(i) {
    h <- 1e-6 * max(1, abs(p[i]))
    up <- p; dn <- p
    up[i] <- up[i] + h; dn[i] <- dn[i] - h
    (f(up) - f(dn)) / (2 * h)
  })
  do.call(cbind, cols)
}

make_model <- function(expr) suppressWarnings(eval(parse(text = expr)))

# parameter vector with plausible values; b gets the requested sign
model_parm <- function(m, bsign = 1) {
  base <- c(b = 1.5, c = 0.5, d = 8, e = 2, f = 0.8)
  nm <- m$names
  v <- vapply(nm, function(n) if (n %in% names(base)) base[[n]] else 1, numeric(1))
  if ("b" %in% nm) v["b"] <- bsign * abs(v["b"])
  unname(v)
}

# fixed, deterministic, positive-definite "covariance matrix"
fake_vcov <- function(np) {
  A <- outer(seq_len(np), seq_len(np), function(i, j) sin(i * j + 1))
  crossprod(A) + diag(np)
}

se_from_grad <- function(g, V) sqrt(drop(t(g) %*% V %*% g))

# ---- deriv1 / derivx versus finite differences --------------------------------

deriv_models <- c(
  "gammadr()",
  "LL.2()", "LL.3()", "LL.3u()", "LL.4()", "LL.5()", "MM.2()", "MM.3()",
  "W1.2()", "W1.3()", "W1.3u()", "W1.4()", "EXD.2()", "EXD.3()",
  "W2.2()", "W2.3()", "W2.3u()", "W2.4()", "AR.2()", "AR.3()",
  "LN.2()", "LN.3()", "LN.3u()", "LN.4()",
  "LL2.2()", "LL2.3()", "LL2.3u()", "LL2.4()", "LL2.5()",
  "G.2()", "G.3()", "G.3u()", "G.4()", "L.3()", "L.4()", "L.5()",
  "BC.4()", "BC.5()",
  "CRS.4a()", "CRS.4b()", "CRS.4c()", "CRS.5a()", "CRS.5b()", "CRS.5c()", "CRS.5()",
  "UCRS.4a()", "UCRS.4b()", "UCRS.4c()", "UCRS.5a()", "UCRS.5b()", "UCRS.5c()",
  "FPL.4(-1, 1)", "gaussian()", "lgaussian()"
)

for (expr in deriv_models) {
  local({
    expr <- expr
    test_that(paste("deriv1/derivx agree with finite differences:", expr), {
      m <- make_model(expr)
      dose <- c(0.1, 0.3, 0.7, 1.5, 3, 6, 12)
      np <- length(m$names)
      parm <- model_parm(m, bsign = if (grepl("^L\\.", expr)) -1 else 1)
      pm <- matrix(parm, length(dose), np, byrow = TRUE)

      if (is.function(m$deriv1)) {
        analytic <- matrix(m$deriv1(dose, pm), length(dose), np)
        numeric <- num_jac(function(q) m$fct(dose, matrix(q, length(dose), np, byrow = TRUE)), parm)
        relerr <- apply(abs(analytic - numeric), 2, max) / pmax(1e-8, apply(abs(numeric), 2, max))
        expect_true(all(relerr < 1e-4), info = paste("parameter(s):", paste(m$names[relerr >= 1e-4], collapse = ", ")))
      }
      if (is.function(m$derivx)) {
        h <- 1e-6
        numx <- (m$fct(dose + h, pm) - m$fct(dose - h, pm)) / (2 * h)
        expect_lt(max(abs(m$derivx(dose, pm) - numx)) / max(1e-8, max(abs(numx))), 1e-4)
      }
    })
  })
}

# ---- ED values: f(ED) must reproduce the requested effect ------------------------
# c and d are read from the parameter vector (c = 0 / d = 1 when not free).
# Decreasing curve: f(ED_p) = d - p/100 (d - c); increasing: c + p/100 (d - c).

ed_value_models <- c(
  "LL.4()", "LL.5()", "W1.4()", "W2.4()", "LN.4()",
  "L.3()", "L.4()", "L.5()",
  "G.2()", "G.3()", "G.4()",
  "LL2.2()", "LL2.3()", "LL2.4()", "LL2.5()",
  "W2x.3()", "W2x.4()"
)
hormesis_models <- c("BC.4()", "BC.5()", "CRS.4a()", "CRS.5a()")

check_ed_values <- function(expr, bsigns, tol) {
  m <- make_model(expr)
  np <- length(m$names)
  logED <- grepl("^LL2", expr)   # LL2 returns log(ED)
  for (bs in bsigns) {
    parm <- model_parm(m, bs)
    cc <- if ("c" %in% m$names) parm[match("c", m$names)] else 0
    dd <- if ("d" %in% m$names) parm[match("d", m$names)] else 1
    fx <- function(x) m$fct(x, matrix(parm, 1, np, byrow = TRUE))
    decreasing <- fx(1e-3) > fx(1e3)
    for (p in c(10, 50, 90)) {
      ed <- m$edfct(parm, p, reference = "control", type = "relative")[[1]]
      if (logED) ed <- exp(ed)
      expected <- if (decreasing) dd - p / 100 * (dd - cc) else cc + p / 100 * (dd - cc)
      expect_lt(abs(fx(ed) - expected), tol * abs(dd - cc),
                label = sprintf("%s b-sign %+d relative p=%d", expr, bs, p))
    }
    respl <- cc + 0.3 * (dd - cc)
    eda <- m$edfct(parm, respl, reference = "control", type = "absolute")[[1]]
    if (logED) eda <- exp(eda)
    expect_lt(abs(fx(eda) - respl), tol * abs(dd - cc),
              label = sprintf("%s b-sign %+d absolute", expr, bs))
  }
}

for (expr in ed_value_models) {
  local({
    expr <- expr
    test_that(paste("ED reproduces the effect level:", expr), {
      m <- make_model(expr)
      bsigns <- if ("b" %in% m$names) c(1, -1) else 1
      check_ed_values(expr, bsigns, tol = 1e-8)
    })
  })
}

for (expr in hormesis_models) {
  local({
    expr <- expr
    test_that(paste("ED reproduces the effect level (hormesis, root finding):", expr), {
      check_ed_values(expr, 1, tol = 1e-3)
    })
  })
}

# ---- ED gradient (-> standard error) versus finite differences ---------------------

se_models <- c(
  "LL.2()", "LL.3()", "LL.3u()", "LL.4()", "LL.5()", "MM.2()", "MM.3()",
  "W1.2()", "W1.3()", "W1.3u()", "W1.4()", "EXD.2()", "EXD.3()",
  "W2.2()", "W2.3()", "W2.3u()", "W2.4()", "AR.2()", "AR.3()", "W2x.3()", "W2x.4()",
  "LN.2()", "LN.3()", "LN.3u()", "LN.4()",
  "LL2.2()", "LL2.3()", "LL2.3u()", "LL2.4()", "LL2.5()",
  "G.2()", "G.3()", "G.3u()", "G.4()", "L.3()", "L.4()", "L.5()",
  "BC.4()", "BC.5()", "CRS.4a()", "CRS.4b()", "CRS.5a()", "CRS.5b()",
  "UCRS.4a()", "UCRS.5a()", "UCRS.5b()", "FPL.4(-1, 1)"
)
root_finding <- "^(BC|CRS|UCRS|FPL)"

for (expr in se_models) {
  local({
    expr <- expr
    test_that(paste("ED standard error matches numerical delta method:", expr), {
      m <- make_model(expr)
      np <- length(m$names)
      V <- fake_vcov(np)
      tol <- if (grepl(root_finding, expr)) 2e-2 else 1e-4
      checked <- 0
      for (bs in c(1, -1)) {
        parm <- model_parm(m, bs)
        for (p in c(20, 50, 80)) {
          res <- tryCatch(suppressWarnings(m$edfct(parm, p, reference = "control", type = "relative")),
                          error = function(e) NULL)
          if (is.null(res) || !is.finite(res[[1]])) next
          expect_length(res[[2]], np)
          edf <- function(q) suppressWarnings(m$edfct(q, p, reference = "control", type = "relative")[[1]])
          g_num <- as.vector(num_jac(edf, parm))
          if (!all(is.finite(g_num))) next
          expect_equal(se_from_grad(as.numeric(res[[2]]), V), se_from_grad(g_num, V), tolerance = tol,
                       info = sprintf("%s b-sign %+d p=%d", expr, bs, p))
          checked <- checked + 1
        }
      }
      expect_gt(checked, 0)
    })
  })
}

# ---- fitted models: ED() output, regression tests for the reported defects -----------

num_ed_se <- function(fit, p) {
  f <- function(q) fit$fct$edfct(q, p, reference = "control", type = "relative")[[1]]
  se_from_grad(as.vector(num_jac(f, coef(fit))), vcov(fit))
}

test_that("L.4 ED follows the log-logistic convention on decreasing data", {
  mLL <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  mL <- drm(rootl ~ conc, data = ryegrass, fct = L.4())
  p <- c(10, 50, 90)
  edL <- ED(mL, p, display = FALSE)[, 1]

  expect_true(all(diff(edL) > 0))   # ED10 < ED50 < ED90 for a decreasing curve
  cf <- coef(mL)
  expected <- cf[3] - p / 100 * (cf[3] - cf[2])
  expect_equal(as.vector(predict(mL, data.frame(conc = edL))), unname(expected), tolerance = 1e-6)
  # logistic and log-logistic differ in shape, but a swapped ED10/ED90 would be far off
  expect_equal(edL, ED(mLL, p, display = FALSE)[, 1], tolerance = 0.5)

  # absolute levels
  edAbs <- ED(mL, 4, type = "absolute", display = FALSE)[, 1]
  expect_equal(as.vector(predict(mL, data.frame(conc = edAbs))), 4, tolerance = 1e-6)

  # the S3 object keeps the usual structure
  expect_identical(names(mL), names(mLL))
  expect_s3_class(mL, "drc")
})

test_that("ED standard errors of fitted models match the numerical delta method", {
  fits <- list(
    LL.4 = drm(rootl ~ conc, data = ryegrass, fct = LL.4()),
    L.4 = drm(rootl ~ conc, data = ryegrass, fct = L.4()),
    G.4 = drm(rootl ~ conc, data = ryegrass, fct = G.4()),
    LL2.5 = drm(rootl ~ conc, data = ryegrass, fct = LL2.5())
  )
  for (nm in names(fits)) {
    for (p in c(10, 50, 90)) {
      expect_equal(ED(fits[[nm]], p, display = FALSE)[1, 2], num_ed_se(fits[[nm]], p),
                   tolerance = 1e-4, info = paste(nm, "p =", p))
    }
  }
})

test_that("W2x ED returns a finite estimate and standard error", {
  x <- rep(c(1.5, 2, 3, 4, 6, 8, 12, 16, 24), each = 3)
  y <- 10 * (1 - exp(-exp(log(x - 1) - log(4)))) + 0.2 * sin(seq_along(x) * 1.7)
  mW <- drm(resp ~ dose, data = data.frame(dose = x, resp = y), fct = W2x.3())
  ed <- mW$fct$edfct(coef(mW), 50, reference = "control", type = "relative")
  expect_length(ed[[2]], length(coef(mW)))
  expect_equal(ed[[2]][length(ed[[2]])], 1)   # d ED / d t0

  res <- ED(mW, 50, display = FALSE)
  expect_true(all(is.finite(res[1, 1:2])))
  expect_equal(res[1, 2], num_ed_se(mW, 50), tolerance = 1e-4)
  expect_false(is.na(ED_robust(mW, 50)$Estimate))
})

test_that("predict() standard errors for CRS models match the numerical delta method", {
  mC <- suppressWarnings(drm(weight ~ conc, data = lettuce, fct = CRS.4a()))
  nd <- data.frame(conc = c(0.5, 2, 10))
  se_fit <- predict(mC, nd, se.fit = TRUE)[, "SE"]
  pc <- coef(mC)
  se_num <- vapply(nd$conc, function(x) {
    g <- as.vector(num_jac(function(q) mC$fct$fct(x, matrix(q, 1, length(q))), pc))
    se_from_grad(g, vcov(mC))
  }, numeric(1))
  expect_equal(unname(se_fit), se_num, tolerance = 1e-4)
})

# ---- separate fits must not depend on an attached 'magic' package -------------------

test_that("blockDiag builds block-diagonal matrices", {
  a <- matrix(1:4, 2, 2); b <- matrix(5, 1, 1); d <- matrix(6:11, 2, 3)
  out <- blockDiag(list(a, b, d))
  expect_equal(dim(out), c(5, 6))
  expect_equal(out[1:2, 1:2], a)
  expect_equal(out[3, 3], 5)
  expect_equal(out[4:5, 4:6], d)
  expect_equal(sum(out != 0), sum(a != 0) + 1 + sum(d != 0))
})

test_that("vcov/summary/confint/ED/predict work for separate fits", {
  ms <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL.4(), separate = TRUE)
  np <- length(coef(ms))

  v <- vcov(ms, pool = FALSE)
  expect_equal(dim(v), c(np, np))
  k <- length(coef(ms$objList[[1]]))
  expect_equal(v[seq_len(k), seq_len(k)], vcov(ms$objList[[1]]), ignore_attr = TRUE)
  expect_equal(v[1:k, (k + 1):(2 * k)], matrix(0, k, k), ignore_attr = TRUE)
  expect_true(all(diag(vcov(ms)) > 0))

  expect_s3_class(summary(ms), "summary.drc")
  expect_equal(nrow(confint(ms)), np)
  expect_equal(nrow(ED(ms, 50, display = FALSE)), length(ms$objList))
  expect_length(predict(ms), nrow(spinach))
})

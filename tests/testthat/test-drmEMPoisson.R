# Tests for the Poisson log-likelihood (R/drmEMPoisson.R)

poisson_data <- function(scale = 1) {
  set.seed(42)
  data.frame(
    dose = rep(c(0, 0.1, 0.5, 1, 2, 5, 10), each = 3),
    count = rpois(21, rep(c(20, 18, 15, 10, 5, 2, 1), each = 3) * scale)
  )
}

test_that("logLik() of a Poisson drm fit equals the Poisson log-density at the fitted values", {
  pd <- poisson_data()
  m <- drm(count ~ dose, data = pd, fct = LL.4(), type = "Poisson")

  ll <- as.numeric(logLik(m))
  expect_equal(ll, sum(dpois(pd$count, fitted(m), log = TRUE)), tolerance = 1e-8)
  expect_lt(ll, 0)
  expect_equal(AIC(m), -2 * ll + 2 * length(coef(m)), tolerance = 1e-8)
})

test_that("Poisson logLik stays finite for large counts (no overflow of gamma())", {
  pd <- poisson_data(scale = 20)
  m <- drm(count ~ dose, data = pd, fct = LL.4(), type = "Poisson")

  expect_true(max(pd$count) > 171)
  expect_true(is.finite(as.numeric(logLik(m))))
  expect_equal(as.numeric(logLik(m)), sum(dpois(pd$count, fitted(m), log = TRUE)), tolerance = 1e-8)
})

test_that("logLik() df is the number of estimated parameters (Poisson, negbin, event)", {
  pd <- poisson_data()
  mp <- drm(count ~ dose, data = pd, fct = LL.4(), type = "Poisson")
  expect_equal(attr(logLik(mp), "df"), length(coef(mp)))
  expect_equal(AIC(mp), -2 * as.numeric(logLik(mp)) + 2 * length(coef(mp)), tolerance = 1e-8)

  mn <- drm(count ~ dose, data = pd, fct = LL.4(), type = "negbin2")
  expect_equal(attr(logLik(mn), "df"), length(coef(mn)))

  me <- drm(count ~ start + end, data = chickweed, fct = LL.3(), type = "event")
  expect_equal(attr(logLik(me), "df"), length(coef(me)))
})

test_that("Poisson estimates maximise the Poisson likelihood", {
  pd <- poisson_data()
  m <- drm(count ~ dose, data = pd, fct = LL.4(), type = "Poisson")
  nll <- function(p) {
    mu <- ifelse(pd$dose == 0, p[3],
                 p[2] + (p[3] - p[2]) / (1 + exp(p[1] * (log(pd$dose) - log(p[4])))))
    if (any(!is.finite(mu)) || any(mu <= 0)) return(1e10)
    -sum(dpois(pd$count, mu, log = TRUE))
  }
  opt <- suppressWarnings(optim(coef(m), nll, method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-14)))
  expect_gte(as.numeric(logLik(m)), -opt$value - 1e-2)
})

# Compare drm() output to drm_legacy() output
# drm() must produce an identically structured and valued output as drm_legacy()

# Helper function to compare two drc fits element by element
compare_drc_fits <- function(fit1, fit2, tol = .Machine$double.eps^0.5) {
  # Check that names match
  expect_identical(names(fit1), names(fit2))

  # Compare all non-function elements
  for (nm in names(fit1)) {
    # 'call' naturally differs because the function names are different
    if (nm == "call") next

    el1 <- fit1[[nm]]
    el2 <- fit2[[nm]]

    if (is.function(el1)) next

    if (is.list(el1) && !is.data.frame(el1)) {
      # Compare list sub-elements
      expect_identical(names(el1), names(el2), info = paste("List names for", nm))
      for (subnm in names(el1)) {
        if (is.function(el1[[subnm]])) next
        expect_equal(el1[[subnm]], el2[[subnm]],
          tolerance = tol,
          info = paste(nm, "$", subnm, sep = "")
        )
      }
    } else {
      expect_equal(el1, el2, tolerance = tol, info = nm)
    }
  }
}

# Create test datasets
ryegrass <- data.frame(
  rootl = c(
    7.58, 8.00, 8.33, 7.25, 7.17, 7.00, 7.17, 7.83, 7.92, 7.58,
    6.17, 5.75, 5.83, 6.00, 5.83, 4.92, 4.50, 4.17, 4.42, 4.00,
    2.67, 2.08, 2.42, 2.50, 2.25, 1.17, 0.75, 0.92, 1.00, 0.58
  ),
  conc = c(
    rep(0, 5), rep(0.94, 5), rep(1.88, 5),
    rep(3.75, 5), rep(7.50, 5), rep(15, 5)
  )
)

set.seed(42)
multi_data <- data.frame(
  dose = rep(c(0, 0.5, 1, 2, 5, 10), each = 5, times = 2),
  resp = c(
    rnorm(5, 100, 5), rnorm(5, 95, 5), rnorm(5, 85, 5),
    rnorm(5, 60, 5), rnorm(5, 20, 5), rnorm(5, 5, 5),
    rnorm(5, 100, 5), rnorm(5, 90, 5), rnorm(5, 70, 5),
    rnorm(5, 40, 5), rnorm(5, 10, 5), rnorm(5, 3, 5)
  ),
  group = rep(c("A", "B"), each = 30)
)

test_that("drm() matches drm_legacy() for single curve LL.4", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4())
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() for multi-curve LL.4", {
  fit1 <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  fit2 <- drc:::drm_legacy(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with Box-Cox transformation", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), bcVal = 0)
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4(), bcVal = 0)
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with logDose", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), logDose = 10)
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4(), logDose = 10)
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with Weibull W1.4", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = W1.4())
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with pmodels", {
  fit1 <- drm(resp ~ dose,
    curveid = group, data = multi_data, fct = LL.4(),
    pmodels = data.frame(multi_data$group, multi_data$group, 1, multi_data$group)
  )
  fit2 <- drc:::drm_legacy(resp ~ dose,
    curveid = group, data = multi_data, fct = LL.4(),
    pmodels = data.frame(multi_data$group, multi_data$group, 1, multi_data$group)
  )
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with constrained optimization", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), lowerl = c(-Inf, 0, -Inf, -Inf))
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4(), lowerl = c(-Inf, 0, -Inf, -Inf))
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with binomial type", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  fit1 <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  fit2 <- drc:::drm_legacy(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with Poisson type", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  fit1 <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")
  fit2 <- drc:::drm_legacy(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with negbin1 type", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  fit1 <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "negbin1")
  fit2 <- drc:::drm_legacy(count ~ dose, data = poisson_data, fct = LL.4(), type = "negbin1")
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with negbin2 type", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  fit1 <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "negbin2")
  fit2 <- drc:::drm_legacy(count ~ dose, data = poisson_data, fct = LL.4(), type = "negbin2")
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with robust estimation (median)", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), robust = "median")
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4(), robust = "median")
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with LL.3 model", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.3())
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with user-supplied start values", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), start = c(2.5, 1, 0.5, 8))
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4(), start = c(2.5, 1, 0.5, 8))
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with varcov matrix", {
  vcov_mat <- diag(30) * 2
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), varcov = vcov_mat)
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4(), varcov = vcov_mat)
  compare_drc_fits(fit1, fit2)
})

test_that("drm() matches drm_legacy() with noMessage=TRUE", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), control = drmc(noMessage = TRUE))
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4(), control = drmc(noMessage = TRUE))
  compare_drc_fits(fit1, fit2)
})

test_that("drm() return value has correct class", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_identical(class(fit1), class(fit2))
  expect_identical(class(fit1), "drc")
})

test_that("drm() return value has all expected named elements", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expected_names <- c(
    "varParm", "fit", "curve", "summary", "start", "parNames",
    "predres", "call", "data",
    "parmMat", "fct", "robust", "estMethod", "df.residual",
    "sumList", "scaleFct", "pmFct", "pfFct", "type", "indexMat",
    "logDose", "cm", "deriv1",
    "curveVarNam", "origData", "weights",
    "dataList", "coefficients", "boxcox", "indexMat2"
  )
  expect_identical(names(fit1), expected_names)
})

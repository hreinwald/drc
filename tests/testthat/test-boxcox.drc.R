# Tests for boxcox.drc(), boxcoxCI(), and anovaFormula()

# =========================================================================
# Tests for boxcox.drc()
# =========================================================================

# --- Method "ml" ---

test_that("boxcox.drc works with method='ml' and plotit=TRUE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  pdf(file = tempfile(fileext = ".pdf"))
  result <- boxcox(m1, lambda = seq(-2, 2, by = 0.5), plotit = TRUE, method = "ml")
  dev.off()

  expect_true(inherits(result, "drc"))
  expect_true(is.list(result$boxcox))
  expect_true(is.numeric(result$boxcox$lambda))
  expect_equal(length(result$boxcox$lambda), 1)
  expect_equal(length(result$boxcox$ci), 2)
  expect_equal(result$boxcox$bcAdd, 0)
  expect_false(is.null(result$call$bcVal))
  expect_false(is.null(result$call$bcAdd))
})

test_that("boxcox.drc works with method='ml' and plotit=FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- boxcox(m1, lambda = seq(-1, 2, by = 0.5), plotit = FALSE, method = "ml")

  expect_true(inherits(result, "drc"))
  expect_true(is.list(result$boxcox))
  expect_true(is.numeric(result$boxcox$lambda))
})

test_that("boxcox.drc ml method with bcAdd", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- boxcox(m1, lambda = seq(0, 2, by = 0.5), plotit = FALSE, bcAdd = 1)

  expect_true(inherits(result, "drc"))
  expect_equal(result$boxcox$bcAdd, 1)
  expect_equal(result$call$bcAdd, 1)
})

test_that("boxcox.drc ml method with custom level", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- boxcox(m1, lambda = seq(-1, 2, by = 0.5), plotit = FALSE, level = 0.90)

  expect_true(inherits(result, "drc"))
  expect_true(is.list(result$boxcox))
})

# --- Method "fixed" ---

test_that("boxcox.drc works with a single lambda (fixed method)", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- boxcox(m1, lambda = 1)

  expect_true(inherits(result, "drc"))
  expect_equal(result$boxcox$lambda, 1)
  expect_true(all(is.na(result$boxcox$ci)))
  expect_equal(result$call$bcVal, 1)
})

test_that("boxcox.drc fixed method with lambda = 0", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- boxcox(m1, lambda = 0)

  expect_true(inherits(result, "drc"))
  expect_equal(result$boxcox$lambda, 0)
  expect_equal(result$boxcox$ci, c(NA, NA))
})

# --- Method "anova" ---

test_that("boxcox.drc works with method='anova' single curve", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  pdf(file = tempfile(fileext = ".pdf"))
  result <- boxcox(m1, lambda = seq(-2, 2, by = 0.5), method = "anova", plotit = TRUE)
  dev.off()

  expect_true(inherits(result, "drc"))
  expect_true(is.list(result$boxcox))
  expect_true(is.numeric(result$boxcox$lambda))
  expect_equal(length(result$boxcox$ci), 2)
})

test_that("boxcox.drc works with method='anova' and plotit=FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- boxcox(m1, lambda = seq(-2, 2, by = 0.5), method = "anova", plotit = FALSE)

  expect_true(inherits(result, "drc"))
  expect_true(is.list(result$boxcox))
})

test_that("boxcox.drc works with method='anova' multi-curve", {
  m1 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = S.alba, fct = LL.4())
  pdf(file = tempfile(fileext = ".pdf"))
  result <- boxcox(m1, lambda = seq(0, 2, by = 0.5), method = "anova", plotit = TRUE)
  dev.off()

  expect_true(inherits(result, "drc"))
  expect_true(is.list(result$boxcox))
})

test_that("boxcox.drc anova method errors without replicates", {
  # Create a dataset with no replicates (each dose has exactly 1 obs)
  no_rep_data <- data.frame(
    dose = c(0, 1, 2, 5, 10, 20),
    resp = c(100, 90, 70, 40, 10, 3)
  )
  m1 <- drm(resp ~ dose, data = no_rep_data, fct = LL.4())

  expect_error(
    boxcox(m1, method = "anova"),
    "ANOVA-based TBS approach requires replicates for each dose value"
  )
})

# --- Return value structure ---

test_that("boxcox.drc returns invisible drc object with boxcox component", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- boxcox(m1, lambda = 1)
  expect_true(inherits(result, "drc"))
  expect_true("boxcox" %in% names(result))
  expect_true(is.list(result$boxcox))
  expect_true(all(c("lambda", "ci", "bcAdd") %in% names(result$boxcox)))
})

# =========================================================================
# Tests for boxcoxCI()
# =========================================================================

test_that("boxcoxCI computes correct confidence interval", {
  # Create a synthetic bell-shaped log-likelihood
  x <- seq(-2, 2, by = 0.25)
  y <- -((x - 0.5)^2) + 10  # peak at x = 0.5

  ci <- drc:::boxcoxCI(x, y, level = 0.95)

  expect_equal(length(ci), 2)
  expect_true(is.numeric(ci))
  # CI should bracket the optimal value of 0.5
  expect_true(ci[1] < 0.5)
  expect_true(ci[2] > 0.5)
})

test_that("boxcoxCI handles NA values in y", {
  x <- seq(-2, 2, by = 0.25)
  y <- -((x - 0.5)^2) + 10
  y[c(1, 2, 16, 17)] <- NA  # add some NAs

  ci <- drc:::boxcoxCI(x, y, level = 0.95)

  expect_equal(length(ci), 2)
  expect_true(is.numeric(ci))
})

test_that("boxcoxCI returns NA for lower bound when max is at left edge", {
  # Log-likelihood that is monotonically decreasing
  x <- seq(-2, 2, by = 0.25)
  y <- -x + 10  # max at x = -2 (left edge), loglik[1] >= lim
  # Here loglik[1] = max, so loglik[1] is NOT < lim => xx[1] stays NA

  ci <- drc:::boxcoxCI(x, y, level = 0.95)

  expect_true(is.na(ci[1]))  # lower bound NA because max is at edge
})

test_that("boxcoxCI returns NA for upper bound when max is at right edge", {
  # Log-likelihood that is monotonically increasing
  x <- seq(-2, 2, by = 0.25)
  y <- x + 10  # max at x = 2 (right edge), loglik[m] >= lim
  # Here loglik[m] = max, so loglik[m] is NOT < lim => xx[2] stays NA

  ci <- drc:::boxcoxCI(x, y, level = 0.95)

  expect_true(is.na(ci[2]))  # upper bound NA because max is at edge
})

test_that("boxcoxCI with different confidence levels", {
  # Use wide grid centered at 0 to ensure both CI bounds are captured
  x <- seq(-5, 5, by = 0.25)
  y <- -((x - 0)^2) + 10

  ci_95 <- drc:::boxcoxCI(x, y, level = 0.95)
  ci_99 <- drc:::boxcoxCI(x, y, level = 0.99)

  # 99% CI should be wider than 95% CI
  expect_true((ci_99[2] - ci_99[1]) >= (ci_95[2] - ci_95[1]))
})

# =========================================================================
# Tests for anovaFormula()
# =========================================================================

test_that("anovaFormula creates formula for single curve", {
  dose <- c(0, 1, 2, 5, 10)
  resp <- c(100, 80, 60, 30, 5)
  curveid <- rep("A", 5)
  bcAdd <- 0

  result <- drc:::anovaFormula(dose, resp, curveid, bcAdd)

  expect_true(is.list(result))
  expect_true("anovaFormula" %in% names(result))
  expect_true("anovaData" %in% names(result))
  expect_true(inherits(result$anovaFormula, "formula"))
  expect_true(is.data.frame(result$anovaData))
  expect_equal(nrow(result$anovaData), 5)
  expect_true(all(c("dose", "resp", "curveid", "bcc") %in% names(result$anovaData)))
  # Single curve: should NOT have interaction term
  formula_str <- deparse(result$anovaFormula)
  expect_false(grepl("\\*", paste(formula_str, collapse = " ")))
})

test_that("anovaFormula creates formula for multiple curves", {
  dose <- c(0, 1, 2, 5, 10, 0, 1, 2, 5, 10)
  resp <- c(100, 80, 60, 30, 5, 95, 75, 55, 25, 3)
  curveid <- rep(c("A", "B"), each = 5)
  bcAdd <- 0

  result <- drc:::anovaFormula(dose, resp, curveid, bcAdd)

  expect_true(is.list(result))
  expect_true(inherits(result$anovaFormula, "formula"))
  expect_equal(nrow(result$anovaData), 10)
  # Multi curve: should have interaction term
  formula_str <- deparse(result$anovaFormula)
  expect_true(grepl("\\*", paste(formula_str, collapse = " ")))
})

test_that("anovaFormula with non-zero bcAdd", {
  dose <- c(0, 1, 2)
  resp <- c(100, 80, 60)
  curveid <- rep("A", 3)
  bcAdd <- 2

  result <- drc:::anovaFormula(dose, resp, curveid, bcAdd)

  expect_true(all(result$anovaData$bcc == 2))
  expect_equal(nrow(result$anovaData), 3)
})

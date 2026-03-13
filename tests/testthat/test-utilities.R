# Test additional utility and model functions

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

# Tests for update.drc()

test_that("update.drc can update model function", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m2 <- update(m1, fct = LL.3())

  expect_true(inherits(m2, "drc"))
  expect_equal(length(coef(m2)), 3)
  expect_equal(length(coef(m1)), 4)
})

test_that("update.drc can update data", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Create modified dataset
  ryegrass_subset <- ryegrass[1:20, ]
  m2 <- update(m1, data = ryegrass_subset)

  expect_true(inherits(m2, "drc"))
  expect_equal(nrow(m2$data), 20)
})

test_that("update.drc preserves original model when nothing changed", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m2 <- update(m1)

  expect_equal(coef(m1), coef(m2))
})

# Tests for logLik.drc()

test_that("logLik.drc returns log-likelihood value", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  ll <- logLik(m1)

  expect_true(is.numeric(ll))
  expect_equal(length(ll), 1)
  expect_true(ll < 0)  # Log-likelihood is typically negative
})

test_that("logLik.drc has df attribute", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  ll <- logLik(m1)

  expect_true("df" %in% names(attributes(ll)))
  # df includes model parameters plus the error variance parameter
  expect_equal(attr(ll, "df"), length(coef(m1)) + 1)
})

test_that("logLik.drc works with different model types", {
  m_ll3 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  m_ll4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  ll3 <- logLik(m_ll3)
  ll4 <- logLik(m_ll4)

  # LL.4 should have higher log-likelihood (better fit with more parameters)
  expect_true(ll4 > ll3)
})

# Tests for anova.drc()

test_that("anova.drc compares two models", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  result <- anova(m1, m2)

  expect_true(is.data.frame(result) || is.matrix(result))
  expect_true(nrow(result) >= 2)
})

test_that("anova.drc returns valid p-values", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  result <- anova(m1, m2, details = FALSE)

  # Find p-value column
  p_col <- which(grepl("p-value|Pr\\(", colnames(result), ignore.case = TRUE))
  if (length(p_col) > 0) {
    p_values <- result[, p_col]
    p_values <- p_values[!is.na(p_values)]
    if (length(p_values) > 0) {
      expect_true(all(p_values >= 0 & p_values <= 1))
    }
  }
})

test_that("anova.drc with single model directs to modelFit", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_error(anova(m1, test = "Chi"), "modelFit")
})

test_that("modelFit performs lack-of-fit test for single model", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  result <- modelFit(m1)

  expect_true(is.data.frame(result) || inherits(result, "anova"))
})

# Tests for df.residual()

test_that("df.residual returns correct degrees of freedom", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  df <- df.residual(m1)

  expect_true(is.numeric(df))
  expect_equal(length(df), 1)
  expect_true(df > 0)
  expect_equal(df, nrow(ryegrass) - length(coef(m1)))
})

# Tests for AIC and BIC

test_that("AIC works with drc models", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  aic_val <- AIC(m1)

  expect_true(is.numeric(aic_val))
  expect_equal(length(aic_val), 1)
})

test_that("BIC works with drc models", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  bic_val <- BIC(m1)

  expect_true(is.numeric(bic_val))
  expect_equal(length(bic_val), 1)
})

test_that("AIC and BIC can compare models", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  aic1 <- AIC(m1)
  aic2 <- AIC(m2)

  bic1 <- BIC(m1)
  bic2 <- BIC(m2)

  # Both should be numeric
  expect_true(is.numeric(aic1) && is.numeric(aic2))
  expect_true(is.numeric(bic1) && is.numeric(bic2))
})

# Tests for boxcox.drc()

test_that("boxcox.drc works with Box-Cox transformation", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), bcVal = 0)

  expect_true(!is.null(m1$boxcox))
  expect_true(is.list(m1$boxcox))
})

test_that("boxcox.drc model can make predictions", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), bcVal = 0)

  predictions <- predict(m1)
  expect_true(length(predictions) == nrow(ryegrass))
  expect_true(all(is.finite(predictions)))
})

# Tests for modelFit() and related functions

test_that("drm creates valid model object structure", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Check essential components
  expect_true("coefficients" %in% names(m1))
  expect_true("fct" %in% names(m1))
  expect_true("fit" %in% names(m1))
  expect_true("data" %in% names(m1))
  expect_true("call" %in% names(m1))
})

test_that("drm with start values", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
            start = c(2, 0.5, 8, 3))

  expect_true(inherits(m1, "drc"))
  expect_equal(length(coef(m1)), 4)
})

test_that("drm with weights", {
  weights <- rep(1, nrow(ryegrass))
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
            weights = weights)

  expect_true(inherits(m1, "drc"))
  expect_true(!is.null(m1$weights))
})

test_that("drm with subset", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
            subset = conc < 10)

  expect_true(inherits(m1, "drc"))
  expect_true(nrow(m1$data) < nrow(ryegrass))
})

# Tests for formula handling

test_that("drm formula with different specifications", {
  # Standard formula
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  expect_true(inherits(m1, "drc"))

  # With curveid
  m2 <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())
  expect_true(inherits(m2, "drc"))
})

# Tests for different dose-response functions

test_that("drm works with log-logistic functions", {
  m_ll2 <- drm(rootl ~ conc, data = ryegrass, fct = LL.2())
  m_ll3 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  m_ll4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m_ll5 <- drm(rootl ~ conc, data = ryegrass, fct = LL.5())

  expect_equal(length(coef(m_ll2)), 2)
  expect_equal(length(coef(m_ll3)), 3)
  expect_equal(length(coef(m_ll4)), 4)
  expect_equal(length(coef(m_ll5)), 5)
})

test_that("drm works with Weibull functions", {
  m_w12 <- drm(rootl ~ conc, data = ryegrass, fct = W1.2())
  m_w13 <- drm(rootl ~ conc, data = ryegrass, fct = W1.3())
  m_w14 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())

  expect_equal(length(coef(m_w12)), 2)
  expect_equal(length(coef(m_w13)), 3)
  expect_equal(length(coef(m_w14)), 4)
})

test_that("drm works with two-parameter Weibull", {
  m_w22 <- drm(rootl ~ conc, data = ryegrass, fct = W2.2())
  m_w23 <- drm(rootl ~ conc, data = ryegrass, fct = W2.3())
  m_w24 <- drm(rootl ~ conc, data = ryegrass, fct = W2.4())

  expect_equal(length(coef(m_w22)), 2)
  expect_equal(length(coef(m_w23)), 3)
  expect_equal(length(coef(m_w24)), 4)
})

# Tests for control parameters

test_that("drm with drmc control", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
            control = drmc(errorm = FALSE))

  expect_true(inherits(m1, "drc"))
})

test_that("drm with maxIt control", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
            control = drmc(maxIt = 500))

  expect_true(inherits(m1, "drc"))
})

# Tests for robust estimation

test_that("drm with robust = 'median'", {
  m_robust <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
                  robust = "median")

  expect_true(inherits(m_robust, "drc"))
  expect_true(!is.null(m_robust$robust))
  expect_equal(m_robust$robust, "median")
})

# Tests for different data types

test_that("drm handles binomial data correctly", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(),
                 type = "binomial", weights = n)

  expect_true(inherits(m_binom, "drc"))
  expect_equal(m_binom$type, "binomial")
})

test_that("drm handles Poisson data correctly", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  m_poisson <- drm(count ~ dose, data = poisson_data, fct = LL.4(),
                   type = "Poisson")

  expect_true(inherits(m_poisson, "drc"))
  expect_equal(m_poisson$type, "Poisson")
})

# Tests for model validation functions

test_that("drm model has proper class", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_true(inherits(m1, "drc"))
  expect_equal(class(m1), "drc")
})

test_that("drm coefficients are named", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  coeffs <- coef(m1)

  expect_true(!is.null(names(coeffs)))
  expect_equal(length(names(coeffs)), length(coeffs))
})

# Tests for hatvalues (leverage)

test_that("hatvalues.drc returns leverage values", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Check if hatvalues function exists for drc
  if (exists("hatvalues.drc") || "hatvalues" %in% methods(class = "drc")) {
    hv <- hatvalues(m1)
    expect_true(is.numeric(hv))
    expect_equal(length(hv), nrow(ryegrass))
    expect_true(all(hv >= 0 & hv <= 1))
  }
})

# Tests for cooks.distance

test_that("cooks.distance.drc returns influence measures", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Check if cooks.distance function exists for drc
  if (exists("cooks.distance.drc") || "cooks.distance" %in% methods(class = "drc")) {
    cd <- cooks.distance(m1)
    expect_true(is.numeric(cd))
    expect_equal(length(cd), nrow(ryegrass))
    expect_true(all(cd >= 0))
  }
})

# Integration tests

test_that("Complete workflow: fit, summarize, predict, plot", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Summarize
  summ <- summary(m1)
  expect_true(is.list(summ))

  # Predict
  newdata <- data.frame(conc = c(1, 2, 3))
  pred <- predict(m1, newdata = newdata)
  expect_equal(length(pred), 3)

  # ED values
  ed50 <- ED(m1, 50, display = FALSE)
  expect_true(is.matrix(ed50))

  # Plot
  expect_no_error({
    pdf(file = tempfile())
    plot(m1)
    dev.off()
  })
})

test_that("Model comparison workflow", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Compare AIC
  aic1 <- AIC(m1)
  aic2 <- AIC(m2)
  expect_true(is.numeric(aic1) && is.numeric(aic2))

  # ANOVA comparison
  comp <- anova(m1, m2)
  expect_true(nrow(comp) >= 2)

  # Log-likelihood comparison
  ll1 <- logLik(m1)
  ll2 <- logLik(m2)
  expect_true(ll2 > ll1)  # More parameters should give better fit
})

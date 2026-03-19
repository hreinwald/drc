# Tests for lin.test() function
# Lack-of-fit test for the mean structure based on cumulated residuals

# ---- Setup: create reusable model objects ----

# Model WITH replicates (ryegrass: 24 obs, 7 unique doses)
ryegrass_data <- data.frame(
  rootl = c(7.58, 8.0, 8.33, 7.25, 7.37, 7.96, 8.36, 6.91, 7.75,
            6.87, 6.45, 5.92, 1.93, 2.89, 4.23, 1.19, 0.86, 1.06,
            0.69, 0.52, 0.82, 0.25, 0.22, 0.44),
  conc = c(0, 0, 0, 0, 0, 0, 0.94, 0.94, 0.94,
           1.88, 1.88, 1.88, 3.75, 3.75, 3.75, 7.5, 7.5, 7.5,
           15, 15, 15, 30, 30, 30)
)

# Model WITHOUT replicates (unique x values only)
norep_data <- data.frame(
  resp = c(7.5, 7.0, 6.5, 5.5, 4.0, 2.5, 1.5, 0.8, 0.4, 0.2),
  dose = c(0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)
)


# ===========================================================
# Test block 1: Basic functionality with replicates
# ===========================================================
test_that("lin.test returns numeric p-value with replicate data", {
  m_rep <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  pval <- lin.test(m_rep, noksSim = 20, plotit = FALSE)
  expect_true(is.numeric(pval))
  expect_length(pval, 1)
  expect_true(pval >= 0 && pval <= 1)
})


# ===========================================================
# Test block 2: Basic functionality without replicates
# ===========================================================
test_that("lin.test returns numeric p-value without replicate data", {
  m_norep <- drm(resp ~ dose, data = norep_data, fct = LL.4())
  pval <- lin.test(m_norep, noksSim = 20, plotit = FALSE)
  expect_true(is.numeric(pval))
  expect_length(pval, 1)
  expect_true(pval >= 0 && pval <= 1)
})


# ===========================================================
# Test block 3: Reproducibility via seed
# ===========================================================
test_that("lin.test produces reproducible results with same seed", {
  m_rep <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  pval1 <- lin.test(m_rep, noksSim = 50, seed = 42, plotit = FALSE)
  pval2 <- lin.test(m_rep, noksSim = 50, seed = 42, plotit = FALSE)
  expect_identical(pval1, pval2)
})


# ===========================================================
# Test block 4: seed = NULL path
# ===========================================================
test_that("lin.test works with seed = NULL", {
  m_rep <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  pval <- lin.test(m_rep, noksSim = 20, seed = NULL, plotit = FALSE)
  expect_true(is.numeric(pval))
  expect_length(pval, 1)
  expect_true(pval >= 0 && pval <= 1)
})


# ===========================================================
# Test block 5: Plotting with default parameters (log="", missing ylim/xlab/ylab)
# Covers: plotit=TRUE, log="" (else branch), missing(ylim), missing(xlab), missing(ylab)
# ===========================================================
test_that("lin.test plots with default parameters (replicates)", {
  m_rep <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  pdf(NULL)  # suppress graphical output
  on.exit(dev.off(), add = TRUE)
  pval <- lin.test(m_rep, noksSim = 10, plotit = TRUE)
  expect_true(is.numeric(pval))
})


# ===========================================================
# Test block 6: Plotting with log="x"
# Covers: if (identical(log, "x")) branch
# ===========================================================
test_that("lin.test plots with log='x' scale (replicates)", {
  m_rep <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  pval <- lin.test(m_rep, noksSim = 10, plotit = TRUE, log = "x")
  expect_true(is.numeric(pval))
})


# ===========================================================
# Test block 7: Plotting with custom ylim
# Covers: else branch of if (missing(ylim))
# ===========================================================
test_that("lin.test plots with custom ylim", {
  m_rep <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  pval <- lin.test(m_rep, noksSim = 10, plotit = TRUE, ylim = c(-5, 5))
  expect_true(is.numeric(pval))
})


# ===========================================================
# Test block 8: Plotting with custom xlab and ylab
# Covers: non-missing xlab and ylab paths in ifelse()
# ===========================================================
test_that("lin.test plots with custom xlab and ylab", {
  m_rep <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  pval <- lin.test(m_rep, noksSim = 10, plotit = TRUE,
                   xlab = "Dose (mg/L)", ylab = "Cum. Residuals")
  expect_true(is.numeric(pval))
})


# ===========================================================
# Test block 9: Plotting without replicates (repAdjust=FALSE + plotit=TRUE)
# ===========================================================
test_that("lin.test plots without replicates (repAdjust=FALSE)", {
  m_norep <- drm(resp ~ dose, data = norep_data, fct = LL.4())
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  pval <- lin.test(m_norep, noksSim = 10, plotit = TRUE)
  expect_true(is.numeric(pval))
})


# ===========================================================
# Test block 10: Plotting without replicates with log="x"
# ===========================================================
test_that("lin.test plots without replicates with log='x'", {
  m_norep <- drm(resp ~ dose, data = norep_data, fct = LL.4())
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  pval <- lin.test(m_norep, noksSim = 10, plotit = TRUE, log = "x")
  expect_true(is.numeric(pval))
})


# ===========================================================
# Test block 11: All optional plot params together
# Covers: custom ylim + xlab + ylab + log="x" combined
# ===========================================================
test_that("lin.test plots with all custom parameters combined", {
  m_rep <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  pval <- lin.test(m_rep, noksSim = 10, plotit = TRUE,
                   log = "x", ylim = c(-3, 3),
                   xlab = "Concentration", ylab = "Residuals")
  expect_true(is.numeric(pval))
})


# ===========================================================
# Test block 12: Different noksSim values
# ===========================================================
test_that("lin.test works with different noksSim values", {
  m_rep <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  pval <- lin.test(m_rep, noksSim = 5, plotit = FALSE)
  expect_true(is.numeric(pval))
  expect_length(pval, 1)
})

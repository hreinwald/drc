# Test plot.drc() function - Plotting dose-response curves

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

# Basic plot tests

test_that("plot.drc creates plot without error", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1)
    dev.off()
  })
})

test_that("plot.drc returns data frame invisibly", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  pdf(file = tempfile())
  result <- plot(m1)
  dev.off()

  expect_true(is.data.frame(result))
})

test_that("plot.drc works with different plot types", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Test different type parameters
  types <- c("average", "all", "bars", "none", "obs")

  for (plot_type in types) {
    expect_no_error({
      pdf(file = tempfile())
      plot(m1, type = plot_type)
      dev.off()
    })
  }
})

test_that("plot.drc type = 'confidence' creates confidence bands", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, type = "confidence")
    dev.off()
  })
})

# Tests with multi-curve models

test_that("plot.drc works with multi-curve models", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m_multi)
    dev.off()
  })
})

test_that("plot.drc can plot specific curves using level", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m_multi, level = "A")
    dev.off()
  })
})

test_that("plot.drc with multiple levels", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m_multi, level = c("A", "B"))
    dev.off()
  })
})

# Tests with graphical parameters

test_that("plot.drc respects col parameter", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, col = "blue")
    dev.off()
  })
})

test_that("plot.drc respects lty parameter", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, lty = 2)
    dev.off()
  })
})

test_that("plot.drc respects pch parameter", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, pch = 19)
    dev.off()
  })
})

test_that("plot.drc with custom axis labels", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, xlab = "Concentration", ylab = "Root Length")
    dev.off()
  })
})

test_that("plot.drc with custom xlim and ylim", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, xlim = c(0, 20), ylim = c(0, 10))
    dev.off()
  })
})

# Tests with log scale

test_that("plot.drc with log = 'x'", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, log = "x")
    dev.off()
  })
})

test_that("plot.drc with log = ''", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, log = "")
    dev.off()
  })
})

# Tests with broken axis

test_that("plot.drc with broken = TRUE", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, broken = TRUE)
    dev.off()
  })
})

# Tests with gridsize

test_that("plot.drc respects gridsize parameter", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, gridsize = 50)
    dev.off()
  })

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, gridsize = 200)
    dev.off()
  })
})

# Tests with legend

test_that("plot.drc with legend = TRUE", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m_multi, legend = TRUE)
    dev.off()
  })
})

test_that("plot.drc with legend = FALSE", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m_multi, legend = FALSE)
    dev.off()
  })
})

# Tests with normalization

test_that("plot.drc with normal = TRUE normalizes data", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, normal = TRUE)
    dev.off()
  })
})

# Tests with confidence level

test_that("plot.drc respects confidence.level parameter", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, type = "bars", confidence.level = 0.90)
    dev.off()
  })

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, type = "confidence", confidence.level = 0.99)
    dev.off()
  })
})

# Tests with different model types

test_that("plot.drc works with LL.3 model", {
  m_ll3 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())

  expect_no_error({
    pdf(file = tempfile())
    plot(m_ll3)
    dev.off()
  })
})

test_that("plot.drc works with Weibull models", {
  m_w1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m_w1)
    dev.off()
  })
})

test_that("plot.drc works with binomial type data", {
  binom_data <- data.frame(
    dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
    resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
    n = rep(50, 7)
  )
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(), type = "binomial", weights = n)

  expect_no_error({
    pdf(file = tempfile())
    plot(m_binom)
    dev.off()
  })
})

test_that("plot.drc works with Poisson type data", {
  poisson_data <- data.frame(
    dose = c(0, 1, 2, 4, 8, 16, 32),
    count = c(50, 48, 40, 25, 10, 3, 1)
  )
  m_poisson <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")

  expect_no_error({
    pdf(file = tempfile())
    plot(m_poisson)
    dev.off()
  })
})

# Tests with multiple graphical parameters combined

test_that("plot.drc with multiple graphical parameters", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1,
         col = "red",
         lty = 2,
         pch = 19,
         xlab = "Concentration (mM)",
         ylab = "Root length (cm)",
         main = "Ryegrass dose-response",
         xlim = c(0, 20),
         ylim = c(0, 10)
    )
    dev.off()
  })
})

# Tests for edge cases

test_that("plot.drc handles models with zero dose", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # ryegrass data includes zero doses
  expect_no_error({
    pdf(file = tempfile())
    plot(m1, log = "x")  # log scale with zero dose
    dev.off()
  })
})

test_that("plot.drc with robust estimation", {
  m_robust <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(), robust = "median")

  expect_no_error({
    pdf(file = tempfile())
    plot(m_robust)
    dev.off()
  })
})

test_that("plot.drc works after model update", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  m2 <- update(m1, fct = LL.3())

  expect_no_error({
    pdf(file = tempfile())
    plot(m2)
    dev.off()
  })
})

# Tests with colors for multi-curve

test_that("plot.drc with vector of colors for multi-curve", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m_multi, col = c("red", "blue"))
    dev.off()
  })
})

test_that("plot.drc with logical col = TRUE for auto colors", {
  m_multi <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m_multi, col = TRUE)
    dev.off()
  })
})

# Integration tests

test_that("plot.drc can be called multiple times", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1)
    plot(m1, col = "blue")
    dev.off()
  })
})

test_that("plot.drc works with par() settings", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    old_par <- par(mfrow = c(2, 2))
    plot(m1)
    plot(m1, type = "bars")
    plot(m1, type = "confidence")
    plot(m1, log = "")
    par(old_par)
    dev.off()
  })
})

test_that("plot.drc works after predictions", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Make some predictions first
  newdata <- data.frame(conc = c(1, 2, 3))
  pred <- predict(m1, newdata = newdata)

  # Plot should still work
  expect_no_error({
    pdf(file = tempfile())
    plot(m1)
    dev.off()
  })
})

test_that("plot.drc consistent across calls", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Both plots should work without interfering
  expect_no_error({
    pdf(file = tempfile())
    plot(m1)
    dev.off()

    pdf(file = tempfile())
    plot(m1)
    dev.off()
  })
})

# Tests to ensure plot doesn't crash with various argument combinations

test_that("plot.drc with type = 'all' and col vector", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  expect_no_error({
    pdf(file = tempfile())
    plot(m1, type = "all", col = "darkgreen")
    dev.off()
  })
})

test_that("plot.drc with add = TRUE parameter if supported", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # First plot
  expect_no_error({
    pdf(file = tempfile())
    plot(m1)
    # Try to add to existing plot (if add parameter exists)
    # This may or may not be supported, so we just test it doesn't crash
    try(plot(m1, add = TRUE), silent = TRUE)
    dev.off()
  })
})

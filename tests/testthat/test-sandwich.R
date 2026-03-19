# Tests for estfun.drc() and bread.drc() from R/sandwich.R

# ============================================================================
# Setup: Create test data and models for different types
# ============================================================================

# --- Continuous type (single curve) ---
test_that("estfun.drc works for continuous type (single curve)", {
  data(ryegrass)
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  ef <- estfun.drc(m1)

  # Should return a matrix
  expect_true(is.matrix(ef))

  # Number of rows should match number of observations
  expect_equal(nrow(ef), nrow(ryegrass))

  # Number of columns should match number of parameters
  expect_equal(ncol(ef), length(coef(m1)))

  # Column names should match coefficient names
  expect_equal(colnames(ef), names(coef(m1)))

  # Values should be numeric and finite
  expect_true(all(is.finite(ef)))
})

# --- Continuous type (multi-curve) ---
test_that("estfun.drc works for continuous type (multi-curve)", {
  data(ryegrass)
  # Create multi-curve data
  multi_data <- data.frame(
    resp = c(ryegrass$rootl, ryegrass$rootl * 0.8),
    dose = c(ryegrass$conc, ryegrass$conc),
    group = factor(rep(c("A", "B"), each = nrow(ryegrass)))
  )
  m2 <- drm(resp ~ dose, curveid = group, data = multi_data, fct = LL.4())

  ef <- estfun.drc(m2)

  expect_true(is.matrix(ef))
  expect_equal(nrow(ef), nrow(multi_data))
  expect_equal(ncol(ef), length(coef(m2)))
  expect_equal(colnames(ef), names(coef(m2)))
})

# --- Binomial type ---
test_that("estfun.drc works for binomial type", {
  data(selenium)
  m3 <- drm(dead / total ~ conc, weights = total, data = selenium,
            fct = LL.2(), type = "binomial")

  ef <- estfun.drc(m3)

  expect_true(is.matrix(ef))
  expect_equal(nrow(ef), nrow(selenium))
  expect_equal(ncol(ef), length(coef(m3)))
  expect_equal(colnames(ef), names(coef(m3)))
})

# --- Poisson type ---
test_that("estfun.drc works for Poisson type", {
  # Create Poisson-distributed data
  set.seed(42)
  poisson_data <- data.frame(
    dose = rep(c(0, 0.1, 0.5, 1, 2, 5, 10), each = 3),
    count = c(rpois(3, 20), rpois(3, 18), rpois(3, 15),
              rpois(3, 10), rpois(3, 5), rpois(3, 2), rpois(3, 1))
  )
  m4 <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")

  ef <- estfun.drc(m4)

  expect_true(is.matrix(ef))
  expect_equal(nrow(ef), nrow(poisson_data))
  expect_equal(ncol(ef), length(coef(m4)))
  expect_equal(colnames(ef), names(coef(m4)))
})

# --- Binomial type (multi-curve, tests is.matrix(indexMat0) == TRUE branch) ---
test_that("estfun.drc works for binomial type (multi-curve)", {
  data(selenium)
  multi_binom <- data.frame(
    dead = c(selenium$dead, selenium$dead),
    total = c(selenium$total, selenium$total),
    conc = c(selenium$conc, selenium$conc),
    group = factor(rep(c("A", "B"), each = nrow(selenium)))
  )
  m5 <- drm(dead / total ~ conc, curveid = group, weights = total,
            data = multi_binom, fct = LL.2(), type = "binomial")

  ef <- estfun.drc(m5)

  expect_true(is.matrix(ef))
  expect_equal(nrow(ef), nrow(multi_binom))
  expect_equal(ncol(ef), length(coef(m5)))
  expect_equal(colnames(ef), names(coef(m5)))
})

# --- Poisson type (multi-curve) ---
test_that("estfun.drc works for Poisson type (multi-curve)", {
  set.seed(42)
  poisson_data_mc <- data.frame(
    dose = rep(rep(c(0, 0.5, 1, 5, 10), each = 3), 2),
    count = c(rpois(15, c(rep(20, 3), rep(10, 3), rep(5, 3), rep(2, 3), rep(1, 3))),
              rpois(15, c(rep(15, 3), rep(8, 3), rep(4, 3), rep(2, 3), rep(1, 3)))),
    group = factor(rep(c("A", "B"), each = 15))
  )
  m6 <- drm(count ~ dose, curveid = group, data = poisson_data_mc,
            fct = LL.4(), type = "Poisson")

  ef <- estfun.drc(m6)

  expect_true(is.matrix(ef))
  expect_equal(nrow(ef), nrow(poisson_data_mc))
  expect_equal(ncol(ef), length(coef(m6)))
  expect_equal(colnames(ef), names(coef(m6)))
})

# ============================================================================
# bread.drc tests
# ============================================================================

test_that("bread.drc works for continuous type", {
  data(ryegrass)
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  br <- bread.drc(m1)

  # Should return a matrix
  expect_true(is.matrix(br))

  # Should be square with dimension = number of parameters
  np <- length(coef(m1))
  expect_equal(dim(br), c(np, np))

  # Should be finite
  expect_true(all(is.finite(br)))
})

test_that("bread.drc works for binomial type (non-continuous path)", {
  data(selenium)
  m3 <- drm(dead / total ~ conc, weights = total, data = selenium,
            fct = LL.2(), type = "binomial")

  br <- bread.drc(m3)

  expect_true(is.matrix(br))
  np <- length(coef(m3))
  expect_equal(dim(br), c(np, np))
  expect_true(all(is.finite(br)))
})

test_that("bread.drc works for Poisson type (non-continuous path)", {
  set.seed(42)
  poisson_data <- data.frame(
    dose = rep(c(0, 0.1, 0.5, 1, 2, 5, 10), each = 3),
    count = c(rpois(3, 20), rpois(3, 18), rpois(3, 15),
              rpois(3, 10), rpois(3, 5), rpois(3, 2), rpois(3, 1))
  )
  m4 <- drm(count ~ dose, data = poisson_data, fct = LL.4(), type = "Poisson")

  br <- bread.drc(m4)

  expect_true(is.matrix(br))
  np <- length(coef(m4))
  expect_equal(dim(br), c(np, np))
})

# ============================================================================
# Integration: bread and estfun work together
# ============================================================================

test_that("bread and estfun are consistent for continuous model", {
  data(ryegrass)
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  ef <- estfun.drc(m1)
  br <- bread.drc(m1)

  # Both should have compatible dimensions
  expect_equal(ncol(ef), nrow(br))
  expect_equal(ncol(ef), ncol(br))
})

test_that("bread and estfun are consistent for binomial model", {
  data(selenium)
  m3 <- drm(dead / total ~ conc, weights = total, data = selenium,
            fct = LL.2(), type = "binomial")

  ef <- estfun.drc(m3)
  br <- bread.drc(m3)

  expect_equal(ncol(ef), nrow(br))
  expect_equal(ncol(ef), ncol(br))
})

# ============================================================================
# Edge case: single curve with indexMat0 not being a matrix (vector branch)
# ============================================================================

test_that("estfun.drc handles non-matrix indexMat0 (single curve)", {
  data(ryegrass)
  # Single curve model - indexMat2 may be a vector
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  # Verify indexMat2 structure for single curve
  # This should exercise the is.matrix() == FALSE branch
  ef <- estfun.drc(m1)
  expect_true(is.matrix(ef))
  expect_equal(ncol(ef), length(coef(m1)))
})

# ============================================================================
# Event type
# ============================================================================

test_that("estfun.drc works for event type (single curve)", {
  data(chickweed)
  m_event <- drm(count ~ start + end, data = chickweed, fct = LL.3(), type = "event")

  ef <- estfun.drc(m_event)

  expect_true(is.matrix(ef))
  expect_equal(nrow(ef), nrow(m_event$data))
  expect_equal(ncol(ef), length(coef(m_event)))
  expect_equal(colnames(ef), names(coef(m_event)))
})

test_that("estfun.drc works for event type (multi-curve)", {
  data(germination)
  germ_sub <- germination[germination$species == "wheat" &
                            germination$temp %in% c(10, 22), ]
  m_event_mc <- drm(germinated ~ start + end, factor(temp),
                    data = germ_sub, fct = LL.3(), type = "event")

  ef <- estfun.drc(m_event_mc)

  expect_true(is.matrix(ef))
  expect_equal(nrow(ef), nrow(m_event_mc$data))
  expect_equal(ncol(ef), length(coef(m_event_mc)))
  expect_equal(colnames(ef), names(coef(m_event_mc)))
})

test_that("bread.drc works for event type (non-continuous path)", {
  data(chickweed)
  m_event <- drm(count ~ start + end, data = chickweed, fct = LL.3(), type = "event")

  br <- bread.drc(m_event)

  expect_true(is.matrix(br))
  np <- length(coef(m_event))
  expect_equal(dim(br), c(np, np))
})

# ============================================================================
# Correctness checks
# ============================================================================

test_that("estfun.drc returns correct values structure for continuous", {
  data(ryegrass)
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  ef <- estfun.drc(m1)

  # Sum of estimating functions should be close to zero at MLE
  col_sums <- colSums(ef)
  expect_true(all(is.numeric(col_sums)))
  expect_true(all(abs(col_sums) < 1))
})

test_that("bread.drc is consistent with vcov for continuous", {
  data(ryegrass)
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  br <- bread.drc(m1)
  vc <- vcov(m1)

  # Both should be square matrices of same size

  expect_equal(dim(br), dim(vc))
})

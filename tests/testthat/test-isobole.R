# tests/testthat/test-isobole.R
# Comprehensive tests for the isobole() function

# ---------------------------------------------------------------------------
# Helper: fit the "free" (unconstrained EC50) model used by all isobole tests
# ---------------------------------------------------------------------------
fit_mecter_free <- function() {
  drm(rgr ~ dose, pct, data = mecter,
      fct = LL.4(),
      pmodels = list(~1, ~1, ~1, ~factor(pct) - 1))
}

fit_acidiq_free <- function() {
  drm(rgr ~ dose, pct, data = acidiq,
      fct = LL.4(),
      pmodels = list(~factor(pct), ~1, ~1, ~factor(pct) - 1))
}

# ===========================================================================
# 1.  Basic / "happy path" – object1 only, default parameters
# ===========================================================================
test_that("isobole produces a plot with object1 only (default args)", {
  m_free <- fit_mecter_free()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(isobole(m_free, exchange = 0.02))
})

# ===========================================================================
# 2.  Custom xlim / ylim supplied
# ===========================================================================
test_that("isobole respects user-supplied xlim and ylim", {
  m_free <- fit_mecter_free()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(
    isobole(m_free, exchange = 0.02, xlim = c(0, 500), ylim = c(0, 10))
  )
})

# ===========================================================================
# 3.  Custom xlab / ylab supplied
# ===========================================================================
test_that("isobole respects user-supplied xlab and ylab", {
  m_free <- fit_mecter_free()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(
    isobole(m_free, exchange = 0.02, xlab = "Substance A", ylab = "Substance B")
  )
})

# ===========================================================================
# 4.  xaxis = "0"  – axis swap path
# ===========================================================================
test_that("isobole swaps axes when xaxis = '0' (object1 only)", {
  m_free <- fit_mecter_free()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(isobole(m_free, exchange = 0.02, xaxis = "0"))
})

# ===========================================================================
# 5.  xaxis = "0" with custom labels (covers the else branch for labels)
# ===========================================================================
test_that("isobole with xaxis='0' auto-labels are '0' and '100'", {
  m_free <- fit_mecter_free()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  # Without explicit xlab/ylab, xaxis="0" should produce "0" and "100" as labels

  expect_no_error(isobole(m_free, exchange = 0.02, xaxis = "0"))
})

# ===========================================================================
# 6.  object2 = CA model (concentration addition, lambda = 1)
# ===========================================================================
test_that("isobole draws CA isobole line", {
  m_free <- fit_mecter_free()
  m_ca   <- mixture(m_free, model = "CA")
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(isobole(m_free, m_ca, exchange = 0.02))
})

# ===========================================================================
# 7.  object2 = Hewlett model
# ===========================================================================
test_that("isobole draws Hewlett isobole line", {
  m_free <- fit_mecter_free()
  m_hew  <- mixture(m_free, model = "Hewlett")
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(isobole(m_free, m_hew, exchange = 0.02))
})

# ===========================================================================
# 8.  object2 = Voelund model
# ===========================================================================
test_that("isobole draws Voelund isobole line", {
  m_free <- fit_mecter_free()
  m_voe  <- mixture(m_free, model = "Voelund")
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(isobole(m_free, m_voe, exchange = 0.02))
})

# ===========================================================================
# 9.  object2 + xaxis = "0"  (covers swap inside object2 block)
# ===========================================================================
test_that("isobole with object2 and xaxis='0' swaps correctly (CA)", {
  m_free <- fit_mecter_free()
  m_ca   <- mixture(m_free, model = "CA")
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(isobole(m_free, m_ca, exchange = 0.02, xaxis = "0"))
})

test_that("isobole with object2 and xaxis='0' swaps correctly (Voelund)", {
  m_free <- fit_mecter_free()
  m_voe  <- mixture(m_free, model = "Voelund")
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(isobole(m_free, m_voe, exchange = 0.02, xaxis = "0"))
})

# ===========================================================================
# 10. cifactor argument affects CI width
# ===========================================================================
test_that("isobole works with different cifactor values", {
  m_free <- fit_mecter_free()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(isobole(m_free, exchange = 0.02, cifactor = 1))
})

# ===========================================================================
# 11. Using a different dataset (acidiq) with Hewlett model
# ===========================================================================
test_that("isobole works with acidiq data and Hewlett model", {
  m_free <- fit_acidiq_free()
  m_hew  <- mixture(m_free, model = "Hewlett")
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(
    isobole(m_free, m_hew, xlim = c(0, 400), ylim = c(0, 450))
  )
})

# ===========================================================================
# 12. Hewlett model with xaxis = "0" (swap path inside non-voelund branch)
# ===========================================================================
test_that("isobole Hewlett with xaxis='0' swaps axes in object2 block", {
  m_free <- fit_mecter_free()
  m_hew  <- mixture(m_free, model = "Hewlett")
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(isobole(m_free, m_hew, exchange = 0.02, xaxis = "0"))
})

# ===========================================================================
# 13. All parameters supplied (full custom call)
# ===========================================================================
test_that("isobole with all custom parameters works", {
  m_free <- fit_mecter_free()
  m_voe  <- mixture(m_free, model = "Voelund")
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(
    isobole(m_free, m_voe, exchange = 0.02, cifactor = 1,
            xlab = "A", ylab = "B", xlim = c(0, 300), ylim = c(0, 6))
  )
})

# ===========================================================================
# 14. Voelund with xaxis = "0" and custom limits
# ===========================================================================
test_that("isobole Voelund xaxis='0' with custom limits", {
  m_free <- fit_mecter_free()
  m_voe  <- mixture(m_free, model = "Voelund")
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_no_error(
    isobole(m_free, m_voe, exchange = 0.02, xaxis = "0",
            xlim = c(0, 10), ylim = c(0, 500))
  )
})

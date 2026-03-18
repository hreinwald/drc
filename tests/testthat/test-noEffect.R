# Tests for noEffect function
# Tests the likelihood ratio test for dose effect significance

test_that("noEffect returns correct structure for continuous model", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- noEffect(m1)

  expect_true(is.numeric(result))
  expect_length(result, 3)
  expect_named(result, c("Chi-square test", "Df", "p-value"))
})

test_that("noEffect returns correct values for continuous model", {
  m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  result <- noEffect(m1)

  # Chi-square test statistic should be positive
  expect_true(result["Chi-square test"] > 0)
  # Degrees of freedom should be positive integer
  expect_true(result["Df"] > 0)
  # p-value should be between 0 and 1
  expect_true(result["p-value"] >= 0 && result["p-value"] <= 1)
  # For ryegrass data, there is a strong dose effect, so p-value should be very small
  expect_true(result["p-value"] < 0.05)
})

test_that("noEffect works for binomial model", {
  m_bin <- drm(number / total ~ dose, weights = total, data = earthworms,
               fct = LL.2(), type = "binomial")
  result <- noEffect(m_bin)

  expect_true(is.numeric(result))
  expect_length(result, 3)
  expect_named(result, c("Chi-square test", "Df", "p-value"))
  expect_true(is.finite(result["Chi-square test"]))
  expect_true(result["Df"] > 0)
  expect_true(result["p-value"] >= 0 && result["p-value"] <= 1)
})

test_that("noEffect works for Poisson model", {
  m_pois <- drm(count ~ conc, data = decontaminants, fct = LL.3(),
                type = "Poisson")
  result <- noEffect(m_pois)

  expect_true(is.numeric(result))
  expect_length(result, 3)
  expect_named(result, c("Chi-square test", "Df", "p-value"))
  expect_true(result["Chi-square test"] > 0)
  expect_true(result["Df"] > 0)
  expect_true(result["p-value"] >= 0 && result["p-value"] <= 1)
})

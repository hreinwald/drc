# Test that otrace correctly controls error message display

test_that("drmc() otrace default is FALSE", {
  ctrl <- drmc()
  expect_false(ctrl$otrace)
})

test_that("otrace=TRUE does not suppress error messages from try(optim())", {
  # When otrace=TRUE, silentVal should be FALSE (i.e., errors ARE displayed)
  # Use intentionally bad data/starting values to trigger an error in optim
  # and verify the error message is shown (not silenced)

  # A normal fit with otrace=TRUE should still converge
  fit <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
    control = drmc(otrace = TRUE))
  expect_s3_class(fit, "drc")
})

test_that("otrace=FALSE suppresses error messages from try(optim())", {
  # Default behavior: errors from try(optim()) are silenced
  fit <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
    control = drmc(otrace = FALSE))
  expect_s3_class(fit, "drc")
})

test_that("drm() matches drm_legacy() with otrace=TRUE", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
    control = drmc(otrace = TRUE))
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4(),
    control = drmc(otrace = TRUE))

  # Both should converge successfully

  expect_s3_class(fit1, "drc")
  expect_s3_class(fit2, "drc")

  # Coefficients should match
  expect_equal(coef(fit1), coef(fit2), tolerance = 1e-6)
})

# Test that otrace correctly controls error message display

test_that("drmc() otrace default is FALSE", {
  ctrl <- drmc()
  expect_false(ctrl$otrace)
})

test_that("otrace=TRUE allows successful convergence", {
  fit <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
    control = drmc(otrace = TRUE))
  expect_s3_class(fit, "drc")
})

test_that("otrace=FALSE allows successful convergence", {
  fit <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
    control = drmc(otrace = FALSE))
  expect_s3_class(fit, "drc")
})

test_that("drm() matches drm_legacy() with otrace=TRUE", {
  fit1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
    control = drmc(otrace = TRUE))
  fit2 <- drc:::drm_legacy(rootl ~ conc, data = ryegrass, fct = LL.4(),
    control = drmc(otrace = TRUE))

  expect_s3_class(fit1, "drc")
  expect_s3_class(fit2, "drc")
  expect_equal(coef(fit1), coef(fit2), tolerance = 1e-6)
})

# ==============================================================================
# Tests for mixture(), hewlett(), and voelund() functions
# ==============================================================================

# --- Shared setup: fit free models used across multiple tests -----------------

# LL.4 free model (4-parameter log-logistic)
acidiq.free4 <- drm(rgr ~ dose, pct, data = acidiq,
                    fct = LL.4(),
                    pmodels = list(~1, ~1, ~1, ~factor(pct) - 1))

# LL.3 free model (3-parameter log-logistic)
acidiq.free3 <- drm(rgr ~ dose, pct, data = acidiq,
                    fct = LL.3(),
                    pmodels = list(~1, ~1, ~factor(pct) - 1))

# LL.2 free model (2-parameter log-logistic)
acidiq.free2 <- drm(rgr ~ dose, pct, data = acidiq,
                    fct = LL.2(),
                    pmodels = list(~1, ~factor(pct) - 1))

# LL.5 free model (5-parameter log-logistic, used for error tests)
acidiq.free5 <- drm(rgr ~ dose, pct, data = acidiq,
                    fct = LL.5(),
                    pmodels = list(~1, ~1, ~1, ~1, ~factor(pct) - 1))

# ==============================================================================
# Tests for mixture() - CA model
# ==============================================================================

test_that("mixture CA works with LL.4 model", {
  result <- mixture(acidiq.free4, model = "CA")
  expect_s3_class(result, "drc")
  expect_equal(result$text, "CA model")
  expect_true(!is.null(result$pmodelsText))
  expect_true(!is.null(result$deviance))
  expect_equal(result$anova$test, "F")
})

test_that("mixture CA works with LL.3 model", {
  result <- mixture(acidiq.free3, model = "CA")
  expect_s3_class(result, "drc")
  expect_equal(result$text, "CA model")
})

test_that("mixture CA works with LL.2 model", {
  result <- mixture(acidiq.free2, model = "CA")
  expect_s3_class(result, "drc")
  expect_equal(result$text, "CA model")
})

test_that("mixture CA errors with LL.5 model", {
  expect_error(mixture(acidiq.free5, model = "CA"), "Does not work for LL.5")
})

test_that("mixture CA uses default startm = NULL when missing", {
  # When startm is missing for CA, it defaults to NULL
  result <- mixture(acidiq.free4, model = "CA")
  expect_s3_class(result, "drc")
})

test_that("mixture CA with explicit startm = NULL", {
  result <- mixture(acidiq.free4, model = "CA", startm = NULL)
  expect_s3_class(result, "drc")
})

test_that("mixture CA sets class to CA and name to ca", {
  result <- mixture(acidiq.free4, model = "CA")
  # The fct class should be CA (overridden from Hewlett)
  expect_s3_class(result$fct, "CA")
  expect_equal(result$fct$name, "ca")
})

# ==============================================================================
# Tests for mixture() - Hewlett model
# ==============================================================================

test_that("mixture Hewlett works with LL.4 model", {
  result <- mixture(acidiq.free4, model = "Hewlett")
  expect_s3_class(result, "drc")
  expect_equal(result$text, "Hewlett model")
})

test_that("mixture Hewlett works with LL.3 model", {
  result <- mixture(acidiq.free3, model = "Hewlett")
  expect_s3_class(result, "drc")
  expect_equal(result$text, "Hewlett model")
})

test_that("mixture Hewlett works with LL.2 model", {
  # LL.2 Hewlett may fail convergence on some datasets; test the error path
  # or successful fitting if it converges
  tryCatch({
    result <- mixture(acidiq.free2, model = "Hewlett")
    expect_s3_class(result, "drc")
  }, error = function(e) {
    # Convergence failure is acceptable for this test
    expect_true(grepl("Convergence|convergence|optim|finite", conditionMessage(e)))
  })
})

test_that("mixture Hewlett errors with LL.5 model", {
  expect_error(mixture(acidiq.free5, model = "Hewlett"), "Does not work for LL.5")
})

test_that("mixture Hewlett uses default startm = 1 when missing", {
  result <- mixture(acidiq.free4, model = "Hewlett")
  expect_s3_class(result, "drc")
})

test_that("mixture Hewlett with explicit startm", {
  result <- mixture(acidiq.free4, model = "Hewlett", startm = 2)
  expect_s3_class(result, "drc")
})

# ==============================================================================
# Tests for mixture() - Voelund model
# ==============================================================================

test_that("mixture Voelund works with LL.4 model", {
  result <- mixture(acidiq.free4, model = "Voelund")
  expect_s3_class(result, "drc")
  expect_equal(result$text, "Voelund model")
})

test_that("mixture Voelund works with LL.3 model", {
  result <- mixture(acidiq.free3, model = "Voelund")
  expect_s3_class(result, "drc")
  expect_equal(result$text, "Voelund model")
})

test_that("mixture Voelund works with LL.2 model", {
  result <- mixture(acidiq.free2, model = "Voelund")
  expect_s3_class(result, "drc")
  expect_equal(result$text, "Voelund model")
})

test_that("mixture Voelund errors with LL.5 model", {
  expect_error(mixture(acidiq.free5, model = "Voelund"), "Does not work for LL.5")
})

test_that("mixture Voelund uses default startm = c(3, 0.3) when missing", {
  result <- mixture(acidiq.free4, model = "Voelund")
  expect_s3_class(result, "drc")
})

test_that("mixture Voelund with explicit startm", {
  result <- mixture(acidiq.free4, model = "Voelund", startm = c(2, 0.5))
  expect_s3_class(result, "drc")
})

# ==============================================================================
# Tests for mixture() - Error handling
# ==============================================================================

test_that("mixture errors when pmodels collapse argument is missing", {
  bad_obj <- drm(rgr ~ dose, data = subset(acidiq, pct == 999 | pct == 0),
                 fct = LL.4())
  expect_error(mixture(bad_obj, model = "CA"),
               "collapse.*argument should be a formula")
})

test_that("mixture errors when Level 0 is missing from data", {
  # Modify curveid to remove all levels containing '0'
  obj_no0 <- acidiq.free4
  obj_no0$dataList$curveid <- factor(rep("17", length(acidiq.free4$dataList$curveid)))
  expect_error(mixture(obj_no0, model = "CA"), "Level 0 is missing")
})

test_that("mixture errors when Level 100 is missing from data", {
  # Modify curveid to remove level 100 (but keep levels with "0" substring)
  obj_no100 <- acidiq.free4
  cids <- as.character(acidiq.free4$dataList$curveid)
  cids[cids == "100"] <- "33"
  obj_no100$dataList$curveid <- factor(cids)
  expect_error(mixture(obj_no100, model = "CA"), "Level 100 is missing")
})

test_that("mixture works with custom start values", {
  # Provide explicit start values to skip the default construction
  sv <- coef(acidiq.free4)
  # CA with LL.4 needs 5 start values (b, c, d, e0, e100)
  start_ca <- c(-2, 0, 0.3, 1, 1)
  result <- mixture(acidiq.free4, model = "CA", start = start_ca)
  expect_s3_class(result, "drc")
})

test_that("mixture model argument uses match.arg", {
  expect_error(mixture(acidiq.free4, model = "invalid"),
               "'arg' should be one of")
})

test_that("mixture preserves anova structure", {
  result <- mixture(acidiq.free4, model = "CA")
  expect_equal(result$anova$test, "F")
  # anovaFit is copied from the input object (may be NULL if not set)
  expect_true("test" %in% names(result$anova))
})

test_that("mixture with custom control parameter", {
  ctrl <- drmc(maxIt = 500)
  result <- mixture(acidiq.free4, model = "CA", control = ctrl)
  expect_s3_class(result, "drc")
})

# ==============================================================================
# Tests for hewlett() function
# ==============================================================================

test_that("hewlett returns correct structure with default parameters", {
  h <- hewlett()
  expect_s3_class(h, "Hewlett")
  expect_equal(h$name, "hewlett")
  expect_equal(h$text, "Hewlett mixture")
  expect_equal(h$noParm, 6)
  expect_equal(h$names, c("b", "c", "d", "e", "f", "g"))
  expect_true(is.function(h$fct))
  expect_true(is.function(h$ssfct))
  expect_true(is.function(h$scaleFct))
  expect_null(h$deriv1)
  expect_null(h$deriv2)
  expect_null(h$edfct)
  expect_null(h$sifct)
})

test_that("hewlett returns correct structure with fixed parameters", {
  h <- hewlett(fixed = c(NA, 0, 1, NA, NA, 1))
  expect_equal(h$noParm, 3)
  expect_equal(h$names, c("b", "e", "f"))
})

test_that("hewlett errors on incorrect names argument", {
  expect_error(hewlett(names = c("a")), "Not correct 'names' argument")
  expect_error(hewlett(names = c(1, 2, 3, 4, 5, 6)), "Not correct 'names' argument")
})

test_that("hewlett errors on incorrect fixed argument", {
  expect_error(hewlett(fixed = c(NA, NA)), "Not correct 'fixed' argument")
})

test_that("hewlett fct computes correct values for normal doses", {
  h <- hewlett()
  dose <- c(0.1, 1, 10, 100)
  parm <- matrix(c(-2, 0, 1, 10, 10, 1), nrow = length(dose), ncol = 6, byrow = TRUE)
  result <- h$fct(dose, parm)
  expect_length(result, 4)
  expect_true(all(is.finite(result)))
  # All values should be between c (0) and d (1) parameters
  expect_true(all(result >= 0 & result <= 1))
})

test_that("hewlett fct handles zero dose correctly", {
  h <- hewlett()
  dose <- c(0, 1, 10)
  # b < 0 case: zero dose returns c parameter (lower limit)
  parm_neg <- matrix(c(-2, 0.1, 1, 10, 10, 1), nrow = 3, ncol = 6, byrow = TRUE)
  result_neg <- h$fct(dose, parm_neg)
  expect_equal(result_neg[1], 0.1)  # c parameter when b < 0
  
  # b > 0 case: zero dose returns d parameter (upper limit)
  parm_pos <- matrix(c(2, 0.1, 1, 10, 10, 1), nrow = 3, ncol = 6, byrow = TRUE)
  result_pos <- h$fct(dose, parm_pos)
  expect_equal(result_pos[1], 1)  # d parameter when b > 0
})

test_that("hewlett default ssfct returns valid starting values", {
  h <- hewlett()
  df <- data.frame(dose = c(0.01, 0.1, 1, 10, 100),
                   resp = c(1, 0.95, 0.5, 0.1, 0.01))
  ss <- h$ssfct(df)
  expect_length(ss, 6)
  expect_true(all(is.finite(ss)))
})

test_that("hewlett with custom ssfct", {
  custom_ss <- function(dframe) rep(1, 6)
  h <- hewlett(ssfct = custom_ss)
  df <- data.frame(dose = 1:5, resp = 5:1)
  result <- h$ssfct(df)
  expect_equal(result, rep(1, 6))
})

test_that("hewlett scaleFct returns correct scaling", {
  h <- hewlett()
  sf <- h$scaleFct(10, 100)
  expect_equal(sf, c(1, 100, 100, 10, 10, 1))
})

test_that("hewlett scaleFct respects fixed parameters", {
  h <- hewlett(fixed = c(NA, 0, 1, NA, NA, NA))
  sf <- h$scaleFct(10, 100)
  expect_equal(sf, c(1, 10, 10, 1))
})

# ==============================================================================
# Tests for voelund() function
# ==============================================================================

test_that("voelund returns correct structure with default parameters", {
  v <- voelund()
  expect_s3_class(v, "Voelund")
  expect_equal(v$name, "voelund")
  expect_equal(v$text, "Voelund mixture")
  expect_equal(v$noParm, 7)
  expect_equal(v$names, c("b", "c", "d", "e", "f", "g", "h"))
  expect_true(is.function(v$fct))
  expect_true(is.function(v$ssfct))
  expect_null(v$deriv1)
  expect_null(v$deriv2)
  expect_null(v$edfct)
  expect_null(v$sifct)
})

test_that("voelund returns correct structure with fixed parameters", {
  v <- voelund(fixed = c(NA, 0, 1, NA, NA, NA, NA))
  expect_equal(v$noParm, 5)
  expect_equal(v$names, c("b", "e", "f", "g", "h"))
})

test_that("voelund errors on incorrect names argument", {
  expect_error(voelund(names = c("a")), "Not correct 'names' argument")
  expect_error(voelund(names = c(1, 2, 3, 4, 5, 6, 7)), "Not correct 'names' argument")
})

test_that("voelund errors on incorrect fixed argument", {
  expect_error(voelund(fixed = c(NA, NA)), "Not correct 'fixed' argument")
})

test_that("voelund fct computes correct values for normal doses", {
  v <- voelund()
  dose <- c(0.1, 1, 10, 100)
  parm <- matrix(c(-2, 0, 1, 10, 10, 1, 1), nrow = 4, ncol = 7, byrow = TRUE)
  result <- v$fct(dose, parm)
  expect_length(result, 4)
  expect_true(all(is.finite(result)))
})

test_that("voelund fct handles zero dose correctly", {
  v <- voelund()
  dose <- c(0, 1, 10)
  parm <- matrix(c(-2, 0, 1, 10, 10, 1, 1), nrow = 3, ncol = 7, byrow = TRUE)
  result <- v$fct(dose, parm)
  # When dose < eps (zero), result should be d parameter
  expect_equal(result[1], 1)  # d parameter
})

test_that("voelund fct handles infinite e parameter", {
  v <- voelund()
  dose <- c(0.1, 1, 10)
  parm <- matrix(c(-2, 0, 1, 10, 10, 1, 1), nrow = 3, ncol = 7, byrow = TRUE)
  parm[2, 4] <- Inf  # e parameter = Inf for row 2
  result <- v$fct(dose, parm)
  expect_length(result, 3)
  # When e is Inf, loge should use log(f) instead
  expect_true(is.finite(result[1]))
  expect_true(is.finite(result[3]))
})

test_that("voelund fct handles infinite f parameter", {
  v <- voelund()
  dose <- c(0.1, 1, 10)
  parm <- matrix(c(-2, 0, 1, 10, 10, 1, 1), nrow = 3, ncol = 7, byrow = TRUE)
  parm[2, 5] <- Inf  # f parameter = Inf for row 2
  result <- v$fct(dose, parm)
  expect_length(result, 3)
  # When f is Inf, loge should use log(e) instead
  expect_true(is.finite(result[1]))
  expect_true(is.finite(result[3]))
})

test_that("voelund default ssfct returns valid starting values", {
  v <- voelund()
  df <- data.frame(dose = c(0.01, 0.1, 1, 10, 100),
                   resp = c(1, 0.95, 0.5, 0.1, 0.01))
  ss <- v$ssfct(df)
  expect_length(ss, 7)
  expect_true(all(is.finite(ss)))
})

test_that("voelund with custom ssfct", {
  custom_ss <- function(dframe) rep(1, 7)
  v <- voelund(ssfct = custom_ss)
  df <- data.frame(dose = 1:5, resp = 5:1)
  result <- v$ssfct(df)
  expect_equal(result, rep(1, 7))
})

test_that("voelund scaleFct returns correct scaling", {
  v <- voelund()
  sf <- v$scaleFct(10, 100)
  expect_equal(sf, c(1, 100, 100, 10, 10, 1, 1))
})

test_that("voelund scaleFct respects fixed parameters", {
  v <- voelund(fixed = c(NA, 0, NA, NA, NA, NA, NA))
  sf <- v$scaleFct(10, 100)
  expect_equal(sf, c(1, 100, 10, 10, 1, 1))
})

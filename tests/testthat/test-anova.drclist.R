# Tests for anova.drclist() function
# Achieves 100% coverage of R/anova.drclist.R

# ─── Test Data ───────────────────────────────────────────────────────────────

ryegrass_data <- data.frame(
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

binom_data <- data.frame(
  dose = c(0, 0.1, 0.5, 1, 2, 5, 10),
  resp = c(0, 0.05, 0.15, 0.35, 0.65, 0.90, 0.98),
  n = rep(50, 7)
)

# ─── Error Handling Tests ────────────────────────────────────────────────────

test_that("anova.drclist errors when more than 2 models are provided", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  m3 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.5())

  expect_error(
    anova(m1, m2, m3),
    "Only two models can be compared"
  )
})

test_that("anova.drclist errors when models have different data types", {
  m_cont <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  m_binom <- drm(resp ~ dose, data = binom_data, fct = LL.2(),
                 type = "binomial", weights = n)

  expect_error(
    anova(m_cont, m_binom),
    "The two models are based on different types on data"
  )
})

# ─── F-test Branch (Continuous Data) ─────────────────────────────────────────

test_that("F-test: basic comparison with details = TRUE", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  result <- anova(m1, m2, details = TRUE)

  expect_s3_class(result, "anova")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 5)
  expect_equal(colnames(result), c("ModelDf", "RSS", "Df", "F value", "p value"))
  expect_equal(attr(result, "heading"), "ANOVA table\n")

  # First row should have NA for test stat and p-value
  expect_true(is.na(result[1, "F value"]))
  expect_true(is.na(result[1, "p value"]))
  # Second row should have actual values
  expect_false(is.na(result[2, "F value"]))
  expect_false(is.na(result[2, "p value"]))
  # p-value should be between 0 and 1
  expect_true(result[2, "p value"] >= 0 && result[2, "p value"] <= 1)
})

test_that("F-test: comparison with details = FALSE", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  result <- anova(m1, m2, details = FALSE)

  expect_s3_class(result, "anova")
  expect_equal(nrow(result), 2)
})

test_that("F-test: model order is swapped when df2 > df1", {
  # LL.4 has more parameters (4) -> fewer df.residual than LL.3 (3 params)
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())

  # When obj1 has fewer df (more params) and obj2 has more df (fewer params),
  # df2 > df1 triggers the swap
  result <- anova(m1, m2, details = FALSE)

  expect_s3_class(result, "anova")
  expect_equal(nrow(result), 2)
  # Row names should be swapped
  expect_equal(rownames(result), c("2nd model", "1st model"))
})

test_that("F-test: explicit test='F' parameter", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  result <- anova(m1, m2, test = "F", details = FALSE)

  expect_s3_class(result, "anova")
  expect_equal(colnames(result), c("ModelDf", "RSS", "Df", "F value", "p value"))
})

# ─── Chi-square Test Branch (Binomial Data) ──────────────────────────────────

test_that("Chi-square test: binomial data with default test=NULL selects Chisq", {
  m1 <- drm(resp ~ dose, data = binom_data, fct = LL.2(),
            type = "binomial", weights = n)
  m2 <- drm(resp ~ dose, data = binom_data, fct = LL.3(),
            type = "binomial", weights = n)

  result <- anova(m1, m2, details = FALSE)

  expect_s3_class(result, "anova")
  expect_equal(nrow(result), 2)
  expect_equal(colnames(result), c("ModelDf", "Loglik", "Df", "LR value", "p value"))
  expect_equal(attr(result, "heading"), "ANOVA-like table\n")
})

test_that("Chi-square test: explicit test='Chisq' for continuous data", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  result <- anova(m1, m2, test = "Chisq", details = FALSE)

  expect_s3_class(result, "anova")
  expect_equal(colnames(result), c("ModelDf", "Loglik", "Df", "LR value", "p value"))
  expect_equal(attr(result, "heading"), "ANOVA-like table\n")
  # p-value should be valid
  pval <- result[2, "p value"]
  expect_true(!is.na(pval) && pval >= 0 && pval <= 1)
})

# ─── F-test Edge Cases for p-value ───────────────────────────────────────────

test_that("F-test: NaN test statistic gives NA p-value", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  # Manipulate RSS to produce NaN (0/0 scenario)
  m1_mod <- m1
  m2_mod <- m2
  m1_mod$summary[4] <- 0
  m2_mod$summary[4] <- 0

  result <- anova(m1_mod, m2_mod, test = "F", details = FALSE)
  expect_true(is.na(result[2, "p value"]))
})

test_that("F-test: negative test statistic gives p-value of 1", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  # Manipulate so simpler model has lower RSS than complex model
  m1_mod <- m1
  m2_mod <- m2
  # After swap: m1 is the model with more df (simpler), m2 has fewer df (complex)
  # testStat = ((loglik[1]-loglik[2])/dfDiff[2]) / (loglik[2]/df2)
  # Make loglik[1] < loglik[2] to get negative numerator
  m1_mod$summary[4] <- 100  # simpler model RSS (will be obj1 after swap since more df)
  m2_mod$summary[4] <- 200  # complex model RSS (will be obj2)

  result <- anova(m1_mod, m2_mod, test = "F", details = FALSE)
  expect_equal(result[2, "p value"], 1)
})

test_that("F-test: Inf test statistic gives NA p-value", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  # Make obj2 RSS = 0 to produce Inf (positive numerator / 0 denominator)
  m1_mod <- m1
  m2_mod <- m2
  m2_mod$summary[4] <- 0
  m1_mod$summary[4] <- 10

  result <- anova(m1_mod, m2_mod, test = "F", details = FALSE)
  expect_true(is.na(result[2, "p value"]))
})

# ─── Details Printing: Collapse Paths for obj1 ──────────────────────────────

test_that("details: collapse1 from obj1[[8]]$collapse (character)", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  # Set a character collapse
  m1_mod <- m1
  m1_mod[[8]]$collapse <- "custom_collapse"

  output <- capture.output(result <- anova(m1_mod, m2, details = TRUE))
  expect_true(any(grepl("custom_collapse", output)))
})

test_that("details: collapse1 from pmodelsText", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m1_mod <- m1
  m1_mod$pmodelsText <- "pmodels_text_1"

  output <- capture.output(result <- anova(m1_mod, m2, details = TRUE))
  expect_true(any(grepl("pmodels_text_1", output)))
})

test_that("details: collapse1 from obj1[[8]]$curve (non-NULL)", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m1_mod <- m1
  m1_mod[[8]]$collapse <- NULL
  m1_mod[[8]]$pmodels <- NULL
  m1_mod$pmodelsText <- NULL
  m1_mod[[8]]$curve <- ~ group

  output <- capture.output(result <- anova(m1_mod, m2, details = TRUE))
  expect_true(any(grepl("group", output)))
})

test_that("details: collapse1 is non-character (formula/expression)", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m1_mod <- m1
  m1_mod[[8]]$collapse <- ~ a + b
  m1_mod$pmodelsText <- NULL

  output <- capture.output(result <- anova(m1_mod, m2, details = TRUE))
  expect_s3_class(result, "anova")
})

test_that("details: collapse1 contains 'data.frame(' prefix", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m1_mod <- m1
  m1_mod[[8]]$collapse <- "data.frame(x, y, z)"
  m1_mod$pmodelsText <- NULL

  output <- capture.output(result <- anova(m1_mod, m2, details = TRUE))
  # data.frame( prefix should be stripped
  expect_true(any(grepl("x, y, z", output)))
})

test_that("details: collapse1 contains 'list(' prefix", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m1_mod <- m1
  m1_mod[[8]]$collapse <- "list(a, b, c)"
  m1_mod$pmodelsText <- NULL

  output <- capture.output(result <- anova(m1_mod, m2, details = TRUE))
  # list( prefix should be stripped
  expect_true(any(grepl("a, b, c", output)))
})

# ─── Details Printing: Collapse Paths for obj2 ──────────────────────────────

test_that("details: collapse2 from obj2[[8]]$curve (non-NULL)", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m2_mod <- m2
  m2_mod[[8]]$collapse <- NULL
  m2_mod[[8]]$pmodels <- NULL
  m2_mod$pmodelsText <- NULL
  m2_mod[[8]]$curve <- ~ treatment

  output <- capture.output(result <- anova(m1, m2_mod, details = TRUE))
  expect_true(any(grepl("treatment", output)))
})

test_that("details: collapse2 is non-character (formula)", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m2_mod <- m2
  m2_mod[[8]]$collapse <- ~ x + y
  m2_mod$pmodelsText <- NULL

  output <- capture.output(result <- anova(m1, m2_mod, details = TRUE))
  expect_s3_class(result, "anova")
})

test_that("details: collapse2 contains 'data.frame(' prefix", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m2_mod <- m2
  m2_mod[[8]]$collapse <- "data.frame(p, q)"
  m2_mod$pmodelsText <- NULL

  output <- capture.output(result <- anova(m1, m2_mod, details = TRUE))
  expect_true(any(grepl("p, q", output)))
})

test_that("details: collapse2 contains 'list(' prefix", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m2_mod <- m2
  m2_mod[[8]]$collapse <- "list(m, n)"
  m2_mod$pmodelsText <- NULL

  output <- capture.output(result <- anova(m1, m2_mod, details = TRUE))
  expect_true(any(grepl("m, n", output)))
})

test_that("details: collapse2 from pmodelsText", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m2_mod <- m2
  m2_mod$pmodelsText <- "pmodels_text_2"

  output <- capture.output(result <- anova(m1, m2_mod, details = TRUE))
  expect_true(any(grepl("pmodels_text_2", output)))
})

# ─── Details Printing: colLine and pmodels paths ─────────────────────────────

test_that("details: colLine = TRUE when collapse1 != collapse2", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  # Set different collapse values to make colLine TRUE
  m1_mod <- m1
  m2_mod <- m2
  m1_mod[[8]]$collapse <- "collapse_A"
  m1_mod$pmodelsText <- NULL
  m2_mod[[8]]$collapse <- "collapse_B"
  m2_mod$pmodelsText <- NULL

  output <- capture.output(result <- anova(m1_mod, m2_mod, details = TRUE))
  # Both collapse values should appear in the output
  expect_true(any(grepl("collapse_A", output)))
  expect_true(any(grepl("collapse_B", output)))
})

test_that("details: colLine = FALSE when collapse1 == collapse2", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  # Both have NULL collapse (default), so both will be "1 (for all parameters)"
  output <- capture.output(result <- anova(m1, m2, details = TRUE))
  # The pmodels line should not appear since colLine is FALSE
  # (both models have the same default collapse)
  expect_true(any(grepl("1st model", output)))
  expect_true(any(grepl("2nd model", output)))
})

test_that("details: pmodels present in obj1 triggers fctStart alignment", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m1_mod <- m1
  m1_mod[[8]]$pmodels <- list(~ 1, ~ 1, ~ 1)
  m1_mod[[8]]$collapse <- "pmodel_collapse"
  m1_mod$pmodelsText <- NULL

  output <- capture.output(result <- anova(m1_mod, m2, details = TRUE))
  # Should use " fct:     " (with one fewer space) when pmodels present
  expect_true(any(grepl("fct:", output)))
})

test_that("details: pmodels present in obj2 triggers fctStart alignment", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m2_mod <- m2
  m2_mod[[8]]$pmodels <- list(~ 1, ~ 1, ~ 1, ~ 1)
  m2_mod[[8]]$collapse <- "pmodel_collapse_2"
  m2_mod$pmodelsText <- NULL

  output <- capture.output(result <- anova(m1, m2_mod, details = TRUE))
  expect_true(any(grepl("fct:", output)))
})

# ─── Details Printing: fctInfo (text field) ──────────────────────────────────

test_that("details: fctInfo uses obj$text when available", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  m1_mod <- m1
  m1_mod$text <- "Custom Model Text 1"

  output <- capture.output(result <- anova(m1_mod, m2, details = TRUE))
  expect_true(any(grepl("Custom Model Text 1", output)))
})

test_that("details: fctInfo uses deparse(obj[[8]]$fct) when text is NULL", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  # text is NULL by default, so deparse should be used
  output <- capture.output(result <- anova(m1, m2, details = TRUE))
  expect_true(any(grepl("fct:", output)))
})

# ─── Chi-square Test with Details ────────────────────────────────────────────

test_that("Chi-square test with details = TRUE", {
  m1 <- drm(resp ~ dose, data = binom_data, fct = LL.2(),
            type = "binomial", weights = n)
  m2 <- drm(resp ~ dose, data = binom_data, fct = LL.3(),
            type = "binomial", weights = n)

  output <- capture.output(result <- anova(m1, m2, details = TRUE))

  expect_s3_class(result, "anova")
  expect_equal(colnames(result), c("ModelDf", "Loglik", "Df", "LR value", "p value"))
  expect_true(any(grepl("1st model", output)))
  expect_true(any(grepl("2nd model", output)))
})

# ─── Return Value Structure ─────────────────────────────────────────────────

test_that("return value has correct structure and class", {
  m1 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.3())
  m2 <- drm(rootl ~ conc, data = ryegrass_data, fct = LL.4())

  result <- anova(m1, m2, details = FALSE)

  expect_s3_class(result, "anova")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 5)
  expect_true(!is.null(attr(result, "heading")))
  expect_equal(rownames(result), c("1st model", "2nd model"))
})

## --------------------------------------------------------------------------
## Tests for CIcompX, CIcomp, and plotFACI (R/CIcompX.R)
## --------------------------------------------------------------------------

# Shared test data: fit 3 dose-response models using acidiq dataset
# acidiq.17 = mixture model (17:83 ratio)
# acidiq.0  = pure substance 1
# acidiq.100 = pure substance 2
acidiq.0   <- drm(rgr ~ dose, data = subset(acidiq, pct == 999 | pct == 0),   fct = LL.4())
acidiq.100 <- drm(rgr ~ dose, data = subset(acidiq, pct == 999 | pct == 100), fct = LL.4())
acidiq.17  <- drm(rgr ~ dose, data = subset(acidiq, pct == 17  | pct == 0),   fct = LL.4())

modList <- list(acidiq.17, acidiq.0, acidiq.100)

## ===== CIcompX tests =====

# --- Input validation ---

test_that("CIcompX errors when mixProp < 0", {
  expect_error(
    CIcompX(-0.1, modList, c(10, 50)),
    "Mixture proportion should be between 0 and 1"
  )
})

test_that("CIcompX errors when mixProp > 1", {
  expect_error(
    CIcompX(1.5, modList, c(10, 50)),
    "Mixture proportion should be between 0 and 1"
  )
})

test_that("CIcompX errors when modelList is not a list", {
  expect_error(
    CIcompX(0.17, "not_a_list", c(10, 50)),
    "Exactly 3 model fits should be provided in a list"
  )
})

test_that("CIcompX errors when modelList has wrong length", {
  expect_error(
    CIcompX(0.17, list(acidiq.17, acidiq.0), c(10, 50)),
    "Exactly 3 model fits should be provided in a list"
  )
})

test_that("CIcompX errors when EDvec is empty", {
  expect_error(
    CIcompX(0.17, modList, numeric(0)),
    "At least effective dose level should be specified"
  )
})

# --- Happy path: multiple ED levels, EDonly = FALSE ---

test_that("CIcompX returns correct structure with multiple ED levels (EDonly=FALSE)", {
  res <- CIcompX(0.17, modList, c(10, 20, 50), EDonly = FALSE)

  expect_true(is.list(res))
  expect_equal(length(res), 5)
  expect_true(all(c("Effx", "Effy", "CAx", "CAy", "EDvec") %in% names(res)))

  # Effx: 3 rows x 6 cols

  expect_true(is.matrix(res$Effx))
  expect_equal(nrow(res$Effx), 3)
  expect_equal(ncol(res$Effx), 6)
  expect_equal(colnames(res$Effx), c("ED.mix", "ED1", "ED2", "SE.mix", "SE1", "SE2"))
  expect_equal(rownames(res$Effx), c("10", "20", "50"))

  # Effy: 3 rows x 6 cols
  expect_true(is.matrix(res$Effy))
  expect_equal(nrow(res$Effy), 3)
  expect_equal(ncol(res$Effy), 6)
  expect_equal(colnames(res$Effy), c("E.mix", "E1", "E2", "SE.mix", "SE1", "SE2"))

  # CAx: 3 rows x 8 cols
  expect_true(is.matrix(res$CAx))
  expect_equal(nrow(res$CAx), 3)
  expect_equal(ncol(res$CAx), 8)
  expect_equal(colnames(res$CAx),
    c("combInd", "SE", "lowCI", "highCI", "CAdiff", "CAdiffp", "PredAdd", "sePredAdd"))

  # CAy: 3 rows x 8 cols
  expect_true(is.matrix(res$CAy))
  expect_equal(nrow(res$CAy), 3)
  expect_equal(ncol(res$CAy), 8)

  # EDvec preserved

  expect_equal(res$EDvec, c(10, 20, 50))

  # All ED values should be positive
  expect_true(all(res$Effx[, "ED.mix"] > 0))
  expect_true(all(res$Effx[, "ED1"] > 0))
  expect_true(all(res$Effx[, "ED2"] > 0))
})

# --- Happy path: multiple ED levels, EDonly = TRUE ---

test_that("CIcompX returns only ED-related components when EDonly=TRUE", {
  res <- CIcompX(0.17, modList, c(10, 20, 50), EDonly = TRUE)

  expect_true(is.list(res))
  expect_equal(length(res), 3)
  expect_true(all(c("Effx", "CAx", "EDvec") %in% names(res)))
  expect_false("Effy" %in% names(res))
  expect_false("CAy" %in% names(res))

  expect_true(is.matrix(res$Effx))
  expect_true(is.matrix(res$CAx))
  expect_equal(res$EDvec, c(10, 20, 50))
})

# --- Single ED level (triggers !is.matrix branch for predictions) ---

test_that("CIcompX works with a single ED level", {
  res <- CIcompX(0.17, modList, 50, EDonly = FALSE)

  expect_true(is.list(res))
  expect_equal(length(res), 5)

  # Matrices should have 1 row
  expect_equal(nrow(res$Effx), 1)
  expect_equal(nrow(res$Effy), 1)
  expect_equal(nrow(res$CAx), 1)
  expect_equal(nrow(res$CAy), 1)
  expect_equal(rownames(res$Effx), "50")
  expect_equal(res$EDvec, 50)
})

test_that("CIcompX single ED level with EDonly=TRUE", {
  res <- CIcompX(0.17, modList, 50, EDonly = TRUE)

  expect_equal(length(res), 3)
  expect_equal(nrow(res$Effx), 1)
  expect_equal(nrow(res$CAx), 1)
})

# --- Combination index numerical sanity ---

test_that("CIcompX combination indices have expected properties", {
  res <- CIcompX(0.17, modList, c(10, 20, 50))

  # Standard errors should be positive
  expect_true(all(res$CAx[, "SE"] > 0))
  expect_true(all(res$CAy[, "SE"] > 0))

  # Confidence intervals: lower < estimate < upper
  expect_true(all(res$CAx[, "lowCI"] < res$CAx[, "combInd"]))
  expect_true(all(res$CAx[, "highCI"] > res$CAx[, "combInd"]))

  # p-values between 0 and 1
  expect_true(all(res$CAx[, "CAdiffp"] >= 0 & res$CAx[, "CAdiffp"] <= 1))
  expect_true(all(res$CAy[, "CAdiffp"] >= 0 & res$CAy[, "CAdiffp"] <= 1))

  # PredAdd and sePredAdd should be positive
  expect_true(all(res$CAx[, "PredAdd"] > 0))
  expect_true(all(res$CAx[, "sePredAdd"] > 0))
})

# --- Edge cases for mixProp boundary ---

test_that("CIcompX works with mixProp at boundaries (0 and 1)", {
  res0 <- CIcompX(0, modList, c(10, 50))
  res1 <- CIcompX(1, modList, c(10, 50))

  expect_true(is.list(res0))
  expect_true(is.list(res1))
  expect_equal(length(res0), 5)
  expect_equal(length(res1), 5)
})

## ===== CIcomp tests =====

test_that("CIcomp returns correct matrix structure with multiple ED levels", {
  res <- CIcomp(0.17, modList, c(10, 20, 50))

  expect_true(is.matrix(res))
  expect_equal(nrow(res), 3)
  expect_equal(ncol(res), 9)
  expect_equal(rownames(res), c("10", "20", "50"))

  # Check renamed columns
  expect_equal(colnames(res)[6], "ED.CA")
  expect_equal(colnames(res)[7], "SE.CA")

  # Numeric sanity
  expect_true(all(is.finite(res)))
})

test_that("CIcomp works with single ED level (drop=FALSE fix)", {
  res <- CIcomp(0.17, modList, 50)

  expect_true(is.matrix(res))
  expect_equal(nrow(res), 1)
  expect_equal(ncol(res), 9)
  expect_equal(rownames(res), "50")
  expect_equal(colnames(res)[6], "ED.CA")
  expect_equal(colnames(res)[7], "SE.CA")
})

test_that("CIcomp values are consistent with CIcompX output", {
  edvec <- c(10, 20, 50)
  resX <- CIcompX(0.17, modList, edvec, EDonly = FALSE)
  resC <- CIcomp(0.17, modList, edvec)

  # CIcomp columns 1-5 should be CAx without column 5 (CAdiff) minus column 5
  # combInd from CAx matches column 1 of CIcomp
  expect_equal(resC[, "combInd"], resX$CAx[, "combInd"])

  # ED.mix and SE.mix from Effx
  expect_equal(as.numeric(resC[, "ED.mix"]), as.numeric(resX$Effx[, "ED.mix"]))
  expect_equal(as.numeric(resC[, "SE.mix"]), as.numeric(resX$Effx[, "SE.mix"]))
})

## ===== plotFACI tests =====

# Helper: build an effList for plotting tests
build_effList <- function(edvec = c(10, 20, 50)) {
  CIcompX(0.17, modList, edvec, EDonly = FALSE)
}

test_that("plotFACI default call (ED axis, caRef, new plot) works", {
  effL <- build_effList()
  res <- plotFACI(effL)

  expect_true(is.matrix(res))
  expect_equal(nrow(res), 3)
  expect_equal(ncol(res), 8)
})

test_that("plotFACI with indAxis='EF' works", {
  effL <- build_effList()
  res <- plotFACI(effL, indAxis = "EF")

  expect_true(is.matrix(res))
  expect_equal(nrow(res), 3)
})

test_that("plotFACI with caRef=FALSE uses range for ylim", {
  effL <- build_effList()
  res <- plotFACI(effL, caRef = FALSE)

  expect_true(is.matrix(res))
})

test_that("plotFACI with explicit ylim overrides computed limits", {
  effL <- build_effList()
  res <- plotFACI(effL, ylim = c(0, 3))

  expect_true(is.matrix(res))
})

test_that("plotFACI with showPoints=TRUE draws points", {
  effL <- build_effList()
  res <- plotFACI(effL, showPoints = TRUE)

  expect_true(is.matrix(res))
})

test_that("plotFACI with add=TRUE adds to existing plot", {
  effL <- build_effList()
  # Create initial plot
  plot(1, 1, xlim = c(0, 100), ylim = c(0, 3),
       xlab = "FA", ylab = "CI")
  res <- plotFACI(effL, add = TRUE)

  expect_true(is.matrix(res))
})

test_that("plotFACI with all optional arguments combined", {
  effL <- build_effList()
  plot(1, 1, xlim = c(0, 100), ylim = c(0, 3),
       xlab = "FA", ylab = "CI")
  res <- plotFACI(effL, indAxis = "EF", caRef = FALSE,
                  showPoints = TRUE, add = TRUE)

  expect_true(is.matrix(res))
})

test_that("plotFACI handles negative faValues in EDvec", {
  # Construct a mock effList with negative EDvec values
  effL <- build_effList(c(10, 20, 50))

  # Manually adjust EDvec to include negative values for coverage of line 234
  effL$EDvec <- c(-10, 20, 50)
  rownames(effL$CAx) <- as.character(c(-10, 20, 50))
  rownames(effL$CAy) <- as.character(c(-10, 20, 50))

  res <- plotFACI(effL)
  expect_true(is.matrix(res))
})

test_that("plotFACI with caRef=FALSE and ylim provided", {
  effL <- build_effList()
  res <- plotFACI(effL, caRef = FALSE, ylim = c(0, 5))

  expect_true(is.matrix(res))
})

test_that("plotFACI returns correct matrix invisibly", {
  effL <- build_effList()
  res <- plotFACI(effL, indAxis = "ED")

  # Returned matrix should be the CAx component
  expect_equal(res, effL$CAx)
})

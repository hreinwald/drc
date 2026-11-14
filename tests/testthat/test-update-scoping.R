# update.drc() must resolve the stored call in the environment the model was fitted in
# (not in the caller's or the global environment), so that refits by update(), mselect()
# and maED() also work for fits made inside functions or sourced into a new.env().

fit_in_function <- function(fct = LL.4()) {
  drm.formula <- rootl ~ conc
  df <- ryegrass
  drm(drm.formula, data = df, fct = fct)
}

test_that("update() works outside the environment the model was fitted in", {
  m <- fit_in_function()
  m2 <- update(m, fct = W1.4())
  expect_s3_class(m2, "drc")
  expect_equal(coef(m2), coef(drm(rootl ~ conc, data = ryegrass, fct = W1.4())))
})

test_that("objects of the same name in the global environment are not used", {
  assign("drm.formula", conc ~ rootl, envir = globalenv())
  assign("df", ryegrass[1:12, ], envir = globalenv())
  on.exit(rm(list = c("drm.formula", "df"), envir = globalenv()))

  m2 <- update(fit_in_function(), fct = W1.4())
  expect_equal(nrow(m2$data), nrow(ryegrass))
  expect_equal(coef(m2), coef(drm(rootl ~ conc, data = ryegrass, fct = W1.4())))
})

test_that("mselect() and maED() refit all candidates of a fit made in a function", {
  m <- fit_in_function()
  ref <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

  refit_warnings <- character()
  collect <- function(expr) withCallingHandlers(expr, warning = function(w) {
    refit_warnings <<- c(refit_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  ms <- collect(mselect(m, list(W1.4(), LL.3())))
  expect_false(anyNA(ms[, "IC"]))
  expect_equal(ms, suppressWarnings(mselect(ref, list(W1.4(), LL.3()))))

  ma <- collect(maED(m, list(W1.4()), 50, display = FALSE))
  expect_equal(ma, suppressWarnings(maED(ref, list(W1.4()), 50, display = FALSE)))
  expect_false(any(grepl("Could not fit|excluded|not found", refit_warnings)))
})

test_that("a script sourced into a fresh environment can be refitted", {
  env <- new.env()
  evalq({
    drm.formula <- rootl ~ conc
    df <- ryegrass
    m <- drm(drm.formula, data = df, fct = LL.4())
  }, env)
  expect_false(anyNA(suppressWarnings(mselect(env$m, list(W1.4())))[, "IC"]))
})

test_that("robust = 'tukey' works with local formula and data symbols", {
  f <- function() {
    fo <- rootl ~ conc
    d <- ryegrass
    drm(fo, data = d, fct = LL.4(), robust = "tukey")
  }
  m <- f()
  expect_s3_class(m, "drc")
  expect_equal(coef(m), coef(drm(rootl ~ conc, data = ryegrass, fct = LL.4(), robust = "tukey")))
})

test_that("new arguments are evaluated in the caller's environment, also when chained", {
  m <- fit_in_function()
  refit <- function(object) {
    myfct <- W1.3()
    update(object, fct = myfct)
  }
  m2 <- refit(m)
  expect_equal(coef(m2), coef(drm(rootl ~ conc, data = ryegrass, fct = W1.3())))

  m3 <- update(m2, fct = LL.3())
  expect_equal(coef(m3), coef(drm(rootl ~ conc, data = ryegrass, fct = LL.3())))

  newdata <- function(object) {
    df <- ryegrass[1:18, ]
    update(object, data = df)
  }
  expect_no_warning(m4 <- newdata(m))
  expect_equal(nrow(m4$data), 18)
})

test_that("data changed after the fit gives a warning and the stored data is used", {
  env <- new.env()
  env$d <- ryegrass
  m <- evalq(drm(rootl ~ conc, data = d, fct = LL.4()), env)
  env$d <- ryegrass[1:12, ]

  expect_warning(m2 <- update(m, fct = W1.4()), "differs from the data")
  expect_equal(nrow(m2$data), nrow(ryegrass))
})

test_that("drm() restores the 'na.action' option, also for embedded function objects", {
  old <- options(na.action = "na.omit")
  on.exit(options(old), add = TRUE)

  drm(rootl ~ conc, data = ryegrass, fct = LL.4(), na.action = na.exclude)
  expect_identical(getOption("na.action"), "na.omit")

  m <- do.call(drm, list(formula = rootl ~ conc, data = ryegrass, fct = LL.4(), na.action = na.omit))
  expect_identical(getOption("na.action"), "na.omit")
  expect_error(drm(rootl ~ conc, data = ryegrass, fct = LL.4(), na.action = na.exclude, type = "foo"))
  expect_identical(getOption("na.action"), "na.omit")

  # mselect()/modelFit() rely on the option being intact
  ms <- expect_no_warning(mselect(m, list(LL.3())), message = "lack-of-fit")
  expect_false(anyNA(ms[, "Lack of fit"]))
})

test_that("fits without the stored environment are evaluated in the caller's environment", {
  m <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  attr(m, "callEnv") <- NULL
  expect_equal(coef(update(m, fct = W1.4())), coef(drm(rootl ~ conc, data = ryegrass, fct = W1.4())))
  expect_identical(update(m, fct = W1.4(), evaluate = FALSE)$fct, quote(W1.4()))
})

# Changelog

## drc (development version)

### Bug Fixes

- Fixed
  [`logistic()`](https://hreinwald.github.io/drc/reference/logistic.md)
  / [`L.3()`](https://hreinwald.github.io/drc/reference/L.3.md) /
  [`L.4()`](https://hreinwald.github.io/drc/reference/L.4.md) /
  [`L.5()`](https://hreinwald.github.io/drc/reference/L.5.md) ED:
  relative levels were mirrored (ED10 and ED90 swapped) for decreasing
  curves (b \> 0) and `type = "absolute"` returned the dose of the
  mirrored response. The edfct now uses
  [`EDhelper()`](https://hreinwald.github.io/drc/reference/EDhelper.md)
  like the log-logistic models.
- Fixed
  [`gompertz()`](https://hreinwald.github.io/drc/reference/gompertz.md)
  ED gradient, which was multiplied by the ED and gave wrong standard
  errors and confidence intervals for all `G.*` models.
- Fixed
  [`llogistic2()`](https://hreinwald.github.io/drc/reference/llogistic2.md)
  ED gradient with respect to `f`, which gave wrong standard errors for
  [`LL2.5()`](https://hreinwald.github.io/drc/reference/LL2.5.md).
- Fixed
  [`weibull2x()`](https://hreinwald.github.io/drc/reference/weibull2x.md)
  ED: the gradient omitted the lag time `t0`, so
  [`ED()`](https://hreinwald.github.io/drc/reference/ED.md) failed for
  [`W2x.3()`](https://hreinwald.github.io/drc/reference/W2x.3.md) /
  [`W2x.4()`](https://hreinwald.github.io/drc/reference/W2x.4.md);
  `type = "absolute"` now delegates to
  [`weibull2()`](https://hreinwald.github.io/drc/reference/weibull2.md)
  instead of returning the mirrored response.
- Fixed
  [`cedergreen()`](https://hreinwald.github.io/drc/reference/cedergreen.md)
  first derivative with respect to `e` (missing factor), which affected
  [`predict()`](https://rdrr.io/r/stats/predict.html) standard errors,
  sandwich estimators and studentised residuals for all `CRS.*` models.
- Fixed [`vcov()`](https://rdrr.io/r/stats/vcov.html) /
  [`summary()`](https://rdrr.io/r/base/summary.html) /
  [`confint()`](https://rdrr.io/r/stats/confint.html) /
  [`ED()`](https://hreinwald.github.io/drc/reference/ED.md) /
  [`predict()`](https://rdrr.io/r/stats/predict.html) for
  `drm(..., separate = TRUE)` fits, which failed with
  `could not find function "adiag"` unless `magic` was attached; the
  block-diagonal assembly is now internal and `magic` was dropped from
  `Suggests`.
- `ED(interval = "fls")` now stops with an explanatory error (reason and
  how to fix) when the model returns ED on the dose scale; `"fls"`
  back-transforms with [`exp()`](https://rdrr.io/r/base/Log.html) and is
  only valid for models with a log(ED) parameter (`LL2.*`,
  `LN.*(loge = TRUE)`, flagged by the new `edLogScale` element of the
  model list) or log-dose fits (`logBase`).
  [`ED_robust()`](https://hreinwald.github.io/drc/reference/ED_robust.md)
  raises the same error instead of returning `NA` rows.
  [`get_ed_interval()`](https://hreinwald.github.io/drc/reference/get_ed_interval.md)
  no longer recommends `"fls"` for dose-scale models (`small_n = FALSE`
  now gives `"delta"`) and returns `"fls"` for log(ED) models.
- [`maED()`](https://hreinwald.github.io/drc/reference/maED.md) and
  [`maED_robust()`](https://hreinwald.github.io/drc/reference/maED_robust.md)
  now stop with an explanatory error if the base model or a candidate
  has a log(ED) parameterisation (`LL2.*`, `LN.*(loge = TRUE)`), since
  averaging mixes log-scale and dose-scale EDs; use the dose-scale
  equivalents `LL.*` / `LN.*`.
- Fixed
  [`ED_robust()`](https://hreinwald.github.io/drc/reference/ED_robust.md)
  for `LL2.*` models, which reported log(ED) as ED and turned negative
  log-EDs into `NA`; estimates and limits are now back-transformed and
  the SE is the delta-method SE of the dose.
- Fixed
  [`gammadr()`](https://hreinwald.github.io/drc/reference/gammadr.md)
  first derivative with respect to the shape parameter (missing digamma
  term).
- [`predict()`](https://rdrr.io/r/stats/predict.html) now computes
  standard errors and intervals with a numerical gradient for models
  without analytic `deriv1` (`NEC.*`, `CRS.6`, `baro5`, …) instead of
  silently returning plain predictions.
- [`confint()`](https://rdrr.io/r/stats/confint.html) and
  `ED(interval = "delta")` now work for `type = "ssd"` fits.
- `drm(separate = TRUE)` now forwards `robust` (except `"trimmed"`,
  `"winsor"`, `"tukey"`), `bcVal`, `bcAdd`, `lowerl`, `upperl`,
  `logDose` and `start` to the per-curve fits instead of silently
  ignoring them; `varcov`, `pshifts` and the three M-estimators raise an
  error.
- Fixed [`maED()`](https://hreinwald.github.io/drc/reference/maED.md)
  Buckland/Kang standard errors, which centred on the unweighted mean of
  the model estimates instead of the model-averaged estimate.
- Fixed `comped(log = TRUE, operator = "/")`, which printed the log
  ratio with a confidence interval on the ratio scale; estimate,
  delta-method SE and interval are now all on the ratio scale.
- Fixed [`simDR()`](https://hreinwald.github.io/drc/reference/simDR.md)
  failing for increasing curves (missing `reference` argument).
- [`noEffect()`](https://hreinwald.github.io/drc/reference/noEffect.md)
  now stops with a clear message for response types without a null model
  and no longer rounds tiny p-values to 0.
- [`drm()`](https://hreinwald.github.io/drc/reference/drm.md) evaluates
  a `pmodels` argument in the calling environment (not only
  `.GlobalEnv`) when `data` is `NULL`.

### Changes

- Added `tests/testthat/test-numerical-derivatives-ed.R` and
  `tests/testthat/test-review-fixes-f07-f17.R`: cross-check analytic
  `deriv1`/`derivx`, ED values and ED standard errors of all model
  functions against finite differences, and cover the fixes listed
  above.
- Documentation: every model-function help page now has a runnable
  [`drm()`](https://hreinwald.github.io/drc/reference/drm.md) +
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) example, and
  most families show the effect of fixing the lower/upper asymptotes.
  The `arandaordaz` documentation was corrected; roxygen descriptions of
  [`drm()`](https://hreinwald.github.io/drc/reference/drm.md)
  (`separate`, `pmodels`),
  [`ED_robust()`](https://hreinwald.github.io/drc/reference/ED_robust.md),
  [`maED_robust()`](https://hreinwald.github.io/drc/reference/maED_robust.md),
  [`predict.drc()`](https://hreinwald.github.io/drc/reference/predict.drc.md),
  [`noEffect()`](https://hreinwald.github.io/drc/reference/noEffect.md)
  and [`comped()`](https://hreinwald.github.io/drc/reference/comped.md)
  were updated.

------------------------------------------------------------------------

## drc 3.3.2

### New Features

- Enhanced
  [`plot.drc()`](https://hreinwald.github.io/drc/reference/plot.drc.md):
  added `errbar.lwd` parameter for independent control of error bar line
  width in `type = "bars"` plots. When `NULL` (default), error bars
  inherit the line width from the `lwd` argument or fall back to
  `par("lwd")`.

### Bug Fixes

- Fixed
  [`update.drc()`](https://hreinwald.github.io/drc/reference/update.drc.md)
  to fall back to stored data (`object$origData`) when `call$data`
  cannot be resolved in the calling frame, enabling use of
  [`update()`](https://rdrr.io/r/stats/update.html) inside
  [`lapply()`](https://rdrr.io/r/base/lapply.html),
  [`purrr::map()`](https://purrr.tidyverse.org/reference/map.html), and
  other functional programming contexts.
- Fixed `vcDisc()` to validate that the inverse Hessian has non-negative
  variances before returning, preventing invalid variance-covariance
  matrices from propagating downstream.

### Changes

- Tightened default `relTol` in
  [`drmc()`](https://hreinwald.github.io/drc/reference/drmc.md) from
  `1e-7` to `1e-10` for improved cross-platform reproducibility of
  optimization results.
- Added comparative analysis vignette (`comparative-analysis.Rmd`)
  documenting differences between `hreinwald/drc` and the original
  `DoseResponse/drc` package.
- Added “Articles” section to pkgdown navigation in `_pkgdown.yml` with
  categorized guides and technical reports.
- Fixed flaky `maED` tests that depended on platform-specific
  convergence behavior.
- Added test suites for `boxcox.drc` (functional programming context),
  `update.drc` (data resolution fallback), and `vcDisc` (singular
  Hessian handling).

------------------------------------------------------------------------

## drc 3.3.1

### New Features

- Enhanced
  [`plot.drc()`](https://hreinwald.github.io/drc/reference/plot.drc.md):
  error bars in `type = "bars"` plots now match curve colors by default.
  Added `errbar.col` parameter to allow manual control of error bar
  colors. Set `errbar.col = "black"` to restore the previous behavior of
  black error bars.

### Bug Fixes

- Fixed [`summary()`](https://rdrr.io/r/base/summary.html) crash for
  binomial `drm` models when the Hessian is singular
  (DoseResponse/drc#36): `vcDisc` now uses a robust fallback chain
  (matching `vcCont`) instead of a bare
  [`solve()`](https://rdrr.io/r/base/solve.html) call. When inversion
  fails, a warning is emitted and standard errors are reported as `NA`.
- Fixed
  [`logistic()`](https://hreinwald.github.io/drc/reference/logistic.md)
  model ED calculation with `type="absolute"`: the `edfct` function now
  correctly handles absolute-to-relative conversion without applying the
  incorrect p-swap from
  [`EDhelper()`](https://hreinwald.github.io/drc/reference/EDhelper.md).
  The logistic model has opposite b-sign convention from log-logistic (b
  \< 0 means increasing, not decreasing), so `EDhelper`’s p-swap for b
  \< 0 would incorrectly swap ED values. The fix uses inline
  absolute-to-relative conversion (`p = 100·(d−respl)/(d−c)`) for
  absolute type and `p = respl` directly for relative type.
- Fixed model-level `edfct` derivatives for absolute ED type in
  [`braincousens()`](https://hreinwald.github.io/drc/reference/braincousens.md),
  [`fplogistic()`](https://hreinwald.github.io/drc/reference/fplogistic.md),
  [`llogistic()`](https://hreinwald.github.io/drc/reference/llogistic.md),
  [`llogistic2()`](https://hreinwald.github.io/drc/reference/llogistic2.md),
  [`lnormal()`](https://hreinwald.github.io/drc/reference/lnormal.md),
  [`weibull1()`](https://hreinwald.github.io/drc/reference/weibull1.md),
  and
  [`weibull2()`](https://hreinwald.github.io/drc/reference/weibull2.md):
  when `type = "absolute"`, the gradient functions previously set ∂ED/∂c
  and ∂ED/∂d to 0, which was incorrect because the absolute-to-relative
  conversion makes p a function of c and d. The chain rule requires
  non-zero partials. Now compute ∂ED/∂c and ∂ED/∂d via central
  differences on a closure that captures the chain-rule contribution
  through the full ED computation path.

------------------------------------------------------------------------

## drc 3.3.0.03

### New Features

- Enhanced
  [`plot.drc()`](https://hreinwald.github.io/drc/reference/plot.drc.md):
  error bars in `type = "bars"` plots now match curve colors by default.
  Added `errbar.col` parameter to allow manual control of error bar
  colors. Set `errbar.col = "black"` to restore the previous behavior of
  black error bars.

### Bug Fixes

- Fixed [`predict()`](https://rdrr.io/r/stats/predict.html) “incorrect
  number of dimensions” error for models with many fixed parameters
  (e.g., `EXD.3(fixed = c(lower, upper, NA))`): when only one parameter
  is estimated, `indexMat` in the fitted model object is a vector rather
  than a matrix, causing
  [`predict.drc()`](https://hreinwald.github.io/drc/reference/predict.drc.md)
  to fail when computing standard errors or confidence intervals.
  Ensured `indexMat` is always coerced to a matrix before column
  subsetting.

### Changes

- Updated package version and date in `DESCRIPTION` and website
  documentation to `3.3.0.03`.
- Updated logo path in `README.md` to point to `man/figures/logo.png`
  for consistency with package structure.
- Added favicon and manifest links to HTML documentation files for
  improved branding and browser integration.
- Added the package website (`https://hreinwald.github.io/drc`) as the
  primary URL in the `DESCRIPTION` file for better discoverability.
- Added the [`rss()`](https://hreinwald.github.io/drc/reference/rss.md)
  function to the reference index in `_pkgdown.yml`.
- Added logo image to the dose-response workflow vignette and updated
  the vignette date.
- Simplified labeling of effective dose (ED) estimates in the workflow
  vignette outputs for clarity, removing the `e:1:` prefix.
- Updated model comparison output in the vignette to include additional
  columns and more precise values.

------------------------------------------------------------------------

## drc 3.3.0.02

### New Features

- Added [`rss()`](https://hreinwald.github.io/drc/reference/rss.md)
  function for computing the residual sum of squares of a fitted `drc`
  model. Refactored
  [`Rsq()`](https://hreinwald.github.io/drc/reference/Rsq.md) to reuse
  [`rss()`](https://hreinwald.github.io/drc/reference/rss.md)
  internally; both functions are now exported.

### Bug Fixes

- Fixed [`ED()`](https://hreinwald.github.io/drc/reference/ED.md) for
  exponential decay models (EXD.2, EXD.3, AR.2, AR.3, W1.x, W2.x) with
  two fixed parameters: when only one parameter is estimated (1×1
  variance-covariance matrix), the function previously failed with
  “incorrect number of dimensions” errors. Enhanced `ED.drc` to
  defensively coerce scalar/vector `vcov` inputs to proper matrices and
  to always strip names from gradients for consistent matrix algebra.
  This fix now allows retrieving ED values from exponential decay models
  with two fixed parameters, which was previously impossible.
- Fixed gradient handling in
  [`ED()`](https://hreinwald.github.io/drc/reference/ED.md) to ensure
  model-specific derivative functions always return unnamed numeric
  vectors, preventing dimension errors in delta-method standard error
  calculations.
- Fixed boundary detection bugs in
  [`MAX()`](https://hreinwald.github.io/drc/reference/MAX.md): used
  [`unname()`](https://rdrr.io/r/base/unname.html) so named return
  values from cedergreen models are compared correctly with unnamed
  lower/upper scalars, and added tolerance in boundary check since
  numerical optimizers return values near but not exactly at boundaries.
- Fixed [`PR()`](https://hreinwald.github.io/drc/reference/PR.md)
  dropping `...` arguments for single-curve models.
- Fixed all 17 issues in
  [`ucedergreen()`](https://hreinwald.github.io/drc/reference/ucedergreen.md)
  function: missing `+c` term in model formula, `edfct` signature
  mismatch with the drc framework, undefined `xlogx` function call in
  `deriv1`, missing
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) validation for
  `method`, vectorized `|` operators in scalar `if()` guards, missing
  `useFixed` flag computation, `maxfct` signature mismatch and unsafe
  parameter indexing, broken self-starter ignoring
  `alpha`/`method`/`useFixed`, missing `fctName`/`fctText` parameters,
  `deriv1` excluded from return list, and documentation issues.
- Fixed SE calculation for absolute type
  [`ED()`](https://hreinwald.github.io/drc/reference/ED.md): the
  model-specific `edfct` gradient functions treated asymptote parameters
  as constants when `type="absolute"`, missing the chain-rule
  contribution from the `absToRel` conversion and underestimating the
  standard error. Now uses numerical central differences with an
  improved adaptive step size. Added internal helpers
  `.centralDiffGradient()`, `.safeConfintBasic()`, and `.computeSE()` to
  make SE computation more robust: `.computeSE()` guards against
  non-positive-definite variance-covariance matrix slices (returning
  `NA` instead of erroring), and `.safeConfintBasic()` validates
  residual degrees of freedom before calling
  [`confint.basic()`](https://hreinwald.github.io/drc/reference/confint.basic.md),
  falling back to a z-distribution when
  [`df.residual()`](https://rdrr.io/r/stats/df.residual.html) returns an
  invalid value.
- Fixed inverted `otrace`/`silentVal` logic in
  [`drmOpt()`](https://hreinwald.github.io/drc/reference/drmOpt.md)
  where `otrace=TRUE` incorrectly caused `silent=TRUE` in
  `try(optim())`, suppressing error messages instead of displaying them.
- Fixed
  [`searchdrc()`](https://hreinwald.github.io/drc/reference/searchdrc.md)
  regex error and convergence failure behavior.
- Fixed citation URL: reordered URLs in DESCRIPTION so `citation('drc')`
  returns the GitHub repository URL instead of r-project.org.
- Fixed [`ED()`](https://hreinwald.github.io/drc/reference/ED.md)
  “incorrect number of dimensions” error for models with few estimated
  parameters (e.g., EXD.3 with fixed c and d): ensured `indexMat` is
  always treated as a matrix before column subsetting.
- Fixed [`ED()`](https://hreinwald.github.io/drc/reference/ED.md)
  returning NaN with warning for LL.5 models with ill-conditioned
  parameters: added validity check to return `Inf` (indicating EC50 is
  outside valid range) instead of NaN when
  `exp(-tempVal/parmVec[5]) - 1` is non-positive. Also fixed NaN
  handling in the check condition to prevent “missing value where
  TRUE/FALSE needed” errors in
  [`backfit()`](https://hreinwald.github.io/drc/reference/backfit.md)
  and other functions.
- Fixed additional robustness issues in
  [`ED()`](https://hreinwald.github.io/drc/reference/ED.md) / `ED.drc`:
  loop now always iterates over all curves and all response levels,
  filtering by `clevel` after computation rather than before;
  `invMatList` is grown dynamically to avoid NULL holes; curve label
  construction uses a single structured object with explicit `match` and
  `display` fields; variance-covariance matrix slices always use
  `drop = FALSE` to remain matrices.
- Fixed
  [`mselect()`](https://hreinwald.github.io/drc/reference/mselect.md)
  missing two closing braces that caused a parse error when the function
  was sourced directly.
- Fixed `ED.lin.R` bugs: removed a duplicate `if`-block (dead code that
  evaluated the same condition twice), removed a stray debug
  [`print()`](https://rdrr.io/r/base/print.html) statement, and added
  the missing `parameterNames = c("b0", "b1", "b2")` argument to the
  `deltaMethod()` call for quadratic models (the omission caused
  incorrect parameter mapping and wrong confidence intervals).
- Fixed
  [`CRS.4b()`](https://hreinwald.github.io/drc/reference/CRS.4b.md)
  display text: `fctText` incorrectly showed `"alpha="` instead of
  `"alpha=0.5"`.
- Fixed
  [`gammadr()`](https://hreinwald.github.io/drc/reference/gammadr.md)
  first-derivative (`deriv1`) calculation: the gradient with respect to
  the dose parameter incorrectly used `parmMat[, 1]` (the rate
  parameter) where `dose` was required, producing wrong gradient values.
- Fixed [`maED()`](https://hreinwald.github.io/drc/reference/maED.md)
  model-averaging: models whose ED estimates are non-finite (`Inf` or
  `NaN`) are now detected and excluded from the weighted average (with a
  warning naming the model and the offending values); models that
  returned a `try-error` during fitting are also excluded. When all
  candidate models are excluded, the function returns `NA` for all
  estimates instead of `0` or `NaN`.
- Added warning to
  [`noEffect()`](https://hreinwald.github.io/drc/reference/noEffect.md)
  when degrees of freedom difference is ≤ 0, clarifying that the
  likelihood ratio test may not be meaningful when the dose-response
  model has no additional parameters compared to the null model (e.g.,
  when most parameters are fixed).

### Changes

- Added `NEWS.md` version control log. Reformatted legacy news file into
  properly formatted `NEWS.md` with categorized sections.
- Improved documentation for Weibull starting value `method` parameter
  across
  [`weibull1()`](https://hreinwald.github.io/drc/reference/weibull1.md),
  [`weibull2()`](https://hreinwald.github.io/drc/reference/weibull2.md),
  and all wrapper functions (`W1.2`, `W1.3`, `W1.4`, `W2.2`, `W2.3`,
  `W2.4`, `AR.2`, `AR.3`, `EXD.2`, `EXD.3`).
- Enhanced roxygen2 documentation for `ED` and `ED.drc` functions with
  improved parameter descriptions and examples.
- Added comprehensive test suites for `anova.drclist`,`summary.drc`,
  `print.summary.drc`, `noEffect`, `searchdrc`, `backfit`, `getInitial`,
  `drmEMeventtime`, `repChar`, `rdrm`, `gompertzd`,
  [`MAX()`](https://hreinwald.github.io/drc/reference/MAX.md), and
  [`PR()`](https://hreinwald.github.io/drc/reference/PR.md) functions.
- Added comprehensive test suites for `llogistic`/LL.x models,
  `weibull1`/W1.x/EXD.x models, `logistic.ssf`, `gammadr`, `EDcomp`,
  `mselect`, `drmOpt`, `modelFunction`, `modelFit`, `anova.drclist`,
  `rss`, and `ED.lin`.
- Large-scale dead code removal across 70+ R source files: removed
  commented-out function implementations, stray
  [`print()`](https://rdrr.io/r/base/print.html) debug statements, old
  code paths, and `if(FALSE){...}` blocks. No logic changes; all
  roxygen2 documentation and meaningful explanatory comments were
  preserved.
- Removed dead code `iband.R` and all associated references.
- Removed unused `inst/citation` file, superseded by `CITATION.cff` at
  repository root.
- Deleted `build_pkgdown.R` build script.
- Added PLoS ONE 2015 article and CRC Press 2019 book references to
  `CITATION.cff`.
- Updated installation instructions and README documentation.
- Added `magic` to Suggests in DESCRIPTION for test dependency.

------------------------------------------------------------------------

## drc 3.3.0.01

### New Features

- Created comprehensive vignettes: `dose-response-workflow.Rmd`
  providing a complete tutorial on dose-response analysis, and
  `nec-models.Rmd` documenting No Effect Concentration modeling with
  `NEC.2`/`NEC.3`/`NEC.4` function variants.
- Set up pkgdown website infrastructure: added `_pkgdown.yml` with
  Bootstrap 5 configuration, created `build_pkgdown.R` script for build
  automation, documented pkgdown build process in README, and generated
  pkgdown documentation site.
- Added computationally robust (stable) wrapper functions in new
  `ED_robust.R` module:
  [`ED_robust()`](https://hreinwald.github.io/drc/reference/ED_robust.md)
  for calculating ED values with proper error handling that returns `NA`
  instead of failing when an ED value is not estimable,
  [`maED_robust()`](https://hreinwald.github.io/drc/reference/maED_robust.md)
  for model-averaged ED estimation with the same graceful error
  handling, and
  [`get_ed_interval()`](https://hreinwald.github.io/drc/reference/get_ed_interval.md)
  for recommending appropriate confidence interval methods based on
  model type.
- Added comprehensive test suite covering ED calculations, predictions,
  plotting, residuals, model selection, and utility functions.
- Added `drm_name()` helper function to `ED_robust.R`.
- Enhanced package startup message with citations and developer credits.
- Added
  [`drm_legacy()`](https://hreinwald.github.io/drc/reference/drm_legacy.md)
  as an internal reference function preserving the original
  [`drm()`](https://hreinwald.github.io/drc/reference/drm.md)
  implementation.
- Added testthat infrastructure with tests verifying
  [`drm()`](https://hreinwald.github.io/drc/reference/drm.md) output
  matches
  [`drm_legacy()`](https://hreinwald.github.io/drc/reference/drm_legacy.md)
  output across continuous, binomial, Poisson, and negative binomial
  data types.
- Added comprehensive anova tests.

### Bug Fixes

- Fixed vignette build by removing vignettes from `.Rbuildignore` and
  correcting incorrect
  [`mselect()`](https://hreinwald.github.io/drc/reference/mselect.md)
  usage in examples.
- Fixed Rd comment warning by escaping the `%*%` operator in
  documentation.
- Fixed all
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  errors and warnings: added roxygen2 `@keywords` and lifecycle
  deprecation notices for deprecated CRS functions, expanded dataset
  documentation files with examples, added missing dataset aliases, and
  fixed Weibull model documentation.
- Fixed division-by-zero in
  [`Rsq()`](https://hreinwald.github.io/drc/reference/Rsq.md) and
  [`absToRel()`](https://hreinwald.github.io/drc/reference/absToRel.md).
- Removed dead `scaleEst()` stub function.
- Fixed [`inherits()`](https://rdrr.io/r/base/class.html) bug in
  `mselect.R`.
- Added edge case handling in `modelFit.R`.
- Added input validation for
  [`comped()`](https://hreinwald.github.io/drc/reference/comped.md) and
  [`compParm()`](https://hreinwald.github.io/drc/reference/compParm.md).
- Fixed unsafe global state modification via `options(warn)`, incorrect
  `compParm` od/pool handling, and residuals division by zero.
- Fixed NaN warning in `summary.drc` for robust estimation methods
  (metric trimming, Winsorizing, Tukey’s biweight).
- Improved `predict.drc` and `vcov.drc` to resolve 23 test failures.
- Fixed
  [`mselect()`](https://hreinwald.github.io/drc/reference/mselect.md) to
  always compute Lack of fit p-values for all models, not only when
  `nested=TRUE`.
- Fixed a bug in `anova.drclist` where negative or non-finite F
  statistics produced NaN p-values; negative F statistics now return
  p-value of 1 and non-finite F statistics return NA.
- Fixed duplicate aliases and unstated dependencies in examples.
- Fixed package dependency warnings: added `data.table` and `dplyr` to
  Imports, updated NAMESPACE with required imports.
- Fixed S3 method consistency issues and changed `confint.basic` roxygen
  tag from `@exportS3Method` to `@export`.
- Fixed escaped LaTeX special characters in roxygen2 documentation and
  Rd files.
- Fixed escaped percent signs in roxygen docs causing Rd parse warnings.

### Changes

- Added vignette access information to README.
- Completed comprehensive roxygen2 documentation audit: added missing
  `@param` tags, removed `dontrun`/`donttest` wrappers to enable
  automated example testing, and fixed broken examples across
  documentation files.
- Enhanced dataset documentation: improved descriptions and fixed typos
  in dataset `.Rd` files, added examples sections to dataset `.Rd` files
  that were missing them.
- Improved `confint.drc` robustness: added
  [`stop()`](https://rdrr.io/r/base/stop.html) fallback to
  [`switch()`](https://rdrr.io/r/base/switch.html) in
  [`confint.basic()`](https://hreinwald.github.io/drc/reference/confint.basic.md)
  to handle unknown `intType` values gracefully instead of returning
  silent NULL.
- Removed `@export` from
  [`confint.basic()`](https://hreinwald.github.io/drc/reference/confint.basic.md)
  as internal helpers should not be part of the public API.
- Enhanced roxygen2 documentation for `CRS.5`, convenience functions,
  and `ED_robust` with improved argument descriptions.
- Updated DESCRIPTION: added Hannes Reinwald as maintainer and
  co-author, updated package version to 3.3.0.01.
- Removed external `drcData` package dependency; example datasets are
  now bundled directly in the package `data/` directory.
- Added `.Rd` documentation files for all bundled datasets.
- Renamed internal variables for clarity: `ndRows` to `nRows` in
  `predict.drc`, `posIdx` to `validVar` in `summary.drc`.
- Added test coverage documentation.
- Renamed all 33 R source files from lowercase `.r` to uppercase `.R`
  extensions for consistency.
- Updated all `.Rd` documentation files to reference the new file names.
- Added GNU General Public License version 2 file and updated license
  version from GPL-2 to GPL-2.0 in DESCRIPTION.
- Added Hannes Reinwald as author in the DESCRIPTION file.
- Updated README with revised installation instructions and bug report
  link.
- Migrated all package documentation to roxygen2-generated Rd files.
- Regenerated NAMESPACE via roxygen2.
- Added `@exportS3Method` tags to S3 methods in `confint.drc.R` and
  `mrdrm.r`.
- Updated package version format from 3.3-0 to 3.3.0.
- Lowered the minimum R version requirement to 4.0.0.

### Breaking Changes

- Removed deprecated developmental `cedergreen2` function.

------------------------------------------------------------------------

## drc 3.3.0

### New Features

- Added new `CRS.5` wrapper function and `CRS.6` six-parameter model
  where the alpha exponent is estimated rather than fixed.

### Bug Fixes

- Fixed a bug where the [`stop()`](https://rdrr.io/r/base/stop.html)
  call for using separate curves with control measurements was inside
  the `if(!noMessage)` block, meaning it would be silently skipped when
  messages were suppressed.
- Fixed a bug in `noEffect.R` where the Poisson null model incorrectly
  referenced `resp` instead of using the response vector from the fitted
  object.

### Changes

- Refactored the Cedergreen-Ritz-Streibig hormesis model: extracted
  `edfct` and `maxfct` into standalone helper functions
  (`cedergreen_edfct`, `cedergreen_maxfct`), refactored the self-starter
  function, and improved documentation.
- Cleaned up [`drm()`](https://hreinwald.github.io/drc/reference/drm.md)
  function by removing approximately 900 lines of commented-out dead
  code, debug print statements, and old experimental implementations.
- Removed unused variable `isfi` and redundant variable `lenData`
  (identical to `numObs`).
- Removed a dead loop over `pmodelsList2` that could never execute.
- Added roxygen2 documentation headers to all R source files across the
  package.
- Fixed typos in source code and manual pages: ‘insted’ to ‘instead’ in
  `gaussian.r` and `lgaussian.R`, duplicate parameter name ‘e1’ to ‘e2’
  in `ursa.r`, ‘contain’ to ‘contents’ in `EDcomp.R`, ‘mising’ to
  ‘missing’ and ‘reponses’ to ‘responses’ in `drm.Rd`, and ‘reponse’ to
  ‘response’ in `CRS.5a.Rd`.
- Updated DESCRIPTION file: added Encoding field (UTF-8), fixed
  `Authors@R` to use proper
  [`person()`](https://rdrr.io/r/utils/person.html) format, removed
  deprecated Maintainer and LazyLoad fields, added missing Imports
  (`graphics`, `utils`), and updated URLs from HTTP to HTTPS.
- Comprehensive repository cleanup and code quality improvements.
- Removed obsolete configuration files (`.travis.yml`, `drc.Rproj`,
  `_pkgdown.yml`, `README.Rmd`) and redundant reference files
  (`_gitignore`, `_Rbuildignore`).
- Removed the `/tests` directory containing outdated development
  artifacts with no testing value.
- Updated `.gitignore` and `.Rbuildignore` with standard R/RStudio
  settings.
- Rewrote `README.md` with comprehensive documentation including
  quick-start examples, available models, key functions, and supported
  data types.
- Removed debug [`print()`](https://rdrr.io/r/base/print.html)
  statements in `drmEMstandard.R` and `findbe.r`.
- Removed dead code block (commented-out experimental code wrapped in
  `if (FALSE)`) in `drmEMstandard.R` and `llogistic.ssf.R`.
- Replaced unsafe `eval(parse(text=...))` calls with
  [`match.fun()`](https://rdrr.io/r/base/match.fun.html) and
  [`do.call()`](https://rdrr.io/r/base/do.call.html) in `rdrm.r`.
- Improved [`options()`](https://rdrr.io/r/base/options.html) handling
  in `searchdrc.R` by saving and restoring the original warn setting
  using `on.exit(add=TRUE)`.

### Deprecated

- Deprecated old CRS function names (`CRS.4a`, `CRS.4b`, `CRS.4c`,
  `CRS.5a`, `CRS.5b`, `CRS.5c`) with lifecycle notices in favor of new
  wrappers.

------------------------------------------------------------------------

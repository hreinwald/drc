# drc 3.3.0.03

## New Features
* Enhanced `plot.drc()`: error bars in `type = "bars"` plots now match curve colors by default. Added `errbar.col` parameter to allow manual control of error bar colors. Set `errbar.col = "black"` to restore the previous behavior of black error bars.

## Bug Fixes
* Fixed `predict()` "incorrect number of dimensions" error for models with many fixed parameters (e.g., `EXD.3(fixed = c(lower, upper, NA))`): when only one parameter is estimated, `indexMat` in the fitted model object is a vector rather than a matrix, causing `predict.drc()` to fail when computing standard errors or confidence intervals. Ensured `indexMat` is always coerced to a matrix before column subsetting.

## Changes
* Updated package version and date in `DESCRIPTION` and website documentation to `3.3.0.03`.
* Updated logo path in `README.md` to point to `man/figures/logo.png` for consistency with package structure.
* Added favicon and manifest links to HTML documentation files for improved branding and browser integration.
* Added the package website (`https://hreinwald.github.io/drc`) as the primary URL in the `DESCRIPTION` file for better discoverability.
* Added the `rss()` function to the reference index in `_pkgdown.yml`.
* Added logo image to the dose-response workflow vignette and updated the vignette date.
* Simplified labeling of effective dose (ED) estimates in the workflow vignette outputs for clarity, removing the `e:1:` prefix.
* Updated model comparison output in the vignette to include additional columns and more precise values.

---

# drc 3.3.0.02

## New Features
* Added `rss()` function for computing the residual sum of squares of a fitted `drc` model. Refactored `Rsq()` to reuse `rss()` internally; both functions are now exported.

## Bug Fixes
* Fixed `ED()` for exponential decay models (EXD.2, EXD.3, AR.2, AR.3, W1.x, W2.x) with two fixed parameters: when only one parameter is estimated (1×1 variance-covariance matrix), the function previously failed with "incorrect number of dimensions" errors. Enhanced `ED.drc` to defensively coerce scalar/vector `vcov` inputs to proper matrices and to always strip names from gradients for consistent matrix algebra. This fix now allows retrieving ED values from exponential decay models with two fixed parameters, which was previously impossible.
* Fixed gradient handling in `ED()` to ensure model-specific derivative functions always return unnamed numeric vectors, preventing dimension errors in delta-method standard error calculations.
* Fixed boundary detection bugs in `MAX()`: used `unname()` so named return values from cedergreen models are compared correctly with unnamed lower/upper scalars, and added tolerance in boundary check since numerical optimizers return values near but not exactly at boundaries.
* Fixed `PR()` dropping `...` arguments for single-curve models.
* Fixed all 17 issues in `ucedergreen()` function: missing `+c` term in model formula, `edfct` signature mismatch with the drc framework, undefined `xlogx` function call in `deriv1`, missing `match.arg()` validation for `method`, vectorized `|` operators in scalar `if()` guards, missing `useFixed` flag computation, `maxfct` signature mismatch and unsafe parameter indexing, broken self-starter ignoring `alpha`/`method`/`useFixed`, missing `fctName`/`fctText` parameters, `deriv1` excluded from return list, and documentation issues.
* Fixed SE calculation for absolute type `ED()`: the model-specific `edfct` gradient functions treated asymptote parameters as constants when `type="absolute"`, missing the chain-rule contribution from the `absToRel` conversion and underestimating the standard error. Now uses numerical central differences with an improved adaptive step size. Added internal helpers `.centralDiffGradient()`, `.safeConfintBasic()`, and `.computeSE()` to make SE computation more robust: `.computeSE()` guards against non-positive-definite variance-covariance matrix slices (returning `NA` instead of erroring), and `.safeConfintBasic()` validates residual degrees of freedom before calling `confint.basic()`, falling back to a z-distribution when `df.residual()` returns an invalid value.
* Fixed inverted `otrace`/`silentVal` logic in `drmOpt()` where `otrace=TRUE` incorrectly caused `silent=TRUE` in `try(optim())`, suppressing error messages instead of displaying them.
* Fixed `searchdrc()` regex error and convergence failure behavior.
* Fixed citation URL: reordered URLs in DESCRIPTION so `citation('drc')` returns the GitHub repository URL instead of r-project.org.
* Fixed `ED()` "incorrect number of dimensions" error for models with few estimated parameters (e.g., EXD.3 with fixed c and d): ensured `indexMat` is always treated as a matrix before column subsetting.
* Fixed `ED()` returning NaN with warning for LL.5 models with ill-conditioned parameters: added validity check to return `Inf` (indicating EC50 is outside valid range) instead of NaN when `exp(-tempVal/parmVec[5]) - 1` is non-positive. Also fixed NaN handling in the check condition to prevent "missing value where TRUE/FALSE needed" errors in `backfit()` and other functions.
* Fixed additional robustness issues in `ED()` / `ED.drc`: loop now always iterates over all curves and all response levels, filtering by `clevel` after computation rather than before; `invMatList` is grown dynamically to avoid NULL holes; curve label construction uses a single structured object with explicit `match` and `display` fields; variance-covariance matrix slices always use `drop = FALSE` to remain matrices.
* Fixed `mselect()` missing two closing braces that caused a parse error when the function was sourced directly.
* Fixed `ED.lin.R` bugs: removed a duplicate `if`-block (dead code that evaluated the same condition twice), removed a stray debug `print()` statement, and added the missing `parameterNames = c("b0", "b1", "b2")` argument to the `deltaMethod()` call for quadratic models (the omission caused incorrect parameter mapping and wrong confidence intervals).
* Fixed `CRS.4b()` display text: `fctText` incorrectly showed `"alpha="` instead of `"alpha=0.5"`.
* Fixed `gammadr()` first-derivative (`deriv1`) calculation: the gradient with respect to the dose parameter incorrectly used `parmMat[, 1]` (the rate parameter) where `dose` was required, producing wrong gradient values.
* Fixed `maED()` model-averaging: models whose ED estimates are non-finite (`Inf` or `NaN`) are now detected and excluded from the weighted average (with a warning naming the model and the offending values); models that returned a `try-error` during fitting are also excluded. When all candidate models are excluded, the function returns `NA` for all estimates instead of `0` or `NaN`.
* Added warning to `noEffect()` when degrees of freedom difference is ≤ 0, clarifying that the likelihood ratio test may not be meaningful when the dose-response model has no additional parameters compared to the null model (e.g., when most parameters are fixed).

## Changes
* Added `NEWS.md` version control log. Reformatted legacy news file into properly formatted `NEWS.md` with categorized sections.
* Improved documentation for Weibull starting value `method` parameter across `weibull1()`, `weibull2()`, and all wrapper functions (`W1.2`, `W1.3`, `W1.4`, `W2.2`, `W2.3`, `W2.4`, `AR.2`, `AR.3`, `EXD.2`, `EXD.3`).
* Enhanced roxygen2 documentation for `ED` and `ED.drc` functions with improved parameter descriptions and examples.
* Added comprehensive test suites for ``anova.drclist``,`summary.drc`, `print.summary.drc`, `noEffect`, `searchdrc`, `backfit`, `getInitial`, `drmEMeventtime`, `repChar`, `rdrm`, `gompertzd`, `MAX()`, and `PR()` functions.
* Added comprehensive test suites for `llogistic`/LL.x models, `weibull1`/W1.x/EXD.x models, `logistic.ssf`, `gammadr`, `EDcomp`, `mselect`, `drmOpt`, `modelFunction`, `modelFit`, `anova.drclist`, `rss`, and `ED.lin`.
* Large-scale dead code removal across 70+ R source files: removed commented-out function implementations, stray `print()` debug statements, old code paths, and `if(FALSE){...}` blocks. No logic changes; all roxygen2 documentation and meaningful explanatory comments were preserved.
* Removed dead code `iband.R` and all associated references.
* Removed unused `inst/citation` file, superseded by `CITATION.cff` at repository root.
* Deleted `build_pkgdown.R` build script.
* Added PLoS ONE 2015 article and CRC Press 2019 book references to `CITATION.cff`.
* Updated installation instructions and README documentation.
* Added `magic` to Suggests in DESCRIPTION for test dependency.

---

# drc 3.3.0.01

## New Features
* Created comprehensive vignettes: `dose-response-workflow.Rmd` providing a complete tutorial on dose-response analysis, and `nec-models.Rmd` documenting No Effect Concentration modeling with `NEC.2`/`NEC.3`/`NEC.4` function variants.
* Set up pkgdown website infrastructure: added `_pkgdown.yml` with Bootstrap 5 configuration, created `build_pkgdown.R` script for build automation, documented pkgdown build process in README, and generated pkgdown documentation site.
* Added computationally robust (stable) wrapper functions in new `ED_robust.R` module: `ED_robust()` for calculating ED values with proper error handling that returns `NA` instead of failing when an ED value is not estimable, `maED_robust()` for model-averaged ED estimation with the same graceful error handling, and `get_ed_interval()` for recommending appropriate confidence interval methods based on model type.
* Added comprehensive test suite covering ED calculations, predictions, plotting, residuals, model selection, and utility functions.
* Added `drm_name()` helper function to `ED_robust.R`.
* Enhanced package startup message with citations and developer credits.
* Added `drm_legacy()` as an internal reference function preserving the original `drm()` implementation.
* Added testthat infrastructure with tests verifying `drm()` output matches `drm_legacy()` output across continuous, binomial, Poisson, and negative binomial data types.
* Added comprehensive anova tests.

## Bug Fixes
* Fixed vignette build by removing vignettes from `.Rbuildignore` and correcting incorrect `mselect()` usage in examples.
* Fixed Rd comment warning by escaping the `%*%` operator in documentation.
* Fixed all `devtools::check()` errors and warnings: added roxygen2 `@keywords` and lifecycle deprecation notices for deprecated CRS functions, expanded dataset documentation files with examples, added missing dataset aliases, and fixed Weibull model documentation.
* Fixed division-by-zero in `Rsq()` and `absToRel()`.
* Removed dead `scaleEst()` stub function.
* Fixed `inherits()` bug in `mselect.R`.
* Added edge case handling in `modelFit.R`.
* Added input validation for `comped()` and `compParm()`.
* Fixed unsafe global state modification via `options(warn)`, incorrect `compParm` od/pool handling, and residuals division by zero.
* Fixed NaN warning in `summary.drc` for robust estimation methods (metric trimming, Winsorizing, Tukey's biweight).
* Improved `predict.drc` and `vcov.drc` to resolve 23 test failures.
* Fixed `mselect()` to always compute Lack of fit p-values for all models, not only when `nested=TRUE`.
* Fixed a bug in `anova.drclist` where negative or non-finite F statistics produced NaN p-values; negative F statistics now return p-value of 1 and non-finite F statistics return NA.
* Fixed duplicate aliases and unstated dependencies in examples.
* Fixed package dependency warnings: added `data.table` and `dplyr` to Imports, updated NAMESPACE with required imports.
* Fixed S3 method consistency issues and changed `confint.basic` roxygen tag from `@exportS3Method` to `@export`.
* Fixed escaped LaTeX special characters in roxygen2 documentation and Rd files.
* Fixed escaped percent signs in roxygen docs causing Rd parse warnings.

## Changes
* Added vignette access information to README.
* Completed comprehensive roxygen2 documentation audit: added missing `@param` tags, removed `dontrun`/`donttest` wrappers to enable automated example testing, and fixed broken examples across documentation files.
* Enhanced dataset documentation: improved descriptions and fixed typos in dataset `.Rd` files, added examples sections to dataset `.Rd` files that were missing them.
* Improved `confint.drc` robustness: added `stop()` fallback to `switch()` in `confint.basic()` to handle unknown `intType` values gracefully instead of returning silent NULL.
* Removed `@export` from `confint.basic()` as internal helpers should not be part of the public API.
* Enhanced roxygen2 documentation for `CRS.5`, convenience functions, and `ED_robust` with improved argument descriptions.
* Updated DESCRIPTION: added Hannes Reinwald as maintainer and co-author, updated package version to 3.3.0.01.
* Removed external `drcData` package dependency; example datasets are now bundled directly in the package `data/` directory.
* Added `.Rd` documentation files for all bundled datasets.
* Renamed internal variables for clarity: `ndRows` to `nRows` in `predict.drc`, `posIdx` to `validVar` in `summary.drc`.
* Added test coverage documentation.
* Renamed all 33 R source files from lowercase `.r` to uppercase `.R` extensions for consistency.
* Updated all `.Rd` documentation files to reference the new file names.
* Added GNU General Public License version 2 file and updated license version from GPL-2 to GPL-2.0 in DESCRIPTION.
* Added Hannes Reinwald as author in the DESCRIPTION file.
* Updated README with revised installation instructions and bug report link.
* Migrated all package documentation to roxygen2-generated Rd files.
* Regenerated NAMESPACE via roxygen2.
* Added `@exportS3Method` tags to S3 methods in `confint.drc.R` and `mrdrm.r`.
* Updated package version format from 3.3-0 to 3.3.0.
* Lowered the minimum R version requirement to 4.0.0.

## Breaking Changes
* Removed deprecated developmental `cedergreen2` function.

---

# drc 3.3.0

## New Features
* Added new `CRS.5` wrapper function and `CRS.6` six-parameter model where the alpha exponent is estimated rather than fixed.

## Bug Fixes
* Fixed a bug where the `stop()` call for using separate curves with control measurements was inside the `if(!noMessage)` block, meaning it would be silently skipped when messages were suppressed.
* Fixed a bug in `noEffect.R` where the Poisson null model incorrectly referenced `resp` instead of using the response vector from the fitted object.

## Changes
* Refactored the Cedergreen-Ritz-Streibig hormesis model: extracted `edfct` and `maxfct` into standalone helper functions (`cedergreen_edfct`, `cedergreen_maxfct`), refactored the self-starter function, and improved documentation.
* Cleaned up `drm()` function by removing approximately 900 lines of commented-out dead code, debug print statements, and old experimental implementations.
* Removed unused variable `isfi` and redundant variable `lenData` (identical to `numObs`).
* Removed a dead loop over `pmodelsList2` that could never execute.
* Added roxygen2 documentation headers to all R source files across the package.
* Fixed typos in source code and manual pages: 'insted' to 'instead' in `gaussian.r` and `lgaussian.R`, duplicate parameter name 'e1' to 'e2' in `ursa.r`, 'contain' to 'contents' in `EDcomp.R`, 'mising' to 'missing' and 'reponses' to 'responses' in `drm.Rd`, and 'reponse' to 'response' in `CRS.5a.Rd`.
* Updated DESCRIPTION file: added Encoding field (UTF-8), fixed `Authors@R` to use proper `person()` format, removed deprecated Maintainer and LazyLoad fields, added missing Imports (`graphics`, `utils`), and updated URLs from HTTP to HTTPS.
* Comprehensive repository cleanup and code quality improvements.
* Removed obsolete configuration files (`.travis.yml`, `drc.Rproj`, `_pkgdown.yml`, `README.Rmd`) and redundant reference files (`_gitignore`, `_Rbuildignore`).
* Removed the `/tests` directory containing outdated development artifacts with no testing value.
* Updated `.gitignore` and `.Rbuildignore` with standard R/RStudio settings.
* Rewrote `README.md` with comprehensive documentation including quick-start examples, available models, key functions, and supported data types.
* Removed debug `print()` statements in `drmEMstandard.R` and `findbe.r`.
* Removed dead code block (commented-out experimental code wrapped in `if (FALSE)`) in `drmEMstandard.R` and `llogistic.ssf.R`.
* Replaced unsafe `eval(parse(text=...))` calls with `match.fun()` and `do.call()` in `rdrm.r`.
* Improved `options()` handling in `searchdrc.R` by saving and restoring the original warn setting using `on.exit(add=TRUE)`.

## Deprecated
* Deprecated old CRS function names (`CRS.4a`, `CRS.4b`, `CRS.4c`, `CRS.5a`, `CRS.5b`, `CRS.5c`) with lifecycle notices in favor of new wrappers.

---

# drc Changes in 2017

## New Features
* The argument `checkND` has been added to the predict method, allowing switching off comparison of variable names in the original data frame and the `newdata` data frame; useful for predicting in mixture models (after a report by Evan Palmer-Young).
* Confidence intervals for ED values may now be obtained using inverse regression (`interval = "inv"`).
* Species sensitivity distributions may now be fitted using `drm()` with `type = "ssd"`. The predict method now constrains predicted values to meaningful ranges by default and allows incorporating standard errors of estimates in confidence bands for fitted SSDs.

## Bug Fixes
* Small bug in `mixture()` resolved (after a report by Andrew Kniss).

## Changes
* Updated the event-time part of drc, in particular in `drm()` (improved code provided by Andrea Onofri).

---

# drc Changes in 2016

## New Features
* Negative binomial distributions may now be fitted using the argument `type` with values `"negbin1"` and `"negbin2"` (after a suggestion from Signe M. Jensen).
* Argument `conCheck` added to `drmc()` to switch on/off handling of control measurements.
* `na.omit()` is now the default in `drm()`.
* New functions `CIcomp()`, `CIcompX()`, and `plotFACI()` for calculating combination indices based on effective doses and effects as described in Martin-Betancor et al. (2015). An accompanying dataset `metals` has also been included.

## Bug Fixes
* Fixed a small bug in printing confidence intervals (after a report from Johannes Ranke).
* Fixed a small error for binomial data in the function `estfun.drc` for use with the package "sandwich" (after a report from Andrew Kniss).

## Changes
* Help page for two- and three-parameter Weibull models updated, removing typos (after a report by Mikael Gustavsson).
* Output text for confidence intervals modified slightly.

---

# drc Changes in 2015

## New Features
* The plot method extended to provide confidence bands (contribution by Gregory Warnes).

## Bug Fixes
* Robust estimation is working again (after reports from Kathy Mutambanengwe, Sten Ilmjärv, and Corina Dueñas Roca).
* Minor labelling issue in the plot method resolved (with help from Bert Oosthuyse).

## Changes
* The predict method has been updated (after a report from Duncan Mackay).
* The function `isobole()` now also propagates graphical arguments to fitted isoboles.
* The `print.summary.drc` method shows a warning message in case of df < 1.
* Calculation of the log likelihood for binomial data has been updated and improved.
* The function `backfit()` has been updated.
* Help page for `EDcomp()` updated.
* Help pages for `ED()` and `isobole()` have been updated. Some code has been tidied up.

## Breaking Changes
* The function `SI()` has been completely replaced by the function `EDcomp()`.
* The function `diagnostics()` has been removed.

---

# drc Changes in 2014

## New Features
* The argument `pshifts` has been added to `drm()` to allow weights on parameters (after a comment from Florent Baty).
* The model functions `gammadr()` and `multi2()` have been added.
* The argument `vcov.` has been added to `compParm()`.
* The argument `vcov.` has been added to `EDcomp()` and `predict()`.
* The argument `vcov.` has been added to `ED()` to allow choosing between the standard vcov method and the sandwich function for robust standard errors.
* Added `fixed` value to output from `logistic()` (after a suggestion from Daniel Gerhard).

## Bug Fixes
* The argument `control` is now correctly propagated in case `separate = TRUE` in `drm()` (after a report from Andy Liaw).
* Small bug in the plot method (not ordering labels in legend text correctly) has been fixed (after a report by Francois Keck).

## Changes
* The help page for `ursa()` has been improved (in particular the example section).
* The bread and estfun methods have been extended to event-time data. The help page for `G.aparine` has been improved.
* A number of functions depending on `ED()` and `SI()` have been updated. The predict method and plot functionality for event-time data have been improved (after a report from Eshagh Keshtkar).
* `ED()` and `SI()` now return (invisibly) a list whose second component can be used directly with the package multcomp (after a suggestion by Daniel Gerhard).

## Breaking Changes
* A number of not fully implemented model functions for mixture data have been removed.
* `SI()` has been renamed to `EDcomp()`.

---

# drc Changes in 2013

## New Features
* An argument for specifying the reference for the normalization has been added (after a suggestion from Sunniva Foerster).
* An argument for showing normalized data and fitted curves has been added to the plot method (after a suggestion from Ludwig A. Hothorn).

## Bug Fixes
* Plot method now works for fits obtained from `drm()` using `separate = TRUE` (after a report from Thomas Kroeber).

## Changes
* Help page for `germination` has been updated.
* `bread.drc()` and `estfun.drc()` have been updated to handle fits for binomial and Poisson data (after a suggestion by Signe M. Jensen).
* Help page for `plot.drc()` has been improved regarding the explanation on the use of error bars (after an enquiry from Julien Delafontaine).
* Help page for `selenium` has been improved (after a comment from Keith Taulbee).
* Help page for `mselect()` has been improved (after a comment from Sona Jesenska).
* `cedergreen()` has been improved to provide meaningful names (after a suggestion from Dave Smithson).

---

# drc Changes in 2012

## New Features
* Model functions `gaussian()` and `lgaussian()` have been included (after a suggestion from Ismael Rodea).
* `hatvalues` and `cooks.distance` methods have been added (after a question from Sunniva Förster).
* `noEffect()` function included for testing the dose-response model against a simpler model with no dose effect (after a suggestion by Ryan Hechinger).
* `backfit()` function added (after a suggestion by Keld Sorensen).
* Poisson models can now be fitted with weights (after a suggestion by Marie Laure Delignette-Muller).
* `ED()` now also works for Gompertz models (after a question from Calvin Odero). The dataset `germination` has been included.
* The dataset `chickweed` has been included.
* The dataset `selenium` has been included.

## Bug Fixes
* Small bug in `drm()` related to starting values for event-time models resolved.
* Small error in `plot.drc` for event times fixed (after a report from Christian Andreasen).
* Small scaling error in fitted method has been removed (after a report from Andreas Betz).
* Limit in dose scaling has been lifted (after a report from Andreas Wernitznig).

## Changes
* The method for residuals has been extended to provide residuals on the transformed scale in case a Box-Cox transformation was applied.
* Help page for `maED()` has been extended.
* `vcov.drc()` has been updated.
* Help page of `lettuce` has been revised.
* Help page for the function `modelFit()` has been revised, removing a typo in the title (after a comment from John Lynch).
* Help page for the function `mr.test()` has been updated.
* Summary output for `separate = TRUE` now shows the original labels (levels in the variable curveid) (after a comment from Radu Slobodeanu).
* Help page of `NEC()` has been updated.

---

# drc Changes in 2011

## New Features
* The model function `W2x.4()` has been added (after a suggestion by Cécile Cornou).
* `display` and `type` arguments have been added to `maED()`.
* The argument `type` can now also take the value `"event"` for fitting event times.
* The model function `weibull2x()` (a model including a sort of "lag time" parameter) has been added (after a suggestion by Cécile Cornou).
* Functions `iceLoewe1()` and `iceLoewe2.1()` have been added.

## Bug Fixes
* Error in calculation of ED values for `fplogistic()` has been fixed.
* Error in label ordering and calculation of standard errors in `SI()` have been fixed (bug report by Andrew Kniss).
* Problem with confidence intervals in the predict method has been solved (reported on R-help 2010-11-28).
* Problem with `logDose` argument in `drm()` and subsequent plotting has been solved (reported by Ralf Schäfer).
* Constrained estimation now uses the actual limits provided (bug report by Andrew Kniss).
* Error in ED calculation for the logistic models (e.g., `L.4()`) has been fixed, caused by an update of `deltaMethod()` in alr3 (after a report from Daniel Gerhard).
* A small bug in `ED()` has been fixed (after a report from Nathan Pace).
* The functions `isobole()` and `mixture()` have been modified to fix an error in the plotting of the isoboles (after a report from Andreas Betz).
* The self starter function for llogistic models can now handle infinite dose values (after a report from Marc Weimer).

## Changes
* Plotting for event-time data has been improved.
* Help page of `drm()` (`type` argument) has been updated (after a comment from Radu Slobodeanu).
* Help page of `mixture()` has been slightly improved.
* The functions `genLoewe()`, `genLoewe2()`, `genBliss()`, `genBliss2()`, and `ursa()` have been improved with respect to starting values. The model specification has also been simplified.
* `genBliss()` and `genBliss2()` have been updated (after suggestions from Hugo Ceulemans).
* Small changes for the anova method; added df for event time and Poisson models.
* Small update in the calculation of confidence intervals (after a suggestion from Scott Ray).
* Small improvement in the calculation of ED values for the log-normal and Weibull type 2 models.
* The help page for `spinach` has been improved.
* Simplified printing in `print.summary.drc()` for model fit summary output.

## Breaking Changes
* `multdrc()` has been completely removed.
* The BIC method has been removed (after a suggestion from Prof. Brian Ripley).

---

# drc Changes in 2010

## New Features
* The model `genursa2()` for fitting the generalized URSA model has been included (after an idea by Hugo Ceulemans).
* The model `genursa()` for fitting the generalized URSA model has been included (after an idea by Hugo Ceulemans).
* The model functions `genBliss()`, `genBliss2()`, and `genLoewe2()` for fitting generalized Bliss independence and Loewe additivity with different maxima have been included (after an idea by Hugo Ceulemans).
* The model function `genLoewe()` for fitting generalized Loewe additivity has been included (after an idea by Hugo Ceulemans).
* `mselect()` has been extended to include an argument for specifying the type of information criterion to use.
* `maED()` has been extended to include simple linear regression.
* `mselect()` has been extended to include a few standard polynomial regression models.
* Studentised residuals are now also available for binomial responses (after an inquiry from Stuart Rosen).
* Studentised residuals are now available (after an inquiry from John).
* The argument `clevel` has been added to `ED()`. The function `maED()` has been extended to handle model fits involving several curves (after a suggestion from Andre Kleensang).
* The `comped()` function has been re-introduced (after a suggestion from Jochen Zubrod).
* `logLik` has been extended with a `nobs` attribute (provided by Tobias Verbeke). An S4 BIC method has been added (also provided by Tobias Verbeke).
* The model function `ursa()` for describing combination effects has been added (after an idea by Hugo Ceulemans).

## Bug Fixes
* A small bug in `ED()` (mismatch of curve names and parameter estimates) has been fixed (after a report from Andreas Betz).
* A small bug in the plot produced by `isobole()` has been fixed (after an inquiry from Andreas Betz).
* A bug in the calculation of standard errors in `SI()` has been fixed (after a report from Andrea Onofri).
* ED values are now correct for `cedergreen()` models (after a bug report by Claire Della Vedova).
* Small bug in `plot.drc()` related to `xt` and `xtlab` arguments has been fixed (after a comment by Anja Coors). A small bug in `comped()` has also been fixed (after a question from Jochen Zubrod).
* Minor bug in the predict method has been fixed, so `type="bars"` in the plot method now works again (after a bug report by Andy Robinson).
* The model function `cedergreen()` has been updated to ensure correct calculation of ED values (reported by Claire Della Vedova).

## Changes
* Small improvement in the summary output; the actual curve names are now used.
* The help for `NEC()` has been improved (after a question from Inés González).
* The help page of `drm()` has been improved with respect to the use of weights (after a question by Xuesong Yu).
* A slight modification in `vcov.drc()` to suppress unnecessary error messages (reported by Xuesong Yu).
* A slightly different bisection method has been implemented in `ursa()`. The corresponding help page has also been extended.
* The help page for `MM.2()` and `MM.3()` has been improved.

## Breaking Changes
* The model function `genursa2()` has been replaced by `actimL()`.
* Dataset `ecvam` has been removed. The `mixdrc()` function has been temporarily removed.

---

# drc Changes in 2009

## New Features
* The functions `NEC.2()`, `NEC.3()`, `NEC.4()` for estimation of no effect concentration have been included (after an idea by Ralf Schaefer).
* The argument `extended` has been added to the function `maED()`.
* The model function `twophase()` based on log-logistic models has been added (after an idea by Ida Katarina Auf der Maur Hindrichsen).
* The functions `lin.test()`, `mr.test()`, and `neill.test()` have been added.
* The data frame `etmotc` has been included.
* The argument `display` (with same functionality as in `ED()`) has been added for `compParm()` (after a suggestion from Scott Ray).
* The dataset `ecvam` has been included (migrated from the package 'mrdrc').
* The function `maED()` for parametric model averaging has been added.
* The functions `fplogistic()` and `FPL.4()` enable fitting dose-response models based on fractional polynomials.
* The function `getInitial()` has been included.
* `drm()` extended to allow fitting models separately for each curve. `confint`, `summary`, `vcov`, `ED`, `MAX`, `SI` have a new argument `pool` to allow pooling of separate fits.
* `modelFit()` now also works for binomial data.
* The dataset `algae` has been included.
* The argument `xsty` has been introduced to control the arrangement of tick marks on the dose axis.
* The function `yieldLoss()` has been included to handle a different parameterization of the Michaelis-Menten model (after a suggestion by Andrew Kniss).
* A new function `modelFit()` has been introduced for assessing the model fit, partly replacing the anova method.
* Argument `fixed` added to `BC.4()` and `BC.5()` (thanks to Nina Cedergreen).
* `mixture()` now also works for binomial data.

## Bug Fixes
* Small bug in `cedergreen()` and `ucedergreen()` related to calculation of ED values has been fixed (reported by Clare Della Vedova).
* Small bug in the function `modelFit()` has been fixed (after feedback from Heike Schmitt).
* Bug in `ED()` concerning `type="absolute"` and `reference="upper"` has been resolved (reported by Yue Zeng-Li).
* Error in likelihood calculation for some binomial models has been fixed.
* Bug in `mselect()` has been fixed (reported by John Lewis).
* Bug in `boxcox.drc` has been fixed.
* Bug in `predict.drc` fixed (thanks to Mario D'Antuono).

## Changes
* The help page of `drm()` has been updated.
* `confint` has been improved to automatically use the appropriate reference distribution for the confidence intervals; internal structure of `ED()` has also been modified (reported by Marc Weimer).
* The self starter for `twophase()` has been slightly improved.
* Help page for `CRS.5a` has been improved (after a comment from Claire Della Vedova).
* The help pages for `anova.drc` and `BC.4()`, `BC.5()`, `CRS.5.()` have been improved.
* Structure of self starter functions has been completely revamped. Four initial value procedures are now available for almost all implemented dose-response models.
* Help page for `ED()` has been improved (after a comment from Claire Della Vedova).
* Help page for `gompertz()` has been improved slightly.
* Help page for `earthworms` dataset has been improved.
* Minor internal changes in `drm()` and in the plot method.
* `compParm()`, `ED()`, `SI()` have been restructured.
* `vcov` method has been re-structured.
* The argument `conLevel` now has a more sensible default (no longer hardcoded at 0.01).
* `boxcox.drc` method extended to include functionality previously available through `drm()`.
* `drm()` has been improved with respect to handling extremely small or large dose or response values.
* `mixture()` has been completely revised with a lot of changes to the arguments.
* Redundant encoding removed.

## Breaking Changes
* The function `comped()` has been removed.
* The convenience functions `b.3()`, `B.3()`, `b.4()`, `B.4()`, `b.5()`, `B.5()`, and `boltzmann` have been removed. Use `L.3()`, `L.4()`, and `L.5()` instead.
* Argument `ci` in `relpot()` and `SI()` has been renamed to `interval`.
* Argument `ci` in `ED()` renamed to `interval`.
* The function `plotraw()` has been removed. Use R's standard plotting functionality instead.
* The argument `fctList` has been removed from `drm()`.
* The model function `richards()` has been removed as it is a different parameterization of the five-parameter log-logistic model. `colFct()` has also been removed.
* `multdrc()` and associated `mdControl` have been taken completely out of use.
* Arguments `lowerc` and `upperc` have been removed in model functions as they were redundant.
* Argument `legendCex` in the plot method renamed to `cex.legend` in line with other cex arguments.

---

# drc Changes in 2008

## New Features
* Dataset `H.virescens` added.
* `lnormal` function has been added plus three new datasets.
* `gompertz` function has been added.
* Datasets `lepidium` and `nasturtium` have been added.
* New function `mrdrm()` for model-robust modelling included, with accompanying `ED` and `predict` methods.

## Bug Fixes
* Error in `level` argument in plot method has been fixed.
* Bug in `level` argument in plot method has been fixed.
* The `lettuce` dataset got correct row numbers.

## Changes
* Asymptotic regression and exponential decay implemented differently.

# Changelog

## drc 3.3.0.02

### Changes

- Added `NEWS.md` version control log.

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
- Added robust estimation methods in new `ED_robust.R` module:
  [`ED_robust()`](https://hreinwald.github.io/drc/reference/ED_robust.md)
  for calculating ED values using robust median-based estimation,
  [`maED_robust()`](https://hreinwald.github.io/drc/reference/maED_robust.md)
  for model-averaged robust ED estimation, and
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
  [`devtools::check()`](https://rdrr.io/pkg/devtools/man/check.html)
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
- Fixed NaN warning in `summary.drc` for robust median estimation.
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

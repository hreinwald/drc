# Release Notes: drc 3.3.0.02

**Release Date:** March 9, 2026

## Overview

Version 3.3.0.02 is a maintenance release focused on bug fixes and code quality improvements. This release addresses 17 critical issues in the `ucedergreen()` function, improves error handling across multiple core functions, and adds comprehensive test coverage.

## Key Highlights

### Critical Bug Fixes
- **Fixed 17 issues in `ucedergreen()` function** - Resolved model formula errors, signature mismatches, unsafe parameter handling, and broken self-starter functionality
- **Enhanced SE calculation in `ED()`** - Corrected standard error estimation for absolute type effective dose calculations
- **Improved boundary detection in `MAX()`** - Fixed numerical comparison issues with named return values and boundary tolerance

### Improvements
- **Added comprehensive test suites** for 12+ core functions including `summary.drc`, `searchdrc`, `backfit`, `MAX()`, and `PR()`
- **Enhanced documentation** for Weibull starting value methods across all wrapper functions
- **Code cleanup** - Removed dead code (`iband.R`), unused files, and build scripts

### Additional Fixes
- Fixed `PR()` argument handling for single-curve models
- Corrected inverted trace/silent logic in `drmOpt()`
- Fixed regex error and convergence behavior in `searchdrc()`
- Updated citation URL configuration

## Installation

Install the latest development version from GitHub:

```r
devtools::install_github("hreinwald/drc")
```

## Full Changelog

For detailed information about all changes, see [NEWS.md](NEWS.md).

## Credits

Maintained by Christian Ritz, Jens C. Streibig, and Hannes Reinwald.

For bug reports and feature requests, visit: https://github.com/hreinwald/drc/issues

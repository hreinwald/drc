[![GitHub dev version](https://img.shields.io/github/r-package/v/hreinwald/drc)](https://github.com/hreinwald/drc)
[![Documentation](https://img.shields.io/static/v1?style=flat-square&message=ReadTheDocs&color=2C4AA8&logo=ReadTheDocs&logoColor=FFFFFF&label=Documentation)](https://hreinwald.github.io/drc/)
<a href="https://github.com/hreinwald/drc/actions/workflows/R-CMD-check.yaml"><img src="https://github.com/hreinwald/drc/actions/workflows/R-CMD-check.yaml/badge.svg" alt="R-CMD-check status"></a>
<a href="https://app.codecov.io/gh/hreinwald/drc"><img src="https://codecov.io/gh/hreinwald/drc/branch/dev/graph/badge.svg" alt="Code coverage status"></a>
<a href="https://lifecycle.r-lib.org/articles/stages.html"><img src="https://img.shields.io/badge/lifecycle-stable-green.svg" alt="Lifecycle: stable"></a>

<a href="https://cran.r-project.org/package=drc"><img src="https://www.r-pkg.org/badges/version/drc" alt="CRAN version"></a>
[![Downloads](https://cranlogs.r-pkg.org/badges/drc)](https://cranlogs.r-pkg.org/)
<a href="https://github.com/hreinwald/drc/blob/dev/LICENSE"><img src="https://img.shields.io/github/license/hreinwald/drc" alt="License: GPL-2.0"></a>
<a href="https://github.com/hreinwald/drc/commits/dev"><img src="https://img.shields.io/github/last-commit/hreinwald/drc" alt="Last commit date"></a>
<a href="https://github.com/hreinwald/drc/issues"><img src="https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat" alt="Contributions welcome"></a>

<p align="center">
  <img src="man/figures/logo.png" alt="drc Logo" width="250">
</p>


# drc — Dose-Response Curve Analysis in R

## Note

This repository contains a refactored development version of the [*drc*](https://github.com/DoseResponse/drc) R package first published by **Christian Ritz, Florent Baty, Jens C. Streibig und Daniel Gerhard** [(2015)](https://doi.org/10.1371/journal.pone.0146021). Their foundational work on dose–response modeling in R is gratefully acknowledged and inspired the present refactoring.

The goal of this project is to modernize the codebase, improve maintainability, and provide a clearer development structure while preserving the core functionality of the original package.

This repository focuses on structural refactoring and development improvements. Behavior and interfaces may change as the codebase is modernized.

## Overview

The **drc** package provides a comprehensive framework for fitting, analyzing, and visualizing dose-response curves in R. It is widely used in bioassay, toxicology, pharmacology, and agricultural research to model the relationship between an exposure (e.g., concentration of a substance) or dose and a biological response.

The package offers:

- **Flexible model fitting** via the central `drm()` function, supporting multiple data types (continuous, binomial, Poisson, negative binomial, event-time, and species sensitivity distributions).
- **40+ built-in parametric models** including log-logistic, Weibull, Gompertz, Brain-Cousens, Cedergreen, and many more, each with self-starting parameter initialization.
- **Effective dose (ED) estimation** with confidence intervals (delta method, Fieller, inverse regression) through `ED()`.
- **Model comparison and diagnostics**: ANOVA, lack-of-fit tests, Neill's test, Box-Cox transformations, R-squared, Cook's distance, and hat values.
- **Multi-curve analysis**: fit and compare dose-response curves across groups, compute relative potency and selectivity indices via `EDcomp()`.
- **Robust inference**: sandwich variance estimators for heteroscedasticity-consistent standard errors.
- **Simulation tools**: generate random dose-response data for power analysis and method comparison.

For more details visit:

:book: **[drc github documentation](https://hreinwald.github.io/drc/)**  
:zap: **[drc example workflow](https://hreinwald.github.io/drc/articles/dose-response-workflow.html)**

Feature requests or ideas?

:bulb: **[Post them here](https://github.com/hreinwald/drc/discussions)**

## Installation

**⚠️ Important:** We **do not recommend** installing the currently heavily outdated CRAN version of this package. Instead, we recommend installing the development (`dev`) or stable beta (`main_beta`) version from GitHub.

### Install from GitHub (Recommended)

``` r
# install.packages("devtools")

# Install the re-factored development version 
devtools::install_github("hreinwald/drc")

# Install the re-factored stable beta version
devtools::install_github("hreinwald/drc@main_beta")
```

### Local Installation from tar.gz

If GitHub installation is failing, you can run the installation from the local tar.gz file. 
[Download the latest release](https://github.com/hreinwald/drc/archive/refs/tags/3.3.1.tar.gz).

After downloading the file, run the following:

``` r
# Specify the path to the directory where you saved the downloaded tar.gz file.
# Make sure to specify the correct file path below.
targz  <- file.path("~/Downloads/drc-3.3.1.tar.gz") # adjust file path here accordingly

# Local installation with base R
install.packages(targz, repos = NULL, type = "source")
```

### Outdated CRAN Version (Not Recommended)

To install the outdated version from CRAN:

``` r
install.packages("drc")
```

## Quick Start

### Fitting a basic dose-response model

``` r
library(drc)

# Fit a four-parameter log-logistic model to the built-in 'ryegrass' dataset
model <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

# View model summary with parameter estimates and standard errors
summary(model)

# Plot the fitted dose-response curve
plot(model, xlab = "Concentration", ylab = "Root length")
```

### Estimating effective doses (ED values)

``` r
# Estimate the ED50 (dose producing 50% effect) with confidence intervals
ED(model, respLev = c(10, 50, 90), interval = "delta")
```

### Comparing curves across groups

``` r
# Fit separate curves for multiple groups
model_multi <- drm(rootl ~ conc, curveid = herbicide,
                   data = ryegrass, fct = LL.4())

# Compare ED50 values between groups
EDcomp(model_multi, percVec = c(50), interval = "delta")
```

### Model selection

``` r
# Compare different dose-response model families
mselect(model, fctList = list(W1.4(), W2.4(), LL.3()))
```

## Vignettes

The package includes detailed vignettes to help you understand specific topics:

``` r
# View available vignettes
vignette(package = "drc")

# Access the NEC models vignette
vignette("nec-models", package = "drc")
```

## Available Models

| Function  | Description                                     |
|-----------|-------------------------------------------------|
| `LL.2()` – `LL.5()` | Log-logistic models (2 to 5 parameters)  |
| `W1.2()` – `W1.4()` | Weibull type 1 models                    |
| `W2.2()` – `W2.4()` | Weibull type 2 models                    |
| `G.3()`, `G.4()`    | Gompertz models                          |
| `LN.2()` – `LN.4()` | Log-normal models                       |
| `BC.4()`, `BC.5()`  | Brain-Cousens models (hormesis)          |
| `CRS.4a()` – `CRS.4c()` | Cedergreen-Ritz-Streibig 4-parameter models (hormesis) |
| `CRS.5()`, `CRS.5a()` – `CRS.5c()` | Cedergreen-Ritz-Streibig 5-parameter models (hormesis) |
| `CRS.6()` | Generalised Cedergreen-Ritz-Streibig model (hormesis) |
| `UCRS.4a()` – `UCRS.4c()` | U-shaped Cedergreen-Ritz-Streibig 4-parameter models (hormesis) |
| `UCRS.5a()` – `UCRS.5c()` | U-shaped Cedergreen-Ritz-Streibig 5-parameter models (hormesis) |
| `NEC.2()` – `NEC.4()` | No-effect-concentration models         |
| `L.3()` – `L.5()`   | Logistic models                          |
| `baro5()`            | Baro five-parameter model                |
| `gammadr()` | Gamma dose-response model                        |

## Key Functions

| Function      | Purpose                                            |
|---------------|----------------------------------------------------|
| `drm()`       | Fit dose-response models                           |
| `ED()`        | Estimate effective doses (ED10, ED50, ...)  |
| `maED()`      | Model averaged estimate effective doses (ED10, ED50, ...)  |
| `EDcomp()`    | Compare ED values between curves                   |
| `compParm()`  | Compare model parameters between curves            |
| `noEffect()`  | Testing if there is a dose effect at all           |
| `plot()`      | Plot fitted dose-response curves                   |
| `summary()`   | Model summary with parameter estimates             |
| `anova()`     | ANOVA and lack-of-fit tests                        |
| `mselect()`   | Model selection among candidate models             |
| `predict()`   | Predictions with confidence/prediction intervals   |
| `modelFit()`  | Goodness-of-fit test                               |
| `Rsq()`       | R-squared calculation                              |
| `rdrm()`      | Simulate dose-response data                        |

## Data Types Supported

The `drm()` function supports multiple response types via the `type` argument:

- **`"continuous"`** (default): Standard continuous dose-response data.
- **`"binomial"`**: Quantal/binary response data (e.g., proportion of individuals affected).
- **`"Poisson"`**: Count data following a Poisson distribution.
- **`"negbin1"`, `"negbin2"`**: Negative binomial count data.
- **`"event"`**: Event-time / time-to-event data (e.g., germination time).
- **`"ssd"`**: Species sensitivity distributions for ecotoxicology.

## Dependencies

**drc** depends on:
- R (≥ 4.0.0), MASS, stats

and imports from: car, graphics, gtools, lifecycle, multcomp, plotrix, sandwich, scales, utils.

## References

- Ritz, C., Baty, F., Streibig, J. C., and Gerhard, D. (2015). Dose-Response Analysis Using R. *PLOS ONE*, 10(12), e0146021.
- Ritz, C. and Streibig, J. C. (2005). Bioassay Analysis using R. *Journal of Statistical Software*, 12(5), 1–22.

## Bug Reports

Please report issues with this re-factory version [here](https://github.com/hreinwald/drc/issues/).

## License

GPL-2.0


[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/drc)](https://cran.r-project.org/package=drc) [![Build Status](https://travis-ci.org/DoseResponse/drc.svg?branch=master)](https://travis-ci.org/DoseResponse/drc) [![Downloads](https://cranlogs.r-pkg.org/badges/drc)](https://cranlogs.r-pkg.org/)

# drc — Dose-Response Curve Analysis in R

## Overview

The **drc** package provides a comprehensive framework for fitting, analyzing, and visualizing dose-response curves in R. It is widely used in bioassay, toxicology, pharmacology, and agricultural research to model the relationship between a dose (e.g., concentration of a substance) and a biological response.

The package offers:

- **Flexible model fitting** via the central `drm()` function, supporting multiple data types (continuous, binomial, Poisson, negative binomial, event-time, and species sensitivity distributions).
- **40+ built-in parametric models** including log-logistic, Weibull, Gompertz, Brain-Cousens, Cedergreen, and many more, each with self-starting parameter initialization.
- **Effective dose (ED) estimation** with confidence intervals (delta method, Fieller, inverse regression) through `ED()`.
- **Model comparison and diagnostics**: ANOVA, lack-of-fit tests, Neill's test, Box-Cox transformations, R-squared, Cook's distance, and hat values.
- **Multi-curve analysis**: fit and compare dose-response curves across groups, compute relative potency and selectivity indices via `EDcomp()`.
- **Robust inference**: sandwich variance estimators for heteroscedasticity-consistent standard errors.
- **Simulation tools**: generate random dose-response data for power analysis and method comparison.

## Installation

To install the stable but outdated version from CRAN:

``` r
install.packages("drc")
```

Install the bug-fix development version from GitHub:

``` r
# install.packages("devtools")

# Install the companion data package first
devtools::install_github("DoseResponse/drcData")

# Then install the development version of drc
# devtools::install_github("DoseResponse/drc")      # <- unmaintained version 3.2-0 
devtools::install_github("hreinwald/drc")           # <- re-factored development version 3.3.0
devtools::install_github("hreinwald/drc@main_beta") # <- re-factored stable beta version 3.3.0
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
| `ED()`        | Estimate effective doses (ED10, ED50, ED90, etc.)  |
| `EDcomp()`    | Compare ED values between curves                   |
| `compParm()`  | Compare model parameters between curves            |
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
- R (≥ 4.0.0), MASS, stats, [drcData](https://github.com/DoseResponse/drcData)

and imports from: car, graphics, gtools, lifecycle, multcomp, plotrix, sandwich, scales, utils.

## References

- Ritz, C., Baty, F., Streibig, J. C., and Gerhard, D. (2015). Dose-Response Analysis Using R. *PLOS ONE*, 10(12), e0146021.
- Ritz, C. and Streibig, J. C. (2005). Bioassay Analysis using R. *Journal of Statistical Software*, 12(5), 1–22.

## Bug Reports

Please report issues with this re-factory version at <https://github.com/hreinwald/drc/issues/>.

## License

GPL-2

# Dose-response profile of degradation of agrochemical using nasturtium

Estimation of the degradation profile of an agrochemical based on soil
samples at depth 0-10cm from a calibration experiment.

## Usage

``` r
data(nasturtium)
```

## Format

A data frame with 42 observations on the following 2 variables.

- `conc`:

  a numeric vector of concentrations (g/ha)

- `wt`:

  a numeric vector of plant weight (mg) after 3 weeks' growth

- `rep`:

  a numeric vector of replicates

## Details

It is an experiment with seven concentrations and six replicates per
concentration. *Nasturtium* is sensitive and its weight reduces
noticeable at low concentrations.

Racine-Poon (1988) suggests using a three-parameter log-logistic model.

## Source

Racine-Poon, A. (1988) A Bayesian Approach to Nonlinear Calibration
Problems, *J. Am. Statist. Ass.*, **83**, 650–656.

## Examples

``` r
library(drc)

nasturtium.m1 <- drm(wt~conc, data=nasturtium, fct = LL.3())

modelFit(nasturtium.m1)
#> Lack-of-fit test
#> 
#>           ModelDf    RSS Df F value p value
#> ANOVA          35 104464                   
#> DRC model      39 120387  4  1.3337  0.2768

plot(nasturtium.m1, type = "all", log = "", xlab = "Concentration (g/ha)", ylab = "Weight (mg)")
```

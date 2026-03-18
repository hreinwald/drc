# Dose-response profile of degradation of agrochemical using lepidium

Estimation of the degradation profile of an agrochemical based on soil
samples at depth 0-10cm from a calibration experiment.

## Usage

``` r
data(lepidium)
```

## Format

A data frame with 42 observations on the following 2 variables.

- `conc`:

  a numeric vector of concentrations (g/ha)

- `weight`:

  a numeric vector of plant weight (g) after 3 weeks' growth

## Details

It is an experiment with seven concentrations and six replicates per
concentration. *Lepidium* is rather robust as it only responds to high
concentrations.

## Source

Racine-Poon, A. (1988) A Bayesian Approach to Nonlinear Calibration
Problems, *J. Am. Statist. Ass.*, **83**, 650–656.

## Examples

``` r
library(drc)

lepidium.m1 <- drm(weight~conc, data=lepidium, fct = LL.4())

modelFit(lepidium.m1)
#> Lack-of-fit test
#> 
#>           ModelDf    RSS Df F value p value
#> ANOVA          35 14.187                   
#> DRC model      38 14.449  3  0.2159  0.8847

plot(lepidium.m1, type = "all", log = "")
```

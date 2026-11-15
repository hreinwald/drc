# ANOVA Model Comparison for Dose-Response Models

Compares two nested dose-response model fits using a likelihood-ratio
test (for binomial data) or an F-test (for continuous data). Two `drc`
objects must be provided. For a lack-of-fit test of a single model, use
[`modelFit`](https://hreinwald.github.io/drc/reference/modelFit.md)
instead.

## Usage

``` r
# S3 method for class 'drc'
anova(object, ..., details = TRUE, test = NULL)
```

## Arguments

- object:

  an object of class ‘drc’.

- ...:

  a second object of class ‘drc’ to compare against `object`. Exactly
  two models must be supplied; passing a single model will result in an
  error directing the user to
  [`modelFit`](https://hreinwald.github.io/drc/reference/modelFit.md).

- details:

  logical indicating whether or not details on the models compared
  should be displayed. Default is `TRUE` (details are displayed).

- test:

  a character string specifying the test statistic to be applied. For
  continuous data the default is `"F"` (F-test); for binomial data the
  default is `"Chisq"` (likelihood-ratio test). Use `"Chisq"` to force a
  likelihood-ratio test for continuous data.

## Value

An object of class ‘anova’ (inheriting from `data.frame`) with columns
for model degrees of freedom, residual sum of squares (or
log-likelihood), the difference in degrees of freedom, the test
statistic, and the p-value.

## Details

Two `drc` objects must be specified. The function performs a test for
reduction from the larger to the smaller model. This only makes
statistical sense if the models are nested, that is: one model is a
submodel of the other model.

For continuous data an F-test is used by default. For binomial data a
likelihood-ratio (chi-square) test is used by default.

If a single model is passed, the function raises an error. To assess the
fit of a single dose-response model (lack-of-fit test comparing the
model to a more general ANOVA model), use
[`modelFit`](https://hreinwald.github.io/drc/reference/modelFit.md)
instead.

## See also

[`modelFit`](https://hreinwald.github.io/drc/reference/modelFit.md) for
lack-of-fit testing of a single model,
[`drm`](https://hreinwald.github.io/drc/reference/drm.md) for fitting
dose-response models,
[`logLik.drc`](https://hreinwald.github.io/drc/reference/logLik.drc.md)
for log-likelihood extraction,
[`summary.drc`](https://hreinwald.github.io/drc/reference/summary.drc.md)
for model summaries.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
## Comparing two nested models (two-model comparison)
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
ryegrass.m2 <- drm(rootl ~ conc, data = ryegrass, fct = W1.3())
anova(ryegrass.m2, ryegrass.m1)
#> 
#> 1st model
#>  fct:      W1.3()
#> 2nd model
#>  fct:      W1.4()
#> 
#> ANOVA table
#> 
#>           ModelDf    RSS Df F value p value
#> 1st model      21 8.9520                   
#> 2nd model      20 6.0242  1  9.7205  0.0054

anova(ryegrass.m2, ryegrass.m1, details = FALSE)  # without details
#> ANOVA table
#> 
#>           ModelDf    RSS Df F value p value
#> 1st model      21 8.9520                   
#> 2nd model      20 6.0242  1  9.7205  0.0054

## For a lack-of-fit test on a single model, use modelFit():
modelFit(ryegrass.m1)
#> Lack-of-fit test
#> 
#>           ModelDf    RSS Df F value p value
#> ANOVA          17 5.1799                   
#> DRC model      20 6.0242  3  0.9236  0.4506
```

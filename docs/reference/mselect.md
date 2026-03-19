# Dose-response model selection

Model selection by comparison of different models using the following
criteria: the log likelihood value, Akaike's information criterion
(AIC), the estimated residual standard error or the p-value from a
lack-of-fit test.

## Usage

``` r
mselect(
  object,
  fctList = NULL,
  nested = FALSE,
  sorted = c("IC", "Res var", "Lack of fit", "no"),
  linreg = FALSE,
  icfct = AIC
)
```

## Arguments

- object:

  an object of class 'drc'.

- fctList:

  a list of dose-response functions to be compared.

- nested:

  logical. TRUE results in F tests between adjacent models (in
  `fctList`). Only sensible for nested models.

- sorted:

  character string determining according to which criterion the model
  fits are ranked.

- linreg:

  logical indicating whether or not additionally polynomial regression
  models (linear, quadratic, and cubic models) should be fitted.

- icfct:

  function for supplying the information criterion to be used.
  [`AIC`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC`](https://rdrr.io/r/stats/AIC.html) are two options.

## Value

A matrix with one row for each model and one column for each criterion.

## Details

For Akaike's information criterion and the residual standard error: the
smaller the better and for lack-of-fit test (against a one-way ANOVA
model): the larger (the p-value) the better. Note that the residual
standard error is only available for continuous dose-response data.

Log likelihood values cannot be used for comparison unless the models
are nested.

## Author

Christian Ritz

## Examples

``` r
### Example with continuous/quantitative data
## Fitting initial four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## Model selection
mselect(ryegrass.m1, list(LL.3(), LL.5(), W1.3(), W1.4(), W2.4(), baro5()))
#>          logLik       IC Lack of fit   Res var
#> W2.4  -15.91352 41.82703           0 0.2646283
#> LL.4  -16.15514 42.31029           0 0.2700107
#> baro5 -15.86422 43.72844           0 0.2774141
#> LL.5  -15.87828 43.75656           0 0.2777393
#> W1.4  -17.46720 44.93439           0 0.3012075
#> LL.3  -18.60413 45.20827           0 0.3153724
#> W1.3  -22.22047 52.44094           0 0.4262881

## Model selection including linear, quadratic, and cubic regression models
mselect(ryegrass.m1, list(LL.3(), LL.5(), W1.3(), W1.4(), W2.4(), baro5()), linreg = TRUE)
#>          logLik        IC Lack of fit   Res var
#> W2.4  -15.91352  41.82703           0 0.2646283
#> LL.4  -16.15514  42.31029           0 0.2700107
#> baro5 -15.86422  43.72844           0 0.2774141
#> LL.5  -15.87828  43.75656           0 0.2777393
#> W1.4  -17.46720  44.93439           0 0.3012075
#> LL.3  -18.60413  45.20827           0 0.3153724
#> W1.3  -22.22047  52.44094           0 0.4262881
#> Cubic -25.53428  61.06856          NA 0.5899609
#> Quad  -35.11558  78.23116          NA 1.2485122
#> Lin   -50.47554 106.95109          NA 4.2863247

## Comparing nested models
mselect(ryegrass.m1, list(LL.5()), nested = TRUE)
#>         logLik       IC Lack of fit   Res var Nested F test
#> LL.4 -16.15514 42.31029   0.8664830 0.2700107            NA
#> LL.5 -15.87828 43.75656   0.8538476 0.2777393     0.5134602
```

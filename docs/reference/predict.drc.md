# Prediction

Predicted values for models of class 'drc'.

## Usage

``` r
# S3 method for class 'drc'
predict(
  object,
  newdata,
  se.fit = FALSE,
  interval = c("none", "confidence", "prediction", "ssd"),
  level = 0.95,
  na.action = na.pass,
  od = FALSE,
  vcov. = vcov,
  ssdSEfct = NULL,
  constrain = TRUE,
  checkND = TRUE,
  ...
)
```

## Arguments

- object:

  an object of class 'drc'.

- newdata:

  an optional data frame in which to look for variables with which to
  predict. If omitted, the fitted values are used.

- se.fit:

  logical. If TRUE standard errors are required.

- interval:

  character string. Type of interval calculation: `"none"`,
  `"confidence"`, `"prediction"`, or `"ssd"`.

- level:

  tolerance/confidence level.

- na.action:

  function determining what should be done with missing values in
  `newdata`. The default is to predict `NA`.

- od:

  logical. If TRUE adjustment for over-dispersion is used.

- vcov.:

  function providing the variance-covariance matrix.
  [`vcov`](https://rdrr.io/r/stats/vcov.html) is the default, but
  `sandwich` is also an option (for obtaining robust standard errors).

- ssdSEfct:

  specifies the function for interpolating standard errors between
  observed standard errors. The default is linear interpolation on
  log-log scale (back-transformed).

- constrain:

  logical. If TRUE (default) predicted values are truncated within
  meaningful limits, i.e., 0 and, possibly, 1.

- checkND:

  logical indicating whether or not names in `newdata` data frame match
  the names in the original data frame used for fitting the model.
  Default is TRUE.

- ...:

  further arguments passed to or from other methods.

## Value

A matrix with as many rows as there are dose values provided in
`newdata` or in the original dataset (in case `newdata` is not
specified) and, at most, 4 columns containing fitted values, standard
errors, lower and upper limits of confidence/prediction intervals. If
only fitted values are requested (`se.fit = FALSE` and
`interval = "none"`) a numeric vector is returned.

## Details

Standard errors and intervals use the delta method with the analytic
first derivatives of the mean function (`deriv1`) when the model
provides them. For models without analytic derivatives (e.g.
[`NEC.4`](https://hreinwald.github.io/drc/reference/NEC.4.md),
[`CRS.6`](https://hreinwald.github.io/drc/reference/CRS.6.md),
[`baro5`](https://hreinwald.github.io/drc/reference/baro5.md)) a
numerical gradient is used instead. If neither `se.fit` nor `interval`
is requested only the predicted values are computed.

## See also

For details see the help page for
[`predict.lm`](https://rdrr.io/r/stats/predict.lm.html).

## Author

Christian Ritz

## Examples

``` r
## Fitting a model
spinach.model1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())

## Predicting values at dose=2 (with standard errors)
predict(spinach.model1, data.frame(dose=2, CURVE=c("1", "2", "3")), se.fit = TRUE)
#>      Prediction         SE
#> [1,]  0.9048494 0.02496159
#> [2,]  0.4208338 0.02924978
#> [3,]  0.5581703 0.03067242

## Getting confidence intervals
predict(spinach.model1, data.frame(dose=2, CURVE=c("1", "2", "3")),
interval = "confidence")
#>      Prediction     Lower     Upper
#> [1,]  0.9048494 0.8552190 0.9544797
#> [2,]  0.4208338 0.3626774 0.4789903
#> [3,]  0.5581703 0.4971853 0.6191553

## Getting prediction intervals
predict(spinach.model1, data.frame(dose=2, CURVE=c("1", "2", "3")),
interval = "prediction")
#>      Prediction     Lower     Upper
#> [1,]  0.9048494 0.7504606 1.0592382
#> [2,]  0.4208338 0.2634969 0.5781707
#> [3,]  0.5581703 0.3997661 0.7165745
```

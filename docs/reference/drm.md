# Fitting dose-response models

A general model fitting function for analysis of various types of
dose-response data.

## Usage

``` r
drm(
  formula,
  curveid,
  pmodels,
  weights,
  data = NULL,
  subset,
  fct,
  type = c("continuous", "binomial", "Poisson", "negbin1", "negbin2", "event", "ssd"),
  bcVal = NULL,
  bcAdd = 0,
  start,
  na.action = na.omit,
  robust = "mean",
  logDose = NULL,
  control = drmc(),
  lowerl = NULL,
  upperl = NULL,
  separate = FALSE,
  pshifts = NULL,
  varcov = NULL
)
```

## Arguments

- formula:

  a symbolic description of the model to be fit. Either of the form
  `response ~ dose` or as a data frame with response values in first
  column and dose values in second column.

- curveid:

  a numeric vector or factor containing the grouping of the data.

- pmodels:

  a data frame with as many columns as there are parameters in the
  non-linear function. Or a list containing a formula for each parameter
  in the nonlinear function.

- weights:

  a numeric vector containing weights. For continuous/quantitative
  responses, inverse weights are multiplied inside the squared errors
  (weights should have the same unit as the response). For binomial
  responses weights provide information about the total number of binary
  observations used to obtain the response.

- data:

  an optional data frame containing the variables in the model.

- subset:

  an optional vector specifying a subset of observations to be used in
  the fitting process.

- fct:

  a list with three or more elements specifying the non-linear function,
  the accompanying self starter function, the names of the parameters in
  the non-linear function and, optionally, the first and second
  derivatives as well as information used for calculation of ED values.
  Use
  [`getMeanFunctions`](https://hreinwald.github.io/drc/reference/getMeanFunctions.md)
  for a full list.

- type:

  a character string specifying the distribution of the data. The
  default is `"continuous"`, corresponding to a normal distribution.
  Other choices include `"binomial"`, `"Poisson"`, `"negbin1"`,
  `"negbin2"`, `"event"`, and `"ssd"`.

- bcVal:

  a numeric value specifying the lambda parameter to be used in the
  Box-Cox transformation.

- bcAdd:

  a numeric value specifying the constant to be added on both sides
  prior to Box-Cox transformation. The default is 0.

- start:

  an optional numeric vector containing starting values for all mean
  parameters in the model. Overrules any self starter function.

- na.action:

  a function for treating missing values (`NA`s). Default is
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html).

- robust:

  a character string specifying the rho function for robust estimation.
  Default is non-robust least squares estimation (`"mean"`). Available
  robust methods are: `"median"`, `"lms"`, `"lts"`, `"trimmed"`,
  `"winsor"`, and `"tukey"`.

- logDose:

  a numeric value or `NULL`. If log dose values are provided the base of
  the logarithm should be specified (e.g., `exp(1)` for natural
  logarithm, `10` for base 10).

- control:

  a list of arguments controlling constrained optimisation, maximum
  iterations, relative tolerance, and warnings. See
  [`drmc`](https://hreinwald.github.io/drc/reference/drmc.md).

- lowerl:

  a numeric vector of lower limits for all parameters in the model (the
  default corresponds to minus infinity for all parameters).

- upperl:

  a numeric vector of upper limits for all parameters in the model (the
  default corresponds to plus infinity for all parameters).

- separate:

  logical value indicating whether curves should be fit separately
  (independent of each other).

- pshifts:

  a matrix of constants to be added to the matrix of parameters. Default
  is no shift for all parameters.

- varcov:

  an optional user-defined known variance-covariance matrix for the
  responses. Default is the identity matrix (`NULL`), corresponding to
  independent response values with a common standard deviation estimated
  from the data.

## Value

An object of (S3) class `"drc"`.

## Details

This function relies on [`optim`](https://rdrr.io/r/stats/optim.html)
for minimisation of the negative log likelihood function. For a
continuous response this reduces to least squares estimation. Response
values are assumed to be mutually independent unless `varcov` is
specified. For robust estimation MAD (median absolute deviance) is used
to estimate the residual variance. Setting `lowerl` and/or `upperl`
automatically invokes constrained optimisation. Control arguments may be
specified using
[`drmc`](https://hreinwald.github.io/drc/reference/drmc.md).

## See also

[`drmc`](https://hreinwald.github.io/drc/reference/drmc.md),
[`LL.4`](https://hreinwald.github.io/drc/reference/LL.4.md),
[`getMeanFunctions`](https://hreinwald.github.io/drc/reference/getMeanFunctions.md)

## Author

Christian Ritz, Jens C. Streibig and Hannes Reinwald

## Examples

``` r
## Fitting a four-parameter log-logistic model to the ryegrass data
model <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
summary(model)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept)  2.98222    0.46506  6.4125 2.960e-06 ***
#> c:(Intercept)  0.48141    0.21219  2.2688   0.03451 *  
#> d:(Intercept)  7.79296    0.18857 41.3272 < 2.2e-16 ***
#> e:(Intercept)  3.05795    0.18573 16.4644 4.268e-13 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.5196256 (20 degrees of freedom)
```

# Legacy dose-response model fitting (internal)

This is the legacy implementation of the dose-response model fitting
function. It is kept only as an internal reference point in case
questions or errors might occur with the current
[`drm()`](https://hreinwald.github.io/drc/reference/drm.md)
implementation.

## Usage

``` r
drm_legacy(
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

## See also

[`drm()`](https://hreinwald.github.io/drc/reference/drm.md) for the
current implementation.

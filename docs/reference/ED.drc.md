# Estimating effective doses

`ED` estimates effective concentration or doses for one or more
specified absolute or relative response levels.

## Usage

``` r
# S3 method for class 'drc'
ED(
  object,
  respLev = c(10, 20, 50),
  interval = c("none", "delta", "fls", "tfls", "inv"),
  clevel = NULL,
  level = 0.95,
  reference = c("control", "upper"),
  type = c("relative", "absolute"),
  lref,
  uref,
  bound = TRUE,
  vcov. = vcov,
  display = TRUE,
  logBase = NULL,
  multcomp = FALSE,
  intType = "confidence",
  ...
)
```

## Arguments

- object:

  an object of class `drc`.

- respLev:

  a numeric vector containing the response levels.

- interval:

  character string specifying the type of confidence intervals to be
  supplied. The default is `"none"`. See Details below for more
  explanation.

- clevel:

  character string specifying the curve id in case estimates for a
  specific curve or compound are requested. By default estimates are
  shown for all curves.

- level:

  numeric. The level for the confidence intervals. Must be a single
  value strictly between 0 and 1. The default is `0.95`.

- reference:

  character string. Is the upper limit or the control level the
  reference?

- type:

  character string. Whether the specified response levels are absolute
  or relative (default).

- lref:

  numeric value specifying the lower limit to serve as reference.

- uref:

  numeric value specifying the upper limit to serve as reference (e.g.,
  100%).

- bound:

  logical. Default is `TRUE`, in which case only ED values between 0 and
  100% are allowed. Set to `FALSE` for hormesis models.

- vcov.:

  function providing the variance-covariance matrix, or a
  variance-covariance matrix directly.
  [`vcov`](https://rdrr.io/r/stats/vcov.html) is the default, but
  `sandwich` is also an option for obtaining robust standard errors.

- display:

  logical. If `TRUE` results are displayed. Otherwise they are not
  (useful in simulations).

- logBase:

  numeric. The base of the logarithm in case logarithm transformed dose
  values are used.

- multcomp:

  logical to switch on output for use with the package multcomp (which
  needs to be activated first). Default is `FALSE`.

- intType:

  string specifying the type of interval to use with the predict method
  in case the type of confidence interval chosen is inverse regression.

- ...:

  additional arguments passed to the ED function in the model.

## Value

An invisible matrix containing the estimates and the corresponding
estimated standard errors and possibly lower and upper confidence
limits. Or, alternatively, a list with elements that may be plugged
directly into `parm` in the package multcomp (when `multcomp = TRUE`).

## Details

There are several options for calculating confidence intervals through
the argument `interval`. The option `"delta"` results in asymptotical
Wald-type confidence intervals (using the delta method and the normal or
t-distribution depending on the type of response). The option `"fls"`
produces (possibly skewed) confidence intervals through
back-transformation from the logarithm scale (only meaningful in case
the parameter in the model is log(ED50) as for the
[`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md)
models). The option `"tfls"` is for transforming back and forth from log
scale (experimental). The option `"inv"` results in confidence intervals
obtained through inverse regression.

For hormesis models
([`braincousens`](https://hreinwald.github.io/drc/reference/braincousens.md)
and
[`cedergreen`](https://hreinwald.github.io/drc/reference/cedergreen.md)),
the additional arguments `lower` and `upper` may be supplied. These
arguments specify the lower and upper limits of the bisection method
used to find the ED values.

## See also

[`EDcomp`](https://hreinwald.github.io/drc/reference/EDcomp.md) for
estimating differences and ratios of ED values,
[`compParm`](https://hreinwald.github.io/drc/reference/compParm.md) for
comparing other model parameters, and
[`backfit`](https://hreinwald.github.io/drc/reference/backfit.md).

## Author

Christian Ritz

## Examples

``` r
## Fitting a 4-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## Calculating EC/ED values
ED(ryegrass.m1, c(10, 50, 90))
#> 
#> Estimated effective doses
#> 
#>        Estimate Std. Error
#> e:1:10  1.46371    0.18677
#> e:1:50  3.05795    0.18573
#> e:1:90  6.38864    0.84510

## Displaying 95% confidence intervals using the delta method
ED(ryegrass.m1, c(10, 50, 90), interval = "delta")
#> 
#> Estimated effective doses
#> 
#>        Estimate Std. Error   Lower   Upper
#> e:1:10  1.46371    0.18677 1.07411 1.85330
#> e:1:50  3.05795    0.18573 2.67053 3.44538
#> e:1:90  6.38864    0.84510 4.62580 8.15148

## Displaying 95% confidence intervals using back-transformation
ED(ryegrass.m1, c(10, 50, 90), interval = "fls")
#> 
#> Estimated effective doses
#> 
#>         Estimate     Lower     Upper
#> e:1:10    4.3219    2.9274    6.3809
#> e:1:50   21.2840   14.4476   31.3553
#> e:1:90  595.0468  102.0842 3468.5164

## Displaying 95% confidence intervals using inverse regression
ED(ryegrass.m1, c(10, 50, 90), interval = "inv")
#> 
#> Estimated effective doses
#> 
#>        Estimate  Lower  Upper
#> e:1:10   1.4637 1.1423 1.8225
#> e:1:50   3.0580 2.7490 3.4017
#> e:1:90   6.3886 5.1514 8.1965
```

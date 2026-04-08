# Estimating effective doses

Default method for class `drc`. `ED.drc` estimates effective
concentrations (EC) or effective doses (ED) for one or more specified
response levels. Response levels may be given as relative percentages of
the response range (e.g. ED50 = 50\\ values. The function computes point
estimates, delta-method standard errors, and optional confidence
intervals for each combination of curve and response level in the fitted
model.

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

The function carries out the following computational steps:

1.  **Input validation.** Arguments are checked for correct types and
    ranges (e.g. `respLev` must be numeric, `level` must be in (0, 1),
    and relative response levels must lie strictly inside the interval
    (0, 100) when `bound = TRUE`).

2.  **Model component extraction.** The model-specific ED function
    (`edfct`), parameter matrix (`parmMat`), and index matrix
    (`indexMat`) are retrieved from the fitted `drc` object. The
    variance-covariance matrix is obtained from `vcov.`, which may be a
    function (e.g. [`vcov`](https://rdrr.io/r/stats/vcov.html) or
    [`sandwich::vcovHC`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html))
    or a pre-computed matrix.

3.  **Curve ordering.** When multiple curves are present, they are
    sorted alphabetically by name, unless the names are purely numeric,
    in which case the original order is preserved.

4.  **ED estimation and delta-method standard errors.** For each curve
    and each requested response level, the model-specific `edfct` is
    called to obtain the ED point estimate and its analytical gradient
    with respect to the model parameters. Standard errors are then
    computed via the delta method: \\SE = \sqrt{g' V g}\\, where \\g\\
    is the gradient vector and \\V\\ is the relevant sub-matrix of the
    variance-covariance matrix.

5.  **Numerical gradient for absolute responses.** When
    `type = "absolute"`, the analytical gradient returned by the model
    may miss the chain-rule contribution from the asymptote parameters
    involved in converting absolute to relative response levels. In that
    case a numerical central-difference gradient is computed to ensure
    correct standard errors.

6.  **Log-base back-transformation.** If `logBase` is specified
    (indicating that dose values were log-transformed prior to model
    fitting), the ED estimates and their derivatives are
    back-transformed via \\ED^\* = b^{ED}\\ (where \\b\\ is the log
    base) so that results are reported on the original dose scale.

7.  **Confidence interval construction.** Depending on `interval`:

    - `"delta"`:

      Asymptotic Wald-type intervals using the delta method, based on
      the normal or t-distribution (depending on the response type).

    - `"fls"`:

      Intervals obtained by back-transforming from the log scale. Only
      meaningful when the model parameterises the ED on the log scale
      (e.g.
      [`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md)).

    - `"tfls"`:

      Experimental: intervals obtained by transforming to the log scale,
      computing Wald intervals there, then back-transforming.

    - `"inv"`:

      Intervals derived from inverse regression via
      [`EDinvreg`](https://hreinwald.github.io/drc/reference/EDinvreg.md),
      where confidence limits on the predicted response are inverted to
      the dose axis.

8.  **Output.** Results are returned as an invisible matrix with columns
    for the estimate, standard error, and (optionally) lower and upper
    confidence limits. When `multcomp = TRUE`, a list compatible with
    [`parm`](https://rdrr.io/pkg/multcomp/man/parm.html) is returned
    instead, enabling multiple-comparison procedures.

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
#>      Estimate Std. Error
#> e:10  1.46371    0.18677
#> e:50  3.05795    0.18573
#> e:90  6.38864    0.84510

## Displaying 95% confidence intervals using the delta method
ED(ryegrass.m1, c(10, 50, 90), interval = "delta")
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error   Lower   Upper
#> e:10  1.46371    0.18677 1.07411 1.85330
#> e:50  3.05795    0.18573 2.67053 3.44538
#> e:90  6.38864    0.84510 4.62580 8.15148

## Displaying 95% confidence intervals using back-transformation
ED(ryegrass.m1, c(10, 50, 90), interval = "fls")
#> 
#> Estimated effective doses
#> 
#>       Estimate     Lower     Upper
#> e:10    4.3219    2.9274    6.3809
#> e:50   21.2840   14.4476   31.3553
#> e:90  595.0468  102.0842 3468.5164

## Displaying 95% confidence intervals using inverse regression
ED(ryegrass.m1, c(10, 50, 90), interval = "inv")
#> 
#> Estimated effective doses
#> 
#>      Estimate  Lower  Upper
#> e:10   1.4637 1.1423 1.8225
#> e:50   3.0580 2.7490 3.4017
#> e:90   6.3886 5.1514 8.1965
```

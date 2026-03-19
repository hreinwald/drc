# Maximum mean response

Estimates the maximum mean response and the dose at which it occurs,
using a bisection method to locate the peak of the fitted dose-response
curve. This function is only implemented for the built-in model
functions of class
[`braincousens`](https://hreinwald.github.io/drc/reference/braincousens.md)
and
[`cedergreen`](https://hreinwald.github.io/drc/reference/cedergreen.md),
which are capable of exhibiting hormesis (i.e., a non-monotone response
with a stimulatory effect at low doses).

## Usage

``` r
MAX(object, lower = 0.001, upper = 1000, pool = TRUE)
```

## Arguments

- object:

  an object of class `drc`, fitted using
  [`drm`](https://hreinwald.github.io/drc/reference/drm.md) with a
  hormesis model such as
  [`CRS.4c`](https://hreinwald.github.io/drc/reference/CRS.4c.md) or
  [`BC.4`](https://hreinwald.github.io/drc/reference/BC.4.md).

- lower:

  numeric. Lower bound of the interval used by the bisection method to
  search for the dose at maximum response. Must be strictly smaller than
  `upper` and should be set below the expected dose at maximum response.
  Defaults to `1e-3`.

- upper:

  numeric. Upper bound of the interval used by the bisection method to
  search for the dose at maximum response. Must be strictly larger than
  `lower` and should be set above the expected dose at maximum response.
  Defaults to `1000`.

- pool:

  logical. If `TRUE` (default), curves are pooled when computing the
  variance-covariance matrix. Otherwise they are not. This argument only
  works for models with independently fitted curves as specified in
  [`drm`](https://hreinwald.github.io/drc/reference/drm.md). Note:
  currently the variance-covariance matrix is retrieved for internal
  consistency but standard errors are not yet reported in the output.

## Value

Invisibly returns a numeric matrix with one row per curve in the data
set and two columns:

- Dose:

  The dose at which the maximum mean response occurs, found via
  bisection within `[lower, upper]`.

- Response:

  The estimated maximum mean response at that dose.

Row names correspond to curve identifiers. If the computation fails for
a given curve, the corresponding row will contain `NA` values and a
warning is issued. The matrix is also printed to the console via
[`printCoefmat`](https://rdrr.io/r/stats/printCoefmat.html).

## Details

The function numerically locates the dose \\d^\*\\ that maximises the
fitted dose-response curve over the search interval `[lower, upper]`:
\$\$d^\* = \arg\max\_{d} f(d, \hat{\theta})\$\$ where \\f\\ is the
fitted dose-response function and \\\hat{\theta}\\ is the vector of
estimated parameters. The search is performed using a bisection approach
defined internally by the model's `maxfct` component.

It is the user's responsibility to ensure that the true maximum lies
within `[lower, upper]`. If the maximum falls outside this interval, the
function will silently return a boundary value and a warning is issued.

## References

Cedergreen, N., Ritz, C., and Streibig, J. C. (2005) Improved empirical
models describing hormesis, *Environmental Toxicology and Chemistry*
**24**, 3166–3172.

## Author

Christian Ritz. Issues fixed and documentation enhanced by Hannes
Reinwald.

## Examples

``` r
## Fitting a Cedergreen-Ritz-Streibig model
lettuce.m1 <- drm(weight ~ conc, data = lettuce, fct = CRS.4c())

## Finding the maximum mean response and the corresponding dose
MAX(lettuce.m1)
#>      Dose Response
#> 1 0.25587    1.178

## Custom search interval
MAX(lettuce.m1, lower = 1e-5, upper = 500)
#>      Dose Response
#> 1 0.25587    1.178

## Capture the result matrix
result <- MAX(lettuce.m1)
#>      Dose Response
#> 1 0.25587    1.178
result["Dose"]
#> [1] NA
```

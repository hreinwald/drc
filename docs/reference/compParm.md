# Comparison of parameters

Compare parameters from different assays, either by means of ratios or
differences.

## Usage

``` r
compParm(
  object,
  strVal,
  operator = "/",
  vcov. = vcov,
  od = FALSE,
  pool = TRUE,
  display = TRUE
)
```

## Arguments

- object:

  an object of class 'drc'.

- strVal:

  a name of parameter to compare.

- operator:

  a character. If equal to `"/"` (default) parameter ratios are
  compared. If equal to `"-"` parameter differences are compared.

- vcov.:

  function providing the variance-covariance matrix.
  [`vcov`](https://rdrr.io/r/stats/vcov.html) is the default, but
  `sandwich` is also an option (for obtaining robust standard errors).

- od:

  logical. If TRUE adjustment for over-dispersion is used.

- pool:

  logical. If TRUE curves are pooled. Otherwise they are not. This
  argument only works for models with independently fitted curves as
  specified in
  [`drm`](https://hreinwald.github.io/drc/reference/drm.md).

- display:

  logical. If TRUE results are displayed. Otherwise they are not (useful
  in simulations).

## Value

A matrix with columns containing the estimates, estimated standard
errors, values of t-statistics and p-values for the null hypothesis that
the ratio equals 1 or that the difference equals 0 (depending on the
`operator` argument).

## See also

[`ED.drc`](https://hreinwald.github.io/drc/reference/ED.drc.md) for
calculating effective doses and
[`EDcomp`](https://hreinwald.github.io/drc/reference/EDcomp.md) for
comparing effective doses.

## Author

Christian Ritz

## Examples

``` r
spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach,
fct = LL.4(names = c("b", "lower", "upper", "ed50")))

## Calculating ratios of parameter estimates for "ed50"
compParm(spinach.m1, "ed50")
#> 
#> Comparison of parameter 'ed50' 
#> 
#>     Estimate Std. Error t-value  p-value   
#> 1/2  1.89836    0.71185  1.2620 0.210398   
#> 1/3  1.30730    0.55416  0.5545 0.580668   
#> 1/4  9.09638    2.46866  3.2797 0.001508 **
#> 1/5  8.51523    2.33645  3.2165 0.001836 **
#> 2/3  0.68865    0.29081 -1.0706 0.287360   
#> 2/4  4.79171    1.28835  2.9431 0.004189 **
#> 2/5  4.48557    1.21960  2.8580 0.005361 **
#> 3/4  6.95813    2.32206  2.5659 0.012045 * 
#> 3/5  6.51359    2.18960  2.5181 0.013674 * 
#> 4/5  0.93611    0.07814 -0.8176 0.415868   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Calculating differences between parameter estimates for "ed50"
compParm(spinach.m1, "ed50", "-")
#> 
#> Comparison of parameter 'ed50' 
#> 
#>      Estimate Std. Error t-value  p-value   
#> 1-2  0.849425   0.539400  1.5748 0.119027   
#> 1-3  0.421932   0.658506  0.6407 0.523414   
#> 1-4  1.597629   0.478341  3.3399 0.001246 **
#> 1-5  1.584161   0.478432  3.3112 0.001365 **
#> 2-3 -0.427493   0.516885 -0.8271 0.410521   
#> 2-4  0.748204   0.249701  2.9964 0.003580 **
#> 2-5  0.734736   0.249876  2.9404 0.004221 **
#> 3-4  1.175696   0.452799  2.5965 0.011095 * 
#> 3-5  1.162229   0.452896  2.5662 0.012034 * 
#> 4-5 -0.013467   0.017174 -0.7842 0.435128   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

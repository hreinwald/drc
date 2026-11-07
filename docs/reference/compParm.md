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
#> 1/2  1.89833    0.71186  1.2620 0.210418   
#> 1/3  1.30726    0.55417  0.5545 0.580720   
#> 1/4  9.09644    2.46873  3.2796 0.001508 **
#> 1/5  8.51520    2.33649  3.2164 0.001837 **
#> 2/3  0.68864    0.29082 -1.0706 0.287363   
#> 2/4  4.79182    1.28842  2.9430 0.004189 **
#> 2/5  4.48563    1.21966  2.8579 0.005362 **
#> 3/4  6.95838    2.32227  2.5658 0.012049 * 
#> 3/5  6.51376    2.18978  2.5179 0.013679 * 
#> 4/5  0.93610    0.07814 -0.8177 0.415799   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Calculating differences between parameter estimates for "ed50"
compParm(spinach.m1, "ed50", "-")
#> 
#> Comparison of parameter 'ed50' 
#> 
#>      Estimate Std. Error t-value  p-value   
#> 1-2  0.849417   0.539420  1.5747 0.119044   
#> 1-3  0.421896   0.658545  0.6406 0.523475   
#> 1-4  1.597643   0.478355  3.3399 0.001246 **
#> 1-5  1.584174   0.478446  3.3111 0.001365 **
#> 2-3 -0.427521   0.516930 -0.8270 0.410531   
#> 2-4  0.748226   0.249715  2.9963 0.003581 **
#> 2-5  0.734757   0.249890  2.9403 0.004222 **
#> 3-4  1.175747   0.452842  2.5964 0.011099 * 
#> 3-5  1.162278   0.452939  2.5661 0.012038 * 
#> 4-5 -0.013469   0.017174 -0.7843 0.435064   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

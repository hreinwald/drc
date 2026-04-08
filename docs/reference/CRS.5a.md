# Cedergreen-Ritz-Streibig Five-Parameter Model with Alpha = 1 (Deprecated)

**\[deprecated\]**

This function is deprecated as of version 3.3.0. Please use
[`CRS.5()`](https://hreinwald.github.io/drc/reference/CRS.5.md) instead,
which provides a more general and flexible interface.

A five-parameter Cedergreen-Ritz-Streibig (CRS) hormesis model where the
alpha parameter controlling the steepness of the hormetic component is
fixed at 1. All five parameters `b`, `c`, `d`, `e`, and `f` are freely
estimated.

## Usage

``` r
CRS.5a(names = c("b", "c", "d", "e", "f"), fixed = c(NA, NA, NA, NA, NA), ...)
```

## Arguments

- names:

  A character vector of length 5 specifying the names of the model
  parameters in the following order:

  `b`

  :   Hill slope (steepness of the dose-response curve).

  `c`

  :   Lower asymptote (freely estimated).

  `d`

  :   Upper asymptote.

  `e`

  :   Effective dose producing a response midway between `c` and `d`
      (ED50).

  `f`

  :   Hormesis parameter controlling the magnitude of the stimulatory
      effect at low doses.

  Defaults to `c("b", "c", "d", "e", "f")`.

- fixed:

  A numeric vector of length 5 specifying fixed (non-estimated)
  parameter values. Use `NA` for parameters that should be estimated
  freely. Defaults to `c(NA, NA, NA, NA, NA)`, meaning all five
  parameters are freely estimated.

- ...:

  Additional arguments passed to
  [`cedergreen()`](https://hreinwald.github.io/drc/reference/cedergreen.md).

## Value

A list of class `"drcMean"` as returned by
[`cedergreen()`](https://hreinwald.github.io/drc/reference/cedergreen.md),
containing the model definition including the mean function, its
gradient, parameter names, and fixed values. This object is intended for
use as the `fct` argument in
[`drm()`](https://hreinwald.github.io/drc/reference/drm.md).

## See also

- [`CRS.5()`](https://hreinwald.github.io/drc/reference/CRS.5.md) — the
  recommended replacement for this deprecated function.

- [`cedergreen()`](https://hreinwald.github.io/drc/reference/cedergreen.md)
  — the underlying model constructor.

- [`CRS.4a()`](https://hreinwald.github.io/drc/reference/CRS.4a.md) —
  the four-parameter CRS model with lower limit fixed at 0 and alpha =
  1.

- [`UCRS.5a()`](https://hreinwald.github.io/drc/reference/UCRS.5a.md) —
  the unconstrained five-parameter CRS model with alpha = 1.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
# NOTE: CRS.5a() is deprecated. Use CRS.5() instead.
# The example below is retained for backward compatibility illustration only.

lettuce.m1 <- drm( lettuce[, c(2, 1)], fct = CRS.5a() )
summary(lettuce.m1)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig (alpha=1) (5 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 1.334173   0.358675  3.7197  0.004773 ** 
#> c:(Intercept) 0.447962   0.080700  5.5510  0.000356 ***
#> d:(Intercept) 1.035658   0.077323 13.3940 3.004e-07 ***
#> e:(Intercept) 1.337869   1.185153  1.1289  0.288148    
#> f:(Intercept) 1.993259   2.017541  0.9880  0.348985    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1305067 (9 degrees of freedom)
ED(lettuce.m1, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50   5.5439     1.9480

# Recommended replacement:
lettuce.crs5 <- drm( lettuce[, c(2, 1)], fct = CRS.5(alpha_type = "a") )
summary(lettuce.crs5)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig (alpha=1) (5 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 1.334173   0.358675  3.7197  0.004773 ** 
#> c:(Intercept) 0.447962   0.080700  5.5510  0.000356 ***
#> d:(Intercept) 1.035658   0.077323 13.3940 3.004e-07 ***
#> e:(Intercept) 1.337869   1.185153  1.1289  0.288148    
#> f:(Intercept) 1.993259   2.017541  0.9880  0.348985    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1305067 (9 degrees of freedom)
ED(lettuce.crs5, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50   5.5439     1.9480
```

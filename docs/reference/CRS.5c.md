# Cedergreen-Ritz-Streibig Five-Parameter Model with Alpha = 0.25 (Deprecated)

**\[deprecated\]**

This function is deprecated as of version 3.3.0. Please use
[`CRS.5()`](https://hreinwald.github.io/drc/reference/CRS.5.md) instead,
which provides a more general and flexible interface.

A five-parameter Cedergreen-Ritz-Streibig (CRS) hormesis model where the
alpha parameter controlling the steepness of the hormetic component is
fixed at 0.25. All five parameters `b`, `c`, `d`, `e`, and `f` are
freely estimated.

## Usage

``` r
CRS.5c(names = c("b", "c", "d", "e", "f"), fixed = c(NA, NA, NA, NA, NA), ...)
```

## Arguments

- names:

  A character vector of length 5 specifying the names of the model
  parameters in the following order:

  `b`

  : Hill slope (steepness of the dose-response curve).

  `c`

  : Lower asymptote (freely estimated).

  `d`

  : Upper asymptote.

  `e`

  : Effective dose producing a response midway between `c` and `d`
    (ED50).

  `f`

  : Hormesis parameter controlling the magnitude of the stimulatory
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

- [`CRS.4c()`](https://hreinwald.github.io/drc/reference/CRS.4c.md) —
  the four-parameter CRS model with lower limit fixed at 0 and alpha =
  0.25.

- [`CRS.5b()`](https://hreinwald.github.io/drc/reference/CRS.5b.md) —
  the five-parameter CRS model with alpha = 0.5.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
# NOTE: CRS.5c() is deprecated. Use CRS.5() instead.
# The example below is retained for backward compatibility illustration only.

lettuce.m3 <- drm( lettuce[, c(2, 1)], fct = CRS.5c() )
plot(lettuce.m3, type = "all", main = "CRS.5c (deprecated)")

summary(lettuce.m3)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig (alpha=0.25) (5 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.975182   0.562159  1.7347    0.1168    
#> c:(Intercept) 0.335014   0.185691  1.8041    0.1047    
#> d:(Intercept) 0.969038   0.088264 10.9788 1.636e-06 ***
#> e:(Intercept) 3.856221   2.477090  1.5568    0.1540    
#> f:(Intercept) 1.040418   0.783943  1.3272    0.2171    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1256823 (9 degrees of freedom)
ED(lettuce.m3, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50  11.5195     8.9099

# Recommended replacement:
lettuce.crs5 <- drm( lettuce[, c(2, 1)], fct = CRS.5(alpha_type = "c") )
summary(lettuce.crs5)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig (alpha=0.25) (5 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.975182   0.562159  1.7347    0.1168    
#> c:(Intercept) 0.335014   0.185691  1.8041    0.1047    
#> d:(Intercept) 0.969038   0.088264 10.9788 1.636e-06 ***
#> e:(Intercept) 3.856221   2.477090  1.5568    0.1540    
#> f:(Intercept) 1.040418   0.783943  1.3272    0.2171    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1256823 (9 degrees of freedom)
ED(lettuce.crs5, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50  11.5195     8.9099
```

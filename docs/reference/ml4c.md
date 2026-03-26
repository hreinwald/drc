# Alias for CRS.5c (Deprecated)

**\[deprecated\]**

This function is a deprecated alias for
[`CRS.5c()`](https://hreinwald.github.io/drc/reference/CRS.5c.md),
itself deprecated as of version 3.3.0. Please use
[`CRS.5()`](https://hreinwald.github.io/drc/reference/CRS.5.md) instead,
which provides a more general and flexible interface.

## Usage

``` r
ml4c(names = c("b", "c", "d", "e", "f"), fixed = c(NA, NA, NA, NA, NA), ...)
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

- [`CRS.5c()`](https://hreinwald.github.io/drc/reference/CRS.5c.md) —
  the function this alias points to.

- [`cedergreen()`](https://hreinwald.github.io/drc/reference/cedergreen.md)
  — the underlying model constructor.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
# NOTE: ml4c() is a deprecated alias for CRS.5c(). Use CRS.5() instead.
# The example below is retained for backward compatibility illustration only.

lettuce.m3 <- drm( lettuce[, c(2, 1)], fct = ml4c() )
summary(lettuce.m3)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig (alpha=0.25) (5 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.981945   0.559334  1.7556   0.11305    
#> c:(Intercept) 0.336670   0.182883  1.8409   0.09877 .  
#> d:(Intercept) 0.969845   0.088261 10.9883 1.624e-06 ***
#> e:(Intercept) 3.883893   2.462313  1.5773   0.14917    
#> f:(Intercept) 1.027934   0.766823  1.3405   0.21293    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1256841 (9 degrees of freedom)
ED(lettuce.m3, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50  11.4243     8.7214

# Recommended replacement:
lettuce.crs5 <- drm( lettuce[, c(2, 1)], fct = CRS.5(alpha_type = "c") )
summary(lettuce.crs5)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig (alpha=0.25) (5 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.981945   0.559334  1.7556   0.11305    
#> c:(Intercept) 0.336670   0.182883  1.8409   0.09877 .  
#> d:(Intercept) 0.969845   0.088261 10.9883 1.624e-06 ***
#> e:(Intercept) 3.883893   2.462313  1.5773   0.14917    
#> f:(Intercept) 1.027934   0.766823  1.3405   0.21293    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1256841 (9 degrees of freedom)
ED(lettuce.crs5, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50  11.4243     8.7214
```

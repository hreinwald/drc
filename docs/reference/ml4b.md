# Alias for CRS.5b (Deprecated)

**\[deprecated\]**

This function is a deprecated alias for
[`CRS.5b()`](https://hreinwald.github.io/drc/reference/CRS.5b.md),
itself deprecated as of version 3.3.0. Please use
[`CRS.5()`](https://hreinwald.github.io/drc/reference/CRS.5.md) instead,
which provides a more general and flexible interface.

## Usage

``` r
ml4b(names = c("b", "c", "d", "e", "f"), fixed = c(NA, NA, NA, NA, NA), ...)
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

- [`CRS.5b()`](https://hreinwald.github.io/drc/reference/CRS.5b.md) —
  the function this alias points to.

- [`cedergreen()`](https://hreinwald.github.io/drc/reference/cedergreen.md)
  — the underlying model constructor.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
# NOTE: ml4b() is a deprecated alias for CRS.5b(). Use CRS.5() instead.
# The example below is retained for backward compatibility illustration only.

lettuce.m2 <- drm( lettuce[, c(2, 1)], fct = ml4b() )
summary(lettuce.m2)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig (alpha=.5) (5 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.806096   0.537800  1.4989    0.1681    
#> c:(Intercept) 0.316586   0.199024  1.5907    0.1461    
#> d:(Intercept) 0.971581   0.081936 11.8577 8.523e-07 ***
#> e:(Intercept) 0.814111   2.969068  0.2742    0.7901    
#> f:(Intercept) 3.288976   8.216399  0.4003    0.6983    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1167711 (9 degrees of freedom)
ED(lettuce.m2, c(50))
#> 
#> Estimated effective doses
#> 
#>        Estimate Std. Error
#> e:1:50   11.550      8.603

# Recommended replacement:
lettuce.crs5 <- drm( lettuce[, c(2, 1)], fct = CRS.5(alpha_type = "b") )
#> Error in CRS.5(alpha_type = "b"): could not find function "CRS.5"
summary(lettuce.crs5)
#> Error: object 'lettuce.crs5' not found
ED(lettuce.crs5, c(50))
#> Error: object 'lettuce.crs5' not found
```

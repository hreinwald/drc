# Cedergreen-Ritz-Streibig Model with Lower Limit Fixed at 0 and Alpha = 0.25 (Deprecated)

**\[deprecated\]**

This function is deprecated as of version 3.3.0. Please use
[`CRS.5()`](https://hreinwald.github.io/drc/reference/CRS.5.md) instead,
which provides a more general and flexible interface.

A four-parameter Cedergreen-Ritz-Streibig (CRS) hormesis model where the
lower asymptote (`c`) is fixed at 0 and the alpha parameter controlling
the steepness of the hormetic component is fixed at 0.25. The four free
parameters are `b`, `d`, `e`, and `f`.

## Usage

``` r
CRS.4c(names = c("b", "c", "d", "e", "f"), fixed = c(NA, 0, NA, NA, NA), ...)
```

## Arguments

- names:

  A character vector of length 5 specifying the names of the model
  parameters in the following order:

  `b`

  :   Hill slope (steepness of the dose-response curve).

  `c`

  :   Lower asymptote (fixed at 0 via the `fixed` argument).

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
  freely. Defaults to `c(NA, 0, NA, NA, NA)`, which fixes the lower
  asymptote `c` at 0.

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

- [`CRS.5c()`](https://hreinwald.github.io/drc/reference/CRS.5c.md) —
  the five-parameter CRS model with alpha = 0.25.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
# NOTE: CRS.4c() is deprecated. Use CRS.5() instead.
# The example below is retained for backward compatibility illustration only.

lettuce.crsm3 <- drm( lettuce[, c(2, 1)], fct = CRS.4c() )
summary(lettuce.crsm3)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig with lower limit 0 (alpha=0.25) (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.488440   0.133643  3.6548  0.004427 ** 
#> d:(Intercept) 0.973666   0.086883 11.2066 5.544e-07 ***
#> e:(Intercept) 1.314657   3.614266  0.3637  0.723624    
#> f:(Intercept) 2.998547   3.626210  0.8269  0.427579    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.123575 (10 degrees of freedom)
ED(lettuce.crsm3, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50   37.033     15.437

# Recommended replacement:
fct_spec <- CRS.5(alpha_type = "c", fixed = c(NA, 0, NA, NA, NA))
lettuce.crs5 <- drm(lettuce[, c(2, 1)], fct = fct_spec)
summary(lettuce.crs5)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig (alpha=0.25) (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.488440   0.133643  3.6548  0.004427 ** 
#> d:(Intercept) 0.973666   0.086883 11.2066 5.544e-07 ***
#> e:(Intercept) 1.314657   3.614266  0.3637  0.723624    
#> f:(Intercept) 2.998547   3.626210  0.8269  0.427579    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.123575 (10 degrees of freedom)
ED(lettuce.crs5, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50   37.033     15.437
```

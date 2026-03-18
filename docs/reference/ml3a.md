# Alias for CRS.4a (Deprecated)

**\[deprecated\]**

This function is a deprecated alias for
[`CRS.4a()`](https://hreinwald.github.io/drc/reference/CRS.4a.md),
itself deprecated as of version 3.3.0. Please use
[`CRS.5()`](https://hreinwald.github.io/drc/reference/CRS.5.md) instead,
which provides a more general and flexible interface.

## Usage

``` r
ml3a(names = c("b", "c", "d", "e", "f"), fixed = c(NA, 0, NA, NA, NA), ...)
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

- [`CRS.4a()`](https://hreinwald.github.io/drc/reference/CRS.4a.md) —
  the function this alias points to.

- [`cedergreen()`](https://hreinwald.github.io/drc/reference/cedergreen.md)
  — the underlying model constructor.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
# NOTE: ml3a() is a deprecated alias for CRS.4a(). Use CRS.5() instead.
# The example below is retained for backward compatibility illustration only.

lettuce.crsm1 <- drm( lettuce[, c(2, 1)], fct = ml3a() )
summary(lettuce.crsm1)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig with lower limit 0 (alpha=1) (4 parms)
#> 
#> Parameter estimates:
#> 
#>                Estimate Std. Error t-value   p-value    
#> b:(Intercept)  0.774519   0.248592  3.1156   0.01096 *  
#> d:(Intercept)  1.108705   0.078481 14.1270 6.212e-08 ***
#> e:(Intercept) 27.620019  30.307666  0.9113   0.38357    
#> f:(Intercept)  0.013090   0.417215  0.0314   0.97559    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1556406 (10 degrees of freedom)
ED(lettuce.crsm1, c(50))
#> 
#> Estimated effective doses
#> 
#>        Estimate Std. Error
#> e:1:50   28.436     11.618

# Recommended replacement:
fct_spec <- CRS.5(alpha_type = "a", fixed = c(NA, 0, NA, NA, NA))
#> Error in CRS.5(alpha_type = "a", fixed = c(NA, 0, NA, NA, NA)): could not find function "CRS.5"
lettuce.crs5 <- drm(lettuce[, c(2, 1)], fct = fct_spec)
#> Error: object 'fct_spec' not found
summary(lettuce.crs5)
#> Error: object 'lettuce.crs5' not found
ED(lettuce.crs5, c(50))
#> Error: object 'lettuce.crs5' not found
```

# Alias for CRS.4c (Deprecated)

**\[deprecated\]**

This function is a deprecated alias for
[`CRS.4c()`](https://hreinwald.github.io/drc/reference/CRS.4c.md),
itself deprecated as of version 3.3.0. Please use
[`CRS.5()`](https://hreinwald.github.io/drc/reference/CRS.5.md) instead,
which provides a more general and flexible interface.

## Usage

``` r
ml3c(names = c("b", "c", "d", "e", "f"), fixed = c(NA, 0, NA, NA, NA), ...)
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

- [`CRS.4c()`](https://hreinwald.github.io/drc/reference/CRS.4c.md) —
  the function this alias points to.

- [`cedergreen()`](https://hreinwald.github.io/drc/reference/cedergreen.md)
  — the underlying model constructor.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
# NOTE: ml3c() is a deprecated alias for CRS.4c(). Use CRS.5() instead.
# The example below is retained for backward compatibility illustration only.

lettuce.crsm3 <- drm( lettuce[, c(2, 1)], fct = ml3c() )
summary(lettuce.crsm3)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig with lower limit 0 (alpha=.25) (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.490601   0.137206  3.5757  0.005048 ** 
#> d:(Intercept) 0.973940   0.086907 11.2067 5.543e-07 ***
#> e:(Intercept) 1.378741   3.812958  0.3616  0.725179    
#> f:(Intercept) 2.935678   3.597993  0.8159  0.433550    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1235846 (10 degrees of freedom)
ED(lettuce.crsm3, c(50))
#> 
#> Estimated effective doses
#> 
#>        Estimate Std. Error
#> e:1:50   36.935     15.377

# Recommended replacement:
fct_spec <- CRS.5(alpha_type = "c", fixed = c(NA, 0, NA, NA, NA))
#> Error in CRS.5(alpha_type = "c", fixed = c(NA, 0, NA, NA, NA)): could not find function "CRS.5"
lettuce.crs5 <- drm(lettuce[, c(2, 1)], fct = fct_spec)
#> Error: object 'fct_spec' not found
summary(lettuce.crs5)
#> Error: object 'lettuce.crs5' not found
ED(lettuce.crs5, c(50))
#> Error: object 'lettuce.crs5' not found
```

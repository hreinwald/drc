# Cedergreen-Ritz-Streibig Model with Lower Limit Fixed at 0 and Alpha = 0.5 (Deprecated)

**\[deprecated\]**

This function is deprecated as of version 3.3.0. Please use
[`CRS.5()`](https://hreinwald.github.io/drc/reference/CRS.5.md) instead,
which provides a more general and flexible interface.

A four-parameter Cedergreen-Ritz-Streibig (CRS) hormesis model where the
lower asymptote (`c`) is fixed at 0 and the alpha parameter controlling
the steepness of the hormetic component is fixed at 0.5. The four free
parameters are `b`, `d`, `e`, and `f`.

## Usage

``` r
CRS.4b(names = c("b", "c", "d", "e", "f"), fixed = c(NA, 0, NA, NA, NA), ...)
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

- [`CRS.5b()`](https://hreinwald.github.io/drc/reference/CRS.5b.md) —
  the five-parameter CRS model with alpha = 0.5.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
# NOTE: CRS.4b() is deprecated. Use CRS.5() instead.
# The example below is retained for backward compatibility illustration only.

lettuce.crsm2 <- drm( lettuce[, c(2, 1)], fct = CRS.4b() )
summary(lettuce.crsm2)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig with lower limit 0 (alpha=0.5) (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.569426   0.068538  8.3081 8.444e-06 ***
#> d:(Intercept) 1.008915   0.094919 10.6292 9.061e-07 ***
#> e:(Intercept) 0.642290   1.533937  0.4187    0.6843    
#> f:(Intercept) 4.446933   5.821389  0.7639    0.4626    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1345066 (10 degrees of freedom)
ED(lettuce.crsm2, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50  26.1252     8.6286

# Recommended replacement:
fct_spec <- CRS.5(alpha_type = "b", fixed = c(NA, 0, NA, NA, NA))
lettuce.crs5 <- drm(lettuce[, c(2, 1)], fct = fct_spec)
summary(lettuce.crs5)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig (alpha=0.5) (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.569426   0.068538  8.3081 8.444e-06 ***
#> d:(Intercept) 1.008915   0.094919 10.6292 9.061e-07 ***
#> e:(Intercept) 0.642290   1.533937  0.4187    0.6843    
#> f:(Intercept) 4.446933   5.821389  0.7639    0.4626    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1345066 (10 degrees of freedom)
ED(lettuce.crs5, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50  26.1252     8.6286
```

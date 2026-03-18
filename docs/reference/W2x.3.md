# Three-parameter Weibull type 2 model with lag time

A three-parameter Weibull type 2 model with lag time, where `b` is fixed
at 1 and `c` is fixed at 0. This is a convenience wrapper around
[`weibull2x`](https://hreinwald.github.io/drc/reference/weibull2x.md).

## Usage

``` r
W2x.3(fixed = c(NA, NA, NA), names = c("d", "e", "t0"), ...)
```

## Arguments

- fixed:

  numeric vector of length 3. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated
  (default is `c(NA, NA, NA)`).

- names:

  character vector of length 3 giving the names of the parameters
  (default is `c("d", "e", "t0")`).

- ...:

  additional arguments passed to
  [`weibull2x`](https://hreinwald.github.io/drc/reference/weibull2x.md).

## Value

A list of class `"Weibull-2"` containing the nonlinear function, self
starter function, and parameter names.

## See also

[`weibull2x`](https://hreinwald.github.io/drc/reference/weibull2x.md),
[`W2x.4`](https://hreinwald.github.io/drc/reference/W2x.4.md),
[`W2.3`](https://hreinwald.github.io/drc/reference/W2.3.md)

## Examples

``` r
spinach.m1 <- drm(SLOPE ~ DOSE, data = spinach, fct = W2x.3())
summary(spinach.m1)
#> 
#> Model fitted: Weibull (type 2) with lower limit at 0 (3 parms)
#> 
#> Parameter estimates:
#> 
#>                 Estimate Std. Error t-value   p-value    
#> d:(Intercept)   0.827617   0.067677  12.229 < 2.2e-16 ***
#> e:(Intercept)   0.008972        NaN     NaN       NaN    
#> t0:(Intercept) -0.135527        NaN     NaN       NaN    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.6934874 (102 degrees of freedom)
```

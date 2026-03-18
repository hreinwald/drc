# Summarising non-linear model fits

`summary` compiles a comprehensive summary for objects of class 'drc'.

## Usage

``` r
# S3 method for class 'drc'
summary(object, od = FALSE, pool = TRUE, ...)
```

## Arguments

- object:

  an object of class 'drc'.

- od:

  logical. If TRUE adjustment for over-dispersion is used.

- pool:

  logical. If TRUE curves are pooled. Otherwise they are not. This
  argument only works for models with independently fitted curves as
  specified in
  [`drm`](https://hreinwald.github.io/drc/reference/drm.md).

- ...:

  additional arguments.

## Value

A list of summary statistics that includes parameter estimates and
estimated standard errors.

## See also

[`drm`](https://hreinwald.github.io/drc/reference/drm.md),
[`coef.drc`](https://hreinwald.github.io/drc/reference/coef.drc.md),
[`confint.drc`](https://hreinwald.github.io/drc/reference/confint.drc.md)

## Author

Christian Ritz

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
summary(ryegrass.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept)  2.98222    0.46506  6.4125 2.960e-06 ***
#> c:(Intercept)  0.48141    0.21219  2.2688   0.03451 *  
#> d:(Intercept)  7.79296    0.18857 41.3272 < 2.2e-16 ***
#> e:(Intercept)  3.05795    0.18573 16.4644 4.268e-13 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.5196256 (20 degrees of freedom)
```

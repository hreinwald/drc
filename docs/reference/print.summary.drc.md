# Printing summary of non-linear model fits

This method produces formatted output of the summary statistics:
parameter estimates, estimated standard errors, z-test statistics and
corresponding p-values.

## Usage

``` r
# S3 method for class 'summary.drc'
print(x, ...)
```

## Arguments

- x:

  an object of class 'drc'.

- ...:

  additional arguments.

## Value

The object (argument `x`) is returned invisibly.

## Author

Christian Ritz

## Examples

``` r
ryegrass.m1 <- drm(rootl~conc, data=ryegrass, fct= LL.4())

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

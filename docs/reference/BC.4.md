# Four-parameter Brain-Cousens hormesis model

`BC.4` provides the Brain-Cousens modified log-logistic model with the
lower limit fixed at 0.

## Usage

``` r
BC.4(fixed = c(NA, NA, NA, NA), names = c("b", "d", "e", "f"), ...)
```

## Arguments

- fixed:

  numeric vector of length 4 specifying fixed parameters (NAs for free
  parameters).

- names:

  a vector of character strings giving the names of the parameters.

- ...:

  additional arguments passed to
  [`braincousens`](https://hreinwald.github.io/drc/reference/braincousens.md).

## Value

A list (see
[`braincousens`](https://hreinwald.github.io/drc/reference/braincousens.md)).

## References

van Ewijk, P. H. and Hoekstra, J. A. (1993) Calculation of the EC50 and
its Confidence Interval When Subtoxic Stimulus Is Present,
*Ecotoxicology and Environmental Safety*, **25**, 25–32.

## See also

[`braincousens`](https://hreinwald.github.io/drc/reference/braincousens.md),
[`BC.5`](https://hreinwald.github.io/drc/reference/BC.5.md)

## Author

Christian Ritz

## Examples

``` r
## Fit with drm() and plot
lettuce.bcm2 <- drm(weight ~ conc, data = lettuce, fct = BC.4())
plot(lettuce.bcm2, type = "all", main = "BC.4 (lower limit fixed at 0)")

summary(lettuce.bcm2)
#> 
#> Model fitted: Brain-Cousens (hormesis) with lower limit fixed at 0 (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 1.282745   0.049347 25.9945 1.633e-10 ***
#> d:(Intercept) 0.967195   0.077127 12.5403 1.929e-07 ***
#> e:(Intercept) 0.846276   0.435726  1.9422   0.08078 .  
#> f:(Intercept) 1.623790   0.982039  1.6535   0.12924    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1117921 (10 degrees of freedom)
ED(lettuce.bcm2, c(50))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:50   35.045     15.446
```

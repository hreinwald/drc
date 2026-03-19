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
lettuce.bcm2 <- drm(weight ~ conc, data = lettuce, fct = BC.4())
summary(lettuce.bcm2)
#> 
#> Model fitted: Brain-Cousens (hormesis) with lower limit fixed at 0 (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 1.282812   0.049346 25.9964 1.632e-10 ***
#> d:(Intercept) 0.967302   0.077123 12.5423 1.926e-07 ***
#> e:(Intercept) 0.847633   0.436093  1.9437   0.08059 .  
#> f:(Intercept) 1.620703   0.979711  1.6543   0.12908    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1117922 (10 degrees of freedom)
ED(lettuce.bcm2, c(50))
#> 
#> Estimated effective doses
#> 
#>        Estimate Std. Error
#> e:1:50   35.023     15.427
```

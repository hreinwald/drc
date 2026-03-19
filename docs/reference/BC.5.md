# Five-parameter Brain-Cousens hormesis model

`BC.5` provides the full five-parameter Brain-Cousens modified
log-logistic model for describing hormesis.

## Usage

``` r
BC.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
```

## Arguments

- fixed:

  numeric vector of length 5 specifying fixed parameters (NAs for free
  parameters).

- names:

  a vector of character strings giving the names of the parameters.

- ...:

  additional arguments passed to
  [`braincousens`](https://hreinwald.github.io/drc/reference/braincousens.md).

## Value

A list (see
[`braincousens`](https://hreinwald.github.io/drc/reference/braincousens.md)).

## See also

[`braincousens`](https://hreinwald.github.io/drc/reference/braincousens.md),
[`BC.4`](https://hreinwald.github.io/drc/reference/BC.4.md)

## Author

Christian Ritz

## Examples

``` r
lettuce.bcm1 <- drm(weight ~ conc, data = lettuce, fct = BC.5())
modelFit(lettuce.bcm1)
#> Lack-of-fit test
#> 
#>           ModelDf      RSS Df F value p value
#> ANOVA           7 0.088237                   
#> DRC model       9 0.118842  2  1.2140  0.3527
plot(lettuce.bcm1)

```

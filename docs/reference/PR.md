# Expected or predicted response

Returns the expected or predicted response for specified dose values.
This is a convenience function for easy access to predicted values.

## Usage

``` r
PR(object, xVec, ...)
```

## Arguments

- object:

  object of class `drc` obtained from fitting a dose-response model.

- xVec:

  numeric vector of dose values.

- ...:

  additional arguments passed to
  [`predict.drc`](https://hreinwald.github.io/drc/reference/predict.drc.md).

## Value

A numeric vector of predicted values or possibly a matrix of predicted
values and corresponding standard errors.

## See also

[`predict.drc`](https://hreinwald.github.io/drc/reference/predict.drc.md)

## Author

Christian Ritz after a suggestion from Andrew Kniss.

## Examples

``` r
ryegrass.m1 <- drm(ryegrass, fct = LL.4())
PR(ryegrass.m1, c(5, 10))
#>         5        10 
#> 1.8523337 0.6888809 
```

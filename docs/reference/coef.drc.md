# Extract Model Coefficients

Extract parameter estimates.

## Usage

``` r
# S3 method for class 'drc'
coef(object, ...)
```

## Arguments

- object:

  an object of class 'drc'.

- ...:

  additional arguments.

## Value

A vector of parameter coefficients which are extracted from the model
object `object`.

## Author

Christian Ritz

## Examples

``` r
## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
coef(ryegrass.m1)
#> b:(Intercept) c:(Intercept) d:(Intercept) e:(Intercept) 
#>     2.9822191     0.4814132     7.7929583     3.0579550 
```

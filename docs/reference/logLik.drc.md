# Extracting the log likelihood

`logLik` extracts the value of the log likelihood function evaluated at
the parameter estimates.

## Usage

``` r
# S3 method for class 'drc'
logLik(object, ...)
```

## Arguments

- object:

  an object of class 'drc'.

- ...:

  additional arguments.

## Value

The evaluated log likelihood as a numeric value and the corresponding
degrees of freedom as well as the number of observations as attributes.

## Author

Christian Ritz

## Examples

``` r
## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
logLik(ryegrass.m1)
#> 'log Lik.' -16.15514 (df=5)
```

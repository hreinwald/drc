# Estimating function for the sandwich estimator

Evaluates the estimating function (the "meat") for the sandwich
estimator of the variance-covariance matrix for objects of class 'drc'.

## Usage

``` r
# S3 method for class 'drc'
estfun(x, ...)
```

## Arguments

- x:

  object of class `drc`.

- ...:

  additional arguments. At the moment none are supported.

## Value

The estimating function evaluated at the data and the parameter
estimates. By default no clustering is assumed, corresponding to robust
standard errors under independence.

## Details

The details are provided by Zeileis (2006).

## References

Zeileis, A. (2006) Object-oriented Computation of Sandwich Estimators,
*J. Statist. Software*, **16**, Issue 9.

## Author

Christian Ritz

## Examples

``` r
## The lines below requires that the packages
## 'lmtest' and 'sandwich' are installed
# library(lmtest)
# library(sandwich)

# ryegrass.m1<-drm(rootl ~ conc, data = ryegrass, fct = LL.4())

# Standard summary output
# coeftest(ryegrass.m1)

# Output with robust standard errors
# coeftest(ryegrass.m1, vcov = sandwich)
```

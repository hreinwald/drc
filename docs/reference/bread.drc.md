# Bread for the sandwich estimator

Computes the "bread" (unscaled hessian) for the sandwich estimator of
the variance-covariance matrix for objects of class 'drc'.

## Usage

``` r
# S3 method for class 'drc'
bread(x, ...)
```

## Arguments

- x:

  object of class `drc`.

- ...:

  additional arguments. At the moment none are supported.

## Value

The unscaled hessian matrix.

## Details

The details are provided by Zeileis (2006).

## References

Zeileis, A. (2006) Object-oriented Computation of Sandwich Estimators,
*J. Statist. Software*, **16**, Issue 9.

## Author

Christian Ritz

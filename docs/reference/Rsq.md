# R-squared for dose-response models

Calculates and displays R-squared values for a fitted dose-response
model. For models with multiple curves, per-curve and total R-squared
values are returned.

## Usage

``` r
Rsq(object)
```

## Arguments

- object:

  an object of class 'drc'.

## Value

Invisibly returns a matrix of R-squared values. For single-curve models,
a 1x1 matrix. For multi-curve models, includes per-curve values and a
total R-squared.

## Author

Christian Ritz

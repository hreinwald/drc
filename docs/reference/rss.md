# Residual sum of squares for dose-response models

Calculates and displays the residual sum of squares (RSS) for a fitted
dose-response model. For models with multiple curves, per-curve and
total RSS values are returned.

## Usage

``` r
rss(object, print = TRUE)
```

## Arguments

- object:

  an object of class 'drc'.

- print:

  logical. If `TRUE` (the default), the RSS values are printed.

## Value

Invisibly returns a matrix of RSS values. For single-curve models, a 1x1
matrix. For multi-curve models, includes per-curve values and a total
RSS.

## See also

[`Rsq()`](https://hreinwald.github.io/drc/reference/Rsq.md) which uses
this function to compute R-squared.

## Author

Christian Ritz

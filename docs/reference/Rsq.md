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

## Details

R-squared is computed as \\1 - RSS / TSS\\ where RSS is the residual sum
of squares (obtained via
[`rss()`](https://hreinwald.github.io/drc/reference/rss.md)) and TSS is
the total sum of squares.

## See also

[`rss()`](https://hreinwald.github.io/drc/reference/rss.md) for the
underlying residual sum of squares.

## Author

Christian Ritz

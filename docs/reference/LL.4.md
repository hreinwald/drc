# Four-parameter log-logistic function

A four-parameter log-logistic function.

## Usage

``` r
LL.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)

l4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 4, specifying fixed parameters (use NA for
  non-fixed parameters).

- names:

  character vector of length 4, specifying the names of the parameters
  (default: b, c, d, e).

- ...:

  additional arguments to
  [`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Value

See
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Details

The four-parameter log-logistic function is given by the expression
\$\$f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}\$\$

## See also

[`LL.3`](https://hreinwald.github.io/drc/reference/LL.3.md),
[`LL.5`](https://hreinwald.github.io/drc/reference/LL.5.md),
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md)

## Author

Christian Ritz and Jens C. Streibig

## Examples

``` r
spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())
```

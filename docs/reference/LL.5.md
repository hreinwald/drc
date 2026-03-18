# Five-parameter log-logistic function

A five-parameter (generalized) log-logistic function. The function is
asymmetric when f differs from 1.

## Usage

``` r
LL.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)

l5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
```

## Arguments

- fixed:

  numeric vector of length 5, specifying fixed parameters (use NA for
  non-fixed parameters).

- names:

  character vector of length 5, specifying the names of the parameters
  (default: b, c, d, e, f).

- ...:

  additional arguments to
  [`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Value

See
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Details

The five-parameter log-logistic function is given by the expression
\$\$f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}\$\$

## See also

[`LL.3`](https://hreinwald.github.io/drc/reference/LL.3.md),
[`LL.4`](https://hreinwald.github.io/drc/reference/LL.4.md),
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md)

## Author

Christian Ritz

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.5())
```

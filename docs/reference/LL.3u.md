# Three-parameter log-logistic function with upper limit fixed

A three-parameter log-logistic function with upper limit fixed (default
1), primarily for use with binomial/quantal dose-response data.

## Usage

``` r
LL.3u(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)

l3u(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
```

## Arguments

- upper:

  numeric value, the fixed upper limit (default 1).

- fixed:

  numeric vector of length 3, specifying fixed parameters (use NA for
  non-fixed parameters).

- names:

  character vector of length 3, specifying the names of the parameters
  (default: b, c, e).

- ...:

  additional arguments to
  [`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Value

See
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Details

The three-parameter log-logistic function with upper limit fixed is
given by \$\$f(x) = c + \frac{upper-c}{1+\exp(b(\log(x)-\log(e)))}\$\$

## See also

[`LL.2`](https://hreinwald.github.io/drc/reference/LL.2.md),
[`LL.3`](https://hreinwald.github.io/drc/reference/LL.3.md),
[`LL.4`](https://hreinwald.github.io/drc/reference/LL.4.md),
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md)

## Author

Christian Ritz

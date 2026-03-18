# Two-parameter log-logistic function

A two-parameter log-logistic function with lower limit fixed at 0 and
upper limit fixed (default 1), primarily for use with binomial/quantal
dose-response data.

## Usage

``` r
LL.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)

l2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
```

## Arguments

- upper:

  numeric value, the fixed upper limit (default 1).

- fixed:

  numeric vector of length 2, specifying fixed parameters (use NA for
  non-fixed parameters).

- names:

  character vector of length 2, specifying the names of the parameters
  (default: b, e).

- ...:

  additional arguments to
  [`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Value

See
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Details

The two-parameter log-logistic function is given by the expression
\$\$f(x) = \frac{upper}{1+\exp(b(\log(x)-\log(e)))}\$\$

## See also

[`LL.3`](https://hreinwald.github.io/drc/reference/LL.3.md),
[`LL.4`](https://hreinwald.github.io/drc/reference/LL.4.md),
[`LL.5`](https://hreinwald.github.io/drc/reference/LL.5.md),
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md)

## Author

Christian Ritz

## Examples

``` r
earthworms.m1 <- drm(number/total~dose, weights=total,
  data = earthworms, fct = LL.2(), type = "binomial")
```

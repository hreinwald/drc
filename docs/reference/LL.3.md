# Three-parameter log-logistic function

A three-parameter log-logistic function with lower limit fixed at 0.

## Usage

``` r
LL.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)

l3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 3, specifying fixed parameters (use NA for
  non-fixed parameters).

- names:

  character vector of length 3, specifying the names of the parameters
  (default: b, d, e).

- ...:

  additional arguments to
  [`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Value

See
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Details

The three-parameter log-logistic function is given by the expression
\$\$f(x) = \frac{d}{1+\exp(b(\log(x)-\log(e)))}\$\$

## See also

[`LL.2`](https://hreinwald.github.io/drc/reference/LL.2.md),
[`LL.4`](https://hreinwald.github.io/drc/reference/LL.4.md),
[`LL.5`](https://hreinwald.github.io/drc/reference/LL.5.md),
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md)

## Author

Christian Ritz

## Examples

``` r
ryegrass.model1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.3())
```

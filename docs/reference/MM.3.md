# Three-parameter Michaelis-Menten function

A three-parameter (shifted) Michaelis-Menten function where b is fixed
at -1 and f at 1.

## Usage

``` r
MM.3(fixed = c(NA, NA, NA), names = c("c", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 3, specifying fixed parameters (use NA for
  non-fixed parameters).

- names:

  character vector of length 3, specifying the names of the parameters
  (default: c, d, e).

- ...:

  additional arguments to
  [`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Value

See
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Details

The three-parameter Michaelis-Menten function is \$\$f(x) = c +
\frac{d-c}{1+(e/x)}\$\$

## See also

[`MM.2`](https://hreinwald.github.io/drc/reference/MM.2.md),
[`AR.2`](https://hreinwald.github.io/drc/reference/AR.2.md),
[`AR.3`](https://hreinwald.github.io/drc/reference/AR.3.md)

## Author

Christian Ritz

## Examples

``` r
met.mm.m1 <- drm(gain~dose, product, data = methionine, fct = MM.3())
#> Control measurements detected for level: control
```

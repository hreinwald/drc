# Two-parameter Michaelis-Menten function

A two-parameter Michaelis-Menten function where b is fixed at -1, c at
0, and f at 1. Commonly used for enzyme kinetics and weed density
studies.

## Usage

``` r
MM.2(fixed = c(NA, NA), names = c("d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 2, specifying fixed parameters (use NA for
  non-fixed parameters).

- names:

  character vector of length 2, specifying the names of the parameters
  (default: d, e).

- ...:

  additional arguments to
  [`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Value

See
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Details

The two-parameter Michaelis-Menten function is \$\$f(x) = \frac{d \cdot
x}{e + x}\$\$ which is equivalent to \\d/(1+(e/x))\\.

## See also

[`MM.3`](https://hreinwald.github.io/drc/reference/MM.3.md),
[`AR.2`](https://hreinwald.github.io/drc/reference/AR.2.md),
[`AR.3`](https://hreinwald.github.io/drc/reference/AR.3.md)

## Author

Christian Ritz

## Examples

``` r
met.mm.m1 <- drm(gain~dose, product, data = methionine, fct = MM.2())
#> Control measurements detected for level: control
```

# Three-parameter shifted asymptotic regression model

A three-parameter asymptotic regression model where `b` is fixed at 1.
The model is given by the equation \$\$f(x) = c + (d - c)(1 - \exp(-x /
e))\$\$

## Usage

``` r
AR.3(fixed = c(NA, NA, NA), names = c("c", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 3, specifying fixed parameters (use `NA` for
  parameters that should be estimated).

- names:

  character vector of length 3 giving the names of the parameters
  (default `c("c", "d", "e")`).

- ...:

  additional arguments passed to
  [`weibull2`](https://hreinwald.github.io/drc/reference/weibull2.md).

## Value

A list of class `"Weibull-2"` as returned by
[`weibull2`](https://hreinwald.github.io/drc/reference/weibull2.md).

## See also

[`AR.2`](https://hreinwald.github.io/drc/reference/AR.2.md),
[`weibull2`](https://hreinwald.github.io/drc/reference/weibull2.md),
[`EXD.3`](https://hreinwald.github.io/drc/reference/EXD.3.md)

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = AR.3())
```

# Two-parameter asymptotic regression model

A two-parameter asymptotic regression model where `b` is fixed at 1 and
the lower limit is fixed at 0. The model is given by the equation
\$\$f(x) = d \cdot (1 - \exp(-x / e))\$\$

## Usage

``` r
AR.2(fixed = c(NA, NA), names = c("d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 2, specifying fixed parameters (use `NA` for
  parameters that should be estimated).

- names:

  character vector of length 2 giving the names of the parameters
  (default `c("d", "e")`).

- ...:

  additional arguments passed to
  [`weibull2`](https://hreinwald.github.io/drc/reference/weibull2.md).

## Value

A list of class `"Weibull-2"` as returned by
[`weibull2`](https://hreinwald.github.io/drc/reference/weibull2.md).

## See also

[`AR.3`](https://hreinwald.github.io/drc/reference/AR.3.md),
[`weibull2`](https://hreinwald.github.io/drc/reference/weibull2.md),
[`EXD.2`](https://hreinwald.github.io/drc/reference/EXD.2.md)

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = AR.2())
```

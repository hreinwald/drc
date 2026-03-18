# Four-Parameter Log-Logistic Model with log(ED50) as Parameter

A four-parameter log-logistic model where the ED50 is parameterised on
the log scale. The asymmetry parameter `f` is fixed at 1. The estimated
parameters are the slope `b`, the lower limit `c`, the upper limit `d`,
and the log(ED50) `e`.

## Usage

``` r
LL2.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 4. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated.

- names:

  character vector of length 4 giving the names of the parameters `b`,
  `c`, `d`, and `e`.

- ...:

  additional arguments passed to
  [`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md).

## Value

A list of class `"llogistic"` with the nonlinear function, self-starter,
and related components.

## See also

[`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md),
[`LL2.3`](https://hreinwald.github.io/drc/reference/LL2.3.md),
[`LL2.5`](https://hreinwald.github.io/drc/reference/LL2.5.md)

## Examples

``` r
spinach.m1 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL2.4())
```

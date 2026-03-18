# Three-Parameter Log-Logistic Model with log(ED50) and Lower Limit at 0

A three-parameter log-logistic model where the lower limit is fixed at
0. The estimated parameters are the slope `b`, the upper limit `d`, and
the log(ED50) `e`.

## Usage

``` r
LL2.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 3. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated.

- names:

  character vector of length 3 giving the names of the parameters `b`,
  `d`, and `e`.

- ...:

  additional arguments passed to
  [`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md).

## Value

A list of class `"llogistic"` with the nonlinear function, self-starter,
and related components.

## See also

[`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md),
[`LL2.2`](https://hreinwald.github.io/drc/reference/LL2.2.md),
[`LL2.4`](https://hreinwald.github.io/drc/reference/LL2.4.md)

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.3())
```

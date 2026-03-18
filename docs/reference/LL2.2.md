# Two-Parameter Log-Logistic Model with log(ED50) as Parameter

A two-parameter log-logistic model where the lower limit is fixed at 0
and the upper limit is fixed at a specified value (default 1). The
estimated parameters are the slope `b` and the log(ED50) `e`.

## Usage

``` r
LL2.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
```

## Arguments

- upper:

  numeric value giving the fixed upper limit. Defaults to 1.

- fixed:

  numeric vector of length 2. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated.

- names:

  character vector of length 2 giving the names of the parameters `b`
  and `e`.

- ...:

  additional arguments passed to
  [`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md).

## Value

A list of class `"llogistic"` with the nonlinear function, self-starter,
and related components.

## See also

[`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md),
[`LL2.3`](https://hreinwald.github.io/drc/reference/LL2.3.md),
[`LL2.4`](https://hreinwald.github.io/drc/reference/LL2.4.md),
[`LL2.5`](https://hreinwald.github.io/drc/reference/LL2.5.md)

## Examples

``` r
earthworms.m1 <- drm(number/total ~ dose, weights = total,
  data = earthworms, fct = LL2.2(), type = "binomial")
```

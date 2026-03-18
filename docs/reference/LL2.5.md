# Five-Parameter Generalised Log-Logistic Model with log(ED50) as Parameter

A five-parameter generalised log-logistic model where the ED50 is
parameterised on the log scale. All five parameters (`b`, `c`, `d`, `e`,
`f`) are estimated.

## Usage

``` r
LL2.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
```

## Arguments

- fixed:

  numeric vector of length 5. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated.

- names:

  character vector of length 5 giving the names of the parameters `b`,
  `c`, `d`, `e`, and `f`.

- ...:

  additional arguments passed to
  [`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md).

## Value

A list of class `"llogistic"` with the nonlinear function, self-starter,
and related components.

## See also

[`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md),
[`LL2.3`](https://hreinwald.github.io/drc/reference/LL2.3.md),
[`LL2.4`](https://hreinwald.github.io/drc/reference/LL2.4.md)

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.5())
```

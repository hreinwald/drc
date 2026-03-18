# Three-Parameter Log-Logistic Model with log(ED50) and Fixed Upper Limit

A three-parameter log-logistic model where the upper limit is fixed at a
specified value (default 1). The estimated parameters are the slope `b`,
the lower limit `c`, and the log(ED50) `e`.

## Usage

``` r
LL2.3u(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
```

## Arguments

- upper:

  numeric value giving the fixed upper limit. Defaults to 1.

- fixed:

  numeric vector of length 3. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated.

- names:

  character vector of length 3 giving the names of the parameters `b`,
  `c`, and `e`.

- ...:

  additional arguments passed to
  [`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md).

## Value

A list of class `"llogistic"` with the nonlinear function, self-starter,
and related components.

## See also

[`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md),
[`LL2.2`](https://hreinwald.github.io/drc/reference/LL2.2.md),
[`LL2.3`](https://hreinwald.github.io/drc/reference/LL2.3.md)

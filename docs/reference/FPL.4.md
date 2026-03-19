# Four-parameter fractional polynomial-logistic model

Convenience function for the four-parameter fractional
polynomial-logistic model.

## Usage

``` r
FPL.4(p1, p2, fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
```

## Arguments

- p1:

  numeric denoting the negative power of log(dose+1) in the fractional
  polynomial.

- p2:

  numeric denoting the positive power of log(dose+1) in the fractional
  polynomial.

- fixed:

  numeric vector of length 4 specifying fixed parameters (NAs for free
  parameters).

- names:

  character vector of parameter names.

- ...:

  additional arguments passed to
  [`fplogistic`](https://hreinwald.github.io/drc/reference/fplogistic.md).

## Value

A list (see
[`fplogistic`](https://hreinwald.github.io/drc/reference/fplogistic.md)).

## See also

[`fplogistic`](https://hreinwald.github.io/drc/reference/fplogistic.md),
[`maED`](https://hreinwald.github.io/drc/reference/maED.md)

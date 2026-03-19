# Four-parameter log-normal dose-response model

`LN.4` is a convenience function for the full four-parameter log-normal
model.

## Usage

``` r
LN.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 4 specifying fixed parameters (NAs for free
  parameters).

- names:

  character vector of parameter names.

- ...:

  additional arguments passed to
  [`lnormal`](https://hreinwald.github.io/drc/reference/lnormal.md).

## Value

A list (see
[`lnormal`](https://hreinwald.github.io/drc/reference/lnormal.md)).

## See also

[`lnormal`](https://hreinwald.github.io/drc/reference/lnormal.md),
[`LN.2`](https://hreinwald.github.io/drc/reference/LN.2.md),
[`LN.3`](https://hreinwald.github.io/drc/reference/LN.3.md)

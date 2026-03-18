# Three-parameter log-normal dose-response model

`LN.3` is a convenience function for the log-normal model with the lower
limit fixed at 0.

## Usage

``` r
LN.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 3 specifying fixed parameters (NAs for free
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
[`LN.4`](https://hreinwald.github.io/drc/reference/LN.4.md)

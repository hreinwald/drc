# Two-parameter log-normal dose-response model

`LN.2` is a convenience function for the log-normal model with lower
limit fixed at 0 and upper limit fixed (default 1), corresponding to the
classic probit model.

## Usage

``` r
LN.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
```

## Arguments

- upper:

  numeric specifying the fixed upper horizontal asymptote. Default is 1.

- fixed:

  numeric vector of length 2 specifying fixed parameters (NAs for free
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
[`LN.3`](https://hreinwald.github.io/drc/reference/LN.3.md),
[`LN.4`](https://hreinwald.github.io/drc/reference/LN.4.md)

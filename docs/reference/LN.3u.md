# Three-parameter log-normal model with upper limit fixed

`LN.3u` is a convenience function for the log-normal model with the
upper limit fixed (default 1).

## Usage

``` r
LN.3u(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
```

## Arguments

- upper:

  numeric specifying the fixed upper horizontal asymptote. Default is 1.

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
[`LN.3`](https://hreinwald.github.io/drc/reference/LN.3.md),
[`LN.4`](https://hreinwald.github.io/drc/reference/LN.4.md)

# Two-parameter Gompertz model

Convenience function for the Gompertz model with lower limit fixed at 0
and upper limit fixed.

## Usage

``` r
G.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
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
  [`gompertz`](https://hreinwald.github.io/drc/reference/gompertz.md).

## Value

A list (see
[`gompertz`](https://hreinwald.github.io/drc/reference/gompertz.md)).

## See also

[`gompertz`](https://hreinwald.github.io/drc/reference/gompertz.md),
[`G.3`](https://hreinwald.github.io/drc/reference/G.3.md),
[`G.4`](https://hreinwald.github.io/drc/reference/G.4.md)

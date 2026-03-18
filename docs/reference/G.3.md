# Three-parameter Gompertz model

Convenience function for the Gompertz model with the lower limit fixed
at 0.

## Usage

``` r
G.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 3 specifying fixed parameters (NAs for free
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
[`G.2`](https://hreinwald.github.io/drc/reference/G.2.md),
[`G.4`](https://hreinwald.github.io/drc/reference/G.4.md)

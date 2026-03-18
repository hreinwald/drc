# Three-parameter NEC model

Convenience function for the NEC model with the lower limit fixed at 0.

## Usage

``` r
NEC.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 3 specifying fixed parameters (NAs for free
  parameters).

- names:

  character vector of parameter names.

- ...:

  additional arguments passed to
  [`NEC`](https://hreinwald.github.io/drc/reference/NEC.md).

## Value

A list (see [`NEC`](https://hreinwald.github.io/drc/reference/NEC.md)).

## See also

[`NEC`](https://hreinwald.github.io/drc/reference/NEC.md),
[`NEC.2`](https://hreinwald.github.io/drc/reference/NEC.2.md),
[`NEC.4`](https://hreinwald.github.io/drc/reference/NEC.4.md)

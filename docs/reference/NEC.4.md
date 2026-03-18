# Four-parameter NEC model

Convenience function for the full four-parameter NEC model.

## Usage

``` r
NEC.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 4 specifying fixed parameters (NAs for free
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
[`NEC.3`](https://hreinwald.github.io/drc/reference/NEC.3.md)

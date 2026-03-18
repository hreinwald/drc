# Four-parameter Gompertz model

Convenience function for the full four-parameter Gompertz model.

## Usage

``` r
G.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 4 specifying fixed parameters (NAs for free
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
[`G.3`](https://hreinwald.github.io/drc/reference/G.3.md)

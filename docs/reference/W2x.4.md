# Four-parameter Weibull type 2 model with lag time

A four-parameter Weibull type 2 model with lag time, where `b` is fixed
at 1. This is a convenience wrapper around
[`weibull2x`](https://hreinwald.github.io/drc/reference/weibull2x.md).

## Usage

``` r
W2x.4(fixed = c(NA, NA, NA, NA), names = c("c", "d", "e", "t0"), ...)
```

## Arguments

- fixed:

  numeric vector of length 4. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated
  (default is `c(NA, NA, NA, NA)`).

- names:

  character vector of length 4 giving the names of the parameters
  (default is `c("c", "d", "e", "t0")`).

- ...:

  additional arguments passed to
  [`weibull2x`](https://hreinwald.github.io/drc/reference/weibull2x.md).

## Value

A list of class `"Weibull-2"` containing the nonlinear function, self
starter function, and parameter names.

## See also

[`weibull2x`](https://hreinwald.github.io/drc/reference/weibull2x.md),
[`W2x.3`](https://hreinwald.github.io/drc/reference/W2x.3.md),
[`W2.4`](https://hreinwald.github.io/drc/reference/W2.4.md)

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W2x.4())
#> Warning: NaNs produced
#> Warning: NaNs produced
```

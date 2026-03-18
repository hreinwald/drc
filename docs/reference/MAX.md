# Maximum mean response

Estimates the maximum mean response and the dose at which it occurs.
This function is only implemented for the built-in functions of class
[`braincousens`](https://hreinwald.github.io/drc/reference/braincousens.md)
and
[`cedergreen`](https://hreinwald.github.io/drc/reference/cedergreen.md).

## Usage

``` r
MAX(object, lower = 0.001, upper = 1000, pool = TRUE)
```

## Arguments

- object:

  an object of class 'drc'.

- lower:

  numeric. Lower limit for bisection method. Must be smaller than the
  EDx level to be calculated.

- upper:

  numeric. Upper limit for bisection method. Must be larger than the EDx
  level to be calculated.

- pool:

  logical. If TRUE curves are pooled. Otherwise they are not. This
  argument only works for models with independently fitted curves as
  specified in
  [`drm`](https://hreinwald.github.io/drc/reference/drm.md).

## Value

A matrix with one row per curve in the data set and two columns: one
containing the dose at which the maximum occurs and one containing the
corresponding maximum response.

## References

Cedergreen, N. and Ritz, C. and Streibig, J. C. (2005) Improved
empirical models describing hormesis, *Environmental Toxicology and
Chemistry* **24**, 3166–3172.

## Author

Christian Ritz

## Examples

``` r
## Fitting a Cedergreen-Ritz-Streibig model
lettuce.m1 <- drm(weight~conc, data = lettuce, fct = CRS.4c())

## Finding maximum average response and the corresponding dose
MAX(lettuce.m1)
#>      Dose Response
#> 1 0.25587    1.178
```

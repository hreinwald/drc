# Generalised Cedergreen-Ritz-Streibig Model for Hormesis

A six-parameter extension of the Cedergreen-Ritz-Streibig model for
describing hormesis, where the alpha parameter is estimated rather than
fixed.

## Usage

``` r
CRS.6(
  fixed = c(NA, NA, NA, NA, NA, NA),
  names = c("b", "c", "d", "e", "f", "g"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL
)
```

## Arguments

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":").

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used (optional).

## Value

A list containing the nonlinear model function, the self starter
function, and the parameter names.

## Details

The model function is:

\$\$f(x) = c + \frac{d-c+f \exp(-1/x^g)}{1+\exp(b(\log(x)-\log(e)))}\$\$

This generalises the five-parameter
[`CRS.5a`](https://hreinwald.github.io/drc/reference/CRS.5a.md) model by
estimating the alpha exponent (parameter \\g\\) instead of fixing it.

## Note

This function is for use with
[`drm`](https://hreinwald.github.io/drc/reference/drm.md).

## See also

[`CRS.5a`](https://hreinwald.github.io/drc/reference/CRS.5a.md),
[`cedergreen`](https://hreinwald.github.io/drc/reference/cedergreen.md)

## Author

Christian Ritz
